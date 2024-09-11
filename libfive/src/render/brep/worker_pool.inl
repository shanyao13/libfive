/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "libfive/render/brep/free_thread_handler.hpp"
#include "libfive/render/brep/settings.hpp"
#include "libfive/render/brep/worker_pool.hpp"
#include "libfive/render/brep/vol/vol_tree.hpp"
#include "libfive/eval/evaluator.hpp"

namespace libfive {

template <typename T, typename Neighbors, unsigned N>
Root<T> WorkerPool<T, Neighbors, N>::build(
        const Tree& t_, const Region<N>& region_,
        const BRepSettings& settings)
{
    // Build evaluators for the pool
    const auto t = t_.optimized();
    std::vector<Evaluator, Eigen::aligned_allocator<Evaluator>> es;
    es.reserve(settings.workers);
    for (unsigned i=0; i < settings.workers; ++i) {
        es.emplace_back(Evaluator(t));
    }
    return build(es.data(), region_, settings);
}

//定义模版，允许该类或函数与不同的类型和参数一起使用。T 树结点类型；neighbors 邻居关系类型，用于管理树节点之间的关系； N 维度，通常为2或3。
// Root<T> 返回类型； 函数名和所属类
template <typename T, typename Neighbors, unsigned N>
Root<T> WorkerPool<T, Neighbors, N>::build(
        Evaluator* eval, const Region<N>& region_,
        const BRepSettings& settings)
{
    if (settings.vol && !settings.vol->contains(region_)) {
        std::cerr << "WorkerPool::build: Invalid region for vol tree\n";
    }

    //调整区域层级
    // region_.withResolution 主要功能是根据给定的最小特征尺寸 (min_feature) 调整区域的层级 (level)
    const auto region = region_.withResolution(settings.min_feature);
    // 创建根节点
    // 使用调整后的区域 region 创建了一个树的根节点 root。此节点通常是递归树结构（如四叉树、八叉树）的起点，为后续的计算和细分操作提供基础
    auto root(new T(nullptr, 0, region));

    //任务初始化
    //创建一个无锁栈 LockFreeStack tasks(settings.workers);
    // 用于存储任务，并将第一个任务压入栈中：tasks.push({root, eval->getDeck()->tape, Neighbors(), settings.vol});。这个任务包括了根节点、评估器的计算图、邻居信息和体积设置
    LockFreeStack tasks(settings.workers);
    tasks.push({root, eval->getDeck()->tape, Neighbors(), settings.vol});

    //多线程任务调度：
    //
    //初始化一个 std::vector<std::future<void>> futures 向量来存储每个线程的执行结果。
    //std::mutex root_lock; 用于线程同步，保护对根节点的访问。
    //std::atomic_bool done(false); 用于指示任务的完成状态。
    std::vector<std::future<void>> futures;
    futures.resize(settings.workers);

    Root<T> out(root);
    std::mutex root_lock;

    // Kick off the progress tracking thread, based on the number of
    // octree levels and a fixed split per level
    //计算一个与树结构层级相关的进度计数 ticks，并用它来启动或更新进度跟踪器 progress_handler。
    // 这种设计让程序在处理递归树（如四叉树或八叉树）的构建过程中，可以实时跟踪进度，并向用户或系统反馈当前的工作状态，有助于管理长时间运行的任务
    uint64_t ticks = 0;
    for (int i=0; i < region.level; ++i) {
        ticks = (ticks + 1) * (1 << N);
    }
    if (settings.progress_handler) {
        settings.progress_handler->nextPhase(ticks + 1);
    }

    //启动线程执行任务：
    //
    //启动多线程执行任务，std::async(std::launch::async, ...) 用于异步执行 run 函数。
    //传递给 run 函数的参数包括评估器、任务栈、输出树、锁、设置以及完成标志。
    std::atomic_bool done(false);
    for (unsigned i=0; i < settings.workers; ++i)
    {
        //启动多个线程并行执行 run 函数，以便在多核环境中高效处理任务。lambda表达式，Lambda 函数是每个线程执行的任务。
        // 通过 & 符号按引用捕获多个外部变量，使得线程在执行时能够访问和修改这些变量
        futures[i] = std::async(std::launch::async,
                [&eval, &tasks, &out, &root_lock, &settings, &done, i](){
                    run(eval + i, tasks, out, root_lock, settings, done);
                });
    }

    // Wait on all of the futures
    for (auto& f : futures)
    {
        f.get();
    }

    assert(done.load() || settings.cancel.load());

    if (settings.cancel.load())
    {
        return Root<T>();
    }
    else
    {
        return out;
    }
}

//这段代码实现了 WorkerPool 类的 run 方法，其主要功能是通过多线程处理任务，以递归方式构建树（如四叉树或八叉树）
template <typename T, typename Neighbors, unsigned N>
void WorkerPool<T, Neighbors, N>::run(
        Evaluator* eval, LockFreeStack& tasks,
        Root<T>& root, std::mutex& root_lock,
        const BRepSettings& settings,
        std::atomic_bool& done)
{
    // Tasks to be evaluated by this thread (populated when the
    // MPMC stack is completely full).
    std::stack<Task, std::vector<Task>> local;

    typename T::Pool object_pool;

    while (!done.load() && !settings.cancel.load())
    {
        // Prioritize picking up a local task before going to
        // the MPMC queue, to keep things in this thread for
        // as long as possible.
        Task task;
        if (local.size())
        {
            task = local.top();
            local.pop();
        }
        else if (!tasks.pop(task))
        {
            task.target = nullptr;
        }

        // If we failed to get a task, keep looping
        // (so that we terminate when either of the flags are set).
        if (task.target == nullptr)
        {
            if (settings.free_thread_handler != nullptr) {
                settings.free_thread_handler->offerWait();
            }
            continue;
        }

        auto tape = task.tape;
        auto t = task.target;

        // Find our local neighbors.  We do this at the last minute to
        // give other threads the chance to populate more pointers.
        Neighbors neighbors;
        if (t->parent)
        };
        up();
        while (t != nullptr && t->collectChildren(eval, tape,
                                                  object_pool,
                                                  settings.max_err))
        {
            // Report the volume of completed trees as we walk back
            // up towards the root of the tree.
            if (settings.progress_handler) {
                settings.progress_handler->tick();
            }
            up();
        }

        // Termination condition:  if we've ended up pointing at the parent
        // of the tree's root (which is nullptr), then we're done and break
        if (t == nullptr)
        {
            break;
        }
    }

    // If we've broken out of the loop, then we should set the done flag
    // so that other worker threads also terminate.
    done.store(true);

    {   // Release the pooled objects to the root
        std::lock_guard<std::mutex> lock(root_lock);
        root.claim(object_pool);
    }
}
        {
            neighbors = task.parent_neighbors.push(
                t->parent_index, t->parent->children);
        }

        // If this tree is larger than the minimum size, then it will either
        // be unambiguously filled/empty, or we'll need to recurse.
        const bool can_subdivide = t->region.level > 0;
        if (can_subdivide)
        {
            Tape::Handle next_tape;
            if (task.vol) {
                auto i = task.vol->check(t->region);
                if (i == Interval::EMPTY || i == Interval::FILLED) {
                    t->setType(i);
                }
            }
            if (t->type == Interval::UNKNOWN) {
                next_tape = t->evalInterval(eval, task.tape, object_pool);
            }
            if (next_tape != nullptr) {
                tape = next_tape;
            }

            // If this Tree is ambiguous, then push the children to the stack
            // and keep going (because all the useful work will be done
            // by collectChildren eventually).
            assert(t->type != Interval::UNKNOWN);
            if (t->type == Interval::AMBIGUOUS)
            {
                auto rs = t->region.subdivide();
                for (unsigned i=0; i < t->children.size(); ++i)
                {
                    // If there are available slots, then pass this work
                    // to the queue; otherwise, undo the decrement and
                    // assign it to be evaluated locally.
                    auto next_tree = object_pool.get(t, i, rs[i]);
                    auto next_vol = task.vol ? task.vol->push(i, rs[i].perp)
                                             : nullptr;
                    Task next{next_tree, tape, neighbors, next_vol};
                    if (!tasks.bounded_push(next))
                    {
                        local.push(next);
                    }
                }

                // If we did an interval evaluation, then we either
                // (a) are done with this tree because it is empty / filled
                // (b) don't do anything until all of its children are done
                //
                // In both cases, we should keep looping; the latter case
                // is handled in collectChildren below.
                continue;
            }
        }
        else
        {
            t->evalLeaf(eval, tape, object_pool, neighbors);
        }

        if (settings.progress_handler)
        {
            if (can_subdivide)
            {
                // Accumulate all of the child XTree cells that would have been
                // included if we continued to subdivide this tree, then pass
                // all of them to the progress tracker
                uint64_t ticks = 0;
                for (int i=0; i < t->region.level; ++i) {
                    ticks = (ticks + 1) * (1 << N);
                }
                settings.progress_handler->tick(ticks + 1);
            }
            else
            {
                settings.progress_handler->tick(1);
            }
        }

        // If all of the children are done, then ask the parent to collect them
        // (recursively, merging the trees on the way up, and reporting
        // completed tree cells to the progress tracker if present).
        auto up = [&]{
            t = t->parent;
            if (t) {
                tape = tape->getBase(t->region.region3());
            }
        };
        up();
        while (t != nullptr && t->collectChildren(eval, tape,
                                                  object_pool,
                                                  settings.max_err))
        {
            // Report the volume of completed trees as we walk back
            // up towards the root of the tree.
            if (settings.progress_handler) {
                settings.progress_handler->tick();
            }
            up();
        }

        // Termination condition:  if we've ended up pointing at the parent
        // of the tree's root (which is nullptr), then we're done and break
        if (t == nullptr)
        {
            break;
        }
    }

    // If we've broken out of the loop, then we should set the done flag
    // so that other worker threads also terminate.
    done.store(true);

    {   // Release the pooled objects to the root
        std::lock_guard<std::mutex> lock(root_lock);
        root.claim(object_pool);
    }
}

}   // namespace libfive
