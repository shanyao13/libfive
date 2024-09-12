/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2017  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <stack>
#include <boost/lockfree/stack.hpp>

#include "libfive/render/brep/per_thread_brep.hpp"
#include "libfive/render/brep/free_thread_handler.hpp"
#include "libfive/render/brep/settings.hpp"
#include "libfive/render/brep/mesh.hpp"
#include "libfive/render/brep/root.hpp"

#include "libfive/render/axes.hpp"
#include "libfive/eval/interval.hpp"

namespace libfive {

/*
 *  Class to walk a dual grid for a quad or octree
 */
template <unsigned N>
class Dual
{
public:
     /*
      *  Basic dual-walking function
      *
      *  The mesher type M needs
      *     load(const std::array<T*, N>& trees)
      *     Input (typename)
      *     Output (typename)
      *  and must have a constructor of the form
      *     M(PerThreadBRep<N>&, A...)
      */
    template<typename M, typename ... A>
    static std::unique_ptr<typename M::Output> walk(
            const Root<typename M::Input>& t,
            const BRepSettings& settings,
            A... args);

     /*
      *  Flexible dual-walking function
      *
      *  The mesher type M needs
      *     load(const std::array<T*, N>& trees)
      *     Input (typename)
      *     Output (typename)
      *
      *  The factory can be anything that spits out valid M objects,
      *  given a PerThreadBRep and worker index.
      */
    template<typename M>
    static std::unique_ptr<typename M::Output> walk_(
            const Root<typename M::Input>& t,
            const BRepSettings& settings,
            std::function<M(PerThreadBRep<N>&, int)> MesherFactory);

protected:
    template<typename T, typename Mesher>
    static void run(Mesher& m,
                    boost::lockfree::stack<const T*,
                                           boost::lockfree::fixed_sized<true>>& tasks,
                    const BRepSettings& settings,
                    std::atomic_bool& done);

    template <typename T, typename Mesher>
    static void work(const T* t, Mesher& m);

    template <typename T, typename Mesher>
    static void handleTopEdges(T* t, Mesher& m);
};

////////////////////////////////////////////////////////////////////////////////
// 2D Implementation
template <typename T, typename V, Axis::Axis A>
void edge2(const std::array<const T*, 2>& ts, V& v)
{
    constexpr uint8_t perp = (Axis::X | Axis::Y) ^ A;

    if (std::any_of(ts.begin(), ts.end(),
        [](const T* t){ return t->isBranch(); }))
    {
        edge2<T, V, A>({{ts[0]->child(perp), ts[1]->child(0)}}, v);
        edge2<T, V, A>({{ts[0]->child(A|perp), ts[1]->child(A)}}, v);
    }
    else if (std::all_of(ts.begin(), ts.end(),
        [](const T* t){ return t->type == Interval::AMBIGUOUS &&
                               !t->isBranch(); }))
    {
        v.template load<A>(ts);
    }
}

template <>
template <typename T, typename V>
void Dual<2>::work(const T* t, V& v)
{
    edge2<T, V, Axis::Y>({{t->child(0), t->child(Axis::X)}}, v);
    edge2<T, V, Axis::Y>({{t->child(Axis::Y), t->child(Axis::Y | Axis::X)}}, v);
    edge2<T, V, Axis::X>({{t->child(0), t->child(Axis::Y)}}, v);
    edge2<T, V, Axis::X>({{t->child(Axis::X), t->child(Axis::X | Axis::Y)}}, v);
}

template <>
template <typename T, typename Mesher>
void Dual<2>::handleTopEdges(T* t, Mesher& m)
{
    (void)t;
    (void)m;

    // TODO
    // No one should be calling this yet, because simplex meshing
    // isn't implemented in 2D.
    assert(false);
}

////////////////////////////////////////////////////////////////////////////////
// 3D Implementation
    //edge3 函数是 DC (Dual Contouring) 算法的一部分，专门用于处理 3D 网格中的边缘生成。此函数采用递归方式遍历并处理八叉树结构的节点，根据需要生成边缘数据。
    //这种边缘数据是 3D 网格的重要组成部分，决定了网格的几何形态
    //edge3 函数是一个模板函数，用于递归地处理网格的边缘，生成相应的几何信息。
template <typename T, typename V, Axis::Axis A>
void edge3(const std::array<const T*, 4> ts, V& v)
{
    constexpr auto Q = Axis::Q(A);
    constexpr auto R = Axis::R(A);

    //edge3 函数的主要功能是递归遍历给定的四个节点（ts），并根据节点是否为分支节点决定下一步操作。
    //如果遇到分支节点（非叶子节点），函数会继续递归处理子节点，直到处理到叶子节点，再由 v 对象加载边缘数据
    if (std::any_of(ts.begin(), ts.end(),
        [](const T* t){ return t->isBranch(); }))
    {
        edge3<T, V, A>({{ts[0]->child(Q|R), ts[1]->child(R), ts[2]->child(Q), ts[3]->child(0)}}, v);
        edge3<T, V, A>({{ts[0]->child(Q|R|A), ts[1]->child(R|A), ts[2]->child(Q|A), ts[3]->child(A)}}, v);
    }
    else
    {
        //ts 长度为 4 的节点数组，包含四个树节点，代表当前边缘的四个顶点。
        v.template load<A>(ts);
    }
}

    //用于处理 3D 网格的面生成
    //在 3D 网格处理中，面是构成体的基本单元，而 face3 函数递归地遍历树结构，确保面和边被正确地处理和生成
template <typename T, typename V, Axis::Axis A>
void face3(const std::array<const T*, 2> ts, V& v)
{
    if (std::any_of(ts.begin(), ts.end(),
        [](const T* t){ return t->isBranch(); }))
    {
        constexpr auto Q = Axis::Q(A);
        constexpr auto R = Axis::R(A);

        for (unsigned k : {0, (int)Q, (int)R, Q|R})
        {
            face3<T, V, A>({{ts[0]->child(k|A), ts[1]->child(k)}}, v);
        }

        edge3<T, V, Q>({{ts[0]->child(A), ts[0]->child(R|A), ts[1]->child(0), ts[1]->child(R)}}, v);
        edge3<T, V, Q>({{ts[0]->child(Q|A), ts[0]->child(Q|R|A), ts[1]->child(Q), ts[1]->child(Q|R)}}, v);

        edge3<T, V, R>({{ts[0]->child(A), ts[1]->child(0), ts[0]->child(A|Q), ts[1]->child(Q)}}, v);
        edge3<T, V, R>({{ts[0]->child(R|A), ts[1]->child(R), ts[0]->child(R|A|Q), ts[1]->child(R|Q)}}, v);
    }
}

template <typename T, typename V, Axis::Axis A>
void call_edge3(const T* t, V& v)
{
    for (auto a : {Axis::Axis(0), A})
    {
        edge3<T, V, A>({{t->child(a),
             t->child(Axis::Q(A) | a),
             t->child(Axis::R(A) | a),
             t->child(Axis::Q(A) | Axis::R(A) | a)}}, v);
    }
}

template <typename T, typename V, Axis::Axis A>
void call_face3(const T* t, V& v)
{
    constexpr auto q = Axis::Q(A);
    constexpr auto r = Axis::R(A);

    face3<T, V, A>({{t->child(0), t->child(A)}}, v);
    face3<T, V, A>({{t->child(q), t->child(q|A)}}, v);
    face3<T, V, A>({{t->child(r), t->child(r|A)}}, v);
    face3<T, V, A>({{t->child(q|r), t->child(q|r|A)}}, v);
}

template <>
template <typename T, typename V>
void Dual<3>::work(const T* t, V& v)
{
    // Call the face procedure on every pair of cells (4x per axis)
    //call_face3 函数用于处理每个坐标轴的面。对于三维空间，每个坐标轴都需要处理 4 个相邻的面，这里分别调用 Axis::X、Axis::Y 和 Axis::Z 方向的处理函数。
    //call_face3 函数的具体实现会根据当前的坐标轴 (Axis::X, Axis::Y, Axis::Z) 来处理面，并生成相应的网格数据
    call_face3<T, V, Axis::X>(t, v);
    call_face3<T, V, Axis::Y>(t, v);
    call_face3<T, V, Axis::Z>(t, v);

    // Call the edge function 6 times (2x per axis)
    //call_edge3 函数用于处理每个坐标轴的边缘。每个坐标轴需要处理 2 个相邻的边缘，分别调用 Axis::X、Axis::Y 和 Axis::Z 方向的处理函数。
    //call_edge3 函数的具体实现会根据当前的坐标轴 (Axis::X, Axis::Y, Axis::Z) 来处理边缘，并生成相应的网格数据
    call_edge3<T, V, Axis::X>(t, v);
    call_edge3<T, V, Axis::Y>(t, v);
    call_edge3<T, V, Axis::Z>(t, v);
}

template <>
template <typename T, typename V>
void Dual<3>::handleTopEdges(T* t, V& v)
{
    auto e = T::empty();

    for (unsigned i=0; i < 4; ++i)
    {
        std::array<T*, 4> ts = {{e.get(), e.get(), e.get(), e.get()}};
        ts[i] = t;
        edge3<T, V, Axis::X>(ts, v);
        edge3<T, V, Axis::Y>(ts, v);
        edge3<T, V, Axis::Z>(ts, v);
    }

    for (unsigned i=0; i < 2; ++i)
    {
        std::array<T*, 2> ts = {{e.get(), e.get()}};
        ts[i] = t;
        face3<T, V, Axis::X>(ts, v);
        face3<T, V, Axis::Y>(ts, v);
        face3<T, V, Axis::Z>(ts, v);
    }
}

////////////////////////////////////////////////////////////////////////////////
//        out = Dual<3>::walk<DCMesher>(t, settings); M:DCMesher
template <unsigned N>
template<typename M, typename ... A>
std::unique_ptr<typename M::Output> Dual<N>::walk(
            const Root<typename M::Input>& t,
            const BRepSettings& settings,
            A... args)
{
    //[&args...](PerThreadBRep<N>& brep, int i)：这是一个lambda表达式，它捕获了参数包args（通过引用捕获），
    // 并接受两个参数：PerThreadBRep<N>& brep和int i。PerThreadBRep<N>是一个与线程相关的B-rep（边界表示）数据结构，i是一个整数索引
    return walk_<M>(
            t, settings,
            [&args...](PerThreadBRep<N>& brep, int i) {
                (void)i;
                return M(brep, args...);
                });

}


//这段代码是 Dual<N>::walk_ 函数的实现，主要功能是并行处理三维（或更高维）数据，利用多线程和任务队列生成几何体或网格。
// 这是一个通用的并行计算框架，结合了模板编程、锁自由数据结构、多线程和回调工厂函数来完成任务
template <unsigned N>
template<typename M>
std::unique_ptr<typename M::Output> Dual<N>::walk_(
            const Root<typename M::Input>& t,
            const BRepSettings& settings,
            std::function<M(PerThreadBRep<N>&, int)> MesherFactory)
{
    boost::lockfree::stack<
        const typename M::Input*,
        boost::lockfree::fixed_sized<true>> tasks(settings.workers);
    tasks.push(t.get());
    t->resetPending();

    //线程设置和进度初始化:
    //
    //std::atomic<uint32_t> global_index(1): 一个全局索引，用于多线程之间的同步。
    //std::vector<PerThreadBRep<N>> breps: 每个线程的网格生成器数据结构，用于存储生成的部分网格。
    //settings.progress_handler->nextPhase(t.size() + 1): 如果有进度处理器，设置进度阶段。
    std::atomic<uint32_t> global_index(1);
    std::vector<PerThreadBRep<N>> breps;
    for (unsigned i=0; i < settings.workers; ++i) {
        breps.emplace_back(PerThreadBRep<N>(global_index));
    }

    if (settings.progress_handler) {
        settings.progress_handler->nextPhase(t.size() + 1);
    }

    //启动多线程任务:
    //
    //std::vector<std::future<void>> futures: 存储每个线程的 future 对象。
    //std::atomic_bool done(false): 用于指示所有线程的完成状态。
    //使用 std::async 启动多个线程，每个线程创建一个 Mesher 实例并运行 Dual<N>::run 方法来处理任务
    std::vector<std::future<void>> futures;
    futures.resize(settings.workers);
    std::atomic_bool done(false);
    for (unsigned i=0; i < settings.workers; ++i) {
        futures[i] = std::async(std::launch::async,
            [&breps, &tasks, &MesherFactory, &settings, &done, i]()
            {
                //auto m = MesherFactory(breps[i], i);:
                //
                //使用 MesherFactory 工厂函数为当前线程创建一个网格生成器对象 m。
                //breps[i] 是当前线程的局部网格生成数据结构，i 是当前线程的索引。
                //Dual<N>::run(m, tasks, settings, done);:
                //
                //调用 Dual<N>::run 函数，该函数使用 Mesher 对象 m 处理任务队列 tasks。
                //settings 包含任务的相关配置，done 用于线程之间的同步，标识任务的完成状态。
                auto m = MesherFactory(breps[i], i);
                Dual<N>::run(m, tasks, settings, done);
            });
    }

    // Wait on all of the futures
    for (auto& f : futures) {
        f.get();
    }

    assert(done.load() || settings.cancel.load());

    // Handle the top tree edges (only used for simplex meshing)
    if (M::needsTopEdges()) {
        auto m = MesherFactory(breps[0], 0);
        Dual<N>::handleTopEdges(t.get(), m);
    }

    auto out = std::make_unique<typename M::Output>();
    out->collect(breps);
    return out;
}


template <unsigned N>
template <typename T, typename V>
void Dual<N>::run(V& v,
                  boost::lockfree::stack<const T*,
                                         boost::lockfree::fixed_sized<true>>& tasks,
                  const BRepSettings& settings,
                  std::atomic_bool& done)

{
    // Tasks to be evaluated by this thread (populated when the
    // MPMC stack is completely full).
    std::stack<const T*, std::vector<const T*>> local;

    while (!done.load() && !settings.cancel.load())
    {
        // Prioritize picking up a local task before going to
        // the MPMC queue, to keep things in this thread for
        // as long as possible.
        const T* t;
        if (local.size())
        {
            t = local.top();
            local.pop();
        }
        else if (!tasks.pop(t))
        {
            t = nullptr;
        }

        // If we failed to get a task, keep looping
        // (so that we terminate when either of the flags are set).
        if (t == nullptr)
        {
            if (settings.free_thread_handler != nullptr) {
                settings.free_thread_handler->offerWait();
            }
            continue;
        }

        if (t->isBranch())
        {
            // Recurse, calling the cell procedure for every child
            for (const auto& c_ : t->children)
            {
                const auto c = c_.load();
                if (!tasks.bounded_push(c)) {
                    local.push(c);
                }
            }
            continue;
        }

        // Special-case for singleton trees, which have null parents
        // (and have already been subtracted from pending)
        if (T::isSingleton(t)) {
            continue;
        }

        if (settings.progress_handler) {
            settings.progress_handler->tick();
        }

        for (t = t->parent; t && t->pending-- == 0; t = t->parent)
        {
            // Do the actual DC work (specialized for N = 2 or 3)
            // v：DCMesher
            Dual<N>::work(t, v);

            // Report trees as completed
            if (settings.progress_handler) {
                settings.progress_handler->tick();
            }
        }

        // Termination condition:  if we've ended up pointing at the parent
        // of the tree's root (which is nullptr), then we're done and break
        if (t == nullptr) {
            break;
        }
    }

    // If we've broken out of the loop, then we should set the done flag
    // so that other worker threads also terminate.
    done.store(true);
}

}   // namespace libfive
