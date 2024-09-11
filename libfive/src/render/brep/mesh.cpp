/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2017  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include <numeric>
#include <fstream>
#include <boost/algorithm/string/predicate.hpp>

#include "libfive/eval/evaluator.hpp"

#include "libfive/render/brep/mesh.hpp"
#include "libfive/render/brep/dual.hpp"
#include "libfive/render/brep/region.hpp"
#include "libfive/render/brep/settings.hpp"

// Dual contouring
#include "libfive/render/brep/dc/dc_worker_pool.hpp"
#include "libfive/render/brep/dc/dc_mesher.hpp"

// Simplex meshing
#include "libfive/render/brep/simplex/simplex_worker_pool.hpp"
#include "libfive/render/brep/simplex/simplex_mesher.hpp"

// Hybrid meshing
#include "libfive/render/brep/hybrid/hybrid_worker_pool.hpp"
#include "libfive/render/brep/hybrid/hybrid_mesher.hpp"

namespace libfive {

std::unique_ptr<Mesh> Mesh::render(const Tree& t_, const Region<3>& r,
                                   const BRepSettings& settings)
{
    // 定义了一个使用 Eigen 对齐分配器的 std::vector，用于存储 Evaluator 类型的对象
    //  Eigen::aligned_allocator 是 Eigen 库提供的一个自定义内存分配器，用于处理对齐要求较高的内存分配
    std::vector<Evaluator, Eigen::aligned_allocator<Evaluator>> es;
    // 预分配足够的内存空间以容纳 settings.workers (8) 个 Evaluator 对象
    // reserve() 是 std::vector 的一个成员函数，用于提前分配内存空间以容纳指定数量的元素。这样做可以避免在向 std::vector 中添加元素时频繁地重新分配内存，提高性能
    es.reserve(settings.workers);
    //优化树。eg：优化后的树还会折叠嵌套的仿射形式，例如：(2*X + 3*Y) + 5*(X - Y) ==> 7*X - 2*Y
    const auto t = t_.optimized();
    for (unsigned i=0; i < settings.workers; ++i) {
        es.emplace_back(Evaluator(t));
    }

    //es.data() 在代码中用于获取 es 容器中存储的数据的指针
    // Evaluator对象 完成deck（tape）block value等数据设置、评估点是否inside、雅可比矩阵计算等方法
    return render(es.data(), r, settings);
}

//函数的返回值是一个指向 Mesh 对象的 std::unique_ptr，这意味着生成的网格对象的所有权被转移给调用者，保证资源的自动管理和释放
std::unique_ptr<Mesh> Mesh::render(
        Evaluator* es,
        const Region<3>& r, const BRepSettings& settings)
{
    //是 C++ 中使用智能指针 std::unique_ptr 定义一个指向 Mesh 对象的独占所有权指针
    // out是输出结果。继承父类属性：verts，branes
    std::unique_ptr<Mesh> out;
    if (settings.alg == DUAL_CONTOURING)
    {
        if (settings.progress_handler) {
            // Pool::build, Dual::walk, t.reset
            settings.progress_handler->start({1, 1, 1});
        }
        //这段代码实现了 WorkerPool::build 函数的功能，主要用于构建一个递归树结构（如四叉树、八叉树），
        // 通过多线程处理任务队列，以高效地处理评估任务并生成所需的树结构
        auto t = DCWorkerPool<3>::build(es, r, settings);

        if (settings.cancel.load() || t.get() == nullptr) {
            if (settings.progress_handler) {
                settings.progress_handler->finish();
            }
            return nullptr;
        }

        // Perform marching squares
        //Dual<3>::walk<DCMesher>(t, settings) 表示：
        //
        //调用 Dual<3> 类中的 walk 静态模板方法。
        //walk 函数被实例化为使用 DCMesher 类型的具体实现。
        //传入 t 和 settings 作为参数，进行具体的操作。
        //walk函数将针对DCMesher这个特定的类型进行特化。
        //Dual<3>::walk<DCMesher>(t, settings); 这段代码中的 <DCMesher> 是模板参数，表示 walk 函数在执行时使用了 DCMesher 类型作为模板参数 M。
        // 在 C++ 中，模板参数的作用是让函数或类在编译时根据具体的类型进行实例化，从而实现泛型编程
        out = Dual<3>::walk<DCMesher>(t, settings);

        // TODO: check for early return here again
        t.reset(settings);
    }
    else if (settings.alg == ISO_SIMPLEX)
    {
        if (settings.progress_handler) {
            // Pool::build, Dual::walk, t->assignIndices, t.reset
            settings.progress_handler->start({1, 1, 1});
        }
        auto t = SimplexWorkerPool<3>::build(es, r, settings);

        if (settings.cancel.load() || t.get() == nullptr) {
            if (settings.progress_handler) {
                settings.progress_handler->finish();
            }
            return nullptr;
        }

        t->assignIndices(settings);

        out = Dual<3>::walk_<SimplexMesher>(t, settings,
                [&](PerThreadBRep<3>& brep, int i) {
                    return SimplexMesher(brep, &es[i]);
                });
        t.reset(settings);
    }
    else if (settings.alg == HYBRID)
    {
        if (settings.progress_handler) {
            // Pool::build, Dual::walk, t->assignIndices, t.reset
            settings.progress_handler->start({1, 1, 1});
        }
        auto t = HybridWorkerPool<3>::build(es, r, settings);

        if (settings.cancel.load() || t.get() == nullptr) {
            if (settings.progress_handler) {
                settings.progress_handler->finish();
            }
            return nullptr;
        }

        t->assignIndices(settings);

        out = Dual<3>::walk_<HybridMesher>(t, settings,
                [&](PerThreadBRep<3>& brep, int i) {
                    return HybridMesher(brep, &es[i]);
                });
        t.reset(settings);
    }

    if (settings.progress_handler) {
        settings.progress_handler->finish();
    }
    return out;
}

void Mesh::line(const Eigen::Vector3f& a, const Eigen::Vector3f& b)
{
    uint32_t a_ = verts.size();
    verts.push_back(a);
    uint32_t b_ = verts.size();
    verts.push_back(b);

    branes.push_back({a_, a_, b_});
}

////////////////////////////////////////////////////////////////////////////////

bool Mesh::saveSTL(const std::string& filename,
                   const std::list<const Mesh*>& meshes)
{
    if (!boost::algorithm::iends_with(filename, ".stl"))
    {
        std::cerr << "Mesh::saveSTL: filename \"" << filename
                  << "\" does not end in .stl" << std::endl;
    }
    std::ofstream file;
    file.open(filename, std::ios::out | std::ios::binary);
    if (!file.is_open())
    {
        std::cout << "Mesh::saveSTL: could not open " << filename
                  << std::endl;
        return false;
    }

    // File header (giving human-readable info about file type)
    std::string header = "This is a binary STL exported from libfive.";
    file.write(header.c_str(), header.length());

    // Pad the rest of the header to 80 bytes
    for (int i=header.length(); i < 80; ++i)
    {
        file.put(' ');
    }

    // Write the triangle count to the file
    uint32_t num = std::accumulate(meshes.begin(), meshes.end(), (uint32_t)0,
            [](uint32_t i, const Mesh* m){ return i + m->branes.size(); });
    file.write(reinterpret_cast<char*>(&num), sizeof(num));

    for (const auto& m : meshes)
    {
        for (const auto& t : m->branes)
        {
            // Write out the normal vector for this face (all zeros)
            float norm[3] = {0, 0, 0};
            file.write(reinterpret_cast<char*>(&norm), sizeof(norm));

            // Iterate over vertices (which are indices into the verts list)
            for (unsigned i=0; i < 3; ++i)
            {
                auto v = m->verts[t[i]];
                float vert[3] = {v.x(), v.y(), v.z()};
                file.write(reinterpret_cast<char*>(&vert), sizeof(vert));
            }

            // Write out this face's attribute short
            uint16_t attrib = 0;
            file.write(reinterpret_cast<char*>(&attrib), sizeof(attrib));
        }
    }

    return true;
}

bool Mesh::saveSTL(const std::string& filename) const
{
    return saveSTL(filename, {this});
}

}   // namespace libfive
