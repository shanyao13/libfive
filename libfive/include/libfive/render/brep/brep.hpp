/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2017  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <future>
#include <vector>

#include <Eigen/Eigen>
#include <Eigen/StdVector>

#include "libfive/render/brep/per_thread_brep.hpp"

namespace libfive {

    //这段代码表明 BRep 是一个模板类，其中 N 是一个无符号整数类型的模板参数。如 N=2 时表示二维几何，N=3 时表示三维几何
template <unsigned N>
class BRep
{
public:
    BRep() { verts.push_back(Eigen::Matrix<float, N, 1>::Zero()); }

    /*  Flat array of point positions
     *  The 0th position is reserved as a marker
     *  顶点数组
     *  定义顶点数组，使用 N 维向量。
     *  Eigen::Matrix<float, N, 1>：表示 N 维的浮点型向量。例如，当 N=3 时，它是一个三维向量 (x, y, z)*/
    std::vector<Eigen::Matrix<float, N, 1>,
                Eigen::aligned_allocator<Eigen::Matrix<float, N, 1>>> verts;

    /*  [N-1]-dimensional objects (line segments, triangles)
     * 几何对象数组
     * 定义 branes 数组，存储 N-1 维的几何对象。
     * Eigen::Matrix<uint32_t, N, 1>：表示 N 维的无符号整数向量，通常用于存储顶点的索引*/
    std::vector<Eigen::Matrix<uint32_t, N, 1>,
                Eigen::aligned_allocator<Eigen::Matrix<uint32_t, N, 1>>> branes;

    //pushVertex 函数接收一个 N 维的向量，并将其添加到 verts 中
    uint32_t pushVertex(const Eigen::Matrix<float, N, 1>& v) {
        uint32_t out = verts.size();
        verts.push_back(v);
        return out;
    }

    /*
     *  Collects a set of PerThreadBRep objects into this BRep.
     *  The children must be a valid set of breps, i.e. generated with
     *  the same atomic index, so that their indices are globally
     *  unique and completely fill the range starting from 1.
     *
     *  If workers is 0, spins up one thread per child, otherwise spins
     *  up the requested number of threads to do the merge in parallel.
     *  collect 方法用于将多个子 BRep 对象（PerThreadBRep）合并到当前对象中。
使用多线程并行处理，以提高合并速度。workers 参数控制并行线程的数量，如果为 0，则默认每个子对象使用一个线程。
方法首先根据所有子对象的大小调整 verts 和 branes 数组的大小。
在多线程中，通过异步操作 std::async 并行处理每个子对象的顶点和几何对象，确保所有数据被正确地合并到当前对象中
     */
    void collect(const std::vector<PerThreadBRep<N>>& children,
                 unsigned workers=0)
    {
        assert(verts.size() == 1);
        assert(branes.size() == 0);

        if (workers == 0) {
            workers = children.size();
        }

        // Build big enough vectors to hold everything, since we're going
        // to be dropping items in through multiple threads.
        size_t num_verts = 1;
        size_t num_branes = 0;
        for (const auto& c : children) {
            num_verts += c.verts.size();
            num_branes += c.branes.size();
        }
        verts.resize(num_verts);
        branes.resize(num_branes);

        std::vector<std::future<void>> futures;
        futures.resize(workers);

        for (unsigned i=0; i < workers; ++i) {
            futures[i] = std::async(std::launch::async,
                [i, workers, this, &children]() {
                    for (unsigned j=i; j < children.size(); j += workers) {
                        const auto& c = children[j];

                        // Unpack vertices, which all have unique indexes into
                        // our collecting vertex array.
                        for (unsigned k=0; k < c.indices.size(); ++k) {
                            verts.at(c.indices.at(k)) = c.verts.at(k);
                        }

                        // Figure out where to position the branes in the
                        // collecting branes array, using a simple offset
                        // from the start.
                        size_t offset = 0;
                        for (unsigned k=0; k < j; ++k) {
                            offset += children[k].branes.size();
                        }

                        // Then save all of the branes
                        for (unsigned k=0; k < c.branes.size(); ++k) {
                            branes[offset + k] = c.branes[k];
                        }
                    }
                }
            );
        }

        for (auto& f : futures) {
            f.wait();
        }
    }
};

}   // namespace libfive
