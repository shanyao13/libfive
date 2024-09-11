/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2017  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include "libfive/tree/tree.hpp"

#include "libfive/render/brep/brep.hpp"

namespace libfive {

// Forward declaration：在 C++ 中，前向声明（Forward Declaration）用于告知编译器某个类或结构体的存在，而不提供其完整定义。这种声明方式在编译时可以减少依赖关系，提高编译速度，并解决循环依赖的问题。
class Evaluator;
struct BRepSettings;

template <unsigned N> class Region;

//在这个 Mesh 类中，有多个静态成员函数和一个保护成员函数，主要用于网格渲染、文件保存和调试功能
class Mesh : public BRep<3> {
public:
    /*
     *  Core render function
     *
     *  Returns nullptr if min_feature is invalid or cancel is set to true
     *  partway through the computation.
     */
    static std::unique_ptr<Mesh> render(
            const Tree& t, const Region<3>& r,
            const BRepSettings& settings);

    /*
     *  Render function that re-uses evaluators
     *  es must be a pointer to at least [settings.workers] Evaluators
     *
     *  Returns nullptr if min_feature is invalid or cancel is set to true
     *  partway through the computation.
     */
    static std::unique_ptr<Mesh> render(
            Evaluator* es, const Region<3>& r,
            const BRepSettings& settings);

    /*
     *  Writes the mesh to a file
     */
    bool saveSTL(const std::string& filename) const;

    /*
     *  Merge multiple bodies and write them to a single file
     */
    static bool saveSTL(const std::string& filename,
                        const std::list<const Mesh*>& meshes);

protected:

    /*
     *  Inserts a line into the mesh as a zero-size triangle
     *  (used for debugging)
     */
    void line(const Eigen::Vector3f& a, const Eigen::Vector3f& b);

};

}   // namespace libfive
