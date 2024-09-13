/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2017  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "util/shapes.hpp"
#include "../stdlib/libfive_stdlib.h"
#include "libfive/render/brep/mesh.hpp"
#include "libfive/render/brep/settings.hpp"

#include <iostream>

#include "/Users/shanyao/zjlWorkSpace/clworkSpace/libfive-scCode-read/libfive/libfive/stdlib/stdlib_impl.cpp"

int main()
{
    auto mySphere = sphere(0.2);
    auto myBox = box({-0.4, -0.4, -0.2}, {0.4, 0.4, 0});

    //1 union
    auto myUnion_sphereBox = min(mySphere, myBox);

    //2 intersection
    auto myIntersection = max(mySphere, myBox);


    //3 inverse
    auto myInverse = inverse(myBox);

    //4 inverse
    auto myDifferance = difference(myBox, mySphere);

    //5 offset
    auto myOffset = difference(myBox, offset(mySphere, 0.1));

    //6 clearance
    auto myClearance = clearance(myBox, mySphere, 0.1);

    //7 shell
    auto myShell = shell(myBox, 0.1);

    //8 blend
    float para = 0.5;
    auto myBlend = blend(myBox, mySphere, para);

    //9 blend_expt  m取值0.5到2，m>1更柔和的过度和光滑表面
    auto myBlend_expt_1 = blend_expt(myBox, mySphere, 1);

    //9 blend_expt_unit
    auto myBlend_expt_unit = blend_expt_unit(myBox, mySphere, 0.1);

    //10 blend_difference
    auto myBlend_difference = blend_difference(myBox, mySphere, 0.1, 0.2);

    // morph
    auto myMorph = morph(myBox, mySphere, 0.01);
    auto myMorph_2 = morph(myBox, mySphere, 0.9);

    int sopnge_level = 2; //     //表示递归的次数或分形的深度。它决定了几何图形被细分和分形的程度，从而影响了结构的复杂性和精细度。
    libfive::Tree sponge = menger(sopnge_level);

    libfive::Region<3> r({-0.5, -0.5, -0.5}, {1, 1, 1});

    libfive::BRepSettings settings ;
    settings.workers = 8;
    // settings.alg = libfive::ISO_SIMPLEX;
//         DUAL_CONTOURING,
// ISO_SIMPLEX,
// HYBRID,

    // auto mesh = libfive::Mesh::render(sponge, r, settings);

    // 开始计时
    auto start = std::chrono::high_resolution_clock::now();

    // 调用目标函数
    auto mesh = libfive::Mesh::render(sponge, r, settings);

    // 结束计时
    auto end = std::chrono::high_resolution_clock::now();

    // 计算耗时
    std::chrono::duration<double> elapsed = end - start;

    // 输出耗时
    std::cout << "Function execution time: " << elapsed.count() << " seconds" << std::endl;

    std::string para_ALgori = "iso-worker16";
    std::string Path = "./sponge_" + para_ALgori + ".stl";
//    mesh->saveSTL("/home/dell/clWorkSpace2/libfive/result/csg/myBlend.stl");
    mesh->saveSTL(Path);

    std::cout << "hi" << std::endl;

    return 0;
}

//int main(int argc, char** argv)
//{
//    return Catch::Session().run(argc, argv);
//}
