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

#include "/home/dell/clWorkSpace2/libfive/libfive/stdlib/stdlib_impl.cpp"

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



    libfive::Region<3> r({-1, -1, -1}, {1, 1, 1});

    libfive::BRepSettings settings;
    auto mesh = libfive::Mesh::render(myBlend, r, settings);

    std::string Path = "/home/dell/clWorkSpace2/libfive/result/csg/myBlend_" + std::to_string(para) + ".stl";
//    mesh->saveSTL("/home/dell/clWorkSpace2/libfive/result/csg/myBlend.stl");
    mesh->saveSTL(Path);

    std::cout << "hi" << std::endl;

    return 0;
}

//int main(int argc, char** argv)
//{
//    return Catch::Session().run(argc, argv);
//}
