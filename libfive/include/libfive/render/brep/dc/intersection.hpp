/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <Eigen/Eigen>

#include "libfive/render/brep/dc/dc_flags.hpp"
#include "libfive/render/brep/default_new_delete.hpp"

namespace libfive {

    //这段代码定义了一个名为 Intersection 的模板结构体，属于 libfive 命名空间，用于处理和存储在离散曲面网格（DCMesh）算法中的交点信息。它涉及到求解最小二乘问题以确定交点的几何特性。
template <unsigned N>
struct Intersection {
    Intersection()
    {
        reset();
    }

    void reset() {
        AtA.array() = 0.0;
        AtB.array() = 0.0;
        BtB = 0.0;
        mass_point.array() = 0.0;
    }

    void push(Eigen::Matrix<double, N, 1> pos,
              Eigen::Matrix<double, N, 1> deriv,
              double value)
    {
        Eigen::Matrix<double, N + 1, 1> mp;
        mp << pos, 1.0;
        mass_point += mp;

#if !LIBFIVE_UNNORMALIZED_DERIVS
        // Find normalized derivatives and distance value
        const double norm = deriv.matrix().norm();
        deriv /= norm;
        value /= norm;

        // If the point has an invalid normal, then skip it
        if (norm <= 1e-12 || !deriv.array().isFinite().all()) {
            return;
        }
#endif
        AtA += deriv * deriv.transpose();
        const double b = deriv.dot(pos) - value;

        AtB += deriv * b;
        BtB += b * b;

        // Reset the recorded rank
        rank = -1;
    }

    Eigen::Matrix<double, N + 1, 1> normalized_mass_point() const {
        Eigen::Matrix<double, N + 1, 1> mp;
        mp << (mass_point.template head<N>() / mass_point(N)), 1;
        return mp;
    }

    //计算和返回交点的伪秩。伪秩用于描述 AtA 矩阵的非满秩情况，通常用于最小二乘问题的条件数检查
    uint8_t get_rank() const {
        if (rank == -1) {
            // Use the pseudo-rank of the AtA matrix to assign
            // a rank to this particular intersection
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, N, N>> es(AtA);
            auto eigenvalues = es.eigenvalues().real();

#if LIBFIVE_UNNORMALIZED_DERIVS
            auto highest_val = eigenvalues.template lpNorm<Eigen::Infinity>();
            const double cutoff = (highest_val > 1e-20)
                ? highest_val * EIGENVALUE_CUTOFF
                : 0;
#else
            const double cutoff = EIGENVALUE_CUTOFF;
#endif
            rank = 0;
            for (unsigned j=0; j < N; ++j) {
                rank += (fabs(eigenvalues[j]) < cutoff);
            }
        }
        return rank;
    }

    // Represents the accumulated position of the intersection
    // (the final item is total number of samples)
    Eigen::Matrix<double, N + 1, 1> mass_point;

    Eigen::Matrix<double, N, N> AtA;
    Eigen::Matrix<double, N, 1> AtB;
    double BtB;
    mutable int8_t rank=-1;

    ALIGNED_OPERATOR_NEW_AND_DELETE(Intersection)
};

}   // namespace libfive
