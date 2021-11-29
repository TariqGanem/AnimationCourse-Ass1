// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2016 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "shortest_edge_and_midpoint.h"
#include <vector>
#include <Eigen/LU>

IGL_INLINE void igl::shortest_edge_and_midpoint(
  const int e,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & /*F*/,
  const Eigen::MatrixXi & E,
  const Eigen::VectorXi & /*EMAP*/,
  const Eigen::MatrixXi & /*EF*/,
  const Eigen::MatrixXi & /*EI*/,
  double & cost,
  Eigen::RowVectorXd & p)
{
  cost = (V.row(E(e,0))-V.row(E(e,1))).norm();
  p = 0.5*(V.row(E(e,0))+V.row(E(e,1)));
}


IGL_INLINE void igl::quadric_error(
    const int e,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi&,
    const Eigen::MatrixXi& E,
    const Eigen::VectorXi&,
    const Eigen::MatrixXi&,
    const Eigen::MatrixXi&,
    double& cost,
    Eigen::RowVectorXd& p,
    std::vector<Eigen::Matrix4d> Q4)
{
    using namespace std;
    // Sum of the two Quadrics and mentioned in chapter 5
    Eigen::Matrix4d quadricRes = (Q4.at(E(e, 0)) + Q4.at(E(e, 1)));
    quadricRes.row(3) = Eigen::RowVector4d(0, 0, 0, 1);
    // if the matrix is invertible then inverse it and then multiply by [0,0,0,1] as the algo showed
    if (quadricRes.determinant() != 0 && !isnan(quadricRes.determinant())) {
        p = quadricRes.inverse() * Eigen::Vector4d(0, 0, 0, 1);
    }
    // if the matrix is not invertible we choose the midpoint
    else {
        p = 0.5 * (V.row(E(e, 0)) + V.row(E(e, 1)));
        Eigen::Vector4d temp = { p[0], p[1], p[2], 1 };
        p = temp;
    }
    // this is the error becomes the cost of the edge
    cost = (p * (Q4.at(E(e, 0)) + Q4.at(E(e, 1))) * p.transpose());
}

