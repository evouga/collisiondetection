#ifndef PENALTYPOTENTIAL_H
#define PENALTYPOTENTIAL_H

#include <Eigen/Core>
#include "Stencils.h"

// V(q) = 0,                             primitiveDist >= eta
//        1/2 k (primitiveDist - eta)^2  primitiveDist < eta

class VertexFacePenaltyPotential
{
 public:
  static bool addForce(const Eigen::VectorXd &q, const Eigen::VectorXd &v, Eigen::VectorXd &F, VertexFaceStencil stencil, double outerEta, double innerEta, double stiffness, double CoR);
};

class EdgeEdgePenaltyPotential
{
 public:
  static bool addForce(const Eigen::VectorXd &q, const Eigen::VectorXd &v, Eigen::VectorXd &F, EdgeEdgeStencil stencil, double outerEta, double innerEta, double stiffness, double CoR);
};


#endif
