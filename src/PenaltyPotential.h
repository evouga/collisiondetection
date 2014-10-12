#ifndef PENALTYPOTENTIAL_H
#define PENALTYPOTENTIAL_H

#include <Eigen/Core>
#include "Stencils.h"

// V(q) = 0,                             primitiveDist >= eta
//        1/2 k (primitiveDist - eta)^2  primitiveDist < eta

class PenaltyPotential
{
 public:
  PenaltyPotential(double outerEta, double innerEta, double stiffness) : outerEta_(outerEta), innerEta_(innerEta), stiffness_(stiffness) {}
  virtual ~PenaltyPotential() {}

  virtual void addForce(const Eigen::VectorXd &q, Eigen::VectorXd &F)=0;

 protected:
  double outerEta_;
  double innerEta_;
  double stiffness_;
  double CoR_;
};

class VertexFacePenaltyPotential : public PenaltyPotential
{
 public:
  VertexFacePenaltyPotential(VertexFaceStencil stencil, double outerEta, double innerEta, double stiffness) : PenaltyPotential(outerEta, innerEta, stiffness), stencil_(stencil) {}

  virtual void addForce(const Eigen::VectorXd &q, Eigen::VectorXd &F);

 private:
  VertexFaceStencil stencil_;
};

class EdgeEdgePenaltyPotential : public PenaltyPotential
{
 public:
  EdgeEdgePenaltyPotential(EdgeEdgeStencil stencil, double outerEta, double innerEta, double stiffness) : PenaltyPotential(outerEta, innerEta, stiffness), stencil_(stencil) {}

  virtual void addForce(const Eigen::VectorXd &q, Eigen::VectorXd &F);

 private:
  EdgeEdgeStencil stencil_;
};


#endif
