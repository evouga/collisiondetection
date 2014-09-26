#ifndef PENALTYPOTENTIAL_H
#define PENALTYPOTENTIAL_H

#include <Eigen/Core>

struct VertexFaceStencil
{
  int p;
  int q0,q1,q2;
};

struct EdgeEdgeStencil
{
  int p0,p1;
  int q0,q1;
};



// V(q) = 0,                             primitiveDist >= eta
//        1/2 k (primitiveDist - eta)^2  primitiveDist < eta

class PenaltyPotential
{
 public:
  PenaltyPotential(double eta, double stiffness) : eta_(eta), stiffness_(stiffness) {}
  virtual ~PenaltyPotential() {}

  virtual void addForce(const Eigen::VectorXd &q, Eigen::VectorXd &F)=0;

 protected:
  double eta_;
  double stiffness_;
};

class VertexFacePenaltyPotential : public PenaltyPotential
{
 public:
  VertexFacePenaltyPotential(VertexFaceStencil stencil, double eta, double stiffness) : PenaltyPotential(eta, stiffness), stencil_(stencil) {}

  virtual void addForce(const Eigen::VectorXd &q, Eigen::VectorXd &F);

 private:
  VertexFaceStencil stencil_;
};

class EdgeEdgePenaltyPotential : public PenaltyPotential
{
 public:
  EdgeEdgePenaltyPotential(EdgeEdgeStencil stencil, double eta, double stiffness) : PenaltyPotential(eta, stiffness), stencil_(stencil) {}

  virtual void addForce(const Eigen::VectorXd &q, Eigen::VectorXd &F);

 private:
  EdgeEdgeStencil stencil_;
};


#endif
