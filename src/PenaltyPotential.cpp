#include "PenaltyPotential.h"
#include "Distance.h"
#include <iostream>

using namespace Eigen;

void VertexFacePenaltyPotential::addForce(const Eigen::VectorXd &q, Eigen::VectorXd &F)
{
  double bary0, bary1, bary2;
  Vector3d closestVec = Distance::vertexFaceDistance(q.segment<3>(3*stencil_.p), q.segment<3>(3*stencil_.q0), q.segment<3>(3*stencil_.q1), q.segment<3>(3*stencil_.q2),
						     bary0, bary1, bary2);
  double dist = closestVec.norm();
  if(dist >= eta_ || dist < 1e-12)
    return;

  Vector3d localF = stiffness_ * (eta_ - dist) * closestVec/dist;
  //std::cout << q.segment<3>(3*stencil_.p).transpose() << " " << q.segment<3>(3*stencil_.q0).transpose() << " " <<  q.segment<3>(3*stencil_.q1).transpose() << " " << q.segment<3>(3*stencil_.q2).transpose() << std::endl;
  F.segment<3>(3*stencil_.p) -= localF;
  F.segment<3>(3*stencil_.q0) += bary0*localF;
  F.segment<3>(3*stencil_.q1) += bary1*localF;
  F.segment<3>(3*stencil_.q2) += bary2*localF;
}

void EdgeEdgePenaltyPotential::addForce(const Eigen::VectorXd &q, Eigen::VectorXd &F)
{
  double baryp0, baryp1, baryq0, baryq1;
  Vector3d closestVec = Distance::edgeEdgeDistance(q.segment<3>(3*stencil_.p0), q.segment<3>(3*stencil_.p1), 
						   q.segment<3>(3*stencil_.q0), q.segment<3>(3*stencil_.q1),
						   baryp0, baryp1, baryq0, baryq1);
  double dist = closestVec.norm();
  if(dist >= eta_ || dist < 1e-12)
    return;

  Vector3d localF = stiffness_ * (eta_ - dist) * closestVec/dist;
  F.segment<3>(3*stencil_.p0) -= baryp0*localF;
  F.segment<3>(3*stencil_.p1) -= baryp1*localF;
  F.segment<3>(3*stencil_.q0) += baryq0*localF;
  F.segment<3>(3*stencil_.q1) += baryq1*localF;
}
