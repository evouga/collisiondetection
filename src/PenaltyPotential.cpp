#include "PenaltyPotential.h"
#include "Distance.h"
#include <iostream>

using namespace Eigen;

bool VertexFacePenaltyPotential::addForce(const Eigen::VectorXd &q, Eigen::VectorXd &F, VertexFaceStencil stencil, double outerEta, double innerEta, double stiffness)
{
  double bary0, bary1, bary2;
  Vector3d closestVec = Distance::vertexFaceDistance(q.segment<3>(3*stencil.p), q.segment<3>(3*stencil.q0), q.segment<3>(3*stencil.q1), q.segment<3>(3*stencil.q2),
						     bary0, bary1, bary2);
  double dist = closestVec.norm();
  if(dist >= outerEta || dist < innerEta)
    return false;

  Vector3d localF = stiffness * (outerEta - dist )/(outerEta - innerEta) * closestVec/dist;
  //std::cout << "VF " << stencil.p << " " << stencil.q0 << " " <<  stencil.q1 << " " << stencil.q2 << std::endl;
  F.segment<3>(3*stencil.p) -= localF;
  F.segment<3>(3*stencil.q0) += bary0*localF;
  F.segment<3>(3*stencil.q1) += bary1*localF;
  F.segment<3>(3*stencil.q2) += bary2*localF;
  return true;
}

bool EdgeEdgePenaltyPotential::addForce(const Eigen::VectorXd &q, Eigen::VectorXd &F, EdgeEdgeStencil stencil, double outerEta, double innerEta, double stiffness)
{
  double baryp0, baryp1, baryq0, baryq1;
  Vector3d closestVec = Distance::edgeEdgeDistance(q.segment<3>(3*stencil.p0), q.segment<3>(3*stencil.p1), 
						   q.segment<3>(3*stencil.q0), q.segment<3>(3*stencil.q1),
						   baryp0, baryp1, baryq0, baryq1);
  double dist = closestVec.norm();
  if(dist >= outerEta || dist < innerEta)
    return false;

  Vector3d localF = stiffness * (outerEta - dist)/(outerEta - innerEta) * closestVec/dist;
  //std::cout << "EE " << stencil.p0 << " " << stencil.p1 << " " <<  stencil.q0 << " " << stencil.q1 << std::endl;
  F.segment<3>(3*stencil.p0) -= baryp0*localF;
  F.segment<3>(3*stencil.p1) -= baryp1*localF;
  F.segment<3>(3*stencil.q0) += baryq0*localF;
  F.segment<3>(3*stencil.q1) += baryq1*localF;
  return true;
}
