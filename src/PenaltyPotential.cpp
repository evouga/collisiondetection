#include "PenaltyPotential.h"
#include "Distance.h"
#include <iostream>

using namespace Eigen;	

bool VertexFacePenaltyPotential::addForce(const Eigen::VectorXd &q, const Eigen::VectorXd &v, Eigen::VectorXd &F, VertexFaceStencil stencil, double outerEta, double innerEta, double stiffness, double CoR)
{
  double bary0, bary1, bary2;
  const Vector3d &ppos = q.segment<3>(3*stencil.p);
  const Vector3d &q0pos = q.segment<3>(3*stencil.q0);
  const Vector3d &q1pos = q.segment<3>(3*stencil.q1);
  const Vector3d &q2pos = q.segment<3>(3*stencil.q2);
  if(!Distance::vertexPlaneDistanceLessThan(ppos, q0pos, q1pos, q2pos, outerEta))
	return false;

  Vector3d closestVec = Distance::vertexFaceDistance(ppos, q0pos, q1pos, q2pos,
						     bary0, bary1, bary2);
  double dist = closestVec.norm();
  if(dist >= outerEta || dist < innerEta)
    return false;

  Vector3d relvel = -v.segment<3>(3*stencil.p) + bary0*v.segment<3>(3*stencil.q0) + bary1*v.segment<3>(3*stencil.q1) + bary2*v.segment<3>(3*stencil.q2);
  if(relvel.dot(closestVec) > 0)
	stiffness *= CoR;

  Vector3d localF = stiffness * (outerEta - dist )/(outerEta - innerEta) * closestVec/dist;
  //std::cout << "VF " << stencil.p << " " << stencil.q0 << " " <<  stencil.q1 << " " << stencil.q2 << " " << stiffness << " " << dist << " " << localF.transpose() << std::endl;
  F.segment<3>(3*stencil.p) -= localF;
  F.segment<3>(3*stencil.q0) += bary0*localF;
  F.segment<3>(3*stencil.q1) += bary1*localF;
  F.segment<3>(3*stencil.q2) += bary2*localF;
  return true;
}

bool EdgeEdgePenaltyPotential::addForce(const Eigen::VectorXd &q, const Eigen::VectorXd &v, Eigen::VectorXd &F, EdgeEdgeStencil stencil, double outerEta, double innerEta, double stiffness, double CoR)
{
  double baryp0, baryp1, baryq0, baryq1;
  const Vector3d &p0pos = q.segment<3>(3*stencil.p0);
  const Vector3d &p1pos = q.segment<3>(3*stencil.p1);
  const Vector3d &q0pos = q.segment<3>(3*stencil.q0);
  const Vector3d &q1pos = q.segment<3>(3*stencil.q1);
  if(!Distance::lineLineDistanceLessThan(p0pos, p1pos, q0pos, q1pos, outerEta))
	return false;

  Vector3d closestVec = Distance::edgeEdgeDistance(p0pos, p1pos, 
						   q0pos, q1pos,
						   baryp0, baryp1, baryq0, baryq1);
  double dist = closestVec.norm();
  if(dist >= outerEta || dist < innerEta)
    return false;

  Vector3d relvel = -baryp0*v.segment<3>(3*stencil.p0) - baryp1*v.segment<3>(3*stencil.p1) + baryq0*v.segment<3>(3*stencil.q0) + baryq1*v.segment<3>(3*stencil.q1);
  if(relvel.dot(closestVec) > 0)
	stiffness *= CoR;

  Vector3d localF = stiffness * (outerEta - dist)/(outerEta - innerEta) * closestVec/dist;
  //std::cout << "EE " << stencil.p0 << " " << stencil.p1 << " " <<  stencil.q0 << " " << stencil.q1 << " " << stiffness << " " << dist << " " << localF.transpose() << " " << baryp0 << " " << baryp1 << " " << baryq0 << " " << baryq1 << std::endl;
  F.segment<3>(3*stencil.p0) -= baryp0*localF;
  F.segment<3>(3*stencil.p1) -= baryp1*localF;
  F.segment<3>(3*stencil.q0) += baryq0*localF;
  F.segment<3>(3*stencil.q1) += baryq1*localF;
  return true;
}
