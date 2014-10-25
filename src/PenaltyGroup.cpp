#include "PenaltyGroup.h"
#include "PenaltyPotential.h"
#include <iostream>

using namespace Eigen;
using namespace std;

PenaltyGroup::PenaltyGroup(double dt, double outerEta, double innerEta, double stiffness, double CoR) : nextstep_(1), dt_(dt), outerEta_(outerEta), innerEta_(innerEta), stiffness_(stiffness), CoR_(CoR)
{	
}

PenaltyGroup::~PenaltyGroup()
{
}

void PenaltyGroup::addVFStencil(VertexFaceStencil vfstencil)
{
	vfstencils_.push_back(vfstencil);
	groupStencil_.insert(vfstencil.p);
	groupStencil_.insert(vfstencil.q0);
	groupStencil_.insert(vfstencil.q1);
	groupStencil_.insert(vfstencil.q2);
}

void PenaltyGroup::addEEStencil(EdgeEdgeStencil eestencil)
{
	eestencils_.push_back(eestencil);
	groupStencil_.insert(eestencil.p0);
	groupStencil_.insert(eestencil.p1);
	groupStencil_.insert(eestencil.q0);
	groupStencil_.insert(eestencil.q1);
}

bool PenaltyGroup::addForce(const Eigen::VectorXd &q, const Eigen::VectorXd &v, Eigen::VectorXd &F)
{
	VectorXd groupforce(q.size());
	groupforce.setZero();
	bool newused = false;
	for(vector<VertexFaceStencil>::iterator it = vfstencils_.begin(); it != vfstencils_.end(); ++it)
	{
		bool fired = VertexFacePenaltyPotential::addForce(q, v, groupforce, *it, outerEta_, innerEta_, stiffness_, CoR_);
		newused |= (fired && it->isnew);
	}
	for(vector<EdgeEdgeStencil>::iterator it = eestencils_.begin(); it != eestencils_.end(); ++it)		
	{
		bool fired = EdgeEdgePenaltyPotential::addForce(q, v, groupforce, *it, outerEta_, innerEta_, stiffness_, CoR_);
		newused |= (fired && it->isnew);
	}

	F += groupforce * dt_;
	return newused;
}

void PenaltyGroup::incrementTimeStep()
{
	nextstep_++;
}

double PenaltyGroup::nextFireTime() const
{
	return nextstep_ * dt_;
}

void PenaltyGroup::rollback()
{
	nextstep_ = 1;
	for(vector<VertexFaceStencil>::iterator it = vfstencils_.begin(); it != vfstencils_.end(); ++it)
		it->isnew = false;
	for(vector<EdgeEdgeStencil>::iterator it = eestencils_.begin(); it != eestencils_.end(); ++it)
		it->isnew = false;
}
