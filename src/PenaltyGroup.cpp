#include "PenaltyGroup.h"
#include "PenaltyPotential.h"
#include <iostream>

using namespace Eigen;
using namespace std;

PenaltyGroup::PenaltyGroup(double dt, double outerEta, double innerEta, double stiffness) : nextstep_(0), dt_(dt), outerEta_(outerEta), innerEta_(innerEta), stiffness_(stiffness)
{	
}

PenaltyGroup::~PenaltyGroup()
{
}

void PenaltyGroup::addVFStencil(VertexFaceStencil vfstencil)
{
	vfstencils_.push_back(vfstencil);
}

void PenaltyGroup::addEEStencil(EdgeEdgeStencil eestencil)
{
	eestencils_.push_back(eestencil);
}

bool PenaltyGroup::addForce(const Eigen::VectorXd &q, Eigen::VectorXd &F)
{
	VectorXd groupforce(q.size());
	groupforce.setZero();
	bool newused = false;
	for(vector<VertexFaceStencil>::iterator it = vfstencils_.begin(); it != vfstencils_.end(); ++it)
	{
		bool fired = VertexFacePenaltyPotential::addForce(q,groupforce, *it, outerEta_, innerEta_, stiffness_);
		newused |= (fired && it->isnew);
	}
	for(vector<EdgeEdgeStencil>::iterator it = eestencils_.begin(); it != eestencils_.end(); ++it)		
	{
		bool fired = EdgeEdgePenaltyPotential::addForce(q,groupforce, *it, outerEta_, innerEta_, stiffness_);
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
	nextstep_ = 0;
	for(vector<VertexFaceStencil>::iterator it = vfstencils_.begin(); it != vfstencils_.end(); ++it)
		it->isnew = false;
	for(vector<EdgeEdgeStencil>::iterator it = eestencils_.begin(); it != eestencils_.end(); ++it)
		it->isnew = false;
}
