#include "PenaltyGroup.h"
#include "PenaltyPotential.h"

using namespace Eigen;
using namespace std;

PenaltyGroup::PenaltyGroup(double curt, double dt, double eta, double stiffness) : dt_(dt), eta_(eta), stiffness_(stiffness)
{
	double steps = curt/dt;
	nextstep_ = ceil(steps);
}

PenaltyGroup::~PenaltyGroup()
{
	for(vector<VertexFacePenaltyPotential *>::iterator it = vfforces_.begin(); it != vfforces_.end(); ++it)
		delete *it;
	for(vector<EdgeEdgePenaltyPotential *>::iterator it = eeforces_.begin(); it != eeforces_.end(); ++it)
		delete *it;
}

void PenaltyGroup::addVFStencil(VertexFaceStencil vfstencil)
{
	vfforces_.push_back(new VertexFacePenaltyPotential(vfstencil, eta_, stiffness_));
}

void PenaltyGroup::addEEStencil(EdgeEdgeStencil eestencil)
{
	eeforces_.push_back(new EdgeEdgePenaltyPotential(eestencil, eta_, stiffness_));
}

void PenaltyGroup::addForce(const Eigen::VectorXd &q, Eigen::VectorXd &F)
{
	for(vector<VertexFacePenaltyPotential *>::iterator it = vfforces_.begin(); it != vfforces_.end(); ++it)
		(*it)->addForce(q,F);
	for(vector<EdgeEdgePenaltyPotential *>::iterator it = eeforces_.begin(); it != eeforces_.end(); ++it)		
		(*it)->addForce(q,F);
}

void PenaltyGroup::incrementTimeStep()
{
	nextstep_++;
}

double PenaltyGroup::nextFireTime() const
{
	return nextstep_ * dt_;
}
