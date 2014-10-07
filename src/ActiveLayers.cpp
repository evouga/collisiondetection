#include "ActiveLayers.h"
#include "PenaltyGroup.h"
#include "SimulationState.h"

using namespace Eigen;

ActiveLayers::ActiveLayers(double eta, double baseDt, double baseStiffness, double terminationTime) : eta_(eta), baseDt_(baseDt), baseStiffness_(baseStiffness), termTime_(terminationTime), deepestLayer_(0)
{
}

ActiveLayers::~ActiveLayers()
{
	for(std::vector<PenaltyGroup *>::iterator it = groups_.begin(); it != groups_.end(); ++it)
		delete *it;
}

void ActiveLayers::addVFStencil(double curt, VertexFaceStencil stencil)
{
	int olddepth = vfdepth_[stencil];
	int newdepth = 1 + olddepth;

	addGroups(curt, newdepth);

	groups_[olddepth]->addVFStencil(stencil);
	vfdepth_[stencil]++;	
}

void ActiveLayers::addEEStencil(double curt, EdgeEdgeStencil stencil)
{
	int olddepth = eedepth_[stencil];
	int newdepth = 1 + olddepth;

	addGroups(curt, newdepth);

	groups_[olddepth]->addEEStencil(stencil);
	eedepth_[stencil]++;	
}

void ActiveLayers::addGroups(double curt, int maxdepth)
{
	while(deepestLayer_ < maxdepth)
	{
		int depth = deepestLayer_+1;
		double ki = baseStiffness_*depth;
		double etai = eta_ / sqrt(double(depth));
		double dti = baseDt_ / sqrt(double(depth));

		PenaltyGroup *newgroup = new PenaltyGroup(curt, dti, etai, ki);
		groups_.push_back(newgroup);
		groupQueue_.push(newgroup);

		deepestLayer_++;
	}
}

bool ActiveLayers::step(SimulationState &s)
{
	if(!groupQueue_.empty())
	{
		PenaltyGroup *group = groupQueue_.top();
		groupQueue_.pop();

		double newtime = group->nextFireTime();
		if(newtime < termTime_)
		{
			double dt = newtime - s.time;
			s.q += dt * s.v;
			VectorXd F(s.q.size());
			group->addForce(s.q, F);
			for(int i=0; i<F.size(); i++)
				F[i] /= s.m[i];
			s.v += F;
			group->incrementTimeStep();
			groupQueue_.push(group);
			return false;
		}
	}

	double dt = termTime_ - s.time;
	s.q += dt * s.v;
	return true;
}
