#include "ActiveLayers.h"
#include "PenaltyGroup.h"
#include "SimulationState.h"
#include "CTCD.h"
#include "Distance.h"
#include "Mesh.h"
#include <set>
#include <iostream>
#include "History.h"
#include "CTCDNarrowPhase.h"
#include "KDOPBroadPhase.h"
#include "TrivialBroadPhase.h"
#include "SeparatingPlaneNarrowPhase.h"

using namespace Eigen;
using namespace std;

bool PenaltyGroupComparator::operator()(const PenaltyGroup *first, const PenaltyGroup *second) const
{
	return second->nextFireTime() < first->nextFireTime();
}

ActiveLayers::ActiveLayers(double outerEta, double innerEta, double baseDt, double baseStiffness, double terminationTime, double CoR, bool verbose) 
	: outerEta_(outerEta), innerEta_(innerEta), baseDt_(baseDt), baseStiffness_(baseStiffness), termTime_(terminationTime), CoR_(CoR), deepestLayer_(0), history_(NULL), verbose_(verbose)
{
	bp_ = new KDOPBroadPhase();
	//bp_ = new TrivialBroadPhase();
	np_ = new SeparatingPlaneNarrowPhase();
}

ActiveLayers::~ActiveLayers()
{
	for(std::vector<PenaltyGroup *>::iterator it = groups_.begin(); it != groups_.end(); ++it)
		delete *it;
	
	delete history_;
	delete bp_;
	delete np_;
}

void ActiveLayers::addVFStencil(VertexFaceStencil stencil)
{
	int olddepth = vfdepth_[stencil];

	for(int i=0; i<5; i++)
	{
		addGroups(olddepth+i+1);

		groups_[olddepth+i]->addVFStencil(stencil);
		vfdepth_[stencil]++;	
	}
}

void ActiveLayers::addEEStencil(EdgeEdgeStencil stencil)
{
	int olddepth = eedepth_[stencil];

	for(int i=0; i<5; i++)
	{
		addGroups(olddepth+i+1);

		groups_[olddepth+i]->addEEStencil(stencil);
		eedepth_[stencil]++;	
	}
}

void ActiveLayers::addGroups(int maxdepth)
{
	while(deepestLayer_ < maxdepth)
	{
		int depth = deepestLayer_+1;
		double fudge = 1e-4; // To stop layer timesteps from exactly coinciding
		double ki = baseStiffness_*depth*depth*depth;
		double etai = layerDepth(depth);
		double dti = baseDt_ / double(depth) / sqrt(double(depth)+fudge);

		PenaltyGroup *newgroup = new PenaltyGroup(dti, etai, innerEta_, ki, CoR_);
		groups_.push_back(newgroup);
		groupQueue_.push(newgroup);

		deepestLayer_++;
	}
}

bool ActiveLayers::step(SimulationState &s)
{
	int nverts = s.q.size()/3;
	VectorXd F(s.q.size());
	VectorXd newq(s.q.size());
	VectorXd newv(s.q.size());

	if(!groupQueue_.empty())
	{
		PenaltyGroup *group = groupQueue_.top();

		double newtime = group->nextFireTime();
		if(newtime < termTime_)
		{
			groupQueue_.pop();
			for(set<int>::iterator it = group->getGroupStencil().begin(); it != group->getGroupStencil().end(); ++it)
			{
				int vert = *it;
				for(int j=0; j<3; j++)
				{
					double dt = newtime - s.lastUpdateTime[3*vert+j];
					newq[3*vert+j] = s.q[3*vert+j] + (newtime-s.lastUpdateTime[3*vert+j])*s.v[3*vert+j];
					newv[3*vert+j] = (newq[3*vert+j]-s.q[3*vert+j])/dt;
				}
			}

			F.setZero();
			group->addForce(newq, newv, F);

			for(set<int>::iterator it = group->getGroupStencil().begin(); it != group->getGroupStencil().end(); ++it)
			{
				int vert = *it;
				bool touched = false;
				for(int j=0; j<3; j++)
				{
					if(F[3*vert+j] != 0.0)
						touched = true;
				}
				if(touched)
				{
					for(int j=0; j<3; j++)
					{
						s.v[3*vert+j] += s.minv[3*vert+j]*F[3*vert+j];
						s.q[3*vert+j] = newq[3*vert+j];
						s.lastUpdateTime[3*vert+j] = newtime;
					}
			
					history_->addHistory(vert, newtime, s.q.segment<3>(3*vert));
				}

			}
			group->incrementTimeStep();
			groupQueue_.push(group);
			return false;
		}
	}

	for(int i=0; i<3*nverts; i++)
	{
		s.q[i] += (termTime_ - s.lastUpdateTime[i])*s.v[i];
		s.lastUpdateTime[i] = termTime_;
	}
	history_->finishHistory(s.q);

	return true;
}

void ActiveLayers::rollback()
{
	for(std::vector<PenaltyGroup *>::iterator it = groups_.begin(); it != groups_.end(); ++it)
		(*it)->rollback();

	std::priority_queue<PenaltyGroup *, std::vector<PenaltyGroup *>, PenaltyGroupComparator> newqueue;
	while(!groupQueue_.empty())
	{
		newqueue.push(groupQueue_.top());
		groupQueue_.pop();
	}	

	groupQueue_ = newqueue;
}

double ActiveLayers::layerDepth(int layer)
{
	return innerEta_ + (outerEta_-innerEta_)/double(layer);
}

double ActiveLayers::VFStencilThickness(VertexFaceStencil stencil)
{
	std::map<VertexFaceStencil, int>::iterator it = vfdepth_.find(stencil);
	if(it != vfdepth_.end())
		return layerDepth(it->second + 1);
	return layerDepth(1);
}

double ActiveLayers::EEStencilThickness(EdgeEdgeStencil stencil)
{
	std::map<EdgeEdgeStencil, int>::iterator it = eedepth_.find(stencil);
	if(it != eedepth_.end())
		return layerDepth(it->second + 1);
	return layerDepth(1);
}

bool ActiveLayers::collisionDetection(const Mesh &m, set<VertexFaceStencil> &vfsToAdd, set<EdgeEdgeStencil> &eesToAdd, const set<int> &fixedVerts)
{
	vfsToAdd.clear();
	eesToAdd.clear();
	
	bp_->findCollisionCandidates(*history_, m, outerEta_, vfsToAdd, eesToAdd, fixedVerts);
	if(verbose_)
	  std::cout << "Broad phase found " << vfsToAdd.size() << " vertex-face and " << eesToAdd.size() << " edge-edge candidates" << std::endl;
	set<pair<VertexFaceStencil, double> > etavfs;
	set<pair<EdgeEdgeStencil, double> > etaees;
	for(set<VertexFaceStencil>::iterator it = vfsToAdd.begin(); it != vfsToAdd.end(); ++it)
	{
		pair<VertexFaceStencil, double> vp(*it, VFStencilThickness(*it));
		etavfs.insert(vp);
	}
	for(set<EdgeEdgeStencil>::iterator it = eesToAdd.begin(); it != eesToAdd.end(); ++it)
	{
		pair<EdgeEdgeStencil, double> vp(*it, EEStencilThickness(*it));
		etaees.insert(vp);
	}

	eesToAdd.clear();
	vfsToAdd.clear();

	np_->findCollisions(*history_, etavfs, etaees, vfsToAdd, eesToAdd);
	
	return(!vfsToAdd.empty() || !eesToAdd.empty());
}

bool ActiveLayers::runOneIteration(const Mesh &m, SimulationState &s)
{
	if(verbose_)
	{
		std::cout << "Taking an outer iteration, deepest layer is currently " << deepestLayer_;
	 	if(deepestLayer_)
			std::cout << " with outer thickness " << groups_[deepestLayer_-1]->getOuterEta() << " and dt " << groups_[deepestLayer_-1]->getDt();
		std::cout << std::endl;
	}
	delete history_;
	history_ = new History(s.q);

	while(!step(s));

	if(verbose_)
		std::cout << "Done simulating, accumulated " << history_->countHistoryEntries() << " history entries" << std::endl;

	set<VertexFaceStencil> vfsToAdd;
	set<EdgeEdgeStencil> eesToAdd;

	set<int> fixedVerts;
	for(int i=0; i<(int)s.minv.size(); i++)
	{
		if(s.minv[i] == 0.0)
			fixedVerts.insert(i/3);
	}

	bool collisions = collisionDetection(m, vfsToAdd, eesToAdd, fixedVerts);
	if(verbose_)
		std::cout << "Found " << vfsToAdd.size() << " vertex-face and " << eesToAdd.size() << " edge-edge collisions" << std::endl;

	rollback();

	for(std::set<VertexFaceStencil>::iterator it = vfsToAdd.begin(); it != vfsToAdd.end(); ++it)
		addVFStencil(*it);
	for(std::set<EdgeEdgeStencil>::iterator it = eesToAdd.begin(); it != eesToAdd.end(); ++it)
		addEEStencil(*it);

	return !collisions;
}
