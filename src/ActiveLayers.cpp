#include "ActiveLayers.h"
#include "PenaltyGroup.h"
#include "SimulationState.h"
#include "CTCD.h"
#include "Distance.h"
#include "Mesh.h"
#include <set>
#include <iostream>
#include "History.h"
#include "RetrospectiveDetection.h"

using namespace Eigen;
using namespace std;

bool PenaltyGroupComparator::operator()(const PenaltyGroup *first, const PenaltyGroup *second) const
{
	return second->nextFireTime() < first->nextFireTime();
}

ActiveLayers::ActiveLayers(double outerEta, double innerEta, double baseDt, double baseStiffness, double terminationTime, bool verbose) 
	: outerEta_(outerEta), innerEta_(innerEta), baseDt_(baseDt), baseStiffness_(baseStiffness), termTime_(terminationTime), deepestLayer_(0), history_(NULL), oldhistory_(NULL), verbose_(verbose), earliestTime_(0)
{
	bp_ = new TrivialBroadPhase();
	np_ = new CTCDNarrowPhase();
}

ActiveLayers::~ActiveLayers()
{
	for(std::vector<PenaltyGroup *>::iterator it = groups_.begin(); it != groups_.end(); ++it)
		delete *it;
	
	delete oldhistory_;
	delete history_;
	delete bp_;
	delete np_;
}

void ActiveLayers::addVFStencil(VertexFaceStencil stencil)
{
	int olddepth = vfdepth_[stencil];
	int newdepth = 1 + olddepth;

	addGroups(newdepth);

	groups_[olddepth]->addVFStencil(stencil);
	vfdepth_[stencil]++;	
}

void ActiveLayers::addEEStencil(EdgeEdgeStencil stencil)
{
	int olddepth = eedepth_[stencil];
	int newdepth = 1 + olddepth;

	addGroups(newdepth);

	groups_[olddepth]->addEEStencil(stencil);
	eedepth_[stencil]++;	
}

void ActiveLayers::addGroups(int maxdepth)
{
	while(deepestLayer_ < maxdepth)
	{
		int depth = deepestLayer_+1;
		double ki = baseStiffness_*depth*depth*depth;
		double etai = layerDepth(depth);
		double dti = baseDt_ / double(depth) / sqrt(double(depth));

		PenaltyGroup *newgroup = new PenaltyGroup(dti, etai, innerEta_, ki);
		groups_.push_back(newgroup);
		groupQueue_.push(newgroup);

		deepestLayer_++;
	}
}

bool ActiveLayers::step(SimulationState &s)
{
	int nverts = s.q.size()/3;

	if(!groupQueue_.empty())
	{
		PenaltyGroup *group = groupQueue_.top();

		double newtime = group->nextFireTime();
		if(newtime < termTime_)
		{
			std::cout << "now at time " << newtime << std::endl;
			groupQueue_.pop();
			VectorXd newq(3*nverts);
			for(int i=0; i<3*nverts; i++)
				newq[i] = s.q[i] + (newtime-s.lastUpdateTime[i])*s.v[i];

			VectorXd F(s.q.size());
			F.setZero();
			bool newtouched = group->addForce(newq, F);
			if(newtouched && newtime < earliestTime_)
			{
				std::cout << newtime << " wasn't suppose to fire before " << earliestTime_ << std::endl;
				exit(0);
			}
			s.v += s.minv.asDiagonal()*F;
			for(int i=0; i<nverts; i++)
			{
				bool touched = false;
				for(int j=0; j<3; j++)
					if(F[3*i+j] != 0.0)
						touched = true;
				if(touched)
				{
					for(int j=0; j<3; j++)
					{
						s.q[3*i+j] = newq[3*i+j];
						s.lastUpdateTime[3*i+j] = newtime;
					}
					history_->addHistory(i, newtime, s.q.segment<3>(3*i), oldhistory_, newtime < earliestTime_);
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

double ActiveLayers::closestDistance(const VectorXd &q, const Mesh &m)
{
	int nverts = m.vertices.size()/3;
	int nfaces = m.faces.cols();
	double closest = std::numeric_limits<double>::infinity();

	// Vertex-face proper	
	for(int i=0; i<nverts; i++)
	{
		for(int j=0; j<nfaces; j++)
		{
			if(m.vertexOfFace(i, j))
				continue;

			double b1, b2, b3;
			closest = min(closest, Distance::vertexFaceDistance(q.segment<3>(3*i), q.segment<3>(3*m.faces.coeff(0, j)), q.segment<3>(3*m.faces.coeff(1, j)), q.segment<3>(3*m.faces.coeff(2, j)), b1, b2, b3).norm());
		}
	}

	// Edge-edge proper
	for(int edge1=0; edge1<3*nfaces; edge1++)
	{
		for(int edge2=edge1+1; edge2<3*nfaces; edge2++)
		{
			int face1 = edge1/3;
			int face2 = edge2/3;
			int e1v1 = m.faces.coeff(edge1%3, face1);
			int e1v2 = m.faces.coeff((edge1+1)%3, face1);
			int e2v1 = m.faces.coeff(edge2%3, face2);
			int e2v2 = m.faces.coeff((edge2+1)%3, face2);
			if(e1v1 == e2v1 || e1v1 == e2v2 || e1v2 == e2v1 || e1v2 == e2v2)
				continue;
			double b1, b2, b3, b4;
			closest = min(closest, Distance::edgeEdgeDistance(q.segment<3>(3*e1v1), q.segment<3>(3*e1v2), q.segment<3>(3*e2v1), q.segment<3>(3*e2v2), b1, b2, b3, b4).norm());			
		}
	}

	return closest;
}

bool ActiveLayers::collisionDetection(const Mesh &m, set<VertexFaceStencil> &vfsToAdd, set<EdgeEdgeStencil> &eesToAdd, double &earliestTime)
{
	vfsToAdd.clear();
	eesToAdd.clear();
	earliestTime = 1.0;
	
	bp_->findCollisionCandidates(*history_, m, outerEta_, vfsToAdd, eesToAdd);
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
	np_->findCollisions(*history_, etavfs, etaees, vfsToAdd, eesToAdd, earliestTime);

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

	delete oldhistory_;
	oldhistory_ = history_;
	history_ = new History(s.q);

	while(!step(s));

	if(verbose_)
		std::cout << "Done simulating, accumulated " << history_->countHistoryEntries() << " history entries" << std::endl;

	set<VertexFaceStencil> vfsToAdd;
	set<EdgeEdgeStencil> eesToAdd;

	double t = 0;
	bool collisions = collisionDetection(m, vfsToAdd, eesToAdd, t);
	if(verbose_)
		std::cout << "Found " << vfsToAdd.size() << " vertex-face and " << eesToAdd.size() << " edge-edge collisions, earliest at t=" << t << std::endl;

	if(t < earliestTime_)
	{
		std::cout << "New earliest time " << t << " earlier than old " << earliestTime_ << std::endl;
		exit(0);
	}
	earliestTime_ = t;

	rollback();

	for(std::set<VertexFaceStencil>::iterator it = vfsToAdd.begin(); it != vfsToAdd.end(); ++it)
		addVFStencil(*it);
	for(std::set<EdgeEdgeStencil>::iterator it = eesToAdd.begin(); it != eesToAdd.end(); ++it)
		addEEStencil(*it);

	return !collisions;
}
