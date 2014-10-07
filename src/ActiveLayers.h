#ifndef ACTIVELAYERS_H
#define ACTIVELAYERS_H

#include "Stencils.h"
#include <vector>
#include <map>
#include <queue>
#include <Eigen/Core>

class PenaltyGroup;
struct SimulationState;

class PenaltyGroupComparator
{
public:
	bool operator()(const PenaltyGroup *first, const PenaltyGroup *second) const;
};

class ActiveLayers
{
public:
	ActiveLayers(double eta, double baseDt, double baseStiffness, double terminationTime);	
	~ActiveLayers();

	void addVFStencil(double curt, VertexFaceStencil stencil);	
	void addEEStencil(double curt, EdgeEdgeStencil stencil);

	bool step(SimulationState &s);

private:
	ActiveLayers(const ActiveLayers &other);
	ActiveLayers &operator=(const ActiveLayers &other);

	void addGroups(double curt, int maxdepth);

	double eta_;
	double baseDt_;
	double baseStiffness_;
	double termTime_;
	
	int deepestLayer_;

	std::vector<PenaltyGroup *> groups_;
	std::priority_queue<PenaltyGroup *, std::vector<PenaltyGroup *>, PenaltyGroupComparator> groupQueue_;
	std::map<VertexFaceStencil, int> vfdepth_;
	std::map<EdgeEdgeStencil, int> eedepth_;
};

#endif
