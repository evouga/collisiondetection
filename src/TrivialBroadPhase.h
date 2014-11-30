#ifndef TRIVIALBROADPHASE_H
#define TRIVIALBROADPHASE_H

#include "RetrospectiveDetection.h"

class TrivialBroadPhase : public BroadPhase
{
public:
	virtual void findCollisionCandidates(const History &h, const Mesh &m, double outerEta, std::set<VertexFaceStencil> &vfs, std::set<EdgeEdgeStencil> &ees, const std::set<int> &fixedVerts);
};

#endif
