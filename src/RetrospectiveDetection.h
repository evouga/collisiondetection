#ifndef RETROSPECTIVEDETECTION_H
#define RETROSPECTIVEDETECTION_H

#include "Stencils.h"
#include "History.h"
#include <set>

class Mesh;

class BroadPhase
{
public:
	virtual void findCollisionCandidates(const History &h, const Mesh &m, double outerEta, std::set<VertexFaceStencil> &vfs, std::set<EdgeEdgeStencil> &ees, const std::set<int> &fixedVerts) = 0;
	virtual ~BroadPhase() {}
};

class NarrowPhase
{
public:
	virtual void findCollisions(const History &h, const std::set<std::pair<VertexFaceStencil, double> > &candidateVFS, const std::set<std::pair<EdgeEdgeStencil, double> > &candidateEES,
		std::set<VertexFaceStencil> &vfs, std::set<EdgeEdgeStencil> &ees) = 0;
	virtual ~NarrowPhase() {}
};

#endif
