#ifndef RETROSPECTIVEDETECTION_H
#define RETROSPECTIVEDETECTION_H

#include "Stencils.h"
#include "History.h"
#include <set>

class Mesh;

class BroadPhase
{
public:
	virtual void findCollisionCandidates(const History &h, const Mesh &m, double outerEta, std::set<VertexFaceStencil> &vfs, std::set<EdgeEdgeStencil> &ees) = 0;
	virtual ~BroadPhase() {}
};

class NarrowPhase
{
public:
	virtual void findCollisions(const History &h, const std::set<std::pair<VertexFaceStencil, double> > &candidateVFS, const std::set<std::pair<EdgeEdgeStencil, double> > &candidateEES,
		std::set<VertexFaceStencil> &vfs, std::set<EdgeEdgeStencil> &ees, double &earliestTime) = 0;
	virtual ~NarrowPhase() {}
};

class TrivialBroadPhase : public BroadPhase
{
public:
	virtual void findCollisionCandidates(const History &h, const Mesh &m, double outerEta, std::set<VertexFaceStencil> &vfs, std::set<EdgeEdgeStencil> &ees);
};

class CTCDNarrowPhase : public NarrowPhase
{
public:
	virtual void findCollisions(const History &h, const std::set<std::pair<VertexFaceStencil, double> > &candidateVFS, const std::set<std::pair<EdgeEdgeStencil, double> > &candidateEES,
		std::set<VertexFaceStencil> &vfs, std::set<EdgeEdgeStencil> &ees, double &earliestTime);

private:
	bool checkVFS(const History &h, VertexFaceStencil vfs, double eta, double &earliestTime);
	bool checkEES(const History &h, EdgeEdgeStencil vfs, double eta, double &earliestTime);
};

#endif
