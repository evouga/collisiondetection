#ifndef SEPARATINGPLANENARROWPHASE_H
#define SEPARATINGPLANENARROWPHASE_H

#include "RetrospectiveDetection.h"

class SeparatingPlaneNarrowPhase : public NarrowPhase
{
public:
	virtual void findCollisions(const History &h, const std::set<std::pair<VertexFaceStencil, double> > &candidateVFS, const std::set<std::pair<EdgeEdgeStencil, double> > &candidateEES,
		std::set<VertexFaceStencil> &vfs, std::set<EdgeEdgeStencil> &ees);

private:
	bool checkVFS(const History &h, VertexFaceStencil vfs, double eta);
	bool checkEES(const History &h, EdgeEdgeStencil vfs, double eta);
	double planeIntersect(const Eigen::Vector3d &planePosOld, const Eigen::Vector3d &planeVel, const Eigen::Vector3d &planeNormal, const Eigen::Vector3d &ptOld, const Eigen::Vector3d &ptNew, double ptdt, double eta);


	bool checkVFSInterval(const History &h, VertexFaceStencil vfs, double eta, double mint, double maxt);
	bool checkEESInterval(const History &h, EdgeEdgeStencil vfs, double eta, double mint, double maxt);
};

#endif
