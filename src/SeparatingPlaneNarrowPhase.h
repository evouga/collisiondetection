#ifndef SEPARATINGPLANENARROWPHASE_H
#define SEPARATINGPLANENARROWPHASE_H

#include "RetrospectiveDetection.h"
#include "History.h"

class SeparatingPlaneNarrowPhase : public NarrowPhase
{
public:
	virtual void findCollisions(const History &h, const std::set<std::pair<VertexFaceStencil, double> > &candidateVFS, const std::set<std::pair<EdgeEdgeStencil, double> > &candidateEES,
		std::set<VertexFaceStencil> &vfs, std::set<EdgeEdgeStencil> &ees);

private:
	enum StencilType {ST_VFS, ST_EES};

	bool checkVFS(const History &h, VertexFaceStencil vfs, double eta, double eps);
	bool checkEES(const History &h, EdgeEdgeStencil vfs, double eta, double eps);
	double planeIntersect(const Eigen::Vector3d &planePosOld, const Eigen::Vector3d &planeVel, const Eigen::Vector3d &planeNormal, const Eigen::Vector3d &ptOld, const Eigen::Vector3d &ptNew, double ptdt, double eta);

	double planeTrajectoryIntersect(const std::vector<HistoryEntry> &hist, int startidx, bool forward, const Eigen::Vector3d &planePos, const Eigen::Vector3d &planeVel, const Eigen::Vector3d &planeNormal, const Eigen::Vector3d &ptstart, double timestart, double eta);
	bool checkInterval(StencilType type, const History &h, const std::vector<int> &verts, double eta, double mint, double maxt, double eps);
};

#endif
