#include "Distance.h"
#include "Mesh.h"
#include "History.h"
#include "Stencils.h"
#include "AABBBroadPhase.h"
#include <iostream>
#include <set>

using namespace Eigen;
using namespace std;

double Distance::meshSelfDistance(const Eigen::VectorXd &verts, const Eigen::Matrix3Xi &faces, const set<int> &fixedVerts)
{
	int nverts = verts.size()/3;
	int nfaces = faces.cols();
	Mesh m;
	m.vertices = verts;
	m.faces = faces;
	double closest = std::numeric_limits<double>::infinity();

	for(int i=0; i<nverts; i++)
	{
		for(int j=0; j<nfaces; j++)
		{
			if(m.vertexOfFace(i,j))
				continue;
			for(int k=0; k<3; k++)
			{
				double dist = (verts.segment<3>(3*i)-verts.segment<3>(3*m.faces.coeff(k,j))).squaredNorm();
				if(dist < closest)
					closest = dist;
			}
		}
	}
	closest = sqrt(closest);

	History h(m.vertices);
	h.finishHistory(m.vertices);

	set<VertexFaceStencil> vfs;
	set<EdgeEdgeStencil> ees;

	AABBBroadPhase bp;

	bp.findCollisionCandidates(h, m, closest, vfs, ees, fixedVerts);	

	std::cout << "Checking " << vfs.size() << " vertex-face and " << ees.size() << " edge-edge stencils" << std::endl;

	closest = std::numeric_limits<double>::infinity();

	for(set<VertexFaceStencil>::iterator it = vfs.begin(); it != vfs.end(); ++it)
	{
		double t;
		double dist = Distance::vertexFaceDistance(verts.segment<3>(3*it->p), verts.segment<3>(3*it->q0), verts.segment<3>(3*it->q1), verts.segment<3>(3*it->q2), t, t, t).norm();
		if(dist < closest)
			closest = dist;
	}
	for(set<EdgeEdgeStencil>::iterator it = ees.begin(); it != ees.end(); ++it)
	{
		double t;
		double dist = Distance::edgeEdgeDistance(verts.segment<3>(3*it->p0), verts.segment<3>(3*it->p1), verts.segment<3>(3*it->q0), verts.segment<3>(3*it->q1), t, t, t, t).norm();
		if(dist < closest)
			closest = dist;
	}
	return closest;
}
