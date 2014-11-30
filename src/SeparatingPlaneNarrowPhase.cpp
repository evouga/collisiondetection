#include "SeparatingPlaneNarrowPhase.h"
#include "History.h"
#include "CTCD.h"
#include <set>
#include <iostream>

using namespace std;

void SeparatingPlaneNarrowPhase::findCollisions(const History &h, const set<pair<VertexFaceStencil, double> > &candidateVFS, const set<pair<EdgeEdgeStencil, double> > &candidateEES,
		set<VertexFaceStencil> &vfs, set<EdgeEdgeStencil> &ees)
{
	for(set<pair<VertexFaceStencil, double> >::const_iterator it = candidateVFS.begin(); it != candidateVFS.end(); ++it)
	{
		if(checkVFS(h, it->first, it->second))
			vfs.insert(it->first);
	}
	for(set<pair<EdgeEdgeStencil, double> >::const_iterator it = candidateEES.begin(); it != candidateEES.end(); ++it)
	{
		if(checkEES(h, it->first, it->second))
			ees.insert(it->first);
	}
}

bool SeparatingPlaneNarrowPhase::checkVFS(const History &h, VertexFaceStencil vfs, double eta)
{
	return checkVFSInterval(h, vfs, eta, 0, 1.0);
}

bool SeparatingPlaneNarrowPhase::checkVFSInterval(const History &h, VertexFaceStencil vfs, double eta, double mint, double maxt)
{	
	const double eps = 1e-8;
	vector<int> verts;
	verts.push_back(vfs.p);
	verts.push_back(vfs.q0);
	verts.push_back(vfs.q1);
	verts.push_back(vfs.q2);
	vector<StitchedEntry> sh;
	h.stitchCommonHistory(verts, sh);

	double midt = 0.5*(mint+maxt);
	vector<Eigen::Vector3d> midpos;
	for(int i=0; i<4; i++)
		midpos.push_back(h.getPosAtTime(verts[i], midt));

	double bary1, bary2, bary3;
	Vector3d closest = Distance::vertexFaceDist(verts[0], verts[1], verts[2], verts[3], bary1, bary2, bary3);
	double distsq = closest.normSquared();
	if(distsq < eta*eta)
		return true;

	// find wrapping history entries

	int numentries = sh.size();

	int maxhist = 0;

	while(sh[maxhist].time < midt)
		maxhist++;
	int minhist = maxhist-1;

	Vector3d planepos = 0.5*(midpos[0] + bary1*midpos[1] + bary2*midpos[2] + bary3*midpos[3]);
	Vector3d planevel = 0.5*(sh[maxhist].pos[0]-sh[minhist].pos[0] + bary1*(sh[maxhist].pos[1]-sh[minhist].pos[1]) + bary2*(sh[maxhist].pos[2]-sh[minhist].pos[2]) + bary3*(sh[maxhist].pos[3]-sh[minhist].pos[3]) );
	
	// find upper bound of lower interval
	double lowert = 0;
	bool recurselower = false;
	for(int vert=0; vert<4; vert++)
	{
		Vector3d posnew = sh[minhist].pos[vert];
		double t;
		if(planeIntersect(planepos, planevel, closest, midpos[vert], posnew, eta, t))
		{
			recurselower = true;
			lowert = min(lowert, t);
		}
	}
	if(recurselower)
	{
		double dt = lowert*(midt - sh[minhist].time);
		double newmaxt = midt - dt + eps;
		checkVFSInterval(h, vfs, eta, mint, newmaxt);
	}
	else
	{
	for(int i=minhist; i >= 0; i--)
	{
		for(int vert=0; vert<4; vert++)
		{
			Vector3d posnew = sh[i].pos[vert];
			double t;
			if(planeIntersect(planepos, planevel, closest, midpos[vert], posnew, t))
			{
				
			}
		}
	}
	
	
	

	for(vector<StitchedEntry>::iterator it = sh.begin(); it != sh.end(); ++it)
	{
		
	}
	return false;
}

bool SeparatingPlaneNarrowPhase::checkEES(const History &h, EdgeEdgeStencil ees, double eta)
{
	vector<int> verts;
	verts.push_back(ees.p0);
	verts.push_back(ees.p1);
	verts.push_back(ees.q0);
	verts.push_back(ees.q1);
	vector<StitchedEntry> sh;
	h.stitchCommonHistory(verts, sh);

	for(vector<StitchedEntry>::iterator it = sh.begin(); it != sh.end(); ++it)
	{
	}
	return false;		
}

bool SeparatingPlaneNarrowPhase::planeIntersect(const Eigen::Vector3d &planePos, const Eigen::Vector3d &planeVel, const Eigen::Vector3d &planeNormal, const Eigen::Vector3d &ptold, const Eigen::Vector3d &ptnew, double eta, double t)
{
	double a = (ptnew-ptold-planevel).dot(planeNormal);
	double b = (ptold - planePos).dot(planeNormal);
	double c = eta*eta*planeNormal.dot(planeNormal);

	double A = a*a;
	double B = 2*a*b;
	double C = b*b-c;	
}
