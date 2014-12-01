#include "SeparatingPlaneNarrowPhase.h"
#include "History.h"
#include "CTCD.h"
#include "Distance.h"
#include <set>
#include <iostream>

using namespace std;
using namespace Eigen;

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
	if(maxt < mint)
		return false;

	const double eps = 1e-8;
	vector<int> verts;
	verts.push_back(vfs.p);
	verts.push_back(vfs.q0);
	verts.push_back(vfs.q1);
	verts.push_back(vfs.q2);
	vector<StitchedEntry> sh;
	h.stitchCommonHistoryInterval(verts, sh, mint, maxt);
	int numentries = sh.size();

	if(maxt-mint < 3*eps && numentries == 2)
	{
		vector<Vector3d> oldpos;
		vector<Vector3d> newpos;
		for(int vert=0; vert<4; vert++)
		{
			oldpos.push_back(h.getPosAtTime(verts[vert], mint));
			newpos.push_back(h.getPosAtTime(verts[vert], maxt));
		}
		double t;
		return CTCD::vertexFaceCTCD(oldpos[0], oldpos[1], oldpos[2], oldpos[3], newpos[0], newpos[1], newpos[2], newpos[3], eta, t);
	}

	double midt = 0.5*(mint+maxt);
	vector<Eigen::Vector3d> midpos;
	for(int i=0; i<4; i++)
		midpos.push_back(h.getPosAtTime(verts[i], midt));

	double bary1, bary2, bary3;
	Vector3d closest = Distance::vertexFaceDistance(midpos[0], midpos[1], midpos[2], midpos[3], bary1, bary2, bary3);
	double distsq = closest.squaredNorm();	
	if(distsq < eta*eta)
	{
		return true;
	}

	// find wrapping history entries


	int maxhist = 0;

	while(sh[maxhist].time < midt)
		maxhist++;
	int minhist = maxhist-1;

	Vector3d planepos = 0.5*(midpos[0] + bary1*midpos[1] + bary2*midpos[2] + bary3*midpos[3]);
	Vector3d planevel = 0.5*(sh[maxhist].pos[0]-sh[minhist].pos[0] + bary1*(sh[maxhist].pos[1]-sh[minhist].pos[1]) + bary2*(sh[maxhist].pos[2]-sh[minhist].pos[2]) + bary3*(sh[maxhist].pos[3]-sh[minhist].pos[3]) )/(sh[maxhist].time - sh[minhist].time);

	closest /= closest.norm();
	
	// find upper bound of lower interval
	bool recurselower = false;
	double dt = midt - sh[minhist].time;
	double lowert = dt;
	for(int vert=0; vert<4; vert++)
	{
		double sign = (vert == 0 ? -1.0 : 1.0);
		Vector3d posnew = sh[minhist].pos[vert];
		double t = planeIntersect(planepos, -planevel, sign*closest, midpos[vert], posnew, dt, eta);
		if(t >= 0 && t <= dt)
		{
			recurselower = true;
			lowert = min(lowert, t);
		}
	}
	if(recurselower)
	{
		double newmaxt = midt - lowert + eps;
		if(checkVFSInterval(h, vfs, eta, mint, newmaxt))
			return true;
	}
	else
	{
		for(int i=minhist-1; i >= 0; i--)
		{
			recurselower = false;
			double dt = sh[i+1].time - sh[i].time;
			lowert = dt;
			for(int vert=0; vert<4; vert++)
			{
				Vector3d posold = sh[i+1].pos[vert];
				Vector3d posnew = sh[i].pos[vert];
				double sign = (vert==0 ? -1.0 : 1.0);
				double t = planeIntersect(planepos, -planevel, sign*closest, posold, posnew, dt, eta);
				if(t >= 0 && t <= dt)
				{
					recurselower = true;
					lowert = min(lowert, t);
				}
			}
			if(recurselower)
			{
				double newmaxt = sh[i+1].time - lowert + eps;
				if(checkVFSInterval(h, vfs, eta, mint, newmaxt))
					return true;
				break;
			}
		}
	}
	
	// find lower bound of upper interval
	bool recurseupper = false;
	dt = sh[maxhist].time - midt;
	double uppert = dt;
	for(int vert=0; vert<4; vert++)
	{
		double sign = (vert == 0 ? -1.0 : 1.0);
		Vector3d posnew = sh[maxhist].pos[vert];
		double t = planeIntersect(planepos, planevel, sign*closest, midpos[vert], posnew, dt, eta);
		if(t >= 0 && t <= dt)
		{
			recurseupper = true;
			uppert = min(uppert, t);
		}
	}
	if(recurseupper)
	{
		double newmint = midt + uppert - eps;
		if(checkVFSInterval(h, vfs, eta, newmint, maxt))
			return true;
	}
	else
	{
		for(int i=maxhist; i < numentries-1; i++)
		{
			recurseupper = false;
			double dt = sh[i+1].time - sh[i].time;
			uppert = dt;
			for(int vert=0; vert<4; vert++)
			{
				Vector3d posold = sh[i+1].pos[vert];
				Vector3d posnew = sh[i].pos[vert];
				double sign = (vert==0 ? -1.0 : 1.0);
				double t = planeIntersect(planepos, planevel, sign*closest, posold, posnew, dt, eta);
				if(t >= 0 && t <= dt)
				{
					recurseupper = true;
					uppert = min(uppert, t);
				}
			}
			if(recurseupper)
			{
				double newmint = sh[i].time + uppert - eps;
				if(checkVFSInterval(h, vfs, eta, newmint, maxt))
					return true;
				break;
			}
		}
	}	

	return false;
}

bool SeparatingPlaneNarrowPhase::checkEESInterval(const History &h, EdgeEdgeStencil ees, double eta, double mint, double maxt)
{	
	if(maxt < mint)
		return false;

	const double eps = 1e-8;
	vector<int> verts;
	verts.push_back(ees.p0);
	verts.push_back(ees.p1);
	verts.push_back(ees.q0);
	verts.push_back(ees.q1);
	vector<StitchedEntry> sh;
	h.stitchCommonHistoryInterval(verts, sh, mint, maxt);
	int numentries = sh.size();

	if(maxt-mint < 3*eps && numentries == 2)
	{
		vector<Vector3d> oldpos;
		vector<Vector3d> newpos;
		for(int vert=0; vert<4; vert++)
		{
			oldpos.push_back(h.getPosAtTime(verts[vert], mint));
			newpos.push_back(h.getPosAtTime(verts[vert], maxt));
		}
		double t;
		return CTCD::edgeEdgeCTCD(oldpos[0], oldpos[1], oldpos[2], oldpos[3], newpos[0], newpos[1], newpos[2], newpos[3], eta, t);
	}

	double midt = 0.5*(mint+maxt);
	vector<Eigen::Vector3d> midpos;
	for(int i=0; i<4; i++)
		midpos.push_back(h.getPosAtTime(verts[i], midt));

	double bary1, bary2, bary3, bary4;
	Vector3d closest = Distance::edgeEdgeDistance(midpos[0], midpos[1], midpos[2], midpos[3], bary1, bary2, bary3, bary4);
	double distsq = closest.squaredNorm();	
	if(distsq < eta*eta)
	{
		return true;
	}

	// find wrapping history entries

	int maxhist = 0;

	while(sh[maxhist].time < midt)
		maxhist++;
	int minhist = maxhist-1;

	Vector3d planepos = 0.5*(bary1*midpos[0] + bary2*midpos[1] + bary3*midpos[2] + bary4*midpos[3]);
	Vector3d planevel = 0.5*(bary1*(sh[maxhist].pos[0]-sh[minhist].pos[0]) + bary2*(sh[maxhist].pos[1]-sh[minhist].pos[1]) + bary3*(sh[maxhist].pos[2]-sh[minhist].pos[2]) + bary4*(sh[maxhist].pos[3]-sh[minhist].pos[3]) )/(sh[maxhist].time - sh[minhist].time);

	closest /= closest.norm();
	
	// find upper bound of lower interval
	bool recurselower = false;
	double dt = midt - sh[minhist].time;
	double lowert = dt;
	for(int vert=0; vert<4; vert++)
	{
		double sign = (vert < 2 ? -1.0 : 1.0);
		Vector3d posnew = sh[minhist].pos[vert];
		double t = planeIntersect(planepos, -planevel, sign*closest, midpos[vert], posnew, dt, eta);
		if(t >= 0 && t <= dt)
		{
			recurselower = true;
			lowert = min(lowert, t);
		}
	}
	if(recurselower)
	{
		double newmaxt = midt - lowert + eps;
		if(checkEESInterval(h, ees, eta, mint, newmaxt))
			return true;
	}
	else
	{
		for(int i=minhist-1; i >= 0; i--)
		{
			recurselower = false;
			double dt = sh[i+1].time - sh[i].time;
			lowert = dt;
			for(int vert=0; vert<4; vert++)
			{
				Vector3d posold = sh[i+1].pos[vert];
				Vector3d posnew = sh[i].pos[vert];
				double sign = (vert < 2 ? -1.0 : 1.0);
				double t = planeIntersect(planepos, -planevel, sign*closest, posold, posnew, dt, eta);
				if(t >= 0 && t <= dt)
				{
					recurselower = true;
					lowert = min(lowert, t);
				}
			}
			if(recurselower)
			{
				double newmaxt = sh[i+1].time - lowert + eps;
				if(checkEESInterval(h, ees, eta, mint, newmaxt))
					return true;
				break;
			}
		}
	}
	
	// find lower bound of upper interval
	bool recurseupper = false;
	dt = sh[maxhist].time - midt;
	double uppert = dt;
	for(int vert=0; vert<4; vert++)
	{
		double sign = (vert < 2 ? -1.0 : 1.0);
		Vector3d posnew = sh[maxhist].pos[vert];
		double t = planeIntersect(planepos, planevel, sign*closest, midpos[vert], posnew, dt, eta);
		if(t >= 0 && t <= dt)
		{
			recurseupper = true;
			uppert = min(uppert, t);
		}
	}
	if(recurseupper)
	{
		double newmint = midt + uppert - eps;
		if(checkEESInterval(h, ees, eta, newmint, maxt))
			return true;
	}
	else
	{
		for(int i=maxhist; i < numentries-1; i++)
		{
			recurseupper = false;
			double dt = sh[i+1].time - sh[i].time;
			uppert = dt;
			for(int vert=0; vert<4; vert++)
			{
				Vector3d posold = sh[i+1].pos[vert];
				Vector3d posnew = sh[i].pos[vert];
				double sign = (vert < 2 ? -1.0 : 1.0);
				double t = planeIntersect(planepos, planevel, sign*closest, posold, posnew, dt, eta);
				if(t >= 0 && t <= dt)
				{
					recurseupper = true;
					uppert = min(uppert, t);
				}
			}
			if(recurseupper)
			{
				double newmint = sh[i].time + uppert - eps;
				if(checkEESInterval(h, ees, eta, newmint, maxt))
					return true;
				break;
			}
		}
	}	

	return false;
}

bool SeparatingPlaneNarrowPhase::checkEES(const History &h, EdgeEdgeStencil ees, double eta)
{
	return checkEESInterval(h, ees, eta, 0.0, 1.0);
}

double SeparatingPlaneNarrowPhase::planeIntersect(const Eigen::Vector3d &planePos, const Eigen::Vector3d &planeVel, const Eigen::Vector3d &planeNormal, const Eigen::Vector3d &ptOld, const Eigen::Vector3d &ptNew, double ptdt, double eta)
{
	double numerator = 0.5*eta - (ptOld-planePos).dot(planeNormal);
	double denom = -planeVel.dot(planeNormal);
	if(ptdt != 0.0)
	{
		denom += (ptNew-ptOld).dot(planeNormal)/ptdt;
	}
	return numerator/denom;
}
