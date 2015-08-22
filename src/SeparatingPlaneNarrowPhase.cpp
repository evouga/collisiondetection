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
	double eps = h.computeMinimumGap()/4.0;
	for(set<pair<VertexFaceStencil, double> >::const_iterator it = candidateVFS.begin(); it != candidateVFS.end(); ++it)
	{
		if(checkVFS(h, it->first, it->second, eps))
			vfs.insert(it->first);
	}
	for(set<pair<EdgeEdgeStencil, double> >::const_iterator it = candidateEES.begin(); it != candidateEES.end(); ++it)
	{
		if(checkEES(h, it->first, it->second, eps))
			ees.insert(it->first);
	}
}

bool SeparatingPlaneNarrowPhase::checkVFS(const History &h, VertexFaceStencil vfs, double eta, double eps)
{
	vector<int> verts;
	verts.push_back(vfs.p);
	verts.push_back(vfs.q0);
	verts.push_back(vfs.q1);
	verts.push_back(vfs.q2);
	return checkInterval(ST_VFS, h, verts, eta, 0, 1.0, eps);
}

bool SeparatingPlaneNarrowPhase::checkEES(const History &h, EdgeEdgeStencil ees, double eta, double eps)
{
	vector<int> verts;
	verts.push_back(ees.p0);
	verts.push_back(ees.p1);
	verts.push_back(ees.q0);
	verts.push_back(ees.q1);
	return checkInterval(ST_EES, h, verts, eta, 0.0, 1.0, eps);
}

bool SeparatingPlaneNarrowPhase::checkInterval(StencilType type, const History &h, const vector<int> &verts, double eta, double mint, double maxt, double eps)
{
	if(maxt < mint)
		return false;

	if(maxt-mint < 3*eps)
	{
		bool ok = true;
		vector<Vector3d> oldpos(4);
		vector<Vector3d> newpos(4);
			
		for(int i=0; i<4; i++)
		{
			Vector3d dummy;
			int oldidx, newidx;
			h.getPosAtTime(verts[i], mint, oldpos[i], oldidx);
			h.getPosAtTime(verts[i], maxt, newpos[i], newidx);
			if(oldidx != newidx)
			{
				ok = false;
				break;
			}
		}

		if(ok)
		{
			double t;
			if(type == ST_VFS)
			{
				if(CTCD::vertexFaceCTCD(oldpos[0], oldpos[1], oldpos[2], oldpos[3], newpos[0], newpos[1], newpos[2], newpos[3], eta, t))
					return true;
				// Vertex-face edges
				for(int edge=0; edge<3; edge++)
				{
					if(CTCD::vertexEdgeCTCD(oldpos[0], oldpos[1+(edge%3)], oldpos[1+ ((edge+1)%3)],
						newpos[0], newpos[1+(edge%3)], newpos[1+ ((edge+1)%3)],
						eta, t))
					{
						return true;
					}
				}
				// Vertex-face vertices
				for(int vert=0; vert<3; vert++)
				{
					if(CTCD::vertexVertexCTCD(oldpos[0], oldpos[1+vert],
						  newpos[0], newpos[1+vert],
						  eta, t))
					{
						return true;
					}
				}
				return false;
			}
			else
			{
				if(CTCD::edgeEdgeCTCD(oldpos[0], oldpos[1], oldpos[2], oldpos[3], newpos[0], newpos[1], newpos[2], newpos[3], eta, t))
					return true;

				// Edge-edge vertices
				if(CTCD::vertexEdgeCTCD(oldpos[0], oldpos[2], oldpos[3], newpos[0], newpos[2], newpos[3], eta, t))
				{
					return true;
				}
				if(CTCD::vertexEdgeCTCD(oldpos[1], oldpos[2], oldpos[3], newpos[1], newpos[2], newpos[3], eta, t))
				{
					return true;
				}
				if(CTCD::vertexEdgeCTCD(oldpos[2], oldpos[0], oldpos[1], newpos[2], newpos[0], newpos[1], eta, t))
				{
					return true;
				}
				if(CTCD::vertexEdgeCTCD(oldpos[3], oldpos[0], oldpos[1], newpos[3], newpos[0], newpos[1], eta, t))
				{
					return true;
				}

				// edge vertex-edge vertex
				if(CTCD::vertexVertexCTCD(oldpos[0], oldpos[2], newpos[0], newpos[2], eta, t))
				{
					return true;
				}
				if(CTCD::vertexVertexCTCD(oldpos[0], oldpos[3], newpos[0], newpos[3], eta, t))
				{
					return true;
				}
				if(CTCD::vertexVertexCTCD(oldpos[1], oldpos[2], newpos[1], newpos[2], eta, t))
				{
					return true;
				}
				if(CTCD::vertexVertexCTCD(oldpos[1], oldpos[3], newpos[1], newpos[3], eta, t))
				{
					return true;
				}

				return false;
			}
		}

	}

	double midt = 0.5*(mint+maxt);
	vector<Eigen::Vector3d> midpos(4);
	vector<int> mididx(4);
	for(int i=0; i<4; i++)
		h.getPosAtTime(verts[i], midt, midpos[i], mididx[i]);

	double bary1, bary2, bary3, bary4;
	Vector3d closest;
	if(type == ST_VFS) 
	{
		closest = Distance::vertexFaceDistance(midpos[0], midpos[1], midpos[2], midpos[3], bary1, bary2, bary3);
		bary4=0; //compiler STFU
	}
	else
		closest = Distance::edgeEdgeDistance(midpos[0], midpos[1], midpos[2], midpos[3], bary1, bary2, bary3, bary4);
	double distsq = closest.squaredNorm();	
	if(distsq < eta*eta)
	{
		return true;
	}

	// find separating plane

	Vector3d planepos, planevel;	

	Vector3d prevpos[4];
	Vector3d nextpos[4];
	double dts[4];
	for(int i=0; i<4; i++)
	{
		prevpos[i] = h.getVertexHistory(verts[i])[mididx[i]].pos;
		nextpos[i] = h.getVertexHistory(verts[i])[mididx[i]+1].pos;
		dts[i] = h.getVertexHistory(verts[i])[mididx[i]+1].time - h.getVertexHistory(verts[i])[mididx[i]].time;
	}

	if(type == ST_VFS)
	{
		planepos = 0.5*(midpos[0] + bary1*midpos[1] + bary2*midpos[2] + bary3*midpos[3]);
		planevel = 0.5*( (nextpos[0]-prevpos[0])/dts[0] + bary1*(nextpos[1]-prevpos[1])/dts[1] + bary2*(nextpos[2]-prevpos[2])/dts[2] + bary3*(nextpos[3]-prevpos[3])/dts[3] );
	}
	else
	{
		planepos = 0.5*(bary1*midpos[0] + bary2*midpos[1] + bary3*midpos[2] + bary4*midpos[3]);
		planevel = 0.5*(bary1*(nextpos[0]-prevpos[0])/dts[0] + bary2*(nextpos[1]-prevpos[1])/dts[1] + bary3*(nextpos[2]-prevpos[2])/dts[2] + bary4*(nextpos[3]-prevpos[3])/dts[3] );
	}
	closest /= closest.norm();
	
	// find upper bound of lower interval

	double t = std::numeric_limits<double>::infinity();
	for(int vert=0; vert<4; vert++)
	{
		double sign = (vert == 0 || (vert == 1 && type == ST_EES)) ? -1.0 : 1.0;
		t = min(t, planeTrajectoryIntersect(h.getVertexHistory(verts[vert]), mididx[vert], false, planepos, planevel, sign*closest, midpos[vert], midt, eta));
	}
	if(t < 1.0)
	{
		double newmaxt = midt - max(0.0, t - eps);
		if(checkInterval(type, h, verts, eta, mint, newmaxt, eps))
			return true;
	}
	
	// find lower bound of upper interval

	t = std::numeric_limits<double>::infinity();
	for(int vert=0; vert<4; vert++)
	{
		double sign = (vert == 0 || (vert == 1 && type == ST_EES)) ? -1.0 : 1.0;
		t = min(t, planeTrajectoryIntersect(h.getVertexHistory(verts[vert]), mididx[vert]+1, true, planepos, planevel, sign*closest, midpos[vert], midt, eta));		
	}
	if(t < 1.0)
	{
		double newmint = midt + max(0.0, t - eps);
		if(checkInterval(type, h, verts, eta, newmint, maxt, eps))
			return true;
	}

	return false;
}

double SeparatingPlaneNarrowPhase::planeTrajectoryIntersect(const vector<HistoryEntry> &hist, int startidx, bool forward, const Eigen::Vector3d &planePos, const Eigen::Vector3d &planeVel, const Eigen::Vector3d &planeNormal, const Eigen::Vector3d &ptstart, double timestart, double eta)
{	
	// check partial window
	double dt = (forward ? hist[startidx].time - timestart : timestart - hist[startidx].time);
	Vector3d newpos = hist[startidx].pos;
	double t = planeIntersect(planePos, (forward ? 1.0 : -1.0) * planeVel, planeNormal, ptstart, newpos, dt, eta);
	if(t >= 0 && t <= dt)
		return t;

	if(forward)
	{
		for(int i=startidx; i<(int)hist.size(); i++)
		{
			if(i == (int)hist.size()-1)
				return numeric_limits<double>::infinity();
			int next = i+1;
			Vector3d newplanepos = planePos + (hist[i].time - timestart)*planeVel;
			double dt = hist[next].time - hist[i].time;
			double t = planeIntersect(newplanepos, planeVel, planeNormal, hist[i].pos, hist[next].pos, dt, eta);
			if(t >= 0 && t <= dt)
				return hist[i].time + t - timestart;
		}
	}
	else
	{
		for(int i=startidx; i>=0; i--)
		{
			if(i == 0)
				return numeric_limits<double>::infinity();
			int prev = i-1;
			Vector3d newplanepos = planePos - (timestart - hist[i].time)*planeVel;
			double dt = hist[i].time - hist[prev].time;
			double t = planeIntersect(newplanepos, -planeVel, planeNormal, hist[i].pos, hist[prev].pos, dt, eta);
			if(t >= 0 && t <= dt)
				return timestart - (hist[i].time - t);
		}
	}

	//unreachable
	return numeric_limits<double>::infinity();
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
