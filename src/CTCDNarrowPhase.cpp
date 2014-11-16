#include "RetrospectiveDetection.h"
#include "History.h"
#include "CTCD.h"
#include <set>
#include <iostream>

using namespace std;

void CTCDNarrowPhase::findCollisions(const History &h, const set<pair<VertexFaceStencil, double> > &candidateVFS, const set<pair<EdgeEdgeStencil, double> > &candidateEES,
		set<VertexFaceStencil> &vfs, set<EdgeEdgeStencil> &ees, double &earliestTime)
{
	earliestTime = 1.0;
	for(set<pair<VertexFaceStencil, double> >::const_iterator it = candidateVFS.begin(); it != candidateVFS.end(); ++it)
	{
		if(checkVFS(h, it->first, it->second, earliestTime))
			vfs.insert(it->first);
	}
	for(set<pair<EdgeEdgeStencil, double> >::const_iterator it = candidateEES.begin(); it != candidateEES.end(); ++it)
	{
		if(checkEES(h, it->first, it->second, earliestTime))
			ees.insert(it->first);
	}
}

bool CTCDNarrowPhase::checkVFS(const History &h, VertexFaceStencil vfs, double eta, double &earliestTime)
{
	vector<int> verts;
	verts.push_back(vfs.p);
	verts.push_back(vfs.q0);
	verts.push_back(vfs.q1);
	verts.push_back(vfs.q2);
	vector<StitchedEntry> sh;
	h.stitchCommonHistory(verts, sh);

	for(vector<StitchedEntry>::iterator it = sh.begin(); it != sh.end(); ++it)
	{
		vector<StitchedEntry>::iterator next = it+1;
		if(next == sh.end())
			break;

		double tinterval = next->time - it->time;

		bool found = false;

		double t;
		if(CTCD::vertexFaceCTCD(it->pos[0], it->pos[1], it->pos[2], it->pos[3],
						next->pos[0], next->pos[1], next->pos[2], next->pos[3],
						eta, t))
		{
			earliestTime = min(earliestTime, it->time + t*tinterval);
			//std::cout << "VF " << it->time + t*tinterval << " " << vfs.p << " " << vfs.q0 << " " << vfs.q1 << " " << vfs.q2 << std::endl;
			found = true;
		}
			
		// Vertex-face edges
		for(int edge=0; edge<3; edge++)
		{
			if(CTCD::vertexEdgeCTCD(it->pos[0], it->pos[1+(edge%3)], it->pos[1+ ((edge+1)%3)],
						next->pos[0], next->pos[1+(edge%3)], next->pos[1+ ((edge+1)%3)],
						eta, t))
			{
				earliestTime = min(earliestTime, it->time + t*tinterval);
				found = true;
			}
		}
		// Vertex-face vertices
		for(int vert=0; vert<3; vert++)
		{
			if(CTCD::vertexVertexCTCD(it->pos[0], it->pos[1+vert],
						  next->pos[0], next->pos[1+vert],
						  eta, t))
			{
				earliestTime = min(earliestTime, it->time + t*tinterval);
				found = true;
			}
		}
		if(found)
			return true;
	}
	return false;
}

bool CTCDNarrowPhase::checkEES(const History &h, EdgeEdgeStencil ees, double eta, double &earliestTime)
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
		vector<StitchedEntry>::iterator next = it+1;
		if(next == sh.end())
			break;

		double tinterval = next->time - it->time;

		double t;
		bool print = false;
		if(ees.p0 == 2243 && ees.p1 == 2297 && ees.q0 == 11804 && ees.q1 == 12008)
		{
			std::cout << "Checking times " << it->time << " to " << next->time << " eta " << eta << std::endl;

			print = true;
		}
		if(CTCD::edgeEdgeCTCD(it->pos[0], it->pos[1], it->pos[2], it->pos[3],
					next->pos[0], next->pos[1], next->pos[2], next->pos[3],
					      eta, t, print))
		{
			
			if(ees.p0 == 2243 && ees.p1 == 2297 && ees.q0 == 11804 && ees.q1 == 12008)				
				std::cout << "Detected " << ees.p0 << " " << ees.p1 << " " << ees.q0 << " " << ees.q1 << " at " << it->time + t*tinterval << std::endl;
			earliestTime = min(earliestTime, it->time + t*tinterval);				
			return true;
		}
	}
	return false;		
}
