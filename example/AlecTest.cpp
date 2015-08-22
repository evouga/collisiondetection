#include "../src/History.h"
#include "../src/KDOPBroadPhase.h"
#include "../src/Mesh.h"
#include "../src/SeparatingPlaneNarrowPhase.h"
#include "../src/CTCDNarrowPhase.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "Distance.h"
#include <Eigen/Core>
#include "CTCD.h"

using namespace std;
using namespace Eigen;

void loadMesh(const char *filename, VectorXd &verts, Matrix3Xi &faces)
{
	ifstream ifs(filename);
	if(!ifs)
		return;

	vector<double> pos;
	vector<int> faceidx;

	while(true)
	{
		char c;
		ifs >> c;
		if(!ifs)
			break;
		if(c == 'v')
		{
			for(int i=0; i<3; i++)
			{
				double p;
				ifs >> p;
				pos.push_back(p);
			}
		}
		if(c == 'f')
		{
			for(int i=0; i<3; i++)
			{
				int vert;
				ifs >> vert;
				faceidx.push_back(vert);
			}
		}
	}

	int nverts = pos.size()/3;
	int nfaces = faceidx.size()/3;

	verts.resize(3*nverts);
	for(int i=0; i<3*nverts; i++)
		verts[i] = pos[i];
	
	faces.resize(3, nfaces);
	for(int i=0; i<nfaces; i++)
		for(int j=0; j<3; j++)
			faces.coeffRef(j, i) = faceidx[3*i+j]-1;
}

int main(int argc, char *argv[])
{
	if(argc != 3)
	{
		std::cerr << "Usage: AlecTest (initial mesh) (final mesh)" << std::endl;
		return -1;
	}
	VectorXd initialv;
	Matrix3Xi faces;
	loadMesh(argv[1], initialv, faces);
	int numcoarseverts = initialv.size()/3;
	std::cout << "Loaded initial mesh with " << numcoarseverts << " vertices and " << faces.cols() << " faces" << std::endl;

	VectorXd finalv;
	loadMesh(argv[2], finalv, faces);
	std::cout << "Loaded final mesh with " << finalv.size()/3 << " vertices and " << faces.cols() << " faces" << std::endl;

  // Build data structures: mesh and history (trajectories)
  Mesh m;
  m.faces = faces;
  // These are needed but ignored
  m.vertices = initialv;
  History h(initialv);
  h.finishHistory( finalv);

  // Broadphase to collect candidate vertex-face and edge-edge collisions
  const double outerEta = 1e-8;
  std::set<VertexFaceStencil> can_vfs;
  std::set<EdgeEdgeStencil> can_ees;
  std::set<std::pair<VertexFaceStencil,double> > can_vfs_eta;
  std::set<std::pair<EdgeEdgeStencil,double> > can_ees_eta;
  // No fixed vertices
  const std::set<int> fixedVerts;
  KDOPBroadPhase().findCollisionCandidates(h,m,outerEta,can_vfs,can_ees,fixedVerts);
  // Narrowphase 
  const double innerEta = 1e-8;
  for(std::set<VertexFaceStencil>::iterator it = can_vfs.begin(); it != can_vfs.end(); ++it)
  {
    can_vfs_eta.insert(std::pair<VertexFaceStencil,double>(*it, innerEta));
  }
  for(std::set<EdgeEdgeStencil>::iterator it = can_ees.begin(); it != can_ees.end(); ++it)
  {
    can_ees_eta.insert(std::pair<EdgeEdgeStencil,double>(*it, innerEta));
  }
  std::set<VertexFaceStencil> ctcdvfs;
  {
    std::set<EdgeEdgeStencil> ees;
    CTCDNarrowPhase().findCollisions(h,can_vfs_eta,can_ees_eta,ctcdvfs,ees);
    std::cout<<"CTCDNarrowPhase:"<<std::endl;
    std::cout<<"  "<<ctcdvfs.size()<<" vf collisions"<<std::endl;
    std::cout<<"  "<<ees.size()<<" ee collisions"<<std::endl;
  }
  {
    std::set<VertexFaceStencil> vfs;
    std::set<EdgeEdgeStencil> ees;
    SeparatingPlaneNarrowPhase().findCollisions(h,can_vfs_eta,can_ees_eta,vfs,ees);
    std::cout<<"SeparatingPlaneNarrowPhase:"<<std::endl;
    std::cout<<"  "<<vfs.size()<<" vf collisions"<<std::endl;
    std::cout<<"  "<<ees.size()<<" ee collisions"<<std::endl;

   for(std::set<VertexFaceStencil>::iterator it = ctcdvfs.begin(); it != ctcdvfs.end(); ++it)
	{
		if(vfs.count(*it) == 0)
			cout << it->p << " " << it->q0 << " " << it->q1 << " " << it->q2 << " " << endl;
	}
  }
}
