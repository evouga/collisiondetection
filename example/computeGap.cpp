#include "VelocityFilter.h"
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
		std::cerr << "Usage: computeGap (coarseMesh) (fineMesh)" << std::endl;
		return -1;
	}
	VectorXd qcoarse;
	Matrix3Xi fcoarse;
	loadMesh(argv[1], qcoarse, fcoarse);
	int numcoarseverts = qcoarse.size()/3;
	std::cout << "Loaded coarse mesh with " << numcoarseverts << " vertices and " << fcoarse.cols() << " faces" << std::endl;

	VectorXd qfine;
	Matrix3Xi ffine;
	loadMesh(argv[2], qfine, ffine);
	std::cout << "Loaded fine mesh with " << qfine.size()/3 << " vertices and " << ffine.cols() << " faces" << std::endl;

	for(int i=0; i<fcoarse.cols(); i++)
	{
		for(int j=0; j<ffine.cols(); j++)
		{
			double t;
			for(int k=0; k<3; k++)
			{
				if(CTCD::vertexFaceCTCD(
qcoarse.segment<3>(3*fcoarse.col(i)[k]), 
qfine.segment<3>(3*ffine.col(j)[0]), 
qfine.segment<3>(3*ffine.col(j)[1]), 
qfine.segment<3>(3*ffine.col(j)[2]), 
qcoarse.segment<3>(3*fcoarse.col(i)[(k+1)%3]),
qfine.segment<3>(3*ffine.col(j)[0]),
qfine.segment<3>(3*ffine.col(j)[1]),
qfine.segment<3>(3*ffine.col(j)[2]), 1e-8, t))
				{
					std::cout << "Meshes intersect!" << 
std::endl;
					return 0;
				}
			}
		}
	}

	VectorXd q(qcoarse.size() + qfine.size());
	Matrix3Xi f;
	f.resize(3, fcoarse.cols() + ffine.cols());
	for(int i=0; i<(int)qcoarse.size(); i++)
		q[i] = qcoarse[i];
	for(int i=0; i<qfine.size(); i++)
		q[qcoarse.size()+i] = qfine[i];
	for(int i=0; i<(int)fcoarse.cols(); i++)
		f.col(i) = fcoarse.col(i);
	for(int i=0; i<(int)ffine.cols(); i++)
		for(int j=0; j<3; j++)
			f.coeffRef(j, fcoarse.cols()+i) = ffine.coeff(j,i) + numcoarseverts;

	set<int> fixedVerts;
	for(int i=0; i<(int)qfine.size()/3; i++)
		fixedVerts.insert(qcoarse.size()/3 + i);

	double distance = Distance::meshSelfDistance(q, f, fixedVerts);
	std::cout << "Min distance between coarse mesh and itself/fine mesh: " << distance << std::endl;
}
