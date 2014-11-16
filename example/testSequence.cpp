#include "VelocityFilter.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "Distance.h"
#include <Eigen/Core>

using namespace std;
using namespace Eigen;

static void writeMesh(const char *filename, const VectorXd &verts, const Matrix3Xi &faces)
{
	ofstream ofs(filename);

	for(int i=0; i<verts.size()/3; i++)
	{
		ofs << "v ";
		for(int j=0; j<3; j++)
			ofs << verts[3*i+j] << " ";
		ofs << endl;
	}

	for(int i=0; i<faces.cols(); i++)
	{
		ofs << "f ";
		for(int j=0; j<3; j++)
			ofs << faces.coeff(j, i)+1 << " ";
		ofs << endl;
	}
}

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
	if(argc < 4)
	{
		std::cerr << "Usage: testSequence (fineMesh) (coarseMesh1) (coarseMesh2) ..." << std::endl;
		return -1;
	}
	VectorXd qfine;
	Matrix3Xi ffine;
	loadMesh(argv[1], qfine, ffine);
	int numfineverts = qfine.size()/3;
	std::cout << "Loaded coarse mesh with " << numfineverts << " vertices and " << ffine.cols() << " faces" << std::endl;

	int numcoarse = argc-2;

	std::cout << "Will inflate using " << numcoarse << " fine meshes" << std::endl;

	std::cout << "Checking initial fine mesh" << std::endl;
	VectorXd qcoarse;
	Matrix3Xi fcoarse;
	loadMesh(argv[1+numcoarse], qcoarse, fcoarse);
	std::cout << "Loaded fine mesh with " << qcoarse.size()/3 << " vertices and " << fcoarse.cols() << " faces" << std::endl;
	VectorXd q(qfine.size() + qcoarse.size());
	Matrix3Xi f;
	f.resize(3, ffine.cols() + fcoarse.cols());
	for(int i=0; i<(int)qfine.size(); i++)
		q[i] = qfine[i];
	for(int i=0; i<qcoarse.size(); i++)
		q[qfine.size()+i] = qcoarse[i];
	for(int i=0; i<(int)ffine.cols(); i++)
		f.col(i) = ffine.col(i);
	for(int i=0; i<(int)fcoarse.cols(); i++)
		for(int j=0; j<3; j++)
			f.coeffRef(j, ffine.cols()+i) = fcoarse.coeff(j,i) + numfineverts;

	VectorXd invmasses(qfine.size() + qcoarse.size());
	for(int i=0; i<(int)qfine.size(); i++)
		invmasses[i] = 1.0;
	for(int i=0; i<(int)qcoarse.size(); i++)
		invmasses[qfine.size()+i] = 0;	

	set<int> fixedVerts;
	for(int i=0; i<(int)qcoarse.size()/3; i++)
		fixedVerts.insert(qfine.size()/3 + i);

	double distance = Distance::meshSelfDistance(q, f, fixedVerts);
	while(distance < 1e-4)
	{
		std::cout << "Distance is " << distance << ". Inflating..." << std::endl;
		VectorXd qnew = q;
		VelocityFilter::velocityFilter(q, qnew, f, invmasses, 2.0*distance, 0.5*distance);
		q = qnew;
		distance = Distance::meshSelfDistance(q, f, fixedVerts);
	}
	for(int i=0; i<(int)qfine.size(); i++)
		qfine[i] = q[i];

	for(int iter=0; iter<numcoarse-1; iter++)
	{
		std::cout << "Now processing " << argv[1+numcoarse-iter] << " -> " << argv[numcoarse-iter] << std::endl;
		VectorXd qcoarse1, qcoarse2;
		Matrix3Xi fcoarse1, fcoarse2;
		loadMesh(argv[1+numcoarse-iter], qcoarse1, fcoarse1);
		loadMesh(argv[numcoarse-iter], qcoarse2, fcoarse2);
		if(qcoarse1.size() != qcoarse2.size() || fcoarse1.size() != fcoarse2.size())
		{
			std::cout << "Dimensions don't match!" << std::endl;
			return -1;
		}
		std::cout << "Loaded fine mesh with " << qcoarse1.size()/3 << " vertices and " << fcoarse1.cols() << " faces" << std::endl;
		VectorXd q1(qfine.size() + qcoarse1.size());
		VectorXd q2(qfine.size() + qcoarse2.size());
		Matrix3Xi f1;
		f1.resize(3, ffine.cols() + fcoarse1.cols());
		Matrix3Xi f2;
		f2.resize(3, ffine.cols() + fcoarse2.cols());
		for(int i=0; i<(int)qfine.size(); i++)
		{
			q1[i] = qfine[i];
			q2[i] = qfine[i];
		}
		for(int i=0; i<(int)ffine.cols(); i++)
		{
			f1.col(i) = ffine.col(i);
			f2.col(i) = ffine.col(i);
		}
		for(int i=0; i<qcoarse1.size(); i++)
			q1[qfine.size() + i] = qcoarse1[i];
		for(int i=0; i<qcoarse2.size(); i++)
			q2[qfine.size() + i] = qcoarse2[i];
		for(int i=0; i<(int)fcoarse1.cols(); i++)
		{
			for(int j=0; j<3; j++)
				f1.coeffRef(j, ffine.cols() + i) = fcoarse1.coeff(j, i) + numfineverts;
		}
		for(int i=0; i<(int)fcoarse2.cols(); i++)
		{
			for(int j=0; j<3; j++)
				f2.coeffRef(j, ffine.cols() + i) = fcoarse2.coeff(j, i) + numfineverts;
		}
		VectorXd invmasses(qfine.size() + qcoarse1.size());
		for(int i=0; i<(int)qfine.size(); i++)
			invmasses[i] = 1.0;
		for(int i=0; i<(int)qcoarse1.size(); i++)
			invmasses[qfine.size()+i] = 0;		

		VelocityFilter::velocityFilter(q1, q2, f1, invmasses, 1e-4, 1e-6);
		for(int i=0; i<(int)qfine.size(); i++)
			qfine[i] = q2[i];
	
		stringstream ss;
		ss << "out-iter" << iter << ".obj";

		writeMesh(ss.str().c_str(), qfine, ffine);
	}		
}
