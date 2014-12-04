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
	if(argc != 5)
	{
		std::cerr << "Usage: testNewSequence (outerRadius) (innerRadius) (coarseMesh) (fineMeshBaseName_)" << std::endl;
		return -1;
	}
	VectorXd qcoarse;
	Matrix3Xi fcoarse;
	double outerRadius = strtod(argv[1], NULL);
	double innerRadius = strtod(argv[2], NULL);
	loadMesh(argv[3], qcoarse, fcoarse);
	int numcoarseverts = qcoarse.size()/3;
	std::cout << "Loaded coarse mesh with " << numcoarseverts << " vertices and " << fcoarse.cols() << " faces" << std::endl;

	std::cout << "Checking initial fine mesh" << std::endl;
	VectorXd qfine;
	Matrix3Xi ffine;
	stringstream initial;
	initial << argv[4];
	initial << "0.obj";
	loadMesh(initial.str().c_str(), qfine, ffine);
	std::cout << "Loaded fine mesh with " << qfine.size()/3 << " vertices and " << ffine.cols() << " faces" << std::endl;
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

	VectorXd invmasses(qcoarse.size() + qfine.size());
	for(int i=0; i<(int)qcoarse.size(); i++)
		invmasses[i] = 1.0;
	for(int i=0; i<(int)qfine.size(); i++)
		invmasses[qcoarse.size()+i] = 0;	

	set<int> fixedVerts;
	for(int i=0; i<(int)qfine.size()/3; i++)
		fixedVerts.insert(qcoarse.size()/3 + i);

	double distance = Distance::meshSelfDistance(q, f, fixedVerts);
	while(distance < outerRadius)
	{
		std::cout << "Distance is " << distance << ". Inflating..." << std::endl;
		VectorXd qnew = q;
		VelocityFilter::velocityFilter(q, qnew, f, invmasses, 2.0*distance, 0.5*distance);
		q = qnew;
		distance = Distance::meshSelfDistance(q, f, fixedVerts);
	}
	for(int i=0; i<(int)qcoarse.size(); i++)
		qcoarse[i] = q[i];

	int curmesh = 0;
	while(true)
	{
		stringstream curname;
		curname << argv[4];
		curname << curmesh;
		curname << ".obj";

		stringstream nextname;
		nextname << argv[4];
		nextname << curmesh+1;
		nextname << ".obj";

		std::cout << "Now processing " << curname.str() << " -> " << nextname.str() << std::endl;
		VectorXd qfine1, qfine2;
		Matrix3Xi ffine1, ffine2;
		loadMesh(curname.str().c_str(), qfine1, ffine1);
		loadMesh(nextname.str().c_str(), qfine2, ffine2);
		if(qfine1.size() != qfine2.size() || ffine1.size() != ffine2.size())
		{
			if(qfine2.size() == 0)
			{
				std::cout << "No " << nextname.str() << ", terminating.";
				return 0;
			}
			std::cout << "Dimensions don't match!" << std::endl;
			return -1;
		}
		std::cout << "Loaded fine mesh with " << qfine1.size()/3 << " vertices and " << ffine1.cols() << " faces" << std::endl;
		VectorXd q1(qcoarse.size() + qfine1.size());
		VectorXd q2(qcoarse.size() + qfine2.size());
		Matrix3Xi f1;
		f1.resize(3, fcoarse.cols() + ffine1.cols());
		Matrix3Xi f2;
		f2.resize(3, fcoarse.cols() + ffine2.cols());
		for(int i=0; i<(int)qcoarse.size(); i++)
		{
			q1[i] = qcoarse[i];
			q2[i] = qcoarse[i];
		}
		for(int i=0; i<(int)fcoarse.cols(); i++)
		{
			f1.col(i) = fcoarse.col(i);
			f2.col(i) = fcoarse.col(i);
		}
		for(int i=0; i<qfine1.size(); i++)
			q1[qcoarse.size() + i] = qfine1[i];
		for(int i=0; i<qfine2.size(); i++)
			q2[qcoarse.size() + i] = qfine2[i];
		for(int i=0; i<(int)ffine1.cols(); i++)
		{
			for(int j=0; j<3; j++)
				f1.coeffRef(j, fcoarse.cols() + i) = ffine1.coeff(j, i) + numcoarseverts;
		}
		for(int i=0; i<(int)ffine2.cols(); i++)
		{
			for(int j=0; j<3; j++)
				f2.coeffRef(j, fcoarse.cols() + i) = ffine2.coeff(j, i) + numcoarseverts;
		}
		VectorXd invmasses(qcoarse.size() + qfine1.size());
		for(int i=0; i<(int)qcoarse.size(); i++)
			invmasses[i] = 1.0;
		for(int i=0; i<(int)qfine1.size(); i++)
			invmasses[qcoarse.size()+i] = 0;		

		int ret = VelocityFilter::velocityFilter(q1, q2, f1, invmasses, outerRadius, innerRadius);
		if(ret < 0)
		{
			std::cout << "Velocity filter failed! Error code: " << ret << std::endl;
			return -1;
		}
		for(int i=0; i<(int)qcoarse.size(); i++)
			qcoarse[i] = q2[i];
	
		stringstream ss;
		ss << "out-iter-" << curmesh << ".obj";

		writeMesh(ss.str().c_str(), qcoarse, fcoarse);
		curmesh++;
	}		
}
