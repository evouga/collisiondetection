#include "VelocityFilter.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "Distance.h"
#include <Eigen/Core>

using namespace std;
using namespace Eigen;

void writeMesh(const char *filename, const VectorXd &verts, const Matrix3Xi &faces)
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

int main()
{
	VectorXd q1, q2;
	Matrix3Xi f1, f2;
	loadMesh("mesh1.obj", q1, f1);
	loadMesh("mesh2.obj", q2, f2);

	std::cout << "Loaded " << q1.size() << " vertices, " << f1.size() << " faces." << std::endl;

	VectorXd masses(q1.size()/3);
	for(int i=0; i<157; i++)
		masses[i] = std::numeric_limits<double>::infinity();
	for(int i=157; i<q1.size()/3; i++)
		masses[i] = 1.0;
	std::cout << VelocityFilter::velocityFilter(q1, q2, f1, masses, 1e-7, 1e-10) << std::endl;

	writeMesh("out.obj", q2, f1);
}
