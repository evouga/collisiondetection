#include "VelocityFilter.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "Distance.h"
#include <Eigen/Core>
#include "CTCD.h"

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

double displaceInward(const VectorXd &q, const Matrix3Xi &f, const Matrix3Xd &vertNormals, double maxDisplacement, VectorXd &displacedq)
{
	int nverts = q.size()/3;
	int nfaces = f.cols();
	displacedq = q;
	double displacement = maxDisplacement;

	for(int i=0; i<nverts; i++)
	{
		Vector3d pos = q.segment<3>(3*i);
		Vector3d normal = vertNormals.col(i);
		Vector3d endpos = pos - maxDisplacement*normal;

		// vertex-face

		for(int j=0; j<nfaces; j++)
		{
			bool skip = false;
			for(int k=0; k<3; k++)
				if(f.col(j)[k] == i)
					skip = true;

			if(skip)
				continue;

			double t;
			if(CTCD::vertexFaceCTCD(pos, 
					q.segment<3>(3*f.col(j)[0]), 
					q.segment<3>(3*f.col(j)[1]), 
					q.segment<3>(3*f.col(j)[2]),
					endpos,					
					q.segment<3>(3*f.col(j)[0]), 
					q.segment<3>(3*f.col(j)[1]), 
					q.segment<3>(3*f.col(j)[2]),
					1e-8, t))
			{
				displacement = min(displacement, maxDisplacement*t/2);
			}
		}
		// edge-edge

		for(int j=0; j<nfaces; j++)
		{
			for(int k=0; k<3; k++)
			{
				if(f.col(j)[k] == i)
				{
					int vert2 = f.col(j)[(k+1)%3];
					Vector3d pos2 = q.segment<3>(3*vert2);
					Vector3d endpos2 = pos2 - maxDisplacement*vertNormals.col(vert2);
					for(int j2=0; j2<nfaces; j2++)
					{
						for(int k2=0; k2<3; k2++)
						{
							if(f.col(j2)[k2] == i || f.col(j2)[(k2+1)%3] == i)
								continue;
							if(f.col(j2)[k2] == vert2 || f.col(j2)[(k2+1)%3] == vert2)
								continue;

							double t;
							if(CTCD::edgeEdgeCTCD(pos, pos2,
								q.segment<3>(3*f.col(j2)[k2]),
								q.segment<3>(3*f.col(j2)[(k2+1)%3]),
								endpos, endpos2,
								q.segment<3>(3*f.col(j2)[k2]),
								q.segment<3>(3*f.col(j2)[(k2+1)%3]),
								1e-8, t))
							{
								displacement = min(displacement, maxDisplacement*t/2);
							}
						}
					}
				}
			}
		}
	}
	
	for(int i=0; i<nverts; i++)
	{
		displacedq.segment<3>(3*i) -= displacement*vertNormals.col(i);
	}

	return displacement;
}

int main(int argc, char *argv[])
{
	if(argc != 4)
	{
		std::cerr << "Usage: displaceInward (inputMesh) (outputMesh) (maxDisplacement)" << std::endl;
		return -1;
	}
	VectorXd q;
	Matrix3Xi f;
	loadMesh(argv[1], q, f);
	int numverts = q.size()/3;
	std::cout << "Loaded mesh with " << numverts << " vertices and " << f.cols() << " faces" << std::endl;

	double maxdisp = strtod(argv[3], NULL);

	Matrix3Xd normals(3, numverts);
	normals.setZero();
	for(int face=0; face<f.cols(); face++)
	{
		Vector3i verts = f.col(face);
		Vector3d n = ( q.segment<3>(3*verts[1]) - q.segment<3>(3*verts[0]) ).cross( q.segment<3>(3*verts[2]) - q.segment<3>(3*verts[0]) );
		n /= n.norm();

		for(int i=0; i<3; i++)
		{
			double sinangle = ( q.segment<3>(3*verts[(i+1)%3]) - q.segment<3>(3*verts[i]) ).cross( q.segment<3>(3*verts[(i+2)%3]) - q.segment<3>(3*verts[i]) ).norm();
			double cosangle = ( q.segment<3>(3*verts[(i+1)%3]) - q.segment<3>(3*verts[i]) ).dot( q.segment<3>(3*verts[(i+2)%3]) - q.segment<3>(3*verts[i]) );
			double angle = atan2(sinangle, cosangle);
			normals.col(verts[i]) += angle*n;
		}
	}	
	for(int i=0; i<numverts; i++)
		normals.col(i) /= normals.col(i).norm();

	VectorXd newq;
	double mindisp = displaceInward(q, f, normals, maxdisp, newq);	
	std::cout << "Actual displacement was " << mindisp << std::endl;
	writeMesh(argv[2], newq, f);
}
