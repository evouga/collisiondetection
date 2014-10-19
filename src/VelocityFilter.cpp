#include "VelocityFilter.h"
#include "Stencils.h"
#include "ActiveLayers.h"
#include "Mesh.h"
#include "SimulationState.h"
#include <iostream>
#include <set>
#include <fstream>
#include <sstream>

using namespace Eigen;
using namespace std;


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


int VelocityFilter::velocityFilter(const VectorXd &qstart, VectorXd &qend, const Matrix3Xi &faces, const VectorXd &invmasses,
			double outerEta, double innerEta, double baseStiffness, double baseSubstepSize, int maxRollbacks)
{
	int nverts = qstart.size()/3;
	assert(qend.size() == 3*nverts);
	assert(invmasses.size() == 3*nverts);

	int nfaces = faces.cols();

	for(int i=0; i<nverts; i++)
	{
		bool found = false;
		for(int j=0; j<nfaces; j++)
		{
			for(int k=0; k<3; k++)
				if(faces.coeff(k, j) == i)
					found = true;
		}
		if(!found)
		{
			return DANGLING_VERTICES;
		}
	}

	double maxmass = 0;
	for(int i=0; i<invmasses.size(); i++)
		maxmass = max(maxmass, invmasses[i]);

	double dt = baseSubstepSize;
	double k = baseStiffness/dt*sqrt(maxmass);

	std::cout << "For outer layer thickness " << outerEta << " and max inverse mass " << maxmass << ", using base timestep " << dt << " and stiffness " << k << std::endl;

	ActiveLayers al(outerEta, innerEta, dt, k, 1.0, true);
	Mesh m;
	m.vertices = qstart;
	m.faces = faces;	

	double closest = al.closestDistance(qstart, m);

	std::cout << "Closest primitive pair is distance " << closest << " apart at start" << std::endl;
	
	if(closest < innerEta)
		return STARTING_STATE_INTERSECTS;

	int neededLayers = (outerEta-innerEta)/(closest-innerEta);
	std::cout << "It will take at least " << neededLayers << " layers to resolve this contact problem" << std::endl;
	

	for(int iter=0; iter<maxRollbacks; iter++)
	{
		SimulationState s;
		s.q = qstart;
		s.v = qend-qstart;
		s.minv = invmasses;
		s.time = 0;
		if(al.runOneIteration(m, s))
		{
			qend = s.q;
			return iter;
		}
		stringstream ss;
		ss << "iter_" << iter << ".obj";
		writeMesh(ss.str().c_str(), s.q, m.faces);
	}
	return MAXROLLBACKS_EXCEEDED;
}
