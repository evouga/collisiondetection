#include "VelocityFilter.h"
#include "Stencils.h"
#include "ActiveLayers.h"
#include "Mesh.h"
#include "SimulationState.h"
#include <iostream>
#include <set>

using namespace Eigen;
using namespace std;

int VelocityFilter::velocityFilter(const VectorXd &qstart, VectorXd &qend, const Matrix3Xi &faces, const VectorXd &masses,
			double outerEta, double innerEta, double baseStiffness, double baseSubstepSize, int maxRollbacks)
{
	int nverts = qstart.size()/3;
	assert(qend.size() == 3*nverts);
	assert(masses.size() == nverts);

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

	double maxvel = 0.0;
	for(int i=0; i<nverts; i++)
	{
		double vel = (qend.segment<3>(3*i) - qstart.segment<3>(3*i) ).norm();
		maxvel = max(maxvel, vel);
	}

	double effmaxval = maxvel;

	if(effmaxval == 0)
		effmaxval = 1.0;

	double dt = baseSubstepSize * outerEta/effmaxval;
	double k = baseStiffness/dt;

	std::cout << "For outer layer thickness " << outerEta << " and maximum velocity " << maxvel << ", using base timestep " << dt << " and stiffness " << k << std::endl;

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
		s.m = masses;
		s.time = 0;
		if(al.runOneIteration(m, s))
		{
			qend = s.q;
			return iter;
		}
	}
	return MAXROLLBACKS_EXCEEDED;
}
