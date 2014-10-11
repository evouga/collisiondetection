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
			double eta, double baseStiffness, double baseSubstepSize, int maxRollbacks)
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

	double dt = baseSubstepSize * eta/effmaxval;
	double k = baseStiffness/dt;

	std::cout << "For layer thickness " << eta << " and maximum velocity " << maxvel << ", using base timestep " << dt << " and stiffness " << k << std::endl;

	ActiveLayers al(eta, dt, k, 1.0, true);
	Mesh m;
	m.vertices = qstart;
	m.faces = faces;

	set<VertexFaceStencil> vfs;
	set<EdgeEdgeStencil> ees;
	if(al.collisionDetection(qstart, m, vfs, ees))
	{
		return STARTING_STATE_INTERSECTS;
	}

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
