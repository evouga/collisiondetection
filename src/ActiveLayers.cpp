#include "ActiveLayers.h"
#include "PenaltyGroup.h"
#include "SimulationState.h"
#include "CTCD.h"
#include "Mesh.h"
#include <set>

using namespace Eigen;

ActiveLayers::ActiveLayers(double eta, double baseDt, double baseStiffness, double terminationTime) : eta_(eta), baseDt_(baseDt), baseStiffness_(baseStiffness), termTime_(terminationTime), deepestLayer_(0)
{
}

ActiveLayers::~ActiveLayers()
{
	for(std::vector<PenaltyGroup *>::iterator it = groups_.begin(); it != groups_.end(); ++it)
		delete *it;
}

void ActiveLayers::addVFStencil(VertexFaceStencil stencil)
{
	int olddepth = vfdepth_[stencil];
	int newdepth = 1 + olddepth;

	addGroups(newdepth);

	groups_[olddepth]->addVFStencil(stencil);
	vfdepth_[stencil]++;	
}

void ActiveLayers::addEEStencil(EdgeEdgeStencil stencil)
{
	int olddepth = eedepth_[stencil];
	int newdepth = 1 + olddepth;

	addGroups(newdepth);

	groups_[olddepth]->addEEStencil(stencil);
	eedepth_[stencil]++;	
}

void ActiveLayers::addGroups(int maxdepth)
{
	while(deepestLayer_ < maxdepth)
	{
		int depth = deepestLayer_+1;
		double ki = baseStiffness_*depth;
		double etai = layerDepth(depth);
		double dti = baseDt_ / sqrt(double(depth));

		PenaltyGroup *newgroup = new PenaltyGroup(dti, etai, ki);
		groups_.push_back(newgroup);
		groupQueue_.push(newgroup);

		deepestLayer_++;
	}
}

bool ActiveLayers::step(SimulationState &s)
{
	if(!groupQueue_.empty())
	{
		PenaltyGroup *group = groupQueue_.top();
		groupQueue_.pop();

		double newtime = group->nextFireTime();
		if(newtime < termTime_)
		{
			double dt = newtime - s.time;
			s.q += dt * s.v;
			VectorXd F(s.q.size());
			group->addForce(s.q, F);
			for(int i=0; i<F.size(); i++)
				F[i] /= s.m[i];
			s.v += F;
			group->incrementTimeStep();
			groupQueue_.push(group);
			return false;
		}
	}

	double dt = termTime_ - s.time;
	s.q += dt * s.v;
	return true;
}

void ActiveLayers::rollback()
{
	for(std::vector<PenaltyGroup *>::iterator it = groups_.begin(); it != groups_.end(); ++it)
		(*it)->rollback();	
}

double ActiveLayers::layerDepth(int layer)
{
	return eta_/sqrt(double(layer));
}

double ActiveLayers::VFStencilThickness(VertexFaceStencil stencil)
{
	std::map<VertexFaceStencil, int>::iterator it = vfdepth_.find(stencil);
	if(it != vfdepth_.end())
		return layerDepth(it->second + 1);
	return layerDepth(1);
}

double ActiveLayers::EEStencilThickness(EdgeEdgeStencil stencil)
{
	std::map<EdgeEdgeStencil, int>::iterator it = eedepth_.find(stencil);
	if(it != eedepth_.end())
		return layerDepth(it->second + 1);
	return layerDepth(1);
}

bool ActiveLayers::collisionDetectionAndResponse(const SimulationState &s, const Mesh &m)
{
	int nverts = m.vertices.size();
	int nfaces = m.faces.size();

	const VectorXd &qold = m.vertices;
	const VectorXd &qnew = s.q;

	double t;

	std::set<VertexFaceStencil> vfsToAdd;
	std::set<EdgeEdgeStencil> eesToAdd;

	// Vertex-face proper	
	for(int i=0; i<nverts; i++)
	{
		for(int j=0; j<nfaces; j++)
		{
			if(m.vertexOfFace(i, j))
				continue;

			VertexFaceStencil vfs(i, m.faces.coeff(0, j), m.faces.coeff(1, j), m.faces.coeff(2, j));
			double depth = VFStencilThickness(vfs);
			bool done = false;
			if(CTCD::vertexFaceCTCD(qold.segment<3>(3*vfs.p), qold.segment<3>(3*vfs.q0), qold.segment<3>(3*vfs.q1), qold.segment<3>(3*vfs.q2),
						qnew.segment<3>(3*vfs.p), qnew.segment<3>(3*vfs.q0), qnew.segment<3>(3*vfs.q1), qnew.segment<3>(3*vfs.q2),
						depth, t))
			{
				vfsToAdd.insert(vfs);
				done = true;
			}
			if(!done)
			{
				// Vertex-face edges
				for(int edge=0; edge<3; edge++)
				{
					if(CTCD::vertexEdgeCTCD(qold.segment<3>(3*vfs.p), qold.segment<3>(3*m.faces.coeff(edge%3, j)), qold.segment<3>(3*m.faces.coeff( (edge+1)%3, j)),
								qnew.segment<3>(3*vfs.p), qnew.segment<3>(3*m.faces.coeff(edge&3, j)), qnew.segment<3>(3*m.faces.coeff( (edge+1)%3, j)),
								depth, t))
					{						
						vfsToAdd.insert(vfs);
						done = true;
						break;
					}
				}
			}
			if(!done)
			{
				// Vertex-face vertices
				for(int vert=0; vert<3; vert++)
				{
					if(CTCD::vertexVertexCTCD(qold.segment<3>(3*vfs.p), qold.segment<3>(3*m.faces.coeff(vert, j)),
								  qnew.segment<3>(3*vfs.p), qnew.segment<3>(3*m.faces.coeff(vert, j)),
								  depth, t))
					{
						vfsToAdd.insert(vfs);
						done = true;
						break;
					}
				}
			}
		}
	}

	// Edge-edge proper
	for(int edge1=0; edge1<3*nfaces; edge1++)
	{
		for(int edge2=edge1+1; edge2<3*nfaces; edge2++)
		{
			int face1 = edge1/3;
			int face2 = edge2/3;
			int e1v1 = m.faces.coeff(edge1%3, face1);
			int e1v2 = m.faces.coeff((edge1+1)%3, face1);
			int e2v1 = m.faces.coeff(edge2%3, face2);
			int e2v2 = m.faces.coeff((edge2+1)%3, face2);
			if(e1v1 == e2v1 || e1v1 == e2v2 || e1v2 == e2v1 || e1v2 == e2v2)
				continue;
			
			EdgeEdgeStencil ees(e1v1, e1v2, e2v1, e2v2);
			double depth = EEStencilThickness(ees);
			if(CTCD::edgeEdgeCTCD(qold.segment<3>(3*ees.p0), qold.segment<3>(3*ees.p1), qold.segment<3>(3*ees.q0), qold.segment<3>(3*ees.q1),
					      qnew.segment<3>(3*ees.p0), qnew.segment<3>(3*ees.p1), qnew.segment<3>(3*ees.q0), qnew.segment<3>(3*ees.q1),
					      depth, t))
			{
				eesToAdd.insert(ees);
			}
			// vertex-edge and vertex-vertex already handled above
		}
	}

	if(vfsToAdd.empty() && eesToAdd.empty())
		return false;

	for(std::set<VertexFaceStencil>::iterator it = vfsToAdd.begin(); it != vfsToAdd.end(); ++it)
		addVFStencil(*it);
	for(std::set<EdgeEdgeStencil>::iterator it = eesToAdd.begin(); it != eesToAdd.end(); ++it)
		addEEStencil(*it);

	return true;
}

bool ActiveLayers::runOneIteration(const Mesh &m, const VectorXd &mass)
{
	SimulationState s;
	s.q = m.vertices;
	s.v.resize(s.q.size());
	s.v.setZero();
	s.m = mass;
	s.time = 0;

	rollback();

	while(!step(s));

	return collisionDetectionAndResponse(s, m);
}
