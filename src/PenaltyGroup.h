#ifndef PENALTYGROUP_H
#define PENALTYGROUP_H

#include <Eigen/Core>
#include <vector>
#include <set>
#include "Stencils.h"

class VertexFacePenaltyPotential;
class EdgeEdgePenaltyPotential;

class PenaltyGroup
{
public:
	PenaltyGroup(double dt, double outerEta, double innerEta, double stiffness, double CoR);
	~PenaltyGroup();

	void addVFStencil(VertexFaceStencil stencil);	
	void addEEStencil(EdgeEdgeStencil stencil);

	bool addForce(const Eigen::VectorXd &q, const Eigen::VectorXd &v, Eigen::VectorXd &F);
	void incrementTimeStep();
	double nextFireTime() const;
	void rollback();
	double getDt() {return dt_;}
	double getOuterEta() {return outerEta_;}
	const std::set<int> &getGroupStencil() {return groupStencil_;}

private:
	PenaltyGroup(const PenaltyGroup &other);
	PenaltyGroup &operator=(const PenaltyGroup &other);

	std::vector<VertexFaceStencil> vfstencils_;
	std::vector<EdgeEdgeStencil> eestencils_;
	std::set<int> groupStencil_;
	int nextstep_;
	double dt_;
	double outerEta_;
	double innerEta_;
	double stiffness_;
	double CoR_;
};

#endif
