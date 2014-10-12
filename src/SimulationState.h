#ifndef SIMULATIONSTATE_H
#define SIMULATIONSTATE_H

#include <Eigen/Core>

struct SimulationState
{
	Eigen::VectorXd q;
	Eigen::VectorXd v;
	Eigen::VectorXd minv;

	double time;
};

#endif
