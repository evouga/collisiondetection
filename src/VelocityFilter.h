#ifndef VELOCITYFILTER_H
#define VELOCITYFILTER_H

#include <Eigen/Core>

class VelocityFilter
{
public:
	const static int DANGLING_VERTICES = -1;
	const static int STARTING_STATE_INTERSECTS = -2;
	const static int MAXROLLBACKS_EXCEEDED = -3;

	static int velocityFilter(const Eigen::VectorXd &qstart, Eigen::VectorXd &qend, const Eigen::Matrix3Xi &faces, const Eigen::VectorXd &masses,
			double outerEta, double innerEta, double baseStiffness = 0.001, double baseSubstepSize = 0.1, int maxRollbacks = 1000);
};

#endif
