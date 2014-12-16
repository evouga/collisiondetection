#ifndef VELOCITYFILTER_H
#define VELOCITYFILTER_H

#include <Eigen/Core>

class VelocityFilter
{
public:
	const static int DANGLING_VERTICES = -1;
	const static int STARTING_STATE_INTERSECTS = -2;
	const static int MAXROLLBACKS_EXCEEDED = -3;

	/* Given a mesh with initial vertex positions qstart, triangle faces faces (each column of 
	   faces specifies one triagle), and candidate ending positions qend, resolves collisions by
	   perturbing the ending positions (a velocity filter). qend is overwritten with new ending
	   positions. The mesh does not need to be manifold, nor must the faces be oriented: however,
	   it is assumed that all vertices are part of at least one face (otherwise the algorithm aborts
	   and returns DANGLING_VERTICES.)
	   The collision forces will respect the lumped inverse masses *per DOF* in the masses vector. 
	   You can give some DOFs infinite mass by setting the corresponding entries to zero, and the
	   algorithm will cope (if possible).
	   outerEta and innerEta control the extent of the contact forces: primitives further than outerEta
	   apart will feel no force, and primitives approaching innerEta apart will feel infinite force.
	   The input mesh qstart must not contain any vertex-face or edge-edge pairs that are already
	   closer than innerEta apart. If it does, the algorithm aborts and returns STARTING_STATE_INTERSECTS.
	   The starting positions should also be free of intersections, but the algorithm will not check for this.
           Any collisions present in the initial configuration will be maintained in the output.
	   The other parameters control the physics of the contact response. The coefficient of restitution CoR
           should be strictly positive and introduces damping into the simulation. Set it to 1.0 for perfectly
           elastic response. 
           With default parameters for the step size and stiffness parameters, the algorithm will attempt 
           to guess a good stiffness and substep size given the maximum velocity and minimum mass in the 
           system. You can override this guess by manually adjusting the stiffness and/or substep size 
           parameters. If the stiffness is too large for the time step, the simulation will blow up 
           and never terminate. If too small, it may take many layers and many iterations to find a solution.
	   maxRollbacks controls how many iterations the algorithm will take before giving up and returning
	   MAXROLLBACKS_EXCEEDED. Even if the simulation does not blow up, it is possible for the algorithm
	   to never terminate if the input contains impossible geometry (e.g. due to conflicting constraints
           caused by infinite-mass regions) or if the problem is so hard that the size of the layers reduces
	   to close to the machine epsilon.
	   If the algorithm succeeds, it returns the number of outer iterations that were needed to resolve
	   all contacts. A negative number indicates a detected error (see above).
	*/

	static int velocityFilter(const Eigen::VectorXd &qstart, Eigen::VectorXd &qend, const Eigen::Matrix3Xi &faces, const Eigen::VectorXd &masses,
			double outerEta, double innerEta, double baseStiffness = .0001, double baseSubstepSize = 0.0001, double CoR = 0.1, int maxRollbacks = 1000);
};

#endif
