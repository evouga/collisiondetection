#ifndef DISTANCE_H
#define DISTANCE_H

#include <Eigen/Core>
#include <algorithm>

class Distance
{
 public:
  // Computes the vector between a point p and the closest point to p on the triangle (q0, q1, q2). Also returns the barycentric coordinates of this closest point on the triangle;
  // q0bary is the barycentric coordinate of q0, etc. (The distance from p to the triangle is the norm of this vector.)
  static Eigen::Vector3d vertexFaceDistance(const Eigen::Vector3d &p, 
					    const Eigen::Vector3d &q0, const Eigen::Vector3d &q1, const Eigen::Vector3d &q2, 
					    double &q0bary, double &q1bary, double &q2bary);

  // Computes the shotest vector between a segment (p0, p1) and segment (q0, q1). Also returns the barycentric coordinates of the closest points on both segments; p0bary is the barycentric
  // coordinate of p0, etc. (The distance between the segments is the norm of this vector).
  static Eigen::Vector3d edgeEdgeDistance(const Eigen::Vector3d &p0, const Eigen::Vector3d &p1,
					  const Eigen::Vector3d &q0, const Eigen::Vector3d &q1,
					  double &p0bary, double &p1bary,
					  double &q0bary, double &q1bary);

 private:
  static double clamp(double u)
  {
    return std::min(1.0, std::max(u, 0.0));
  }
};

#endif
