#ifndef MESH_H
#define MESH_H

#include <Eigen/Core>

struct Mesh
{
	Eigen::VectorXd vertices;
	Eigen::Matrix3Xi faces;

	bool neighboringFaces(int face1, int face2) const;
	bool vertexOfFace(int vert, int face) const;
};

#endif
