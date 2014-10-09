#include "Mesh.h"

using namespace std;

bool Mesh::vertexOfFace(int vert, int face) const
{
	for(int i=0; i<3; i++)
		if(vert == faces.coeff(i, face))
			return true;
	return false;
}

bool Mesh::neighboringFaces(int face1, int face2) const
{
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			if(faces.coeff(i, face1) == faces.coeff(j, face2))
				return true;
	return false;
}
