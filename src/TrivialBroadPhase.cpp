#include "TrivialBroadPhase.h"
#include "Mesh.h"

void TrivialBroadPhase::findCollisionCandidates(const History &h, const Mesh &m, double outerEta, std::set<VertexFaceStencil> &vfs, std::set<EdgeEdgeStencil> &ees, const std::set<int> &fixedVerts)
{
	int nverts = m.vertices.size()/3;
	int nfaces = m.faces.cols();

	for(int i=0; i<nverts; i++)
	{
		for(int j=0; j<nfaces; j++)
		{
			if(m.vertexOfFace(i, j))
				continue;

			bool allfixed = true;
			allfixed = allfixed && fixedVerts.count(i);
			for(int k=0; k<3; k++)
				allfixed = allfixed && fixedVerts.count(m.faces.coeff(k,j));
			if(allfixed)
				continue;

			VertexFaceStencil newvfs(i, m.faces.coeff(0, j), m.faces.coeff(1, j), m.faces.coeff(2, j));
			vfs.insert(newvfs);
		}
	}

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

			bool allfixed = true;
			allfixed = allfixed && fixedVerts.count(e1v1);
			allfixed = allfixed && fixedVerts.count(e1v2);
			allfixed = allfixed && fixedVerts.count(e2v1);
			allfixed = allfixed && fixedVerts.count(e2v2);
			if(allfixed)
				continue;
			
			EdgeEdgeStencil newees(e1v1, e1v2, e2v1, e2v2);
			ees.insert(newees);
		}
	}
}
