#include "AABBBroadPhase.h"
#include <set>
#include "History.h"
#include "Mesh.h"
#include <iostream>

using namespace std;

void AABBBroadPhase::findCollisionCandidates(const History &h, const Mesh &m, double outerEta, set<VertexFaceStencil> &vfs, set<EdgeEdgeStencil> &ees, const std::set<int> &fixedVerts)
{
	vfs.clear();
	ees.clear();
	AABBNode *tree = buildAABBTree(h, m, outerEta);
	intersect(tree, tree, m, vfs, ees, fixedVerts);
	delete tree;
}

AABBNode *AABBBroadPhase::buildAABBTree(const History &h, const Mesh &m, double outerEta)
{
	vector<AABBNode *> leaves;
	for(int i=0; i<(int)m.faces.cols(); i++)
	{
		AABBLeafNode *node = new AABBLeafNode;
		node->face = i;
		int verts[3];
		for(int j=0; j<3; j++)
		{
			verts[j] = m.faces.coeff(j, i);
			node->mincorner[j] = std::numeric_limits<double>::infinity();
			node->maxcorner[j] = -std::numeric_limits<double>::infinity();
		}
		for(int j=0; j<3; j++)
		{
			for(vector<HistoryEntry>::const_iterator it = h.getVertexHistory(verts[j]).begin(); it != h.getVertexHistory(verts[j]).end(); ++it)
			{
				for(int k=0; k<3; k++)
				{
					node->mincorner[k] = min(it->pos[k] - outerEta, node->mincorner[k]);
					node->maxcorner[k] = max(it->pos[k] + outerEta, node->maxcorner[k]);
				}
			}
		}
		leaves.push_back(node);		
	}
	return buildAABBInterior(leaves);
}

AABBNode *AABBBroadPhase::buildAABBInterior(vector<AABBNode *> &children)
{
	int nchildren = children.size();
	assert(nchildren > 0);
	if(nchildren == 1)
		return children[0];

	AABBInteriorNode *node = new AABBInteriorNode;
	
	for(int i=0; i<3; i++)	
	{
		node->mincorner[i] = std::numeric_limits<double>::infinity();
		node->maxcorner[i] = -std::numeric_limits<double>::infinity();
	}

	for(vector<AABBNode *>::iterator it = children.begin(); it != children.end(); ++it)
	{
		for(int j=0; j<3; j++)
		{
			node->mincorner[j] = min((*it)->mincorner[j], node->mincorner[j]);
			node->maxcorner[j] = max((*it)->maxcorner[j], node->maxcorner[j]);
		}
	}
	double lengths[3];
	for(int i=0; i<3; i++)
	{
		lengths[i] = node->maxcorner[i] - node->mincorner[i];
	}
	if(lengths[0] >= lengths[1] && lengths[0] >= lengths[2])
		node->splitaxis = 0;
	else if(lengths[1] >= lengths[0] && lengths[1] >= lengths[2])
		node->splitaxis = 1;
	else
		node->splitaxis = 2;

	sort(children.begin(), children.end(), NodeComparator(node->splitaxis));

	vector<AABBNode *> left;
	int child=0;
	for(; child<nchildren/2; child++)
		left.push_back(children[child]);
	node->left = buildAABBInterior(left);
	vector<AABBNode *> right;
	for(; child<nchildren; child++)
		right.push_back(children[child]);
	node->right = buildAABBInterior(right);
	return node;
}

void AABBBroadPhase::intersect(AABBNode *left, AABBNode *right, const Mesh &m, std::set<VertexFaceStencil> &vfs, std::set<EdgeEdgeStencil> &ees, const std::set<int> &fixedVerts)
{
	for(int axis=0; axis<3; axis++)
	{
		if(left->maxcorner[axis] < right->mincorner[axis] ||
			right->maxcorner[axis] < left->mincorner[axis])
			return;		
	}
	if(!left->isLeaf())
	{
		AABBInteriorNode *ileft = (AABBInteriorNode *)left;
		intersect(ileft->left, right, m, vfs, ees, fixedVerts);
		intersect(ileft->right, right, m, vfs, ees, fixedVerts);	
	}
	else if(!right->isLeaf())
	{
		AABBInteriorNode *iright = (AABBInteriorNode *)right;
		intersect(left, iright->left, m, vfs, ees, fixedVerts);
		intersect(left, iright->right, m, vfs, ees, fixedVerts);
	}
	else
	{
		AABBLeafNode *lleft = (AABBLeafNode *)left;
		AABBLeafNode *lright = (AABBLeafNode *)right;
		if(m.neighboringFaces(lleft->face, lright->face))
			return;

		// 6 vertex-face and 9 edge-edge
		for(int i=0; i<3; i++)
		{
			bool alllfixed = true;
			bool allrfixed = true;
			alllfixed = alllfixed && fixedVerts.count(m.faces.coeff(i, lleft->face));
			allrfixed = allrfixed && fixedVerts.count(m.faces.coeff(i, lright->face));
			for(int j=0; j<3; j++)
			{
				alllfixed = alllfixed && fixedVerts.count(m.faces.coeff(j, lright->face));
				allrfixed = allrfixed && fixedVerts.count(m.faces.coeff(j, lleft->face));
			}
			if(!alllfixed)
				vfs.insert(VertexFaceStencil(m.faces.coeff(i, lleft->face), m.faces.coeff(0, lright->face), m.faces.coeff(1, lright->face), m.faces.coeff(2, lright->face)));
			if(!allrfixed)
				vfs.insert(VertexFaceStencil(m.faces.coeff(i, lright->face), m.faces.coeff(0, lleft->face), m.faces.coeff(1, lleft->face), m.faces.coeff(2, lleft->face)));
			for(int j=0; j<3; j++)
			{
				bool allefixed = true;
				allefixed = allefixed && fixedVerts.count(m.faces.coeff(i, lleft->face));
				allefixed = allefixed && fixedVerts.count(m.faces.coeff((i+1)%3, lleft->face));
				allefixed = allefixed && fixedVerts.count(m.faces.coeff(j, lright->face));
				allefixed = allefixed && fixedVerts.count(m.faces.coeff((j+1)%3, lright->face));
				if(!allefixed)
					ees.insert(EdgeEdgeStencil(m.faces.coeff(i, lleft->face), m.faces.coeff((i+1)%3, lleft->face), m.faces.coeff(j, lright->face), m.faces.coeff((j+1)%3, lright->face)));
			}
		}		
	}
}
