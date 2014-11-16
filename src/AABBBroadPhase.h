#include "RetrospectiveDetection.h"

struct AABBNode
{
	virtual ~AABBNode() {}
	virtual bool isLeaf() = 0;

	double mincorner[3];
	double maxcorner[3];
};

class NodeComparator
{
public:
	NodeComparator(int axis) : axis(axis) {}

	bool operator()(AABBNode *left, AABBNode *right)
	{
		return left->mincorner[axis] < right->mincorner[axis];
	}

private:
	int axis;
};

struct AABBInteriorNode : public AABBNode
{
	virtual ~AABBInteriorNode() 
	{
		delete left; 
		delete right;
	}

	virtual bool isLeaf() {return false;}

	int splitaxis;
	AABBNode *left, *right;
};

struct AABBLeafNode : public AABBNode
{
	virtual bool isLeaf() {return true;}

	int face;
};

class AABBBroadPhase : public BroadPhase
{
public:
	virtual void findCollisionCandidates(const History &h, const Mesh &m, double outerEta, std::set<VertexFaceStencil> &vfs, std::set<EdgeEdgeStencil> &ees, const std::set<int> &fixedVerts);
private:
	AABBNode *buildAABBTree(const History &h, const Mesh &m, double outerEta);
	AABBNode *buildAABBInterior(std::vector<AABBNode *> &children);
	void intersect(AABBNode *left, AABBNode *right, const Mesh &m, std::set<VertexFaceStencil> &vfs, std::set<EdgeEdgeStencil> &ees, const std::set<int> &fixedVerts);
};
