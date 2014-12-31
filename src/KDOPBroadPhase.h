#include "RetrospectiveDetection.h"

#define K 13

struct KDOPNode
{
  virtual ~KDOPNode() {}
  virtual bool isLeaf() = 0;

  double mins[K];
  double maxs[K];
};

class NodeComparator
{
public:
	NodeComparator(int axis) : axis(axis) {}

	bool operator()(KDOPNode *left, KDOPNode *right)
	{
		return left->mins[axis] < right->mins[axis];
	}

private:
	int axis;
};

struct KDOPInteriorNode : public KDOPNode
{
	virtual ~KDOPInteriorNode() 
	{
		delete left; 
		delete right;
	}

	virtual bool isLeaf() {return false;}

	int splitaxis;
	KDOPNode *left, *right;
};

struct KDOPLeafNode : public KDOPNode
{
	virtual bool isLeaf() {return true;}

	int face;
};

class KDOPBroadPhase : public BroadPhase
{
public:
  KDOPBroadPhase();

  virtual void findCollisionCandidates(const History &h, const Mesh &m, double outerEta, std::set<VertexFaceStencil> &vfs, std::set<EdgeEdgeStencil> &ees, const std::set<int> &fixedVerts);
 private:
  KDOPNode *buildKDOPTree(const History &h, const Mesh &m, double outerEta);
  KDOPNode *buildKDOPInterior(std::vector<KDOPNode *> &children);
  void intersect(KDOPNode *left, KDOPNode *right, const Mesh &m, std::set<VertexFaceStencil> &vfs, std::set<EdgeEdgeStencil> &ees, const std::set<int> &fixedVerts);
  
  std::vector<Eigen::Vector3d> DOPaxis;
};
