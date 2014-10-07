#ifndef STENCILS_H
#define STENCILS_H

struct VertexFaceStencil
{
  int p;
  int q0,q1,q2;

  bool operator<(const VertexFaceStencil &other) const
  {
    if(p < other.p)
	return true;
    else if(p > other.p)
	return false;

    if(q0 < other.q0)
	return true;
    else if(q0 > other.q0)
	return false;

    if(q1 < other.q1)
	return true;
    else if(q1 > other.q1)
	return false;

    return q2 < other.q2;
  }
};

struct EdgeEdgeStencil
{
  int p0,p1;
  int q0,q1;

  bool operator<(const EdgeEdgeStencil &other) const
  {
    if(p0 < other.p0)
	return true;
    else if(p0 > other.p0)
	return false;
   
    if(p1 < other.p1)
	return true;
    else if(p1 > other.p1)
	return false;

    if(q0 < other.q0)
	return true;
    else if(q0 > other.q0)
	return false;

    return q1 < other.q1;
  }
};

#endif
