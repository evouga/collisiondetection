#include "CTCD.h"
#include <iostream>
#include <Eigen/Core>

using namespace std;
using namespace Eigen;

int main()
{
    {
        Vector3d q0(0,0,0);
        Vector3d p0(1,0,0);
        Vector3d q1(2,0,0);
        Vector3d p1(2,1,0);

        Vector3d q0end(0,0,0);
        Vector3d p0end(1,0,0);
        Vector3d q1end(2,1,0);
        Vector3d p1end(-1,-1,0);

        double eta = 1e-6;
        double t = 0;

        cout << CTCD::edgeEdgeCTCD(q0, p0, q1, p1, q0end, p0end, q1end, p1end, eta, t);
        cout << " " << t << endl;
    }

    {
        Vector3d q0(0.5, 0.5, 1);
        Vector3d q1(0,0,0);
        Vector3d q2(1,0,0);
        Vector3d q3(0,1,0);
        Vector3d q0end(0.5, 0.5, 1);
        Vector3d q1end(0,0,0);
        Vector3d q2end(1,0,3);
        Vector3d q3end(0,1,0);

        double eta = 1e-6;
        double t = 0;
        cout << CTCD::vertexFaceCTCD(q0, q1, q2, q3, q0end, q1end, q2end, q3end, eta, t);
        cout << " " << t << endl;
    }

    {
        Vector3d q0(0.5,0,1);
        Vector3d q1(0, 0, 0);
        Vector3d q2(1,0,0);
        Vector3d q0end(0.5,0,-3);
        Vector3d q1end(0,0,0);
        Vector3d q2end(1,0,0);

        double eta = 1e-6;
        double t = 0;
        cout << CTCD::vertexEdgeCTCD(q0, q1, q2, q0end, q1end, q2end, eta, t);
        cout << " " << t << endl;
    }

    {
        Vector3d q0(0,0,0);
        Vector3d q1(1,0,0);
        Vector3d q0end(6,0,0);
        Vector3d q1end(1,0,0);
        double eta = 1e-6;
        double t = 0;

        cout << CTCD::vertexVertexCTCD(q0, q1, q0end, q1end, eta, t);
        cout << " " << t << endl;
    }
}
