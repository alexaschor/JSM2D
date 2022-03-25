#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <random>
#include <chrono>
#include <thread>

#include "SETTINGS.h"
#include "field.h"
#include "julia.h"

using namespace std;

Real clamp(Real x, Real a, Real b) {
    return (x < a)? a : ((x > b)? b : x);
}

class PolySDF: public FieldFunction2D {
public:
    vector<VEC2F> v;

    PolySDF(vector<VEC2F> points): v(points) {};

    virtual Real getFieldValue(const VEC2F &p) const override {
        Real d = (p-v[0]).dot(p-v[0]);
        Real s = 1.0;
        for(int i=0, j=v.size()-1; i<v.size(); j=i, i++) {
            VEC2F e = v[j] - v[i];
            VEC2F w =  p - v[i];
            VEC2F b = w - e*clamp( w.dot(e)/e.dot(e), 0.0, 1.0 );
            d = min(d, b.dot(b));

            bool ba = p.y()>=v[i].y();
            bool bb = p.y()<v[j].y();
            bool bc = e.x()*w.y()>e.y()*w.x();

            if(ba==bb && bb==bc) s*=-1.0;
        }
        return s*sqrt(d);
    }

};


int main(int argc, char *argv[]) {

    if(argc != 4) {
        cout << "USAGE: " << endl;
        cout << "To generate a 2D SDF from a *.poly file" << endl;
        cout << " " << argv[0] << " <*.poly input> <SDF res> <*.f2d output>" << endl;
        exit(0);
    }

    FILE* file = fopen(argv[1], "r");
    if (file == NULL) {
        PRINT("Failed to read *.poly file: file open failed!");
        exit(0);
    }

    int numVerts;
    fscanf(file, "%d 2 0 0\n", &numVerts);

    vector<VEC2F> verts;

    for (int i = 0; i < numVerts; ++i) {
        Real x, y;
        fscanf(file, "%*d %lf %lf\n", &x, &y);
        verts.push_back(VEC2F(x, y));
    }

    PolySDF sdf(verts);

    int res = atoi(argv[2]);

    VirtualGrid2D vg(res, res, VEC2F(-1, -1), VEC2F(1,1), &sdf);
    vg.writeF2D(string(argv[3]), true);

    return 0;
}

