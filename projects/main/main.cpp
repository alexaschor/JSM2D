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

int main(int argc, char *argv[]) {

    if(argc != 5) {
        cout << "USAGE: " << endl;
        cout << "To create a shaped Julia set from a distance field:" << endl;
        cout << " " << argv[0] << " <distance field> <output resolution> <fill level> <output ppm>" << endl;
        exit(0);
    }

    ArrayGrid2D distanceField((string(argv[1])));
    // distanceField.writePPM("sdf.ppm");

    // Shift it up so that (0,0) is inside the object
    distanceField.mapBox.min() += VEC2F(0, 0.1);
    distanceField.mapBox.max() += VEC2F(0, 0.1);

    POLYNOMIAL_2D map = POLYNOMIAL_2D::randPolynomialInBox(AABB_2D(VEC2F(-1,-1), VEC2F(1,1)), 1, 1, 30, false);

    Real fillLevel = atof(argv[3]);
    DistanceGuidedMap m2(&distanceField, &map, fillLevel);

    JuliaSet julia(&m2, 5);
    // julia.mode = JuliaSet::SET_MEMBERSHIP;
    julia.mode = JuliaSet::LOG_MAGNITUDE;

    int res = atoi(argv[2]);

    VirtualGrid2D vg(res, res, VEC2F(-1,-1), VEC2F(1,1), &julia);
    vg.writePPM(string(argv[4]));

    return 0;
}

