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

    if(argc != 11) {
        cout << "USAGE: " << endl;
        cout << "To create a shaped Julia set from a distance field:" << endl;
        cout << " " << argv[0] << " <distance field> <rotation field> <output resolution> <fill level> <offset> <min x> <min y> <max x> <max y> <output ppm>" << endl;
        exit(0);
    }

    AABB_2D bounds = AABB_2D(VEC2F(-1,-1), VEC2F(1,1));

    ArrayGrid2D distanceFieldCoarse((string(argv[1])));
    distanceFieldCoarse.setMapBox(bounds);
    InterpolationGrid2D distanceField(&distanceFieldCoarse);

    // ArrayGrid2D rotFieldCoarse((string(argv[2])));
    // rotFieldCoarse.setMapBox(bounds);
    // InterpolationGrid2D rotField(&rotFieldCoarse);
    // RotationMap2D rm(&rotField);

    POLYNOMIAL_2D poly = POLYNOMIAL_2D::randPolynomialInBox(AABB_2D(VEC2F(-0.5,-0.5), VEC2F(0.5,0.5)),
        8,
        12,
        15,
        true);

    int res = atoi(argv[3]);

    Real fillLevel = atof(argv[4]);
    Real offset = atof(argv[5]);

    VEC2F minView(atof(argv[6]), atof(argv[7]));
    VEC2F maxView(atof(argv[8]), atof(argv[9]));

    // DistanceGuidedMap m2(&distanceField, &rm, fillLevel, offset);
    DistanceGuidedMap m2(&distanceField, &poly, fillLevel, offset);
    JuliaSet julia(&m2, 30);

    julia.mode = JuliaSet::SET_MEMBERSHIP;
    // julia.mode = JuliaSet::LOG_MAGNITUDE;


    VirtualGrid2D vg(res, res, minView, maxView, &julia);
    // vg.writePPM(string(argv[10]), true);
    vg.writePPM(string(argv[10]), false);

    // VirtualGrid2D rfvg(res, res, minView, maxView, &rotField);
    // rfvg.writePPM(string(argv[10]), false);

    return 0;
}

