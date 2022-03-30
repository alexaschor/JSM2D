#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include "SETTINGS.h"
#include "field.h"
#include "julia.h"

#define GL_SILENCE_DEPRECATION

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/gl.h> // OpenGL itself.
#include <GL/glu.h> // GLU support library.
#include <GL/glut.h> // GLUT support library.
#endif

#include <iostream>
using namespace std;

typedef struct JuliaParams {
    POLYNOMIAL_2D poly;
    Real C = 40;
    Real B = 0;
    int maxIterations = 3;
    Real escape = 20;
} JuliaParams;

JuliaParams GLOBAL_params;

VEC2F findDistanceFieldCenter(const Grid2D& distField) {
    // Descending the gradient would be a ton faster, but
    // let's just exhaustively search for now

    VEC2I minPt;
    Real min = numeric_limits<Real>::max();

    for (uint x=0; x<distField.xRes; x++) {
        for (uint y=0; y<distField.yRes; y++) {
            Real val = distField.get(x, y);
            if (val < min){
                minPt = VEC2I(x,y);
                min = val;
            }
        }
    }

    return distField.getCellCenter(minPt);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {
    if(argc < 7 || argc > 8) {
        cout << "USAGE: " << endl;
        cout << "To render an animation script file generated using the interactive viewer" << endl;
        cout << " " << argv[0] << " <params file> <anim file> <output res> <speed float value> <primary distance field> <output path> <swap distance field (optional)>" << endl;
        exit(0);
    }

    // Read params
    FILE * params;
    params = fopen(argv[1], "r");

    // FIXME all this will break if we later change Real to be float,
    // not sure why we'd do that but watch out just in case

    fscanf(params, "%lf\n", &GLOBAL_params.C);
    fscanf(params, "%lf\n", &GLOBAL_params.B);
    fscanf(params, "%d\n", &GLOBAL_params.maxIterations);
    fscanf(params, "%lf\n", &GLOBAL_params.escape);

    // Read roots from params file
    int numRoots;
    vector<COMPLEX> roots;
    vector<Real> powers;
    fscanf(params, "%d\n", &numRoots);

    for (int i = 0; i < numRoots; ++i) {
        Real real, imag, power;
        fscanf(params, "%lf, %lf, %lf\n", &real, &imag, &power);
        roots.push_back(COMPLEX(real, imag));
        powers.push_back(power);
    }

    GLOBAL_params.poly = POLYNOMIAL_2D(roots, powers);

    // Read frame bounds from script file
    FILE * script;
    script = fopen(argv[2], "r");

    vector<AABB_2D> keypoints;
    Real minX, minY, maxX, maxY;
    for ( ; fscanf(script, "%lf, %lf, %lf, %lf", &minX, &minY, &maxX, &maxY) != EOF; ) {
        keypoints.push_back(AABB_2D(VEC2F(minX, minY), VEC2F(maxX, maxY)));
    }

    // Interpolate given max speed (% of view width zoom per frame)
    Real maxSpeed = atof(argv[4]);
    vector<AABB_2D> frames;
    for (int i = 0; i < keypoints.size()-1; ++i) {
        AABB_2D f1 = keypoints[i];
        AABB_2D f2 = keypoints[i+1];

        Real dist = (f2.max() - f1.max()).norm() + (f2.min() - f1.min()).norm();
        Real screenDist = dist / (f1.span()[0]);

        // And subframe to limit that max speed
        int subSteps = (screenDist / maxSpeed);
        for (int j = 0; j < subSteps; ++j) {
            VEC2F min = f1.min() + (f2.min() - f1.min()) * ((Real) j / subSteps);
            VEC2F max = f1.max() + (f2.max() - f1.max()) * ((Real) j / subSteps);
            frames.push_back(AABB_2D(min, max));
        }

    }

    frames.push_back(keypoints[keypoints.size() - 1]);

    PRINTF("Got %lu keyframes, generated %lu frames", keypoints.size(), frames.size());

    // Read distance field
    ArrayGrid2D *distFieldCoarse = new ArrayGrid2D((string(argv[5])));
    distFieldCoarse->setMapBox(AABB_2D(VEC2F(-0.5,-0.5), VEC2F(0.5,0.5))); // Set to size 1, centered at origin
    InterpolationGrid2D* distField = new InterpolationGrid2D(distFieldCoarse, InterpolationGrid2D::LINEAR); // Interpolate

    // Shift the distance field so that the minimum distance is at the origin
    VEC2F dfCenter = findDistanceFieldCenter(*distField);
    distField->mapBox.min() -= dfCenter;
    distField->mapBox.max() -= dfCenter;

    // Compute julia set
    DistanceGuidedMap dMap(distField, &GLOBAL_params.poly, GLOBAL_params.C, GLOBAL_params.B);
    JuliaSet julia(&dMap, GLOBAL_params.maxIterations, GLOBAL_params.escape);
    julia.mode = JuliaSet::SET_MEMBERSHIP;

    int res = atoi(argv[3]);

    PB_START("Rendering zoom-in");
    for (int i = 0; i < frames.size(); ++i) {
        char buf[256];
        sprintf(buf, "%s/out_%04d.ppm", argv[6], i);

        VirtualGrid2D viewingField(res, res, frames[i].min(), frames[i].max(), &julia);
        viewingField.writePPM(string(buf));
        PB_PROGRESS((Real) i / frames.size());
    }
    PB_END();


    return 0;
}
