#include "SETTINGS.h"
#include <cmath>
#include <iostream>
#include <vector>

#define NANOSVG_IMPLEMENTATION
#include "nanosvg.h"


using namespace std;

// the actual triangle data
vector<VEC2F> nodes;
vector<int> indices;
vector<VEC3I> triangles;

VEC2F mins, maxs;

///////////////////////////////////////////////////////////////////////
// Get the bounding box of the nodes
///////////////////////////////////////////////////////////////////////
void getBoundingBox(VEC2F& localMins, VEC2F& maxs)
{
    mins = nodes[0];
    maxs = nodes[0];

    for (unsigned int x = 1; x < nodes.size(); x++)
    {
        if (nodes[x][0] < mins[0])
            mins[0] = nodes[x][0];
        if (nodes[x][1] < mins[1])
            mins[1] = nodes[x][1];
        if (nodes[x][0] > maxs[0])
            maxs[0] = nodes[x][0];
        if (nodes[x][1] > maxs[1])
            maxs[1] = nodes[x][1];
    }

    cout << " Bounding box mins: " << mins[0] << " " << mins[1] << endl;
    cout << " Bounding box maxs: " << maxs[0] << " " << maxs[1] << endl;
}

///////////////////////////////////////////////////////////////////////
// write the SVG path out the a POLY file for Triangle
///////////////////////////////////////////////////////////////////////
void writePoly(const char* filename)
{
    FILE* file = NULL;
    file = fopen(filename, "w");
    if (file == NULL)
    {
        cout << " Could not open file " << filename << "!!!" << endl;
        exit(0);
    }
    cout << " Writing out " << filename << endl;

    // output vertices
    fprintf(file, "%i 2 0 0\n", (int)nodes.size());
    for (unsigned int x = 0; x < nodes.size(); x++)
        fprintf(file, "%i %f %f\n", x + 1, (float)nodes[x][0], (float)nodes[x][1]);

    // output edges
    fprintf(file, "%i 0\n", (int)nodes.size());
    for (unsigned int x = 0; x < nodes.size() - 1; x++)
        fprintf(file, "%i %i %i\n", x + 1, x + 1, x + 2);
    fprintf(file, "%i %i 1\n", (int)nodes.size(), (int)nodes.size());

    // output holes
    fprintf(file, "0\n");

    fclose(file);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
    if (argc < 3)
    {
        cout << " USAGE: " << argv[0] << " <SVG filename> <POLY filename>" << endl;
        return 0;
    }

    NSVGimage* image;
    image = nsvgParseFromFile(argv[1], "px", 96);
    printf("size: %f x %f\n", image->width, image->height);

    NSVGshape* shape;
    NSVGpath* path;
    for (shape = image->shapes; shape != NULL; shape = shape->next) {
        for (path = shape->paths; path != NULL; path = path->next) {
            for (int i = 0; i < path->npts-1; i += 3) {
                float* p = &path->pts[i*2];
                cout << " point " << i << ": " << p[0] << " " << p[1] << endl;

                VEC2F node(p[0], p[1]);
                nodes.push_back(node);
            }
        }
    }

    // Delete
    nsvgDelete(image);
    getBoundingBox(mins, maxs);

    VEC2F center = (mins + maxs) * 0.5;
    VEC2F diff = maxs - mins;
    Real scale = (diff[0] > diff[1]) ? diff[0] : diff[1];
    for (unsigned int x = 0; x < nodes.size(); x++)
    {
        nodes[x] -= center;
        nodes[x] *= 1.0 / scale;

        // Triangle seems to flip things
        nodes[x][1] *= -1;
    }

    writePoly(argv[2]);
    return 1;
}
