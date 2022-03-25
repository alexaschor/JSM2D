#ifndef FIELD_H
#define FIELD_H

#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <queue>

#include "SETTINGS.h"

using namespace std;

class AABB_2D: public AlignedBox<Real, 2> {
public:
    using AlignedBox<Real, 2>::AlignedBox;
    AABB_2D(VEC2F c1, VEC2F c2): AlignedBox<Real, 2>(c1.cwiseMin(c2), c1.cwiseMax(c2)) {}
    AABB_2D(): AABB_2D(VEC2F(0,0), VEC2F(0,0)){}

    VEC2F span() const {
        return max()-min();
    }

    VEC2F clamp(VEC2F pos) const {
        return pos.cwiseMax(min()).cwiseMin(max());
    }

    void setCenter(VEC2F newCenter) {
        VEC2F offset = newCenter - this->center();
        max() += offset;
        min() += offset;
    }
};

class FieldFunction2D {
private:
    Real (*fieldFunction)(VEC2F pos);
public:
    FieldFunction2D(Real (*fieldFunction)(VEC2F pos)):fieldFunction(fieldFunction) {}

    FieldFunction2D():fieldFunction(nullptr) {} // Only to be used by subclasses

    virtual Real getFieldValue(const VEC2F& pos) const {
        return fieldFunction(pos);
    }

    virtual Real operator()(const VEC2F& pos) const {
        return getFieldValue(pos);
    }
};

class Grid2D: public FieldFunction2D {
public:
    uint xRes, yRes;
    bool supportsNonIntegerIndices = false;

    // Bounds (aka center + lengths) for mapping into this grid
    // as a field function
    AABB_2D mapBox;
    bool hasMapBox = false;

    virtual uint totalCells() const {
        return xRes * yRes;
    }

    virtual Real get(uint x, uint y) const = 0;

    virtual Real operator[](size_t x) { // Access as C-style array;
        uint Cx = (x % yRes);
        uint Cy = (x - Cx)/yRes;
        return get(Cx, Cy);
    }

    virtual Real getf(Real x, Real y) const {
        (void) x; (void) y; // Suppress unused argument warning
        printf("This grid doesn't support non-integer indices!\n");
        exit(1);
    }

    virtual Real get(VEC2I pos) const {
        return get(pos[0], pos[1]);
    }

    virtual Real getf(VEC2F pos) const {
        return getf(pos[0], pos[1]);
    }

    virtual void setMapBox(AABB_2D box) {
        mapBox = box;
        hasMapBox = true;
    }

    virtual Real getFieldValue(const VEC2F& pos) const override {
        if (!hasMapBox) {
            printf("Attempting getFieldValue on a Grid without a mapBox!\n");
            exit(1);
        }

        VEC2F samplePoint = (pos - mapBox.min()).cwiseQuotient(mapBox.span());
        samplePoint = samplePoint.cwiseMax(VEC2F(0,0)).cwiseMin(VEC2F(1,1));

        const VEC2F indices = samplePoint.cwiseProduct(VEC2F(xRes-1, yRes-1));

        if (supportsNonIntegerIndices) {
            return getf(indices);
        } else {
            return get(indices.cast<int>());
        }
    }

    virtual void writeCSV(string filename) {
        ofstream out;
        out.open(filename);
        if (out.is_open() == false)
            return;

        for (uint i = 0; i < xRes; ++i) {
            for (uint j = 0; j < yRes; ++j) {
                out << i << ", " << j << ", " << get(i, j) << endl;
            }
        }

        printf("Wrote %d x %d field (%d values) to %s\n", xRes, yRes, (xRes * yRes), filename.c_str());

        out.close();
    }

    void writePPM(const string& filename, bool verbose=false) {
        int totalCells = xRes * yRes;

        Real min = numeric_limits<Real>::max();
        Real max = -min;

        Real* values = new Real[3 * totalCells];
        PB_DECL();
        if (verbose) PB_STARTD("Fetching values for PPM write");
#pragma omp parallel for
        for(uint y=0; y<yRes; ++y) {
            for(uint x=0; x<xRes; ++x) {
                int idx = ((y*xRes) + x) * 3;

                Real val = get(x,y);

                values[idx] = val;
                values[idx + 1] = val;
                values[idx + 2] = val;

                if (val < min) min = val;
                if (val > max) max = val;
            }
            if (verbose && y%30==0) PB_PROGRESS((float)y / yRes);
        }
        if (verbose) PB_END();

        unsigned char* pixels = new unsigned char[3 * totalCells];
#pragma omp parallel for
        for(uint y=0; y<yRes; ++y) {
            for(uint x=0; x<xRes; ++x) {
                int idx = ((y*xRes) + x) * 3;
                Real val = 255 * ((values[idx] - min) / (max - min));

                pixels[idx] = val;
                pixels[idx + 1] = val;
                pixels[idx + 2] = val;
            }
        }

        FILE *fp;
        fp = fopen(filename.c_str(), "wb");
        if (fp == NULL) {
            PRINTF("Could not open file '%s' for PPM write!", filename.c_str());
            exit(1);
        }

        fprintf(fp, "P6\n%d %d\n255\n", xRes, yRes);
        fwrite(pixels, 1, totalCells * 3, fp);
        fclose(fp);
        delete[] pixels;

        printf("Wrote %d x %d field (%d values) to %s\n", xRes, yRes, (xRes * yRes), filename.c_str());
    }

    // Writes to F2D file using the field bounds if the field has them, otherwise using
    // the resolution ofthe grid.
    void writeF2D(string filename, bool verbose = false) const {
        if (hasMapBox) {
            writeF2D(filename, mapBox, verbose);
        } else {
            writeF2D(filename, AABB_2D(VEC2F(0,0), VEC2F(xRes, yRes)), verbose);
        }
    }

    void writeF2D(string filename, AABB_2D bounds, bool verbose = false) const {
        FILE* file = fopen(filename.c_str(), "wb");

        if (file == NULL) {
            PRINT("Failed to write F2D: file open failed!");
            exit(0);
        }

        PB_DECL();
        if (verbose) {
            PB_STARTD("Writing %dx%d field to %s", xRes, yRes, filename.c_str());
        }

        // write dimensions
        fwrite((void*)&xRes, sizeof(int), 1, file);
        fwrite((void*)&yRes, sizeof(int), 1, file);

        MyEigen::write_vec2f(file, bounds.center());
        MyEigen::write_vec2f(file, bounds.span());

        const int totalCells = xRes*yRes;

        if (totalCells <= 0)
            return;

        // write data
        for (uint i = 0; i < yRes; ++i) {
            for (uint j = 0; j < xRes; ++j) {
                const Real val = get(j, i);
                double out;
                if (sizeof(Real) != sizeof(double)) {
                    out = (double) val;
                } else {
                    out = val;
                }

                fwrite((void*) (&out), sizeof(double), 1, file);
            }

            if (verbose && i % 10 == 0) {
                PB_PROGRESS((Real) i / xRes);
            }
        }

        if (verbose) {
            PB_END();
        }
    }
};

class ArrayGrid2D: public Grid2D {
private:
    Real* values;
public:

    // Create empty (not zeroed) field with given resolution
    ArrayGrid2D(uint xRes, uint yRes) {
        this->xRes = xRes;
        this->yRes = yRes;
        values = new Real[xRes * yRes];
    }

    // Create empty (not zeroed) field with given resolution
    ArrayGrid2D(VEC3I resolution): ArrayGrid2D(resolution[0], resolution[1]) {}

    // Create deep copy of another ArrayGrid2D
    ArrayGrid2D(const ArrayGrid2D& other) {
        this->xRes = other.xRes;
        this->yRes = other.yRes;

        values = new Real[xRes * yRes];
        memcpy(values, other.values, xRes * yRes * sizeof(Real));
    }

    // Read ArrayGrid2D from F2D
    ArrayGrid2D(string filename, string format = "f2d", bool verbose = false) {

        if (format == "f2d") {
            FILE* file = fopen(filename.c_str(), "rb");
            if (file == NULL) {
                PRINT("Failed to read F2D: file open failed!");
                exit(0);
            }

            int xRes, yRes;
            VEC2F center, lengths;

            // read dimensions
            fread((void*)&xRes, sizeof(int), 1, file);
            fread((void*)&yRes, sizeof(int), 1, file);

            MyEigen::read_vec2f(file, center);
            MyEigen::read_vec2f(file, lengths);

            this->xRes = xRes;
            this->yRes = yRes;

            try {
                values = new Real[xRes * yRes];
            } catch(bad_alloc& exc) {
                printf("Failed to allocate %.2f MB for ArrayGrid2D read from file!\n", (xRes * yRes * sizeof(Real)) / pow(2.0,20.0));
                exit(0);
            }

            if (verbose) {
                printf("Reading %d x %d field from %s... ", xRes, yRes, filename.c_str());
                fflush(stdout);
            }

            AABB_2D bounds((center - lengths/2), (center + lengths/2));
            setMapBox(bounds);

            const int totalCells = xRes * yRes;
            // always read in as a double
            if (sizeof(Real) != sizeof(double)) {
                double* dataDouble = new double[totalCells];
                fread((void*)dataDouble, sizeof(double), totalCells, file);

                for (int x = 0; x < totalCells; x++)
                    values[x] = dataDouble[x];

                delete[] dataDouble;
            } else fread((void*)values, sizeof(Real), totalCells, file);

        } else {
            PRINT("CSV import not implemented yet!");
            exit(1);
        }

    }
    // Destructor
    virtual ~ArrayGrid2D() {
        delete values;
    }

    // Access value based on integer indices
    Real get(uint x, uint y) const override {
        return values[y * xRes + x];
    }

    // Access value directly (allows setting)
    Real& at(uint x, uint y) {
        return values[y * xRes + x];
    }

    // Access value directly in C-style array (allows setting)
    Real& at(size_t x) {
        return values[x];
    }

    Real& atFieldPos(VEC2F pos) {
        if (!hasMapBox) {
            printf("Attempting atFieldPos on an ArrayGrid without a mapBox!\n");
            exit(1);
        }

        VEC2F samplePoint = (pos - mapBox.min()).cwiseQuotient(mapBox.span());
        samplePoint = samplePoint.cwiseMax(VEC2F(0,0)).cwiseMin(VEC2F(1,1));

        const VEC2F indices = samplePoint.cwiseProduct(VEC2F(xRes-1, yRes-1));

        return at(indices[0], indices[1]);
    }

    Real& atFieldPos(Real x, Real y) {
        return atFieldPos(VEC2F(x,y));
    }


    Real& operator()(VEC2F pos) {
        return atFieldPos(pos);
    }

    void operator+=(Real x){
        for (unsigned i = 0; i < (xRes*yRes); ++i) {
            values[i] += x;
        }
    }

    void operator-=(Real x){
        for (unsigned i = 0; i < (xRes*yRes); ++i) {
            values[i] -= x;
        }
    }

    void operator*=(Real x){
        for (unsigned i = 0; i < (xRes*yRes); ++i) {
            values[i] *= x;
        }
    }

    void operator/=(Real x){
        for (unsigned i = 0; i < (xRes*yRes); ++i) {
            values[i] /= x;
        }
    }

    void abs(){
        for (unsigned i = 0; i < (xRes*yRes); ++i) {
            values[i] = ::abs(values[i]);
        }
    }

    void log(Real x){
        for (unsigned i = 0; i < (xRes*yRes); ++i) {
            values[i] = ::log(values[i]) / ::log(x);
        }
    }

    // Create field from scalar function by sampling it on a regular grid
    ArrayGrid2D(uint xRes, uint yRes, VEC2F functionMin, VEC2F functionMax, FieldFunction2D *function):ArrayGrid2D(xRes, yRes){

        VEC2F gridResF(xRes, yRes);

#pragma omp parallel for
        for (uint i = 0; i < xRes; i++) {
            for (uint j = 0; j < yRes; j++) {
                VEC2F gridPointF(i, j);
                VEC2F fieldDelta = functionMax - functionMin;

                VEC2F samplePoint = functionMin + (gridPointF.cwiseQuotient(gridResF - VEC2F(1,1)).cwiseProduct(fieldDelta));

                this->at(i, j) = function->getFieldValue(samplePoint);
            }
        }

    }


};

class VirtualGrid2D: public Grid2D {
private:
    FieldFunction2D *fieldFunction;
    VEC2F functionMin, functionMax;

    VEC2F getSamplePoint(Real x, Real y) const {
        VEC2F gridPointF(x, y);
        VEC2F gridResF(xRes, yRes);
        VEC2F fieldDelta = functionMax - functionMin;
        VEC2F samplePoint = functionMin + (gridPointF.cwiseQuotient(gridResF - VEC2F(1,1)).cwiseProduct(fieldDelta));

        return samplePoint;
    }
public:

    VirtualGrid2D(uint xRes, uint yRes, VEC2F functionMin, VEC2F functionMax, FieldFunction2D *fieldFunction):
        fieldFunction(fieldFunction),
        functionMin(functionMin),
        functionMax(functionMax) {
            this->xRes = xRes;
            this->yRes = yRes;

            this->setMapBox(AABB_2D(functionMin, functionMax));

            this->supportsNonIntegerIndices = true;
        }

    // Virtually (shallow) downsample an existing grid
    // VirtualGrid2D(uint xRes, uint yRes, uint zRes, Grid2D *other) {
    // TODO implement
    // }

    virtual Real get(uint x, uint y) const override {
        return getf(x, y);
    }

    virtual Real getf(Real x, Real y) const override {
        return fieldFunction->getFieldValue(getSamplePoint(x, y));
    }
};

// Hash function for Eigen matrix and vector.
// From https://wjngkoh.wordpress.com/2015/03/04/c-hash-function-for-eigen-matrix-and-vector/
template<typename T> struct matrix_hash : unary_function<T, size_t> {
    size_t operator()(T const& matrix) const {
        // Note that it is oblivious to the storage order of Eigen matrix (column- or
        // row-major). It will give you the same hash value for two different matrices if they
        // are the transpose of each other in different storage order.
        size_t seed = 0;
        for (long i = 0; i < matrix.size(); ++i) { // For some reason Eigen size is not unsigned?!
            auto elem = *(matrix.data() + i);
            seed ^= hash<typename T::Scalar>()(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};


class VirtualGrid2DCached: public VirtualGrid2D {
protected:
    mutable unordered_map<VEC2F, Real, matrix_hash<VEC2F>> map;

public:
    mutable int numQueries = 0;
    mutable int numHits = 0;
    mutable int numMisses = 0;

    using VirtualGrid2D::VirtualGrid2D;

    virtual Real get(uint x, uint y) const override {
        return getf(x,y);
    }

    virtual Real getf(Real x, Real y) const override {
        VEC2F key(x,y);
        numQueries++;

        auto search = map.find(key);
        if (search != map.end()) {
            numHits++;
            return search->second;
        }

        Real result = VirtualGrid2D::get(x,y);
        map[key] = result;
        numMisses++;
        return result;
    }

};

class VirtualGrid2DLimitedCache: public VirtualGrid2DCached {
private:
    size_t maxSize;
    mutable queue<VEC2F> cacheQueue;

public:
    // Instantiates a VirtualGrid2D with a limited-size cache. When additional
    // items are inserted into the cache (beyond the capacity), the cache will
    // forget the item that was least recently inserted. If capacity -1 is
    // specified (default), it defaults to a size equal to three Y slices
    // through the field.
    VirtualGrid2DLimitedCache(uint xRes, uint yRes, VEC2F functionMin, VEC2F functionMax,  FieldFunction2D *fieldFunction, int capacity = -1):
        VirtualGrid2DCached(xRes, yRes, functionMin, functionMax, fieldFunction) {
            PRINTV2(functionMin);
            PRINTV2(functionMax);
            if (capacity == -1) {
                maxSize = xRes * 3;
            } else {
                maxSize = capacity;
            }
        }

    virtual Real get(uint x, uint y) const override {
        return getf(x,y);
    }

    virtual Real getf(Real x, Real y) const override {
        VEC2F key(x,y);
        numQueries++;

        auto search = map.find(key);
        if (search != map.end()) {
            numHits++;
            return search->second;
        }

        Real result = VirtualGrid2D::getf(x,y);
        map[key] = result;

        if (cacheQueue.size() >= maxSize) {
            map.erase(cacheQueue.front());
            cacheQueue.pop();
        }

        cacheQueue.push(key);

        numMisses++;
        return result;
    }
};


class InterpolationGrid2D: public Grid2D {
private:
    Real interpolate(Real x0, Real x1, Real d) const {
        switch (mode) {
        case LINEAR:
            return ((1 - d) * x0) + (d * x1);
        case SMOOTHSTEP:
            d = (3 * d * d) - (2 * d * d * d);
            return ((1 - d) * x0) + (d * x1);
        }
        assert(false);
        return -1;
    }


public:
    Grid2D* baseGrid;

    enum INTERPOLATION_MODE {
        LINEAR,
        SMOOTHSTEP
    };

    INTERPOLATION_MODE mode;

    InterpolationGrid2D(Grid2D* baseGrid, INTERPOLATION_MODE mode = LINEAR) {
        this->baseGrid = baseGrid;
        xRes = baseGrid->xRes;
        yRes = baseGrid->yRes;
        this->mode = mode;
        this->supportsNonIntegerIndices = true;

        if (baseGrid->hasMapBox) this->setMapBox(baseGrid->mapBox);

        if (baseGrid->supportsNonIntegerIndices) {
            PRINT("Warning: laying an InterpolationGrid2D over a grid which already supports non-integer indices!");
        }
    }

    virtual Real get(uint x, uint y) const override {
        return baseGrid->get(x, y);
    }

    virtual Real getf(Real x, Real y) const override {
        // "Bilinear" interpolation with whatever technique we select

        uint x0 = floor(x);
        uint y0 = floor(y);

        uint x1 = x0 + 1;
        uint y1 = y0 + 1;

        // Clamp if out of bounds
        x0 = (x0 > xRes - 1) ? xRes - 1 : x0;
        y0 = (y0 > yRes - 1) ? yRes - 1 : y0;

        x1 = (x1 > xRes - 1) ? xRes - 1 : x1;
        y1 = (y1 > yRes - 1) ? yRes - 1 : y1;


        const Real xd = min(1.0, max(0.0, (x - x0) / ((Real) x1 - x0)));
        const Real yd = min(1.0, max(0.0, (y - y0) / ((Real) y1 - y0)));

        // First grab 2D surroundings...
        const Real c00 = baseGrid->get(x0, y0);
        const Real c01 = baseGrid->get(x0, y1);
        const Real c10 = baseGrid->get(x1, y0);
        const Real c11 = baseGrid->get(x1, y1);

        // Extract 1D interpolated line...
        const Real c0 = interpolate(c00, c10, xd);
        const Real c1 = interpolate(c01, c11, xd);

        // And grab 0D point.
        const Real output = interpolate(c0, c1, yd);

        return output;
    }



};



#endif
