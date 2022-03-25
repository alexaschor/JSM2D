#ifndef JULIA_H
#define JULIA_H

#include "field.h"
#include <complex>
#include <limits>
#include <random>

typedef std::complex<Real> COMPLEX;

inline COMPLEX normalizedComplex(COMPLEX q){ return (q / (sqrt(real(q)*real(q) + imag(q)*imag(q)))); }
inline Real complexMagnitude(COMPLEX q){ return sqrt(real(q)*real(q) + imag(q)*imag(q)); }

class ComplexMap {
public:
    virtual COMPLEX getFieldValue(COMPLEX q) const = 0;

    virtual COMPLEX operator()(COMPLEX q) const {
        return getFieldValue(q);
    }
};

class SimpleJuliaQuat: public ComplexMap {
public:
    COMPLEX c;
    SimpleJuliaQuat(COMPLEX c): c(c) {}

    virtual COMPLEX getFieldValue(COMPLEX q) const override {
        return (q * q) + c;
    }
};

class POLYNOMIAL_2D: public ComplexMap {
public:
    vector<COMPLEX> roots;
    vector<Real> powers;

    static POLYNOMIAL_2D randPolynomialInBox(AABB_2D box, Real minPower, Real maxPower, uint nRoots, bool allowNonIntegerPowers = true) {
        random_device rd;
        default_random_engine eng(rd());
        // default_random_engine eng(12321); // XXX For consistent results

        Real minX = box.min()[0];
        Real minY = box.min()[1];

        Real maxX = box.max()[0];
        Real maxY = box.max()[1];

        uniform_real_distribution<> spaceXD(minX, maxX);
        uniform_real_distribution<> spaceYD(minY, maxY);

        uniform_real_distribution<> powerD(minPower, maxPower);

        vector<COMPLEX> roots{};
        vector<Real> powers{};
        for (uint i = 0; i < nRoots; ++i) {
            roots.push_back(COMPLEX(spaceXD(eng), spaceYD(eng)));
            powers.push_back( allowNonIntegerPowers? powerD(eng) : int(powerD(eng)) );
        }

        POLYNOMIAL_2D poly(roots, powers);
        return poly;
    }

    bool operator==(const POLYNOMIAL_2D& other) {
        return (roots == other.roots && powers == other.powers); //XXX won't work for out of order duplicates
    }



    POLYNOMIAL_2D(){}

    POLYNOMIAL_2D(vector<COMPLEX> roots, vector<Real> powers): roots(roots), powers(powers) {}

    COMPLEX evaluate(COMPLEX q) const {
        COMPLEX sum = 0;
        for (unsigned i=0; i < roots.size(); ++i) {
            sum += pow(q - roots[i], powers[i]);
        }
        return sum;
    }

    virtual COMPLEX getFieldValue(COMPLEX q) const override {
        return evaluate(q);
    }
};

class DistanceGuidedMap: public ComplexMap {
public:
    Grid2D* distanceField;
    ComplexMap* p;

    Real c;

public:
    DistanceGuidedMap(Grid2D* distanceField,  ComplexMap* p, Real c = 300):
        distanceField(distanceField), p(p), c(c){}

    virtual COMPLEX getFieldValue(COMPLEX q) const override {
        VEC2F qV2(real(q), imag(q));

        const Real distance = (*distanceField)(qV2);
        Real radius = exp(c * distance);

        // Evaluate polynomial
        q = p->getFieldValue(q);
        q = normalizedComplex(q) * radius;

        return q;
    }

};

class JuliaSet: public FieldFunction2D {
public:
    ComplexMap* map;
    int maxIterations;
    Real escape;

    typedef enum {
        LOG_MAGNITUDE,
        SET_MEMBERSHIP
    } JuliaOutputMode;

    JuliaOutputMode mode = LOG_MAGNITUDE;

    JuliaSet(ComplexMap* map, int maxIterations=3, Real escape=20): map(map), maxIterations(maxIterations), escape(escape) {}

    virtual Real getFieldValue(const VEC2F& pos) const override {
        COMPLEX iterate(pos[0], pos[1]);

        iterate = map->getFieldValue(iterate);
        Real magnitude = complexMagnitude(iterate);

        int totalIterations = 1;

        while (magnitude < escape && totalIterations < maxIterations) {
            // Evaluate polynomial
            iterate = map->getFieldValue(iterate);
            magnitude = complexMagnitude(iterate);

            totalIterations++;
        }

        if (mode == LOG_MAGNITUDE) {
            return log(magnitude);
        }
        return (magnitude < escape)? 1 : 0;
    }
};


// =============== INSPECTION TOOLS =======================

class ComplexMapRotField: public FieldFunction2D {
public:
    ComplexMap *func;
    int i;

    ComplexMapRotField(ComplexMap* func, int i): func(func), i(i) {}

    virtual Real getFieldValue(const VEC2F& pos) const override {
        COMPLEX in(pos[0], pos[1]);
        COMPLEX out = normalizedComplex(func->getFieldValue(in));

        return (i==1)? real(out) : imag(out);
    }
};

class DistanceMapInspection {
public:

    static void sampleAttractingPercent(FieldFunction2D& ff, int gridRes, Real min, Real max, Real radDelta) {

        PB_START("Sampling radii");

        if (min == 0) min = radDelta;

        for (Real rad = min; rad < max; rad += radDelta) {
            uint numPts = 0, numAttracting = 0;

            Real gridDelta = (2 * rad)/gridRes;
            for (Real x = -rad; x <= rad; x+=gridDelta) {
                for (Real y = -rad; y <= rad; y+=gridDelta) {
                    Real ptRad = sqrt(x*x + y*y);
                    if (ptRad <= rad) {
                        if (ff(VEC2F(x,y)) > rad) {
                            numAttracting++;
                        }
                        numPts++;
                    }
                }
            }

            printf("%f, %f\n", rad, (Real) numAttracting / numPts);
            PB_PROGRESS((rad - min) / (max - min));

        }

        PB_END();

    }

    static bool isAttractingRegion(FieldFunction2D& ff, Real radius, Real delta) {
        // Do it the stupid way
        for (Real x = -radius; x <= radius; x+=delta) {
            for (Real y = -radius; y <= radius; y+=delta) {
                if (sqrt(x*x + y*y) <= radius) {
                    if (ff(VEC2F(x,y)) > radius) {
                        // PRINTF("Found value %.5f inside radius %.5f", ff(VEC2F(x,y)), radius);
                        return false;
                    }
                }
            }
        }
        return true;
    }

    static Real getMaxAttractingRadius(FieldFunction2D& ff, Real gridDelta) {
        Real radMin = 0;
        Real radDelta = 1e-1;
        int maxIterations = 10;

        Real maxFound = 0;

        for (int i = 0; i < maxIterations; ++i) {
            Real radius = radMin + radDelta;
            while (isAttractingRegion(ff, radius, gridDelta)) {
                // PRINTD(radius);
                if (radius > maxFound) maxFound = radius;
                radius += radDelta;
            }
            radMin = radius - radDelta;
            radDelta *= 1e-1;
        }

        return maxFound;
    }

};


#endif
