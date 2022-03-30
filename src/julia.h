#ifndef JULIA_H
#define JULIA_H

#include "field.h"
#include <cmath>
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
    Real b;

public:
    DistanceGuidedMap(Grid2D* distanceField,  ComplexMap* p, Real c = 300, Real b=0):
        distanceField(distanceField), p(p), c(c), b(b){}

    virtual COMPLEX getFieldValue(COMPLEX q) const override {
        VEC2F qV2(real(q), imag(q));

        const Real distance = (*distanceField)(qV2);
        Real radius = exp((c * distance) + b);

        // Evaluate polynomial
        q = p->getFieldValue(q);
        q = normalizedComplex(q) * radius;

        return q;
    }

};

class JuliaSet: public FieldFunction2D {
public:
    ComplexMap* map;
    ComplexMap* secondaryMap = 0;
    int maxIterations;
    Real escape;

    typedef enum {
        LOG_MAGNITUDE,
        SET_MEMBERSHIP
    } JuliaOutputMode;

    JuliaOutputMode mode = LOG_MAGNITUDE;

    JuliaSet(ComplexMap* map, int maxIterations=3, Real escape=20): map(map), secondaryMap(0), maxIterations(maxIterations), escape(escape) {}
    JuliaSet(ComplexMap* map1, ComplexMap* map2, int maxIterations=3, Real escape=20): map(map1), secondaryMap(map2), maxIterations(maxIterations), escape(escape) {}

    virtual Real getFieldValue(const VEC2F& pos) const override {
        COMPLEX iterate(pos[0], pos[1]);

        // First iteration
        iterate = map->getFieldValue(iterate);
        Real magnitude = complexMagnitude(iterate);

        int totalIterations = 1;

        while (magnitude < escape && totalIterations < maxIterations) {
            // Evaluate polynomial
            iterate = (secondaryMap == 0? map : secondaryMap)->getFieldValue(iterate);
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

    static Real findMaxValInShell(const FieldFunction2D& ff, int numRings, int numAngles, Real minRad, Real maxRad) {
        Real maxFound = numeric_limits<Real>::lowest();

        for (Real theta = 0; theta <= M_2_PI; theta+=(M_2_PI/numAngles)) {
            for (Real rad = minRad; rad <= maxRad; rad+=(maxRad-minRad)/numRings) {
                VEC2F samplePoint = rad * VEC2F(cos(theta), sin(theta));
                Real val = ff(samplePoint);
                if (val  > maxFound) maxFound = val;
            }
        }

        return maxFound;
    }

    static Real getPercentAttractingInRadius(const FieldFunction2D& ff, int numRings, int numAngles, Real minRad, Real maxRad) {
        int totalSamples = 0, attractingSamples = 0;

        for (Real theta = 0; theta <= M_2_PI; theta+=(M_2_PI/numAngles)) {
            for (Real rad = minRad; rad < maxRad; rad+=(maxRad-minRad)/numRings) {
                VEC2F samplePoint = rad * VEC2F(cos(theta), sin(theta));
                Real val = ff(samplePoint);
                if (val  < rad) attractingSamples++;
                totalSamples++;
            }
        }

        return (Real) attractingSamples / totalSamples;
    }

    static vector<pair<Real, Real>> samplePercentAttracting(const FieldFunction2D& ff, int numRings, int numAngles, Real minRad, Real maxRad, Real radDelta) {
        vector<pair<Real, Real>> samples;
        for (Real rad = minRad; rad < maxRad; rad += radDelta) {
            Real val = getPercentAttractingInRadius(ff, numRings, numAngles, 0, rad);
            samples.push_back(pair<Real, Real>(rad,val));
            // PRINTF("%.4f\t=>\t%.4f", rad, val);
        }
        return samples;
    }

    static Real findMaxAttractingRadius(const FieldFunction2D& ff, int angularRes, int shellRes, Real min, Real escape, Real radDelta) {
        // PB_START("Sampling radii");

        Real maxFound = numeric_limits<Real>::lowest();
        Real maxRadius = 0;

        for (Real rad = min+radDelta; rad < escape; rad += radDelta) {
            Real val = findMaxValInShell(ff, shellRes, angularRes, rad-radDelta, rad);
            if (val > maxFound) maxFound = val;


            if (val <= rad) { // Is attracting radius
                if (rad > maxRadius) maxRadius = rad;
                // PRINTF("Max value in radius %.4f:\t %.4f\t***", rad, maxFound);
            } else {
                // PRINTF("Max value in radius %.4f:\t %.4f", rad, maxFound);
            }

            // PB_PROGRESS((rad - min) / (escape - min));
        }

        // PB_END();

        return maxRadius;
    }

    static vector<Real> samplePercentInside(const FieldFunction2D& distField, int angularRes, Real maxRad, Real radDelta) {
        vector<Real> out{};
        for(Real rad = 0; rad < maxRad; rad+=radDelta) {
            int numSamples=0, numInside=0;
            for (Real theta = 0; theta <= M_2_PI; theta+=(M_2_PI/angularRes)) {
                VEC2F samplePoint = rad * VEC2F(cos(theta), sin(theta));
                if (distField(samplePoint) < 0) {
                    numInside++;
                }
                numSamples++;
            }

            out.push_back((Real) numInside / numSamples);
        }

        return out;
    }

    static vector<Real> sampleAverageDistance(const FieldFunction2D& distField, int angularRes, Real maxRad, Real radDelta) {
        vector<Real> out{};
        for(Real rad = 0; rad < maxRad; rad+=radDelta) {
            Real totalDist=0, numSamples=0;
            for (Real theta = 0; theta <= M_2_PI; theta+=(M_2_PI/angularRes)) {
                VEC2F samplePoint = rad * VEC2F(cos(theta), sin(theta));
                totalDist += distField(samplePoint);
                numSamples++;
            }

            out.push_back(totalDist / numSamples);
        }

        return out;
    }

};


#endif
