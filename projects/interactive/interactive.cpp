#include <cmath>
#include <cstring>
#include <limits>
#include <iomanip>
#include "SETTINGS.h"
#include "field.h"
#include "julia.h"

#include "AntTweakBar/include/AntTweakBar.h"

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

string windowLabel("FIELD_3D Viewer");

// underlying field
Grid2D *originalField;

// cached current field
ArrayGrid2D *cachedViewingField;
// scaled, biased and manipulated field
ArrayGrid2D *viewingField;

int xScreenRes = 800;
int yScreenRes = 800;
int xMouse, yMouse;
int mouseButton;
int mouseState;
uint xField = -1;
uint yField = -1;
float zoom = 1.0;
float fieldZoom = 1.0; // When resetting field bounds, zoom resets but fieldZoom does not;

float scale = 1.0;
float bias = 0.0;

bool drawingGrid = false;
bool useAbsolute = false;
bool useLog = false;
bool normalized = true;
float oldScale = 1.0;
float oldBias = 0.0;

VEC2F eyeCenter(0.5, 0.5);
VEC2F fieldEyeCenter = eyeCenter;

VEC2F fieldMin(-1,-1);
VEC2F fieldMax(1, 1);

Grid2D *distField;

// Kinda hacky to make this global but avoids having to use a capturing lambda std::function below
vector<pair<Real, Real>> radiusAttractionPercents;

// List of values sampled at consistent spacing from zero up,
// so we can have fast lookup. Really I should just make a
// Field1D class but oh well
vector<Real> GLOBAL_radInsidePercents;
Real GLOBAL_radInsidePercentSpacing;

// Same now for average distance
vector<Real> GLOBAL_radAvgDistances;
Real GLOBAL_radAvgDistanceSpacing;

// And for swap field
vector<Real> GLOBAL_radAvgDistancesSwap;
Real GLOBAL_radAvgDistanceSwapSpacing;

Real getPercentInsideForRad(Real rad) { return GLOBAL_radInsidePercents[min((uint) GLOBAL_radInsidePercents.size()-1, (uint) (rad / GLOBAL_radInsidePercentSpacing))]; }
Real getAvgDistanceForRad(Real rad) { return GLOBAL_radAvgDistances[min((uint) GLOBAL_radAvgDistances.size()-1, (uint) (rad / GLOBAL_radAvgDistanceSpacing))]; }
Real getAvgDistanceForRadSwap(Real rad) { return GLOBAL_radAvgDistancesSwap[min((uint) GLOBAL_radAvgDistancesSwap.size()-1, (uint) (rad / GLOBAL_radAvgDistanceSwapSpacing))]; }
void dumpEvenlySpacedArrToCSV(vector<Real> vec, string filename, Real spacing) {
    ofstream out;
    out.open(filename);
    if (out.is_open() == false)
        return;

    Real key = 0;
    for (Real r : vec) {
        out << key << ", " << r << endl;
        key+=spacing;
    }

    out.close();
}

typedef enum ViewMode {
    MEMBERSHIP,
    LOG_MAGNITUDE,
    DISTANCE,
    SWAP_DISTANCE,
    FIRST_ITER_RAD,
    IN_OR_OUT,
    CAPTIVITY_RATE,
    INSIDE_RATE,
    AVG_DISTANCE,
    ROTATION,
} ViewMode;

bool GLOBAL_liveView = true;
bool GLOBAL_liveZoom = true;

bool GLOBAL_drawRoots = false;
bool GLOBAL_drawOrigin = false;

typedef struct JuliaParams {
    POLYNOMIAL_2D poly;
    Real C = 40;
    Real B = 0;
    int maxIterations = 3;
    Real escape = 20;
    int res = 300;
    ViewMode mode = LOG_MAGNITUDE;
    bool dfSwapOn = false;
} JuliaParams;

typedef struct PolyParams {
    uint numRoots = 15;
    Real minPower = 5;
    Real maxPower = 10;
    bool allowNonIntegerPowers = false;
} PolyParams;

bool operator==(JuliaParams a, JuliaParams b) {
    return (
        a.poly == b.poly &&
        a.C == b.C &&
        a.B == b.B &&
        a.maxIterations == b.maxIterations &&
        a.escape == b.escape &&
        a.res == b.res &&
        a.mode == b.mode);
}

bool operator==(PolyParams a, PolyParams b) {
    return (
        a.numRoots == b.numRoots &&
        a.minPower == b.minPower &&
        a.maxPower == b.maxPower &&
        a.allowNonIntegerPowers == b.allowNonIntegerPowers);
}

PolyParams GLOBAL_polyParams;

JuliaParams GLOBAL_params;
JuliaParams GLOBAL_oldParams;


// Globals for ATW distfieldswap control
char GLOBAL_dfSwapFilename[128] = "star.sdf";
ArrayGrid2D *swapDistFieldCoarse = NULL;
Grid2D *swapDistField = NULL;

// Globals for animation export
char GLOBAL_animationScriptFilename[128] = "anim.txt";
char GLOBAL_juliaParamsFilename[128] = "params.txt";

///////////////////////////////////////////////////////////////////////
// Print a string to the GL window
///////////////////////////////////////////////////////////////////////
void printGlString(string output)
{
    //glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    glColor4f(0.0f, 1.0f, 0.0f, 1.0f);
    for (unsigned int x = 0; x < output.size(); x++)
        glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, output[x]);
}

///////////////////////////////////////////////////////////////////////
// initialize the GL texture
///////////////////////////////////////////////////////////////////////
void initTexture(Grid2D* texture)
{
    // do the dumb thing here. GL_DOUBLE doesn't seem well-supported.
    int xRes = texture->xRes;
    int yRes = texture->yRes;
    float* data = new float[xRes * yRes];

    for (int x = 0; x < xRes * yRes; x++)
        data[x] = (float)((*texture)[x]);

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, 3, texture->xRes, texture->yRes, 0, GL_LUMINANCE, GL_FLOAT, data);

    delete[] data;

    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
    glEnable(GL_TEXTURE_2D);
}

///////////////////////////////////////////////////////////////////////
// update the texture being viewed
///////////////////////////////////////////////////////////////////////
void updateViewingTexture() {

    if (viewingField != 0) delete viewingField;

    viewingField = new ArrayGrid2D(*cachedViewingField);

    if (useAbsolute)
        viewingField->abs();

    *viewingField += bias;
    *viewingField *= scale;

    if (useLog)
        viewingField->log(10.0);

    initTexture(viewingField);
}


VEC2F fieldSpaceToScreenSpace(VEC2F p) {
    VEC2F halfInvFieldZoom(1/fieldZoom/2, 1/fieldZoom/2);

    VEC2F fieldOffset = ((fieldMax + fieldMin)/2.0) / fieldZoom / 2;

    // These two points are (0, 0) and (1, 1) respectively
    VEC2F fieldScreenMin = VEC2F(0.5, 0.5) - fieldOffset;
    VEC2F fieldScreenMax = VEC2F(0.5, 0.5) + halfInvFieldZoom - fieldOffset;

    return fieldScreenMin + p.cwiseProduct(fieldScreenMax - fieldScreenMin);
}


///////////////////////////////////////////////////////////////////////
// GL and GLUT callbacks
///////////////////////////////////////////////////////////////////////
//Forward declare
void normalize();
void recomputeJuliaSet(int res=-1);
void glutDisplay()
{

    // Check if params have updated
    if (GLOBAL_liveView && !(GLOBAL_params == GLOBAL_oldParams)) {
        recomputeJuliaSet();
        updateViewingTexture();
        if (normalized) normalize();
        GLOBAL_oldParams = GLOBAL_params;
    }

    // Make ensuing transforms affect the projection matrix
    glMatrixMode(GL_PROJECTION);

    // set the projection matrix to an orthographic view
    glLoadIdentity();
    float halfZoom = zoom * 0.5;

    glOrtho(-halfZoom, halfZoom, -halfZoom, halfZoom, -10, 10);

    // set the matric mode back to modelview
    glMatrixMode(GL_MODELVIEW);

    // set the lookat transform
    glLoadIdentity();
    gluLookAt(eyeCenter[0], eyeCenter[1], 1,  // eye
        eyeCenter[0], eyeCenter[1], 0,  // center
        0, 1, 0);   // up

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    Real xLength = 1.0;
    Real yLength = 1.0;

    int xRes = originalField->xRes;
    int yRes = originalField->yRes;
    if (xRes < yRes)
        xLength = (Real)xRes / yRes;
    if (yRes < xRes)
        yLength = (Real)yRes / xRes;

    glEnable(GL_TEXTURE_2D);
    glBegin(GL_QUADS);
    glTexCoord2f(0.0, 0.0); glVertex3f(0.0, 0.0, 0.0);
    glTexCoord2f(0.0, 1.0); glVertex3f(0.0, yLength, 0.0);
    glTexCoord2f(1.0, 1.0); glVertex3f(xLength, yLength, 0.0);
    glTexCoord2f(1.0, 0.0); glVertex3f(xLength, 0.0, 0.0);
    glEnd();
    glDisable(GL_TEXTURE_2D);

    if (GLOBAL_drawRoots) {

        for(COMPLEX c : GLOBAL_params.poly.roots) {
            VEC2F v(real(c), imag(c));
            VEC2F vs = fieldSpaceToScreenSpace(v);

            glColor3f(1, 0, 0);
            glPointSize(10);
            glBegin(GL_POINTS);
            glVertex2f(vs[0], vs[1]);
            glEnd();
        }
    }

    if (GLOBAL_drawOrigin) {
        VEC2F os = fieldSpaceToScreenSpace(VEC2F(0,0));

        glColor3f(0, 1, 0);
        glPointSize(10);
        glBegin(GL_POINTS);
        glVertex2f(os[0], os[1]);
        glEnd();
    }



    // if there's a valid field index, print it
    if (xField >= 0 && yField >= 0 && xField < originalField->xRes && yField < originalField->yRes) {
        glLoadIdentity();

        // must set color before setting raster position, otherwise it won't take
        //glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        glColor4f(0.0f, 1.0f, 0.0f, 1.0f);

        // normalized screen coordinates (-0.5, 0.5), due to the glLoadIdentity
        float halfZoom = 0.5 * zoom;
        glRasterPos3f(-halfZoom* 0.95, -halfZoom* 0.95, 0);

        // build the field value string
        char buffer[256];

        // add the integer index
        int index = xField + yField * originalField->xRes;
        string fieldValue;
        sprintf(buffer, "%i = (", index);
        fieldValue = fieldValue + string(buffer);
        sprintf(buffer, "%i", xField);
        fieldValue = fieldValue + string(buffer);
        sprintf(buffer, "%i", yField);
        fieldValue = fieldValue + string(", ") + string(buffer) + string(") = ");

        // add the global position
        VEC2F position = fieldMin + (VEC2F(xField, yField).cwiseQuotient(VEC2F(originalField->xRes, originalField->yRes))).cwiseProduct(fieldMax - fieldMin);

        sprintf(buffer, "%f", position[0]);
        fieldValue = fieldValue + string("(") + string(buffer);
        sprintf(buffer, "%f", position[1]);
        fieldValue = fieldValue + string(",") + string(buffer) + string(") = ");

        Real value = originalField->get(xField, yField);

        if (isnan(value))
            sprintf(buffer, "nan");
        else
            sprintf(buffer, "%.10f", value);
        fieldValue = fieldValue + string(buffer);

        printGlString(fieldValue);
    }

    TwDraw();

    glutSwapBuffers();
}

///////////////////////////////////////////////////////////////////////
// animate and display new result
///////////////////////////////////////////////////////////////////////
void glutIdle()
{
    glutDisplay();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void printCommands() {
    cout << " q           - quits" << endl;
    cout << " left mouse  - pan around" << endl;
    cout << " right mouse - zoom in and out " << endl;
    cout << " right arrow - increase bias " << endl;
    cout << " left arrow  - decrease bias " << endl;
    cout << " up arrow    - increase scale " << endl;
    cout << " down arrow  - decrease scale " << endl;
    cout << " n           - normalize (auto-set scale and bias) " << endl;
    cout << " g           - throw a grid over the pixels " << endl;
    cout << " a           - take absolute value of cells " << endl;
    cout << " l           - take log of cells " << endl;
    cout << " r           - look at real component of FFT" << endl;
    cout << " i           - look at imaginary component of FFT" << endl;
    cout << " s           - look at spatial (non-FFT)" << endl;
    cout << " m           - print min and max of field" << endl;
    cout << " M           - print min and max of current z slice" << endl;
}

///////////////////////////////////////////////////////////////////////
// normalize the texture
///////////////////////////////////////////////////////////////////////
void normalize() {

    Real minFound = numeric_limits<Real>::max();
    Real maxFound = numeric_limits<Real>::min();

    for (uint i = 0; i < cachedViewingField->totalCells(); ++i) {
        if ((*cachedViewingField)[i] < minFound) minFound = (*cachedViewingField)[i];
        if ((*cachedViewingField)[i] > maxFound) maxFound = (*cachedViewingField)[i];
    }

    // cache the values in case we want to undo the normalization
    oldScale = scale;
    oldBias = bias;

    // go ahead and compute the normalized version
    bias = -minFound;
    scale = 1.0 / (maxFound - minFound);

    updateViewingTexture();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutSpecial(int key, int x, int y) {
    if (!TwEventSpecialGLUT(key, x, y)) {
        switch (key)
        {
        case GLUT_KEY_LEFT:
            bias -= 0.01;
            break;
        case GLUT_KEY_RIGHT:
            bias += 0.01;
            break;
        case GLUT_KEY_UP:
            scale += 0.01;
            break;
        case GLUT_KEY_DOWN:
            scale -= 0.01;
            break;
        }
        cout << " scale: " << scale << " bias: " << bias << endl;

        updateViewingTexture();
    }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutKeyboard(unsigned char key, int x, int y)
{
    if (!TwEventKeyboardGLUT(key, x, y)) {

        switch (key)
        {
        case 'a':
            useAbsolute = !useAbsolute;
            updateViewingTexture();
            break;
        case 'l':
            useLog = !useLog;
            updateViewingTexture();
            break;
        case 'g':
            drawingGrid = !drawingGrid;
            break;
        case 'n':
            if (normalized) {
                scale = oldScale;
                bias = oldBias;
                updateViewingTexture();
            }
            else
                normalize();

            normalized = !normalized;
            break;
        case '?':
            printCommands();
            break;
        case 'q':
            exit(0);
            break;
        default:
            break;
        }

    }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TW_CALL setFieldBoundsCallback(void *clientData);
void glutMouseClick(int button, int state, int x, int y) {
    if (!TwEventMouseButtonGLUT(button, state, x, y)) {
        xMouse = x;
        yMouse = y;

        mouseButton = button;
        mouseState = state;

        if (mouseState == GLUT_UP && GLOBAL_liveZoom) {
            setFieldBoundsCallback(0);
        }
    }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutMouseMotion(int x, int y) {
    if (!TwEventMouseMotionGLUT(x, y)) {
        float xDiff = x - xMouse;
        float yDiff = y - yMouse;
        float speed = 0.001;

        if (mouseButton == GLUT_LEFT_BUTTON)
        {
            eyeCenter[0] -= xDiff * speed * zoom;
            eyeCenter[1] += yDiff * speed * zoom;
        }
        if (mouseButton == GLUT_RIGHT_BUTTON)
        {
            zoom -= yDiff * speed;
        }

        xMouse = x;
        yMouse = y;

    }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutPassiveMouseMotion(int x, int y) {

    if (!TwEventMouseMotionGLUT(x, y)) {
        // make the lower left the origin
        y = yScreenRes - y;

        float xNorm = (float)x / xScreenRes;
        float yNorm = (float)y / yScreenRes;

        float halfZoom = 0.5 * zoom;
        float xWorldMin = eyeCenter[0] - halfZoom;
        float xWorldMax = eyeCenter[0] + halfZoom;

        // get the bounds of the field in screen coordinates
        // if non-square textures are ever supported, change the 0.0 and 1.0 below
        float xMin = (0.0 - xWorldMin) / (xWorldMax - xWorldMin);
        float xMax = (1.0 - xWorldMin) / (xWorldMax - xWorldMin);

        float yWorldMin = eyeCenter[1] - halfZoom;
        float yWorldMax = eyeCenter[1] + halfZoom;

        float yMin = (0.0 - yWorldMin) / (yWorldMax - yWorldMin);
        float yMax = (1.0 - yWorldMin) / (yWorldMax - yWorldMin);

        int xRes = originalField->xRes;
        int yRes = originalField->yRes;

        Real xScale = 1.0;
        Real yScale = 1.0;

        if (xRes < yRes)
            xScale = (Real)yRes / xRes;
        if (xRes > yRes)
            yScale = (Real)xRes / yRes;

        // index into the field after normalizing according to screen
        // coordinates
        xField = xScale * xRes * ((xNorm - xMin) / (xMax - xMin));
        yField = yScale * yRes * ((yNorm - yMin) / (yMax - yMin));
    }
}

Real closestLookupInRealRealPairs(Real key, vector<pair<Real, Real>> pairs) {
    for (uint i = 0; i < pairs.size(); ++i) {
        if (pairs[i].first > key) {
            return pairs[i].second;
        }
    }
    return pairs[pairs.size()-1].second;
}

// Radius-distance relation is used in a few different viewing modes so
// we split it off here
Real radFn(VEC2F pt) {
    return exp(((*distField)(pt) + GLOBAL_params.B) * GLOBAL_params.C);
}
// And a second one for the swap field.
Real radFnSwap(VEC2F pt) {
    return exp(((*swapDistField)(pt) + GLOBAL_params.B) * GLOBAL_params.C);
}



void recomputeJuliaSet(int res) { //XXX default param value was given above in forward decl
    if (cachedViewingField != 0) delete cachedViewingField;

    if (res == -1) res = GLOBAL_params.res;

    // Julia set viewing
    if (GLOBAL_params.mode == ViewMode::LOG_MAGNITUDE || GLOBAL_params.mode == ViewMode::MEMBERSHIP) {
        if (GLOBAL_params.dfSwapOn) {
            // WITH SWAPPING
            DistanceGuidedMap dMap(distField, &GLOBAL_params.poly, GLOBAL_params.C, GLOBAL_params.B);
            DistanceGuidedMap dMapSwap(swapDistField, &GLOBAL_params.poly, GLOBAL_params.C, GLOBAL_params.B);

            JuliaSet julia(&dMap, &dMapSwap, GLOBAL_params.maxIterations, GLOBAL_params.escape);
            julia.mode = (GLOBAL_params.mode == MEMBERSHIP)? JuliaSet::SET_MEMBERSHIP: JuliaSet::LOG_MAGNITUDE;

            cachedViewingField = new ArrayGrid2D(res, res, fieldMin, fieldMax, &julia);
            originalField = cachedViewingField;

        } else {
            // NORMAL
            DistanceGuidedMap dMap(distField, &GLOBAL_params.poly, GLOBAL_params.C, GLOBAL_params.B);

            JuliaSet julia(&dMap, GLOBAL_params.maxIterations, GLOBAL_params.escape);
            julia.mode = (GLOBAL_params.mode == MEMBERSHIP)? JuliaSet::SET_MEMBERSHIP: JuliaSet::LOG_MAGNITUDE;

            cachedViewingField = new ArrayGrid2D(res, res, fieldMin, fieldMax, &julia);
            originalField = cachedViewingField;

        }
    }

    // View distance field directly
    else if (GLOBAL_params.mode == DISTANCE) {
        cachedViewingField = new ArrayGrid2D(res, res, fieldMin, fieldMax, distField);
        originalField = cachedViewingField;
    }

    // View (secondary) distance field directly, but only if it exists to avoid crash
    else if (GLOBAL_params.mode == SWAP_DISTANCE) {
        cachedViewingField = new ArrayGrid2D(res, res, fieldMin, fieldMax, swapDistField==0? distField:swapDistField);
        originalField = cachedViewingField;
    }

    // Plot radius after one iteration
    else if (GLOBAL_params.mode == FIRST_ITER_RAD) {
        FieldFunction2D radField(radFn);

        cachedViewingField = new ArrayGrid2D(res, res, fieldMin, fieldMax, &radField);
        originalField = cachedViewingField;
    }

    // Plot % of pts inside the circle of R(i+1), so polynomial agnostic, whose next-iteration radius is <R(i+1).
    // If this is =1.0, then the iterate is guaranteed to never leave that circle
    else if (GLOBAL_params.mode == CAPTIVITY_RATE) {
        FieldFunction2D radField(radFn);
        radiusAttractionPercents = DistanceMapInspection::samplePercentAttracting(radField, 10, 20, 0, 4, 0.01); //XXX really should be adjustable
        FieldFunction2D captivityField( [](VEC2F pt) { return closestLookupInRealRealPairs(radFn(pt), radiusAttractionPercents); } );

        cachedViewingField = new ArrayGrid2D(res, res, fieldMin, fieldMax, &captivityField);
        originalField = cachedViewingField;
    }

    // Plot whether the point moves towards the origin or away from the origin on the next iteration.
    else if (GLOBAL_params.mode == IN_OR_OUT) {
        FieldFunction2D inOutField( [](VEC2F pt) { return (radFn(pt) < pt.norm()) ? 1.0:0.0; } );

        cachedViewingField = new ArrayGrid2D(res, res, fieldMin, fieldMax, &inOutField);
        originalField = cachedViewingField;
    }

    // Plot the % of pts along the circle of R(i+1), so polynomial agnostic, which lie inside the target shape.
    // The intuition here is that points inside the target shape are likely to stay inside the target shape and vice versa,
    // so this gives us somewhat of a view of the "shell of chaos" along the target surface.
    else if (GLOBAL_params.mode == INSIDE_RATE) {
        FieldFunction2D insideTargetField( [](VEC2F pt) { return getPercentInsideForRad(radFn(pt)); } );

        cachedViewingField = new ArrayGrid2D(res, res, fieldMin, fieldMax, &insideTargetField);
        originalField = cachedViewingField;
    }

    // Plot the avg signed distance along the circle of R(i+1), so polynomial agnostic.
    else if (GLOBAL_params.mode == AVG_DISTANCE) {
        if (GLOBAL_params.dfSwapOn && swapDistField != 0) {
            FieldFunction2D avgDistanceField( [](VEC2F pt) { return getAvgDistanceForRadSwap(radFn(pt)); } );

            cachedViewingField = new ArrayGrid2D(res, res, fieldMin, fieldMax, &avgDistanceField);
            originalField = cachedViewingField;
        } else {
            FieldFunction2D avgDistanceField( [](VEC2F pt) { return getAvgDistanceForRad(radFn(pt)); } );

            cachedViewingField = new ArrayGrid2D(res, res, fieldMin, fieldMax, &avgDistanceField);
            originalField = cachedViewingField;
        }
    }

    // Plot the rotation caused by the polynomial at each point in space, as a scalar corresponding to a CCW angle.
    else if (GLOBAL_params.mode == ROTATION) {
        FieldFunction2D rotationField( [](VEC2F pt) {
            COMPLEX out = GLOBAL_params.poly.evaluate(COMPLEX(pt[0], pt[1]));
            VEC2F outV2 = VEC2F(real(out), imag(out));
            return acos(outV2.dot(pt)/(outV2.norm() * pt.norm()));
            });

        cachedViewingField = new ArrayGrid2D(res, res, fieldMin, fieldMax, &rotationField);
        originalField = cachedViewingField;
    }

    // for (pair<Real, Real> p : samples) {
    //     PRINTF("%.4f:\t%.4f", p.first, p.second);
    // }

    // PRINTD(DistanceMapInspection::findMaxAttractingRadius(radField, 300, 5, 0, 1, 1e-3));

    // DistanceMapInspection::getMaxAttractingRadius(radField, 300);
    // DistanceMapInspection::sampleAttractingPercent(radField, 300, 0, 0.25, 1e-3);

}

void TW_CALL setFieldBoundsCallback(void *clientData) {

    fieldZoom *= zoom;

    VEC2F fieldCenter = fieldMin + eyeCenter.cwiseProduct((fieldMax - fieldMin));

    fieldMin = fieldCenter - VEC2F(fieldZoom, fieldZoom);
    fieldMax = fieldCenter + VEC2F(fieldZoom, fieldZoom);

    zoom = 1;
    eyeCenter = VEC2F(0.5,0.5);

    recomputeJuliaSet();
    updateViewingTexture();
    if (normalized) normalize();
}

void TW_CALL recomputeCallback(void *clientData) {
    recomputeJuliaSet();
    updateViewingTexture();
    if (normalized) normalize();
}

void TW_CALL resetFieldBoundsCallback(void *clientData) {

    fieldZoom = 1;
    zoom = 1;

    fieldMin = VEC2F(-1, -1);
    fieldMax = VEC2F(1, 1);

    recomputeJuliaSet();
    updateViewingTexture();
    if (normalized) normalize();
}

void TW_CALL newPolynomialCallback(void *clientData) {
    GLOBAL_params.poly = POLYNOMIAL_2D::randPolynomialInBox(AABB_2D(VEC2F(-0.5,-0.5), VEC2F(0.5,0.5)),
        GLOBAL_polyParams.minPower,
        GLOBAL_polyParams.maxPower,
        GLOBAL_polyParams.numRoots,
        GLOBAL_polyParams.allowNonIntegerPowers);
    // recomputeJuliaSet();
    // updateViewingTexture();
    // if (normalized) normalize();
}

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

void TW_CALL loadSwapFieldCallback(void *clientData) {
    if (swapDistField != 0) delete swapDistField;
    if (swapDistFieldCoarse != 0) delete swapDistFieldCoarse;
    swapDistFieldCoarse = new ArrayGrid2D(GLOBAL_dfSwapFilename);
    swapDistField = new InterpolationGrid2D(swapDistFieldCoarse);

    // Shift the distance field so that the minimum distance is at the origin
    VEC2F dfCenter = findDistanceFieldCenter(*swapDistField);
    swapDistField->mapBox.setCenter(-dfCenter); // XXX I'm not sure why this has to be negative. Something might be afoot.

    // Compute the radius-distance relation
    GLOBAL_radAvgDistanceSwapSpacing = 0.0005;
    GLOBAL_radAvgDistancesSwap = DistanceMapInspection::sampleAverageDistance(*swapDistField, 360, 0.25, GLOBAL_radAvgDistanceSpacing);
}

void TW_CALL writeAnimFrameCallback(void *clientData) {
    ofstream out;
    out.open(GLOBAL_animationScriptFilename, ios_base::app);
    if (out.is_open() == false)
        return;

    out << setprecision(17) << fieldMin[0] << ", " << fieldMin[1] << ", " << fieldMax[0] << ", " << fieldMax[1] << endl;

    out.close();
}

void TW_CALL writeParamsCallback(void *clientData) {
    ofstream out;
    out.open(GLOBAL_juliaParamsFilename);
    if (out.is_open() == false)
        return;

    out << GLOBAL_params.C << endl;
    out << GLOBAL_params.B << endl;
    out << GLOBAL_params.maxIterations << endl;
    out << GLOBAL_params.escape << endl;

    out << GLOBAL_params.poly.roots.size() << endl;
    for (int i = 0; i < GLOBAL_params.poly.roots.size(); ++i) {
        out << GLOBAL_params.poly.roots[i].real() << ", " << GLOBAL_params.poly.roots[i].imag() << ", " << GLOBAL_params.poly.powers[i] << endl;
    }

    out.close();
}

//////////////////////////////////////////////////////////////////////////////
// open the GLVU window
//////////////////////////////////////////////////////////////////////////////
int glvuWindow() {
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE| GLUT_RGBA);
    glutInitWindowSize(xScreenRes, yScreenRes);
    glutInitWindowPosition(10, 10);
    //glutCreateWindow("FIELD_3D Viewer");
    glutCreateWindow(windowLabel.c_str());

    // set the viewport resolution (w x h)
    glViewport(0, 0, (GLsizei) xScreenRes, (GLsizei) yScreenRes);

    glClearColor(0.1, 0.1, 0.1, 0);
    glutDisplayFunc(&glutDisplay);
    glutIdleFunc(&glutIdle);
    glutKeyboardFunc(&glutKeyboard);
    glutSpecialFunc(&glutSpecial);
    glutMouseFunc(&glutMouseClick);
    glutMotionFunc(&glutMouseMotion);
    glutPassiveMotionFunc(&glutPassiveMouseMotion);

    updateViewingTexture();

    TwInit(TW_OPENGL, NULL);

    TwGLUTModifiersFunc(glutGetModifiers);

    TwWindowSize(xScreenRes, yScreenRes);

    TwBar *juliaBar;
    juliaBar = TwNewBar("params");
    TwDefine("params label='Julia set parameters' size='250 270' valueswidth='100'");

    TwEnumVal ViewModeEV[] = {
        {MEMBERSHIP, "Set membership"},
        {LOG_MAGNITUDE, "Log magnitude"},
        {DISTANCE, "Signed distance field"},
        {SWAP_DISTANCE, "Secondary SDF"},
        {FIRST_ITER_RAD, "First iteration radius"},
        {CAPTIVITY_RATE, "R(i+1) captivity rate"},
        {INSIDE_RATE, "R(i+1) \% inside target"},
        {AVG_DISTANCE, "R(i+1) average distance from surface"},
        {IN_OR_OUT, "R(i+1) < R(i)?"},
        {ROTATION, "Rotation field"}
    };
    TwType viewEnumType = TwDefineEnum("View Mode", ViewModeEV, 10);

    TwAddVarRW(juliaBar, "View mode", viewEnumType, &GLOBAL_params.mode, NULL);
    TwAddSeparator(juliaBar, "sep1", NULL);
    TwAddVarRW(juliaBar, "C", TW_TYPE_DOUBLE, &GLOBAL_params.C, "min=1 max=1000 step=1");
    TwAddVarRW(juliaBar, "B", TW_TYPE_DOUBLE, &GLOBAL_params.B, "min=-5 max=5 step=0.01");
    TwAddVarRW(juliaBar, "Max iterations", TW_TYPE_INT32, &GLOBAL_params.maxIterations, "min=1 max=30 step=1");
    TwAddVarRW(juliaBar, "Resolution", TW_TYPE_INT32, &GLOBAL_params.res, "min=20 max=1000 step=10");
    TwAddVarRW(juliaBar, "Escape radius", TW_TYPE_DOUBLE, &GLOBAL_params.escape, "min=1 step=0.5");
    TwAddSeparator(juliaBar, "sep2", NULL);
    TwAddVarRW(juliaBar, "Live field bounds", TW_TYPE_BOOL8, &GLOBAL_liveZoom, NULL);
    TwAddButton(juliaBar, "Zoom field bounds", setFieldBoundsCallback, NULL, NULL);
    TwAddButton(juliaBar, "Reset field bounds", resetFieldBoundsCallback, NULL, NULL);
    TwAddVarRW(juliaBar, "Live update", TW_TYPE_BOOL8, &GLOBAL_liveView, NULL);
    TwAddButton(juliaBar, "Recompute", recomputeCallback, NULL, NULL);
    TwAddVarRW(juliaBar, "Draw origin", TW_TYPE_BOOL8, &GLOBAL_drawOrigin, NULL);

    TwBar *polyBar;
    polyBar = TwNewBar("poly");
    TwDefine("poly label='Polynomial parameters' size='300 120' valueswidth='50' iconified=true");
    TwAddVarRW(polyBar, "Num roots", TW_TYPE_INT32, &GLOBAL_polyParams.numRoots, "min=1 max=30 step=1");
    TwAddVarRW(polyBar, "Min power", TW_TYPE_DOUBLE, &GLOBAL_polyParams.minPower, "min=1 max=100");
    TwAddVarRW(polyBar, "Max power", TW_TYPE_DOUBLE, &GLOBAL_polyParams.maxPower, "min=1 max=100");
    TwAddVarRW(polyBar, "Allow noninteger powers", TW_TYPE_BOOL8, &GLOBAL_polyParams.allowNonIntegerPowers, NULL);
    TwAddButton(polyBar, "Generate new polynomial", newPolynomialCallback, NULL, NULL);
    TwAddVarRW(polyBar, "Draw root positions", TW_TYPE_BOOL8, &GLOBAL_drawRoots, NULL);

    TwBar *dfSubBar;
    dfSubBar = TwNewBar("dfsub");
    TwDefine("dfsub label='Distance field swap' size='200 100' valueswidth='100' iconified=true");
    TwAddVarRW(dfSubBar, "Use alternate distance field for i>1", TW_TYPE_BOOL8, &GLOBAL_params.dfSwapOn, NULL);
    TwAddVarRW(dfSubBar, "Alternate distance field", TW_TYPE_CSSTRING(128), &GLOBAL_dfSwapFilename, "");
    TwAddButton(dfSubBar, "Load alternate distance field", loadSwapFieldCallback, NULL, NULL);

    TwBar *animBar;
    animBar = TwNewBar("anim");
    TwDefine("anim label='Export for render' size='200 100' valueswidth='100' iconified=true");
    TwAddVarRW(animBar, "Animation filename", TW_TYPE_CSSTRING(128), &GLOBAL_animationScriptFilename, "");
    TwAddButton(animBar, "Append frame to file", writeAnimFrameCallback, NULL, NULL);
    TwAddSeparator(animBar, "sep1", NULL);
    TwAddVarRW(animBar, "Params filename", TW_TYPE_CSSTRING(128), &GLOBAL_juliaParamsFilename, "");
    TwAddButton(animBar, "Write params to file", writeParamsCallback, NULL, NULL);

    glutMainLoop();

    // Control flow will never reach here
    return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {
    if(argc < 2 || argc > 5) {
        cout << "USAGE: " << endl;
        cout << "To open an interactive visualization of a shaped Julia set from a distance field:" << endl;
        cout << " " << argv[0] << " <distance field> <output res (optional)> <fill level (optional)> <second distance field (optional)>" << endl;
        exit(0);
    }

    ArrayGrid2D *distFieldCoarse = new ArrayGrid2D((string(argv[1])));
    distFieldCoarse->setMapBox(AABB_2D(VEC2F(-0.5,-0.5), VEC2F(0.5,0.5))); // Set to size 1, centered at origin
    distField = new InterpolationGrid2D(distFieldCoarse, InterpolationGrid2D::LINEAR); // Interpolate

    // Shift the distance field so that the minimum distance is at the origin
    VEC2F dfCenter = findDistanceFieldCenter(*distField);
    PRINTV2(dfCenter);
    PRINTI(distField->hasMapBox);
    PRINTV2(distField->mapBox.min());
    PRINTV2(distField->mapBox.max());

    distField->mapBox.min() -= dfCenter;
    distField->mapBox.max() -= dfCenter;

    // Compute the radius-inside relations once, at the beginning
    GLOBAL_radInsidePercentSpacing = 0.0005;
    GLOBAL_radInsidePercents = DistanceMapInspection::samplePercentInside(*distField, 360, 0.25, GLOBAL_radInsidePercentSpacing);

    GLOBAL_radAvgDistanceSpacing = 0.0005;
    GLOBAL_radAvgDistances = DistanceMapInspection::sampleAverageDistance(*distField, 360, 0.25, GLOBAL_radAvgDistanceSpacing);
    // dumpEvenlySpacedArrToCSV(GLOBAL_radAvgDistances, "distances.csv", GLOBAL_radAvgDistanceSpacing);

    // distField->writePPM("sdf.ppm");

    newPolynomialCallback(NULL);

    if (argc > 2) GLOBAL_params.res = atoi(argv[2]);
    if (argc > 3) GLOBAL_params.C = atof(argv[3]);

    if (argc > 4) {
        strncpy(GLOBAL_dfSwapFilename, argv[4], 128);
        GLOBAL_params.dfSwapOn = true;
        loadSwapFieldCallback(NULL);
    }

    recomputeJuliaSet();

    windowLabel = string(argv[1]);
    glutInit(&argc, argv);
    glvuWindow();

    return 1;
}
