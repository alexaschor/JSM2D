#include <cmath>
#include <limits>
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
int xField = -1;
int yField = -1;
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

VEC2F fieldMin(-1,-1);
VEC2F fieldMax(1, 1);

Grid2D *distField;

typedef enum ViewMode {
    MEMBERSHIP,
    LOG_MAGNITUDE,
    DISTANCE,
    FIRST_ITER_RAD
} ViewMode;

bool GLOBAL_liveView = true;

typedef struct JuliaParams {
    POLYNOMIAL_2D poly;
    Real C = 2;
    int maxIterations = 3;
    Real escape = 20;
    int res = 100;
    ViewMode mode = LOG_MAGNITUDE;
} JuliaParams;

bool operator==(JuliaParams a, JuliaParams b) {
    return (
        a.poly == b.poly &&
        a.C == b.C &&
        a.maxIterations == b.maxIterations &&
        a.escape == b.escape &&
        a.res == b.res &&
        a.mode == b.mode);
}

JuliaParams GLOBAL_params;
JuliaParams GLOBAL_oldParams;

///////////////////////////////////////////////////////////////////////
// Print a string to the GL window
///////////////////////////////////////////////////////////////////////
void printGlString(string output)
{
    //glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
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

    // if there's a valid field index, print it
    if (xField >= 0 && yField >= 0 && xField < originalField->xRes && yField < originalField->yRes) {
        glLoadIdentity();

        // must set color before setting raster position, otherwise it won't take
        //glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        glColor4f(1.0f, 0.0f, 0.0f, 1.0f);

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

    for (int i = 0; i < viewingField->totalCells(); ++i) {
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
void glutMouseClick(int button, int state, int x, int y) {
    if (!TwEventMouseButtonGLUT(button, state, x, y)) {
        xMouse = x;
        yMouse = y;

        mouseButton = button;
        mouseState = state;
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
            eyeCenter[0] -= xDiff * speed;
            eyeCenter[1] += yDiff * speed;
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

void recomputeJuliaSet(int res) { //XXX default param value was given above in forward decl
    if (cachedViewingField != 0) delete cachedViewingField;

    if (res == -1) res = GLOBAL_params.res;

    if (GLOBAL_params.mode == ViewMode::LOG_MAGNITUDE || GLOBAL_params.mode == ViewMode::MEMBERSHIP) {
        DistanceGuidedMap m2(distField, &GLOBAL_params.poly, GLOBAL_params.C);

        JuliaSet julia(&m2, GLOBAL_params.maxIterations, GLOBAL_params.escape);
        julia.mode = (GLOBAL_params.mode == MEMBERSHIP)? JuliaSet::SET_MEMBERSHIP: JuliaSet::LOG_MAGNITUDE;


        cachedViewingField = new ArrayGrid2D(res, res, fieldMin, fieldMax, &julia);
        originalField = cachedViewingField;
    }

    if (GLOBAL_params.mode == DISTANCE) {
        cachedViewingField = new ArrayGrid2D(res, res, fieldMin, fieldMax, distField);
        originalField = cachedViewingField;
    }

    if (GLOBAL_params.mode == FIRST_ITER_RAD) {
        FieldFunction2D radField( [](VEC2F pt) { return exp((*distField)(pt) * GLOBAL_params.C); } );

        cachedViewingField = new ArrayGrid2D(res, res, fieldMin, fieldMax, &radField);
        originalField = cachedViewingField;
    }

    // Compute min radius
    // FieldFunction2D radField( [](VEC2F pt) { return exp((*distField)(pt) * GLOBAL_params.C); } );
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

    TwBar *myBar;
    myBar = TwNewBar("params");
    TwDefine("params label='Julia set parameters' size='320 200' valueswidth='fit'");

    TwEnumVal ViewModeEV[] = { {MEMBERSHIP, "Set membership"}, {LOG_MAGNITUDE, "Log magnitude"}, {DISTANCE, "Signed distance"}, {FIRST_ITER_RAD, "First iteration radius"} };
    TwType seasonType = TwDefineEnum("View Mode", ViewModeEV, 4);

    TwAddVarRW(myBar, "View mode", seasonType, &GLOBAL_params.mode, NULL);
    TwAddVarRW(myBar, "Live view", TW_TYPE_BOOL32, &GLOBAL_liveView, NULL);
    TwAddVarRW(myBar, "C", TW_TYPE_DOUBLE, &GLOBAL_params.C, "min=1 max=1000 step=1");
    TwAddVarRW(myBar, "Max iterations", TW_TYPE_INT32, &GLOBAL_params.maxIterations, "min=1 max=30 step=1");
    TwAddVarRW(myBar, "Resolution", TW_TYPE_INT32, &GLOBAL_params.res, "min=20 max=1000 step=1");
    TwAddVarRW(myBar, "Escape radius", TW_TYPE_DOUBLE, &GLOBAL_params.escape, "min=1 max=1000 step=0.5");
    TwAddButton(myBar, "Zoom field bounds", setFieldBoundsCallback, NULL, NULL);
    TwAddButton(myBar, "Reset field bounds", resetFieldBoundsCallback, NULL, NULL);
    TwAddButton(myBar, "Recompute", recomputeCallback, NULL, NULL);

    glutMainLoop();

    // Control flow will never reach here
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {
    if(argc != 4) {
        cout << "USAGE: " << endl;
        cout << "To create a shaped Julia set from a distance field:" << endl;
        cout << " " << argv[0] << " <distance field> <output res> <fill level>" << endl;
        exit(0);
    }

    ArrayGrid2D *distFieldCoarse = new ArrayGrid2D((string(argv[1])));
    distField = new InterpolationGrid2D(distFieldCoarse, InterpolationGrid2D::LINEAR);

    // distField->writePPM("sdf.ppm");

    // Shift it up so that (0,0) is inside the object
    distField->mapBox.min() += VEC2F(0, 0.1);
    distField->mapBox.max() += VEC2F(0, 0.1);

    GLOBAL_params.poly = POLYNOMIAL_2D::randPolynomialInBox(AABB_2D(VEC2F(-1,-1), VEC2F(1,1)), 1, 1, 30, false);

    GLOBAL_params.C = atof(argv[3]);
    GLOBAL_params.res = atoi(argv[2]);

    GLOBAL_oldParams = GLOBAL_params;

    recomputeJuliaSet();

    windowLabel = string(argv[2]);
    glutInit(&argc, argv);
    glvuWindow();

    return 1;
}
