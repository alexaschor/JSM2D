# JSM2D

<img width="800" alt="image" src="https://github.com/alexaschor/JSM2D/assets/77604978/62d6f342-eede-4cbe-a933-80e3af296a61">

## Note
This is a 2D demo implementation of [alexaschor/JuliaShapeModulus](https://github.com/alexaschor/JuliaShapeModulus). For the code used to generate the examples in the paper, check out that repo.  

## Building:
Requires OpenGL and GLUT. `./bin/interactive` uses [ANTTweakBar](http://anttweakbar.sourceforge.net/doc/), but that's bundled with this code and should just build with `make`.
1. `git clone https://github.com/alexaschor/DQM2D`
2. `cd DQM2D`
3. `mkdir bin`
4. If on Linux, go to `include.mk`, comment out the Mac part and uncomment the Linux part. For either platform you may need to change CXX to whichever `g++` you want to use.
5. `make`

## SDF from SVG (only supports single polygons):
1. `./bin/svgToPoly data/XXX.svg data/XXX.poly`
2. `./bin/polyToSDF data/XXX.poly <RESOLUTION> data/XXX.sdf`


## SDF from PPM (will flip image vertically, should be white inside shape and black everywhere else):
Uses a ["dead reckoning"](https://www.sciencedirect.com/science/article/pii/S1077314204000682) transform, so I don't think it's exact but it works fine for this purpose.
1. `./bin/ppmToSDF data/XXX.ppm <OUTPUT RESOLUTION> data/XXX.sdf`

## Generate/view a Julia set:
1. To run `./bin/interactive`, you'll need to set your `DYLD_LIBRARY_PATH` on Mac, or `LD_LIBRARY_PATH` on Linux. `make setPath` in `projects/interactive` will generate the command you need to run before executing `./bin/interactive`
2. `./bin/interactive data/XXX.sdf`
3. Most of the controls are on-screen and the rest are from `fieldViewer3D`. LMB to pan, RMB and drag to zoom. N togges auto-normalizing, and arrow keys adjust scale/bias. By default there's just the one config window open, but there's a few more minimized to the bottom left of the viewing window.
