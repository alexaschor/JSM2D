# DQM2D

2D version of https://github.com/alexaschor/DistQuatMap

Building:
1. `git clone https://github.com/alexaschor/DQM2D`
2. `mkdir bin`
3. If on Linux, go to `include.mk`, comment out the Mac part and uncomment the Linux part. For either platform you may need to change CXX to whichever `g++` you want to use.
4. `make`

SDF from SVG (only supports single polygons):
1. `./bin/svgToPoly data/XXX.svg data/XXX.poly`
2. `./bin/polyToSDF data/XXX.poly <RESOLUTION> data/XXX.sdf`


SDF from PPM (will flip image vertically, should be white inside shape and black everywhere else):
1. `./bin/ppmToSDF data/XXX.ppm <RESOLUTION> data/XXX.sdf`

Generate a Julia set:
1. To run `./bin/interactive`, you'll need to set your `DYLD_LIBRARY_PATH` on Mac, or `LD_LIBRARY_PATH` on Linux. `make setPath` in `projects/interactive` will generate the command you need to run before executing `./bin/interactive`
2. `./bin/interactive data/XXX.sdf`
3. Most of the controls are on-screen and the rest are from `fieldViewer3D`. LMB to pan, RMB and drag to zoom. N togges auto-normalizing, and arrow keys adjust scale/bias.
