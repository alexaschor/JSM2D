# Settings for all project Makefiles

# For OSX:
CXX=g++-9
CXXFLAGS=-Wall -MMD -g -std=c++11 -fopenmp -framework GLUT -framework OpenGL -I../../lib -I../../src/ -I../../

# For Linux:
# CXX=g++
# CXXFLAGS=-Wall -MMD -g -std=c++11 -fopenmp -lGL -lGLU -lglut -I../../lib -I../../src/ -I../../
