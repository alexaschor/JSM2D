CXX=g++-9
CXXFLAGS=-Wall -MMD -g -std=c++11 -fopenmp -framework GLUT -framework OpenGL -I../../lib -I../../src/ -I../../

# Additional flags for AntTweakBar library
CXXFLAGS+= -L../../lib/AntTweakBar/lib -lAntTweakBar


EXECUTABLE = ../../bin/interactive

SOURCES    = interactive.cpp

OBJECTS = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(WARNING) $(CXXFLAGS) $^ -o $@

.cpp.o:
	$(CXX) $(WARNING) $(CXXFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -f *.o

DYLIBPATH=$(CURDIR)/../../lib/AntTweakBar/lib
setPath:
	export DYLD_LIBRARY_PATH=$(DYLIBPATH)

