CXX=g++-9
CXXFLAGS=-Wall -MMD -g -std=c++11 -fopenmp -I../../lib -I../../src/ -I../../

EXECUTABLE = ../../bin/ppmToSDF

SOURCES = ./ppmToSDF.cpp

OBJECTS = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(WARNING) $(CXXFLAGS) $^ -o $@

.cpp.o:
	$(CXX) $(WARNING) $(CXXFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -f *.o
