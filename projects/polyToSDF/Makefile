CXX=g++-9
CXXFLAGS=-Wall -MMD -std=c++11 -fopenmp -I../../lib -I../../src/ -I../../

EXECUTABLE = ../../bin/polyToSDF

SOURCES = ./polyToSDF.cpp

OBJECTS = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(WARNING) $(CXXFLAGS) $^ -o $@

.cpp.o:
	$(CXX) $(WARNING) $(CXXFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -f *.o
