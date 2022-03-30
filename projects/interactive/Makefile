include ../include.mk

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
	echo "Run the following command before running bin/interactive:"
	export $(LD_PATH_VAR)=$(DYLIBPATH)

