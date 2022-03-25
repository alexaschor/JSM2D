SHELL := /bin/bash -e
LIBS =

all: main

svgToPoly:
	cd projects/svgToPoly; make

polyToSDF:
	cd projects/polyToSDF; make

main:
	cd projects/main; make

interactive:
	cd projects/interactive; make

clean:
	@-for d in $(LIBS); do (echo -e "cd ./lib/$$d; rm *.o";cd ./lib/$$d; rm *.o; cd ../..); done
	cd projects/main; make clean
	cd projects/interactive; make clean
	cd projects/polyToSDF; make clean
	cd projects/svgToPoly; make clean
