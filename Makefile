SHELL := /bin/bash -e
LIBS =

include ./projects/include.mk

all: svgToPoly polyToSDF ppmToSDF main render ANTTweakBar interactive

svgToPoly:
	cd projects/svgToPoly; make

polyToSDF:
	cd projects/polyToSDF; make

ppmToSDF:
	cd projects/ppmToSDF; make

main:
	cd projects/main; make

interactive:
	cd projects/interactive; make

render:
	cd projects/render; make

ANTTweakBar:
	cd lib/AntTweakBar/src; make -f $(ANT_MAKEFILE_NAME)

clean:
	@-for d in $(LIBS); do (echo -e "cd ./lib/$$d; rm *.o";cd ./lib/$$d; rm *.o; cd ../..); done
	cd projects/main; make clean
	cd projects/interactive; make clean
	cd projects/render; make clean
	cd projects/polyToSDF; make clean
	cd projects/ppmToSDF; make clean
	cd projects/svgToPoly; make clean
	cd lib/AntTweakBar/src; make -f $(ANT_MAKEFILE_NAME) clean
