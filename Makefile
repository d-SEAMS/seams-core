
all:
	cd src && g++ -g -lm -Wall -o runme transition.cpp chunk.cpp density.cpp geometry.cpp output.cpp structure_factor.cpp rdf2D.cpp rdf3D.cpp molecular_system.cpp molecule.cpp parameter.cpp main.cpp

molecule.o: molecule.cpp molecule.h
	cd src && g++ -c molecule.cpp

molecular_system.o: molecular_system.cpp molecular_system.h
	cd src && g++ -c molecular_system.cpp

parameter.o: parameter.cpp
	cd src && g++ -c parameter.cpp

rdf3D.o: rdf3D.cpp
	cd src && g++ -c rdf3D.cpp

rdf2D.o: rdf2D.cpp
	cd src && g++ -c rdf2D.cpp

structure_factor.o: structure_factor.cpp
	cd src && g++ -c structure_factor.cpp

output.o: output.cpp
	cd src && g++ -c output.cpp

geometry.o: geometry.cpp
	cd src && g++ -c geometry.cpp

density.o: density.cpp
	cd src && g++ -c density.cpp

transition.o: transition.cpp
	cd src && g++ -c transition.cpp

chunk.o: chunk.cpp
	cd src && g++ -c chunk.cpp

clean:
	cd src && rm -f *.o runme
