
all:
	g++ -g -lm -Wall -o runme output.cpp rdf2D.cpp rdf3D.cpp molecular_system.cpp molecule.cpp parameter.cpp main.cpp 

molecule.o: molecule.cpp molecule.h
	g++ -c molecule.cpp

molecular_system.o: molecular_system.cpp molecular_system.h
	g++ -c molecular_system.cpp

parameter.o: parameter.cpp
	g++ -c parameter.cpp

rdf3D.o: rdf3D.cpp
	g++ -c rdf3D.cpp 

rdf2D.o: rdf2D.cpp
	g++ -c rdf2D.cpp 

output.o: output.cpp
	g++ -c output.cpp 
	
clean:
	rm -f *.o runme 
