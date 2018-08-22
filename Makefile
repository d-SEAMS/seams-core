
all:
	g++ -g -lm -Wall -o runme molecular_system.cpp molecule.cpp parameter.cpp main.cpp 


molecule.o: molecule.cpp molecule.h
	g++ -c molecule.cpp

molecular_system.o: molecular_system.cpp molecular_system.h
	g++ -c molecular_system.cpp

parameter.o: parameter.cpp
	g++ -c parameter.cpp 
	
clean:
	rm -f *.o runme 
