
all:
	g++ -g -lm -Wall -o runme output.cpp analysis.cpp molecular_system.cpp molecule.cpp parameter.cpp main.cpp -std=c++17

# Try -std=c++11  

molecule.o: molecule.cpp molecule.h
	g++ -c molecule.cpp

molecular_system.o: molecular_system.cpp molecular_system.h
	g++ -c molecular_system.cpp

parameter.o: parameter.cpp
	g++ -c parameter.cpp

analysis.o: analysis.cpp
	g++ -c analysis.cpp  

output.o: output.cpp
	g++ -c output.cpp 
	
clean:
	rm -f *.o runme 
