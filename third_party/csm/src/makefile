# the compiler
CC = gcc
CXX = g++

# the compiler flags
CFLAGS = -Wall -I/usr/local/include/openbabel-2.0 -g 
CXXFLAGS = -Wall -I/usr/local/include/openbabel-2.0 -g

LFLAGS = -lm -L. -L/usr/local/lib -lopenbabel -g

BIN = bin

# object files for executables
SHAPEOBJ = babelAdapter.o shape.o matrix.o Molecule.o svdcmp.o parseFunctions.o permuter.o

# the first target. All executables will be made when 'make' is
# invoked with no target
all: shape

# general rule how to compile a source file and produce an object file
%.c,%.o:
	$(CXX) $(CFLAGS) -c $<

%.cpp,%.o:
	$(CXX) $(CFLAGS) -c $<

# linking rule for the first executable
shape: $(SHAPEOBJ)
	$(CXX) $(LFLAGS)  $(SHAPEOBJ) -o $(BIN)/$@	

# make .o file dependencies
depend:
	echo -e "###\n" >> makefile
	$(CXX) -MM *.c >> makefile

# remove object files and backups
clean:
	rm -f *.o *~ $(BIN)/*

