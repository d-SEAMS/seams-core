#include "output.h"

//Constructor
COutput::COutput()
{
}
// Destructor
COutput::~COutput()
{
}

// This function creates the output directory inside the
// base directory, if it does not already exist
void COutput::createOutputDir(const char *path)
{
	struct stat info;

	int statRC = stat( path, &info );
    if( statRC != 0 )
    {
		cout<< "The output directory does not yet exist.\n"; // something in path prefix is not a dir
    }

	if( info.st_mode & S_IFDIR )  
    	cout << "Output directory exists \n";
	else
	{
		// The output directory does not exist
		// Create the folder
		int status = mkdir("./output", 0077);
		switch	(status) {
			case -1:
			cout << "Folder exists or you can't access it.\n";
			break;
			case 0:
			cout << "I created something.\n";
			break;
			default:
			cout << "I have died. Avenge me.\n";
		}
	}
 }


// Function for printing out a file to the output directory 
// The function prints arrays to file, with the size, and two x,y arrays, and
// desired filename as arguments
void COutput::printToFile(int nbin, double* x, double*, const char *filename = "output")
{
	// First check if the output directory exists or not. If it does not exist create it.
	this->createOutputDir("output");

	// Write out the arrays x and y to a file in the output folder
}
