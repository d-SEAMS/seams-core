#include "output.h"
#include <sys/stat.h>

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

	if( info.st_mode & S_IFDIR )  
    	std::cout << "Output directory exists \n";
	else
	{
		// The output directory does not exist
		// Create the folder
		int status = mkdir("./output", 0777);
		switch	(status) {
			case -1:
			std::cout << "Folder exists or you can't access it.\n";
			break;
			case 0:
			std::cout << "I created something.\n";
			break;
			default:
			std::cout << "I have died. Avenge me.\n";
		}
	}
 }


// Function for printing out a file to the output directory 
// The function prints arrays to file, with the size, and two x,y arrays, and
// desired filename as arguments
void COutput::printToFile(int nbin, double* x, double* y, const std::string& filename, const std::string& xName, const std::string& yName)
{
	std::ofstream outputFile; 

	// First check if the output directory exists or not. If it does not exist create it.
	this->createOutputDir("output");

	// Create a new file in the output directory 
	outputFile.open (("./output/" + filename + ".dat").c_str());

	if (outputFile.is_open())
	{
		// First line
		outputFile << "# " << xName << "\t" <<  yName << "\n";
		// Write out the arrays x and y to a file in the output folder
		for (int ibin=0; ibin<nbin; ibin++)
		{
			outputFile << x[ibin] << "\t" <<y[ibin] <<"\n";
		}
		// Close the file
		outputFile.close();
	}
}
