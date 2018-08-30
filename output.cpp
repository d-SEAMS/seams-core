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

	// if( stat( path, &info ) == 0 )
 //    	cerr << "Cannot access " << path << "\n";
	if( info.st_mode & S_IFDIR )  
    	cout << "Output directory exists \n";
	else
	{
		// The output directory does not exist
		// Create the folder
		int status = mkdir("/root/output", 0077);
		switch	(status) {
			case -1:
			cout << "Folder exists or you can't access it.\n"<<status<<endl;
			break;
			case 0:
			cout << "I created something.\n";
			break;
			default:
			cout << "I have died. Avenge me.\n"<<status<<endl;
		}
	}
 }