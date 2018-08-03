#include <iostream>
#include <string.h>

#include "Globals.h"
#include "ElementList.h"
#include "Pore.h"
#include "Throat.h"


int __cdecl main (int argc, char* argv[]) {
	char FPath[MAX_PATH_LENGTH], FName[MAX_PATH_LENGTH];
	std::ofstream DrainageOutFile, ImbibitionOutFile;
	
	
#ifdef PETSC_ENABLED
	PetscErrorCode ierr;
	omp_set_num_threads(4);
	PetscInitialize(&argc, &argv, (char*)0, NULL);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &wsize);
	VectorMap = new int[wsize];
	DisplacementMap = new int[wsize];
	if (!rank) {
#else
	init_paralution();
	//set_omp_threads_paralution(3);

#endif
		std::cout<<"Enter the file path:\n";
		std::cin.getline(FPath, MAX_PATH_LENGTH);
		std::cout<<"Enter the file prefix:\n";
		std::cin.getline(FName, MAX_PATH_LENGTH);
		//strcpy(FPath, "data\\Berea\\");
		//strcpy(FName, "Berea");			
#ifdef PETSC_ENABLED
	}	
#endif
	
	ReadStatoilFormat(FPath, FName, DrainageOutFile, ImbibitionOutFile);
	Drainage(DrainageOutFile);
	//StoreToBinFile();
	//ReadFromBinFile();
	Imbibition(ImbibitionOutFile);

	
#ifdef PETSC_ENABLED		
	if (!rank) {
#endif
		DrainageOutFile.close();
		ImbibitionOutFile.close();
		PrintQuotes();
#ifdef PETSC_ENABLED
	}
#endif
	NormalFinish();

#ifdef PETSC_ENABLED
	ierr = PetscFinalize();
#else
	stop_paralution();
#endif
	
	return 0;
}