#include <cstdlib>
#include <iostream>
#include "Globals.h"
#include "Pore.h"
#include "Throat.h"



void TerM(char *ErrorMessage) {
	NormalFinish();

	std::cout<<ErrorMessage<<std::endl;
	std::exit(EXIT_FAILURE);
}

void NormalFinish(void) {
	delete[] pores;
	delete[] throats;
#ifdef PETSC_ENABLED	
	ierr = KSPDestroy(&ksp); CHKERRXX(ierr);
	ierr = PCDestroy(&pc); CHKERRXX(ierr);
	ierr = VecDestroy(&Petscx); CHKERRXX(ierr);
	ierr = VecDestroy(&Petscb); CHKERRXX(ierr);
	ierr = MatDestroy(&PetscA); CHKERRXX(ierr);
	delete[] VectorMap;
	delete[] DisplacementMap;
#else
	delete[] CoeffMatrix;
	delete[] Ans;
	delete[] Row;
	delete[] Col;
#endif	
}