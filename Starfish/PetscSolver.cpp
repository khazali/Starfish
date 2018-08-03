#include "Globals.h"
#ifdef PETSC_ENABLED
#include <petscksp.h>
#include <mpi.h>


void PetscSolve(FloatType *xx) {
	double rtol;
	double *pans;
	PetscInt its;

	ierr = MatAssemblyBegin(PetscA, MAT_FINAL_ASSEMBLY); CHKERRXX(ierr);
	ierr = MatAssemblyEnd(PetscA, MAT_FINAL_ASSEMBLY); CHKERRXX(ierr);
	ierr = VecAssemblyBegin(Petscb); CHKERRXX(ierr);
	ierr = VecAssemblyEnd(Petscb); CHKERRXX(ierr);
	ierr = VecAssemblyBegin(Petscx); CHKERRXX(ierr);
	ierr = VecAssemblyEnd(Petscx); CHKERRXX(ierr);

	ierr = MatSetOption(PetscA, MAT_SYMMETRIC, PETSC_TRUE); CHKERRXX(ierr);
	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRXX(ierr);
	ierr = PCCreate(PETSC_COMM_WORLD, &pc); CHKERRXX(ierr);

	ierr = KSPSetOperators(ksp, PetscA, PetscA); CHKERRXX(ierr);
	ierr = PCSetOperators(pc, PetscA, PetscA); CHKERRXX(ierr);
	ierr = PCSetType(pc, PCJACOBI); CHKERRXX(ierr);
	ierr = KSPSetType(ksp, KSPCG); CHKERRXX(ierr);

	//rtol = 1e-15 / (((double)(PoreNO))*((double)(PoreNO)));
	rtol = 1e-15;
	ierr = KSPSetTolerances(ksp, rtol, 1.e-50, PETSC_DEFAULT, 10000); CHKERRXX(ierr);
	//ierr = KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED); CHKERRXX(ierr);

	KSPSetUp(ksp);
	PCSetUp(pc);
	ierr = KSPSolve(ksp, Petscb, Petscx); CHKERRXX(ierr);
	//ierr = PCView(pc, PETSC_VIEWER_STDOUT_WORLD); CHKERRXX(ierr);
	//ierr = KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD); CHKERRXX(ierr);
	//ierr = KSPGetIterationNumber(ksp, &its); CHKERRXX(ierr);
	//ierr = PetscPrintf(PETSC_COMM_WORLD, "iterations %D\n", its); CHKERRXX(ierr);
	ierr = VecGetArray(Petscx, &pans); CHKERRXX(ierr);
	MPI_Gatherv(pans, VectorMap[rank], MPI_DOUBLE, xx, VectorMap, DisplacementMap, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
	ierr = VecRestoreArray(Petscx, &pans); CHKERRXX(ierr);
	//for (int i = 0; i < PoreNO; i++) PetscPrintf(PETSC_COMM_WORLD, "%f\n", xx[i]);
	ierr = MatSetOption(PetscA, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE); CHKERRXX(ierr);
}
#endif