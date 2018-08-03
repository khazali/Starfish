#include "Globals.h"
#include "ElementList.h"
#ifdef PETSC_ENABLED
#include <petscksp.h>
#include <mpi.h>
#endif

Pore *pores;
Throat *throats;

FloatType Dx, Dy, Dz;				//Dimensions of the porous media
FloatType RhoW, RhoO, DeltaRho;		//Water, Oil Densities and their difference
FloatType IFT;						//Oil-Water interfacial tension
unsigned int PoreNO, ThroatNO;
FloatType TotalVolume;

FloatType RecedingContactAngle;
FloatType AdvancingContactAngle;

FloatType Rmin;

ElementList OutletThroats, InletThroats;
FloatType TotalOilSaturation;
FloatType TotalWaterSaturation;
FloatType OilViscosity, WaterViscosity;
#ifdef PETSC_ENABLED
Vec Petscx, Petscb;
Mat PetscA;
KSP ksp;
PetscErrorCode ierr;
PC pc;
PetscInt Istart, Iend;
int *VectorMap, *DisplacementMap;
int rank, wsize;
#else
FloatType *CoeffMatrix, *Ans;
//unsigned int *Col, *Row;
int *Col, *Row;
unsigned int TotalMSize;
#endif


FloatType WaterRelPerm, OilRelPerm, AbsPerm;
unsigned int *PorePointer, *ThroatPointer;

