#ifndef GLOBALS_H
#define GLOBALS_H

#include <omp.h>
#include <fstream>
#define PETSC_ENABLED

#ifdef PETSC_ENABLED
#define PETSC_CLANGUAGE_CXX
#include <petscksp.h>
#include <mpi.h>
#else
#include <paralution.hpp>
using namespace paralution;
#endif

#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286
#define SQRT3 1.73205080756887729352744634150587236694280525381038062805580
#define G_ACC 9.80665
//#define G_ACC 0
#define MAX_STRING_LENGTH 256
#define MAX_PATH_LENGTH 1000
#define FLOATCOMPARETOL 1e-15
//#define CAPILLARYINCREMENT 1000
#define CAPILLARYINCREMENT 100
#define CAPILLARYLIMIT	40000
#define RTOL	1e-15
#define ATMOSPHERP 10.1325
//#define ATMOSPHERP 0
#define SOLVERTOL 1e-60
#define PERMCONVF 1.01325e15	//mili Darcy
#define CHANGEFACT 1e3
#define PORELENGTHFACTOR 1
#define DELTAP 100
#define MAXSOLVERITERATIONS 40000
#define SHAPETOL 1e-2
#define MAXNITERS 1000

typedef double FloatType;



enum ElementType {PORE, THROAT};
enum ImbibitionType {SPONTANEOUS, FORCED};
enum SortType {DESCENDING, ASCENDING};
enum ProcessType {DRAINAGE, IMBIBITION};
enum FluidType {WATER, OIL};
enum ExitType {INLET, OUTLET};

class Pore;
class Throat;
class ElementList;

extern Pore *pores;
extern Throat *throats;

extern FloatType Dx, Dy, Dz;				//Dimensions of the porous media
extern FloatType RhoW, RhoO, DeltaRho;		//Water, Oil Densities and their difference
extern FloatType IFT;						//Oil-Water interfacial tension
extern unsigned int PoreNO, ThroatNO;
extern FloatType TotalVolume;
extern FloatType RecedingContactAngle;
extern FloatType AdvancingContactAngle;
extern FloatType Rmin;						//Smallest radius of curvature obtained during Drainage
extern ElementList OutletThroats, InletThroats;
extern FloatType TotalOilSaturation;
extern FloatType TotalWaterSaturation;
extern FloatType OilViscosity, WaterViscosity;
#ifdef PETSC_ENABLED
extern Vec Petscx, Petscb;
extern Mat PetscA;
extern KSP ksp;
extern PetscErrorCode ierr;
extern PC pc;
extern PetscInt Istart, Iend;
extern int *VectorMap, *DisplacementMap;
extern int rank, wsize;
#else
extern FloatType *CoeffMatrix, *Ans;
//extern unsigned int *Col, *Row;
extern int *Col, *Row;
extern unsigned int TotalMSize;
#endif


extern FloatType WaterRelPerm, OilRelPerm, AbsPerm;
extern unsigned int *PorePointer, *ThroatPointer;
extern int rank;

void TerM(char *);
void NormalFinish(void);
void ReadStatoilFormat(char *, char *, std::ofstream&, std::ofstream&);
void PrintQuotes(void);
void RecursiveSweepForConnection(FluidType);
void Drainage(std::ofstream&);
void Imbibition(std::ofstream&);
void CalcRelPerm(void);
void CalcAbsPerm(void);

void StoreToBinFile(void);
void ReadFromBinFile(void);
#ifdef PETSC_ENABLED
void PetscSolve(FloatType *);
#else
FloatType ParaSolver(FloatType *, int *, int *, FloatType *, FloatType *);
#endif

#endif