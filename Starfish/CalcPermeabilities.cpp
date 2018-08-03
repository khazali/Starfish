#include "Globals.h"
#include "Pore.h"
#include "ElementList.h"
#include "Throat.h"
//#include "AMG_prototypes.h"
#include <fstream>

void CalcAbsPerm(void) {
	register unsigned int i, OutletThroatIndex, OutletThroatsListLengh;
	FloatType PorePressure;
	FloatType QwT, Qww;
	FloatType *x;

	x=new FloatType[PoreNO];
#ifdef PETSC_ENABLED	
	if (!rank) {
#endif
		RecursiveSweepForConnection(WATER);
#pragma omp parallel for
		for (i = 0; i < PoreNO; i++) {
			pores[i].CalculateConductanceForAbs();
		}
#pragma omp parallel for
		for (i = 0; i < ThroatNO; i++) {
			throats[i].CalculateConductanceForAbs();
		}
#pragma omp parallel for
		for (i = 0; i < PoreNO; i++) {
			pores[i].BuildMyRowForAbs();
		}
		
#ifdef PETSC_ENABLED	
	}
#else
	for (i = 0; i < PoreNO; i++) {
		x[i] = ATMOSPHERP / CHANGEFACT;
	}
	ParaSolver(CoeffMatrix, Row, Col, Ans, x);
#endif
	//PreconditionedConjugateGradient(CoeffMatrix, Row, Col, Ans, x);
	//ParaSolver(CoeffMatrix, Row, Col, Ans, x);
	//LisSolver(CoeffMatrix, Row, Col, Ans, x);
	//ierr = VecSet(Petscx, ATMOSPHERP / CHANGEFACT); CHKERRXX(ierr);
	
	
#ifdef PETSC_ENABLED
	PetscSolve(x);
	if (!rank) {
#endif
#pragma omp parallel for
		for (i = 0; i < PoreNO; i++) {
			pores[i].SetPressure(x[i] * CHANGEFACT);
		}

		QwT = 0;
		OutletThroatsListLengh = OutletThroats.GetListLength();		
#pragma omp parallel for reduction(+:QwT)
			for (i = 0; i < OutletThroatsListLengh; i++) {
			OutletThroatIndex = OutletThroats.GetListContent(i);
			PorePressure = throats[OutletThroatIndex].GetOutletPorePressure();

			Qww = throats[OutletThroatIndex].GetWaterConductancePerLength()*(PorePressure - ATMOSPHERP);
			if (Qww < 0) Qww = 0;
			QwT += Qww;
		}
		AbsPerm = PERMCONVF*QwT*WaterViscosity*Dx / (Dy*Dz*(DELTAP));
		

		if ((AbsPerm - AbsPerm) != 0) AbsPerm = 0;
#pragma omp parallel for
		for (i = 0; i < PoreNO; i++) {
			pores[i].SetPressure(ATMOSPHERP);
		}
#ifdef PETSC_ENABLED	
	}
#endif
	delete[] x;	
}

void CalcRelPerm(void) {
	register unsigned int i, OutletThroatIndex, OutletThroatsListLengh;
	FloatType PorePressure;
	FloatType QwT, QoT, Qt;

	FloatType *x;

	x = new FloatType[PoreNO];
#ifdef PETSC_ENABLED	
	if (!rank) {
#endif
#pragma omp parallel for
		for (i = 0; i < PoreNO; i++) {
			pores[i].CalculateConductance();
		}
#pragma omp parallel for
		for (i = 0; i < ThroatNO; i++) {
			throats[i].CalculateConductance();
		}
#pragma omp parallel for
		for (i = 0; i < PoreNO; i++) {
			pores[i].BuildMyRowForOil();
		}

		

		//CGOilResidual=PreconditionedConjugateGradient(CoeffMatrix, Row, Col, Ans, x);
		//ParaSolver(CoeffMatrix, Row, Col, Ans, x);
		//LisSolver(CoeffMatrix, Row, Col, Ans, x);
#ifdef PETSC_ENABLED	
	}
#else
		for (i = 0; i < PoreNO; i++) {
			x[i] = pores[i].GetPressure() / CHANGEFACT;
			if ((x[i] - x[i]) != 0) x[i] = ATMOSPHERP / CHANGEFACT;
		}
		ParaSolver(CoeffMatrix, Row, Col, Ans, x);
#endif
	

#ifdef PETSC_ENABLED
	PetscSolve(x);
	if (!rank) {
#endif
#pragma omp parallel for
		for (i = 0; i < PoreNO; i++) {
			pores[i].SetPressure(x[i] * CHANGEFACT);
		}

		QoT = 0;
		OutletThroatsListLengh = OutletThroats.GetListLength();
#pragma omp parallel for reduction(+:QoT)
		for (i = 0; i < OutletThroatsListLengh; i++) {
			OutletThroatIndex = OutletThroats.GetListContent(i);
			PorePressure = throats[OutletThroatIndex].GetOutletPorePressure();

			Qt = throats[OutletThroatIndex].GetOilConductancePerLength()*(PorePressure - ATMOSPHERP);
			if (Qt < 0) Qt = 0;

			QoT += Qt;
		}
		//std::cout << QoT << std::endl;
		OilRelPerm = PERMCONVF*QoT*OilViscosity*Dx / (AbsPerm*Dy*Dz*(DELTAP));
		if (((OilRelPerm - OilRelPerm) != 0) || (OilRelPerm < 0)) OilRelPerm = 0;
#pragma omp parallel for
		for (i = 0; i < ThroatNO; i++) {
			throats[i].CalcOilVelocity();
		}
#pragma omp parallel for
		for (i = 0; i < PoreNO; i++) {
			pores[i].BuildMyRowForWater();
		}

		
#ifdef PETSC_ENABLED	
	}
#else
		for (i = 0; i < PoreNO; i++) {
			x[i] = pores[i].GetPressure() / CHANGEFACT;
			if ((x[i] - x[i]) != 0) x[i] = ATMOSPHERP / CHANGEFACT;		
		}
		ParaSolver(CoeffMatrix, Row, Col, Ans, x);
#endif
	


#ifdef PETSC_ENABLED
	PetscSolve(x);
	if (!rank) {
#endif
#pragma omp parallel for
		for (i = 0; i < PoreNO; i++) {
			pores[i].SetPressure(x[i] * CHANGEFACT);
		}

		QwT = 0;
		OutletThroatsListLengh = OutletThroats.GetListLength();
#pragma omp parallel for reduction(+:QwT)
		for (i = 0; i < OutletThroatsListLengh; i++) {
			OutletThroatIndex = OutletThroats.GetListContent(i);
			PorePressure = throats[OutletThroatIndex].GetOutletPorePressure();
			Qt = throats[OutletThroatIndex].GetWaterConductancePerLength()*(PorePressure - ATMOSPHERP);
			if (Qt < 0) Qt = 0;
			QwT += Qt;
		}

		//std::cout << QwT << std:: endl;

		WaterRelPerm = PERMCONVF*QwT*WaterViscosity*Dx / (AbsPerm*Dy*Dz*(DELTAP));
		if (((WaterRelPerm - WaterRelPerm) != 0) || (WaterRelPerm < 0)) WaterRelPerm = 0;
#pragma omp parallel for
		for (i = 0; i < ThroatNO; i++) {
			throats[i].CalcWaterVelocity();
		}
#ifdef PETSC_ENABLED	
	}
#endif

	delete[] x;
}