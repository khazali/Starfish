#ifndef THROAT_H
#define THROAT_H
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "Globals.h"
#include "NetworkElement.h"
//#include "Pore.h"
#include "MIfstream.h"
//#include "GarbageCollection.h"


//extern Pore *pores;

class Throat: public NetworkElement {
private:
	Pore *ConnectingPores[2];
	FloatType TotalLength;		//Throat Total Length (pore center to pore center)

	FloatType WaterConductancePerLength, OilConductancePerLength, TotalConductancePerLength;

	FloatType OilVelocity, WaterVelocity;

public:
	Throat();
	void ReadLink1(MIfstream&);
	void ReadLink2(MIfstream&);
	unsigned int CalculateImbibitionPc(FloatType);	
	unsigned int TestFillReadiness(FloatType);
	void SweepAdjacentPores(FluidType);
	FloatType GetTotalLength(void);
	FloatType GetOutletPorePressure(void);
	void SetWaterConductancePerLength(FloatType);
	void SetOilConductancePerLength(FloatType);
	FloatType GetWaterConductancePerLength(void);
	FloatType GetOilConductancePerLength(void);
	void SetTotalConductancePerLength(FloatType);
	FloatType GetTotalConductancePerLength(void);
	void SweepAdjacentPoresForDeletion(ExitType);
	int GetIOStat(void);
	void operator=(const Throat &);
	void CalcWaterVelocity(void);
	void CalcOilVelocity(void);
};

 
#endif