#ifndef PORE_H
#define PORE_H
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "Globals.h"
#include "NetworkElement.h"
#include "MIfstream.h"
#include "Throat.h"
//#include "GarbageCollection.h"

//class Throat;
//extern Throat* throats;

class Pore: public NetworkElement {
private:

	unsigned int CoordinationNumber;
	Pore** AdjacentPores;
	Throat** ConnectingThroats;

	FloatType X, Y;			//Pore location coordinates
	
	unsigned int MatrixPlace;
	FloatType Pressure;
	
public:
	Pore();
	~Pore();
	void SetLength(FloatType);
	unsigned int ReadNode1(MIfstream&);
	void ReadNode2(MIfstream&);
	
	unsigned int TestFillReadiness(FloatType);
	unsigned int CalculateImbibitionPc(FloatType);	
	void SweepAdjacentThroats(FluidType);
	void SetMatrixPlace(unsigned int);
	unsigned int GetCoordinationNumber(void);
	FloatType GetPressure(void);
	void SetPressure(FloatType);
	void BuildMyRowForAbs(void);
	void SweepAdjacentThroatsForDeletion(ExitType);
	bool IsConnectedThroatNull(unsigned int);	
	unsigned int NumberOfConnections(void);
	void BuildMyRowForWater(void);
	void BuildMyRowForOil(void);	
};
#endif