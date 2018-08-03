#ifndef NETWORKELEMENT_H
#define NETWORKELEMENT_H

#include "Globals.h"

class NetworkElement {
protected:
	FloatType OilSaturation;	//Oil OilSaturation
	FloatType ShapeFactor;
	FloatType InscribedRadius;
	FloatType EntryCapillaryPressure;
	FloatType Fd, SumForR;

	
	//bool IsCircle;		
	//unsigned int CoordinationNumber;
	FloatType ClayVolume;
	FloatType Volume;
	FloatType Length;
	FloatType SurfaceArea;

	FloatType Z;		//Height coordinate
	FloatType Beta[3];
	int IOStat;					//Inlet=-1, Outlet=0, Normal=1;
	bool IsOilFilled;			//Element is Oil filled?
	bool IsWaterFilled;			//Element is Water Filled?
	//bool ReadyToFill;			//Any adjacent elements is filled
	unsigned int Index;
	FloatType R;

	bool IsConnectedToOutletbyWater, IsConnectedToInlet, IsConnectedToOutletbyOil;
	bool IsVisited;

	bool IsSnapOff;
	FloatType MaxAdvancingContactAngle;

	bool IsOilLayerExist[3];
	FloatType Hinge_b[3], HingeAngle[3];

	FloatType AbsoluteConductance, WaterConductance, OilConductance;

	bool NegativePCFilledPistonLike;

	bool WaterCoatingExist;

	unsigned int NumberOfSides;		//0:Circle, 3:Triangle, n:Polygon, 1&2:Invalid
	//FloatType Angles;		//For equal angles

	//FloatType Aw;
	//bool IsCaculatungPc;
public:
	void CalcProps(void);
	bool GetIsOilFilled(void);
	FloatType GetPc(void);
	void SetIndex(unsigned int);
	NetworkElement(void);
	void FillWithOil(void);
	void SetIsConnectedToOutlet(FluidType, bool);
	bool GetIsConnectedToOutlet(FluidType);
	void SetIsConnectedToInlet(bool);
	bool GetIsConnectedToInlet(void);
	void SetIsVisited(bool);
	bool GetIsVisited(void);
	void SetIsWaterFilled(bool);
	bool GetIsWaterFilled(void);
	void OilLayerExist(FloatType);
	void CalculateMaxAdvancingContactAngle(void);
	void FillWithWater(FloatType);
	FloatType GetTotalWaterSaturation(void);
	void CalculateOilSaturation(FloatType);
	void CalculateDrainagePc(void);
	void CalculateConductance(void);
	FloatType GetLength(void);
	FloatType GetAbsoluteConductance(void);
	FloatType GetOilConductance(void);
	FloatType GetWaterConductance(void);
	bool GetWaterCoatingExist (void);
	FloatType GetZ(void);
	unsigned int GetIndex(void);
	FloatType GetReducedVolume(void);
	bool IsIsolated(void);
	void CalculateConductanceForAbs(void);

	FloatType temp_GetOilSaturation(void);
	void temp_SetOilSaturation(FloatType);
	void temp_SetIsOilFilled(bool);
	FloatType temp_GetVolume(void);	
};
#endif