#include "NetworkElement.h"
#include <math.h>
#include <time.h>
#include <cstdlib>
#include <iostream>

void NetworkElement::CalcProps(void){
	FloatType Beta2Min, Beta2Max;
	FloatType tanBeta2;
	register unsigned int i, j;
	FloatType xx, ff, fp, cotx;
	
	//FloatType La;
	xx = 5;
	i = 0;
	do {
		cotx = cos(PI / xx) / sin(PI / xx);
		fp = PI*(cotx*cotx + 1) / (4 * xx*xx*xx) - cotx / (4 * xx*xx);
		ff = cotx / (4 * xx) - ShapeFactor;
		xx -= ff / fp;
		if (((xx - xx) != 0) || (xx <= 3.5)) xx = 4;
		i++;
		if (i > 50) break;
	} while (fabs(ff / fp) > SHAPETOL);
	if ((i > 50) || (xx > 20)) {		
		NumberOfSides = 0;		
	}
	else {
		NumberOfSides = (unsigned int)(xx + 0.5);		
	}

	SurfaceArea = InscribedRadius*InscribedRadius / (4 * ShapeFactor);
	if ((ShapeFactor <= (SQRT3 / 36)) || (NumberOfSides == 3)) {
		//IsCircle=false;
		NumberOfSides = 3;
		//Angles = PI / 3;

		tanBeta2=InscribedRadius*InscribedRadius/(4*ShapeFactor);

		//Ali Zamani's Conductance
		//La = 4 * sqrt(WaterViscosity*9.17e-11 / (InscribedRadius*InscribedRadius*InscribedRadius));
		//AbsoluteConductance = PI*InscribedRadius*InscribedRadius*InscribedRadius*InscribedRadius*La*Length*(2 + exp(La*Length) + exp(-La*Length)) / (16 * (exp(La*Length) + exp(-La*Length)));
		AbsoluteConductance=0.6*tanBeta2*tanBeta2*ShapeFactor;


		Beta2Min=atan((2.0/SQRT3)*cos((4.0*PI+acos(-12.0*SQRT3*ShapeFactor))/3.0));
		Beta2Max=atan((2.0/SQRT3)*cos(acos(-12.0*SQRT3*ShapeFactor)/3.0));
		Beta[1]=Beta2Min+(((FloatType) rand())/RAND_MAX)*(Beta2Max-Beta2Min);
		

		tanBeta2=tan(Beta[1]);
		Beta[0]=0.5*(-Beta[1]+asin(sin(Beta[1])*(tanBeta2+4.0*ShapeFactor)/(tanBeta2-4.0*ShapeFactor)));
		Beta[2]=0.5*PI-Beta[0]-Beta[1];	

		for (i=0; i<2; i++) {				//Remove Candidate
			for (j=i+1; j<3; j++) {
				if (Beta[i]>Beta[j]) {					
					Beta2Max=Beta[i];
					Beta[i]=Beta[j];
					Beta[j]=Beta2Max;
				}
			}
		}
	}
	else if ((ShapeFactor >= (1.0 / (4 * PI))) || (!NumberOfSides)) {
		//IsCircle = true;
		NumberOfSides = 0;
		//Angles = PI;

		tanBeta2 = PI*InscribedRadius*InscribedRadius;

		//Ali Zamani's Conductance
		//La = 4 * sqrt(WaterViscosity*9.17e-11 / (InscribedRadius*InscribedRadius*InscribedRadius));
		//AbsoluteConductance = PI*InscribedRadius*InscribedRadius*InscribedRadius*InscribedRadius*La*Length*(2 + exp(La*Length) + exp(-La*Length)) / (16 * (exp(La*Length) + exp(-La*Length)));
		AbsoluteConductance = 0.5*tanBeta2*tanBeta2*ShapeFactor;
	}
	else {		
		Beta[0] = (PI / (2.0*NumberOfSides))*(NumberOfSides - 2);
		ff = 14.18 + 0.5 / NumberOfSides - 26.4 / (NumberOfSides*NumberOfSides) + 102.18 / (NumberOfSides*NumberOfSides*NumberOfSides);
		AbsoluteConductance = 2 * SurfaceArea*SurfaceArea*sqrt(ShapeFactor) / ff;
	}
}
bool NetworkElement::GetIsOilFilled(void) {
	return IsOilFilled;
}

FloatType NetworkElement::GetPc(void) {
	return EntryCapillaryPressure;
}

void NetworkElement::SetIndex(unsigned int ElementIndex) {
	Index=ElementIndex;
}

void NetworkElement::FillWithOil(void) {
	IsOilFilled=true;
	IsWaterFilled=false;
	if (!NumberOfSides) WaterCoatingExist=false;
}

NetworkElement::NetworkElement(void) {
	IsOilFilled=false;	
	IsWaterFilled=true;

	IsConnectedToOutletbyWater=true;	
	OilSaturation=0;
	IsSnapOff=false;
	IsOilLayerExist[0]=false;
	IsOilLayerExist[1]=false;
	IsOilLayerExist[2]=false;
	NegativePCFilledPistonLike=false;
	WaterCoatingExist=true;

	HingeAngle[0]=RecedingContactAngle;
	HingeAngle[1]=RecedingContactAngle;
	HingeAngle[2]=RecedingContactAngle;
}

void NetworkElement::SetIsConnectedToOutlet(FluidType Fluid, bool ConnectionStatus) {
	if (Fluid == WATER) IsConnectedToOutletbyWater = ConnectionStatus;
	else IsConnectedToOutletbyOil = ConnectionStatus;
}
bool NetworkElement::GetIsConnectedToOutlet(FluidType Fluid) {
	if (Fluid == WATER) return IsConnectedToOutletbyWater;
	else return IsConnectedToOutletbyOil;	
}

void NetworkElement::SetIsConnectedToInlet(bool ConnectionStatus) {
	IsConnectedToInlet = ConnectionStatus;
}
bool NetworkElement::GetIsConnectedToInlet(void) {
	return IsConnectedToInlet;
}

void NetworkElement::SetIsVisited(bool VisitStatus) {
	IsVisited=VisitStatus;
}
bool NetworkElement::GetIsVisited(void) {	
	return IsVisited;
}

void NetworkElement::SetIsWaterFilled(bool WaterFilledStatus) {
	IsWaterFilled=WaterFilledStatus;
}
bool NetworkElement::GetIsWaterFilled(void) {
	return IsWaterFilled;
}

void NetworkElement::CalculateMaxAdvancingContactAngle(void) {			//Remove Candidate
	FloatType S1;	
	register unsigned int i;

	if (NumberOfSides==3) {
		S1=0;
		for (i=0; i<3; i++) {
			if (Beta[i]<(PI/2-RecedingContactAngle)) {
				S1+=cos(RecedingContactAngle+Beta[i]);
			}
		}
		MaxAdvancingContactAngle=acos(-4*ShapeFactor*S1/(InscribedRadius/Rmin-cos(RecedingContactAngle)+12*ShapeFactor*sin(RecedingContactAngle)));	
	}
	else if (NumberOfSides > 3) {		
		S1 = 0;
		if (Beta[0] < (PI / 2 - RecedingContactAngle)) {
			S1 = NumberOfSides*cos(RecedingContactAngle + Beta[0]);
		}	
		MaxAdvancingContactAngle = acos(-4 * ShapeFactor*S1 / (InscribedRadius / Rmin - cos(RecedingContactAngle) + 12 * ShapeFactor*sin(RecedingContactAngle)));
	}
}
void NetworkElement::FillWithWater(FloatType WorkingPc) {
	register unsigned int i;

	IsWaterFilled=true;
	WaterCoatingExist=true;
	if (IsSnapOff || (WorkingPc>=0) || (!NumberOfSides)) IsOilFilled=false;
	//else if ((!IsSnapOff) && (WorkingPc<0) && (NumberOfSides)) {
	else if ((!IsSnapOff) && (WorkingPc<0)) {
		if (NumberOfSides == 3) {
			for (i = 0; i<3; i++) if (AdvancingContactAngle>(PI / 2 + Beta[i])) {
				IsOilFilled = true;
				IsOilLayerExist[i] = true;
			}
		}
		else if (NumberOfSides > 3) {
			if (AdvancingContactAngle>(PI / 2 + Beta[0])) {
				IsOilFilled = true;
				IsOilLayerExist[0] = true;
			}
		}
		
	}
}
void NetworkElement::OilLayerExist(FloatType WorkingPc) {
	register unsigned int i, n;
	FloatType tPc;

	FloatType S1, S2, S3, S4, S5, S6, S7, S8;
	FloatType costr, A, Aeff, Rold;
	FloatType HingeAngle[3], Alpha[3];
	bool MinCond, NRCond;
	FloatType dHingedR[3], dbdR[3], dAlphadR[3];
	FloatType AA, BB, CC, DD, dAAdR, dBBdR, dCCdR, dDDdR, FF, dFFdR;


	if (NegativePCFilledPistonLike) {
		if (NumberOfSides == 3) {
			R = Rmin;
			n = 0;
			NRCond = false;
			do {
				S1 = 0;
				S2 = 0;
				S3 = 0;
				S4 = 0;
				S5 = 0;
				S6 = 0;
				S7 = 0;
				S8 = 0;

				for (i = 0; i < 3; i++) {
					if (Beta[i] < (PI / 2 - RecedingContactAngle)) {
						costr = acos(Rmin*cos(RecedingContactAngle + Beta[i]) / R) - Beta[i];
						MinCond = (costr < AdvancingContactAngle);
						HingeAngle[i] = MinCond ? costr : AdvancingContactAngle;
						if (HingeAngle[i] < RecedingContactAngle) HingeAngle[i] = RecedingContactAngle;
						costr = Rmin*cos(RecedingContactAngle + Beta[i]) / R;
						dHingedR[i] = MinCond ? cos(RecedingContactAngle + Beta[i])*Rmin / (R*R*sqrt(1 - costr*costr)) : 0;
						if (HingeAngle[i] < RecedingContactAngle) dHingedR[i] = 0;

						Hinge_b[i] = (HingeAngle[i] > AdvancingContactAngle) ? R*cos(AdvancingContactAngle + Beta[i]) / sin(Beta[i]) : Rmin*cos(RecedingContactAngle + Beta[i]) / sin(Beta[i]);
						dbdR[i] = (HingeAngle[i] > AdvancingContactAngle) ? cos(AdvancingContactAngle + Beta[i]) / sin(Beta[i]) : 0;
						if (Beta[i] >= (PI / 2 - RecedingContactAngle)) {
							Hinge_b[i] = 0;
						}

						costr = Hinge_b[i] * sin(Beta[i]) / R;
						Alpha[i] = (HingeAngle[i] > AdvancingContactAngle) ? (PI / 2 - AdvancingContactAngle - Beta[i]) : asin(costr);
						dAlphadR[i] = (HingeAngle[i] > AdvancingContactAngle) ? 0 : sin(Beta[i])*(dbdR[i] * R - Hinge_b[i]) / (R*R*sqrt(1 - costr*costr));


						S1 += Hinge_b[i] * cos(HingeAngle[i]);
						S2 += dbdR[i] * cos(HingeAngle[i]) - Hinge_b[i] * dHingedR[i] * sin(HingeAngle[i]);
						S3 += PI / 2 - HingeAngle[i] - Beta[i];
						S4 += dHingedR[i];
						S5 += Alpha[i];
						S6 += dAlphadR[i];
						S7 += dbdR[i];
						S8 += Hinge_b[i];
					}
				}
				A = InscribedRadius*InscribedRadius / (4 * ShapeFactor);
				AA = -R*S1;
				BB = R*R*S3;
				CC = 2 * R*S5;
				DD = (InscribedRadius / (2 * ShapeFactor) - 2 * S8)*cos(AdvancingContactAngle);
				FF = -R + (A + AA + BB) / (CC + DD);

				dAAdR = -(S1 + R*S2);
				dBBdR = 2 * R*S3 - R*R*S4;
				dCCdR = 2 * S5 + 2 * R*S6;
				dDDdR = -2 * S7*cos(AdvancingContactAngle);
				dFFdR = -1 + ((dAAdR + dBBdR)*(CC + DD) - (dCCdR + dDDdR)*(A + AA + BB)) / ((CC + DD)*(CC + DD));

				Rold = R;
				R = R - FF / dFFdR;
				Aeff = A - R*S1 + R*R*S3;

				n++;
				if (n > MAXNITERS) {
					NRCond = true;
					break;
				}
			} while (fabs(R - Rold) > fabs(R*RTOL));

			if (NRCond) {
				R = Rmin;
				n = 0;
				do {
					S1 = 0;
					S2 = 0;
					S3 = 0;
					S4 = 0;
					for (i = 0; i < 3; i++) {
						if (Beta[i] < (PI / 2 - RecedingContactAngle)) {
							costr = acos(Rmin*cos(RecedingContactAngle + Beta[i]) / R) - Beta[i];
							HingeAngle[i] = (costr < AdvancingContactAngle) ? costr : AdvancingContactAngle;
							if (HingeAngle[i] < RecedingContactAngle) HingeAngle[i] = RecedingContactAngle;

							Hinge_b[i] = (HingeAngle[i] < AdvancingContactAngle) ? Rmin*cos(RecedingContactAngle + Beta[i]) / sin(Beta[i]) : R*cos(AdvancingContactAngle + Beta[i]) / sin(Beta[i]);
							if (Beta[i] >= (PI / 2 - RecedingContactAngle)) {
								Hinge_b[i] = 0;
							}
							Alpha[i] = (HingeAngle[i] > AdvancingContactAngle) ? PI / 2 - AdvancingContactAngle - Beta[i] : asin(Hinge_b[i] * sin(Beta[i]) / R);
							S1 += Hinge_b[i] * cos(HingeAngle[i]);
							S2 += PI / 2 - HingeAngle[i] - Beta[i];
							S3 += Alpha[i];
							S4 += Hinge_b[i];
						}
					}
					A = InscribedRadius*InscribedRadius / (4 * ShapeFactor);
					Aeff = A - R*S1 + R*R*S2;
					Rold = R;
					R = Aeff / (2 * R*S3 + (InscribedRadius / (2 * ShapeFactor) - 2 * S4)*cos(AdvancingContactAngle));

					n++;
					if (n > MAXNITERS) {
						TerM("Divergence Error!\n");
					}
				} while (fabs(R - Rold) > (R*RTOL));
			}
			OilSaturation = Aeff / A;
		}
		else if (NumberOfSides > 3) {
			R = Rmin;
			n = 0;
			NRCond = false;
			do {
				S1 = 0;
				S2 = 0;
				S3 = 0;
				S4 = 0;
				S5 = 0;
				S6 = 0;
				S7 = 0;
				S8 = 0;

				if (Beta[0] < (PI / 2 - RecedingContactAngle)) {
					costr = acos(Rmin*cos(RecedingContactAngle + Beta[0]) / R) - Beta[0];
					MinCond = (costr < AdvancingContactAngle);
					HingeAngle[0] = MinCond ? costr : AdvancingContactAngle;
					if (HingeAngle[0] < RecedingContactAngle) HingeAngle[0] = RecedingContactAngle;
					costr = Rmin*cos(RecedingContactAngle + Beta[0]) / R;
					dHingedR[0] = MinCond ? cos(RecedingContactAngle + Beta[0])*Rmin / (R*R*sqrt(1 - costr*costr)) : 0;
					if (HingeAngle[0] < RecedingContactAngle) dHingedR[0] = 0;

					Hinge_b[0] = (HingeAngle[0] > AdvancingContactAngle) ? R*cos(AdvancingContactAngle + Beta[0]) / sin(Beta[0]) : Rmin*cos(RecedingContactAngle + Beta[0]) / sin(Beta[0]);
					dbdR[0] = (HingeAngle[0] > AdvancingContactAngle) ? cos(AdvancingContactAngle + Beta[0]) / sin(Beta[0]) : 0;
					if (Beta[0] >= (PI / 2 - RecedingContactAngle)) {
						Hinge_b[0] = 0;
					}

					costr = Hinge_b[0] * sin(Beta[0]) / R;
					Alpha[0] = (HingeAngle[0] > AdvancingContactAngle) ? (PI / 2 - AdvancingContactAngle - Beta[0]) : asin(costr);
					dAlphadR[0] = (HingeAngle[0] > AdvancingContactAngle) ? 0 : sin(Beta[0])*(dbdR[0] * R - Hinge_b[0]) / (R*R*sqrt(1 - costr*costr));


					S1 = NumberOfSides*(Hinge_b[0] * cos(HingeAngle[0]));
					S2 = NumberOfSides*(dbdR[0] * cos(HingeAngle[0]) - Hinge_b[0] * dHingedR[0] * sin(HingeAngle[0]));
					S3 = NumberOfSides*(PI / 2 - HingeAngle[0] - Beta[0]);
					S4 = NumberOfSides*dHingedR[0];
					S5 = NumberOfSides*Alpha[0];
					S6 = NumberOfSides*dAlphadR[0];
					S7 = NumberOfSides*dbdR[0];
					S8 = NumberOfSides*Hinge_b[0];
				}

				A = InscribedRadius*InscribedRadius / (4 * ShapeFactor);
				AA = -R*S1;
				BB = R*R*S3;
				CC = 2 * R*S5;
				DD = (InscribedRadius / (2 * ShapeFactor) - 2 * S8)*cos(AdvancingContactAngle);
				FF = -R + (A + AA + BB) / (CC + DD);

				dAAdR = -(S1 + R*S2);
				dBBdR = 2 * R*S3 - R*R*S4;
				dCCdR = 2 * S5 + 2 * R*S6;
				dDDdR = -2 * S7*cos(AdvancingContactAngle);
				dFFdR = -1 + ((dAAdR + dBBdR)*(CC + DD) - (dCCdR + dDDdR)*(A + AA + BB)) / ((CC + DD)*(CC + DD));

				Rold = R;
				R = R - FF / dFFdR;
				Aeff = A - R*S1 + R*R*S3;

				n++;
				if (n > MAXNITERS) {
					NRCond = true;
					break;
				}
			} while (fabs(R - Rold) > fabs(R*RTOL));

			if (NRCond) {
				R = Rmin;
				n = 0;
				do {
					S1 = 0;
					S2 = 0;
					S3 = 0;
					S4 = 0;
					if (Beta[0] < (PI / 2 - RecedingContactAngle)) {
						costr = acos(Rmin*cos(RecedingContactAngle + Beta[0]) / R) - Beta[0];
						HingeAngle[0] = (costr < AdvancingContactAngle) ? costr : AdvancingContactAngle;
						if (HingeAngle[0] < RecedingContactAngle) HingeAngle[0] = RecedingContactAngle;

						Hinge_b[0] = (HingeAngle[0] < AdvancingContactAngle) ? Rmin*cos(RecedingContactAngle + Beta[0]) / sin(Beta[0]) : R*cos(AdvancingContactAngle + Beta[0]) / sin(Beta[0]);
						if (Beta[0] >= (PI / 2 - RecedingContactAngle)) {
							Hinge_b[0] = 0;
						}
						Alpha[0] = (HingeAngle[0] > AdvancingContactAngle) ? PI / 2 - AdvancingContactAngle - Beta[0] : asin(Hinge_b[0] * sin(Beta[0]) / R);
						S1 = NumberOfSides*Hinge_b[0] * cos(HingeAngle[0]);
						S2 = NumberOfSides*(PI / 2 - HingeAngle[0] - Beta[0]);
						S3 = NumberOfSides*Alpha[0];
						S4 = NumberOfSides*Hinge_b[0];
					}
					A = InscribedRadius*InscribedRadius / (4 * ShapeFactor);
					Aeff = A - R*S1 + R*R*S2;
					Rold = R;
					R = Aeff / (2 * R*S3 + (InscribedRadius / (2 * ShapeFactor) - 2 * S4)*cos(AdvancingContactAngle));

					n++;
					if (n > MAXNITERS) {
						TerM("Divergence Error!\n");
					}
				} while (fabs(R - Rold) > (R*RTOL));
			}
			OilSaturation = Aeff / A;
		}




		if (NumberOfSides == 3) {
			for (i = 0; i < 3; i++) {
				if (IsOilLayerExist[i]) {
					tPc = IFT*cos(Beta[i] + acos(2 * sin(Beta[i]) + cos(AdvancingContactAngle))) / (Hinge_b[i] * sin(Beta[i]));
					if (WorkingPc <= tPc) IsOilLayerExist[i] = false;
				}
			}
			if ((!IsOilLayerExist[0]) && (!IsOilLayerExist[1]) && (!IsOilLayerExist[2])) {
				IsOilFilled = false;
				NegativePCFilledPistonLike = false;
				OilSaturation = 0;
			}
		}
		else if (NumberOfSides > 3) {
			if (IsOilLayerExist[0]) {
				tPc = IFT*cos(Beta[0] + acos(2 * sin(Beta[0]) + cos(AdvancingContactAngle))) / (Hinge_b[0] * sin(Beta[0]));
				if (WorkingPc <= tPc) IsOilLayerExist[0] = false;
			}
			if (!IsOilLayerExist[0]) {
				IsOilFilled = false;
				NegativePCFilledPistonLike = false;
				OilSaturation = 0;
			}
		}


		/*
		if (NumberOfSides==3) {
			R=Rmin;
			do {
				S1=0;
				S2=0;
				S3=0;
				S4=0;
				for (i=0; i<3; i++) {
					if (IsOilLayerExist[i]) {
						costr=acos(Rmin*cos(RecedingContactAngle+Beta[i])/R)-Beta[i];
						HingeAngle[i]=(costr<AdvancingContactAngle) ? costr : AdvancingContactAngle;
						if (HingeAngle[i]<RecedingContactAngle) HingeAngle[i]=RecedingContactAngle;

						Hinge_b[i]=(HingeAngle[i]<AdvancingContactAngle) ? Rmin*cos(RecedingContactAngle+Beta[i])/sin(Beta[i]) : R*cos(AdvancingContactAngle+Beta[i])/sin(Beta[i]);
						Alpha[i]=(HingeAngle[i]>AdvancingContactAngle) ? PI/2-AdvancingContactAngle-Beta[i] : asin(Hinge_b[i]*sin(Beta[i])/R);

						S1+=Hinge_b[i]*cos(HingeAngle[i]);
						S2+=PI/2-HingeAngle[i]-Beta[i];
						S3+=Alpha[i];
						S4+=Hinge_b[i];
					}
				}
				A=InscribedRadius*InscribedRadius/(4*ShapeFactor);
				Aeff=A-R*S1+R*R*S2;
				Rold=R;
				R=Aeff/(2*R*S3+(InscribedRadius/(2*ShapeFactor)-2*S4)*cos(AdvancingContactAngle));
			} while (fabs(R-Rold)>(R*RTOL));
			OilSaturation=Aeff/A;
		}
		else if (NumberOfSides > 3) {
			R = Rmin;
			S1=0;
			S2=0;
			S3=0;
			S4=0;
			do {
				if (IsOilLayerExist[0]) {
					costr = acos(Rmin*cos(RecedingContactAngle + Beta[0]) / R) - Beta[0];
					HingeAngle[0] = (costr < AdvancingContactAngle) ? costr : AdvancingContactAngle;
					if (HingeAngle[0]<RecedingContactAngle) HingeAngle[0] = RecedingContactAngle;

					Hinge_b[0] = (HingeAngle[0] < AdvancingContactAngle) ? Rmin*cos(RecedingContactAngle + Beta[0]) / sin(Beta[0]) : R*cos(AdvancingContactAngle + Beta[0]) / sin(Beta[0]);
					Alpha[0] = (HingeAngle[0] > AdvancingContactAngle) ? PI / 2 - AdvancingContactAngle - Beta[0] : asin(Hinge_b[0] * sin(Beta[0]) / R);

					S1 = NumberOfSides*(Hinge_b[0] * cos(HingeAngle[0]));
					S2 = NumberOfSides*(PI / 2 - HingeAngle[0] - Beta[0]);
					S3 = NumberOfSides*Alpha[0];
					S4 = NumberOfSides*Hinge_b[0];
				}
				A = InscribedRadius*InscribedRadius / (4 * ShapeFactor);
				Aeff = A - R*S1 + R*R*S2;
				Rold = R;
				R = Aeff / (2 * R*S3 + (InscribedRadius / (2 * ShapeFactor) - 2 * S4)*cos(AdvancingContactAngle));
			} while (fabs(R - Rold)>(R*RTOL));
			OilSaturation = Aeff / A;
		}*/
	}
}


FloatType NetworkElement::GetTotalWaterSaturation(void) {
	return ((1-OilSaturation)*Volume+ClayVolume)/TotalVolume;
}

void NetworkElement::CalculateDrainagePc(void) {
	register unsigned int i;
	FloatType costr;

	
	//Primary Oil Flooding
	
	costr=cos(RecedingContactAngle);
	if (!NumberOfSides) {
		Fd=1;		
	}
	else {
		SumForR = 0;
		if (NumberOfSides == 3) {			
			for (i = 0; i < 3; i++) {
				if (Beta[i] < (PI / 2 - RecedingContactAngle)) {
					SumForR += costr*cos(RecedingContactAngle + Beta[i]) / sin(Beta[i]) + RecedingContactAngle + Beta[i] - PI / 2;
				}
			}
		}
		else if (NumberOfSides > 3) {
			SumForR=0;
			if (Beta[0] < (PI / 2 - RecedingContactAngle)) {
				SumForR = NumberOfSides*(costr*cos(RecedingContactAngle + Beta[0]) / sin(Beta[0]) + RecedingContactAngle + Beta[0] - PI / 2);
			}
		}
		Fd=(1+sqrt(1-4*ShapeFactor*SumForR/(costr*costr)))/(1+2*sqrt(PI*ShapeFactor));		
	}	
	EntryCapillaryPressure=IFT*costr*(1+2*sqrt(PI*ShapeFactor))*Fd/InscribedRadius-DeltaRho*G_ACC*Z;
}

void NetworkElement::CalculateOilSaturation(FloatType WorkingPc) {
	register unsigned int i; 
	//Primary Oil Flooding

	if (IsConnectedToOutletbyWater) {
		if (IsOilFilled) {
			if (!NumberOfSides) {		
				OilSaturation=1;
			}
			else {
				if (WorkingPc){
					R=IFT/WorkingPc;
					OilSaturation=1-4*ShapeFactor*R*R*SumForR/(InscribedRadius*InscribedRadius);
					if (NumberOfSides == 3) for (i = 0; i < 3; i++) Hinge_b[i] = R*cos(RecedingContactAngle + Beta[i]) / sin(Beta[i]);
					else if (NumberOfSides > 3) Hinge_b[0] = R*cos(RecedingContactAngle + Beta[0]) / sin(Beta[0]);
				}
				else {
					OilSaturation=0;
					if (NumberOfSides == 3) for (i = 0; i < 3; i++) Hinge_b[i] = 0;
					else if (NumberOfSides > 3) Hinge_b[0] = 0;					
				}			
			}	
		}
		else {
			OilSaturation=0;
		}	
	}	
}


void NetworkElement::CalculateConductance(void) {
	FloatType Ac, Gc, Gs, C, TempVar;
	register unsigned int i;

	if (IsWaterFilled && (!IsOilFilled)){				
		WaterConductance=AbsoluteConductance/WaterViscosity;
		OilConductance=0;
	}
	else if (((!IsWaterFilled) && IsOilFilled && (!WaterCoatingExist)) || ((!NumberOfSides) && IsOilFilled)){					
		WaterConductance=0;
		OilConductance=AbsoluteConductance/OilViscosity;
	}
	else if ((IsWaterFilled && IsOilFilled) || (IsOilFilled && WaterCoatingExist)) {				
		//WaterConductance = 0.5* AbsoluteConductance*(1 - OilSaturation)*(1 - OilSaturation) / WaterViscosity;
		WaterConductance = 0;
		if (NumberOfSides == 3) {
			for (i = 0; i < 3; i++) {
				if (Beta[i] < (PI / 2 - RecedingContactAngle)) {
					TempVar = (Hinge_b[i] * sin(Beta[i]) / cos(HingeAngle[i] + Beta[i]));
					Ac = TempVar*TempVar*(cos(HingeAngle[i])*cos(HingeAngle[i] + Beta[i]) / sin(Beta[i]) + HingeAngle[i] + Beta[i] - PI / 2);

					TempVar = 2 * Hinge_b[i] * (1 - sin(Beta[i])*(HingeAngle[i] + Beta[i] - PI / 2) / cos(HingeAngle[i] + Beta[i]));
					Gc = Ac / (TempVar*TempVar);

					TempVar = 1 + sin(Beta[i]);
					Gs = sin(Beta[i])*cos(Beta[i]) / (4 * TempVar*TempVar);

					C = 0.364 + 0.28*Gs / Gc;

					WaterConductance += C*Ac*Ac*Gc / WaterViscosity;
				}
			}
		}
		else if (NumberOfSides > 3) {
			if (Beta[0] < (PI / 2 - RecedingContactAngle)) {
				TempVar = (Hinge_b[0] * sin(Beta[0]) / cos(HingeAngle[0] + Beta[0]));
				Ac = TempVar*TempVar*(cos(HingeAngle[0])*cos(HingeAngle[0] + Beta[0]) / sin(Beta[0]) + HingeAngle[0] + Beta[0] - PI / 2);

				TempVar = 2 * Hinge_b[0] * (1 - sin(Beta[0])*(HingeAngle[0] + Beta[0] - PI / 2) / cos(HingeAngle[0] + Beta[0]));
				Gc = Ac / (TempVar*TempVar);

				TempVar = 1 + sin(Beta[0]);
				Gs = sin(Beta[0])*cos(Beta[0]) / (4 * TempVar*TempVar);

				C = 0.364 + 0.28*Gs / Gc;

				WaterConductance = NumberOfSides*(C*Ac*Ac*Gc / WaterViscosity);
			}
		}
		OilConductance=AbsoluteConductance*OilSaturation*OilSaturation/OilViscosity;
	}
}

FloatType NetworkElement::GetLength(void){
	return Length;
}

FloatType NetworkElement::GetAbsoluteConductance(void) {
	return AbsoluteConductance;
}

FloatType NetworkElement::GetWaterConductance(void) {
	return WaterConductance;
}

FloatType NetworkElement::GetOilConductance(void) {
	return OilConductance;
}

bool NetworkElement::GetWaterCoatingExist (void) {
	return WaterCoatingExist;
}

FloatType NetworkElement::GetZ(void) {
	return Z;
}
unsigned int NetworkElement::GetIndex(void) {
	return Index;
}

FloatType NetworkElement::GetReducedVolume(void) {
	return ClayVolume + Volume;
}

bool NetworkElement::IsIsolated(void) {
	return !(IsConnectedToInlet && IsConnectedToOutletbyWater);
}

void NetworkElement::CalculateConductanceForAbs(void) {
	WaterConductance = AbsoluteConductance / WaterViscosity;
	OilConductance = 0;
}



FloatType NetworkElement::temp_GetOilSaturation(void){
	return OilSaturation;
}
void NetworkElement::temp_SetOilSaturation(FloatType InSat){
	OilSaturation=InSat;
}
void NetworkElement::temp_SetIsOilFilled(bool inbool) {
	IsOilFilled=inbool;
}
FloatType NetworkElement::temp_GetVolume(void) {
	return (Volume-ClayVolume)/TotalVolume;
}

