#include <math.h>
#include "Pore.h"
#include "ElementList.h"
#include "MIfstream.h"


Pore::Pore(void) {

}

Pore::~Pore(void) {
	delete[] AdjacentPores;	
	delete[] ConnectingThroats;		
}

void Pore::SetLength(FloatType InLength) {
	Length=InLength;
}

unsigned int Pore::ReadNode1(MIfstream& InputFile){
	char str[MAX_STRING_LENGTH];
	register int j;
	register unsigned int i, Nulls;

	Nulls=0;

	if (!InputFile.ReadWord(str)) TerM("Incorrect node1 file format!");
	if (!InputFile.ReadWord(str)) TerM("Incorrect node1 file format!");
	X=atof(str);

	if (!InputFile.ReadWord(str)) TerM("Incorrect node1 file format!");
	Y=atof(str);

	if (!InputFile.ReadWord(str)) TerM("Incorrect node1 file format!");
	Z=atof(str);

	if (!InputFile.ReadWord(str)) TerM("Incorrect node1 file format!");
	CoordinationNumber=atoi(str);

	AdjacentPores=new Pore* [CoordinationNumber];
	ConnectingThroats=new Throat* [CoordinationNumber];

	for (i=0; i<CoordinationNumber; i++) {
		if (!InputFile.ReadWord(str)) TerM("Incorrect node1 file format!");
		j=atoi(str);
		if (j>0) AdjacentPores[i]=&pores[j-1];
		else {
			AdjacentPores[i]=NULL;
			Nulls++;
		}
	}

	if (!InputFile.ReadWord(str)) TerM("Incorrect node1 file format!");
	i=atoi(str);
	if (!InputFile.ReadWord(str)) TerM("Incorrect node1 file format!");
	j=atoi(str);
	if (i==1) IOStat=-1;
	else if (j==1) IOStat=0;
	else IOStat=1;

	for (i=0; i<CoordinationNumber; i++) {
		if (!InputFile.ReadWord(str)) TerM("Incorrect node1 file format!");
		j=atoi(str);
		if (j>0) ConnectingThroats[i]=&throats[j-1];
		else ConnectingThroats[i]=NULL;		
	}

	return CoordinationNumber-Nulls;
}
void Pore::ReadNode2(MIfstream& InputFile){
	char str[MAX_STRING_LENGTH];

	if (!InputFile.ReadWord(str)) TerM("Incorrect node2 file format!");
	if (!InputFile.ReadWord(str)) TerM("Incorrect node2 file format!");
	Volume=atof(str);
	TotalVolume+=Volume;

	if (!InputFile.ReadWord(str)) TerM("Incorrect node2 file format!");
	InscribedRadius=atof(str);

	if (!InputFile.ReadWord(str)) TerM("Incorrect node2 file format!");
	ShapeFactor=atof(str);

	if (!InputFile.ReadWord(str)) TerM("Incorrect node2 file format!");
	ClayVolume=atof(str);

	TotalVolume+=ClayVolume;
}

unsigned int Pore::CalculateImbibitionPc(FloatType WorkingPc) {
	FloatType S1, S2, S3, S4, S5, S6, S7, S8;
	register unsigned int i, n;
	FloatType costr, sqrtr;
	FloatType tPc1, tPc2;
	FloatType Alpha[3];
	FloatType Aeff, A, Rold;
	bool IsReadyToFill;
	FloatType PcWithZ1, PcWithZ2;
	unsigned int ReturnVal = 0;
	bool MinCond, NRCond;
	FloatType dHingedR[3], dbdR[3], dAlphadR[3];
	FloatType AA, BB, CC, DD, dAAdR, dBBdR, dCCdR, dDDdR, FF, dFFdR;


	if (IsConnectedToOutletbyWater) {
		if (WorkingPc>=0) {					//Secondary Water Flooding- Spontaneous Imbibition
			n=0;			
			IsReadyToFill=false;
			for (i=0; i<CoordinationNumber; i++) {
				if ((ConnectingThroats[i]!=NULL) && (ConnectingThroats[i]->GetIsWaterFilled())) {					
					IsReadyToFill=true;
				}
				if ((ConnectingThroats[i]!=NULL) && (ConnectingThroats[i]->GetIsOilFilled())) n++;				
			}
			
			//Piston Like
			S1=0;
			sqrtr=0.03/sqrt(AbsPerm*9.869233e-16);
			for (i=2; i<=n; i++) {
				S1+=(((FloatType) rand())/RAND_MAX)*sqrtr;
			}
			tPc1=(2*IFT*cos(AdvancingContactAngle)/InscribedRadius)-(IFT*S1);

			//Snap Off
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

				
				
				
				SumForR=0;					
				for (i=0; i<3; i++) {
					if (Beta[i]<(PI/2-RecedingContactAngle)) {
						SumForR+=cos(HingeAngle[i])*cos(HingeAngle[i]+Beta[i])/sin(Beta[i])+HingeAngle[i]+Beta[i]-PI/2;				
					}
				}			
				


				n=0;
				for (i=0; i<3; i++) {
					if ((AdvancingContactAngle<(PI/2-Beta[i])) && (HingeAngle[i]>=AdvancingContactAngle)) n++;
				}
				if (n==0) {
					tPc2=-CAPILLARYLIMIT - 1;
				}
				if (n==1) {
					tPc2=(IFT/InscribedRadius)*((cos(AdvancingContactAngle)/tan(Beta[0])-sin(AdvancingContactAngle)+cos(HingeAngle[2])/tan(Beta[2])-sin(HingeAngle[2]))/(1.0/tan(Beta[0])+1/tan(Beta[1])));
				}
				else {
					tPc2=(IFT/InscribedRadius)*(cos(AdvancingContactAngle)-2*sin(AdvancingContactAngle)/(1.0/tan(Beta[0])+1/tan(Beta[1])));
				}				
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

				SumForR = 0;
				if (Beta[0] < (PI / 2 - RecedingContactAngle)) {
					SumForR = NumberOfSides*(cos(HingeAngle[0])*cos(HingeAngle[0] + Beta[0]) / sin(Beta[0]) + HingeAngle[0] + Beta[0] - PI / 2);
				}

				
				n = 0;
				if ((AdvancingContactAngle < (PI / 2 - Beta[0])) && (HingeAngle[0] >= AdvancingContactAngle)) n = NumberOfSides;
				
				
				if (n == 0) {
					tPc2 = -CAPILLARYLIMIT - 1;
				}
				else {
					tPc2 = (IFT / InscribedRadius)*(cos(AdvancingContactAngle) - 2 * sin(AdvancingContactAngle) / (1.0 / tan(Beta[0]) + 1 / tan(Beta[0])));
				}
			}
			else tPc2 = -CAPILLARYLIMIT - 1;

			PcWithZ1=tPc1-DeltaRho*G_ACC*Z;
			PcWithZ2=tPc2-DeltaRho*G_ACC*Z;
			if (((!NumberOfSides) || (tPc1>=tPc2)) && (IsReadyToFill) && (!IsWaterFilled) && (PcWithZ1>=WorkingPc)){
				EntryCapillaryPressure=PcWithZ1;				
				OilSaturation=0;
				//OilSaturation=tSat;
				FillWithWater(WorkingPc);
				ReturnVal=1;
				IsSnapOff=false;				
			}
			else if ((!IsReadyToFill) && (!IsWaterFilled) && (PcWithZ2>=WorkingPc)) {
				EntryCapillaryPressure=PcWithZ2;			
				OilSaturation=0;				
				FillWithWater(WorkingPc);
				ReturnVal=1;
				IsOilFilled=false;				
				IsSnapOff=true;				
			}
			else {
				ReturnVal=0;
				if (IsOilFilled) {
					if (!NumberOfSides) {		
						OilSaturation=1;
					}
					else {
						if (WorkingPc){
							R=IFT/WorkingPc;
							OilSaturation=1-4*ShapeFactor*R*R*SumForR/(InscribedRadius*InscribedRadius);
						}
						else {
							OilSaturation=0;
						}			
					}	
				}
				else {
					OilSaturation=0;
				}	
			}
		}
		else {					//Secondary Water Flooding- Forced Imbibition
			IsReadyToFill=false;
			for (i=0; i<CoordinationNumber; i++) {
				if ((ConnectingThroats[i]!=NULL) && (ConnectingThroats[i]->GetIsWaterFilled())) {					
					IsReadyToFill=true;
				}							
			}			

			//Piston Like
			costr=cos(AdvancingContactAngle);
			if (!NumberOfSides) {
				Fd=1;		
			}			
			else if (NumberOfSides==3){
				SumForR=0;		
				for (i=0; i<3; i++) {
					if (Beta[i]<(PI/2-AdvancingContactAngle)) {
						SumForR+=costr*cos(AdvancingContactAngle+Beta[i])/sin(Beta[i])+AdvancingContactAngle+Beta[i]-PI/2;				
					}
				}		
				Fd=(1+sqrt(1-4*ShapeFactor*SumForR/(costr*costr)))/(1+2*sqrt(PI*ShapeFactor));
			}
			else if (NumberOfSides>3) {
				SumForR=0;				
				if (Beta[0]<(PI/2-AdvancingContactAngle)) {
					SumForR=NumberOfSides*(costr*cos(AdvancingContactAngle+Beta[0])/sin(Beta[0])+AdvancingContactAngle+Beta[0]-PI/2);				
				}
				Fd=(1+sqrt(1-4*ShapeFactor*SumForR/(costr*costr)))/(1+2*sqrt(PI*ShapeFactor));		
			}
			tPc1=-IFT*costr*(1+2*sqrt(PI*ShapeFactor))*Fd/InscribedRadius;

			//Snap Off
			if (NumberOfSides==3) {
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
			}
			else if (NumberOfSides>3) {
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
			}
			else {
				tPc2 = -CAPILLARYLIMIT - 1;
			}

			if (NumberOfSides) {
				/*if ((AdvancingContactAngle<=(PI-Beta[0])) && (HingeAngle[0]>=AdvancingContactAngle)) {
					tPc2=IFT*cos(AdvancingContactAngle+Beta[0])/(Rmin*cos(RecedingContactAngle+Beta[0]));
				}
				else if ((AdvancingContactAngle>(PI-Beta[0])) && (HingeAngle[0]>=(AdvancingContactAngle-Beta[0]))) {
					tPc2=-IFT/(Rmin*cos(RecedingContactAngle+Beta[0]));
				}
				else {
					tPc2=-CAPILLARYLIMIT-1;				
				}*/

				if (AdvancingContactAngle<=(PI-Beta[0])) {
					tPc2=IFT*cos(AdvancingContactAngle+Beta[0])/(Rmin*cos(RecedingContactAngle+Beta[0]));
				}
				else {
					tPc2=-IFT/(Rmin*cos(RecedingContactAngle+Beta[0]));
				}
			}
			

			PcWithZ1=tPc1-DeltaRho*G_ACC*Z;
			PcWithZ2=tPc2-DeltaRho*G_ACC*Z;
			if (((!NumberOfSides) || (tPc1>=tPc2)) && (IsReadyToFill) && (!IsWaterFilled) && (PcWithZ1>=WorkingPc)){
				EntryCapillaryPressure=PcWithZ1;				
				OilSaturation=0;				
				FillWithWater(WorkingPc);
				ReturnVal=1;
				IsSnapOff=false;
				NegativePCFilledPistonLike=true;
			}
			else if ((!IsReadyToFill) && (!IsWaterFilled) && (PcWithZ2>=WorkingPc)) {
				EntryCapillaryPressure=PcWithZ2;				
				OilSaturation=0;				
				FillWithWater(WorkingPc);
				ReturnVal=1;
				IsOilFilled=false;
				IsSnapOff=true;				
			}
			else {
				ReturnVal=0;
				EntryCapillaryPressure=-CAPILLARYLIMIT-1;
				IsSnapOff=false;
				if (IsOilFilled) {
					if (!NumberOfSides) {		
						OilSaturation=1;
					}
					else {
						if (WorkingPc){
							R=IFT/WorkingPc;
							OilSaturation=1-4*ShapeFactor*R*R*SumForR/(InscribedRadius*InscribedRadius);
						}
						else {
							OilSaturation=0;
						}			
					}	
				}
				else {
					OilSaturation=0;
				}	
			}
		}
	}
	return ReturnVal;
}

unsigned int Pore::TestFillReadiness(FloatType WorkingPc) {
	bool InvestigateThroats=false;
	register unsigned int i;

	for (i=0; i<CoordinationNumber; i++) {
		InvestigateThroats=InvestigateThroats || ((ConnectingThroats[i]!=NULL) && (ConnectingThroats[i]->GetIsOilFilled()));
	}

	if ((!IsOilFilled) && (InvestigateThroats) && (WorkingPc>=EntryCapillaryPressure) && IsConnectedToOutletbyWater) {		
		FillWithOil();
		return 1;
	}
	else return 0;
}

void Pore::SweepAdjacentThroats(FluidType Fluid) {
	register unsigned int i;
	for (i=0; i<CoordinationNumber; i++) {
		if ((ConnectingThroats[i]!=NULL) && (!ConnectingThroats[i]->GetIsVisited()) && (((Fluid==OIL) || (ConnectingThroats[i]->GetIsOilFilled())) || ((Fluid==WATER) || (ConnectingThroats[i]->GetWaterCoatingExist())))) {
			ConnectingThroats[i]->SetIsVisited(true);
			ConnectingThroats[i]->SetIsConnectedToOutlet(Fluid, true);
			ConnectingThroats[i]->SweepAdjacentPores(Fluid);
		}
	}
}

void Pore::SetMatrixPlace(unsigned int MPlace){
	MatrixPlace=MPlace;
}

FloatType Pore::GetPressure(void) {
	return Pressure;
}

void Pore::SetPressure(FloatType PoreP) {
	Pressure=PoreP;
}

void Pore::BuildMyRowForAbs(void){
	register unsigned int i;
	FloatType WaterConductancePerLength;
	//FloatType TotalConductancePerLength;
	FloatType AnsEntry, MEntry;
	//FloatType XX;
#ifdef PETSC_ENABLED
	PetscInt irow, icol;
	PetscScalar ival;
#else
	register unsigned int j=1;
#endif

	AnsEntry=0;
	MEntry=0;
	
	for (i=0; i<CoordinationNumber; i++) {
		if ((AdjacentPores[i]==NULL) && (IOStat==(-1))) {		//Inlet
			WaterConductancePerLength=1.0/(Length/(PORELENGTHFACTOR*WaterConductance)+ConnectingThroats[i]->GetLength()/ConnectingThroats[i]->GetWaterConductance());
			
			if (((WaterConductancePerLength-WaterConductancePerLength)!=0) || (!IsConnectedToOutletbyWater)) WaterConductancePerLength=0;
			
			AnsEntry+=WaterConductancePerLength*(DELTAP+ATMOSPHERP);
			MEntry+=WaterConductancePerLength;			
		}
		else if ((AdjacentPores[i]==NULL) && (IOStat==0)) {		//Outlet
			WaterConductancePerLength=1.0/(Length/(PORELENGTHFACTOR*WaterConductance)+ConnectingThroats[i]->GetLength()/ConnectingThroats[i]->GetWaterConductance());
			
			if (((WaterConductancePerLength-WaterConductancePerLength)!=0) || (!IsConnectedToOutletbyWater)) WaterConductancePerLength=0;
			
			AnsEntry+=WaterConductancePerLength*ATMOSPHERP;
			MEntry+=WaterConductancePerLength;
		}
		else {
			WaterConductancePerLength=1.0/(Length/(PORELENGTHFACTOR*WaterConductance)+ConnectingThroats[i]->GetLength()/ConnectingThroats[i]->GetWaterConductance()+AdjacentPores[i]->GetLength()/(PORELENGTHFACTOR*AdjacentPores[i]->GetWaterConductance()));
			
			if (((WaterConductancePerLength-WaterConductancePerLength)!=0) || (!IsConnectedToOutletbyWater)) WaterConductancePerLength=0;
			

			AnsEntry+=WaterConductancePerLength*RhoW*G_ACC*(AdjacentPores[i]->GetZ()-Z);
			
#ifdef PETSC_ENABLED			
			irow = Index;
			icol = AdjacentPores[i]->GetIndex();
			ival = -WaterConductancePerLength*CHANGEFACT;
			ierr = MatSetValues(PetscA, 1, &irow, 1, &icol, &ival, INSERT_VALUES); CHKERRXX(ierr);			
#else
			CoeffMatrix[MatrixPlace + j] = -WaterConductancePerLength*CHANGEFACT;
			Col[MatrixPlace + j] = AdjacentPores[i]->GetIndex();
			j++;
#endif
			MEntry+=WaterConductancePerLength;			
		}
		ConnectingThroats[i]->SetWaterConductancePerLength(WaterConductancePerLength);
		//ConnectingThroats[i]->SetOilConductancePerLength(OilConductancePerLength);
		//ConnectingThroats[i]->SetTotalConductancePerLength(TotalConductancePerLength);
	}
	if (!MEntry) {
		MEntry = 1.0 / CHANGEFACT;
		AnsEntry = ATMOSPHERP / CHANGEFACT;
	}
#ifdef PETSC_ENABLED	
	irow = Index;
	icol = Index;
	ival = MEntry*CHANGEFACT;
	ierr = MatSetValues(PetscA, 1, &irow, 1, &icol, &ival, INSERT_VALUES); CHKERRXX(ierr);
	ierr = VecSetValue(Petscb, irow, AnsEntry, INSERT_VALUES); CHKERRXX(ierr);
	//ival = ATMOSPHERP / CHANGEFACT;
	//ierr = VecSetValue(Petscx, irow, ival, INSERT_VALUES); CHKERRXX(ierr);
#else
	Ans[Index] = AnsEntry;
	CoeffMatrix[MatrixPlace] = MEntry*CHANGEFACT;
	Col[MatrixPlace] = Index;
#endif
}


void Pore::BuildMyRowForWater(void) {
	register unsigned int i;
	FloatType WaterConductancePerLength;
	FloatType AnsEntry, MEntry;
#ifdef PETSC_ENABLED
	PetscInt irow, icol;
	PetscScalar ival;
#else
	register unsigned int j = 1;
#endif

	AnsEntry = 0;
	MEntry = 0;
	
	for (i = 0; i < CoordinationNumber; i++) {
		if ((AdjacentPores[i] == NULL) && (IOStat == (-1))) {		//Inlet
			WaterConductancePerLength = 1.0 / (Length / (PORELENGTHFACTOR * WaterConductance) + ConnectingThroats[i]->GetLength() / ConnectingThroats[i]->GetWaterConductance());
			if (((WaterConductancePerLength - WaterConductancePerLength) != 0) || (!IsConnectedToOutletbyWater)) WaterConductancePerLength = 0;
			//if (((WaterConductancePerLength - WaterConductancePerLength) != 0)) WaterConductancePerLength = 0;

			AnsEntry += WaterConductancePerLength*(DELTAP+ATMOSPHERP);
			MEntry += WaterConductancePerLength;
		}
		else if ((AdjacentPores[i] == NULL) && (IOStat == 0)) {		//Outlet
			WaterConductancePerLength = 1.0 / (Length / (PORELENGTHFACTOR * WaterConductance) + ConnectingThroats[i]->GetLength() / ConnectingThroats[i]->GetWaterConductance());
			if (((WaterConductancePerLength - WaterConductancePerLength) != 0) || (!IsConnectedToOutletbyWater)) WaterConductancePerLength = 0;
			//if (((WaterConductancePerLength - WaterConductancePerLength) != 0) ) WaterConductancePerLength = 0;

			AnsEntry += WaterConductancePerLength*ATMOSPHERP;			
			MEntry += WaterConductancePerLength;

			//ConnectingThroats[i]->SetWaterConductancePerLength(WaterConductancePerLength);
		}
		else {
			WaterConductancePerLength = 1.0 / (Length / (PORELENGTHFACTOR * WaterConductance) + ConnectingThroats[i]->GetLength() / ConnectingThroats[i]->GetWaterConductance() + AdjacentPores[i]->GetLength() / (PORELENGTHFACTOR * AdjacentPores[i]->GetWaterConductance()));
			if (((WaterConductancePerLength - WaterConductancePerLength) != 0) || (!IsConnectedToOutletbyWater)) WaterConductancePerLength = 0;
			//if (((WaterConductancePerLength - WaterConductancePerLength) != 0) ) WaterConductancePerLength = 0;

			AnsEntry += WaterConductancePerLength*RhoW*G_ACC*(AdjacentPores[i]->GetZ() - Z);
#ifdef PETSC_ENABLED			
			irow = Index;
			icol = AdjacentPores[i]->GetIndex();
			ival = -WaterConductancePerLength*CHANGEFACT;
			ierr = MatSetValues(PetscA, 1, &irow, 1, &icol, &ival, INSERT_VALUES); CHKERRXX(ierr);			
#else
			CoeffMatrix[MatrixPlace + j] = -WaterConductancePerLength*CHANGEFACT;
			Col[MatrixPlace + j] = AdjacentPores[i]->GetIndex();
			j++;
#endif
			MEntry += WaterConductancePerLength;					
		}
		ConnectingThroats[i]->SetWaterConductancePerLength(WaterConductancePerLength);
	}
	

	if (!MEntry) {
		MEntry = 1.0 / CHANGEFACT;
		AnsEntry = ATMOSPHERP / CHANGEFACT;
	}
#ifdef PETSC_ENABLED	
	irow = Index;
	icol = Index;
	ival = MEntry*CHANGEFACT;
	ierr = MatSetValues(PetscA, 1, &irow, 1, &icol, &ival, INSERT_VALUES); CHKERRXX(ierr);
	ierr = VecSetValue(Petscb, irow, AnsEntry, INSERT_VALUES); CHKERRXX(ierr);	
#else
	Ans[Index] = AnsEntry;
	CoeffMatrix[MatrixPlace] = MEntry*CHANGEFACT;
	Col[MatrixPlace] = Index;
#endif
}
	

void Pore::BuildMyRowForOil(void) {
	register unsigned int i;
	FloatType OilConductancePerLength;
	FloatType AnsEntry, MEntry;
#ifdef PETSC_ENABLED
	PetscInt irow, icol;
	PetscScalar ival;
#else
	register unsigned int j = 1;
#endif

	
	AnsEntry = 0;
	MEntry = 0;
	for (i = 0; i < CoordinationNumber; i++) {
		if ((AdjacentPores[i] == NULL) && (IOStat == (-1))) {		//Inlet
			OilConductancePerLength = 1.0 / (Length / (PORELENGTHFACTOR * OilConductance) + ConnectingThroats[i]->GetLength() / ConnectingThroats[i]->GetOilConductance());
			if (((OilConductancePerLength - OilConductancePerLength) != 0) || (!IsConnectedToOutletbyOil)) OilConductancePerLength = 0;
			//if (((OilConductancePerLength - OilConductancePerLength) != 0) ) OilConductancePerLength = 0;

			AnsEntry += OilConductancePerLength*(DELTAP + ATMOSPHERP);
			MEntry += OilConductancePerLength;
		}
		else if ((AdjacentPores[i] == NULL) && (IOStat == 0)) {		//Outlet
			OilConductancePerLength = 1.0 / (Length / (PORELENGTHFACTOR * OilConductance) + ConnectingThroats[i]->GetLength() / ConnectingThroats[i]->GetOilConductance());
			if (((OilConductancePerLength - OilConductancePerLength) != 0) || (!IsConnectedToOutletbyOil)) OilConductancePerLength = 0;
			//if (((OilConductancePerLength - OilConductancePerLength) != 0) ) OilConductancePerLength = 0;

			AnsEntry += OilConductancePerLength*ATMOSPHERP;
			MEntry += OilConductancePerLength;

			//ConnectingThroats[i]->SetOilConductancePerLength(OilConductancePerLength);
		}
		else {
			OilConductancePerLength = 1.0 / (Length / (PORELENGTHFACTOR * OilConductance) + ConnectingThroats[i]->GetLength() / ConnectingThroats[i]->GetOilConductance() + AdjacentPores[i]->GetLength() / (PORELENGTHFACTOR * AdjacentPores[i]->GetOilConductance()));
			if (((OilConductancePerLength - OilConductancePerLength) != 0) || (!IsConnectedToOutletbyOil)) OilConductancePerLength = 0;
			//if (((OilConductancePerLength - OilConductancePerLength) != 0) ) OilConductancePerLength = 0;
						
			AnsEntry += OilConductancePerLength*RhoO*G_ACC*(AdjacentPores[i]->GetZ() - Z);
#ifdef PETSC_ENABLED			
			irow = Index;
			icol = AdjacentPores[i]->GetIndex();
			ival = -(OilConductancePerLength)*CHANGEFACT;
			ierr = MatSetValues(PetscA, 1, &irow, 1, &icol, &ival, INSERT_VALUES); CHKERRXX(ierr);			
#else
			CoeffMatrix[MatrixPlace + j] = -(OilConductancePerLength)*CHANGEFACT;
			Col[MatrixPlace + j] = AdjacentPores[i]->GetIndex();
			j++;
#endif		
			MEntry += OilConductancePerLength;	
		}
		//if (OilConductancePerLength) Zeros = false;
		ConnectingThroats[i]->SetOilConductancePerLength(OilConductancePerLength);
	}
	

	if (!MEntry) {
		MEntry = 1.0 / CHANGEFACT;
		AnsEntry = ATMOSPHERP / CHANGEFACT;
	}

#ifdef PETSC_ENABLED	
	irow = Index;
	icol = Index;
	ival = MEntry*CHANGEFACT;
	ierr = MatSetValues(PetscA, 1, &irow, 1, &icol, &ival, INSERT_VALUES); CHKERRXX(ierr);
	ierr = VecSetValue(Petscb, irow, AnsEntry, INSERT_VALUES); CHKERRXX(ierr);	
#else
	Ans[Index] = AnsEntry;
	CoeffMatrix[MatrixPlace] = MEntry*CHANGEFACT;
	Col[MatrixPlace] = Index;
#endif
}

void Pore::SweepAdjacentThroatsForDeletion(ExitType mExit) {
	register unsigned int i;

	if (mExit == OUTLET) {
		for (i = 0; i < CoordinationNumber; i++) {
			if ((ConnectingThroats[i] != NULL) && (!ConnectingThroats[i]->GetIsVisited())) {
				ConnectingThroats[i]->SetIsVisited(true);
				ConnectingThroats[i]->SetIsConnectedToOutlet(WATER, true);
				ConnectingThroats[i]->SweepAdjacentPoresForDeletion(OUTLET);
			}
		}
	}
	else {
		for (i = 0; i < CoordinationNumber; i++) {
			if ((ConnectingThroats[i] != NULL) && (!ConnectingThroats[i]->GetIsVisited())) {
				ConnectingThroats[i]->SetIsVisited(true);
				ConnectingThroats[i]->SetIsConnectedToInlet(true);
				ConnectingThroats[i]->SweepAdjacentPoresForDeletion(INLET);
			}
		}
	}
	
}


unsigned int Pore::GetCoordinationNumber(void) {
	return CoordinationNumber;
}

bool Pore::IsConnectedThroatNull(unsigned int AdjacentPoresIndex) {
	return (AdjacentPores[AdjacentPoresIndex] == NULL);
}

unsigned int Pore::NumberOfConnections(void) {
	register unsigned int i, j;

	j = 0;
	for (i = 0; i < CoordinationNumber; i++) if (AdjacentPores[i] != NULL) j++;
	return j;
}