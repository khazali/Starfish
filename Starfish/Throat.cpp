#include <math.h>
#include "Throat.h"
#include "Pore.h"
#include "ElementList.h"
#include "MIfstream.h"

Throat::Throat(void) {

}
void Throat::ReadLink1(MIfstream& InputFile) {
	char str[MAX_STRING_LENGTH];
	register int i, j;
	FloatType Z0, Z1;

	if (!InputFile.ReadWord(str)) TerM("Incorrect link1 file format!");
	if (!InputFile.ReadWord(str)) TerM("Incorrect link1 file format!");
	i=atoi(str);
	switch (i) {
		case -1:
			//IOStat=-1;
			ConnectingPores[0]=NULL;
			Z0=0;
			InletThroats.AddElement(Index);
			break;
		case 0:
			//IOStat=0;
			ConnectingPores[0]=NULL;
			Z0=0;
			OutletThroats.AddElement(Index);
			break;
		default:
			//IOStat=1;
			ConnectingPores[0]=&pores[i-1];
			Z0=pores[i-1].GetZ();			
	}	
	if (!InputFile.ReadWord(str)) TerM("Incorrect link1 file format!");
	j=atoi(str);
	switch (j) {
		case -1:
			//IOStat=-1;
			ConnectingPores[1]=NULL;
			Z1=0;
			InletThroats.AddElement(Index);
			break;
		case 0:
			//IOStat=0;
			ConnectingPores[1]=NULL;
			Z1=0;
			OutletThroats.AddElement(Index);
			break;
		default:
			//IOStat=1;
			ConnectingPores[1]=&pores[j-1];
			Z1=pores[j-1].GetZ();
	}
	if (Z0 && Z1) Z=(Z0+Z1)/2;
	else Z=Z0+Z1;

	if ((i==(-1)) || (j==(-1))) IOStat=-1;
	else if ((i==0) || (j==0)) IOStat=0;
	else IOStat=1;

	if (!InputFile.ReadWord(str)) TerM("Incorrect link1 file format!");
	InscribedRadius=atof(str);

	if (!InputFile.ReadWord(str)) TerM("Incorrect link1 file format!");
	ShapeFactor=atof(str);

	if (!InputFile.ReadWord(str)) TerM("Incorrect link1 file format!");
	TotalLength=atof(str);
}

void Throat::ReadLink2(MIfstream& InputFile) {
	char str[MAX_STRING_LENGTH];

	if (!InputFile.ReadWord(str)) TerM("Incorrect link2 file format!");
	if (!InputFile.ReadWord(str)) TerM("Incorrect link2 file format!");
	if (!InputFile.ReadWord(str)) TerM("Incorrect link2 file format!");

	if (!InputFile.ReadWord(str)) TerM("Incorrect link2 file format!");
	if (ConnectingPores[0]!=NULL) ConnectingPores[0]->SetLength(atof(str));

	if (!InputFile.ReadWord(str)) TerM("Incorrect link2 file format!");
	if (ConnectingPores[1]!=NULL) ConnectingPores[1]->SetLength(atof(str));

	if (!InputFile.ReadWord(str)) TerM("Incorrect link2 file format!");
	Length=atof(str);

	if (!InputFile.ReadWord(str)) TerM("Incorrect link2 file format!");
	Volume=atof(str);
	TotalVolume+=Volume;

	if (!InputFile.ReadWord(str)) TerM("Incorrect link2 file format!");
	ClayVolume=atof(str);

	TotalVolume+=ClayVolume;
}

unsigned int Throat::CalculateImbibitionPc(FloatType WorkingPc) {
	FloatType S1, S2, S3, S4, S5, S6, S7, S8;
	register unsigned int i, n;
	FloatType costr;
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
		if (WorkingPc>=0) {						//Secondary Water Flooding- Spontaneous Imbibition
			IsReadyToFill=(((IOStat==-1) ||
				((ConnectingPores[0]!=NULL) && (ConnectingPores[0]->GetIsWaterFilled())) ||
				((ConnectingPores[1]!=NULL) && (ConnectingPores[1]->GetIsWaterFilled()))));
			
			
			//Piston Like
			if (!NumberOfSides) {				
				tPc1=2*IFT*cos(AdvancingContactAngle)/InscribedRadius;					//?????????????????
			}
			else if (NumberOfSides == 3) {
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
				tPc1=IFT/R;
				//tSat=Aeff/A;
				SumForR = 0;
				for (i = 0; i < 3; i++) {
					if (Beta[i] < (PI / 2 - RecedingContactAngle)) {
						SumForR += cos(HingeAngle[i])*cos(HingeAngle[i] + Beta[i]) / sin(Beta[i]) + HingeAngle[i] + Beta[i] - PI / 2;
					}
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
				tPc1 = IFT / R;
				//tSat=Aeff/A;
				SumForR = 0;
				if (Beta[0] < (PI / 2 - RecedingContactAngle)) {
					SumForR = NumberOfSides*(cos(HingeAngle[0])*cos(HingeAngle[0] + Beta[0]) / sin(Beta[0]) + HingeAngle[0] + Beta[0] - PI / 2);
				}				
			}
			

			//Snap Off
			if (NumberOfSides==3) {
				n=0;
				for (i=0; i<3; i++) {
					if ((AdvancingContactAngle<(PI/2-Beta[i])) && (HingeAngle[i]>=AdvancingContactAngle)) n++;
				}

				if (n==1) {
					tPc2=(IFT/InscribedRadius)*((cos(AdvancingContactAngle)/tan(Beta[0])-sin(AdvancingContactAngle)+cos(HingeAngle[2])/tan(Beta[2])-sin(HingeAngle[2]))/(1.0/tan(Beta[0])+1/tan(Beta[1])));
				}
				else if (n>1) {
					tPc2=(IFT/InscribedRadius)*(cos(AdvancingContactAngle)-2*sin(AdvancingContactAngle)/(1.0/tan(Beta[0])+1/tan(Beta[1])));
				}
				else {
					tPc2 = -CAPILLARYLIMIT - 1;
				}
			}
			else if (NumberOfSides > 3) {
				n = 0;				
				if ((AdvancingContactAngle < (PI / 2 - Beta[0])) && (HingeAngle[0] >= AdvancingContactAngle)) n=NumberOfSides;
				

				if (n == 0) {
					tPc2 = -CAPILLARYLIMIT - 1;
					//tPc2 = (IFT / InscribedRadius)*((cos(AdvancingContactAngle) / tan(Beta[0]) - sin(AdvancingContactAngle) + cos(HingeAngle[2]) / tan(Beta[2]) - sin(HingeAngle[2])) / (1.0 / tan(Beta[0]) + 1 / tan(Beta[1])));
				}
				else if (n > 0) {
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
		else {				//Secondary Water Flooding- Forced Imbibition

			IsReadyToFill=(((IOStat==-1) ||
				((ConnectingPores[0]!=NULL) && (ConnectingPores[0]->GetIsWaterFilled())) ||
				((ConnectingPores[1]!=NULL) && (ConnectingPores[1]->GetIsWaterFilled()))));


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
			else if (NumberOfSides > 3) {
				SumForR = 0;
				if (Beta[0] < (PI / 2 - AdvancingContactAngle)) {
					SumForR = NumberOfSides*(costr*cos(AdvancingContactAngle + Beta[0]) / sin(Beta[0]) + AdvancingContactAngle + Beta[0] - PI / 2);
				}				
				Fd = (1 + sqrt(1 - 4 * ShapeFactor*SumForR / (costr*costr))) / (1 + 2 * sqrt(PI*ShapeFactor));
			}
			tPc1 = -IFT*costr*(1 + 2 * sqrt(PI*ShapeFactor))*Fd / InscribedRadius;
			


			
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

				if (AdvancingContactAngle <= (PI - Beta[0])) {
					tPc2 = IFT*cos(AdvancingContactAngle + Beta[0]) / (Rmin*cos(RecedingContactAngle + Beta[0]));
				}
				else {
					tPc2 = -IFT / (Rmin*cos(RecedingContactAngle + Beta[0]));
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
			

unsigned int Throat::TestFillReadiness(FloatType WorkingPc) {
	if ((!IsOilFilled) && (WorkingPc>=EntryCapillaryPressure) && (IsConnectedToOutletbyWater) && ((IOStat==-1) ||
		((ConnectingPores[0]!=NULL) && (ConnectingPores[0]->GetIsOilFilled())) ||
		((ConnectingPores[1]!=NULL) && (ConnectingPores[1]->GetIsOilFilled())))) {			
			FillWithOil();
			return 1;
	}
	else return 0;
}

void Throat::SweepAdjacentPores(FluidType Fluid) {
	if ((ConnectingPores[0]!=NULL) && (!ConnectingPores[0]->GetIsVisited()) && (((Fluid==OIL) && (ConnectingPores[0]->GetIsOilFilled())) || ((Fluid==WATER) && (ConnectingPores[0]->GetWaterCoatingExist())))) {
		ConnectingPores[0]->SetIsVisited(true);
		ConnectingPores[0]->SetIsConnectedToOutlet(Fluid, true);
		ConnectingPores[0]->SweepAdjacentThroats(Fluid);
	}
	if ((ConnectingPores[1]!=NULL) && (!ConnectingPores[1]->GetIsVisited()) && (((Fluid==OIL) && (ConnectingPores[1]->GetIsOilFilled())) || ((Fluid==WATER) && (ConnectingPores[1]->GetWaterCoatingExist())))) {
		ConnectingPores[1]->SetIsVisited(true);
		ConnectingPores[1]->SetIsConnectedToOutlet(Fluid, true);
		ConnectingPores[1]->SweepAdjacentThroats(Fluid);
	}	
}

FloatType Throat::GetTotalLength(void){
	return TotalLength;
}

FloatType Throat::GetOutletPorePressure(void){
	if (ConnectingPores[0]!=NULL) {
		return ConnectingPores[0]->GetPressure();
	}
	else {
		return ConnectingPores[1]->GetPressure();
	}
}

void Throat::SetWaterConductancePerLength(FloatType WCondPerL) {
	WaterConductancePerLength=WCondPerL;
}

void Throat::SetOilConductancePerLength(FloatType OCondPerL) {
	OilConductancePerLength=OCondPerL;
}

void Throat::SetTotalConductancePerLength(FloatType TCondPerL) {
	TotalConductancePerLength=TCondPerL;
}

FloatType Throat::GetWaterConductancePerLength(void) {
	return WaterConductancePerLength;
}

FloatType Throat::GetOilConductancePerLength(void) {
	return OilConductancePerLength;
}

FloatType Throat::GetTotalConductancePerLength(void) {
	return TotalConductancePerLength;
}

void Throat::SweepAdjacentPoresForDeletion(ExitType mExit) {
	if (mExit == OUTLET) {
		if ((ConnectingPores[0] != NULL) && (!ConnectingPores[0]->GetIsVisited())) {
			ConnectingPores[0]->SetIsVisited(true);
			ConnectingPores[0]->SetIsConnectedToOutlet(WATER, true);
			ConnectingPores[0]->SweepAdjacentThroatsForDeletion(OUTLET);
		}
		if ((ConnectingPores[1] != NULL) && (!ConnectingPores[1]->GetIsVisited())) {
			ConnectingPores[1]->SetIsVisited(true);
			ConnectingPores[1]->SetIsConnectedToOutlet(WATER, true);
			ConnectingPores[1]->SweepAdjacentThroatsForDeletion(OUTLET);
		}
	}
	else {
		if ((ConnectingPores[0] != NULL) && (!ConnectingPores[0]->GetIsVisited())) {
			ConnectingPores[0]->SetIsVisited(true);
			ConnectingPores[0]->SetIsConnectedToInlet(true);
			ConnectingPores[0]->SweepAdjacentThroatsForDeletion(INLET);
		}
		if ((ConnectingPores[1] != NULL) && (!ConnectingPores[1]->GetIsVisited())) {
			ConnectingPores[1]->SetIsVisited(true);
			ConnectingPores[1]->SetIsConnectedToInlet(true);
			ConnectingPores[1]->SweepAdjacentThroatsForDeletion(INLET);
		}
	}
	
}

void Throat::operator=(const Throat &AnotherThroat) {
	register unsigned int j;

	if (Index != AnotherThroat.Index) {
		if (AnotherThroat.ConnectingPores[0] == NULL) {
			ConnectingPores[0] = NULL;
		}
		else {
			j = AnotherThroat.ConnectingPores[0]->GetIndex();
			ConnectingPores[0] = &pores[j - PorePointer[j]];
		}

		if (AnotherThroat.ConnectingPores[1] == NULL) {
			ConnectingPores[1] = NULL;
		}
		else {
			j = AnotherThroat.ConnectingPores[1]->GetIndex();
			ConnectingPores[1] = &pores[j - PorePointer[j]];
		}
		//What to do about exit throats added to list?

		IOStat = AnotherThroat.IOStat;
		Z = AnotherThroat.Z;
		InscribedRadius = AnotherThroat.InscribedRadius;
		ShapeFactor = AnotherThroat.ShapeFactor;
		TotalLength = AnotherThroat.TotalLength;
		Length = AnotherThroat.Length;
		TotalVolume -= Volume + ClayVolume;
		Volume = AnotherThroat.Volume;
		ClayVolume = AnotherThroat.ClayVolume;
		//what to do about total volume?
		
		/*AbsoluteConductance = AnotherThroat.AbsoluteConductance;
		Beta[0] = AnotherThroat.Beta[0];
		Beta[1] = AnotherThroat.Beta[1];
		Beta[2] = AnotherThroat.Beta[2];
		IsCircle = AnotherThroat.IsCircle;

		Fd = AnotherThroat.Fd;
		SumForR = AnotherThroat.SumForR;
		EntryCapillaryPressure = AnotherThroat.EntryCapillaryPressure;*/
	}	
}

int Throat::GetIOStat(void) {
	return IOStat;
}

void Throat::CalcWaterVelocity(void) {
	FloatType P1, P2;

	if ((IOStat == (-1)) && (ConnectingPores[0] == NULL)) P1 = DELTAP + ATMOSPHERP;
	else if ((IOStat == 0) && (ConnectingPores[0] == NULL)) P1 = ATMOSPHERP;
	else P1 = ConnectingPores[0]->GetPressure();
	if ((IOStat == (-1)) && (ConnectingPores[1] == NULL)) P2 = DELTAP + ATMOSPHERP;
	else if ((IOStat == 0) && (ConnectingPores[1] == NULL)) P2 = ATMOSPHERP;
	else P2 = ConnectingPores[1]->GetPressure();

	WaterVelocity = WaterConductancePerLength*fabs(P1-P2) / SurfaceArea;	
}

void Throat::CalcOilVelocity(void) {
	FloatType P1, P2;

	if ((IOStat == (-1)) && (ConnectingPores[0] == NULL)) P1 = DELTAP + ATMOSPHERP;
	else if ((IOStat == 0) && (ConnectingPores[0] == NULL)) P1 = ATMOSPHERP;
	else P1 = ConnectingPores[0]->GetPressure();
	if ((IOStat == (-1)) && (ConnectingPores[1] == NULL)) P2 = DELTAP + ATMOSPHERP;
	else if ((IOStat == 0) && (ConnectingPores[1] == NULL)) P2 = ATMOSPHERP;
	else P2 = ConnectingPores[1]->GetPressure();

	OilVelocity = OilConductancePerLength*fabs(P1 - P2) / SurfaceArea;
}