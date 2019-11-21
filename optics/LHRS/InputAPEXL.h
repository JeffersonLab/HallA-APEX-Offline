// Description of APEX optics with PREX target



#ifndef ROOT_Input
#define ROOT_Input


#include "TROOT.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"


const TString InputID = "APEX_PREX";

/////////////////////////////////////////////////////////////////////////
// Input Sections
/////////////////////////////////////////////////////////////////////////
// HRS Positioning
// The central ray of the spectrometer is at -15.993 degrees.
// It is missing the defined target center by 1.63 mm upstream,
// and 0.98 mm vertically (positive = up).

const Double_t D2R = TMath::Pi() / 180.; // degrees to radians conversion

// const Double_t HRSAngle = 12.511/180.*TMath::Pi(); // this is HRS angle
const Double_t HRSAngle = 5.366 * D2R; // this is Sieve angle (in front of septum)

/* const Double_t MissPointZ = -2.22e-3; // APEX = 1.690 */
/* const Double_t MissPointY = -2.87e-3; // APEX = -1.790 */

const Double_t MissPointZ = 1.690e-3;
const Double_t MissPointY = -1.790e-3;



TVector3 fPointingOffset;

//fPointingOffset.SetXYZ(-MissPointZ*TMath::Sin(HRSAngle)*TMath::Cos(HRSAngle),(Double_t)MissPointY,MissPointZ*TMath::Sin(HRSAngle)*TMath::Sin(HRSAngle));





//const Double_t BeamShiftX = 8.8e-3*2;
const Double_t BeamShiftX = 0;


// Used to illustrate where sieve holes should be in theta- phi plot

/* // rastered beam */
const Double_t BeamX_average_V1 = 0.001234;
const Double_t BeamY_average_V1 = 0.002553;

const Double_t BeamX_average_V2 = 0.001618;
const Double_t BeamY_average_V2 = 0.002543;

const Double_t BeamX_average_V3 = 0.001666;
const Double_t BeamY_average_V3 = 0.002538;




//unrastered beam
/* const Double_t BeamX_average = 0.001239; */
/* const Double_t BeamY_average = 0.002545; */

//test (0,0) beam position
/* const Double_t BeamX_average = 0.0; */
/* const Double_t BeamY_average = 0.0; */



//#define DISABLED_SIEVE_HOLE	0
#define DISABLED_SIEVE_HOLE	(Col==7 && Row==9 && FoilID>=7)
// JW: think this just means that these particular holes are skipped in the optimisation process

/////////////////////////////////////////////////////////////////////////
// for Calculating real sieve positions;
/////////////////////////////////////////////////////////////////////////


// -------------------------------------------
//	From PREX people
// -------------------------------------------

//const Double_t SieveYbyCol[] ={
//		0.8435 	* 25.400e-3,
//		0.723 	* 25.400e-3,
//		0.6025 	* 25.400e-3,
//		0.482 	* 25.400e-3,
//		0.241 	* 25.400e-3,
//		0.1205 	* 25.400e-3,
//		0 		* 25.400e-3,
//		-0.188 	* 25.400e-3,
//		-0.282 	* 25.400e-3,
//		-0.376 	* 25.400e-3,
//		-0.564 	* 25.400e-3,
//		-0.658 	* 25.400e-3,
//		1e36
//};
//
//const UInt_t NSieveCol = 12;/*WARNING:12*/
//
//
//const Double_t SieveXbyRow[] ={
//		1.834 * 25.400e-3,
//		1.572 * 25.400e-3,
//		1.310 * 25.400e-3,
//		1.048 * 25.400e-3,
//		0.786 * 25.400e-3,
//		0.524 * 25.400e-3,
//		0.262 * 25.400e-3,
//		0 * 25.400e-3,
//		-0.262 * 25.400e-3,
//		-0.524 * 25.400e-3,
//		-0.786 * 25.400e-3,
//		-1.048 * 25.400e-3,
//		-1.310 * 25.400e-3,
//		-1.572 * 25.400e-3,
//		-1.834 * 25.400e-3,
//		1e36
//};

// -------------------------------------------
//	From Eric Jensen
// -------------------------------------------

/* const Double_t SieveYbyCol[] ={ */
/* 		0.837 	* 25.400e-3,	//0 */
/* 		0.716 	* 25.400e-3,	//1 */
/* 		0.599 	* 25.400e-3,	//2 */
/* 		0.476 	* 25.400e-3,	//3 */
/* 		0.237 	* 25.400e-3,	//4 */
/* 		0.095 	* 25.400e-3,	//5 */
/* 		0 	* 25.400e-3,	//6 */
/* 		-0.189 	* 25.400e-3,	//7 */
/* 		-0.283 	* 25.400e-3,	//8 */
/* 		-0.368 	* 25.400e-3,	//9 */
/* 		-0.555 	* 25.400e-3,	//10 */
/* 		-0.650 	* 25.400e-3,	//11 */
/* 		1e36,1e36,1e36,1e36,1e36 */
/* }; */

/* const UInt_t NSieveCol = 12;/\*WARNING:12*\/ */


/* const Double_t SieveXbyRow[] ={ */
/* 		1.803 * 25.400e-3,	//0 */
/* 		1.544 * 25.400e-3,	//1 */
/* 		1.291 * 25.400e-3,	//2 */
/* 		1.031 * 25.400e-3,	//3 */
/* 		0.772 * 25.400e-3,	//4 */
/* 		0.516 * 25.400e-3,	//5 */
/* 		0.257 * 25.400e-3,	//6 */
/* 		0 * 25.400e-3,	//7 */
/* 		-0.265 * 25.400e-3,	//8 */
/* 		-0.522 * 25.400e-3,	//9 */
/* 		-0.781 * 25.400e-3,	//10 */
/* 		-1.046 * 25.400e-3,	//11 */
/* 		-1.303 * 25.400e-3,	//12 */
/* 		-1.556 * 25.400e-3,	//13 */
/* 		-1.821 * 25.400e-3,	//14 */
/* 		1e36,1e36,1e36,1e36,1e36 */
/* }; */



// -------------------------------------------
//	From John Williamson
// -------------------------------------------
//--------------------------------------------



const Double_t YbyCol = 0.094 * 25.4e-3;

const Double_t SieveYbyCol[] ={ -14*YbyCol, -13*YbyCol, -12*YbyCol, -11*YbyCol, -10*YbyCol, -9*YbyCol, -8*YbyCol, -7*YbyCol, -6*YbyCol,  -5*YbyCol, -4*YbyCol, -3*YbyCol, -2*YbyCol, -1*YbyCol, 0, 1*YbyCol, 2*YbyCol, 3*YbyCol, 4*YbyCol, 5*YbyCol, 6*YbyCol, 7*YbyCol, 8*YbyCol, 9*YbyCol, 10*YbyCol, 12*YbyCol, 14*YbyCol,1e36,1e36};

const UInt_t NSieveCol = 27; 
const UInt_t NSieveCol_even = 15; 


const Int_t NHoles = 225; 
//const UInt_t NHoles = 225;

/* const Double_t SieveYbyColEven[]= {14*YbyCol, 12*YbyCol, 10*YbyCol, 8*YbyCol, 6*YbyCol, 4*YbyCol, 2*YbyCol, 0, -2*YbyCol, -4*YbyCol, -6*YbyCol, -8*YbyCol, -10*YbyCol, -12*YbyCol, -14*YbyCol}; */
/* const Double_t SieveYbyColOdd[]= {13*YbyCol, 11*YbyCol,9*YbyCol, 7*YbyCol, 5*YbyCol, 3*YbyCol, 1*YbyCol, -1*YbyCol, 0, -3*YbyCol, -5*YbyCol, -7*YbyCol, -9*YbyCol}; */




const Double_t XbyRow = 0.23 * 25.4e-3;
const Double_t SieveXbyRow[] = { -8 * XbyRow, -7 * XbyRow, -6 * XbyRow, -5 * XbyRow, -4 * XbyRow, -3 * XbyRow, -2 * XbyRow, -1 * XbyRow, 0.0, 1 * XbyRow, 2 * XbyRow, 3 * XbyRow, 4 * XbyRow, 5 * XbyRow, 6 * XbyRow, 7  * XbyRow, 8 * XbyRow, 1e36};


// List number of holes in each row

const UInt_t NoinEachRow[] = {15, 12, 15, 11, 15, 11, 15, 11, 15, 11, 15, 11, 15, 11, 15, 12, 15};



const UInt_t NSieveRow = 17;/*WARNING:15*/
/////////////////////////////////////////////////////////////////////////


// 							Z X Y Angle
// 		B.L. Sieve Center 798.02 69.91 -1.50 5.007

// For Previous 
/* const Double_t SieveOffX = 1.4*1e-3; */
/* const Double_t SieveOffY = (81.0/TMath::Cos(HRSAngle) - (TMath::Tan(HRSAngle)*81.0+791.8)*TMath::Sin(HRSAngle))*1e-3; */
/* const Double_t ZPos = 1157*1e-3; */


//JW: for APEX
/*
                      Z       X        Y     Angle
B.L. Sieve Centre     791.8   81.0    1.4    5.366

Calculations:

const Double_t SieveOffX = Y;

const Double_t SieveOffY = ( (X/TMath::Cos(HRSAngle)) - ( (X*TMath::Tan(HRSAngle)) + Z )*TMath::Sin(HRSAngle)); 

Double_t ZPos = (Z + (X * TMath::Tan(HRSAngle))) * TMath::Cos(HRSAngle)*1e-3;

Work out as:

SieveOffX = 1.4
SieveOffY = 6.59786
ZPos      = 795.905
(in mm)

 */


const Double_t SieveOffX = 1.4e-3;
const Double_t SieveOffY = 6.59786e-3;
const Double_t ZPos = 795.905e-3;





//  Double_t ZPos = (Z + (X * TMath::Tan(HRSAngle))) * TMath::Cos(HRSAngle)*1e-3;




// 



/////////////////////////////////////////////////////////////////////////
// Vertex Position Inputs

//  N Carbon foil settings could be different !!!

// static const UInt_t NFiols = 5;//WARNING: check
// const Double_t targetfoils[]={-0.15,-0.075,0,0.075,0.15};

//for Sieve and dp calibration
// static const UInt_t NFiols = 2*5;//WARNING: check
// const Double_t targetfoils[]={0,0.075,0,0.075,0,0.075,0,0.075,0,0.075,0,0.075,0,0.075,0,0.075};

////for vertex
//static const UInt_t NFiols = 2;//WARNING: check
////const Double_t targetfoils[]={-0.15,-0.075/*,0*//*,0.075*/,0.15,-0.0044};
//const Double_t targetfoils[]={-0.0044, 0.15};

//for Central Foil
/* static const UInt_t NFiols = 7+5;//WARNING: check */
/* const Double_t targetfoils[] = { -0.0044, -0.0044, -0.0044, -0.0044, -0.0044, */
/* 		-0.0044, -0.0044, -0.0044, -0.0044, -0.0044, -0.0044, -0.0044, -0.0044, */
/* 		-0.0044, -0.0044 }; */


static const UInt_t NFiols = 3;//WARNING: check
//const Double_t targetfoils[] = { -0.20, 0.0, 0.2 };
//const Double_t targetfoils[] = { -0.195, 0.005, 0.205 };
const Double_t targetfoils[] = { -0.205, -0.005, 0.195 };


// Optics 1 foils

/* static const UInt_t NFiols = 4;//WARNING: check */
/* const Double_t targetfoils[] = { -0.295, -0.145, 0.08, 0.224 }; */


// Optics 3 foils

/* static const UInt_t NFiols = 4;//WARNING: check */
/* const Double_t targetfoils[] = { -0.214, -0.07, 0.155, 0.305 }; */


// Optics 2 foils

/* static const UInt_t NFiols = 8;//WARNING: check */
/* const Double_t targetfoils[] = { -0.295,  -0.214, -0.145,  -0.07, 0.08, 0.155, 0.224, 0.305 }; */

// Optics 2 foils

/* static const UInt_t NFiols = 4;//WARNING: check */
/* const Double_t targetfoils[] = { -0.214, -0.07, 0.155, 0.305 }; */



/////////////////////////////////////////////////////////////////////////
// Exciation state Inputs

const UInt_t NKine = 5; //N Delta Scans

#define DIPOLE_MAG2MOM(Mag) (2.702*Mag-1.6e-03*Mag*Mag*Mag)
const Double_t HRSCentralMom[]=
{
		DIPOLE_MAG2MOM((0.42943+0.429431)/2), // %0 run 1170
		DIPOLE_MAG2MOM((0.425724+0.425724)/2), // %1 run 1172
		DIPOLE_MAG2MOM((0.421387+0.421387)/2), // %2 run 1181
		DIPOLE_MAG2MOM((0.417064+0.417064)/2), // %3 run 1183
		DIPOLE_MAG2MOM((0.412774+0.412774)/2), // %4 run 1190 + 1191
		0
};

const Double_t GroundNuclearMass = 180.94788*.931494028-.511e-3*73; //GeV/c^2  // Ta Target
const Double_t ExcitationEnergy[]= //selected exciation states for each kinematics
{0.,0,0,0,0};

const UInt_t NExcitationStates = 1; //Ta Excitation States
const Double_t ExcitationEnergyList[]
                                    ={0};

/////////////////////////////////////////////////////////////////////////
// Radiation loss Inputs

// Xiaodong's version
// const Double_t AllLossExceptFoil = 1e-3*(
// 										 .075//Be Window
// 										+0.374//BeO
// 										+0.0168+0.0228//He4 in target enclosure
// 										+0.109//air, Target Enclosure to HRS
// 										+.1//target enclosure
// 										+.1//HRS Entrance
// 										);
// const Double_t LossEachFoil = 1e-3*(0.08775);

//revised Version
const Double_t AllLossExceptFoil
= 1e-3*(//in MeV
		0
);
const Double_t LossEachFoil = 1e-3*(0.1*1.598);

//Array of FoilID
const Double_t RadiationLossByFoil[]={
		AllLossExceptFoil+LossEachFoil*1,
		AllLossExceptFoil+LossEachFoil*2,
		AllLossExceptFoil+LossEachFoil*3,
		AllLossExceptFoil+LossEachFoil*4,
		AllLossExceptFoil+LossEachFoil*5,
		AllLossExceptFoil+LossEachFoil*6,
		AllLossExceptFoil+LossEachFoil*7
};
//Warning: these numbers are calculated with small angle approximation

/////////////////////////////////////////////////////////////////////////
// Disable Extended Target Correction

//const Double_t ExtTarCor_ThetaCorr = 0;
//const Double_t ExtTarCor_DeltaCorr = 1e36;

const Double_t ExtTarCor_ThetaCorr = 0.0;
//const Double_t ExtTarCor_ThetaCorr = 0.61;
const Double_t ExtTarCor_DeltaCorr = 0;


/////////////////////////////////////////////////////////////////////////

const char * Prefix = "L.vdc.";

const char * DatabaseHeader = "\
[ L.global ]\n\
0.3327 1 0.0 270.2 0.0 -1.6e-03        VDC Angle, Plane Spacing, Gamma Coefficents\n\
L.vdc.matrixelem = \n\
t 0 0 0  -1.001135e+00 -3.313373e-01 -4.290819e-02  4.470852e-03  0.000000e+00  0.000000e+00  0.000000e+00  0\n\
y 0 0 0  -8.060915e-03  1.071977e-03  9.019102e-04 -3.239615e-04  0.000000e+00  0.000000e+00  0.000000e+00  0\n\
p 0 0 0  -2.861912e-03 -2.469069e-03  8.427172e-03  2.274635e-03  0.000000e+00  0.000000e+00  0.000000e+00  0\n\
";

const char * DatabaseFooter = "\
L 0 0 0 0  24.216\n\
L 0 1 0 0 -13.1041\n\
L 0 2 0 0  20.8672\n\
L 0 1 1 0  -1.8728\n\
L 1 0 0 0   0.0482\n\
L 2 0 0 0   0.0675\n\
L 1 0 0 1   0.0281\n\
L 1 1 0 0   2.4314\n\
L 0 0 0 1   0.1408\n\
L 1 0 1 0  -0.0845\n\
L 0 0 0 2   9.3063\n\
L 0 0 2 0   6.3451\n\
L 0 0 1 1  -7.7394\n\
L 0 0 1 0  -0.1443\n\
L 0 1 0 1   1.7292\n\
";

/////////////////////////////////////////////////////////////////////////


#endif


