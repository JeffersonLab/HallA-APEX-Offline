// Description of GMP optics

#ifndef ROOT_Input
#define ROOT_Input

#include "TROOT.h"
#include "TMath.h"

using namespace std;
const TString InputID = "GMp_RHRS";

/////////////////////////////////////////////////////////////////////////
// Input Sections
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
//const Double_t Ebeam = 2.3004;
// HRS Position Inputs
const Double_t D2R = TMath::Pi() / 180.;

//set as central sieve hole angle
const Double_t HRSAngle = -5 * D2R; 

//LH2 target information
const Double_t LH2_TargetLength = 15*1.e-2; //unit m
const Double_t LH2_Target_Tip_Radius = 1.5*2.54*1.e-2; //target width and tip radius
const Double_t LH2_Thickness_Entance = 0.175*1.e-3;//Al 7075, aluminum thickness for the entrance window
const Double_t LH2_Thickness_Side = 0.18*1.e-3;  //Al 7075, aluminum thickness for the side wall
const Double_t LH2_Thickness_Tip = 0.11*1.e-3;  //Al 7075, aluminum thickness for tip

// MissPoint* are in HCS
const Double_t MissPointZ = 0;//1.69*1e-3 //RHRS offsets 
const Double_t MissPointY = 0;//-1.79*1e-3 

const Double_t BeamShiftX = 0;

const Double_t SieveRadius = 0.055*25.4/2.0*1e-3;
const Double_t SieveRadius_c = 0.106*25.4/2.0*1e-3;

//const Double_t SieveRadius = 0;       
//const Double_t SieveRadius_c = 0;   

// average beam positions
const Double_t BeamX_average[] = {1.76e-3, -0.9646e-3, -2.137e-3};
const Double_t BeamY_average[] = {2.726e-3, 2.845*1e-3, 2.879e-3};


/////////////////////////////////////////////////////////////////////////
// Sieve Position Inputs
const Double_t YbyCol = .19 * 25.4e-3;
const Double_t y_off = -31.23 * tan(0.8 *D2R) * 25.4e-3; 
const Double_t SieveYbyCol[]= {7*YbyCol + y_off, 6.5*YbyCol + y_off, 6*YbyCol + y_off, 5.5*YbyCol + y_off, 5*YbyCol + y_off, 4.5*YbyCol + y_off, 4*YbyCol + y_off, 3.5*YbyCol + y_off, 3*YbyCol + y_off, 2.5*YbyCol + y_off, 2*YbyCol + y_off, 1.5*YbyCol + y_off, 1.0*YbyCol + y_off, 0.5*YbyCol + y_off, 0.0*YbyCol + y_off, -0.5*YbyCol + y_off, -1*YbyCol + y_off, -1.5*YbyCol + y_off, -2*YbyCol + y_off, -2.5*YbyCol + y_off, -3*YbyCol + y_off, -3.5*YbyCol + y_off, -4*YbyCol + y_off, -4.5*YbyCol + y_off, -5*YbyCol + y_off, -6*YbyCol + y_off, -7*YbyCol + y_off, 1e36};
const UInt_t NSieveCol = 27; 

const Double_t XbyRow = .46 * 25.4e-3;


const Double_t SieveXbyRow[] = {-4 * XbyRow, -3.5 * XbyRow, -3 * XbyRow, -2.5 * XbyRow, -2 * XbyRow, -1.5 * XbyRow, -1 * XbyRow, -0.5 * XbyRow, 0*XbyRow, 0.5 * XbyRow, 1*XbyRow, 1.5 * XbyRow, 2 * XbyRow, 2.5 * XbyRow, 3 * XbyRow, 3.5 * XbyRow, 4 * XbyRow, 1e36}; // vertical positions of the sieve holes when the column number is odd, column number starts with 0
const UInt_t NSieveRow = 17; 

//////Sieve Survey Inputs////
double yaw = -5.372 * D2R;     //Degree of yaw
double pitch = 90.019 * D2R;  //Degree of pitch 

// SieveOff* are in TCS
const Double_t SieveOffY = 2.07*1.e-3;// -(3.314-0.8)*1.e-3;
const Double_t SieveOffX = -1.65*1.e-3;//-(1.07+1.42)*1.e-3;
const Double_t SieveOffZ = 1.49*1.e-3;//-(1.07+1.42)*1.e-3;

const Double_t ZPos = 31.23 * 25.4e-3;



/////////////////////////////////////////////////////////////////////////
// Vertex Position Inputs

//static const UInt_t NFoils = 9; 
static const UInt_t NFoils = 3; 
const Double_t targetoffset = -5*1e-3;
//const Double_t targetfoils[] = {-10e-2+targetoffset, targetoffset-7.5e-2, targetoffset-5e-2, targetoffset-2.e-2, targetoffset+0.0, targetoffset+2.e-2, targetoffset+5e-2, targetoffset+7.5e-2,10e-2+targetoffset, 1e36};
const Double_t targetfoils[] = {-0.2 + targetoffset, 0.0 + targetoffset, 0.2 + targetoffset, 1e36};

 ///////////////////////////////////////////////////////////////////////// 
 // Excitation State Inputs 
const UInt_t NKine = 5; //N Delta Scans */

#define DIPOLE_MAG2MOM(Mag) (2.702*(Mag)-1.6e-03*(Mag)*(Mag)*(Mag)) 

const Double_t Ebeam[] = { 
  2.222, // -4% run 22771
  2.222, // -2% run 22772
  2.222, //  0% run 22775
  2.222, //  2% run 22776 
  2.222, //  4% run 22778
  0 
  }; 
const Double_t HRSCentralMom[] = { 
  1.27701, // -4% run 22771
  1.25303, // -2% run 22772
  1.22827, //  0% run 22775
  1.20423, //  2% run 22776 
  1.17922, //  4% run 22778
  0 
  }; 

const Double_t GroundNuclearMass = 1*0.938272046 -.511e-3*1; //GeV/c^2  //H2 Target
const Double_t ExcitationEnergy[] = {0.,0.,0.,0.,0.};//selected excitation states for each kinematics
//{0.,0.00443891,0.00443891,0.00443891,0.00443891};

const UInt_t NExcitationStates = 1; // C Excitation States
const Double_t ExcitationEnergyList[] = {0};

// gmp elastic delta scan goes through LH2 target, the AllLossExceptFoil doesn't consider the radiation loss in LH2 target
const Double_t AllLossExceptFoil
        = 1e-3*(//in MeV //rho*dE/dX*X, X in cm
	  2.78*2.035*16*1.e-3*2.54 //Al 2024, rho 2.78, scattering exit window Al
	  +1.205E-03*2.6922*14.79*2.54//air, Target Enclosure to HRS
	  +1.42*2.138*12*1.e-3*2.54 //kapton window on spectrometer entrance
        ); //AllLossExceptFoil==1e-3*0.444;

const Double_t LossEntranceWindow = 1e-3*2.81*2.0795*LH2_Thickness_Entance*100;//eloss at entrance window
const Double_t LossEachUnitB = 0.0723*4.7516*100*1e-3; // Radiation loss in 1m LH2 before scattering, 2.3004GeV
const Double_t LossEachUnitA = 0.0723*4.66154*100*1e-3; // Radiation loss in 1m LH2 after scattering, 1.25397
const Double_t LossEachUnitA_Al7075 = 2.81*2.036*100*1e-3; // Al Eloss in 1m Al 7075
const Double_t LossEachUnitA_Al2024 = 2.78*2.036*100*1e-3; // Al Eloss in 1m Al 2024


/////////////////////////////////////////////////////////////////////////
// Disable Extended Target Correction

//const Double_t ExtTarCor_ThetaCorr = 0.61;//0.00;//
//const Double_t ExtTarCor_DeltaCorr = 5.18;//1e36;//

const Double_t ExtTarCor_ThetaCorr = 0.00;
const Double_t ExtTarCor_DeltaCorr = 1e36;

/////////////////////////////////////////////////////////////////////////
// Database header

const char * Prefix = "R.vdc.";

const char * DatabaseHeader = "\
[ R.global ]   \n\
0.3327 1 0.0 270.2 0.0 -1.6e-03        VDC Angle, Plane Spacing, Gamma Coefficents  \n\
matrix elements   \n\
t 0 0 0  -1.003100e+00 -3.500200e-01 -4.878100e-02  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0\n\
y 0 0 0  -6.725700e-03 -1.141500e-02 -3.353800e-03  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0\n\
p 0 0 0  -2.492000e-03 -2.014700e-02 -1.291400e-02  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0\n\
";

const char * DatabaseFooter = "\
R 0 0 0 0  25.713\n\
R 1 0 0 0  0.1650\n\
R 2 0 0 0 -0.05\n\
R 0 1 0 0 -11.6554\n\
R 0 2 0 0 -9.4951\n\
R 0 0 1 0  0.0\n\
R 0 0 2 0  0.0\n\
R 0 0 0 1  0.0\n\
R 0 0 0 2  0.0\n\
XF 1 0 0 0 0 -2.181E+00\n\
XF 0 1 0 0 0 -1.980E-01\n\
XF 0 0 0 0 1  1.191E+01\n\
TF 1 0 0 0 0 -1.000E-01\n\
TF 0 1 0 0 0 -4.690E-01\n\
TF 0 0 0 0 1  1.967E+00\n\
PF 0 0 1 0 0  3.630E-01\n\
PF 0 0 0 1 0 -0.902E+00\n\
YF 0 0 1 0 0 -5.950E-01\n\
YF 0 0 0 1 0 -1.274E+00\n\
";
/////////////////////////////////////////////////////////////////////////

#endif
