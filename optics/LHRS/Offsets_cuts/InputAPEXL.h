// Description of GMP optics

#ifndef ROOT_Input
#define ROOT_Input

#include "TROOT.h"
#include "TMath.h"
#include "TVector3.h"
#include "TCut.h"
#include <map>
#include <iostream>


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
//const Double_t HRSAngle = -5 * D2R; 
//const Double_t HRSAngle = 5.366 * D2R;
const Double_t HRSAngle = 5.0 * D2R; 

//LH2 target information
const Double_t LH2_TargetLength = 15*1.e-2; //unit m
const Double_t LH2_Target_Tip_Radius = 1.5*2.54*1.e-2; //target width and tip radius
const Double_t LH2_Thickness_Entance = 0.175*1.e-3;//Al 7075, aluminum thickness for the entrance window
const Double_t LH2_Thickness_Side = 0.18*1.e-3;  //Al 7075, aluminum thickness for the side wall
const Double_t LH2_Thickness_Tip = 0.11*1.e-3;  //Al 7075, aluminum thickness for tip

// MissPoint* are in HCS
/* const Double_t MissPointZ =0;//  */
/* const Double_t MissPointY = 0;// */

// APEX LHRS // Hole_pos = [0.010284,-0.00305414,0.795905]
// Hole_pos = [0.010284,0.00659786,0.795905

/* const Double_t MissPointZ = 1.690e-3; */
/* const Double_t MissPointY = -1.790e-3; */

const Double_t MissPointZ = 0;
const Double_t MissPointY = -0;


const Double_t BeamShiftX = 0;

const Double_t SieveRadius = 0.055*25.4/2.0*1e-3;
const Double_t SieveRadius_c = 0.106*25.4/2.0*1e-3;

//const Double_t SieveRadius = 0;       

//const Double_t SieveRadius_c = 0;   

// average beam positions

// case of Optics Foils + vertical foils

Double_t BeamX_average[] = {-0.0006391,-0.000636,-0.0006391,-0.000636,-0.0006391,-0.000636,-0.0006391,-0.000636,0.002064,-0.0005847,-0.002499};

Double_t BeamY_average[] = {0.002405,0.002419,0.002405,0.002419,0.002405,0.002419,0.002405,0.002419,0.002284,0.002424,0.002498};

Double_t BeamXDir_average[] = {0.0002036,0.0001817,0.0002036,0.0001817,0.0002036,0.0001817,0.0002036,0.0001817,0.002392,0.0002429,-0.001303};
Double_t BeamYDir_average[] = {-0.0005503,-0.0005584,-0.0005503,-0.0005584,-0.0005503,-0.0005584,-0.0005503,-0.0005584,-0.0006598,-0.000541,-0.0004732};

// Double_t BeamXDir_average[] = {0};
// Double_t BeamYDir_average[] = {0};

// Double_t BeamX_average[] = {0};

// Double_t BeamY_average[] = {0};
//const Double_t BeamY_average = 0.008;


// average beam positions for Optics targets runs 4771 & 4774
//4771
// BeamX_average = -0.0006361;
// BeamY_average = 0.002419; 
// 4774
/* const Double_t BeamX_average = -0.0006391; */
/* const Double_t BeamY_average = 0.002405; */


// run 4766
/* const Double_t BeamX_average = -0.00268; */
/* const Double_t BeamY_average = 0.00251; */

// run 4768
// const Double_t BeamX_average = -0.000585;  //-0.00268;
// const Double_t BeamY_average = 0.002425;

 
// average beam directions used as test for correction

// const Double_t BeamXDir_average =   0.0002428; //0.0001816;
// const Double_t BeamYDir_average = -0.0005409;
const Double_t BeamZDir_average = 5.131;


/* const Double_t BeamXDir_average = 0.0; */
/* const Double_t BeamYDir_average = 0; */
/* const Double_t BeamZDir_average = 1; */


// create map for runnumber and average beam position
// First int is run_number next is beamz and beamy 

std::map<int, std::pair<int,int>> Beam_info;






/////////////////////////////////////////////////////////////////////////
// Sieve Position Inputs
const Double_t YbyCol = .19 * 25.4e-3;
// 16/1/20
const Double_t y_off = 31.23 * tan(0.8 *D2R) * 25.4e-3; 
//const Double_t y_off = 0.0; 
//const Double_t SieveYbyCol[]= {7*YbyCol + y_off, 6.5*YbyCol + y_off, 6*YbyCol + y_off, 5.5*YbyCol + y_off, 5*YbyCol + y_off, 4.5*YbyCol + y_off, 4*YbyCol + y_off, 3.5*YbyCol + y_off, 3*YbyCol + y_off, 2.5*YbyCol + y_off, 2*YbyCol + y_off, 1.5*YbyCol + y_off, 1.0*YbyCol + y_off, 0.5*YbyCol + y_off, 0.0*YbyCol + y_off, -0.5*YbyCol + y_off, -1*YbyCol + y_off, -1.5*YbyCol + y_off, -2*YbyCol + y_off, -2.5*YbyCol + y_off, -3*YbyCol + y_off, -3.5*YbyCol + y_off, -4*YbyCol + y_off, -4.5*YbyCol + y_off, -5*YbyCol + y_off, -6*YbyCol + y_off, -7*YbyCol + y_off, 1e36};


/* const Double_t SieveYbyCol[]= {7*YbyCol + y_off, 6*YbyCol + y_off, 5*YbyCol + y_off, 4.5*YbyCol + y_off, 4*YbyCol + y_off, 3.5*YbyCol + y_off, 3*YbyCol + y_off, 2.5*YbyCol + y_off, 2*YbyCol + y_off, 1.5*YbyCol + y_off, 1.0*YbyCol + y_off, 0.5*YbyCol + y_off, 0.0*YbyCol + y_off, -0.5*YbyCol + y_off, -1*YbyCol + y_off, -1.5*YbyCol + y_off, -2*YbyCol + y_off, -2.5*YbyCol + y_off, -3*YbyCol + y_off, -3.5*YbyCol + y_off, -4*YbyCol + y_off, -4.5*YbyCol + y_off, -5*YbyCol + y_off, -5.5*YbyCol + y_off, -6*YbyCol + y_off, -6.5*YbyCol + y_off, -7*YbyCol + y_off, 1e36}; */
const Double_t SieveYbyCol[]= { -7*YbyCol + y_off, -6.5*YbyCol + y_off, -6*YbyCol + y_off, -5.5*YbyCol + y_off, -5*YbyCol + y_off, -4.5*YbyCol + y_off, -4*YbyCol + y_off, -3.5*YbyCol + y_off, -3*YbyCol + y_off, -2.5*YbyCol + y_off, -2*YbyCol + y_off, -1.5*YbyCol + y_off, -1*YbyCol + y_off, -0.5*YbyCol + y_off, 0.0*YbyCol + y_off, 0.5*YbyCol + y_off, 1.0*YbyCol + y_off, 1.5*YbyCol + y_off, 2*YbyCol + y_off, 2.5*YbyCol + y_off, 3*YbyCol + y_off, 3.5*YbyCol + y_off, 4*YbyCol + y_off, 4.5*YbyCol + y_off, 5*YbyCol + y_off, 6*YbyCol + y_off, 7*YbyCol + y_off, 1e36}; 
const UInt_t NSieveCol = 27; 

const Double_t XbyRow = .46 * 25.4e-3;


const Double_t SieveXbyRow[] = {-4 * XbyRow, -3.5 * XbyRow, -3 * XbyRow, -2.5 * XbyRow, -2 * XbyRow, -1.5 * XbyRow, -1 * XbyRow, -0.5 * XbyRow, 0*XbyRow, 0.5 * XbyRow, 1*XbyRow, 1.5 * XbyRow, 2 * XbyRow, 2.5 * XbyRow, 3 * XbyRow, 3.5 * XbyRow, 4 * XbyRow, 1e36}; // vertical positions of the sieve holes when the column number is odd, column number starts with 0
const UInt_t NSieveRow = 17; 


/*JW: for APEX

                      Z       X        Y     Angle
B.L. Sieve Centre     791.8   81.0    1.4    5.366
*/

// SieveOff* are in TCS
//const Double_t SieveOffY = 0;
const Double_t SieveOffY = 0.665315e-3; // obtained from technical drawing of set-up (31.23 inches from V2 foil to centre of sieve slit, 0.8 degree angle from centre of sieve slit to largest hole, last number is conversion to metres)
/* const Double_t SieveOffX = 0; */
const Double_t SieveOffX = -1.334e-3; // X points down in TCS, hence negative of Y term in 'Hall-like' system where Survey results are taken from

const Double_t SieveOffZ = 2.6756e-3;

// experiment with SieveOff in HCS and using fTCSinHCS to convert between

/* const Double_t SieveOffX_HCS = 81.0e-3; */
/* const Double_t SieveOffY_HCS = 1.4e-3; */


//Double_t foil_d = 0;
//const Double_t ZPos_HCS = 791.8e-3;
/* const Double_t ZPos_HCS = 791.8e-3; */

//const Double_t ZPos =1059.61e-3+3.314e-3/TMath::Tan(-HRSAngle);//1059.61 * 1e-3;
const Double_t ZPos = 31.23 * 25.4e-3;

/////////////////////////////////////////////////////////////////////////
// Vertex Position Inputs

//static const UInt_t NFoils = 8; 
/* static const UInt_t NFoils = 3;  */
const Double_t targetoffset = 0;
//const Double_t targetfoils[] = {-10e-2+targetoffset, targetoffset-7.5e-2, targetoffset-5e-2, targetoffset-2.e-2, targetoffset+0.0, targetoffset+2.e-2, targetoffset+5e-2, targetoffset+7.5e-2,10e-2+targetoffset, 1e36};
//const Double_t targetfoils[] = {-0.2, 0.0, 0.2, 1e36};

//const Double_t targetfoils[] = { -0.188, 0.012, 0.212 }; 


/* static const UInt_t NFoils = 8; */
//const Double_t targetfoils[] = { -0.288,  -0.207, -0.138,  -0.063, 0.087, 0.162, 0.231, 0.312 };

//const Double_t targetfoils[] = { -0.3,  -0.219, -0.15,  -0.075, 0.075, 0.15, 0.219, 0.3};


// scheme where foils Optics foils are numbered 0-7 and Vertical foils 8-10

static const UInt_t NFoils = 11;
// z-corrections for foils from survey

// general foil z-correction to nearest mm

static const Double_t foil_zcorr = -5e-3;
// static const Double_t foil_zcorr = 0;

// W_wires correction from survey
static const Double_t foil_zcorr_W = -5.01e-3;

// Optics_up
static const Double_t foil_zcorr_up = -5.06e-3;

// Optics_middle/down (middle and down have same offset in survey)
static const Double_t foil_zcorr_down = -5.01e-3;


// x-corrections for foils from survey

// W_wires correction from survey
static const Double_t foil_xcorr_W = -0.01e-3;

// Optics_up
static const Double_t foil_xcorr_up = -0.07e-3;

// Optics_middle/down (middle and down have same offset in survey)
static const Double_t foil_xcorr_down = -0.04e-3;


//const Double_t targetfoils[] = { -0.3,  -0.219, -0.15,  -0.075, 0.075, 0.15, 0.219, 0.3, -0.2, 0, 0.2, 1e36};

// foil z-positions
const Double_t targetfoils[] = { -0.3 + foil_zcorr,  -0.219 + foil_zcorr, -0.15 + foil_zcorr,  -0.075 + foil_zcorr, 0.075 + foil_zcorr, 0.15 + foil_zcorr, 0.219 + foil_zcorr, 0.3 + foil_zcorr, -0.2 + foil_zcorr, 0 + foil_zcorr, 0.2 + foil_zcorr, 1e36};

// foil x-positions
const Double_t orig_targetfoilsX[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0.0025 + foil_xcorr_W, 0 + foil_xcorr_W, -0.0025 + foil_xcorr_W, 1e36};

const Double_t targetfoilsX[] = { 0, 0, 0, 0, 0, 0, 0, 0, 2.124e-3, -0.019e-3, -2.162e-3, 1e36};






std::vector<TString> Foil_names;

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



TVector3 fPointingOffset;

const UInt_t NoinEachRow[] = {15, 12, 15, 11, 15, 11, 15, 11, 15, 11, 15, 11, 15, 11, 15, 12, 15};


const Int_t NHoles = 225; 



////////////////////////////////////////////////
// PID and general track cuts
///////////////////////////////////////////////



TCut GeneralSieveCut = "L.tr.n==1 && L.tr.chi2<0.003 && L.vdc.u1.nclust==1 && L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.v2.nclust==1";

// TCut GeneralSieveCut = "L.tr.n==1 && L.tr.chi2<0.003 && abs(L.gold.th)<0.08 && L.gold.ph>-0.07 && L.gold.ph<0.025 && abs(L.tr.r_x)<0.5 && L.vdc.u1.nclust==1 && L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.v2.nclust==1"; 


TCut FP_cuts = "L.tr.r_x>-0.6 && L.tr.r_x<0.56 && L.tr.r_y>-0.04 && L.tr.r_y<0.037 && L.tr.r_th>-0.029 && L.tr.r_th<0.02 && L.tr.r_ph>-0.05 && L.tr.r_ph<0.04 && abs(L.tr.r_x)<0.1";

// TCut FP_cuts = "L.tr.r_x>-0.6 && L.tr.r_x<0.56 && L.tr.r_y>-0.04 && L.tr.r_y<0.037 && L.tr.r_th>-0.029 && L.tr.r_th<0.02 && L.tr.r_ph>-0.05 && L.tr.r_ph<0.04";

TCut PID_cuts = "(L.prl1.e/(L.gold.p*1000))>0.3 && ((L.prl1.e+L.prl2.e)/(L.gold.p*1000))>0.625 && ((L.prl1.e+L.prl2.e)/(L.gold.p*1000))<1.11 &&  L.cer.asum_c >650 && DL.evtype == 1";

// && Lrb.x>-0.0029 && Lrb.x<-0.0022";


// list of runs for LHRS optics based on target


std::vector<int> V1_runs{4179,4766,4767};

std::vector<int> V2_runs{4181,4768};

std::vector<int> V3_runs{4180,4769,4770};

std::vector<int> Opt1_runs{4771,4772};

std::vector<int> Opt3_runs{4773,4774};

std::vector<int> horizontal_runs{4775,4776,4777};

std::vector<TString> Single_foil{"V1","V2","V3"};
std::vector<TString> Multi_foil{"Optics1","Optics3","horizontal"};



Double_t dp_split[6][2] = {{-0.05,-0.03},{-0.03,-0.01},{-0.01,0.01},{0.01,0.03},{0.03,0.05},{0,0}};



//////Sieve Survey Inputs////
double yaw = 5.366 * D2R;     //Abs value of yaw
double pitch = 89.988 * D2R;  //Degree of pitch 


double target_yaw =  0.1022 * D2R;
double target_pitch = -0.0219 * D2R;  //Degree of pitch 


// function that gets beam cut from run number
TString Return_target(Int_t runnumber);

TString get_Beamcut(Int_t runnumber){

  TString beam_cut = "";

  if(runnumber == 4766){
    // V1
    cout << "V1 wire beam cut!" << endl;
    beam_cut = "0.00170<Lrb.x && Lrb.x<0.0024";
    //    beam_cut = "0.00175<Lrb.x && Lrb.x<0.0021";
  }
  else if(runnumber == 4768){
    // V2
    cout << "V2 wire beam cut!" << endl;
    beam_cut = "-0.001<Lrb.x && Lrb.x<-0.0003";
    // beam_cut = "-0.0008<Lrb.x && Lrb.x<-0.0004";
  }
  else if(runnumber == 4769){
    // V3
    cout << "V3 wire beam cut!" << endl;
    beam_cut = "-0.0029<Lrb.x && Lrb.x<-0.00225";
    // beam_cut = "-0.0027<Lrb.x && Lrb.x<-0.0023";
  }
  else if( !Return_target(runnumber).CompareTo(TString("horizontal"))){
    cout << "Horizontal beam cut!" << endl;
    beam_cut = "(0.00191<Lrb.y && Lrb.y<0.00209)"; 
      //" (0.00240<Lrb.y && Lrb.y<0.00258)";
      //"(0.00191<Lrb.y && Lrb.y<0.00209)"; 
    //||
		   //    (0.00191<Lrb.y && Lrb.y<0.00209)
      //      && //(0.00175<Lrb.y&&Lrb.y<0.00215)";
 //(0.0022<Lrb.y&&Lrb.y<0.0027)";
      
  }
  

  return beam_cut;
}




// functions that can return run type (Foil) from run number

bool Contains(const std::vector<int> &list, int x){

  return std::find(list.begin(), list.end(), x) != list.end();

}

bool Contains(const std::vector<TString> &list, TString x){

  return std::find(list.begin(), list.end(), x) != list.end();

}



TString Return_target(Int_t runnumber){

  TString run_type; 


  if( Contains(V1_runs,runnumber)){
    run_type = "V1";
  }
  else if( Contains(V2_runs,runnumber)){
    run_type = "V2";
  }
  else if( Contains(V3_runs,runnumber)){
    run_type = "V3";
  }
  else if( Contains(Opt1_runs,runnumber)){
    run_type = "Optics1";
  }
  else if( Contains(Opt3_runs,runnumber)){
    run_type = "Optics3";
  }
  else if( Contains(horizontal_runs,runnumber)){
    run_type = "horizontal";
  }

  
    
  return run_type;
  
}

// function returns Foil number based on run number
// for multifoil runs it returns the first foil number (ie 0 for Optics3 and 1 for Optics1)

Int_t GetFoilID(Int_t runnumber){

  Int_t FoilID = -1;
  
  TString run_type = Return_target(runnumber);

  if(run_type == "V1"){
    FoilID = 8;
  }
  else if( run_type == "V2"){
    FoilID = 9;
  }
  else if( run_type == "V3"){
    FoilID = 10;
  }
  else if( run_type == "Optics3"){
    FoilID = 0;    
  }
  else if( run_type == "Optics1"){
    FoilID = 1;    
  }
  else{
    cout << "Run " << runnumber << " does not have recognised type" << endl;
  }
  //    cout << "Not a single foil run" << endl;
  
    
  return FoilID;
}


Bool_t IsMultiFoil(Int_t runnumber){

  Bool_t ismultifoil;
  
  TString run_type = Return_target(runnumber);

  if( Contains(Multi_foil,run_type)){
    ismultifoil = true;
  }
  else if(Contains(Single_foil,run_type)){
    ismultifoil = false;
  }
  else{
    cout << "Chosen run_number not recognised as single or multifoil: defaulted to single foil" << endl;
    ismultifoil = true;
  }
    

  return ismultifoil;
  


}




#endif











