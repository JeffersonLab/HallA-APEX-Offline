///////////////////////////////////////////////////////////////////////////////
//
// LOpticsOpt
//
// HRS optics matrix optimization class
// Based on THaVDC
//
// Units used:
//   For X, Y, and Z coordinates of track    -  meters
//   For Theta and Phi angles of track       -  tan(angle)
//   For Momentums, Masses                   -  GeV, GeV/c^2
//
// Author: Jin Huang <jinhuang@jlab.org>
//
// Modification:
//   Jun 25, 2010 Updated for APEX optics calibration
//   Aug 01, 2013 Updated for G2P optics calibration (Chao Gu)
//
///////////////////////////////////////////////////////////////////////////////


#include <cstdio>
#include <cstdlib>
#include <map>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TDatime.h"
#include "TGraphErrors.h"
#include "TClonesArray.h"
#include "TList.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"
#include "TLatex.h"
#include "TVector3.h"
#include "TLine.h"
#include "TArrow.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TGTableHeader.h"
#include "TGTable.h"
#include "TGWidget.h"
#include "TVirtualTableInterface.h"
#include "TGSimpleTableInterface.h"
#include "TPaveText.h"
#include "TText.h"
#include "TAttLine.h"
#include "TLegend.h"

#include "THaGlobals.h"
#include "THaEvData.h"
#include "THaDetMap.h"
#include "THaTrack.h"
#include "THaScintillator.h"
#include "THaSpectrometer.h"
#include "VarDef.h"

#include "LOpticsOpt.h"

//#include "APEX_Sieve.h"

#ifdef WITH_DEBUG
#include <iostream>
#endif

#define tgy_hole_info false

using namespace std;
using THaString::Split;

///////////////////////////////////////////////////////////////////////////////
// Input Sections
///////////////////////////////////////////////////////////////////////////////

#include "InputAPEXL.h"

///////////////////////////////////////////////////////////////////////////////
// Constructors
///////////////////////////////////////////////////////////////////////////////

LOpticsOpt::LOpticsOpt(const char* name, const char* description, THaApparatus* apparatus) :
THaTrackingDetector(name, description, apparatus)
{
    fPrefix = new char[1000];
    sprintf(fPrefix, "%s", Prefix);

    fCurrentMatrixElems = NULL;

    TVector3 TCSX(0, -1, 0);
    TVector3 TCSZ(TMath::Sin(HRSAngle), 0, TMath::Cos(HRSAngle));
    TVector3 TCSY = TCSZ.Cross(TCSX);
    fTCSInHCS.RotateAxes(TCSX, TCSY, TCSZ);
 
    const Int_t a = (HRSAngle > 0) ? 1 : -1;
    //fPointingOffset, HCS
    fPointingOffset.SetXYZ(-a*MissPointZ*TMath::Cos(HRSAngle), MissPointY, -MissPointZ * TMath::Sin(HRSAngle)); 
    //for left arm, x misspoint is -MissPointZ*cos(theta), for r-arm, x misspoint is MissPointZ*cos(theta) 


    DEBUG_INFO("LOpticsOpt", "Read in configuration " + InputID);
    DEBUG_INFO("LOpticsOpt", "HRS @ %f Degree, PointingOffset = (%f,%f,%f), SievePos = (%f,%f,%f)", HRSAngle / TMath::Pi()*180, fPointingOffset.X(), fPointingOffset.Y(), fPointingOffset.Z(), SieveOffX, SieveOffY, ZPos);

    fNRawData = 0;

    for (UInt_t i = 0; i < 100; i++)
        fArbitaryDpKinShift[i] = fArbitaryVertexShift[i] = 0;


    // give names to foils used
    
    // vertical foils
    if(NFoils == 3){     
      for(UInt_t i = 0; i<NFoils; i++){
	TString foil(Form("V%d",i+1));
	Foil_names.push_back(foil);
      }      
    }
    // Optics foils
    else if(NFoils == 8){
      for(UInt_t i = 0; i<NFoils; i++){
	TString foil(Form("Opt %d",i+1));
	Foil_names.push_back(foil);
      }      
    }
    // Optics and Vertical foils (0-7 wires are optics, 8-10 are Vertical)
    else if(NFoils == 11){
      UInt_t i = 0;
      for(; i<8; i++){
	TString foil(Form("Opt %d",i+1));
	Foil_names.push_back(foil);
      }      
      for(; i<11; i++){
	TString foil(Form("V%d",i-7));
	Foil_names.push_back(foil);
      }      
    }


}

LOpticsOpt::~LOpticsOpt()
{
    // Destructor.
}

///////////////////////////////////////////////////////////////////////////////
// Database I/O
///////////////////////////////////////////////////////////////////////////////
Double_t LOpticsOpt::TravelLength(TVector3 ReactionVertex, TVector3 MomDirectionHCS){

  //code to calculate the travel length in the LH2
  //need to pass the vertex point postion, the scattered electron momentum direction, the beam position
  Double_t pi = 3.14159265358979323846;
  Double_t TL = LH2_TargetLength; //Unit m, Target Length
  Double_t R = LH2_Target_Tip_Radius; //target width and tip radius 

//assume we know the vetex point, beam position X, beam position Y and CalcReacZ

  Double_t vx=ReactionVertex.X(),vy=ReactionVertex.Y(),vz=ReactionVertex.Z(); // vertex position
  Double_t tx=0,ty=0,tz=TL/2.-R;//center postion of the tip
  TVector3 v_tipcenter; //vector from vertex to tip center, the direction is important
  Double_t distance_tip=0;//distance between the vertex and the tip center
  
  //we know the momentum direction from MomDirectionHCS
  Double_t thp=0,php=0; //thp=dx/dz,php=dy/dz
  Double_t px=0,py=0,pz=0;

  //  vx=-4E-3;
  //  vy=-4E-3;
  //  MomDirectionHCS(0)=-1;
  //  MomDirectionHCS(1)=-2;
  //  MomDirectionHCS(2)=TMath::Sqrt(7);
  
  distance_tip = TMath::Sqrt(TMath::Power(vx-tx,2)+TMath::Power(vy-ty,2)+TMath::Power(vz-tz,2));
  px= MomDirectionHCS(0);
  py= MomDirectionHCS(1);
  pz= MomDirectionHCS(2);


  v_tipcenter(0)=tx-vx;
  v_tipcenter(1)=ty-vy;
  v_tipcenter(2)=tz-vz;
  thp=MomDirectionHCS.X()/MomDirectionHCS.Z();
  php=MomDirectionHCS.Y()/MomDirectionHCS.Z(); //thp=dx/dz,php=dy/dz


  Double_t angleR=0,sin_angleR,cos_angleR=0; 
  //angleR is the angle between the scattering momentum and the vertex-tipcenter line
  cos_angleR =  (v_tipcenter.Dot(MomDirectionHCS))/( v_tipcenter.Mag()*MomDirectionHCS.Mag());
  sin_angleR=TMath::Sqrt(1-cos_angleR*cos_angleR);
  angleR = TMath::ACos( cos_angleR );

  //the center of the target is located at (0,0,0)
  Double_t anglex=0,sin_anglex=0;
  Double_t D=0,angleD=0;
  Double_t x,y,z;//the calculated position at the cylindral
  Double_t m1,m2,a,b,c,result1,result2;
  Double_t sin_theta=0;
 
  a=(px*px/py/py+1);
  m1=vx/R;
  m2=px/py/R*vy;
  b=2*px/py*(m1-m2);
  c=(m1-m2)*(m1-m2)-1;
  result1 = (-b+TMath::Sqrt(b*b-4*a*c))/2/a; 
  result2=(-b-TMath::Sqrt(b*b-4*a*c))/2/a;

  sin_theta = (R*result1-vy)*py>=0?result1:sin_theta ;
  sin_theta = (R*result2-vy)*py>=0?result2:sin_theta ;

  y=R*sin_theta;
  x=px/py*y-px/py*vy+vx;
  z=pz/py*(y-vy)+vz;

  if( z<(TL/2. -R) )
    {
      D = TMath::Sqrt((x-vx)*(x-vx)+(y-vy)*(y-vy)+(z-vz)*(z-vz));

    }

  else
    {
      sin_anglex=distance_tip/R*sin(angleR);
      anglex = TMath::ASin(sin_anglex);
      angleD = pi - angleR - anglex;
      D=R/sin(angleR)*sin(angleD);
      
    }

  return D;
}

Double_t LOpticsOpt::ElossTargetBefore(TVector3 ReactionVertex, TVector3 MomDirectionHCS){

  //code to calculate the energy loss in the LH2, including the target shell

  Double_t ElossBefore=0;	
  Double_t TL = LH2_TargetLength; //Unit m, Target Length
  Double_t vz=ReactionVertex.Z(); // vertex position
  ElossBefore = (vz + TL/2.) * LossEachUnitB + LossEntranceWindow; //1e-3*2.699*2.074*0.0175, LH2 Aluminum cell at entrance ; 

  return ElossBefore;
}

Double_t LOpticsOpt::ElossTargetAfter(TVector3 ReactionVertex, TVector3 MomDirectionHCS){

  //code to calculate the energy loss in the LH2, including the target shell

  Double_t ElossAfter=0;	

  Double_t pi = 3.14159265358979323846;
  Double_t TL = LH2_TargetLength; //Unit m, Target Length
  Double_t R = LH2_Target_Tip_Radius; //target width and tip radius 
  Double_t vx=ReactionVertex.X(),vy=ReactionVertex.Y(),vz=ReactionVertex.Z(); // vertex position
  Double_t tx=0,ty=0,tz=TL/2.-R;//center postion of the tip
  TVector3 v_tipcenter; //vector from vertex to tip center, the direction is important
  Double_t distance_tip=0;//distance between the vertex and the tip center
  
  //we know the momentum direction from MomDirectionHCS
  Double_t thp=0,php=0; //thp=dx/dz,php=dy/dz
  Double_t px,py,pz;

  distance_tip = TMath::Sqrt(TMath::Power(vx-tx,2)+TMath::Power(vy-ty,2)+TMath::Power(vz-tz,2));
  px= MomDirectionHCS(0);
  py= MomDirectionHCS(1);
  pz= MomDirectionHCS(2);


  v_tipcenter(0)=tx-vx;
  v_tipcenter(1)=ty-vy;
  v_tipcenter(2)=tz-vz;
  thp=MomDirectionHCS.X()/MomDirectionHCS.Z();
  php=MomDirectionHCS.Y()/MomDirectionHCS.Z(); //thp=dx/dz,php=dy/dz


  Double_t angleR=0,sin_angleR,cos_angleR=0; 
  //angleR is the angle between the scattering momentum and the vertex-tipcenter line
  cos_angleR =  (v_tipcenter.Dot(MomDirectionHCS))/( v_tipcenter.Mag()*MomDirectionHCS.Mag());
  sin_angleR=TMath::Sqrt(1-cos_angleR*cos_angleR);
  angleR = TMath::ACos( cos_angleR );

  //the center of the target is located at (0,0,0)
  Double_t anglex=0,sin_anglex=0;//anglex corresponds to distance_tip
  Double_t D=0,angleD=0;//angleD corresponds to the traveling length
  Double_t x,y,z;
  Double_t m1,m2,a,b,c,result1,result2;
  Double_t sin_theta=0;
 
  a=(px*px/py/py+1);
  m1=vx/R;
  m2=px/py/R*vy;
  b=2*px/py*(m1-m2);
  c=(m1-m2)*(m1-m2)-1;
  result1 = (-b+TMath::Sqrt(b*b-4*a*c))/2/a; 
  result2=(-b-TMath::Sqrt(b*b-4*a*c))/2/a;

  sin_theta = (R*result1-vy)*py>=0?result1:sin_theta ;
  sin_theta = (R*result2-vy)*py>=0?result2:sin_theta ;

  y=R*sin_theta;
  x=px/py*y-px/py*vy+vx;
  z=pz/py*(y-vy)+vz;

  if( z<(TL/2. -R) )
    {
      D = TMath::Sqrt((x-vx)*(x-vx)+(y-vy)*(y-vy)+(z-vz)*(z-vz));
      ElossAfter = D*LossEachUnitA+ LH2_Thickness_Side/TMath::Sqrt(thp*thp+php*php)*TMath::Sqrt(thp*thp+php*php+1)*LossEachUnitA_Al7075;

    }

  else
    {
      sin_anglex=distance_tip/R*sin(angleR);
      anglex = TMath::ASin(sin_anglex);
      angleD = pi - angleR - anglex;
      D=R/sin(angleR)*sin(angleD);
      ElossAfter = D*LossEachUnitA+LH2_Thickness_Tip/cos(anglex)*LossEachUnitA_Al7075;
      
    }

  return (ElossAfter + AllLossExceptFoil);
}

Int_t LOpticsOpt::LoadDataBase(TString DataBaseName)
{
    static const char* const here = "LoadDataBase";
    OldComments = "";

    FILE* file = fopen(DataBaseName, "r");
    if (!file) {
        Error("LoadDataBase", "%s can not be opened", DataBaseName.Data());
        assert(0); //
        return kFileError;
    } else DEBUG_INFO("LoadDataBase", "Parsing Database %s", DataBaseName.Data());

    const int LEN = 200;
    char buff[LEN];

    // Look for the section [<prefix>.global] in the file, e.g. [ R.global ]
    TString tag(fPrefix);
    Ssiz_t tagpos = tag.Index(".");
    if (tagpos != kNPOS)
        tag = tag(0, tagpos + 1);
    else
        tag.Append(".");
    tag.Prepend("[");
    tag.Append("global]");
    TString tag2(tag);
    tag.ToLower();

    bool found = false;
    while (!found && fgets(buff, LEN, file) != NULL) {
        // read in comments
        TString line = buff;
        if (line.BeginsWith("#")) {
            OldComments += line;
            // OldComments += "\n";
        }

        line = ::Compress(buff); // strip blanks
        if (line.EndsWith("\n")) line.Chop();

        line.ToLower();
        if (tag == line) found = true;
    }
    if (!found) {
        Error(Here(here), "Database entry %s not found!", tag2.Data());
        fclose(file);
        assert(0); //
        return kInitError;
    }

    // We found the section, now read the data
    fgets(buff, LEN, file); // Skip constant line
    fgets(buff, LEN, file); // Skip comment line

    fTMatrixElems.clear();
    fDMatrixElems.clear();
    fPMatrixElems.clear();
    fPTAMatrixElems.clear();
    fYMatrixElems.clear();
    fYTAMatrixElems.clear();
    fRMatrixElems.clear();

    fFPMatrixElems.clear();
    fFPMatrixElems.resize(3);

    typedef vector<string>::size_type vsiz_t;
    map<string, vsiz_t> power;
    power["t"] = 3; // transport to focal-plane tensors
    power["y"] = 3;
    power["p"] = 3;
    power["D"] = 3; // focal-plane to target tensors
    power["T"] = 3;
    power["Y"] = 3;
    power["YTA"] = 4;
    power["P"] = 3;
    power["PTA"] = 4;
    power["R"] = 4; // pathlength from z=0 (target) to focal plane (meters)
    power["XF"] = 5; // forward: target to focal-plane (I think)
    power["TF"] = 5;
    power["PF"] = 5;
    power["YF"] = 5;

    map<string, vector<THaMatrixElement>*> matrix_map;
    matrix_map["t"] = &fFPMatrixElems;
    matrix_map["y"] = &fFPMatrixElems;
    matrix_map["p"] = &fFPMatrixElems;
    matrix_map["D"] = &fDMatrixElems;
    matrix_map["T"] = &fTMatrixElems;
    matrix_map["Y"] = &fYMatrixElems;
    matrix_map["YTA"] = &fYTAMatrixElems;
    matrix_map["P"] = &fPMatrixElems;
    matrix_map["PTA"] = &fPTAMatrixElems;
    matrix_map["R"] = &fRMatrixElems;

    map<string, int> fp_map;
    fp_map["t"] = 0;
    fp_map["y"] = 1;
    fp_map["p"] = 2;

    // Read in as many of the matrix elements as there are.
    // Read in line-by-line, so as to be able to handle tensors of
    // different orders.
    while (fgets(buff, LEN, file)) {
        string line(buff);
        // Erase trailing newline
        if (line.size() > 0 && line[line.size() - 1] == '\n') {
            buff[line.size() - 1] = 0;
            line.erase(line.size() - 1, 1);
        }
        // Split the line into whitespace-separated fields
        vector<string> line_spl = Split(line);

        // Stop if the line does not start with a string referring to
        // a known type of matrix element. In particular, this will
        // stop on a subsequent timestamp or configuration tag starting with "["
        if (line_spl.empty()) continue; //ignore empty lines
        const char* w = line_spl[0].c_str();
        vsiz_t npow = power[w];
        if (npow == 0) break;

#if DEBUG_LEVEL>=4
        cout << "Matrix Line = ";
        for (Ssiz_t i = 1; (UInt_t) i < (UInt_t) line_spl.size(); i++) {
            cout << i << "(" << line_spl[i].c_str() << "), ";
        }
        cout << endl;
#endif

        // Looks like a good line, go parse it.
        THaMatrixElement ME;
        ME.pw.resize(npow);
        ME.iszero = true;
        ME.order = 0;
        vsiz_t pos;
        for (pos = 1; pos <= npow && pos < line_spl.size(); pos++) {
            ME.pw[pos - 1] = atoi(line_spl[pos].c_str());
        }
        vsiz_t p_cnt;
        for (p_cnt = 0; pos < line_spl.size() && p_cnt < kPORDER && pos <= npow + kPORDER; pos++, p_cnt++) {
            ME.poly[p_cnt] = atof(line_spl[pos].c_str());
            if (ME.poly[p_cnt] != 0.0) {
                ME.iszero = false;
                ME.order = p_cnt + 1;
            }
        }
        if (p_cnt < 1) {
            Error(Here(here), "Could not read in Matrix Element %s%d%d%d!", w, ME.pw[0], ME.pw[1], ME.pw[2]);
            Error(Here(here), "Line looks like: %s", line.c_str());
            fclose(file);
            return kInitError;
        }

        // order optimize to
        ME.OptOrder = atoi(line_spl[line_spl.size() - 1].c_str());

        // Don't bother with all-zero matrix elements
        if (ME.iszero) continue;

        // Add this matrix element to the appropriate array
        vector<THaMatrixElement> *mat = matrix_map[w];
        if (mat) {
            // Special checks for focal plane matrix elements
            if (mat == &fFPMatrixElems) {
                if (ME.pw[0] == 0 && ME.pw[1] == 0 && ME.pw[2] == 0) {
                    THaMatrixElement& m = (*mat)[fp_map[w]];
                    if (m.order > 0) {
                        Warning(Here(here), "Duplicate definition of focal plane matrix element: %s. Using first definition.", buff);
                    } else
                        m = ME;
                } else
                    Warning(Here(here), "Bad coefficients of focal plane matrix element %s", buff);
            } else {
                // All other matrix elements are just appended to the respective array
                // but ensure that they are defined only once!
                bool match = false;
                for (vector<THaMatrixElement>::iterator it = mat->begin(); it != mat->end() && !(match = it->match(ME)); it++) {
                }
                if (match) {
                    Warning(Here(here), "Duplicate definition of matrix element: %s. Using first definition.", buff);
                } else
                    mat->push_back(ME);
            }
        } else if (fDebug > 0)
            Warning(Here(here), "Not storing matrix for: %s !", w);
    }

    CalcMatrix(1., fRMatrixElems); // tensor without explicit polynomial in x_fp

    fIsInit = true;
    fclose(file);
    return kOK;
}

Int_t LOpticsOpt::SaveDataBase(TString DataBaseName)
{
    // Output database in memory to new database file
    // WARNING: Hard coded text included

    DEBUG_INFO("SaveDataBase", "Saving to %s", DataBaseName.Data());

    typedef vector<THaMatrixElement>::size_type vsiz_t;

    FILE* file = fopen(DataBaseName, "w");
    if (!file) {
        Info("SaveDataBase", "Error Open in %s", DataBaseName.Data());
        return kFileError;
    }
    TDatime dt;
    fprintf(file, "# -------------------------------------------------------------");
    fprintf(file, "\n");
    fprintf(file, "# Optimized by Chao Gu @ %s", dt.AsString());
    fprintf(file, "\n");
    fprintf(file, "# Saved to %s", DataBaseName.Data());
    fprintf(file, "\n");
    fprintf(file, OldComments);

    fprintf(file, DatabaseHeader);

    DEBUG_INFO("SaveDataBase", "Transport Matrix:  D-terms");
    for (vsiz_t i = 0; i < fDMatrixElems.size(); i++) {
        fprintf(file, "D ");
        const THaMatrixElement& m = fDMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            fprintf(file, "%d ", m.pw[j]);
        }
        int j;
        for (j = 0; j < m.order; j++) {
            fprintf(file, " %13.6e", m.poly[j]);
        }
        for (; j < kPORDER; j++) {
            fprintf(file, " %13.6e", 0.0);
        }
        fprintf(file, "  %d", m.OptOrder);
        fprintf(file, "\n");
    }

    DEBUG_INFO("SaveDataBase", "Transport Matrix:  T-terms");
    for (vsiz_t i = 0; i < fTMatrixElems.size(); i++) {
        fprintf(file, "T ");
        const THaMatrixElement& m = fTMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            fprintf(file, "%d ", m.pw[j]);
        }
        int j;
        for (j = 0; j < m.order; j++) {
            fprintf(file, " %13.6e", m.poly[j]);
        }
        for (; j < kPORDER; j++) {
            fprintf(file, " %13.6e", 0.0);
        }
        fprintf(file, "  %d", m.OptOrder);
        fprintf(file, "\n");
    }

    DEBUG_INFO("SaveDataBase", "Transport Matrix:  P-terms");
    for (vsiz_t i = 0; i < fPMatrixElems.size(); i++) {
        fprintf(file, "P ");
        const THaMatrixElement& m = fPMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            fprintf(file, "%d ", m.pw[j]);
        }
        int j;
        for (j = 0; j < m.order; j++) {
            fprintf(file, " %13.6e", m.poly[j]);
        }
        for (; j < kPORDER; j++) {
            fprintf(file, " %13.6e", 0.0);
        }
        fprintf(file, "  %d", m.OptOrder);
        fprintf(file, "\n");
    }

    DEBUG_INFO("SaveDataBase", "Transport Matrix:  Y-terms");
    for (vsiz_t i = 0; i < fYMatrixElems.size(); i++) {
        fprintf(file, "Y ");
        const THaMatrixElement& m = fYMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            fprintf(file, "%d ", m.pw[j]);
        }
        int j;
        for (j = 0; j < m.order; j++) {
            fprintf(file, " %13.6e", m.poly[j]);
        }
        for (; j < kPORDER; j++) {
            fprintf(file, " %13.6e", 0.0);
        }
        fprintf(file, "  %d", m.OptOrder);
        fprintf(file, "\n");
    }

     DEBUG_INFO("SaveDataBase", "Transport Matrix:  PTA-terms");
    for (vsiz_t i = 0; i < fPTAMatrixElems.size(); i++) {
        fprintf(file, "PTA ");
        const THaMatrixElement& m = fPTAMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            fprintf(file, "%d ", m.pw[j]);
        }
        int j;
        for (j = 0; j < m.order; j++) {
            fprintf(file, " %13.6e", m.poly[j]);
        }
        for (; j < kPORDER; j++) {
            fprintf(file, " %13.6e", 0.0);
        }
        fprintf(file, "  %d", m.OptOrder);
        fprintf(file, "\n");
    }

    fprintf(file, DatabaseFooter);

    fclose(file);

    return kOK;
}

void LOpticsOpt::Print(const Option_t* opt) const
{
    // Print current matrix

    THaTrackingDetector::Print(opt);
    typedef vector<THaMatrixElement>::size_type vsiz_t;

    // Print out the optics matrices, to verify they make sense
    //    printf("Matrix FP (t000, y000, p000)\n");
    //    for (vsiz_t i = 0; i < fFPMatrixElems.size(); i++) {
    //        const THaMatrixElement& m = fFPMatrixElems[i];
    //        for (vsiz_t j = 0; j < m.pw.size(); j++) {
    //            printf("  %2d", m.pw[j]);
    //        }
    //        for (int j = 0; j < m.order; j++) {
    //            printf("  %g", m.poly[j]);
    //        }
    //        printf(" : Opt -> %d", m.OptOrder);
    //        printf("\n");
    //    }

    printf("LOpticsOpt::Print: Transport Matrix:  D-terms\n");
    for (vsiz_t i = 0; i < fDMatrixElems.size(); i++) {
        const THaMatrixElement& m = fDMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            printf("  %2d", m.pw[j]);
        }
        for (int j = 0; j < m.order; j++) {
            printf("\t%g", m.poly[j]);
        }
        printf(" : Opt -> %d", m.OptOrder);
        if ((UInt_t) m.order != m.OptOrder) printf(" != Matrix Order !!");
        printf("\n");
    }

    printf("LOpticsOpt::Print: Transport Matrix:  T-terms\n");
    for (vsiz_t i = 0; i < fTMatrixElems.size(); i++) {
        const THaMatrixElement& m = fTMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            printf("  %2d", m.pw[j]);
        }
        for (int j = 0; j < m.order; j++) {
            printf("\t%g", m.poly[j]);
        }
        printf(" : Opt -> %d", m.OptOrder);
        if ((UInt_t) m.order != m.OptOrder) printf(" != Matrix Order !!");
        printf("\n");
    }

    printf("LOpticsOpt::Print: Transport Matrix:  Y-terms\n");
    for (vsiz_t i = 0; i < fYMatrixElems.size(); i++) {
        const THaMatrixElement& m = fYMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            printf("  %2d", m.pw[j]);
        }
        for (int j = 0; j < m.order; j++) {
            printf("\t%g", m.poly[j]);
        }
        printf(" : Opt -> %d", m.OptOrder);
        if ((UInt_t) m.order != m.OptOrder) printf(" != Matrix Order !!");
        printf("\n");
    }

    //    printf("Transport Matrix:  YTA-terms (abs(theta))\n");
    //    for (vsiz_t i = 0; i < fYTAMatrixElems.size(); i++) {
    //        const THaMatrixElement& m = fYTAMatrixElems[i];
    //        for (vsiz_t j = 0; j < m.pw.size(); j++) {
    //            printf("  %2d", m.pw[j]);
    //        }
    //        for (int j = 0; j < m.order; j++) {
    //            printf("\t%g", m.poly[j]);
    //        }
    //        printf(" : Opt -> %d", m.OptOrder);
    //        printf("\n");
    //    }

    printf("LOpticsOpt::Print: Transport Matrix:  P-terms\n");
    for (vsiz_t i = 0; i < fPMatrixElems.size(); i++) {
        const THaMatrixElement& m = fPMatrixElems[i];
        for (vsiz_t j = 0; j < m.pw.size(); j++) {
            printf("  %2d", m.pw[j]);
        }
        for (int j = 0; j < m.order; j++) {
            printf("\t%g", m.poly[j]);
        }
        printf(" : Opt -> %d", m.OptOrder);
        if ((UInt_t) m.order != m.OptOrder) printf(" != Matrix Order !!");
        printf("\n");
    }

       printf("Transport Matrix:  PTA-terms\n");
       for (vsiz_t i = 0; i < fPTAMatrixElems.size(); i++) {
           const THaMatrixElement& m = fPTAMatrixElems[i];
           for (vsiz_t j = 0; j < m.pw.size(); j++) {
               printf("  %2d", m.pw[j]);
           }
           for (int j = 0; j < m.order; j++) {
               printf("\t%g", m.poly[j]);
           }
           printf(" : Opt -> %d", m.OptOrder);
           printf("\n");
       }
    
       printf("Matrix L\n");
       for (vsiz_t i = 0; i < fRMatrixElems.size(); i++) {
           const THaMatrixElement& m = fRMatrixElems[i];
           for (vsiz_t j = 0; j < m.pw.size(); j++) {
               printf("  %2d", m.pw[j]);
           }
           for (int j = 0; j < m.order; j++) {
               printf("\t%g", m.poly[j]);
           }
           printf(" : Opt -> %d", m.OptOrder);
           printf("\n");
       }

    printf("fArbitaryVertexShift[%d] = {", NFoils);
    for (UInt_t FoilID = 0; FoilID < NFoils; FoilID++)
        printf("%f  ", fArbitaryVertexShift[FoilID]);
    printf("}\n");
    printf("fArbitaryDpKinShift[%d] = {", NKine);
    for (UInt_t KineID = 0; KineID < NKine; KineID++)
        printf("%f  ", fArbitaryDpKinShift[KineID]);
    printf("}\n");

    return;
}

UInt_t LOpticsOpt::Matrix2Array(Double_t Array[], const std::vector<THaMatrixElement> &Matrix, Bool_t FreeParaFlag[])
{
    // Matrix -> Array

    typedef vector<THaMatrixElement>::size_type vsiz_t;

    UInt_t idx = 0;

    for (vsiz_t i = 0; i < Matrix.size(); i++) {
        const THaMatrixElement& m = Matrix[i];
        UInt_t j;
        for (j = 0; (int) j < m.order; j++) {
            if (FreeParaFlag) FreeParaFlag[idx] = j < m.OptOrder ? kTRUE : kFALSE;
            Array[idx++] = m.poly[j];
        }
        for (; j < kPORDER; j++) {
            if (FreeParaFlag) FreeParaFlag[idx] = j < m.OptOrder ? kTRUE : kFALSE;
            Array[idx++] = 0;
        }
    }

    DEBUG_INFO("Matrix2Array", "Fill Size = %d", idx);

    return idx;
}

UInt_t LOpticsOpt::Array2Matrix(const Double_t Array[], std::vector<THaMatrixElement> &Matrix)
{
    // Array -> fCurrentMatrixElems

    typedef vector<THaMatrixElement>::size_type vsiz_t;

    UInt_t idx = 0;
    Double_t maxele = 0, sumele = 0;
    vsiz_t max_i = -1;
    int max_j = -1;

    for (vsiz_t i = 0; i < Matrix.size(); i++) {
        THaMatrixElement& m = Matrix[i];
        int j;
        m.order = kPORDER;
        for (j = 0; j < m.order; j++) {
            if (TMath::Abs(Array[idx]) > maxele) {
                max_i = i;
                max_j = j;
            }
            sumele += TMath::Abs(Array[idx]);
            m.poly[j] = Array[idx];
            idx++;
        }
        m.SkimPoly();
    }

    assert(max_i < Matrix.size());
    assert(max_j >= 0);
    THaMatrixElement& m = Matrix[max_i];
    DEBUG_INFO("Array2Matrix", "Load Size = %d, max ele = (%d %d %d %d) = %f, average = %f", idx, max_j, m.pw[0], m.pw[1], m.pw[2], m.poly[max_j], sumele / idx);

    return idx;
}

///////////////////////////////////////////////////////////////////////////////
// Data storage
///////////////////////////////////////////////////////////////////////////////

UInt_t LOpticsOpt::LoadRawData(TString DataFileName, UInt_t NLoad, UInt_t MaxDataPerGroup, TString Old_DataFile_name)
{
  




  
  // ignore/ select holes based on other optimisation

  if(hole_ignore){
    
    hole_selector(Old_DataFile_name);
    
  // if(Old_DataFile_name.CompareTo("")){   
  //   hole_select_copy(Old_DataFile_name);
  // }
  // else{
  //   cout << "Old_DataFile_name = " << Old_DataFile_name << endl;
  //   hole_select_copy("");
  // }
  

  // // ignore prescibed holes (based on file hole_opt_select.C)
  // //  hole_opt_select();

  }


  
    DEBUG_INFO("LoadRawData", "Loading %s", DataFileName.Data());

    if (BeamShiftX != 0)
        DEBUG_WARNING("LoadRawData", "Shift Beam X = %f", BeamShiftX);

    UInt_t datagrpcnt[kMaxDataGroup] = {0};

    FILE* file = fopen(DataFileName, "r");
    if (!file) return kFileError;

    UInt_t NRead = 0;
    const int LEN = 2000;
    char buff[LEN];

    Double_t NDataRead = 0;
    int NLineRead = 0;

    while (fgets(buff, LEN, file)) {
        NLineRead++;

        if (NLineRead % 100000 == 0) DEBUG_INFO("LoadRawData", "%d/%d Entries Loaded", NRead, NLineRead);

        assert(NRead < MaxNRawData); //too much data if fails

        if (NRead >= NLoad) break;

        Double_t * eventdata = fRawData[NRead].Data;

        string line(buff);
        // Erase trailing newline
        if (line.size() > 0 && line[line.size() - 1] == '\n') {
            buff[line.size() - 1] = 0;
            line.erase(line.size() - 1, 1);
        }
        // Split the line into whitespace-separated fields
        vector<string> line_spl = Split(line);

        assert(line_spl.size() <= MaxNEventData); // array size check
        for (UInt_t idx = 0; idx < line_spl.size(); idx++)
            eventdata[idx] = atof(line_spl[idx].c_str());

        // WARNING : shift beam x
        if (BeamShiftX != 0)
            eventdata[kBeamX] += BeamShiftX;

        // // determine whether to save this data
        // UInt_t cutid = (UInt_t) eventdata[kCutID];
        // assert(cutid < kMaxDataGroup); // too many cuts
        // UInt_t & grpcnt = datagrpcnt[cutid];
        // grpcnt++;
        // if (grpcnt > MaxDataPerGroup) {
        //     DEBUG_MASSINFO("LoadRawData", "ignore data %d from cutid %d (%d ev total)", NLineRead, cutid, grpcnt);
        //     continue;
        // }

	UInt_t foilid = eventdata[kFoilID];
	UInt_t colid  = eventdata[kColID];
	UInt_t rowid  = eventdata[kRowID];



	if(hole_ignore){

	  if(!hole_select[foilid][colid][rowid]){
	    continue;
	  }	   
	  if(!foil_select[foilid]){
	    continue;
	  }	
	  if(!col_select[foilid][colid]){
	    continue;
	  }
	  if(!row_select[foilid][rowid]){
	    continue;
	  }
	
	}



        NDataRead += line_spl.size();

        Double_t(*powers)[5] = fRawData[NRead].powers;
        Double_t x_fp = eventdata[kX];
        Double_t th_fp = eventdata[kTh];
        Double_t y_fp = eventdata[kY];
        Double_t ph_fp = eventdata[kPhi];

	
	// Int_t count = 0;
	// if(foilid == 6 && count%1000 == 0){

	//   cout << "Event in Foil "  << foilid << " has col: " << colid << " and row: " << rowid << endl;
	// }
	
	

        // calculate the powers we need
        for (int i = 0; i < kNUM_PRECOMP_POW; i++) {
            powers[i][0] = pow(x_fp, i);
            powers[i][1] = pow(th_fp, i);
            powers[i][2] = pow(y_fp, i);
            powers[i][3] = pow(ph_fp, i);
            powers[i][4] = pow(TMath::Abs(th_fp), i);
        }

        NRead++;
    }

    fclose(file);
    fNRawData = NRead;
    fNCalibData = NRead; // fNCalibData shall be updated later if only part of data read in is for calibration use

    UInt_t goodstatcut = 0, actcutcnt = 0;
    for (int i = 0; i < kMaxDataGroup; i++) {
        if (datagrpcnt[i] > 0) actcutcnt++;

        if (datagrpcnt[i] > MaxDataPerGroup) goodstatcut++;
    }

    DEBUG_INFO("LoadRawData", "Event Limit/Cut = %d, %d / %d ev read, %d / %d cut have enough ev", MaxDataPerGroup, NRead, NLineRead, goodstatcut, actcutcnt);
    DEBUG_INFO("LoadRawData", "%d events x %f record/event read from %s", fNRawData, NDataRead / fNRawData, DataFileName.Data());

    return NRead;
}

///////////////////////////////////////////////////////////////////////////////
// Optimization related Commands
///////////////////////////////////////////////////////////////////////////////

const TVector3 LOpticsOpt::GetSieveHoleTCS(UInt_t Col, UInt_t Row)
{
    assert(Col < NSieveCol);
    assert(Row < NSieveRow);

    Double_t XbyRow =0;
    /*
    if((Col%2) == 0){
      XbyRow = SieveXbyRowEven[Row];
    }
    if((Col%2) == 1){
      XbyRow = SieveXbyRowOdd[Row];
    }
    */

    XbyRow = SieveXbyRow[Row];

    //    TVector3 Sieve_offset_HCS(SieveOffX_HCS,SieveOffY_HCS,ZPos_HCS);

    //    TVector3 Sieve_offset_TCS = fTCSInHCS.Inverse()*(Sieve_offset_HCS);
    
    //    TVector3 SieveHoleTCS( Sieve_offset_TCS.X() + XbyRow, Sieve_offset_TCS.Y() + SieveYbyCol[Col], Sieve_offset_TCS.Z());

    
    TVector3 SieveHoleTCS(XbyRow, SieveYbyCol[Col], 0);

    
    SieveHoleTCS.RotateX(-(yaw - HRSAngle));
    SieveHoleTCS.RotateY(pitch - TMath::Pi()/2);
    //    SieveHoleTCS.SetZ(SieveHoleTCS.Z();
    /*
    cout<<"Col%2:"<<Col%2<<endl;
    cout<<"Col:"<<Col<<endl;
    cout<<"XbyRow:"<<XbyRow<<endl;
    cout<<"Row:"<<Row<<endl;
    cout<<"YbyCol:"<<SieveYbyCol[Col]<<endl;
    cout<<"*******************"<<endl;
    */


    SieveHoleTCS.SetXYZ( SieveHoleTCS.X() + SieveOffX,  SieveHoleTCS.Y() + SieveOffY, SieveHoleTCS.Z() + ZPos + SieveOffZ);
		      

		      
    return SieveHoleTCS;
}

const TVector3 LOpticsOpt::GetSieveHoleCorrectionTCS(UInt_t nfoil, UInt_t Col, UInt_t Row)
{
    assert(nfoil <NFoils);
    assert(Col < NSieveCol);
    assert(Row < NSieveRow);
     
    //consider the difference of real distribution and hole center
    Double_t Y_p=0, Y_m=0, Yback_p=0, Yback_m=0; // Y the sieve hole Y position
    Double_t X_p=0, X_m=0, Xback_p=0, Xback_m=0; // X the sieve hole X position
    Double_t Yreal_p=0, Yreal_m=0, Xreal_p=0, Xreal_m=0; // real distrituion limit
    Double_t Z_distance=0;
    Double_t BeamYHCS=0,BeamXHCS=0;
    //    Double_t SieveY_Correction[NFoils][NSieveCol][NSieveRow] ={{{0}}};
    //    Double_t SieveX_Correction[NFoils][NSieveCol][NSieveRow] ={{{0}}};

    //const TVector3 BeamSpotH CS_average(BeamX_average, BeamY_average, targetfoils[nfoil]); 
   //    const TVector3 BeamSpotHCS_average(BeamX_average[nfoil], BeamY_average, targetfoils[nfoil]);


    //    const TVector3 BeamSpotHCS_average(BeamX_average, BeamY_average, targetfoils[nfoil]);
    const TVector3 BeamSpotHCS_average(BeamX_average[nfoil] + (targetfoils[nfoil]/BeamZDir_average)*BeamXDir_average[nfoil], BeamY_average[nfoil] + (targetfoils[nfoil]/BeamZDir_average)*BeamYDir_average[nfoil], targetfoils[nfoil]);




    const TVector3 BeamSpotTCS_average = fTCSInHCS.Inverse()*(BeamSpotHCS_average - fPointingOffset);        
    BeamXHCS= BeamSpotTCS_average.X();
    BeamYHCS= BeamSpotTCS_average.Y();

    const TVector3 SieveHoleTCS = GetSieveHoleTCS(Col, Row);

    Z_distance = SieveHoleTCS.Z() - BeamSpotTCS_average.Z();
    Y_p = SieveHoleTCS.Y() +  SieveRadius/2.* 25.4e-3; //  // .157/2. * 25.4e-3 is the sieve hole radius
    Y_m = SieveHoleTCS.Y() - SieveRadius/2. * 25.4e-3;
    Yback_p = Z_distance /(Z_distance + 25.4e-3) * (Y_p-BeamYHCS) + BeamYHCS;
    Yback_m = Z_distance /(Z_distance + 25.4e-3) * (Y_m-BeamYHCS) + BeamYHCS;
    Yreal_p = (Y_p >= Yback_p) ? Yback_p : Y_p;
    Yreal_m = (Y_m >= Yback_m) ? Y_m : Yback_m;
    X_p = SieveHoleTCS.X() +  SieveRadius/2. * 25.4e-3; 
    X_m = SieveHoleTCS.X() -  SieveRadius/2. * 25.4e-3;
    Xback_p = Z_distance /(Z_distance + 25.4e-3) * (X_p-BeamXHCS) + BeamXHCS;
    Xback_m = Z_distance /(Z_distance + 25.4e-3) * (X_m-BeamXHCS) + BeamXHCS;
    Xreal_p = (X_p >= Xback_p) ? Xback_p : X_p;
    Xreal_m = (X_m >= Xback_m) ? X_m : Xback_m;

    const TVector3 SieveHoleCorrectionTCS((Xreal_p + Xreal_m)/2, (Yreal_p + Yreal_m)/2 ,SieveHoleTCS.Z());
    //const TVector3 SieveHoleCorrectionTCS( SieveHoleTCS.X(), SieveHoleTCS.Y() ,SieveHoleTCS.Z() );
    return SieveHoleCorrectionTCS;

   
}


void LOpticsOpt::SieveCheck(Int_t FoilID){

TRotation fTCSInHCS;

   TVector3 TCSX(0,-1,0);
  TVector3 TCSZ(TMath::Sin(HRSAngle),0,TMath::Cos(HRSAngle));
  TVector3 TCSY = TCSZ.Cross(TCSX);
  fTCSInHCS.RotateAxes(TCSX,TCSY,TCSZ);

  //  TVector3 fPointingOffset;
  fPointingOffset.SetXYZ(-MissPointZ*TMath::Sin(HRSAngle)*TMath::Cos(HRSAngle),(Double_t)MissPointY,MissPointZ*TMath::Sin(HRSAngle)*TMath::Sin(HRSAngle));


  cout << "fPointingOffset [= " << fPointingOffset.X() << "," << fPointingOffset.Y() << "," << fPointingOffset.Z() << "]" << endl;

  //TVector3 BeamSpotHCS_average(BeamX_average,BeamY_average,targetfoils[FoilID]);
  TVector3 BeamSpotHCS_average(BeamX_average[FoilID], BeamY_average[FoilID], targetfoils[FoilID]);
  TVector3 BeamSpotTCS_average = fTCSInHCS.Inverse()*(BeamSpotHCS_average-fPointingOffset);


  // vector gather all row, column values so they can be sorted later

  // corrected values
  vector<Double_t> Col_cor_vals;
  vector<Double_t> Row_cor_vals;

  Double_t Hole_cor_x;
  Double_t Hole_cor_y;

  
  // uncorrected values
  vector<Double_t> Col_vals;
  vector<Double_t> Row_vals;

  Double_t Hole_x;
  Double_t Hole_y;


  cout << "targetfoils[" << FoilID << "] = " <<  targetfoils[FoilID] << endl;
  
  for( Int_t row = 0; row<NSieveRow; row++){

    
    for( Int_t col = 0; col<NSieveCol; col++){


      //      Int_t Hole = Get_Hole(col, row);
      

      // ROpticsOpt *opt = new ROpticsOpt();

      //      TVector3 Hole_pos = opt->GetSieveHoleTCS(col,row);


      TVector3 Hole_pos_nocorrect = GetSieveHoleTCS(col,row);

      TVector3 Hole_pos = GetSieveHoleCorrectionTCS(FoilID,col,row);
      
      TVector3 MomDirectionTCS_hole = Hole_pos - BeamSpotTCS_average;
      
      Double_t theta_hole = MomDirectionTCS_hole.X()/MomDirectionTCS_hole.Z();
      Double_t phi_hole = MomDirectionTCS_hole.Y()/MomDirectionTCS_hole.Z();
      
      cout << "For row " << row << ", column " << col << " theta = " << theta_hole << ", phi = " << phi_hole << ", Hole_pos_nocorrect.Y() = " << Hole_pos_nocorrect.Y() << ", Hole_pos.Y() = " << Hole_pos.Y() << endl;

      //<< ", BeamSpotTCS_average = " << BeamSpotTCS_average.Y() << endl;

      Hole_cor_x = Hole_pos.X();
      Hole_cor_y = Hole_pos.Y();

      Hole_x = Hole_pos_nocorrect.X();
      Hole_y = Hole_pos_nocorrect.Y();


      if(row == (NSieveRow - 1)){
	Col_vals.push_back(Hole_y);
	Col_cor_vals.push_back(Hole_cor_y);
      }
      
    }

    Row_vals.push_back(Hole_x);
    Row_cor_vals.push_back(Hole_cor_x);
    
  }

  cout << "wefg" << endl;

  // print rows and columns in ascending ord
  sort(Row_vals.begin(),Row_vals.end());
  sort(Col_vals.begin(),Col_vals.end());

  sort(Row_cor_vals.begin(),Row_cor_vals.end());
  sort(Col_cor_vals.begin(),Col_cor_vals.end());


  cout << "Rows (x_values):" << endl;
  for(auto const& value: Row_vals){
    cout << value << " "; 
  }
  cout << endl;

  cout << "Rows corrected (x_values):" << endl;
  for(auto const& value: Row_cor_vals){
    cout << value << " "; 
  }
  cout << endl << endl;

  
  cout << "Cols (y_values):" << endl;
  for(auto const& value: Col_vals){
    cout << value << " "; 
  }
  cout << endl;

  cout << "Cols corrected (y_values):" << endl;
  for(auto const& value: Col_cor_vals){
    cout << value << " "; 
  }
  cout << endl;
  
}


void LOpticsOpt::PrepareSieve(void)
{
    // Calculate kRealTh, kRealPhi

    // DEBUG_INFO("PrepareSieve","Entry Point");

     Double_t exttargcorr_th = 0, rms_exttargcorr_th = 0;
 
     for (UInt_t idx = 0; idx < fNRawData; idx++) {
        EventData &eventdata = fRawData[idx];

        // UInt_t res = (UInt_t) eventdata.Data[kCutID];
        // res = res % (NSieveRow * NSieveCol * NFoils);
        // const UInt_t FoilID = res / (NSieveRow * NSieveCol); //starting 0!
        // res = res % (NSieveRow * NSieveCol);
        // const UInt_t Col = res / (NSieveRow); //starting 0!
        // const UInt_t Row = res % (NSieveRow); //starting 0!

	const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID]; //starting 0!
	// res = res%(NSieveRow*NSieveCol);
	const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
	const UInt_t Row = (UInt_t)eventdata.Data[kRowID]; //starting 0!

	
	
        assert(FoilID < NFoils); //check array index size
	const TVector3 SieveHoleCorrectionTCS = GetSieveHoleCorrectionTCS(FoilID, Col, Row);
        eventdata.Data[kSieveX] = SieveHoleCorrectionTCS.X();
        eventdata.Data[kSieveY] = SieveHoleCorrectionTCS.Y();
        eventdata.Data[kSieveZ] = SieveHoleCorrectionTCS.Z();

	//	const TVector3 BeamSpotHCS(eventdata.Data[kBeamX], eventdata.Data[kBeamY], targetfoils[FoilID]);
	const TVector3 BeamSpotHCS( eventdata.Data[kBeamX] + (targetfoils[FoilID]/BeamZDir_average)*eventdata.Data[kBeamDirX], eventdata.Data[kBeamY] + (targetfoils[FoilID]/BeamZDir_average)*eventdata.Data[kBeamDirY], targetfoils[FoilID]);



	//const TVector3 BeamSpotHCS(BeamX_average[FoilID], BeamY_average, targetfoils[FoilID]);
        eventdata.Data[kBeamZ] = targetfoils[FoilID];

        const TVector3 BeamSpotTCS = fTCSInHCS.Inverse()*(BeamSpotHCS - fPointingOffset);        

	const TVector3 MomDirectionTCS = SieveHoleCorrectionTCS - BeamSpotTCS;

        eventdata.Data[kRealTh] = MomDirectionTCS.X() / MomDirectionTCS.Z();
        eventdata.Data[kRealPhi] = MomDirectionTCS.Y() / MomDirectionTCS.Z();
 
	
	
        const Double_t x_tg = BeamSpotTCS.X() - BeamSpotTCS.Z() * eventdata.Data[kRealTh];
        const Double_t y_tg = BeamSpotTCS.Y() - BeamSpotTCS.Z() * eventdata.Data[kRealPhi];
	
        eventdata.Data[kRealTgX] = x_tg;
        eventdata.Data[kRealTgY] = y_tg;

        // Expected th ph before ext. target correction
        // fDeltaTh = fThetaCorr * x_tg;
        // Double_t theta = trkifo->GetTheta() + fDeltaTh;
        eventdata.Data[kRealThMatrix] = eventdata.Data[kRealTh] - x_tg * ExtTarCor_ThetaCorr;

        exttargcorr_th += x_tg * ExtTarCor_ThetaCorr;
        rms_exttargcorr_th += x_tg * ExtTarCor_ThetaCorr * x_tg * ExtTarCor_ThetaCorr;
	
        DEBUG_MASSINFO("PrepareSieve", "Real_Th_Matrix = %f,\t Real_Phi = %f", eventdata.Data[kRealThMatrix], eventdata.Data[kRealPhi]);
    }

    DEBUG_INFO("PrepareSieve", "Average Extended Target Correction: th = %f,\t rms_th = %f", exttargcorr_th / fNRawData, TMath::Sqrt(rms_exttargcorr_th / fNRawData));

    // Make sure kCalcTh, kCalcPh is filled
    SumSquareDTh();
    SumSquareDPhi();

    DEBUG_INFO("PrepareSieve", "Done!");
}


// function designed to ignore tail of distribs
// can be used after
// preparesieve
void LOpticsOpt::Ignore_tail(){
  
  
  // set events to ingore from other foils


  // check rate of pass/ fail of ignore cut
  Int_t Npass = 0;
  Int_t Nfail = 0;

  Double_t pass_rate = 0;

  cout << "gaus_start[0] = " << gaus_start[0] << ", gaus_end[0] = " << gaus_end[0] << endl;
  cout << "gaus_start[1] = " << gaus_start[1] << ", gaus_end[1] = " << gaus_end[1] << endl;
  cout << "gaus_start[2] = " << gaus_start[2] << ", gaus_end[2] = " << gaus_end[2] << endl;
    
  for (UInt_t idx = 0; idx < fNRawData; idx++) {

      EventData &eventdata = fRawData[idx];
      

      UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID]; //starting 0!


      const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
      const UInt_t Row = (UInt_t)eventdata.Data[kRowID]; //starting 0!
      
      const TVector3 SieveHoleCorrectionTCS = GetSieveHoleCorrectionTCS(FoilID, Col, Row);

      Double_t ProjectionY = eventdata.Data[kRealTgY] + eventdata.Data[kCalcPh] * (SieveHoleCorrectionTCS.Z());
      
      Double_t dY = 1000*(ProjectionY - SieveHoleCorrectionTCS.Y());
      
      UInt_t PlotID = FoilID;

      if(NFoils<11){

	PlotID -= 8;
      }
      
      if(dY < gaus_end[PlotID] && dY > gaus_start[PlotID]){

	
	eventdata.Data[kIgnore] = 1;
	Npass++;

      }
      else{

	eventdata.Data[kIgnore] = 0;
	Nfail++;
	
      }



    }

    pass_rate = Double_t(Npass/(Npass+Nfail));

    cout << "pass_rate = " << pass_rate << endl;
    cout << "Npass = " << Npass << " and Nfail = " << Nfail << endl;

  

}

// function designed to print holes numbers used in optimisation to file

void LOpticsOpt::Print_holes(TString DataSource){


  Bool_t hole_used[NFoils][NSieveCol][NSieveRow] = {kFALSE};

  for (UInt_t idx = 0; idx < fNRawData; idx++) {
    EventData &eventdata = fRawData[idx];


    const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID]; //starting 0!
    const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
    const UInt_t Row = (UInt_t)eventdata.Data[kRowID]; //starting 0!
    
    Int_t PlotID = FoilID;

    if(NFoils<11){
      PlotID -= 8;
    }

    hole_used[PlotID][Col][Row] = kTRUE;
    
  }

  // create output file to print results

  ofstream ofile;
  //  ofile.open(Form("plots/holes_used/"));

  DataSource(0,24) = "";
  TString file_name(TString("plots/holes_used/") + DataSource);
  cout << "File_name = " << file_name << endl;

  ofile.open(TString("plots/holes_used/") + DataSource);
  ofile << DataSource << endl;
  ofile << "Foil Column  Row" << endl;

  for(Int_t iFoil = 0; iFoil<NFoils; iFoil++){

    
    // Int_t PlotID = iFoil;

    // if(NFoils<11){
    //   PlotID -= 8;
    // }

    for (Int_t iCol = 0; iCol<NSieveCol; iCol++){
      for (Int_t iRow = 0; iRow<NSieveRow; iRow++){

	if(hole_used[iFoil][iCol][iRow]){
	  cout << "Foil = " << iFoil << ", Col = " << iCol << ", Row = " << iRow << endl;
	  ofile << iFoil << " " <<  iCol << " " << iRow << endl;
	}
	
      }
    }       
  }

  ofile.close();


}

TCanvas * LOpticsOpt::CheckSieve(Int_t PlotFoilID)
{
    // Visualize Sieve Plane
    DEBUG_INFO("CheckSieve", "Entry Point");
    
    //const UInt_t nplot = (PlotFoilID == -1) ? NFoils : 9;
    const UInt_t nplot = PlotFoilID ;


    // theta and phi diffs
    Double_t dTh = 0, dPh = 0;

    // sieve x and y diffs
    Double_t dX = 0, dY = 0;

    enum {
      kEventID, kRealSieveX, kRealSieveY, kCalcSieveX, kCalcSieveY, kRealSieveTh, kRealSievePh, KCalcSieveTh, kCalcSievePh
    };
    
    Double_t SieveEventID[NFoils][NSieveCol][NSieveRow][5] = {
      {
            {
	      {0}
            }
      }
    };




    TH2D * HSieveAngle[NFoils] = {0};
    Double_t th_lim[NFoils] = {0};
    Double_t ph_lim[NFoils] = {0};

    
    TH2D * HSievePlane[NFoils] = {0};
    Double_t x_lim[NFoils] = {0};
    Double_t y_lim[NFoils] = {0};

    for (UInt_t idx = 0; idx < NFoils; idx++) {
      x_lim[idx] = 1.2 * TMath::Max(TMath::Abs(SieveYbyCol[0]), TMath::Abs(SieveYbyCol[NSieveCol - 1]));
        //y_lim[idx] = 1.5 * TMath::Max(TMath::Abs(SieveXbyRowOdd[0]), TMath::Abs(SieveXbyRowEven[NSieveRow - 1]));
      y_lim[idx] = 1.2 * TMath::Max(TMath::Abs(SieveXbyRow[0]), TMath::Abs(SieveXbyRow[NSieveRow - 1]));
      
      HSievePlane[idx] = new TH2D(Form("Sieve_Foil%d", idx), "Sieve Plane Proj. (tg_X vs tg_Y) for " + Foil_names[idx], 500, -x_lim[idx], x_lim[idx], 500, -y_lim[idx], y_lim[idx]);

      HSievePlane[idx]->SetXTitle("Sieve H [m]");
      HSievePlane[idx]->SetYTitle("Sieve V [m]");
      assert(HSievePlane[idx]); // assure memory allocation


      //      HSieveAngle[idx] = new TH2D(Form("Angle_Sieve_Foil%d",idx),Form("theta_target vs. phi_target for Data set #%d",idx), 500, -0.04, 0.04, 500, -0.08, 0.08);
      HSieveAngle[idx] = new TH2D(Form("Angle_Sieve_Foil%d",idx),Form("theta_target vs. phi_target for Data set #%d",idx), 500, -0.04, 0.04, 500, -0.08, 0.08);
      

      HSieveAngle[idx]->SetXTitle("Phi");
      HSieveAngle[idx]->SetYTitle("Theta");
      assert(HSieveAngle[idx]); // assure memory allocation
    }
    
    Int_t PlotID; // useful for case when only Vertical foils are used
                    // (FoilID based on scheme with optics and vertical foils where optics foils
                    // are 0-7 and vertical foils 8-10) so they can have PlotID 0,1,2

    
    for (UInt_t idx = 0; idx < fNRawData; idx++) {
      const EventData &eventdata = fRawData[idx];
      
      // UInt_t res = (UInt_t) eventdata.Data[kCutID];
      // res = res % (NSieveRow * NSieveCol * NFoils);
      // const UInt_t FoilID = res / (NSieveRow * NSieveCol); //starting 0!
      // res = res % (NSieveRow * NSieveCol);
      // const UInt_t Col = res / (NSieveRow); //starting 0!
      // const UInt_t Row = res % (NSieveRow); //starting 0!
      
      const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID]; //starting 0!
      // res = res%(NSieveRow*NSieveCol);
      const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
      const UInt_t Row = (UInt_t)eventdata.Data[kRowID]; //starting 0!
      
      
      assert(FoilID < NFoils); //array index check
      
      const TVector3 SieveHoleCorrectionTCS = GetSieveHoleCorrectionTCS(FoilID, Col, Row);


      PlotID = FoilID; // useful for case when only Vertical foils are used
                    // (FoilID based on scheme with optics and vertical foils where optics foils
                    // are 0-7 and vertical foils 8-10) so they can have PlotID 0,1,2
      if(NFoils<11){

	PlotID -= 8;
      }
      
      // calculate sieve_x & sieve_y
      
      Double_t ProjectionX = eventdata.Data[kRealTgX] + (eventdata.Data[kCalcTh] + eventdata.Data[kRealTgX] * ExtTarCor_ThetaCorr) * (SieveHoleCorrectionTCS.Z());
      
      
      //Double_t ProjectionX = eventdata.Data[kCalcTh];
      Double_t ProjectionY = eventdata.Data[kRealTgY] + eventdata.Data[kCalcPh] * (SieveHoleCorrectionTCS.Z());
      
	
      //Double_t ProjectionY = eventdata.Data[kCalcPh];
	

      // calculate theta and phi 'real' positions of holes

      //      TVector3 BeamSpotHCS_average(BeamX_average,BeamY_average,targetfoils[FoilID]);
      const TVector3 BeamSpotHCS_average(BeamX_average[FoilID] + (targetfoils[FoilID]/BeamZDir_average)*BeamXDir_average[FoilID], BeamY_average[FoilID] + (targetfoils[FoilID]/BeamZDir_average)*BeamYDir_average[FoilID], targetfoils[FoilID]);
    
      
      TVector3 BeamSpotTCS_average = fTCSInHCS.Inverse()*(BeamSpotHCS_average-fPointingOffset);

      TVector3 MomDirectionTCS_hole = SieveHoleCorrectionTCS - BeamSpotTCS_average;
      
      Double_t theta_hole = MomDirectionTCS_hole.X()/MomDirectionTCS_hole.Z();
      Double_t phi_hole = MomDirectionTCS_hole.Y()/MomDirectionTCS_hole.Z();
      


      
      // fill Sieve x,y and theta,phi plots


      HSievePlane[FoilID]->Fill(ProjectionY, ProjectionX);
      HSieveAngle[FoilID]->Fill(eventdata.Data[kCalcPh],eventdata.Data[kCalcTh]);


      /*
	if(FoilID==1 && Row==6) {
	cout<<"*********"<<endl;
	cout<<"Row:"<<Row<<endl;
	cout<<"XRow:"<<GetSieveHoleTCS(Col,Row).X()<<endl;
	cout<<"ProjectionX:"<<ProjectionX<<endl;
	cout<<"dx:"<<ProjectionX - SieveHoleCorrectionTCS.X()<<endl;
	HSievePlane[FoilID]->Fill(ProjectionY, ProjectionX);
	}
	*/
      dX += ProjectionX - SieveHoleCorrectionTCS.X();
      dY += ProjectionY - SieveHoleCorrectionTCS.Y();
	//	cout<<"dX:"<<dX<<endl;

	
	

        SieveEventID[PlotID][Col][Row][kEventID] = idx;
        SieveEventID[PlotID][Col][Row][kRealSieveX] = SieveHoleCorrectionTCS.X();
        SieveEventID[PlotID][Col][Row][kRealSieveY] = SieveHoleCorrectionTCS.Y();
        SieveEventID[PlotID][Col][Row][kCalcSieveX] = ProjectionX;
        SieveEventID[PlotID][Col][Row][kCalcSieveY] = ProjectionY;
    }



    
    DEBUG_INFO("CheckSieve", "Average : D_X = %f,\t D_Y = %f", dX / fNRawData, dY / fNRawData);

    //    TCanvas * c1 = new TCanvas("SieveCheck", "SieveCheck", 1800, 1100);
    TCanvas * c1;
    
    if (nplot <= 1) {
        c1 = new TCanvas("SieveCheck", "SieveCheck", 1000, 1000);
        c1->Divide(1, 1);
    } else if (nplot <= 3) {
        c1 = new TCanvas("SieveCheck", "SieveCheck", 1800, 1100);
        c1->Divide(3, 1);
    } else if (nplot <= 6) {
        c1 = new TCanvas("SieveCheck", "SieveCheck", 1800, 1100);
        c1->Divide(3, 2);
    } else if (nplot <= 9){
        c1 = new TCanvas("SieveCheck", "SieveCheck", 1800, 1100);
        c1->Divide(3, 3);
    }
    else{
        c1 = new TCanvas("SieveCheck", "SieveCheck", 1800, 1100);
	c1->Divide(4, 3);
    }
    

   
    for (UInt_t idx = 0; idx < nplot; idx++) {

      const EventData &eventdata = fRawData[idx];
      UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID]; //starting 0!
      // UInt_t FoilID = idx;
        // if (PlotFoilID >= 0)
	//   FoilID = PlotFoilID;
	//	idx = FoilID;  //Added this line, does not work with >1 foil


        c1->cd(idx + 1);
        assert(HSievePlane[idx]); //pointer check

        HSievePlane[idx]->Draw("COLZ");

        // Draw Sieve
	const Double_t plotwidth = 0.004;
        for (UInt_t Row = 0; Row < NSieveRow; Row++) {
	  for (UInt_t Col = 0; Col < NSieveCol; Col++) {
	    
	    const TVector3 SieveHoleCorrectionTCS = GetSieveHoleCorrectionTCS(FoilID, Col, Row);
	    
	    const Double_t posx = SieveHoleCorrectionTCS.Y();
	    const Double_t posy = SieveHoleCorrectionTCS.X();
	    
	    TLine *lh = new TLine(posx - plotwidth, posy, posx + plotwidth, posy);
	    TLine *lv = new TLine(posx, posy - plotwidth, posx, posy + plotwidth);


	    
	    if(Col==(NSieveCol-1)/2)   
	      { 
		lh -> SetLineColor(kRed);
		lv -> SetLineColor(kRed);
	      }
	    else{
	      lh -> SetLineColor(kBlack);
	      lv -> SetLineColor(kBlack);
	    }
	    lh -> Draw();
	    lv -> Draw();
	  }
        }

        // Draw arrows
        for (UInt_t Col = 0; Col < NSieveCol; Col++) {
	  for (UInt_t Row = 0; Row < NSieveRow; Row++) {
	    if (SieveEventID[idx][Col][Row][kEventID] > 0) {
	      assert(SieveEventID[idx][Col][Row][kEventID] < fNRawData); //array index bondary check
	      TArrow * ar2 = new TArrow(SieveEventID[idx][Col][Row][kCalcSieveY], SieveEventID[idx][Col][Row][kCalcSieveX], SieveEventID[idx][Col][Row][kRealSieveY], SieveEventID[idx][Col][Row][kRealSieveX], 0.008, "|>");
	      ar2->SetAngle(40);
	      ar2->SetLineColor(kMagenta);
	      ar2->SetFillColor(kMagenta);
	    }
	  }
        }
    }


    // plot equivalent sieve plot for theta-phi distribution

    //    TCanvas * c2 = new TCanvas("AngleCheck", "AngleCheck", 1800, 1100);
    TCanvas * c2;
    
    if (nplot <= 1) {
        c2 = new TCanvas("AngleCheck", "AngleCheck", 1000, 1000);
        c2->Divide(1, 1);
    } else if (nplot <= 3) {
        c2 = new TCanvas("AngleCheck", "AngleCheck", 1800, 1100);
        c2->Divide(3, 1);
    } else if (nplot <= 6) {
        c2 = new TCanvas("AngleCheck", "AngleCheck", 1800, 1100);
        c2->Divide(3, 2);
    } else if (nplot <= 9){
        c2 = new TCanvas("AngleCheck", "AngleCheck", 1800, 1100);
        c2->Divide(3, 3);
    }
    else{
        c2 = new TCanvas("AngleCheck", "AngleCheck", 1800, 1100);
	c2->Divide(4, 3);
    }




   
    for (UInt_t idx = 0; idx < nplot; idx++) {

      PlotID = idx;
      
      if (NFoils <11){
	PlotID += 8;
      }
      
      // if (PlotFoilID >= 0)
      // 	FoilID = PlotFoilID;
      //	idx = FoilID;  //Added this line, does not work with >1 foil
      c2->cd(idx + 1);
      assert(HSieveAngle[idx]); //pointer check
      
      HSieveAngle[idx]->Draw("COLZ");
      
      // Draw Sieve
      const Double_t plotwidth = 0.008;
      for (UInt_t Row = 0; Row < NSieveRow; Row++) {
	for (UInt_t Col = 0; Col < NSieveCol; Col++) {
	  

	  const TVector3 SieveHoleCorrectionTCS = GetSieveHoleCorrectionTCS(PlotID, Col, Row);
	  
	  //	  TVector3 BeamSpotHCS_average(BeamX_average,BeamY_average,targetfoils[PlotID]);
	  const TVector3 BeamSpotHCS_average(BeamX_average[PlotID] + (targetfoils[PlotID]/BeamZDir_average)*BeamXDir_average[PlotID], BeamY_average[PlotID] + (targetfoils[PlotID]/BeamZDir_average)*BeamYDir_average[PlotID], targetfoils[PlotID]);
	      
	  
	  TVector3 BeamSpotTCS_average = fTCSInHCS.Inverse()*(BeamSpotHCS_average-fPointingOffset);

	  TVector3 MomDirectionTCS_hole = SieveHoleCorrectionTCS - BeamSpotTCS_average;
	  
	  Double_t theta_hole = MomDirectionTCS_hole.X()/MomDirectionTCS_hole.Z();
	  Double_t phi_hole = MomDirectionTCS_hole.Y()/MomDirectionTCS_hole.Z();
	  
	  
	  TLine *lh = new TLine(phi_hole - plotwidth, theta_hole, phi_hole + plotwidth, theta_hole);
	  TLine *lv = new TLine(phi_hole, theta_hole - plotwidth, phi_hole, theta_hole + plotwidth);
	    
	    
	  
	  if(Col==(NSieveCol-1)/2)   
	    { 
	      lh -> SetLineColor(kRed);
	      lv -> SetLineColor(kRed);
	    }
	  else{
	    lh -> SetLineColor(kBlack);
	    lv -> SetLineColor(kBlack);
	  }
	  lh -> Draw();
	  lv -> Draw();
	}
      }
      
      // // Draw arrows
      // for (UInt_t Col = 0; Col < NSieveCol; Col++) {
      // 	for (UInt_t Row = 0; Row < NSieveRow; Row++) {
      // 	  if (SieveEventID[FoilID][Col][Row][kEventID] > 0) {
      // 	    assert(SieveEventID[FoilID][Col][Row][kEventID] < fNRawData); //array index bondary check
      // 	    TArrow * ar2 = new TArrow(SieveEventID[FoilID][Col][Row][kCalcSieveY], SieveEventID[FoilID][Col][Row][kCalcSieveX], SieveEventID[FoilID][Col][Row][kRealSieveY], SieveEventID[FoilID][Col][Row][kRealSieveX], 0.008, "|>");
      // 	    ar2->SetAngle(40);
      // 	      ar2->SetLineColor(kMagenta);
      // 	      ar2->SetFillColor(kMagenta);
      // 	    }
      // 	  }
      //   }

    }

    //    c1->cd(12);

    c2->SetEditable(false);
    
    return c1;
}



TCanvas * LOpticsOpt::CheckSieveAccu(Int_t PlotFoilID)
{
 
      gStyle->SetTitleW(0.3); 
      gStyle->SetTitleH(0.18); 
      gStyle->SetTitleX(0.2); 

      DEBUG_INFO("CheckSieveAccu", "Entry Point");

      const UInt_t nplot = (PlotFoilID == -1) ? NFoils : 9;

      Double_t x_lim[NFoils] = {0};
      Double_t y_lim[NFoils] = {0};

      TH1D * HSievePX[NFoils][NSieveCol] = {{0}};
      TH1D * HSievePY[NFoils][NSieveRow] = {{0}};


    for (UInt_t idx = 0; idx < NFoils; idx++) {
        x_lim[idx] = 1.3 * TMath::Max(TMath::Abs(SieveYbyCol[0]), TMath::Abs(SieveYbyCol[NSieveCol - 1]));
        //y_lim[idx] = 1.5 * TMath::Max(TMath::Abs(SieveXbyRowOdd[0]), TMath::Abs(SieveXbyRowEven[NSieveRow - 1]));
	y_lim[idx] = 1.5 * TMath::Max(TMath::Abs(SieveXbyRow[0]), TMath::Abs(SieveXbyRow[NSieveRow - 1]));
    }
    
     for (UInt_t idx = 0; idx < NFoils; idx++) {
      for (UInt_t idy = 0; idy < NSieveCol; idy++) {
	HSievePX[idx][idy] = new TH1D(Form("X_Foil%d_Col%d", idx, idy), Form("X_Foil%d_Col%d", idx, idy), 400, -y_lim[0],y_lim[0] );
	assert(HSievePX[idx][idy]); // assure memory allocation
      }
    }

    for (UInt_t idx = 0; idx < NFoils; idx++) {
      for (UInt_t idy = 0; idy < NSieveRow; idy++) {
	HSievePY[idx][idy] = new TH1D(Form("Y_Foil%d_Row%d", idx, idy), Form("Y_Foil%d_Row%d", idx, idy), 400, -x_lim[0],x_lim[0] );
	assert(HSievePY[idx][idy]); // assure memory allocation
      }
    }

    Double_t dX = 0, dY = 0;

    enum {
        kEventID, kRealSieveX, kRealSieveY, kCalcSieveX, kCalcSieveY
    };

    Double_t SieveEventID[NFoils][NSieveCol][NSieveRow][5] = {
        {
            {
                {0}
            }
        }
    };

    //save the information for each hole, including number of events and fitting information
    enum {
      kNumHole, kYmean, kYmeanErr, kYamp, kYampErr, kXmean, kXmeanErr, kXamp, kXampErr, kYndf,kYchi2,kXndf,kXchi2,kYdiff,kXdiff
    };
    Double_t SieveInfo[NFoils][NSieveCol][NSieveRow][15] = {
        {
            {
                {0}
            }
        }
    };


    cout<<"fNRawData:"<< fNRawData<<endl;

    for (UInt_t idx = 0; idx < fNRawData; idx++) {
        const EventData &eventdata = fRawData[idx];

        // UInt_t res = (UInt_t) eventdata.Data[kCutID];
        // res = res % (NSieveRow * NSieveCol * NFoils);
        // const UInt_t FoilID = res / (NSieveRow * NSieveCol); //starting 0!
        // res = res % (NSieveRow * NSieveCol);
        // const UInt_t Col = res / (NSieveRow); //starting 0!
        // const UInt_t Row = res % (NSieveRow); //starting 0!


	const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID]; //starting 0!
	// res = res%(NSieveRow*NSieveCol);
	const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
	const UInt_t Row = (UInt_t)eventdata.Data[kRowID]; //starting 0!
	
        assert(FoilID < NFoils); //array index check

        const TVector3 SieveHoleCorrectionTCS = GetSieveHoleCorrectionTCS(FoilID, Col, Row);

	Double_t ProjectionX = eventdata.Data[kRealTgX] + (eventdata.Data[kCalcTh] + eventdata.Data[kRealTgX] * ExtTarCor_ThetaCorr) * (SieveHoleCorrectionTCS.Z()); 
        Double_t ProjectionY = eventdata.Data[kRealTgY] + eventdata.Data[kCalcPh] * (SieveHoleCorrectionTCS.Z());

	HSievePX[FoilID][Col]->Fill(ProjectionX);
	HSievePY[FoilID][Row]->Fill(ProjectionY);
		
        dX += ProjectionX - SieveHoleCorrectionTCS.X();
        dY += ProjectionY - SieveHoleCorrectionTCS.Y();

        SieveEventID[FoilID][Col][Row][kEventID] = idx;
        SieveEventID[FoilID][Col][Row][kRealSieveX] = SieveHoleCorrectionTCS.X();
        SieveEventID[FoilID][Col][Row][kRealSieveY] = SieveHoleCorrectionTCS.Y();
        SieveEventID[FoilID][Col][Row][kCalcSieveX] = ProjectionX;
        SieveEventID[FoilID][Col][Row][kCalcSieveY] = ProjectionY;

	SieveInfo[FoilID][Col][Row][0] +=1;

    }

    DEBUG_INFO("CheckSieveAccu", "Average : D_X = %f,\t D_Y = %f", dX / fNRawData, dY / fNRawData);

    //define the gaussian fitting for each hole
    Double_t gXmin[NSieveCol][NSieveRow] ={{0}}; // gXmin, gXman is the Vertical fitting range for each hole
    Double_t gXmax[NSieveCol][NSieveRow] ={{0}};
    Double_t gYmin[NSieveCol][NSieveRow] ={{0}}; // gYmin, gYman is the horizontal fitting range for each hole
    Double_t gYmax[NSieveCol][NSieveRow] ={{0}};

    TF1 * gX[NFoils][NSieveCol][NSieveRow] = {{{0}}}; // define the gaussian fitting function for each hole
    TF1 * gY[NFoils][NSieveCol][NSieveRow] = {{{0}}};


    Double_t xrange,yrange,xoffset,yoffset;
    xrange = 3.3;
    yrange = 1.4;
    yoffset = .0e-3;
    for(UInt_t nfoil=0;nfoil<NFoils;nfoil++){
      for(UInt_t idx=0;idx<NSieveCol;idx++){
	for(UInt_t idy=0;idy<NSieveRow;idy++){

	  //	  if(((idx==4) && (idy==3))||((idx==5) && (idy==1)))
	  //	    {  xrange = 3.5;	  yrange = 1.4;}
	  //	  else { xrange = 2;	  yrange = 1;}

	  gXmin[idx][idy] = SieveEventID[nfoil][idx][idy][kRealSieveX] - xrange*SieveRadius;
	  gXmax[idx][idy] = SieveEventID[nfoil][idx][idy][kRealSieveX] + xrange*SieveRadius;
	  gYmin[idx][idy] = SieveEventID[nfoil][idx][idy][kRealSieveY] - yrange*SieveRadius+ yoffset;
	  gYmax[idx][idy] = SieveEventID[nfoil][idx][idy][kRealSieveY] + yrange*SieveRadius + yoffset;
	  gX[nfoil][idx][idy]=new TF1(Form("gX_F%d_C%d_R%d",nfoil,idx,idy), "gaus", gXmin[idx][idy],gXmax[idx][idy]);
	  gY[nfoil][idx][idy]=new TF1(Form("gY_F%d_C%d_R%d",nfoil,idx,idy), "gaus", gYmin[idx][idy],gYmax[idx][idy]);
	     //	  gX[nfoil[idx][idy]->SetParameter(3,2.0);
	  assert(gX[nfoil][idx][idy]);
	  assert(gY[nfoil][idx][idy]);
	  gX[nfoil][idx][idy]->SetLineColor(kRed);
	  gY[nfoil][idx][idy]->SetLineColor(kRed);
	}
      }
    }

    TCanvas * c1 = new TCanvas("SieveCheck Y","SieveCheck Y", 1800, 1100);
    c1->SetFillColor(kYellow-8);
    c1->SetFillStyle(3013);
    c1->Divide(3,3);
    for(UInt_t idx=0;idx<NFoils;idx++){
      c1->cd(idx+1)->SetFillColor(0);
      c1->cd(idx+1)->Divide(1,NSieveRow);
      for(UInt_t idy=0;idy<NSieveRow;idy++)
	{
	  c1->cd(idx+1)->cd(idy+1);
	  HSievePY[idx][idy]->Draw();
	}
    }
    c1->cd(6)->SetFillColor(0);
       
  //fit one-dimension projection
   for(UInt_t nfoil=0;nfoil<NFoils;nfoil++){
     for(UInt_t idx=0;idx<NSieveCol;idx++){
       for(UInt_t idy=0;idy<NSieveRow;idy++){
	 if (SieveEventID[nfoil][idx][idy][kEventID] > 0) {
	   c1->cd(nfoil+1)->cd(idy+1);
	   HSievePY[nfoil][idy]->Fit(Form("gY_F%d_C%d_R%d",nfoil,idx,idy),"R+");
	   SieveInfo[nfoil][idx][idy][kYmean] = HSievePY[nfoil][idy]->GetFunction(Form("gY_F%d_C%d_R%d",nfoil,idx,idy))->GetParameter(1);
	   SieveInfo[nfoil][idx][idy][kYmeanErr] = HSievePY[nfoil][idy]->GetFunction(Form("gY_F%d_C%d_R%d",nfoil,idx,idy))->GetParError(1);
	   SieveInfo[nfoil][idx][idy][kYamp] = HSievePY[nfoil][idy]->GetFunction(Form("gY_F%d_C%d_R%d",nfoil,idx,idy))->GetParameter(2);
	   SieveInfo[nfoil][idx][idy][kYampErr] = HSievePY[nfoil][idy]->GetFunction(Form("gY_F%d_C%d_R%d",nfoil,idx,idy))->GetParError(2);
	   SieveInfo[nfoil][idx][idy][kYndf] = HSievePY[nfoil][idy]->GetFunction(Form("gY_F%d_C%d_R%d",nfoil,idx,idy))->GetNDF();
	   SieveInfo[nfoil][idx][idy][kYchi2] = HSievePY[nfoil][idy]->GetFunction(Form("gY_F%d_C%d_R%d",nfoil,idx,idy))->GetChisquare();
	   //	   SieveInfo[nfoil][idx][idy][13] = SieveInfo[nfoil][idx][idy][1]-(SieveEventID[nfoil][idx][idy][kRealSieveY]);
	   SieveInfo[nfoil][idx][idy][kYdiff] = SieveInfo[nfoil][idx][idy][1]-(GetSieveHoleTCS(idx,idy).Y());
  	  	  
	 }
       }
     }
   }
    
   Double_t ymax=0;
   c1->Update();
   for(UInt_t nfoil=0;nfoil<NFoils;nfoil++){		  
     for(UInt_t idx=0;idx<NSieveRow;idx++)
       {
	 c1->cd(nfoil+1)->cd(idx+1);
	 ymax=c1->cd(nfoil+1)->cd(idx+1)->GetUymax();
	 //draw line at the position of sieves
	 for(UInt_t idy=0;idy<NSieveCol;idy++)
	   {
	     if (SieveEventID[nfoil][idy][idx][kEventID] > 0) {
	       TLine *line = new TLine( SieveEventID[nfoil][idy][idx][kRealSieveY],0, SieveEventID[nfoil][idy][idx][kRealSieveY],ymax);
	       line->SetLineColor(kBlue);
	       line->Draw("same");
	     }
	   }
       }
   }    

   /////////////////c2checkx////////////////
   
    TCanvas * c2 = new TCanvas("SieveCheck X","SieveCheck X", 1800, 1100);
    c2->SetFillColor(kYellow-8);
    c2->SetFillStyle(3013);
    c2->Divide(3,3);
    for(UInt_t idx=0;idx<NFoils;idx++){
      c2->cd(idx+1)->SetFillColor(0);
      c2->cd(idx+1)->Divide(1,NSieveCol);
      for(UInt_t idy=0;idy<NSieveCol;idy++)
	{
	  c2->cd(idx+1)->cd(idy+1);
	  HSievePX[idx][idy]->Draw();
	}
    }
    c2->cd(6)->SetFillColor(0);
    
   //fit one-dimension projection
    for(UInt_t nfoil=0;nfoil<NFoils;nfoil++){
      for(UInt_t idx=0;idx<NSieveCol;idx++){
	for(UInt_t idy=0;idy<NSieveRow;idy++){
	  if (SieveEventID[nfoil][idx][idy][kEventID] > 0) {
	    c2->cd(nfoil+1)->cd(idx+1);
	    HSievePX[nfoil][idx]->Fit(Form("gX_F%d_C%d_R%d",nfoil,idx,idy),"R+");
	    SieveInfo[nfoil][idx][idy][kXmean] = HSievePX[nfoil][idx]->GetFunction(Form("gX_F%d_C%d_R%d",nfoil,idx,idy))->GetParameter(1);
	    SieveInfo[nfoil][idx][idy][kXmeanErr] = HSievePX[nfoil][idx]->GetFunction(Form("gX_F%d_C%d_R%d",nfoil,idx,idy))->GetParError(1);
	    SieveInfo[nfoil][idx][idy][kXamp] = HSievePX[nfoil][idx]->GetFunction(Form("gX_F%d_C%d_R%d",nfoil,idx,idy))->GetParameter(2);
	    SieveInfo[nfoil][idx][idy][kXampErr] = HSievePX[nfoil][idx]->GetFunction(Form("gX_F%d_C%d_R%d",nfoil,idx,idy))->GetParError(2);
	    SieveInfo[nfoil][idx][idy][kXndf] = HSievePX[nfoil][idx]->GetFunction(Form("gX_F%d_C%d_R%d",nfoil,idx,idy))->GetNDF();
	    SieveInfo[nfoil][idx][idy][kXchi2] = HSievePX[nfoil][idx]->GetFunction(Form("gX_F%d_C%d_R%d",nfoil,idx,idy))->GetChisquare();
	    //	    SieveInfo[nfoil][idx][idy][14] = SieveInfo[nfoil][idx][idy][5]-(SieveEventID[nfoil][idx][idy][kRealSieveX]);
	    SieveInfo[nfoil][idx][idy][kXdiff] = SieveInfo[nfoil][idx][idy][5]-(GetSieveHoleTCS(idx,idy).X());
	  }
	}
      }
    }   
    
    c2->Update();
    for(UInt_t nfoil=0;nfoil<NFoils;nfoil++){
      for(UInt_t idx=0;idx<NSieveCol;idx++)
	{
	  c2->cd(nfoil+1)->cd(idx+1);
	  ymax= c2->cd(nfoil+1)->cd(idx+1)->GetUymax();
	  //draw line at the position of sieves
	  for(UInt_t idy=0;idy<NSieveRow;idy++)
	    {
	      if (SieveEventID[nfoil][idx][idy][kEventID] > 0) {
		TLine *line = new TLine( SieveEventID[nfoil][idx][idy][kRealSieveX],0, SieveEventID[nfoil][idx][idy][kRealSieveX],ymax);
		line->SetLineColor(kBlue);
		line->Draw();
	      }
	    }
	}
    }
     


     /////////////////////////////hole information drawing//////////////////////////////////
     TPaveText *tt;
     TCanvas *c3;
     for(UInt_t nfoil = 0; nfoil <NFoils; nfoil++) {
       if (nfoil ==0)
	 {   
	   gStyle->SetTitleW(0.3); 
	   gStyle->SetTitleH(0.18); 
	   gStyle->SetTitleX(0.7); 
		  }
       else if (nfoil ==1)
	 {   
	   gStyle->SetTitleW(0.3); 
	   gStyle->SetTitleH(0.18); 
	   gStyle->SetTitleX(0.7);
	 }
       else if (nfoil ==2)
	 {  
	   gStyle->SetTitleW(0.3); 
	   gStyle->SetTitleH(0.18); 
	   gStyle->SetTitleX(0.2);
	 }
       else if (nfoil ==3)
	 {  
	   gStyle->SetTitleW(0.3); 
	   gStyle->SetTitleH(0.18); 
	   gStyle->SetTitleX(0.2); 
	 }
       else 
	 { 
	   gStyle->SetTitleW(0.3); 
	   gStyle->SetTitleH(0.18); 
	   gStyle->SetTitleX(0.2); 
	 }

       c3 = new TCanvas(Form("SieveCheck Y #foil %d",nfoil),Form("SieveCheck Y #foil %d",nfoil), 1800, 1100);
       c3->Divide(1,NSieveRow);
       for(UInt_t idy=0;idy<NSieveRow;idy++)
	 {
	   c3->cd(idy+1);
	   HSievePY[nfoil][idy]->Draw();
	   HSievePY[nfoil][idy]->GetXaxis()->SetLabelSize(0.12);
	   HSievePY[nfoil][idy]->GetYaxis()->SetLabelSize(0.1);
	 }
       
       //fit one-dimension projection
       for(UInt_t idx=0;idx<NSieveCol;idx++){
	 for(UInt_t idy=0;idy<NSieveRow;idy++){
	   if (SieveEventID[nfoil][idx][idy][kEventID] > 0) {
	     c3->cd(idy+1);
	     HSievePY[nfoil][idy]->Fit(Form("gY_F%d_C%d_R%d",nfoil,idx,idy),"R+");  	  	  
	   }
	 }
       }
    
       ymax=0;
       c3->Update();
       for(UInt_t idx=0;idx<NSieveRow;idx++)
	 {
	   c3->cd(idx+1);
	   ymax=c3->cd(idx+1)->GetUymax();
	   //draw line at the position of sieves
	   for(UInt_t idy=0;idy<NSieveCol;idy++)
	     {
	       if (SieveEventID[nfoil][idy][idx][kEventID] > 0) {
		 TLine *line = new TLine( SieveEventID[nfoil][idy][idx][kRealSieveY],0, SieveEventID[nfoil][idy][idx][kRealSieveY],ymax);
		 line->SetLineColor(kBlue);
		 line->Draw("same");
		 tt = new TPaveText(SieveEventID[nfoil][idy][idx][kRealSieveY]-0.003, 0.95*ymax,SieveEventID[nfoil][idy][idx][kRealSieveY]+0.0025,ymax*1.13);
		 //	       tt->SetLineColor(0);
		 tt->SetFillColor(0);
		 tt->SetShadowColor(0);
		 tt->AddText( Form( "%1.2f",1000*(SieveEventID[nfoil][idy][idx][kRealSieveY]- SieveInfo[nfoil][idy][idx][kYmean]) ) );
		 tt->Draw();
	       }
	     }
	 }
     }

     ///////////////////////////////////////////////////////////////
     UInt_t y1,y2;
     for(UInt_t nfoil = 0; nfoil <NFoils; nfoil++) {
       if (nfoil ==0)
	 { y1=0+17;y2=9+17; }
       else if(nfoil ==1)
	 { y1=0+14;y2=9+14; }
       else if(nfoil ==2)
	 { y1=0+11;y2=9+11; }
       else if(nfoil ==3)
	 { y1=0+9;y2=9+9; }
       else if(nfoil ==4)
	 { y1=0+8;y2=9+8; }
       else if (nfoil ==5)
	 { y1=0+7;y2=9+7;}
       else if (nfoil ==6)
	 { y1=0+5;y2=9+5;}
       else if (nfoil ==7)
	 { y1=0+3;y2=9+3;}
       else
	 { y1=0;y2=9;}
 
       c3 = new TCanvas(Form("SieveCheck X #foil %d",nfoil),Form("SieveCheck X #foil %d",nfoil), 1800, 1100);
       //       c3->Divide(1,NSieveCol);
       c3->Divide(1,10);
       //       for(UInt_t idy=0;idy<NSieveCol;idy++)
       for(UInt_t idy=y1;idy<y2;idy++)
	 {
	   c3->cd(idy+1-y1);
	   HSievePX[nfoil][idy]->Draw();
	   HSievePX[nfoil][idy]->GetXaxis()->SetLabelSize(0.12);
	   HSievePX[nfoil][idy]->GetYaxis()->SetLabelSize(0.1);
	 }
       
       //fit one-dimension projection
       //       for(UInt_t idx=0;idx<NSieveCol;idx++){
       for(UInt_t idx=y1;idx<y2;idx++){
	 for(UInt_t idy=0;idy<NSieveRow;idy++){
	   if (SieveEventID[nfoil][idx][idy][kEventID] > 0) {
	     c3->cd(idx+1-y1);
	     HSievePX[nfoil][idx]->Fit(Form("gX_F%d_C%d_R%d",nfoil,idx,idy),"R+");  	  	  
	   }
	 }
       }
    
       ymax=0;
       c3->Update();
       //       for(UInt_t idx=0;idx<NSieveRow;idx++)
       //       for(UInt_t idy=0;idy<NSieveCol;idy++)       
       for(UInt_t idy=y1;idy<y2;idy++)
	 {
	   c3->cd(idy+1-y1);
	   ymax=c3->cd(idy+1-y1)->GetUymax();
	   //draw line at the position of sieves
	   //	   for(UInt_t idy=0;idy<NSieveCol;idy++)
	   for(UInt_t idx=0;idx<NSieveRow;idx++)
	     {
	       if (SieveEventID[nfoil][idy][idx][kEventID] > 0) {
		 TLine *line = new TLine( SieveEventID[nfoil][idy][idx][kRealSieveX],0, SieveEventID[nfoil][idy][idx][kRealSieveX],ymax);
		 line->SetLineColor(kBlue);
		 line->Draw("same");
		 tt = new TPaveText(SieveEventID[nfoil][idy][idx][kRealSieveX]-0.007, 0.95*ymax,SieveEventID[nfoil][idy][idx][kRealSieveX]+0.007,ymax*1.13);
		 //	       tt->SetLineColor(0);
		 tt->SetFillColor(0);
		 tt->SetShadowColor(0);
		 tt->AddText( Form( "%1.2f",1000*(SieveEventID[nfoil][idy][idx][kRealSieveX]- SieveInfo[nfoil][idy][idx][kXmean]) ) );
		 tt->Draw();
	       }
	     }
	 }
     }
      //////////////////////////////////hole information drawing////////////////////////////     

     //////////////////////////   
     
    FILE *fp2;
    fp2=fopen("HoleCorrection.dat","w");

    fprintf(fp2,"\n\n");

    for(UInt_t nfoil=0;nfoil<NFoils;nfoil++){
      fprintf(fp2,Form(" The horizontal correction of each hole in Foil#%d.\n",nfoil+1));
      for (UInt_t idx = 0; idx < NSieveRow; idx++) {
	for (UInt_t idy = 0; idy < NSieveCol; idy++) {
	  fprintf(fp2,"%*.2f", 6,1000*(GetSieveHoleCorrectionTCS(nfoil, NSieveCol-1-idy, NSieveRow-1-idx).Y()- GetSieveHoleTCS(NSieveCol-1-idy, NSieveRow-1-idx).Y()) );
	}
	fprintf(fp2,"\n");
      }
      fprintf(fp2,"\n");
    }
    fprintf(fp2,"\n\n");
    for(UInt_t nfoil=0;nfoil<NFoils;nfoil++){
      fprintf(fp2,Form(" The vertical correction of each hole in Foil#%d.\n",nfoil+1));
      for (UInt_t idx = 0; idx < NSieveRow; idx++) {
	for (UInt_t idy = 0; idy < NSieveCol; idy++) {
	  fprintf(fp2,"%*.2f  ", 6,1000*(GetSieveHoleCorrectionTCS(nfoil, NSieveCol-1-idy, NSieveRow-1-idx).X()-GetSieveHoleTCS(NSieveCol-1-idy, NSieveRow-1-idx).X()) );
	}
	fprintf(fp2,"\n");
      }
      fprintf(fp2,"\n");
    }
    fclose(fp2);
     
    ///////////////////////////////////

     ////////////////////////////////////
         
    UInt_t NHole[NFoils] = {0}; // array to store the number of hole for each foil target
    UInt_t NHOLE =0; // the number of total holes
    FILE *fp1;
    fp1=fopen("Opt1\HoleInfo.dat","w");
    fprintf(fp1,"\n");

    for(UInt_t nfoil=0; nfoil < NFoils; nfoil++){
      fprintf(fp1,Form("      The number of events in each hole in Foil#%d.\n",nfoil+1));
      for (UInt_t idx = 0; idx < NSieveRow; idx++) {
	for (UInt_t idy = 0; idy < NSieveCol; idy++) {
	  fprintf(fp1,"%5.0f  ", SieveInfo[nfoil][NSieveCol-1-idy][NSieveRow-1-idx][0]);
	  if(SieveInfo[nfoil][NSieveCol-1-idy][NSieveRow-1-idx][0]>0){
	    NHole[nfoil]++;
	    NHOLE++;
	  }
	}
	fprintf(fp1,"\n");
      }
      fprintf(fp1,"\n");
    }
    cout<<"NHOLE:"<<NHOLE<<endl;
    for(UInt_t i=0;i<NFoils;i++){
      cout<<"NHole["<<i<<"]:"<<NHole[i]<<endl;
    }

    Double_t rms_hole_v[NFoils] = {.0}, rms_hole_h[NFoils] = {.0}; // the square of (fitting center - corrected hole position) for each target foil
    Double_t RMS_hole_v=0, RMS_hole_h=0; // the square of (fitting center - corrected hole position) for all target foils
    Double_t a1=0,a2=0;    
    fprintf(fp1,"\n\n");
    
    for(UInt_t nfoil=0;nfoil<NFoils;nfoil++){
      fprintf(fp1,Form(" The difference of horizontal fitting mean and real position for each hole in Foil#%d.\n",nfoil+1));
      //SieveYbyCol[]
      for(UInt_t ii=0; ii< NSieveCol; ii++){
	  fprintf(fp1,"%-15.2f",1000*SieveYbyCol[NSieveCol-1-ii]);
	  fprintf(fp1," ");
      }
      fprintf(fp1,"\n");
      fprintf(fp1,"\n");
      
      for (UInt_t idx = 0; idx < NSieveRow; idx++) {
	for (UInt_t idy = 0; idy < NSieveCol; idy++) {
	  a1 = GetSieveHoleTCS(NSieveCol-1-idy,NSieveRow-1-idx).Y();
	  a2 = GetSieveHoleCorrectionTCS(nfoil,NSieveCol-1-idy,NSieveRow-1-idx).Y();
	  if(SieveInfo[nfoil][NSieveCol-1-idy][NSieveRow-1-idx][0] > 0){
	    fprintf(fp1,"%-6.2f/%-6.2f  ",1000*(SieveInfo[nfoil][NSieveCol-1-idy][NSieveRow-1-idx][13]),1000*(a2-a1));
	    rms_hole_h[nfoil] += TMath::Power(SieveInfo[nfoil][NSieveCol-1-idy][NSieveRow-1-idx][13] - a2 + a1,2);
	    RMS_hole_h += TMath::Power(SieveInfo[nfoil][NSieveCol-1-idy][NSieveRow-1-idx][13] - a2 + a1,2);
	  }
	  else {fprintf(fp1,"x.xx  /x.xx    ");}
	  //	  fprintf(fp1,"%-6.2f/%-6.2f  ",1000*SieveInfo[nfoil][8-idy][6-idx][13],1000*SieveInfo[nfoil][8-idy][6-idx][1]);
	  fprintf(fp1," ");
	}
	fprintf(fp1,"\n");
      }
      fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n\n");
      
    for(UInt_t nfoil=0;nfoil<NFoils;nfoil++){
      fprintf(fp1,Form(" The difference of vertical fitting mean and real position for each hole in Foil#%d.\n",nfoil+1));
      for (UInt_t idx = 0; idx < NSieveRow; idx++) {
	//	if(Col%2 ==0) fprintf(fp1,"%-6.2f",1000*SieveXbyRowEven[NSieveRow-1-idx]);
	//	else fprintf(fp1,"%-6.2f",1000*SieveXbyRowOdd[NSieveRow-1-idx]);
	fprintf(fp1,"  ");
	for (UInt_t idy = 0; idy < NSieveCol; idy++) {
	  //	  fprintf(fp1,"%-6.2f/%-6.2f  ", 1000*SieveInfo[nfoil][NSieveCol-1-idy][NSieveRow-1-idx][14],1000*SieveInfo[nfoil][NSieveCol-1-idy][NSieveRow-1-idx][5]);
	  a1 = GetSieveHoleTCS(NSieveCol-1-idy,NSieveRow-1-idx).X();
	  a2 = GetSieveHoleCorrectionTCS(nfoil,NSieveCol-1-idy,NSieveRow-1-idx).X();
	  if(SieveInfo[nfoil][NSieveCol-1-idy][NSieveRow-1-idx][0] >0){
	    fprintf(fp1,"%-6.2f/%-6.2f  ", 1000*SieveInfo[nfoil][NSieveCol-1-idy][NSieveRow-1-idx][14],1000*(a2 - a1));
	    rms_hole_v[nfoil] += TMath::Power(SieveInfo[nfoil][NSieveCol-1-idy][NSieveRow-1-idx][14] - a2 + a1,2);
	    RMS_hole_v +=TMath::Power(SieveInfo[nfoil][NSieveCol-1-idy][NSieveRow-1-idx][14] - a2 + a1,2);
	  }
	  else
	    {fprintf(fp1,"x.xx  /x.xx    ");}
	  fprintf(fp1," ");
	}
	fprintf(fp1,"\n");
      }
      fprintf(fp1,"\n");
    }
    
    for(UInt_t i=0;i<NFoils;i++){
      rms_hole_h[i] = 1000*TMath::Sqrt(rms_hole_h[i]/NHole[i]);
      cout<<"rms_hole_h["<<i<<"]:"<<rms_hole_h[i]<<endl;
    }
    RMS_hole_h = 1000*TMath::Sqrt(RMS_hole_h/NHOLE);
    cout<<"RMS_hole_h:"<<RMS_hole_h<<endl;
    cout<<"*******************"<<endl;
    for(UInt_t i=0;i<NFoils;i++){
      rms_hole_v[i] = 1000*TMath::Sqrt(rms_hole_v[i]/NHole[i]);
      cout<<"rms_hole_v["<<i<<"]:"<<rms_hole_v[i]<<endl;
    }
    RMS_hole_v = 1000*TMath::Sqrt(RMS_hole_v/NHOLE);
    cout<<"RMS_hole_v:"<<RMS_hole_v<<endl;
    
    for(UInt_t i=0;i<NFoils;i++){
      fprintf(fp1,"NHOle[%d]:%d \n",i+1,NHole[i]);
    }
    fprintf(fp1,"The total number of visible holes: %d \n", NHOLE);
    fprintf(fp1,"******************\n");
    for(UInt_t i=0;i<NFoils;i++){
      fprintf(fp1,"rms_hole_h[%d]:%.2f \n",i+1,rms_hole_h[i]);
    }
    fprintf(fp1,"RMS_hole_h:%.2f\n",RMS_hole_h);
    fprintf(fp1,"******************\n");
    for(UInt_t i=0;i<NFoils;i++){
      fprintf(fp1,"rms_hole_v[%d]:%.2f \n",i+1,rms_hole_v[i]);
    }
    fprintf(fp1,"RMS_hole_v:%.2f\n",RMS_hole_v);
    
    fclose(fp1);
     
    ////////////////////////////////////////////////////    
    return c1;

}



TCanvas * LOpticsOpt::Sieve_hole_diff(Int_t PlotFoilID){

  // create plots for each foil of difference between expected and real sieve x and y 


   const UInt_t nplot = PlotFoilID;


   
    // theta and phi diffs for each foil
   Double_t dTh[NFoils] = {0}, dPh[NFoils] = {0};

    // sieve x and y diffs for each foil
   Double_t dX[NFoils] = {0}, dY[NFoils] = {0};

    enum {
      kEventID, kRealSieveX, kRealSieveY, kCalcSieveX, kCalcSieveY, kRealSieveTh, kRealSievePh, KCalcSieveTh, kCalcSievePh
    };


    TH1D * HSieveThDiff[NFoils] = {0};
    TH1D * HSievePhDiff[NFoils] = {0};
    Double_t th_lim[NFoils] = {0};
    Double_t ph_lim[NFoils] = {0};

    
    TH1D * HSieveXDiff[NFoils] = {0};
    TH1D * HSieveYDiff[NFoils] = {0};
    TH1D * HSieveYDiff_rem[NFoils] = {0};


    
    Double_t x_lim[NFoils] = {0};
    Double_t y_lim[NFoils] = {0};


    Double_t th_lim_h = 0.08;
    Double_t th_lim_l = -0.08;
    Double_t ph_lim_h = 0.00;
    Double_t ph_lim_l = -0.08;


    for (UInt_t idx = 0; idx < NFoils; idx++) {
      //      x_lim[idx] = 6.5;
      x_lim[idx] = 16.5;
        //y_lim[idx] = 1.5 * TMath::Max(TMath::Abs(SieveXbyRowOdd[0]), TMath::Abs(SieveXbyRowEven[NSieveRow - 1]));
      // y_lim[idx] = 4;
      y_lim[idx] = 14;

      
      HSieveThDiff[idx] = new TH1D(Form("Th Diff Sieve_Foil%d", idx), "\\theta Diff for " + Foil_names[idx], 200, th_lim_l, th_lim_h);
      HSieveThDiff[idx]->SetXTitle("Sieve Th [mrad]");
      HSievePhDiff[idx] = new TH1D(Form("Ph Diff Sieve_Foil%d", idx), Form("Ph Diff #%d", idx), 200, ph_lim_l, ph_lim_h);
      HSievePhDiff[idx]->SetXTitle("Sieve Ph [mrad]");
      

      HSieveXDiff[idx] = new TH1D(Form("X Diff Sieve_Foil%d", idx), "X Diff for " + Foil_names[idx], 150, -x_lim[idx], x_lim[idx]);
      HSieveXDiff[idx]->SetXTitle("Reconstructed - survey sieve X [mm]");
      HSieveXDiff[idx]->GetXaxis()->CenterTitle();
      HSieveYDiff[idx]  = new TH1D(Form("Y Diff Sieve_Foil%d", idx), "Y Diff for " + Foil_names[idx], 200, -y_lim[idx], y_lim[idx]);
      HSieveYDiff[idx]->SetXTitle("Reconstructed - survey sieve Y [mm]");
      HSieveYDiff[idx]->GetXaxis()->CenterTitle();


      // assert(HSievePlaneDiff[idx]); // assure memory     
      // assert(HSieveAngleDiff[idx]); // assure memory allocation
    }


    Int_t PlotID; // useful for case when only Vertical foils are used
                    // (FoilID based on scheme with optics and vertical foils where optics foils
                    // are 0-7 and vertical foils 8-10) so they can have PlotID 0,1,2
    

    // fill difference histograms

    for (UInt_t idx = 0; idx < fNRawData; idx++) {
      const EventData &eventdata = fRawData[idx];
    

      UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID]; //starting 0!
      // res = res%(NSieveRow*NSieveCol);


      

      const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
      const UInt_t Row = (UInt_t)eventdata.Data[kRowID]; //starting 0!
      
      
      assert(FoilID < NFoils); //array index check
      
      const TVector3 SieveHoleCorrectionTCS = GetSieveHoleCorrectionTCS(FoilID, Col, Row);

      PlotID = FoilID; // useful for case when only Vertical foils are used
                    // (FoilID based on scheme with optics and vertical foils where optics foils
                    // are 0-7 and vertical foils 8-10) so they can have PlotID 0,1,2
      if(NFoils<11){

	PlotID -= 8;
      }

      
      // calculate sieve_x & sieve_y
      
      Double_t ProjectionX = eventdata.Data[kRealTgX] + (eventdata.Data[kCalcTh] + eventdata.Data[kRealTgX] * ExtTarCor_ThetaCorr) * (SieveHoleCorrectionTCS.Z());
      
      
      //Double_t ProjectionX = eventdata.Data[kCalcTh];
      Double_t ProjectionY = eventdata.Data[kRealTgY] + eventdata.Data[kCalcPh] * (SieveHoleCorrectionTCS.Z());
      
	

      // calculate theta and phi 'real' positions of holes

      //      TVector3 BeamSpotHCS_average(BeamX_average,BeamY_average,targetfoils[FoilID]);
      const TVector3 BeamSpotHCS_average(BeamX_average[FoilID] + (targetfoils[FoilID]/BeamZDir_average)*BeamXDir_average[FoilID], BeamY_average[FoilID] + (targetfoils[FoilID]/BeamZDir_average)*BeamYDir_average[FoilID], targetfoils[FoilID]);
    
      
      TVector3 BeamSpotTCS_average = fTCSInHCS.Inverse()*(BeamSpotHCS_average-fPointingOffset);





      TVector3 MomDirectionTCS_hole = SieveHoleCorrectionTCS - BeamSpotTCS_average;
      
      Double_t theta_hole = MomDirectionTCS_hole.X()/MomDirectionTCS_hole.Z();
      Double_t phi_hole = MomDirectionTCS_hole.Y()/MomDirectionTCS_hole.Z();
      

      dX[PlotID] = ProjectionX - SieveHoleCorrectionTCS.X();
      dY[PlotID] = ProjectionY - SieveHoleCorrectionTCS.Y();

      dTh[PlotID] = eventdata.Data[kCalcTh] - theta_hole;
      dPh[PlotID] = eventdata.Data[kCalcPh] - phi_hole;
      

      HSieveXDiff[PlotID]->Fill(1000*dX[PlotID]);
      HSieveYDiff[PlotID]->Fill(1000*dY[PlotID]);


    }

    
    // fit histograms to get widths of holes

    Int_t max_bin;
    Double_t max, min;

    Double_t hist_mean;
    Double_t mean_offset[NFoils]; // parameter to judge differnce between gaus and histogram mean

    Double_t width_guess = 3.0; // intial range of first fit


    // varaibles to collect output parameters of fit
    Double_t x_mean,x_wid;
    Double_t x_means[NFoils],x_widths[NFoils];

    Double_t y_mean,y_wid;
    Double_t y_means[NFoils],y_widths[NFoils];
    

    // parameters of fit of Gauss and linear tail for y distrib
    Double_t lin_offset, lin_coeff;
    Double_t lin_start, lin_end;
    Double_t lin_offsets[NFoils],lin_coeffs[NFoils];


    // cutoff for events (based on y diff) being in main gaus of y diff or linear tail
    // Double_t gaus_start[NFoils], gaus_end[NFoils];
    // Double_t lin_starts[NFoils], lin_ends[NFoils];

    
    for (UInt_t idx = 0; idx < NFoils; idx++) {


      max_bin = HSieveXDiff[idx]->GetMaximumBin();
      min = HSieveXDiff[idx]->GetBinCenter(max_bin) - width_guess;
      max = HSieveXDiff[idx]->GetBinCenter(max_bin) + width_guess;
      
      cout << "pre-first-fitting stage" << endl;
      HSieveXDiff[idx]->Fit("gaus","Q","",min,max);     

      if(HSieveXDiff[idx]->GetFunction("gaus")){
	cout << "Fit worked!" << endl;
            
	TF1* xgaus1 = HSieveXDiff[idx]->GetFunction("gaus");
      
	min = xgaus1->GetParameter(1) - 3*xgaus1->GetParameter(2);
	max = xgaus1->GetParameter(1) + 3*xgaus1->GetParameter(2);
      
	HSieveXDiff[idx]->Fit("gaus","Q","",min,max);
	TF1* xgaus2 = HSieveXDiff[idx]->GetFunction("gaus");
                 
	x_mean = xgaus2->GetParameter(1);
	x_wid = xgaus2->GetParameter(2);                
	x_means[idx] = x_mean;
	x_widths[idx] = x_wid;

       	

      }


      max_bin = HSieveYDiff[idx]->GetMaximumBin();
      min = HSieveYDiff[idx]->GetBinCenter(max_bin) - width_guess;
      max = HSieveYDiff[idx]->GetBinCenter(max_bin) + width_guess;

      HSieveYDiff[idx]->Fit("gaus","Q","",min,max);

      
      if(HSieveYDiff[idx]->GetFunction("gaus")){
	cout << "Y fit worked" << endl;


	TF1* ygaus1 = HSieveYDiff[idx]->GetFunction("gaus");
      
	min = ygaus1->GetParameter(1) - 3*ygaus1->GetParameter(2);
	max = ygaus1->GetParameter(1) + 3*ygaus1->GetParameter(2);
      
	HSieveYDiff[idx]->Fit("gaus","Q","",min,max);
	TF1* ygaus2 = HSieveYDiff[idx]->GetFunction("gaus");
                 
	y_mean = ygaus2->GetParameter(1);
	y_wid = ygaus2->GetParameter(2);                
	y_means[idx] = y_mean;
	y_widths[idx] = y_wid;

	hist_mean = HSieveYDiff[idx]->GetMean();

	cout << "For Foil " << idx << ": " <<  endl;
	cout << "gaus_mean = " << y_mean << endl;
	cout << "gaus_width = " << y_wid << endl;
	cout << "hist_mean =  " << hist_mean << endl;
	mean_offset[idx] = TMath::Abs(y_mean-hist_mean)/y_wid;
	cout << "|(gaus_mean-hist_mean)|/gaus_width = " << mean_offset[idx] << endl;

	
      }

    }


        // create canvas for x sieve difference
    TCanvas * c1;
    
    if (nplot <= 1) {
      c1 = new TCanvas("XDiff", "XDiff", 1000, 1000);
      c1->Divide(1, 1);
    } else if (nplot <= 3) {
      c1 = new TCanvas("XDiff", "XDiff", 1800, 1100);
      c1->Divide(3, 1);
    } else if (nplot <= 6) {
      c1 = new TCanvas("XDiff", "XDiff", 1800, 1100);
      c1->Divide(3, 2);
    } else if (nplot <= 9){
      c1 = new TCanvas("XDiff", "XDiff", 1800, 1100);
      c1->Divide(3, 3);
    }
    else{
      c1 = new TCanvas("XDiff", "XDiff", 1800, 1100);
      c1->Divide(4, 3);
    }


    // create canvas for y sieve difference
    TCanvas * c2;
    
    if (nplot <= 1) {
      c2 = new TCanvas("YDiff", "YDiff", 1000, 1000);
      c2->Divide(1, 1);
    } else if (nplot <= 3) {
      c2 = new TCanvas("YDiff", "YDiff", 1800, 1100);
      c2->Divide(3, 1);
    } else if (nplot <= 6) {
      c2 = new TCanvas("YDiff", "YDiff", 1800, 1100);
      c2->Divide(3, 2);
    } else if (nplot <= 9){
      c2 = new TCanvas("YDiff", "YDiff", 1800, 1100);
      c2->Divide(3, 3);
    }
    else{
      c2 = new TCanvas("YDiff", "YDiff", 1800, 1100);
      c2->Divide(4, 3);
    }



      

    for (UInt_t idx = 0; idx < nplot; idx++) {


      c1->cd(idx+1);
     
      HSieveXDiff[idx]->Draw();

      // TPaveText *ttx = new TPaveText(0.6,0.7,0.88,0.82,"NDC");
      // ttx->AddText(Form("\\Delta = %2.3f mm",  x_means[idx]));
      // ttx->AddText(Form("\\sigma = %2.3f mm", x_widths[idx]));
      // ttx->SetShadowColor(0);
      // ttx->SetFillColor(0);
      // ttx->Draw("same");


      TPaveText *ttx_1 = new TPaveText(0.6,0.75,0.85,0.82,"NDC");
      ttx_1->AddText("\\sigma = ");
      ttx_1->AddText(Form("%2.3f mm", x_widths[idx]));


      ttx_1->SetShadowColor(0);
      ttx_1->SetFillColor(0);
      ttx_1->SetTextSize(0.08);
      ttx_1->Draw("same");


      TPaveText *ttx_2 = new TPaveText(0.25,0.75,0.28,0.82,"NDC");

      ttx_2->AddText("\\Delta = ");
      ttx_2->AddText(Form("%2.3f mm",  x_means[idx]));

      ttx_2->SetShadowColor(0);
      ttx_2->SetFillColor(0);
      ttx_2->SetTextSize(0.08);
      ttx_2->Draw("same");


      c2->cd(idx+1);
     
      HSieveYDiff[idx]->Draw();
           

      TPaveText *tty_1 = new TPaveText(0.6,0.75,0.85,0.82,"NDC");
      //      tty_1->AddText(Form("\\Delta = %2.3f mm",  y_means[idx]));
      tty_1->AddText("\\sigma = ");
      tty_1->AddText(Form("%2.3f mm", y_widths[idx]));


      tty_1->SetShadowColor(0);
      tty_1->SetFillColor(0);
      tty_1->SetTextSize(0.08);
      tty_1->Draw("same");


      TPaveText *tty_2 = new TPaveText(0.25,0.75,0.28,0.82,"NDC");

      //      tty_2->AddText(Form("\\Delta = %2.3f mm",  y_means[idx]));
      tty_2->AddText("\\Delta = ");
      tty_2->AddText(Form("%2.3f mm",  y_means[idx]));

      
      //tty_2->AddText(Form("\\sigma = %2.3f mm", y_widths[idx]));
      tty_2->SetShadowColor(0);
      tty_2->SetFillColor(0);
      tty_2->SetTextSize(0.08);
      tty_2->Draw("same");

      cout << "loops or c1 and c2: " << idx << endl;
      
    }       
      
    cout << "c1 returning..." << endl;
    return c1;
    cout << "c1 returned" << endl;


}



TCanvas * LOpticsOpt::Sieve_hole_diff_tail(Int_t PlotFoilID){

  // create plots for each foil of difference between expected and real sieve x and y 


   const UInt_t nplot = PlotFoilID;


   
    // theta and phi diffs for each foil
   Double_t dTh[NFoils] = {0}, dPh[NFoils] = {0};

    // sieve x and y diffs for each foil
   Double_t dX[NFoils] = {0}, dY[NFoils] = {0};

    enum {
      kEventID, kRealSieveX, kRealSieveY, kCalcSieveX, kCalcSieveY, kRealSieveTh, kRealSievePh, KCalcSieveTh, kCalcSievePh
    };


    TH1D * HSieveThDiff[NFoils] = {0};
    TH1D * HSievePhDiff[NFoils] = {0};
    Double_t th_lim[NFoils] = {0};
    Double_t ph_lim[NFoils] = {0};

    
    TH1D * HSieveXDiff[NFoils] = {0};
    TH1D * HSieveYDiff[NFoils] = {0};
    TH1D * HSieveYDiff_rem[NFoils] = {0};


    // series of histogram used to show difference between
    // events in Gaus and tail of y diff 

    TH1D* H_FP_X_tail[NFoils] = {0};
    TH1D* H_FP_X_gaus[NFoils] = {0};

    TH1D* H_FP_Y_tail[NFoils] = {0};
    TH1D* H_FP_Y_gaus[NFoils] = {0};

    TH1D* H_FP_P_tail[NFoils] = {0};
    TH1D* H_FP_P_gaus[NFoils] = {0};

    TH1D* H_FP_T_tail[NFoils] = {0};
    TH1D* H_FP_T_gaus[NFoils] = {0};

    TH2D* H_FP_TP_tail[NFoils] = {0}; // t vs p
    TH2D* H_FP_TP_gaus[NFoils] = {0}; // t vs p

    TH2D* H_FP_YP_tail[NFoils] = {0}; // y vs p
    TH2D* H_FP_YP_gaus[NFoils] = {0}; // y vs p

    // beam parameters

    TH2D* H_beam_xy_tail[NFoils] = {0}; // beam y vs x
    TH2D* H_beam_xy_gaus[NFoils] = {0}; // beam y vs x
    
    
    TF1  * y_bg[NFoils] = {0};
    TF1  * y_gaus[NFoils] = {0};

    
    Double_t x_lim[NFoils] = {0};
    Double_t y_lim[NFoils] = {0};


    Double_t th_lim_h = 0.08;
    Double_t th_lim_l = -0.08;
    Double_t ph_lim_h = 0.00;
    Double_t ph_lim_l = -0.08;


    for (UInt_t idx = 0; idx < NFoils; idx++) {
      //      x_lim[idx] = 6.5;
      x_lim[idx] = 16.5;
        //y_lim[idx] = 1.5 * TMath::Max(TMath::Abs(SieveXbyRowOdd[0]), TMath::Abs(SieveXbyRowEven[NSieveRow - 1]));
      // y_lim[idx] = 4;
      y_lim[idx] = 14;

      
      HSieveThDiff[idx] = new TH1D(Form("Th Diff Sieve_Foil%d", idx), "\\theta Diff for " + Foil_names[idx], 200, th_lim_l, th_lim_h);
      HSieveThDiff[idx]->SetXTitle("Sieve Th [mrad]");
      HSievePhDiff[idx] = new TH1D(Form("Ph Diff Sieve_Foil%d", idx), Form("Ph Diff #%d", idx), 200, ph_lim_l, ph_lim_h);
      HSievePhDiff[idx]->SetXTitle("Sieve Ph [mrad]");
      

      HSieveXDiff[idx] = new TH1D(Form("X Diff Sieve_Foil%d", idx), "X Diff for " + Foil_names[idx], 150, -x_lim[idx], x_lim[idx]);
      HSieveXDiff[idx]->SetXTitle("Reconstructed - survey sieve X [mm]");
      HSieveXDiff[idx]->GetXaxis()->CenterTitle();
      HSieveYDiff[idx]  = new TH1D(Form("Y Diff Sieve_Foil%d", idx), "Y Diff for " + Foil_names[idx], 200, -y_lim[idx], y_lim[idx]);
      HSieveYDiff[idx]->SetXTitle("Reconstructed - survey sieve Y [mm]");
      HSieveYDiff[idx]->GetXaxis()->CenterTitle();


      // assert(HSievePlaneDiff[idx]); // assure memory     
      // assert(HSieveAngleDiff[idx]); // assure memory allocation
    }


    Int_t PlotID; // useful for case when only Vertical foils are used
                    // (FoilID based on scheme with optics and vertical foils where optics foils
                    // are 0-7 and vertical foils 8-10) so they can have PlotID 0,1,2
    

    // fill difference histograms

    for (UInt_t idx = 0; idx < fNRawData; idx++) {
      const EventData &eventdata = fRawData[idx];
    

      UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID]; //starting 0!
      // res = res%(NSieveRow*NSieveCol);


      

      const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
      const UInt_t Row = (UInt_t)eventdata.Data[kRowID]; //starting 0!
      
      
      assert(FoilID < NFoils); //array index check
      
      const TVector3 SieveHoleCorrectionTCS = GetSieveHoleCorrectionTCS(FoilID, Col, Row);

      PlotID = FoilID; // useful for case when only Vertical foils are used
                    // (FoilID based on scheme with optics and vertical foils where optics foils
                    // are 0-7 and vertical foils 8-10) so they can have PlotID 0,1,2
      if(NFoils<11){

	PlotID -= 8;
      }

      
      // calculate sieve_x & sieve_y
      
      Double_t ProjectionX = eventdata.Data[kRealTgX] + (eventdata.Data[kCalcTh] + eventdata.Data[kRealTgX] * ExtTarCor_ThetaCorr) * (SieveHoleCorrectionTCS.Z());
      
      
      //Double_t ProjectionX = eventdata.Data[kCalcTh];
      Double_t ProjectionY = eventdata.Data[kRealTgY] + eventdata.Data[kCalcPh] * (SieveHoleCorrectionTCS.Z());
      
	
      //Double_t ProjectionY = eventdata.Data[kCalcPh];
	

      // calculate theta and phi 'real' positions of holes

      //      TVector3 BeamSpotHCS_average(BeamX_average,BeamY_average,targetfoils[FoilID]);
      const TVector3 BeamSpotHCS_average(BeamX_average[FoilID] + (targetfoils[FoilID]/BeamZDir_average)*BeamXDir_average[FoilID], BeamY_average[FoilID] + (targetfoils[FoilID]/BeamZDir_average)*BeamYDir_average[FoilID], targetfoils[FoilID]);
    
      
      TVector3 BeamSpotTCS_average = fTCSInHCS.Inverse()*(BeamSpotHCS_average-fPointingOffset);





      TVector3 MomDirectionTCS_hole = SieveHoleCorrectionTCS - BeamSpotTCS_average;
      
      Double_t theta_hole = MomDirectionTCS_hole.X()/MomDirectionTCS_hole.Z();
      Double_t phi_hole = MomDirectionTCS_hole.Y()/MomDirectionTCS_hole.Z();
      

      dX[PlotID] = ProjectionX - SieveHoleCorrectionTCS.X();
      dY[PlotID] = ProjectionY - SieveHoleCorrectionTCS.Y();

      dTh[PlotID] = eventdata.Data[kCalcTh] - theta_hole;
      dPh[PlotID] = eventdata.Data[kCalcPh] - phi_hole;
      

      HSieveXDiff[PlotID]->Fill(1000*dX[PlotID]);
      HSieveYDiff[PlotID]->Fill(1000*dY[PlotID]);


    }

    
    // fit histograms to get widths of holes

    Int_t max_bin;
    Double_t max, min;

    Double_t hist_mean;
    Double_t mean_offset[NFoils]; // parameter to judge differnce between gaus and histogram mean

    Double_t width_guess = 3.0; // intial range of first fit


    // varaibles to collect output parameters of fit
    Double_t x_mean,x_wid;
    Double_t x_means[NFoils],x_widths[NFoils];

    Double_t y_mean,y_wid;
    Double_t y_means[NFoils],y_widths[NFoils];
    

    // parameters of fit of Gauss and linear tail for y distrib
    Double_t lin_offset, lin_coeff;
    Double_t lin_start, lin_end;
    Double_t lin_offsets[NFoils],lin_coeffs[NFoils];


    // cutoff for events (based on y diff) being in main gaus of y diff or linear tail
    // Double_t gaus_start[NFoils], gaus_end[NFoils];
    // Double_t lin_starts[NFoils], lin_ends[NFoils];

    
    for (UInt_t idx = 0; idx < NFoils; idx++) {


      max_bin = HSieveXDiff[idx]->GetMaximumBin();
      min = HSieveXDiff[idx]->GetBinCenter(max_bin) - width_guess;
      max = HSieveXDiff[idx]->GetBinCenter(max_bin) + width_guess;
      
      cout << "pre-first-fitting stage" << endl;
      HSieveXDiff[idx]->Fit("gaus","Q","",min,max);     

      if(HSieveXDiff[idx]->GetFunction("gaus")){
	cout << "Fit worked!" << endl;
            
	TF1* xgaus1 = HSieveXDiff[idx]->GetFunction("gaus");
      
	min = xgaus1->GetParameter(1) - 3*xgaus1->GetParameter(2);
	max = xgaus1->GetParameter(1) + 3*xgaus1->GetParameter(2);
      
	HSieveXDiff[idx]->Fit("gaus","Q","",min,max);
	TF1* xgaus2 = HSieveXDiff[idx]->GetFunction("gaus");
                 
	x_mean = xgaus2->GetParameter(1);
	x_wid = xgaus2->GetParameter(2);                
	x_means[idx] = x_mean;
	x_widths[idx] = x_wid;

       	

      }


      max_bin = HSieveYDiff[idx]->GetMaximumBin();
      min = HSieveYDiff[idx]->GetBinCenter(max_bin) - width_guess;
      max = HSieveYDiff[idx]->GetBinCenter(max_bin) + width_guess;

      HSieveYDiff[idx]->Fit("gaus","Q","",min,max);

      
      if(HSieveYDiff[idx]->GetFunction("gaus")){
	cout << "Y fit worked" << endl;


	TF1* ygaus1 = HSieveYDiff[idx]->GetFunction("gaus");
      
	min = ygaus1->GetParameter(1) - 3*ygaus1->GetParameter(2);
	max = ygaus1->GetParameter(1) + 3*ygaus1->GetParameter(2);
      
	HSieveYDiff[idx]->Fit("gaus","Q","",min,max);
	TF1* ygaus2 = HSieveYDiff[idx]->GetFunction("gaus");
                 
	y_mean = ygaus2->GetParameter(1);
	y_wid = ygaus2->GetParameter(2);                
	y_means[idx] = y_mean;
	y_widths[idx] = y_wid;

	hist_mean = HSieveYDiff[idx]->GetMean();

	cout << "For Foil " << idx << ": " <<  endl;
	cout << "gaus_mean = " << y_mean << endl;
	cout << "gaus_width = " << y_wid << endl;
	cout << "hist_mean =  " << hist_mean << endl;
	mean_offset[idx] = TMath::Abs(y_mean-hist_mean)/y_wid;
	cout << "|(gaus_mean-hist_mean)|/gaus_width = " << mean_offset[idx] << endl;

	// condition to check if difference between gaus and histogram means (relative to gaus width) is large enough to indicate that there is backgroud in the plot
	if(mean_offset[idx] > 0.01){

	  if( (y_mean-hist_mean)< 0){
	    // hist mean more positive than gaus so linear background is positive wrt Gauss
	    lin_start = y_mean + 3*y_wid;
	    lin_end = HSieveYDiff[idx]->GetBinCenter(HSieveYDiff[idx]->FindLastBinAbove(2));


	  }
	  else if((y_mean-hist_mean)> 0){
	    // hist mean negative wrt gaus so linear background is negative wrt Gauss
	    lin_end = y_mean - 3*y_wid;
	    lin_start = HSieveYDiff[idx]->GetBinCenter(HSieveYDiff[idx]->FindFirstBinAbove(5));
	    
	  }


	  
	  // TF1* overall_fit = new TF1("overall_fit",overall,std::min(lin_start,y_mean-4*y_wid),std::max(lin_end,y_mean+4*y_wid),5);

	  // //	  TF1* overall_fit = new TF1("overall_fit",overall,-15,15,5);

	  // overall_fit->SetParameter(0,max);
	  // overall_fit->SetParameter(1,y_mean);
	  // overall_fit->SetParameter(2,y_wid);
	  // overall_fit->SetParameter(3,max*0.1);
	  // //	  overall_fit->SetParameter(4,0.5);

	  // HSieveYDiff[idx]->Fit(overall_fit,"R0+");
	  

	  // subtract gaus from histogram

	  HSieveYDiff_rem[idx] =  (TH1D*)HSieveYDiff[idx]->Clone("sieve clone");
	  
	  HSieveYDiff_rem[idx]->Add(ygaus2,-1);
	  HSieveYDiff_rem[idx]->SetLineColor(kMagenta);
	  
	  
	  Double_t fba = HSieveYDiff_rem[idx]->FindFirstBinAbove(max*0.1,y_mean+y_wid);
	  
	  cout << "Found first bin = " << HSieveYDiff_rem[idx]->GetBinCenter(fba) << endl;
	  
	  //	  y_bg[idx] = new TF1(Form("y_bg_%i",idx),bg,HSieveYDiff_rem[idx]->GetBinCenter(fba),lin_end,2);
	  //	  y_bg[idx] = new TF1(Form("y_bg_%i",idx),bg,y_mean+3*y_wid,lin_end,2);
	  y_bg[idx] = new TF1(Form("y_bg_%i",idx),bg,lin_start,lin_end,2);

	  y_bg[idx]->SetLineColor(kBlack);
	  HSieveYDiff_rem[idx]->Fit(Form("y_bg_%i",idx),"R");
	  

	  // fill in cut offs for tail or Gaus

	  gaus_start[idx] = y_mean - 3*y_wid;
	  gaus_end[idx] = y_mean + 3*y_wid;

	  lin_starts[idx] = lin_start;
	  lin_ends[idx] = lin_end;


	  cout << "For foil " << idx << "gaus start -> end = " << gaus_start[idx] << " -> " << gaus_end[idx] << endl;
	  
	  y_bg[idx] = new TF1(Form("y_bg_%i",idx),bg,std::min(lin_start,y_mean-4*y_wid),std::max(lin_end,y_mean+4*y_wid),2);

	  // y_bg[idx]->SetParameter(0,overall_fit->GetParameter(3)); // set heigh of linear bg to 1/10th of Gauss
	  // y_bg[idx]->SetParameter(1,overall_fit->GetParameter(4));

	  // y_bg[idx]->SetLineColor(kBlue);
	  // y_bg[idx]->SetLineStyle(9);

	  // cout << "background parameters: " << endl;
	  // cout << "y = " << y_bg[idx]->GetParameter(1) << "x + " << y_bg[idx]->GetParameter(0) << endl << endl;

	  y_gaus[idx] = new TF1(*ygaus2);

	  
	  // //	  y_gaus[idx] = new TF1(Form("y_gaus_%i",idx),peak,y_mean-3*y_wid,y_mean+3*y_wid,3);
	  // y_gaus[idx] = new TF1(Form("y_gaus_%i",idx),peak,std::min(lin_start,y_mean-4*y_wid),std::max(lin_end,y_mean+4*y_wid),3);

	  // y_gaus[idx]->SetParameter(0,overall_fit->GetParameter(0));
	  // y_gaus[idx]->SetParameter(1,overall_fit->GetParameter(1));
	  // y_gaus[idx]->SetParameter(2,overall_fit->GetParameter(2));

	  // y_gaus[idx]->SetLineColor(kGreen+4);
	  // y_gaus[idx]->SetLineStyle(9);
	  
	  // cout << "Gauss parameters: " << endl;
	  // cout << "mean = " << y_gaus[idx]->GetParameter(1) << ", width = " << y_gaus[idx]->GetParameter(2) << ", peak = " << y_gaus[idx]->GetParameter(0) << endl;
	  
	}
       
	
      }

    }


    // fill histograms for main gaus and tail of y_sieve distrib
    // don't fill histos if no substantial deviation from gaus



    // establish limits of histograms


    Double_t xfp_h = 0.8;
    Double_t xfp_l = -0.8;

    Double_t yfp_h = 0.05;
    Double_t yfp_l = -0.05;

    Double_t thfp_h = 30;
    Double_t thfp_l = -35;

    Double_t phfp_h = 50;
    Double_t phfp_l = -50;

    
    for (UInt_t idx = 0; idx < NFoils; idx++) {

      H_FP_X_tail[idx] = new TH1D(Form("FP X Tail Foil %d", idx), "FP X tail for " + Foil_names[idx], 200, xfp_l,xfp_h);
      H_FP_X_tail[idx]->GetXaxis()->SetTitle("FP x [m]");
      H_FP_Y_tail[idx] = new TH1D(Form("FP Y Tail Foil %d", idx), "FP Y tail for " + Foil_names[idx], 200, yfp_l,yfp_h);
      H_FP_Y_tail[idx]->GetXaxis()->SetTitle("FP y [m]");
      H_FP_P_tail[idx] = new TH1D(Form("FP P Tail Foil %d", idx), "FP P tail for " + Foil_names[idx], 200, phfp_l,phfp_h);
      H_FP_P_tail[idx]->GetXaxis()->SetTitle("FP #phi [mrad]");
      H_FP_T_tail[idx] = new TH1D(Form("FP T Tail Foil %d", idx), "FP T tail for " + Foil_names[idx], 200, thfp_l,thfp_h);
      H_FP_T_tail[idx]->GetXaxis()->SetTitle("FP #theta [mrad]");
      H_FP_TP_tail[idx] = new TH2D(Form("FP T vs P Tail Foil %d", idx), "FP T vs P tail for " + Foil_names[idx], 200, phfp_l,phfp_h, 200, thfp_l,thfp_h);
      H_FP_TP_tail[idx]->GetXaxis()->SetTitle("FP #phi [mrad]");
      H_FP_TP_tail[idx]->GetYaxis()->SetTitle("FP #theta [mrad]");
      H_FP_YP_tail[idx] = new TH2D(Form("FP Y vs P Tail Foil %d", idx), "FP Y vs P tail for " + Foil_names[idx], 200, phfp_l,phfp_h, 200, yfp_l,yfp_h);
      H_FP_YP_tail[idx]->GetXaxis()->SetTitle("FP #phi [mrad]");
      H_FP_YP_tail[idx]->GetYaxis()->SetTitle("FP y [m]");

      H_beam_xy_tail[idx] = new TH2D(Form("Beam Tail Foil %d", idx), "Beam tail for " + Foil_names[idx], 200, -0.005,0.005, 200, -0.005,0.005);
      H_beam_xy_tail[idx]->GetXaxis()->SetTitle("Beam x [m]");
      H_beam_xy_tail[idx]->GetYaxis()->SetTitle("Beam y [m]");



      
      H_FP_X_gaus[idx] = new TH1D(Form("FP X Gaus Foil %d", idx), "FP X gaus for " + Foil_names[idx], 200, xfp_l,xfp_h);
      H_FP_X_gaus[idx]->GetXaxis()->SetTitle("FP x [m]");
      H_FP_Y_gaus[idx] = new TH1D(Form("FP Y Gaus Foil %d", idx), "FP Y gaus for " + Foil_names[idx], 200, yfp_l,yfp_h);
      H_FP_Y_gaus[idx]->GetXaxis()->SetTitle("FP y [m]");
      H_FP_P_gaus[idx] = new TH1D(Form("FP P Gaus Foil %d", idx), "FP P gaus for " + Foil_names[idx], 200, phfp_l,phfp_h);
      H_FP_P_gaus[idx]->GetXaxis()->SetTitle("FP #phi [mrad]");
      H_FP_T_gaus[idx] = new TH1D(Form("FP T Gaus Foil %d", idx), "FP T gaus for " + Foil_names[idx], 200, thfp_l,thfp_h);
      H_FP_T_gaus[idx]->GetXaxis()->SetTitle("FP #theta [mrad]");
      H_FP_TP_gaus[idx] = new TH2D(Form("FP T vs P Gaus Foil %d", idx), "FP T vs P gaus for " + Foil_names[idx], 200, phfp_l,phfp_h, 200, thfp_l,thfp_h);
      H_FP_TP_gaus[idx]->GetXaxis()->SetTitle("FP #phi [mrad]");
      H_FP_TP_gaus[idx]->GetYaxis()->SetTitle("FP #theta [mrad]");
      H_FP_YP_gaus[idx] = new TH2D(Form("FP Y vs P Gaus Foil %d", idx), "FP Y vs P gaus for " + Foil_names[idx], 200, phfp_l,phfp_h, 200, yfp_l,yfp_h);
      H_FP_YP_gaus[idx]->GetXaxis()->SetTitle("FP #phi [mrad]");
      H_FP_YP_gaus[idx]->GetYaxis()->SetTitle("FP y [m]");

      H_beam_xy_gaus[idx] = new TH2D(Form("Beam Gaus Foil %d", idx), "Beam gaus for " + Foil_names[idx], 200, -0.005,0.005, 200, -0.005,0.005);
      H_beam_xy_gaus[idx]->GetXaxis()->SetTitle("Beam x [m]");
      H_beam_xy_gaus[idx]->GetYaxis()->SetTitle("Beam y [m]");
      
      
    }



    UInt_t check_gaus[NFoils] = {0};
    UInt_t check_tail[NFoils] = {0};
    
    for (UInt_t idx = 0; idx < fNRawData; idx++) {
      const EventData &eventdata = fRawData[idx];
      
      UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID]; //starting 0!


      const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
      const UInt_t Row = (UInt_t)eventdata.Data[kRowID]; //starting 0!
      const TVector3 SieveHoleCorrectionTCS = GetSieveHoleCorrectionTCS(FoilID, Col, Row);
      Double_t ProjectionY = eventdata.Data[kRealTgY] + eventdata.Data[kCalcPh] * (SieveHoleCorrectionTCS.Z());

      Double_t fp_x = eventdata.Data[kX];
      Double_t fp_y = eventdata.Data[kY];
      Double_t fp_th = eventdata.Data[kTh];
      Double_t fp_phi = eventdata.Data[kPhi];
      Double_t beam_x = eventdata.Data[kBeamX];
      Double_t beam_y = eventdata.Data[kBeamY];

      
      PlotID = FoilID; // useful for case when only Vertical foils are used
                    // (FoilID based on scheme with optics and vertical foils where optics foils
                    // are 0-7 and vertical foils 8-10) so they can have PlotID 0,1,2
      if(NFoils<11){

	PlotID -= 8;
      }
      

      dY[PlotID] = 1000*(ProjectionY - SieveHoleCorrectionTCS.Y());



      // gaus or tail?
      if(dY[PlotID]>gaus_start[PlotID] && dY[PlotID]< gaus_end[PlotID]){
	H_FP_X_gaus[PlotID]->Fill(fp_x);
	H_FP_Y_gaus[PlotID]->Fill(fp_y);
	H_FP_P_gaus[PlotID]->Fill(1000*fp_phi);
	H_FP_T_gaus[PlotID]->Fill(1000*fp_th);
	H_FP_TP_gaus[PlotID]->Fill(1000*fp_phi,1000*fp_th);
	H_FP_YP_gaus[PlotID]->Fill(1000*fp_phi,fp_y);
	H_beam_xy_gaus[PlotID]->Fill(beam_x,beam_y);
	

	check_gaus[PlotID]++;
      }
      else if(dY[PlotID]>lin_starts[PlotID] && dY[PlotID]< lin_ends[PlotID]){
	H_FP_X_tail[PlotID]->Fill(fp_x);
	H_FP_Y_tail[PlotID]->Fill(fp_y);
	H_FP_P_tail[PlotID]->Fill(1000*fp_phi);
	H_FP_T_tail[PlotID]->Fill(1000*fp_th);
	H_FP_TP_tail[PlotID]->Fill(1000*fp_phi,1000*fp_th);
	H_FP_YP_tail[PlotID]->Fill(1000*fp_phi,fp_y);
	H_beam_xy_tail[PlotID]->Fill(beam_x,beam_y);

	
	check_tail[PlotID]++;
      }
      
      
    }
    
    cout << endl;
    //    cout << "ratio of gaus to tail = " << check_gaus/check_tail << endl;
    cout << "V1: % of tail (" << check_tail[0] <<") to overall (" << HSieveYDiff[0]->GetEntries() << ") = " << 100.*check_tail[0]/(Double_t)HSieveYDiff[0]->GetEntries() << endl;
    cout << "V2: % of tail (" << check_tail[1] <<") to overall (" << HSieveYDiff[1]->GetEntries() << ") = " << 100.*check_tail[1]/(Double_t)HSieveYDiff[1]->GetEntries() << endl;


    

    cout << endl;

    // create canvas for x sieve difference
    TCanvas * c1;
    
    if (nplot <= 1) {
      c1 = new TCanvas("XDiff", "XDiff", 1000, 1000);
      c1->Divide(1, 1);
    } else if (nplot <= 3) {
      c1 = new TCanvas("XDiff", "XDiff", 1800, 1100);
      c1->Divide(3, 1);
    } else if (nplot <= 6) {
      c1 = new TCanvas("XDiff", "XDiff", 1800, 1100);
      c1->Divide(3, 2);
    } else if (nplot <= 9){
      c1 = new TCanvas("XDiff", "XDiff", 1800, 1100);
      c1->Divide(3, 3);
    }
    else{
      c1 = new TCanvas("XDiff", "XDiff", 1800, 1100);
      c1->Divide(4, 3);
    }


    // create canvas for y sieve difference
    TCanvas * c2;
    
    if (nplot <= 1) {
      c2 = new TCanvas("YDiff", "YDiff", 1000, 1000);
      c2->Divide(1, 1);
    } else if (nplot <= 3) {
      c2 = new TCanvas("YDiff", "YDiff", 1800, 1100);
      c2->Divide(3, 1);
    } else if (nplot <= 6) {
      c2 = new TCanvas("YDiff", "YDiff", 1800, 1100);
      c2->Divide(3, 2);
    } else if (nplot <= 9){
      c2 = new TCanvas("YDiff", "YDiff", 1800, 1100);
      c2->Divide(3, 3);
    }
    else{
      c2 = new TCanvas("YDiff", "YDiff", 1800, 1100);
      c2->Divide(4, 3);
    }


    // create canvases for FP distribs
    TCanvas * c3; //x_fp
    TCanvas * c4; //y_fp
    TCanvas * c5; //th_fp
    TCanvas * c6; //ph_fp
    TCanvas * c7; //th_fp vs ph_fp
    TCanvas * c8; //y_fp vs ph_fp
    TCanvas * c9; // beam y vs x
		
    gStyle->SetOptStat(1);
    gPad->Modified();
    
    c3 = new TCanvas("FP x","FP x",1800,1100);
    c3->Divide(3,2);

    c4 = new TCanvas("FP y","FP y",1800,1100);
    c4->Divide(3,2);

    c5 = new TCanvas("FP th","FP th",1800,1100);
    c5->Divide(3,2);

    c6 = new TCanvas("FP ph","FP ph",1800,1100);
    c6->Divide(3,2);

    c7 = new TCanvas("FP th vs ph","FP th vs ph",1800,1100);
    c7->Divide(3,2);

    c8 = new TCanvas("FP y vs ph","FP y vs ph",1800,1100);
    c8->Divide(3,2);

    c9 = new TCanvas("Beam y vs x","Beam y vs x",1800,1100);
    c9->Divide(3,2);

    for (UInt_t idx = 0; idx < nplot; idx++) {

      c3->cd(idx+1);
      
      H_FP_X_gaus[idx]->Draw();
      H_FP_X_gaus[idx]->SetStats(1);
      c3->Update();
      
      c3->cd(idx+4);
      H_FP_X_tail[idx]->Draw();
      H_FP_X_tail[idx]->SetStats(1);
      c3->Update();


      c4->cd(idx+1);
      H_FP_Y_gaus[idx]->Draw();
      H_FP_Y_gaus[idx]->SetStats(1);
      c4->cd(idx+4);
      H_FP_Y_tail[idx]->Draw();
      H_FP_Y_tail[idx]->SetStats(1);

      c5->cd(idx+1);
      H_FP_T_gaus[idx]->Draw();
      H_FP_T_gaus[idx]->SetStats(1);
      c5->cd(idx+4);
      H_FP_T_tail[idx]->Draw();
      H_FP_T_tail[idx]->SetStats(1);

      c6->cd(idx+1);
      H_FP_P_gaus[idx]->Draw();
      H_FP_P_gaus[idx]->SetStats(1);
      c6->cd(idx+4);
      H_FP_P_tail[idx]->Draw();
      H_FP_P_tail[idx]->SetStats(1);

      c7->cd(idx+1);
      H_FP_TP_gaus[idx]->Draw("colz");
      H_FP_TP_gaus[idx]->SetStats(1);
      c7->cd(idx+4);
      H_FP_TP_tail[idx]->Draw("colz");
      H_FP_TP_tail[idx]->SetStats(1);

      c8->cd(idx+1);
      H_FP_YP_gaus[idx]->Draw("colz");
      H_FP_YP_gaus[idx]->SetStats(1);
      c8->cd(idx+4);
      H_FP_YP_tail[idx]->Draw("colz");
      H_FP_YP_tail[idx]->SetStats(1);

      c9->cd(idx+1);
      H_beam_xy_gaus[idx]->Draw("colz");
      H_beam_xy_gaus[idx]->SetStats(1);
      c9->cd(idx+4);
      H_beam_xy_tail[idx]->Draw("colz");
      H_beam_xy_tail[idx]->SetStats(1);
      
    }
    
    cout << "completed c3" << endl;

      

    for (UInt_t idx = 0; idx < nplot; idx++) {


      c1->cd(idx+1);
     
      HSieveXDiff[idx]->Draw();

      // TPaveText *ttx = new TPaveText(0.6,0.7,0.88,0.82,"NDC");
      // ttx->AddText(Form("\\Delta = %2.3f mm",  x_means[idx]));
      // ttx->AddText(Form("\\sigma = %2.3f mm", x_widths[idx]));
      // ttx->SetShadowColor(0);
      // ttx->SetFillColor(0);
      // ttx->Draw("same");


      TPaveText *ttx_1 = new TPaveText(0.6,0.75,0.85,0.82,"NDC");
      ttx_1->AddText("\\sigma = ");
      ttx_1->AddText(Form("%2.3f mm", x_widths[idx]));


      ttx_1->SetShadowColor(0);
      ttx_1->SetFillColor(0);
      ttx_1->SetTextSize(0.08);
      ttx_1->Draw("same");


      TPaveText *ttx_2 = new TPaveText(0.25,0.75,0.28,0.82,"NDC");

      ttx_2->AddText("\\Delta = ");
      ttx_2->AddText(Form("%2.3f mm",  x_means[idx]));

      ttx_2->SetShadowColor(0);
      ttx_2->SetFillColor(0);
      ttx_2->SetTextSize(0.08);
      ttx_2->Draw("same");


      c2->cd(idx+1);
     
      HSieveYDiff[idx]->Draw();

      
      
      y_bg[idx]->Draw("same");
      y_gaus[idx]->Draw("same");

      cout << "Y diff bg drawing for foil " << idx+1 << " : pad name is " << gPad->GetName() << ", title is " << gPad->GetTitle() << endl;
      
      HSieveYDiff_rem[idx]->Draw("same");

      TPaveText *tty_1 = new TPaveText(0.6,0.75,0.85,0.82,"NDC");
      //      tty_1->AddText(Form("\\Delta = %2.3f mm",  y_means[idx]));
      tty_1->AddText("\\sigma = ");
      tty_1->AddText(Form("%2.3f mm", y_widths[idx]));


      tty_1->SetShadowColor(0);
      tty_1->SetFillColor(0);
      tty_1->SetTextSize(0.08);
      tty_1->Draw("same");


      TPaveText *tty_2 = new TPaveText(0.25,0.75,0.28,0.82,"NDC");

      //      tty_2->AddText(Form("\\Delta = %2.3f mm",  y_means[idx]));
      tty_2->AddText("\\Delta = ");
      tty_2->AddText(Form("%2.3f mm",  y_means[idx]));

      
      //tty_2->AddText(Form("\\sigma = %2.3f mm", y_widths[idx]));
      tty_2->SetShadowColor(0);
      tty_2->SetFillColor(0);
      tty_2->SetTextSize(0.08);
      tty_2->Draw("same");

      cout << "loops or c1 and c2: " << idx << endl;
      
    }       
      
    cout << "c1 returning..." << endl;
    return c1;
    cout << "c1 returned" << endl;


}





/*
void LOpticsOpt::check_fit_qual_Th()
{
    // return square sum of diff between calculated tg_th and expected tg_th

    TH1D *h_theta_data = new TH1D("Vertical angle (target)", "Vertical angle (target)", 200, -0.2, 0.2);
    TH1D *h_theta_fit = new TH1D("htheta_fit", "htheta_fit", 200, -0.2, 0.2);
    TH1D *h_dtheta = new TH1D("dtheta", "Vertical angle difference (fit-data)", 2000, -0.05, 0.05);
    h_theta_data->SetLineColor(2);
    h_theta_fit->SetLineColor(4);
    h_theta_data->GetXaxis()->SetTitle("vert angle [rad]");
    h_dtheta->GetXaxis()->SetTitle("difference [rad]");
    h_theta_data->GetYaxis()->SetTitle("Number of events");
    h_theta_data->GetYaxis()->SetTitleOffset(1.35);

    TH2D *h_phi_theta = new TH2D("Horizontal vs vertical angle (target)", "Horizontal vs vertical angle (target)", 200, -0.07, 0.07, 200, -0.07, 0.07);

    TH1D *h_phi_data = new TH1D("Horizontal angle (target)", "Horizontal angle (target)", 500, -0.035, 0.035);
    TH1D *h_phi_fit = new TH1D("hphi_fit", "hphi_fit", 500, -0.035,0.035);
    TH1D *h_dphi = new TH1D("dhphi", "Horizontal angle difference (fit-data)", 2000, -0.02, 0.02);
    h_phi_data->SetLineColor(2);
    h_phi_fit->SetLineColor(4);

    h_phi_data->GetXaxis()->SetTitle("horiz angle [rad]");
    h_dphi->GetXaxis()->SetTitle("difference [rad]");
    h_phi_data->GetYaxis()->SetTitle("Number of events");
    h_phi_data->GetYaxis()->SetTitleOffset(1.35);
    h_phi_data->SetLineColor(2);
    h_phi_fit->SetLineColor(4);
    TLegend *leg = new TLegend(0.7, 0.8, 0.99, 0.95);
    leg->AddEntry(h_phi_data, "Column Positions", "l");
    leg->AddEntry(h_phi_fit, "fit", "l");

    TCanvas *c_th = new TCanvas("theta", "theta", 1600, 900);
    TCanvas *c_ph = new TCanvas("phi", "phi", 1600, 900);
    TCanvas *c_th_ph = new TCanvas("thetaphi", "thetaphi", 1600, 900);
    c_th->Divide(2,1);
    c_ph->Divide(2,1);
    Double_t dth = 0; //Difference
    Double_t dph = 0; //Difference
    Double_t rmsth = 0; //mean square

    static UInt_t NCall = 0;
    NCall++;

    Double_t theta;
    Double_t phi;

    for (UInt_t idx = 0; idx < fNRawData; idx++) 
    {
        EventData &eventdata = fRawData[idx];

        Double_t x_fp = eventdata.Data[kX];
        const Double_t(*powers)[5] = eventdata.powers;

        CalcMatrix(x_fp, fTMatrixElems);
        theta = CalcTargetVar(fTMatrixElems, powers);
//        if (fabs(theta)<1.e+10)
        dth = theta - eventdata.Data[ktarXp];
        h_theta_data->Fill(eventdata.Data[ktarXp]);
        h_theta_fit->Fill(theta);
        h_dtheta->Fill(dth);

        CalcMatrix(x_fp, fPMatrixElems);
        phi = CalcTargetVar(fPMatrixElems, powers);
        dph = phi - eventdata.Data[ktarYp];
        h_phi_data->Fill(eventdata.Data[ktarYp]);
        h_phi_fit->Fill(phi);
        h_dphi->Fill(dph);
        h_phi_theta->Fill(phi,theta);
    }
    c_th->cd(1);
    h_theta_data->Draw();
    h_theta_fit->Draw("Same");
    leg->Draw();

    c_th->cd(2);
    h_dtheta->Draw();

    c_ph->cd(1);
    h_phi_data->Draw();
//    h_theta_data->Draw("Same");
    h_phi_fit->Draw("Same");
    leg->Draw();

    c_ph->cd(2);
    h_dphi->Draw();

    c_th_ph->cd();
    h_phi_theta->Draw();
}

void LOpticsOpt::check_fit_qual_y()
{
    // return square sum of diff between calculated tg_y and expected tg_y

    Double_t dtg_y = 0; //Difference
    Double_t dtg_y_rms = 0; //mean square

    static UInt_t NCall = 0;
    NCall++;

    Double_t tg_y;
    gStyle->SetOptStat(1);

    TH1D *h_y_data = new TH1D("Z coordinate distribution", "Z coordinate distribution", 330, -0.66, 0.66);
    TH1D *h_y_fit = new TH1D("Z vertex", "Z vertex", 330, -0.66, 0.66);
    TH1D *h_dy = new TH1D("Difference (fit-data)", "Difference (fit-data)", 330, -0.66, 0.66);
    h_dy->GetXaxis()->SetTitle("#Deltaz coord (m)");
    h_y_fit->GetXaxis()->SetTitle("z coord (m)");
    h_y_fit->GetYaxis()->SetTitle("Number of events");
    h_y_fit->GetYaxis()->SetTitleOffset(1.27);
    h_y_data->SetLineColor(2);
    h_y_fit->SetLineColor(4);
    TLegend *leg = new TLegend (0.7, 0.8, 0.99, 0.99);
    leg->AddEntry(h_y_data, "Wire position", "l");
    leg->AddEntry(h_y_fit, "Fit", "l");
  

    gStyle->SetStatY(0.8);
    gStyle->SetStatX(0.9);
    gStyle->SetStatW(0.15);
    gStyle->SetStatH(0.25); 




    cout<<"lalallalaa --------------------------------------->"<<fNRawData<<endl;
    for (UInt_t idx = 0; idx < fNRawData; idx++) {
        EventData &eventdata = fRawData[idx];

        Double_t x_fp = eventdata.Data[kX];
        const Double_t(*powers)[5] = eventdata.powers;

        // calculate the matrices we need
        // CalcMatrix(x_fp, fDMatrixElems);
        // CalcMatrix(x_fp, fTMatrixElems);
        CalcMatrix(x_fp, fYMatrixElems);
        CalcMatrix(x_fp, fYTAMatrixElems);
        // CalcMatrix(x_fp, fPMatrixElems);
        // CalcMatrix(x_fp, fPTAMatrixElems);

        // calculate the coordinates at the target
        tg_y = CalcTargetVar(fYMatrixElems, powers) + CalcTargetVar(fYTAMatrixElems, powers);
//        cout<<idx<<"  opt "<<tg_y<<endl;

        if (fabs(tg_y)<1.e+10)
        {
//        cout<<idx<<"           "<<CalcTargetVar(fYMatrixElems, powers)<<",           "<<CalcTargetVar(fYTAMatrixElems, powers)<<endl;

//        int aaa=system("sleep 10m");


        // UInt_t res = (UInt_t) eventdata.Data[kCutID];
        // res = res % (NSieveRow * NSieveCol * NFoils);
        // const UInt_t FoilID = res / (NSieveRow * NSieveCol); //starting 0!

	const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID]; //starting 0!
	// res = res%(NSieveRow*NSieveCol);
	const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
	const UInt_t Row = (UInt_t)eventdata.Data[kRowID]; //starting 0!


        const Double_t ArbitaryVertexShift = fArbitaryVertexShift[FoilID] * (-TMath::Sin(HRSAngle));

        dtg_y = tg_y - eventdata.Data[ktarY] + 0.*ArbitaryVertexShift;
        dtg_y_rms += (tg_y - eventdata.Data[ktarY] + 0.*ArbitaryVertexShift)*(tg_y - eventdata.Data[ktarY] + 0.*ArbitaryVertexShift);

        h_y_data->Fill(eventdata.Data[ktarY]/sin(5.*3.1415927/180.));
        h_y_fit->Fill(tg_y/sin(5.*3.1415927/180.));
        h_dy->Fill(dtg_y/sin(5.*3.1415927/180.));

	*/
/*
    {
      // calculates target coordinates from focal plane coordinates

      const Int_t kNUM_PRECOMP_POW = 10;

      Double_t x_fp, y_fp, th_fp, ph_fp;
      Double_t powers[kNUM_PRECOMP_POW][5];  // {(x), th, y, ph, abs(th) }
      Double_t x, y, theta, phi, dp, p, pathl;

      // first select the coords to use
      {
        x_fp = eventdata.Data[kX];
        y_fp = eventdata.Data[kY];
        th_fp = eventdata.Data[kTh];
        ph_fp = eventdata.Data[kPhi];
      }

      // calculate the powers we need
      for(int i=0; i<kNUM_PRECOMP_POW; i++) 
      {
        powers[i][0] = pow(x_fp, i);
        powers[i][1] = pow(th_fp, i);
        powers[i][2] = pow(y_fp, i);
        powers[i][3] = pow(ph_fp, i);
        powers[i][4] = pow(TMath::Abs(th_fp),i);
      }

      // calculate the matrices we need
      CalcMatrix(x_fp, fYMatrixElems);
      CalcMatrix(x_fp, fYTAMatrixElems);

      // calculate the coordinates at the target
      y = CalcTargetVar(fYMatrixElems, powers)+CalcTargetVar(fYTAMatrixElems,powers);
      cout<<"official "<<y<<endl<<endl;
    }
    /*


	
//	if (fabs(tg_y)<1.e+10) cout<<"dtg_y="<<dtg_y<<"      tg_y="<<tg_y<<"      eventdata.Data["<<kRealTgY<<"]="<<eventdata.Data[kRealTgY]<<"           ArbitaryVertexShift="<<ArbitaryVertexShift<<endl;
//	if (fabs(tg_y)>1.e+10) cout<<"         dtg_y="<<dtg_y<<"      tg_y="<<tg_y<<"      eventdata.Data["<<kRealTgY<<"]="<<eventdata.Data[kRealTgY]<<"           ArbitaryVertexShift="<<ArbitaryVertexShift<<endl;
      

        // save the results
        eventdata.Data[kCalcTgY] = tg_y;
        }

    }
    TCanvas *c_y = new TCanvas("y", "y", 1600, 900);
    c_y->Divide(2,1);

    c_y->cd(1);
    h_y_fit->Draw();
//    h_theta_data->Draw("Same");
    h_y_data->Draw("Same");
    leg->Draw();

    c_y->cd(2);
    h_dy->Draw();
    c_y->SaveAs("vz_calib.gif");
    //    cout<<"fNRawData:"<<fNRawData<<endl;

    //DEBUG_INFO("SumSquareDTgY", "#%d : dtg_y = %f, dtg_y_rms = %f", NCall, dtg_y / fNRawData, TMath::Sqrt(dtg_y_rms / fNRawData));
//    printf("SumSquareDTgY: #%d : dtg_y = %f, dtg_y_rms = %f\n", NCall, dtg_y / fNRawData, TMath::Sqrt(dtg_y_rms / fNRawData));
    printf("SumSquareDTgY: #%d : dtg_y = %f, dtg_y_rms = %f\n, %d rows", NCall, dtg_y / fNRawData, (dtg_y_rms / fNRawData), fNRawData);
//    cout<<fNRawData<<endl;
//    int aa=system("sleep 1m");

}
*/


Double_t LOpticsOpt::SumSquareDTh()
{
    // return square sum of diff between calculated tg_th and expected tg_th

    Double_t dth = 0; //Difference
    Double_t rmsth = 0; //mean square
    Double_t dx = 0;
    Double_t rmsdx = 0;
    
    static UInt_t NCall = 0;
    NCall++;

    Double_t theta;


    for (UInt_t idx = 0; idx < fNRawData; idx++) {
        EventData &eventdata = fRawData[idx];


	// UInt_t res = (UInt_t) eventdata.Data[kCutID];
        // res = res % (NSieveRow * NSieveCol * NFoils);
        // const UInt_t FoilID = res / (NSieveRow * NSieveCol); //starting 0!
        // res = res % (NSieveRow * NSieveCol);
        // const UInt_t Col = res / (NSieveRow); //starting 0!
        // const UInt_t Row = res % (NSieveRow); //starting 0!


      // condition to remove events from outside main gaussian of phi_tg distribution (from initial matrix)
      if(eventdata.Data[kIgnore] == 0){

	continue;
      }
	

	
	const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID]; //starting 0!
	// res = res%(NSieveRow*NSieveCol);
	const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
	const UInt_t Row = (UInt_t)eventdata.Data[kRowID]; //starting 0!
	
	Double_t x_fp = eventdata.Data[kX];
        const Double_t(*powers)[5] = eventdata.powers;

        // calculate the matrices we need
        // CalcMatrix(x_fp, fDMatrixElems);
        CalcMatrix(x_fp, fTMatrixElems);
        // CalcMatrix(x_fp, fYMatrixElems);
        // CalcMatrix(x_fp, fYTAMatrixElems);
        // CalcMatrix(x_fp, fPMatrixElems);
        // CalcMatrix(x_fp, fPTAMatrixElems);

        // calculate the coordinates at the target
        theta = CalcTargetVar(fTMatrixElems, powers);

	
	
        dth += theta - eventdata.Data[kRealThMatrix];
        rmsth += (theta - eventdata.Data[kRealThMatrix])*(theta - eventdata.Data[kRealThMatrix]);

        DEBUG_MASSINFO("SumSquareDTh", "D_Th = %f = \t%f - \t%f", theta - eventdata.Data[kRealThMatrix], theta, eventdata.Data[kRealThMatrix]);

        // save the results
        eventdata.Data[kCalcTh] = theta;

	const TVector3 SieveHoleCorrectionTCS = GetSieveHoleCorrectionTCS(FoilID, Col, Row);
	//	const TVector3 BeamSpotHCS(eventdata.Data[kBeamX], eventdata.Data[kBeamY], targetfoils[FoilID]);
	const TVector3 BeamSpotHCS( eventdata.Data[kBeamX] + (targetfoils[FoilID]/BeamZDir_average)*eventdata.Data[kBeamDirX], eventdata.Data[kBeamY] + (targetfoils[FoilID]/BeamZDir_average)*eventdata.Data[kBeamDirY], targetfoils[FoilID]);

	const TVector3 BeamSpotTCS = fTCSInHCS.Inverse()*(BeamSpotHCS - fPointingOffset);
	const Double_t x_tg = BeamSpotTCS.X() - BeamSpotTCS.Z() * tan(eventdata.Data[kCalcTh]);
	//	double ProjectionX = x_tg + (tan(eventdata.Data[kCalcTh]) + x_tg * ExtTarCor_ThetaCorr) * (ZPos);
	double ProjectionX = x_tg + (tan(eventdata.Data[kCalcTh])*(ZPos));
	
	const Double_t posx = SieveHoleCorrectionTCS.X();
	
	dx += ProjectionX - posx;
	rmsdx += (ProjectionX - posx)*(ProjectionX - posx);
    }
    

    DEBUG_INFO("SumSquareDTh", "#%d : dth = %f,\t rmsth = %f", NCall, dth / fNRawData, TMath::Sqrt(rmsth / fNRawData));
    printf("SumSquareDTh: #%d : dth = %f,\t rmsth = %f\t, orig rmsth = %f\n", NCall, dth / fNRawData, TMath::Sqrt(rmsth / fNRawData), rmsth);
    printf("SumSquareDX: #%d : dx = %f,\t rmsdx = %f\t, orig rmsdx = %f\n", NCall, dx / fNRawData, TMath::Sqrt(rmsdx / fNRawData), rmsdx);

    //return rmsth; //Use for opt to angles
    return rmsdx; //Use for opt to sieve plane
}

 std::pair<Double_t,Double_t> LOpticsOpt::SumSquareDPhi()
{
    // return square sum of diff between calculated tg_ph and expected tg_ph

    Double_t dphi = 0; //Difference
    Double_t rmsphi = 0; //mean square
    Double_t dy = 0;
    Double_t rmsdy = 0;

    static UInt_t NCall = 0;
    NCall++;

    Double_t phi;

    //ignore-related parameters
    Int_t Nignored = 0;
    
    for (UInt_t idx = 0; idx < fNRawData; idx++) {
      

      
      EventData &eventdata = fRawData[idx];


      // condition to remove events from outside main gaussian (from initial matrix)
      if(eventdata.Data[kIgnore] == 0){
	Nignored++;

	continue;
      }
	

	// UInt_t res = (UInt_t) eventdata.Data[kCutID];
        // res = res % (NSieveRow * NSieveCol * NFoils);
        // const UInt_t FoilID = res / (NSieveRow * NSieveCol); //starting 0!
        // res = res % (NSieveRow * NSieveCol);
        // const UInt_t Col = res / (NSieveRow); //starting 0!
        // const UInt_t Row = res % (NSieveRow); //starting 0!

	const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID]; //starting 0!
	// res = res%(NSieveRow*NSieveCol);
	const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
	const UInt_t Row = (UInt_t)eventdata.Data[kRowID]; //starting 0!

	
        Double_t x_fp = eventdata.Data[kX];
        const Double_t(*powers)[5] = eventdata.powers;

        // calculate the matrices we need
        // CalcMatrix(x_fp, fDMatrixElems);
        // CalcMatrix(x_fp, fTMatrixElems);
        // CalcMatrix(x_fp, fYMatrixElems);
        // CalcMatrix(x_fp, fYTAMatrixElems);
        CalcMatrix(x_fp, fPMatrixElems);
        CalcMatrix(x_fp, fPTAMatrixElems);

        // calculate the coordinates at the target
	//        phi = CalcTargetVar(fPMatrixElems, powers) + CalcTargetVar(fPTAMatrixElems, powers);
	phi = CalcTargetVar(fPMatrixElems, powers);
       
        dphi += phi - eventdata.Data[kRealPhi];
        rmsphi += (phi - eventdata.Data[kRealPhi])*(phi - eventdata.Data[kRealPhi]);

        DEBUG_MASSINFO("SumSquareDPhi", "D_Phi = %f = \t%f - \t%f", phi - eventdata.Data[kRealPhi], phi, eventdata.Data[kRealPhi]);

        //save the results
        eventdata.Data[kCalcPh] = phi;

	const TVector3 SieveHoleCorrectionTCS = GetSieveHoleCorrectionTCS(FoilID, Col, Row);
	//	const TVector3 BeamSpotHCS(eventdata.Data[kBeamX], eventdata.Data[kBeamY], targetfoils[FoilID]);
	const TVector3 BeamSpotHCS( eventdata.Data[kBeamX] + (targetfoils[FoilID]/BeamZDir_average)*eventdata.Data[kBeamDirX], eventdata.Data[kBeamY] + (targetfoils[FoilID]/BeamZDir_average)*eventdata.Data[kBeamDirY], targetfoils[FoilID]);

	const TVector3 BeamSpotTCS = fTCSInHCS.Inverse()*(BeamSpotHCS - fPointingOffset);
	const Double_t y_tg = BeamSpotTCS.Y() - BeamSpotTCS.Z() * tan(eventdata.Data[kCalcPh]);

	double ProjectionY = y_tg + tan(eventdata.Data[kCalcPh]) * (ZPos);
  

	
	const Double_t posy = SieveHoleCorrectionTCS.Y();
	
	dy += ProjectionY - posy;
	rmsdy += (ProjectionY - posy)*(ProjectionY - posy);


	// if(FoilID == 6){
	//   cout << "y_tg = " << y_tg << ", ProjectionY = " << ProjectionY << ", posy = " << posy << ", FoilID = " << FoilID << endl;
	// }

	
    }

    //    cout << "Nignored = " << Nignored << endl;
    
    // DEBUG_INFO("SumSquareDPhi", "#%d : dphi = %f,\t rmsphi = %f", NCall, dphi / fNRawData, TMath::Sqrt(rmsphi / fNRawData));
    
    printf("SumSquareDPhi: #%d : dphi = %f,\t rmsphi = %f,\t overall rmsphi = %f\n", NCall, dphi / fNRawData, TMath::Sqrt(rmsphi / fNRawData), rmsphi);
    printf("SumSquareDY: #%d : dy = %f,\t rmsdy = %f,\t overall rmsdy = %f\n", NCall, dy / fNRawData, TMath::Sqrt(rmsdy / fNRawData), rmsdy);

    //return rmsphi;  //Use for opt to angles
    // return rmsdy;   //Use for opt to sieve plane
    return std::make_pair(rmsphi,rmsdy);   //Use for opt to sieve plane
}



// function designed to display hole cuts on Focal Plane
// c1: shows 

void LOpticsOpt::Display_holes_FP(){


  // create 2D plots for each foil
  TH2F* H_FP_TY[NFoils] = {0}; // Theta vs Y
  TH2F* H_FP_PY[NFoils] = {0}; // Phi vs Y


  Double_t yfp_h = 0.05;
  Double_t yfp_l = -0.05;

  Double_t thfp_h = 30;
  Double_t thfp_l = -35;

  Double_t phfp_h = 50;
  Double_t phfp_l = -50;
  
  
  // 2D FP plots for combined foils
  TH2F* H_FP_TY_all = new TH2F("FP T vs Y","FP T vs Y",200, yfp_l,yfp_h, 200, thfp_l,thfp_h); // Theta vs Y
  TH2F* H_FP_PY_all = new TH2F("FP P vs Y", "FP P vs Y", 200, yfp_l,yfp_h, 200, phfp_l,phfp_h); // Phi vs Y
  
  
  // create arrays for each invidual hole
  TH1F* H_FP_T[NFoils][NSieveCol][NSieveRow] = {0}; // FP theta
  TH1F* H_FP_P[NFoils][NSieveCol][NSieveRow] = {0}; // FP phi
  TH1F* H_FP_Y[NFoils][NSieveCol][NSieveRow] = {0}; // FP Y

  // create arrays for line connecting holes
  TLine* Cols[NFoils][NSieveCol][NHoles-1] = {NULL};
  

  for (UInt_t idx = 0; idx < NFoils; idx++) {

    H_FP_TY[idx] = new TH2F(Form("FP T vs Y Foil  %d", idx), "FP T vs Y for " + Foil_names[idx], 200, yfp_l,yfp_h, 200, thfp_l,thfp_h);

    H_FP_TY[idx]->GetXaxis()->SetTitle("FP y [m]");
    H_FP_TY[idx]->GetYaxis()->SetTitle("FP #theta [mrad]");

    H_FP_PY[idx] = new TH2F(Form("FP P vs Y Foil  %d", idx), "FP P vs Y for " + Foil_names[idx], 200, yfp_l,yfp_h, 200, phfp_l,phfp_h);

    H_FP_PY[idx]->GetXaxis()->SetTitle("FP y [m]");
    H_FP_PY[idx]->GetYaxis()->SetTitle("FP #phi [mrad]");

    for(Int_t col = 0; col < NSieveCol; col++){
      for(Int_t row = 0; row < NSieveRow; row++)
	{
	  H_FP_T[idx][col][row] = new TH1F(Form("T Hole c:%d,r:%d",col,row),Form("T Hole c:%d,r:%d",col,row),300,thfp_l,thfp_h);
	  H_FP_Y[idx][col][row] = new TH1F(Form("Y Hole c:%d,r:%d",col,row),Form("Y Hole c:%d,r:%d",col,row),300,yfp_l,yfp_h);
	}
    }


    
  }

    
  

  for (UInt_t idx = 0; idx < fNRawData; idx++) {

    const EventData &eventdata = fRawData[idx];
      
    UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID]; //starting 0!


    const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
    const UInt_t Row = (UInt_t)eventdata.Data[kRowID]; //starting 0!

    Double_t fp_x = eventdata.Data[kX];
    Double_t fp_y = eventdata.Data[kY];
    Double_t fp_th = eventdata.Data[kTh];
    Double_t fp_phi = eventdata.Data[kPhi];
    
    

    Int_t PlotID = FoilID; // useful for case when only Vertical foils are used
                    // (FoilID based on scheme with optics and vertical foils where optics foils
                    // are 0-7 and vertical foils 8-10) so they can have PlotID 0,1,2
    if(NFoils<11){

      PlotID -= 8;
    }


    H_FP_T[PlotID][Col][Row]->Fill(1000*fp_th);
    
    H_FP_Y[PlotID][Col][Row]->Fill(fp_y);
    
    H_FP_TY[PlotID]->Fill(fp_y,1000*fp_th);
    H_FP_PY[PlotID]->Fill(fp_y,1000*fp_phi);

    H_FP_TY_all->Fill(fp_y,1000*fp_th);
    H_FP_PY_all->Fill(fp_y,1000*fp_phi);    
    

  }

  
  const UInt_t nplot = NFoils;

  // create canvas for 2D plots for each foil
  TCanvas * c1_FP_TY;
    
  if (nplot <= 1) {
    c1_FP_TY = new TCanvas( "#theta_{FP} v y_{FP}","theta Y FP", 1000, 1000);
    c1_FP_TY->Divide(1, 1);
  } else if (nplot <= 3) {
    c1_FP_TY = new TCanvas( "#theta_{FP} v y_{FP}","theta Y FP", 1800, 1100);
    c1_FP_TY->Divide(3, 1);
  } else if (nplot <= 6) {
    c1_FP_TY = new TCanvas( "#theta_{FP} v y_{FP}","theta Y FP", 1800, 1100);
    c1_FP_TY->Divide(3, 2);
  } else if (nplot <= 9){
    c1_FP_TY = new TCanvas( "#theta_{FP} v y_{FP}","theta Y FP", 1800, 1100);
    c1_FP_TY->Divide(3, 3);
  }
  else{
    c1_FP_TY = new TCanvas( "#theta_{FP} v y_{FP}","theta Y FP", 1800, 1100);
    c1_FP_TY->Divide(4, 3);
  }

  // create canvas for FP P vs Y
  TCanvas * c2_FP_PY;
    
  if (nplot <= 1) {
    c2_FP_PY = new TCanvas( "#phi_{FP} v y_{FP}","phi Y FP", 1000, 1000);
    c2_FP_PY->Divide(1, 1);
  } else if (nplot <= 3) {
    c2_FP_PY = new TCanvas( "#phi_{FP} v y_{FP}","phi Y FP", 1800, 1100);
    c2_FP_PY->Divide(3, 1);
  } else if (nplot <= 6) {
    c2_FP_PY = new TCanvas( "#phi_{FP} v y_{FP}","phi Y FP", 1800, 1100);
    c2_FP_PY->Divide(3, 2);
  } else if (nplot <= 9){
    c2_FP_PY = new TCanvas( "#phi_{FP} v y_{FP}","phi Y FP", 1800, 1100);
    c2_FP_PY->Divide(3, 3);
  }
  else{
    c2_FP_PY = new TCanvas( "#phi_{FP} v y_{FP}","phi Y FP", 1800, 1100);
    c2_FP_PY->Divide(4, 3);
  }


  // sort through holes for each foil and find lines between

  
  for (UInt_t idx = 0; idx < NFoils; idx++) {

    c1_FP_TY->cd(idx+1);
    cout << "TY drawing for foil " << idx+1 << " : pad name is " << gPad->GetName() << ", title is " << gPad->GetTitle() << endl;
    H_FP_TY[idx]->Draw("colz");  

    
    // check all column rows for that foil and draw position of holes

    for(Int_t col = 0; col < NSieveCol; col++){
      for(Int_t row = 0; row < NSieveRow; row++)
	{

	  if(!H_FP_T[idx][col][row]->GetEntries()==0 && !H_FP_Y[idx][col][row]->GetEntries()==0){

	    Double_t FP_Y_av = H_FP_Y[idx][col][row]->GetMean();
	    Double_t FP_Y_rms = H_FP_Y[idx][col][row]->GetRMS();
	   
	    Double_t FP_T_av = H_FP_T[idx][col][row]->GetMean();
	    Double_t FP_T_rms = H_FP_T[idx][col][row]->GetRMS();
		    
	    
	    TLine* lT = new TLine(FP_Y_av,FP_T_av-FP_T_rms,FP_Y_av,FP_T_av+FP_T_rms);
	    TLine* lP = new TLine(FP_Y_av-FP_Y_rms,FP_T_av,FP_Y_av+FP_Y_rms,FP_T_av);

	    Int_t Hole = Get_Hole(col,row);

	    auto color = kBlack;
	    auto width = 1;
	    
	    if (Hole == 112){
	      color = kRed;
	      width = 2;
	    }
	    if( Hole == 160){
	      color = kRed;
	      width = 2;
	    }

	    lT->SetLineColor(color);
	    lT->SetLineWidth(width);

	    lP->SetLineColor(color);
	    lP->SetLineWidth(width);
	    
	    lT->Draw("same");
	    lP->Draw("same");

	    
	  }

	  // draw lines between holes in same column
	  // checks if previous and next holes have entries first
	  if(row >1 && !H_FP_T[idx][col][row-2]->GetEntries()==0 && !H_FP_Y[idx][col][row-2]->GetEntries()==0 && !H_FP_T[idx][col][row]->GetEntries()==0){


	    Double_t FP_Y_av_p = H_FP_Y[idx][col][row-2]->GetMean();
	    Double_t FP_Y_rms_p = H_FP_Y[idx][col][row-2]->GetRMS();
	   
	    Double_t FP_T_av_p = H_FP_T[idx][col][row-2]->GetMean();
	    Double_t FP_T_rms_p = H_FP_T[idx][col][row-2]->GetRMS();

	    Double_t FP_Y_av = H_FP_Y[idx][col][row]->GetMean();
	    Double_t FP_Y_rms = H_FP_Y[idx][col][row]->GetRMS();
	   
	    Double_t FP_T_av = H_FP_T[idx][col][row]->GetMean();
	    Double_t FP_T_rms = H_FP_T[idx][col][row]->GetRMS();

	    // Double_t FP_Y_av_n = H_FP_Y[idx][col][row+2]->GetMean();
	    // Double_t FP_Y_rms_n = H_FP_Y[idx][col][row+2]->GetRMS();
	   
	    // Double_t FP_T_av_n = H_FP_T[idx][col][row+2]->GetMean();
	    // Double_t FP_T_rms_n = H_FP_T[idx][col][row+2]->GetRMS();


	    TLine* lP = new TLine(FP_Y_av_p,FP_T_av_p,FP_Y_av,FP_T_av);
	    lP->SetLineColor(kOrange);
	    lP->Draw("same");


	    

	  }



	  // draw lines between holes in same row
	  // checks if previous hole has entries first
	  if(col >1 && !H_FP_T[idx][col-2][row]->GetEntries()==0 && !H_FP_Y[idx][col-2][row]->GetEntries()==0 && !H_FP_T[idx][col][row]->GetEntries()==0){


	    Double_t FP_Y_av_p = H_FP_Y[idx][col-2][row]->GetMean();
	    Double_t FP_Y_rms_p = H_FP_Y[idx][col-2][row]->GetRMS();
	   
	    Double_t FP_T_av_p = H_FP_T[idx][col-2][row]->GetMean();
	    Double_t FP_T_rms_p = H_FP_T[idx][col-2][row]->GetRMS();

	    Double_t FP_Y_av = H_FP_Y[idx][col][row]->GetMean();
	    Double_t FP_Y_rms = H_FP_Y[idx][col][row]->GetRMS();
	   
	    Double_t FP_T_av = H_FP_T[idx][col][row]->GetMean();
	    Double_t FP_T_rms = H_FP_T[idx][col][row]->GetRMS();

	    // Double_t FP_Y_av_n = H_FP_Y[idx][col][row+2]->GetMean();
	    // Double_t FP_Y_rms_n = H_FP_Y[idx][col][row+2]->GetRMS();
	   
	    // Double_t FP_T_av_n = H_FP_T[idx][col][row+2]->GetMean();
	    // Double_t FP_T_rms_n = H_FP_T[idx][col][row+2]->GetRMS();


	    TLine* lP = new TLine(FP_Y_av_p,FP_T_av_p,FP_Y_av,FP_T_av);
	    lP->SetLineColor(kRed);
	    lP->Draw("same");


	    

	  }

	  
	}
    }
   
    gPad->SetEditable(kFALSE);
    
  }

  // draw individual foil histos for FP phi vs Y
  for (UInt_t idx = 0; idx < NFoils; idx++) {
    
    c2_FP_PY->cd(idx+1);
    H_FP_PY[idx]->Draw("colz");
    gPad->SetEditable(kFALSE);

  }


  // draw combined foil histos for FP theta vs Y and phi vs Y


    
  TCanvas* c1_FP_TY_all = new TCanvas( "#theta_{FP} v y_{FP} combined","theta Y FP combined", 1800, 1100);

  c1_FP_TY_all->cd(1);
  H_FP_TY_all->Draw("colz");
  

  TCanvas*  c2_FP_PY_all = new TCanvas( "#phi_{FP} v y_{FP} combined","phi Y FP combined", 1800, 1100);
  c2_FP_PY_all->cd(1);
  H_FP_PY_all->Draw("colz");
  gPad->SetEditable(kFALSE);
    



}

void LOpticsOpt::PrepareVertex(void)
{
    // calculate kRealTgY, kRealReactZ

    for (UInt_t idx = 0; idx < fNRawData; idx++) {
        EventData &eventdata = fRawData[idx];

	const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID]; //starting 0!
	// res = res%(NSieveRow*NSieveCol);
	const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
	const UInt_t Row = (UInt_t)eventdata.Data[kRowID]; //starting 0!
	
	
        assert(FoilID < NFoils); //check array index size

        const TVector3 SieveHoleCorrectionTCS = GetSieveHoleCorrectionTCS(FoilID, Col, Row);
        eventdata.Data[kSieveX] = SieveHoleCorrectionTCS.X();
        eventdata.Data[kSieveY] = SieveHoleCorrectionTCS.Y();
        eventdata.Data[kSieveZ] = SieveHoleCorrectionTCS.Z();

	//	TVector3 BeamSpotHCS(eventdata.Data[kBeamX], eventdata.Data[kBeamY], targetfoils[FoilID]);
	const TVector3 BeamSpotHCS( eventdata.Data[kBeamX] + (targetfoils[FoilID]/BeamZDir_average)*eventdata.Data[kBeamDirX], eventdata.Data[kBeamY] + (targetfoils[FoilID]/BeamZDir_average)*eventdata.Data[kBeamDirY], targetfoils[FoilID]);


	//	TVector3 BeamSpotHCS(BeamX_average, BeamY_average , targetfoils[FoilID]);
        eventdata.Data[kBeamZ] = targetfoils[FoilID];
        TVector3 BeamSpotTCS = fTCSInHCS.Inverse()*(BeamSpotHCS - fPointingOffset);

        const TVector3 MomDirectionTCS = SieveHoleCorrectionTCS - BeamSpotTCS;

        Double_t Real_Tg_Phi = MomDirectionTCS.Y() / MomDirectionTCS.Z();

        Double_t Real_Tg_Y = BeamSpotTCS.Y() + Real_Tg_Phi * (0 - BeamSpotTCS.Z());
	//Real_Tg_Y is in TCS, corresponding Tg_X == 0

	eventdata.Data[kRealTgY] = Real_Tg_Y;
        eventdata.Data[kRealReactZ] = targetfoils[FoilID];

        DEBUG_MASSINFO("PrepareVertex", "Real_Tg_Y = %f,\t Real_ReactZ = %f", Real_Tg_Y, eventdata.Data[kRealReactZ]);


	// const Int_t a = (HRSAngle > 0) ? 1 : -1;
	// CalcReacZ = - ( eventdata.Data[kCalcTgY] -a*MissPointZ)*TMath::Cos(TMath::ATan(Real_Tg_Phi))/TMath::Sin(HRSAngle + TMath::ATan(Real_Tg_Phi)) + BeamSpotHCS.X()*TMath::Cos(HRSAngle+TMath::ATan(Real_Tg_Phi))/TMath::Sin(HRSAngle+TMath::ATan(Real_Tg_Phi));


	//	CalcReacZ = - ( eventdata.Data[kCalcTgY] -a*MissPointZ)*TMath::Cos(Real_Tg_Phi)/TMath::Sin(HRSAngle + Real_Tg_Phi) + BeamSpotHCS.X()*TMath::Cos(HRSAngle+Real_Tg_Phi)/TMath::Sin(HRSAngle+Real_Tg_Phi);
	
	// ignore events that are assumed to be from wrong foil

	//	if(event.Data[k

	  
    }

    // make sure kCalcTh, kCalcPh is filled
    SumSquareDTgY();

    DEBUG_INFO("PrepareVertex", "Done!");
}

TCanvas * LOpticsOpt::CheckVertex()
{
  // Visualize Vertex spectrum
  // gStyle->SetOptTitle();
    
    DEBUG_INFO("CheckVertex", "Entry Point");

    Int_t PlotID;
    
    
    const UInt_t nplot = NFoils;
    TH1D * HTgY[NFoils] = {0};
    TH1D * HTgYReal[NFoils] = {0};
    TH1D *HTgYHole[NFoils][NSieveCol][NSieveRow] = {{{0}}};
    TH1D *HTgYHoleReal[NFoils][NSieveCol][NSieveRow] = {{{0}}};
    UInt_t EventHole[NFoils][NSieveCol][NSieveRow] = {{{0}}};
    const Double_t YRange = 10*10e-3; 

    for (UInt_t idx = 0; idx < NFoils; idx++) {
        HTgY[idx] = new TH1D(Form("Target_Y%d", idx), Foil_names[idx], 500, -YRange, YRange);
        HTgYReal[idx] = new TH1D(Form("Target_Y%d", idx), Foil_names[idx], 500, -YRange, YRange);

        HTgY[idx]->SetXTitle("Target Y [m]");
	HTgY[idx]->GetXaxis()->CenterTitle();

        HTgYReal[idx]->SetXTitle("Target Y [m]");
	HTgYReal[idx]->GetXaxis()->CenterTitle();
        assert(HTgY[idx]); // assure memory allocation

	for(UInt_t ncol =0; ncol<NSieveCol; ncol++){
	  for(UInt_t nrow=0; nrow<NSieveRow; nrow++){
	    HTgYHole[idx][ncol][nrow] = new TH1D(Form("TgY%d_%d_%d", idx, ncol, nrow), Form("TgY%d_%d_%d", idx, ncol, nrow), 400, -YRange, YRange);
	    HTgYHoleReal[idx][ncol][nrow] = new TH1D(Form("TgY%d_%d_%d", idx, ncol, nrow),Form("TgY%d_%d_%d", idx, ncol, nrow) , 400, -YRange, YRange);
	    assert(HTgYHole[idx][ncol][nrow]);
	    assert(HTgYHoleReal[idx][ncol][nrow]);
	  }
	}
    }

    Double_t dtg_y = 0;
    Double_t dtg_y_rms = 0;

    for (UInt_t idx = 0; idx < fNRawData; idx++) {
        EventData &eventdata = fRawData[idx];

        // UInt_t res = (UInt_t) eventdata.Data[kCutID];
        // // const UInt_t KineID = res / (NSieveRow * NSieveCol * NFoils); //starting 0!
        // res = res % (NSieveRow * NSieveCol * NFoils);
        // const UInt_t FoilID = res / (NSieveRow * NSieveCol); //starting 0!

        // res = res % (NSieveRow * NSieveCol);
        // const UInt_t Col = res / (NSieveRow); //starting 0!
        // const UInt_t Row = res % (NSieveRow); //starting 0!

	const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID]; //starting 0!
	PlotID = FoilID;	
	
	if(NFoils<11){
	  PlotID -= 8;
      }
	// res = res%(NSieveRow*NSieveCol);
	const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
	const UInt_t Row = (UInt_t)eventdata.Data[kRowID]; //starting 0!


	HTgY[PlotID]->Fill(eventdata.Data[kCalcTgY]);
	HTgYReal[PlotID]->Fill(eventdata.Data[kRealTgY]);

	HTgYHole[PlotID][Col][Row]->Fill(eventdata.Data[kCalcTgY]);
	HTgYHoleReal[PlotID][Col][Row]->Fill(eventdata.Data[kRealTgY]);
	EventHole[PlotID][Col][Row]++;
	
	dtg_y += eventdata.Data[kCalcTgY] - eventdata.Data[kRealTgY];
	dtg_y_rms += (eventdata.Data[kCalcTgY] - eventdata.Data[kRealTgY])*(eventdata.Data[kCalcTgY] - eventdata.Data[kRealTgY]);
    }

    DEBUG_INFO("CheckTgY", "dtg_v = %f,\t dtg_v_rms = %f", dtg_y / fNRawData, TMath::Sqrt(dtg_y_rms / fNRawData));

    TCanvas * c1;
    if (nplot <= 3) {
        c1 = new TCanvas("CheckTgY", "Target Y Check", 1800, 450);
        c1->Divide(3, 1);
    } else if (nplot <= 6) {
        c1 = new TCanvas("CheckTgY", "Target Y Check", 1800, 900);
        c1->Divide(3, 2);
    } else if (nplot <= 9){
        c1 = new TCanvas("CheckTgY", "Target Y Check", 1800, 1100);
        c1->Divide(3, 3);
    }
    else{
        c1 = new TCanvas("CheckTgY", "Target Y Check", 1800, 1100);
	c1->Divide(4, 3);
    }
    

    Double_t YMean[NFoils]={0};
    Double_t MaxPlot = 2000.0;
    Double_t DefResolution = 0.5e-3;
    Double_t DefResolutionZ = 5e-2; //0.5e-2
    Double_t FitRangeMultiply = 6;
    Double_t ymean =0;
    for (UInt_t idx = 0; idx < nplot; idx++) {
        // UInt_t FoilID = idx;

        c1->cd(idx + 1);
        assert(HTgY[idx]);

	YMean[idx] = HTgY[idx]->GetMean();
	HTgYReal[idx]->Draw("same");
	HTgYReal[idx]->GetXaxis()->SetLabelSize(0.05);
	HTgYReal[idx]->GetXaxis()->SetTitleSize(0.05);
	HTgYReal[idx]->GetYaxis()->SetLabelSize(0.05);
        HTgY[idx]->Draw("same");
	HTgY[idx]->GetXaxis()->SetLabelSize(0.05);
	HTgY[idx]->GetXaxis()->SetTitleSize(0.05);
	HTgY[idx]->GetYaxis()->SetLabelSize(0.05);
	HTgY[idx]->GetXaxis()->SetRangeUser(YMean[idx]-0.015,YMean[idx]+0.015 );
	HTgYReal[idx]->GetXaxis()->SetRangeUser(YMean[idx]-0.015,YMean[idx]+0.015 );
	HTgYReal[idx]->SetLineColor(kGreen);
	
	ymean = HTgY[idx]->GetMean();
 
				
        TString FitFunc = Form("YtPeak%d", idx);
        TF1 *f = new TF1(FitFunc, "gaus", ymean - DefResolution*FitRangeMultiply, ymean + DefResolution * FitRangeMultiply);
        f->SetParameter(1, ymean);
	HTgY[idx] -> Fit(FitFunc, "R");
	f->SetLineColor(2);
       	f->Draw("SAME");
	
	ymean = HTgYReal[idx]->GetMean();
	cout<<"ymean:"<<ymean<<endl;
	c1->cd(idx+1)->Update();
	MaxPlot = c1->cd(idx + 1) ->GetUymax();
        TLine *l = new TLine(ymean, 0, ymean, MaxPlot);
        l->SetLineColor(kBlue);
        l->SetLineWidth(2);
	l->Draw("same");

	TPaveText *t1 = new TPaveText(0.65,0.6,0.9,0.82,"NDC");

	if(ymean >0.05){
	  
	  t1->SetBBoxCenterX(0.1);
	}



	/// TPaveText(x1,,y1,x2,y2)

	t1->AddText(Form("\\Delta = %2.1f mm", 1000 * (f->GetParameter(1) - ymean)));
	t1->AddText(Form("\\sigma = %2.1f mm", 1000 * (f->GetParameter(2))));
	t1->SetShadowColor(0);
	t1->SetFillColor(0);
	t1->Draw("same");					
    }
    
    ///////////////////////////check z direction
    TH1D * HTgZ[NFoils] = {0};         // z position from calculated y, calculated y is from matrix
    TH1D * HTgZReal[NFoils] = {0};     //real z postion
    TH1D * HTgZfromTgY[NFoils] = {0};  // z position from 'real' tg_y
    TH2D * HTgZ_v_phi[NFoils] = {0};  // 2d plots reactz vs 'real' phi
    
    //    const Double_t ZRange = 15*10e-3; 
    const Double_t ZRange = 35*10e-3; 
    TH1F* Check_Tong=new TH1F("Check_Tong","React Z",600,-0.12,0.12);
    TH1F* Check_Tong2=new TH1F("Check_Tong","React Z",600,-0.12,0.12);
    Check_Tong2->SetLineColor(2);
    Check_Tong2->SetFillColor(2);
    Check_Tong->SetLineWidth(2);
    TCanvas * c_Tong=new TCanvas("c_Tong", "c_Tong", 1800, 1350);
    
    Double_t Z_offset = 0;

    
    Check_Tong->SetXTitle("Target Z [m]");
    for (UInt_t idx = 0; idx < nplot; idx++) {
      HTgZ[idx] = new TH1D(Form("Target_Z%d", idx),Foil_names[idx], 500, -ZRange+Z_offset, ZRange+Z_offset);
      HTgZReal[idx] = new TH1D(Form("Target_Z%d", idx), Form("Target Z for Foil Target #%d", idx), 500, -ZRange+Z_offset, ZRange+Z_offset);
      HTgZfromTgY[idx] = new TH1D(Form("Target_Z%d", idx), Form("Target Z for Foil Target #%d", idx), 500, -ZRange+Z_offset, ZRange+Z_offset);

      HTgZ_v_phi[idx] =  new TH2D(Form("Z react vs #phi_{tg} %d",idx),Form("Z react vs #phi_{tg} %d",idx),500, -0.04,0.04,500, -ZRange+Z_offset, ZRange+Z_offset);
      HTgZ_v_phi[idx]->GetXaxis()->SetTitle("#phi_{tg} [rad]");
      HTgZ_v_phi[idx]->GetYaxis()->SetTitle("Z vertex [m]");
      
      
      HTgZ[idx]->GetXaxis()->SetTitle("Z vertex [m]");
      HTgZ[idx]->GetXaxis()->CenterTitle();
      assert(HTgZ[idx]); // assure memory allocation
      assert(HTgZReal[idx]);
    }

    Double_t dtg_z = 0;
    Double_t dtg_z_rms = 0;

    Double_t Z_from_tgy; // variable to get 'real' reactz from 'real' tg_y

    for (UInt_t idx = 0; idx < fNRawData; idx++) {
        EventData &eventdata = fRawData[idx];

        // UInt_t res = (UInt_t) eventdata.Data[kCutID];
        // // const UInt_t KineID = res / (NSieveRow * NSieveCol * NFoils); //starting 0!
        // res = res % (NSieveRow * NSieveCol * NFoils);
        // const UInt_t FoilID = res / (NSieveRow * NSieveCol); //starting 0!

        // res = res % (NSieveRow * NSieveCol);
        // const UInt_t Col = res / (NSieveRow); //starting 0!
        // const UInt_t Row = res % (NSieveRow); //starting 0!

	const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID]; //starting 0!
	// res = res%(NSieveRow*NSieveCol);
	const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
	const UInt_t Row = (UInt_t)eventdata.Data[kRowID]; //starting 0!


	PlotID = FoilID;

	if(NFoils<11){
	  PlotID -= 8;
	}

	//calculate ReactZ with CalcTgY
	eventdata.Data[kRealReactZ] = targetfoils[FoilID];
	Double_t CalcReacZ =0;
        const TVector3 SieveHoleCorrectionTCS = GetSieveHoleCorrectionTCS(FoilID, Col, Row);
        eventdata.Data[kSieveX] = SieveHoleCorrectionTCS.X();
        eventdata.Data[kSieveY] = SieveHoleCorrectionTCS.Y();
        eventdata.Data[kSieveZ] = SieveHoleCorrectionTCS.Z();

	//	TVector3 BeamSpotHCS(eventdata.Data[kBeamX], eventdata.Data[kBeamY], targetfoils[FoilID]);
	const TVector3 BeamSpotHCS( eventdata.Data[kBeamX] + (targetfoils[FoilID]/BeamZDir_average)*eventdata.Data[kBeamDirX], eventdata.Data[kBeamY] + (targetfoils[FoilID]/BeamZDir_average)*eventdata.Data[kBeamDirY], targetfoils[FoilID]);


	// 	TVector3 BeamSpotHCS(BeamX_average, BeamY_average, targetfoils[FoilID]);
       eventdata.Data[kBeamZ] = targetfoils[FoilID];
        TVector3 BeamSpotTCS = fTCSInHCS.Inverse()*(BeamSpotHCS - fPointingOffset);

        const TVector3 MomDirectionTCS = SieveHoleCorrectionTCS - BeamSpotTCS;

        Double_t Real_Tg_Phi = MomDirectionTCS.Y() / MomDirectionTCS.Z();

	//        Double_t Real_Tg_Y = BeamSpotTCS.Y() + Real_Tg_Phi * (0 - BeamSpotTCS.Z());
        const Int_t a = (HRSAngle > 0) ? 1 : -1;
	//	CalcReacZ = - ( eventdata.Data[kCalcTgY] -a*MissPointZ)*TMath::Cos(TMath::ATan(Real_Tg_Phi))/TMath::Sin(HRSAngle + TMath::ATan(Real_Tg_Phi)) + BeamSpotHCS.X()*TMath::Cos(HRSAngle+TMath::ATan(Real_Tg_Phi))/TMath::Sin(HRSAngle+TMath::ATan(Real_Tg_Phi));
	CalcReacZ = - ( eventdata.Data[kCalcTgY] -a*MissPointZ)*TMath::Cos(Real_Tg_Phi)/TMath::Sin(HRSAngle + Real_Tg_Phi) + BeamSpotHCS.X()*TMath::Cos(HRSAngle + Real_Tg_Phi)/TMath::Sin(HRSAngle+Real_Tg_Phi);

	Z_from_tgy =  - ( eventdata.Data[kRealTgY] -a*MissPointZ)*TMath::Cos(Real_Tg_Phi)/TMath::Sin(HRSAngle + Real_Tg_Phi) + BeamSpotHCS.X()*TMath::Cos(HRSAngle + Real_Tg_Phi)/TMath::Sin(HRSAngle+Real_Tg_Phi);
	
	//for right arm,eventdata.Data[kCalcTgY] + MissPointZ, for left arm,eventdata.Data[kCalcTgY] - MissPointZ 
	

	HTgZfromTgY[PlotID]->Fill(Z_from_tgy);	

	HTgZ_v_phi[PlotID]->Fill(Real_Tg_Phi,CalcReacZ);

	HTgZ[PlotID]->Fill(CalcReacZ);
	HTgZReal[PlotID]->Fill(eventdata.Data[kRealReactZ]);
	Check_Tong->Fill(CalcReacZ);
	Check_Tong2->Fill(eventdata.Data[kRealReactZ]);
	dtg_z += CalcReacZ - eventdata.Data[kRealReactZ];
	dtg_z_rms += (CalcReacZ - eventdata.Data[kRealReactZ])*(CalcReacZ - eventdata.Data[kRealReactZ]);
    }

    DEBUG_INFO("CheckTgZ", "dtg_z = %f,\t dtg_z_rms = %f", dtg_z / fNRawData, TMath::Sqrt(dtg_z_rms / fNRawData));

    TCanvas * c2;
    if (nplot <= 3) {
        c2 = new TCanvas("CheckTgZ", "Target Z Check", 1800, 450);
        c2->Divide(3, 1);
    } else if (nplot <= 6) {
        c2 = new TCanvas("CheckTgZ", "Target Z Check", 1800, 900);
        c2->Divide(3, 2);
    } else if (nplot <= 9){
        c2 = new TCanvas("CheckTgZ", "Target Z Check", 1800, 1100);
        c2->Divide(3, 3);
    }
    else{
        c2 = new TCanvas("CheckTgZ", "Target Z Check", 1800, 1100);
	c2->Divide(4, 3);
    }


   // c_Tong->cd();
   // Check_Tong->Draw();
    
    
    Double_t ZMean[NFoils] = {0};
    Double_t zmean=0;
    for (UInt_t idx = 0; idx < nplot; idx++) {
        c2->cd(idx + 1);
	//       assert(HTgZ[idx]);

	zmean = HTgZ[idx]->GetMean();
	ZMean[idx] = HTgZ[idx]->GetMean();

	HTgZ[idx]->GetXaxis()->SetTitleSize(0.05);

	//	HTgZReal[idx]->Draw();
	HTgZ[idx]->Draw();
	HTgZReal[idx]->SetLineColor(2);
	//	mean = HTgZReal[idx]->GetMean();
	//	cout<<"mean:"<<mean<<endl;
	// zmean = HTgZ[idx]->GetMean();
	// HTgZ[idx]->GetXaxis()->SetRangeUser(ZMean[idx] - 0.02, ZMean[idx] + 0.02);
	// HTgZReal[idx]->GetXaxis()->SetRangeUser(ZMean[idx] - 0.02, ZMean[idx] + 0.02);	
	
        TString FitFunc1 = Form("ZtPeak%d", idx);
        TF1 *f1 = new TF1(FitFunc1, "gaus", ZMean[idx] - DefResolutionZ*FitRangeMultiply, ZMean[idx] + DefResolutionZ * FitRangeMultiply);
        f1->SetParameter(1, ZMean[idx]);
	HTgZ[idx] -> Fit(FitFunc1, "R0");
	f1->SetLineColor(kPink);
	f1->Draw("SAME");

	//	mean = HTgZReal[idx]->GetMean();
	//	mean=targetfoils[idx];
	//	cout<<"mean:"<<mean<<endl;
	c2->cd(idx+1)->Update();
	MaxPlot = c2->cd(idx + 1) ->GetUymax();


	TLine *l1;
	
	if(NFoils <11)
	  {
	    // case in which only vertical foils are displayed
	    l1 = new TLine(targetfoils[idx+8], 0, targetfoils[idx+8], MaxPlot);
	  }
	else
	  {
	    l1 = new TLine(targetfoils[idx], 0, targetfoils[idx], MaxPlot);
	  }
	
        l1->SetLineColor(kOrange);
        l1->SetLineWidth(2);
	l1->Draw("same");
       
	TPaveText *tt = new TPaveText(0.6,0.6,0.88,0.82,"NDC");

	if(NFoils <11)
	  {
	    // case in which only vertical foils are displayed
	    tt->AddText(Form("\\Delta = %2.1f mm", 1000 * (f1->GetParameter(1) - targetfoils[idx+8]))); //
	  }
	else
	  {
	    tt->AddText(Form("\\Delta = %2.1f mm", 1000 * (f1->GetParameter(1) - targetfoils[idx])));	    
	  }

	tt->AddText(Form("\\sigma = %2.1f mm", 1000 * (f1->GetParameter(2))));
	tt->SetShadowColor(0);
	tt->SetFillColor(0);
	tt->Draw("same");

	HTgZfromTgY[idx]->SetLineColor(kGreen+2);
	HTgZfromTgY[idx]->Draw("same");
	
	

       // c_Tong->cd()->Update();
     }
    c_Tong->cd();
    Check_Tong->Draw("same");
    Check_Tong2->Draw("same");



    // create new plot to display all foils simultaneously

    TCanvas* c_vertices = new TCanvas("c_vertices","All vertices check",1800,1800);
    
    c_vertices->cd(0);

    TLegend* leg = new TLegend(.1,.65,.37,.9,"Key");
    leg->SetFillColor(0);


    // Find larget y value in plots

    Double_t max_val = 0;

    for (Int_t i = 0; i<NFoils; i++){

      Double_t max_bin = HTgZ[i]->GetMaximumBin();
      max_bin = HTgZ[i]->GetBinContent(max_bin);

      if(max_bin > max_val){
	max_val = max_bin;
      }

    }


    // draw lines where 'real' foils are
    //    Double_t MaxPlot_f = c_vertices->cd()->GetUymax();
    TLine* foils[NFoils];

    //    gStyle->SetPallete(

    
    // create small palette of colors to use in drawing 

    // Int_t No_colors = 10;
    // Int_t palette[No_colors];
    // palette[0] = kBlue;
    // palette[1] = kRed;
    // palette[2] = kGreen;
    // palette[3] = kOrange;
    
    
    // gStyle->SetPalette(No_colors,palette);



    for (Int_t i = 0; i<NFoils; i++){

      
      HTgZ[i]->GetYaxis()->SetRangeUser(0,max_val);

      foils[i] = new TLine(targetfoils[i], 0, targetfoils[i], MaxPlot);
      foils[i]->SetLineColor(kBlack);

      

      // HTgZ[i]->SetLineColor(3+i);
      HTgZ[i]->SetLineWidth(2);

      if(i==0){
	//	HTgZ[i]->Draw("L PLC"); // PLC option to do automatic colouring
		HTgZ[i]->Draw("L"); 
	foils[i]->Draw("same");

      }
      else{
	// HTgZ[i]->Draw("same L PLC");
	HTgZ[i]->Draw("same L");
	foils[i]->Draw("same");
      }          

      //      leg->AddEntry(HTgZ[i],Form("V%d",i),"l");
      leg->AddEntry(HTgZ[i],Foil_names[i],"l");
    }
    
    leg->Draw("same");
    // TLegend* leg = new TLegend(.1,.65,.37.9,"Key");
    // leg->SetFillColor(0);
    // for (Int_t i = 0; i<NFoils; i++){
    // leg->AddEntry(


    // draw reactz vs phi plots for each foil
    
    TCanvas *c4;
   
    if (nplot <= 3) {
        c4 = new TCanvas("Z_v_phi","React Z v #phi_{tg}", 1800, 450);
        c4->Divide(3, 1);
    } else if (nplot <= 6) {
        c4 = new TCanvas("Z_v_phi","React Z v #phi_{tg}", 1800, 900);
        c4->Divide(3, 2);
    } else if (nplot <= 9){
        c4 = new TCanvas("Z_v_phi","React Z v #phi_{tg}", 1800, 1100);
        c4->Divide(3, 3);
    }
    else{
        c4 = new TCanvas("Z_v_phi","React Z v #phi_{tg}", 1800, 1100);
	c4->Divide(4, 3);
    }

    for (UInt_t idx = 0; idx < nplot; idx++) {

      c4->cd(idx+1);

      HTgZ_v_phi[idx]->Draw("colz");

      
    }

     
    

    
#if tgy_hole_info
    Double_t YHole_Peak[NFoils][NSieveCol][NSieveRow] = {{{0}}};
    Double_t YHole_Mean[NFoils][NSieveCol][NSieveRow] = {{{0}}};

    for(UInt_t nfoil =0; nfoil <5; nfoil++){

      TCanvas *c3 = new TCanvas(Form("c3_%d",nfoil+1),Form("TgY check for foil target #%d of each sieve hole",nfoil+1),1800,1350);
      c3->Divide(NSieveCol,NSieveRow);

      for(UInt_t col=0; col<NSieveCol; col++){
	for(UInt_t row =0; row<NSieveRow; row++){
	  c3->cd(row*NSieveCol+col+1);
	  HTgYHole[nfoil][8-col][6-row]->GetXaxis()->SetRangeUser(YMean[nfoil]-0.015,YMean[nfoil]+0.015);
	  HTgYHoleReal[nfoil][8-col][6-row]->GetXaxis()->SetRangeUser(YMean[nfoil]-0.015,YMean[nfoil]+0.015);
	  HTgYHole[nfoil][8-col][6-row]->Draw();
	  if(EventHole[nfoil][8-col][6-row]>0){
	    HTgYHoleReal[nfoil][8-col][6-row]->SetLineColor(6);
	    HTgYHoleReal[nfoil][8-col][6-row]->Draw("same");
	    mean = HTgYHoleReal[nfoil][8-col][6-row] ->GetMean();
	    YHole_Mean[nfoil][col][row] = mean;

	    TString FitFunc = Form("YHole_Peak%d_%d_%d", nfoil,8-col,6-row);
	    TF1 *f = new TF1(FitFunc, "gaus", mean - DefResolution*FitRangeMultiply, mean + DefResolution * FitRangeMultiply);
	    f->SetParameter(1, mean);
	    HTgYHole[nfoil][8-col][6-row] -> Fit(FitFunc, "R");
	    f->SetLineColor(2);
	    f->Draw("SAME");
	    YHole_Peak[nfoil][8-col][6-row] = f->GetParameter(1);
	
	    c3->cd(row*NSieveCol+col+1)->Update();
	    MaxPlot = c3->cd(row*NSieveCol+col+1) ->GetUymax();
	    TLine *l = new TLine(mean, 0, mean, MaxPlot);
	    l->SetLineColor(kBlue);
	    l->SetLineWidth(2);
	    l->Draw();

	    TPaveText *tt = new TPaveText(0.6,0.7,0.95,0.82,"NDC");
	    tt->AddText(Form("\\Delta: %2.1f mm", 1000 * (YHole_Peak[nfoil][8-col][6-row] - mean)));
	    tt->SetShadowColor(0);
	    tt->SetFillColor(0);
	    tt->Draw("same");					
 
	  }
	}
      }
    }
#endif  
    return c1;
}

Double_t LOpticsOpt::SumSquareDTgY(void)
{
    // return square sum of diff between calculated tg_y and expected tg_y

    Double_t dtg_y = 0; //Difference
    Double_t dtg_y_rms = 0; //mean square

    // differnces in z-vertex

    Double_t dz = 0; //Difference
    Double_t dz_rms = 0; //mean square

    static UInt_t NCall = 0;
    NCall++;

    Double_t tg_y;

    for (UInt_t idx = 0; idx < fNRawData; idx++) {
        EventData &eventdata = fRawData[idx];


	// condition to remove events from outside main gaussian (from initial matrix)
      if(eventdata.Data[kIgnore] == 0){
	continue;
      }
	

	
        Double_t x_fp = eventdata.Data[kX];
        const Double_t(*powers)[5] = eventdata.powers;

        // calculate the matrices we need
        // CalcMatrix(x_fp, fDMatrixElems);
        // CalcMatrix(x_fp, fTMatrixElems);
        CalcMatrix(x_fp, fYMatrixElems);
	CalcMatrix(x_fp, fYTAMatrixElems);
        // CalcMatrix(x_fp, fPMatrixElems);
        // CalcMatrix(x_fp, fPTAMatrixElems);

        // calculate the coordinates at the target
        tg_y = CalcTargetVar(fYMatrixElems, powers) + CalcTargetVar(fYTAMatrixElems, powers);

	// UInt_t res = (UInt_t) eventdata.Data[kCutID];
        // res = res % (NSieveRow * NSieveCol * NFoils);
        // const UInt_t FoilID = res / (NSieveRow * NSieveCol); //starting 0!

	const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID]; //starting 0!
	// res = res%(NSieveRow*NSieveCol);
	const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
	const UInt_t Row = (UInt_t)eventdata.Data[kRowID]; //starting 0!

  
	//        const Double_t ArbitaryVertexShift = fArbitaryVertexShift[FoilID] * (-TMath::Sin(HRSAngle));

const Double_t ArbitaryVertexShift = 0;

        dtg_y += tg_y - eventdata.Data[kRealTgY] + ArbitaryVertexShift;
        dtg_y_rms += (tg_y - eventdata.Data[kRealTgY] + ArbitaryVertexShift)*(tg_y - eventdata.Data[kRealTgY] + ArbitaryVertexShift);

        // save the results
        eventdata.Data[kCalcTgY] = tg_y;


	
	eventdata.Data[kRealReactZ] = targetfoils[FoilID];
	Double_t CalcReacZ =0;
        const TVector3 SieveHoleCorrectionTCS = GetSieveHoleCorrectionTCS(FoilID, Col, Row);
	//	TVector3 BeamSpotHCS(eventdata.Data[kBeamX], eventdata.Data[kBeamY], targetfoils[FoilID]);
	const TVector3 BeamSpotHCS( eventdata.Data[kBeamX] + (targetfoils[FoilID]/BeamZDir_average)*eventdata.Data[kBeamDirX], eventdata.Data[kBeamY] + (targetfoils[FoilID]/BeamZDir_average)*eventdata.Data[kBeamDirY], targetfoils[FoilID]);
	

	// 	TVector3 BeamSpotHCS(BeamX_average, BeamY_average, targetfoils[FoilID]);
       eventdata.Data[kBeamZ] = targetfoils[FoilID];
        TVector3 BeamSpotTCS = fTCSInHCS.Inverse()*(BeamSpotHCS - fPointingOffset);

        const TVector3 MomDirectionTCS = SieveHoleCorrectionTCS - BeamSpotTCS;
	Double_t Real_Tg_Phi = MomDirectionTCS.Y() / MomDirectionTCS.Z();

	//        Double_t Real_Tg_Y = BeamSpotTCS.Y() + Real_Tg_Phi * (0 - BeamSpotTCS.Z());
        const Int_t a = (HRSAngle > 0) ? 1 : -1;
	//	CalcReacZ = - ( eventdata.Data[kCalcTgY] -a*MissPointZ)*TMath::Cos(TMath::ATan(Real_Tg_Phi))/TMath::Sin(HRSAngle + TMath::ATan(Real_Tg_Phi)) + BeamSpotHCS.X()*TMath::Cos(HRSAngle+TMath::ATan(Real_Tg_Phi))/TMath::Sin(HRSAngle+TMath::ATan(Real_Tg_Phi));
	CalcReacZ = - ( eventdata.Data[kCalcTgY] -a*MissPointZ)*TMath::Cos(Real_Tg_Phi)/TMath::Sin(HRSAngle + Real_Tg_Phi) + BeamSpotHCS.X()*TMath::Cos(HRSAngle+Real_Tg_Phi)/TMath::Sin(HRSAngle+Real_Tg_Phi);

	
	dz += CalcReacZ - eventdata.Data[kRealReactZ];

	dz_rms += (CalcReacZ - eventdata.Data[kRealReactZ])*(CalcReacZ - eventdata.Data[kRealReactZ]);


	



    }
    //    cout<<"fNRawData:"<<fNRawData<<endl;

    DEBUG_INFO("SumSquareDTgY", "#%d : dtg_y = %f, dtg_y_rms = %f", NCall, dtg_y / fNRawData, TMath::Sqrt(dtg_y_rms / fNRawData));
    printf("SumSquareDTgY: #%d : dtg_y = %f, dtg_y_rms = %f\n", NCall, dtg_y / fNRawData, TMath::Sqrt(dtg_y_rms / fNRawData));
    printf("SumSquareDTgY: #%d : dz = %f, dz_rms=%f\n",NCall,dz/fNRawData,TMath::Sqrt(dz_rms/fNRawData));


    //    return dtg_y_rms;
    return dz_rms;
    
}

// void LOpticsOpt::PrepareDp(void)
// {
//     // calate expected dp_kin, dp_kin offsets ....
//     // Fill up fRawData[].Data[] kKineID thr kRealDpKin

//     // print Central Momentums
//     printf("HRSCentralMom[%d] (GeV) = {", NKine);
//     for (UInt_t KineID = 0; KineID < NKine; KineID++)
//         printf("%f  ", HRSCentralMom[KineID]);
//     printf("}\n");

//     // parameter LossEachUnitB is the Radiation loss in 1m LH2 before scattering
//     // parameter LossEachUnitA is the Radiation loss in 1m LH2 after scattering
//     Double_t scatang = 0;
//     Double_t dpkinoff = 0, dpkinoff_rms = 0;
//     fNCalibData = 0;
//     Double_t exttargcorr_dp = 0, rms_exttargcorr_dp = 0;

//     Double_t D1_sum=0, D2_sum=0; //D1, travel length in target before scattering, D2 travel length in target after scattering 
//     //    Double_t ElossAfterTg_sum=0, ElossAfterTotal_sum=0;
//     Double_t ElossTg_sum=0, ElossTg =0;
//     Double_t ElossAfterTg=0, ElossBeforeTg=0, ElossAfter_sum =0;
//     for (UInt_t idx = 0; idx < fNRawData; idx++) {
//         DEBUG_MASSINFO("PrepareDp", "=========== Event %d ===========", idx);

//         EventData &eventdata = fRawData[idx];
// 	//eventdata.Data[kCutID] is the kine number
//         // const UInt_t ExtraDataFlag =(UInt_t) eventdata.Data[kCutID];
// 	// const UInt_t KineID = (UInt_t) eventdata.Data[kCutID];
// 	const UInt_t ExtraDataFlag =(UInt_t) eventdata.Data[kFoilID];
// 	const UInt_t KineID = (UInt_t) eventdata.Data[kFoilID];

//         assert(ExtraDataFlag < 5); //Check flag range. beyond 5 is not used
//         assert(KineID < 5); //check array index size

//         fNCalibData++;
//         eventdata.Data[kKineID] = KineID;
//         eventdata.Data[kCentralp] = HRSCentralMom[KineID];

// 	// calculate th, ph
// 	Double_t theta,phi,tgy,dp;
// 	Double_t x_fp = eventdata.Data[kX];// R.tr.r_x
//         const Double_t(*powers)[5] = eventdata.powers;

//         CalcMatrix(x_fp, fTMatrixElems);
//         theta = CalcTargetVar(fTMatrixElems, powers);

//         CalcMatrix(x_fp, fPMatrixElems);
//         CalcMatrix(x_fp, fPTAMatrixElems);
//         phi = CalcTargetVar(fPMatrixElems, powers) + CalcTargetVar(fPTAMatrixElems, powers);

// 	CalcMatrix(x_fp, fYMatrixElems);
// 	CalcMatrix(x_fp, fYTAMatrixElems);
// 	tgy = CalcTargetVar(fYMatrixElems, powers) + CalcTargetVar(fYTAMatrixElems, powers);
// 	//	eventdata.Data[kCalcTgY] = tgy;
//         CalcMatrix(x_fp, fDMatrixElems);
//         dp = CalcTargetVar(fDMatrixElems, powers);


// 	TVector3 BeamSpotHCS(eventdata.Data[kBeamX], eventdata.Data[kBeamY], eventdata.Data[kBeamZ]);
// 	//        TVector3 BeamSpotHCS(BeamX_average, BeamY_average, eventdata.Data[kBeamVZ]);
//         TVector3 BeamSpotTCS = fTCSInHCS.Inverse()*(BeamSpotHCS - fPointingOffset);
// 	TVector3 MomDirectionTCS(theta,phi,1);

//         eventdata.Data[kRealTh] = theta;
//         eventdata.Data[kRealPhi] = phi;
//         eventdata.Data[kCalcTh] = theta;
//         eventdata.Data[kCalcPh] = phi;

//         const Double_t x_tg = BeamSpotTCS.X() - BeamSpotTCS.Z() * eventdata.Data[kRealTh];
//         eventdata.Data[kRealTgX] = x_tg;
//         eventdata.Data[kRealThMatrix] = eventdata.Data[kRealTh] - x_tg * ExtTarCor_ThetaCorr;

//         DEBUG_MASSINFO("PrepareDp", "RealTh = %f,\t RealPh = %f", eventdata.Data[kRealThMatrix], eventdata.Data[kRealPhi]);
// 	//        DEBUG_MASSINFO("PrepareDp", "SieveHoleY = %f,\t Mom_Y = %f,\t Mom_Z = %f", SieveHoleTCS.Y(), MomDirectionTCS.Y(), MomDirectionTCS.Z());

//         TVector3 MomDirectionHCS = fTCSInHCS*MomDirectionTCS;
//         TVector3 BeamDirection(0, 0, 1);
//         const Double_t ScatteringAngle = BeamDirection.Angle(MomDirectionHCS);
//         eventdata.Data[kScatterAngle] = ScatteringAngle;
//         scatang += ScatteringAngle;

// 	//calculate ReactZ with CalcTgY
// 	Double_t CalcReacZ =0;
//         Double_t Tg_Phi = eventdata.Data[kCalcPh];

// 	CalcReacZ = - ( tgy + MissPointZ)*TMath::Cos(TMath::ATan(Tg_Phi))/TMath::Sin(HRSAngle + TMath::ATan(Tg_Phi)) + BeamSpotHCS.X()*TMath::Cos(HRSAngle+TMath::ATan(Tg_Phi))/TMath::Sin(HRSAngle+TMath::ATan(Tg_Phi));
// 	if(CalcReacZ>0.15/2.) CalcReacZ = 0.15/2.;
// 	if(CalcReacZ<-0.15/2.) CalcReacZ = -0.15/2.;

// 	eventdata.Data[kRealReactZ] = CalcReacZ;
// 	//	cout<<"CalcReacZ:"<<CalcReacZ<<endl;
// 	//	if( CalcReacZ>7.5*1e-2 || CalcReacZ<-7.5*1e-2) cout<<"CalcReacZ:"<<CalcReacZ<<endl; 
// 	//************************************//
// 	TVector3 ReactionVertex(eventdata.Data[kBeamX], eventdata.Data[kBeamY], CalcReacZ);
// 	Double_t D1=0;
// 	Double_t D2=0;
// 	D1 = LH2_TargetLength/2.+CalcReacZ;
// 	D2 = TravelLength(ReactionVertex, MomDirectionHCS);
// 	D1_sum+=D1;
// 	D2_sum+=D2;
// 	ElossBeforeTg = ElossTargetBefore(ReactionVertex, MomDirectionHCS);
// 	ElossAfterTg = ElossTargetAfter(ReactionVertex, MomDirectionHCS);
// 	ElossAfter_sum += ElossAfterTg;
// 	ElossTg = ElossBeforeTg+ElossAfterTg;
// 	ElossTg_sum +=ElossTg;

//         // eventdata.Data[kTravelLengthBefore] = D1;
//         // eventdata.Data[kTravelLengthAfter] = D2;
// 	// eventdata.Data[kElossTgBefore] = ElossTargetBefore(ReactionVertex, MomDirectionHCS);
// 	// eventdata.Data[kElossTgAfter] = ElossTargetAfter(ReactionVertex, MomDirectionHCS);
// 	// eventdata.Data[kElossDp] = ElossAfterTg / eventdata.Data[kCentralp];


//         // calculate difference between dp_kin and dp
//         // dp_kin + kDpKinOffsets = dp
//         const Double_t DM = ExcitationEnergy[KineID];
//         const Double_t Ma = GroundNuclearMass;

//         const Double_t DpKinOffsets = (ScatMom(DM, Ma, Ebeam[KineID] - ElossBeforeTg, ScatteringAngle) - ScatMom(DM, Ma, Ebeam[KineID] - ElossBeforeTg, TMath::Abs(HRSAngle))) / eventdata.Data[kCentralp];
//         eventdata.Data[kDpKinOffsets] = DpKinOffsets;

//         dpkinoff += DpKinOffsets;
//         dpkinoff_rms += DpKinOffsets*DpKinOffsets;

//         // calculate kRealDpKin, should be same for same kine settings
// 	//        eventdata.Data[kElossDp] = ElossByFoil[FoilID] / eventdata.Data[kCentralp];
// 	//********calculate vertex z position to get the ElossDp*******//


// 	eventdata.Data[kRealDpKin] = ScatMom(DM, Ma, Ebeam[KineID] - ElossBeforeTg, TMath::Abs(HRSAngle)) / eventdata.Data[kCentralp] - 1;
// 	//******************//
//         eventdata.Data[kDpKinOffsets] = DpKinOffsets - eventdata.Data[kElossDp] ;
// 	//********************//
//         DEBUG_MASSINFO("PrepareDp", "ScatterMom = %f,\t Centralp = %f,\t radloss = %f,\t ebeam = %f", ScatMom(DM, Ma, Ebeam[KineID], TMath::Abs(HRSAngle)), eventdata.Data[kCentralp], eventdata.Data[kElossDp], Ebeam[KineID]);

// 	eventdata.Data[kRealDpKinMatrix] = eventdata.Data[kRealDpKin] - x_tg / ExtTarCor_DeltaCorr;

//         exttargcorr_dp += x_tg / ExtTarCor_DeltaCorr;
//         rms_exttargcorr_dp += (x_tg / ExtTarCor_DeltaCorr)*(x_tg / ExtTarCor_DeltaCorr);

//         // calculate expected dp_kin for all other excitation states
//         for (UInt_t ExcitID = 0; ExcitID < NExcitationStates; ExcitID++) {
//             assert(kRealDpKinExcitations + ExcitID < kRealTh); //check array index size
//             eventdata.Data[kRealDpKinExcitations + ExcitID] = ScatMom(ExcitationEnergyList[ExcitID], Ma, Ebeam[KineID] - ElossBeforeTg, TMath::Abs(HRSAngle)) / eventdata.Data[kCentralp] - 1;
//         }

//         DEBUG_MASSINFO("PrepareDp", "ScatterAngle = %f,\t DpKinOffsets = %f", ScatteringAngle / TMath::Pi()*180, DpKinOffsets);
//     }

//     DEBUG_INFO("PrepareDp", "%d out of %d data is for calibration", fNCalibData, fNRawData);
//     DEBUG_INFO("PrepareDp", "Average : ScatteringAngle = %f", scatang / fNRawData / TMath::Pi()*180);
//     DEBUG_INFO("PrepareDp", "Average : Travel Length before scattering = %f", D1_sum / fNRawData);
//     DEBUG_INFO("PrepareDp", "Average : Travel Length after scattering = %f", D2_sum / fNRawData);
//     DEBUG_INFO("PrepareDp", "Average : ElossTg = %f", ElossTg_sum / fNRawData);
//     DEBUG_INFO("PrepareDp", "Average : ElossAfter = %f", ElossAfter_sum / fNRawData);
//     DEBUG_INFO("PrepareDp", "Average DpKinOffsets = %f, RMS DpKinOffsets = %f", dpkinoff / fNRawData, TMath::Sqrt(dpkinoff_rms / fNRawData));
//     DEBUG_INFO("PrepareDp", "Average Extended Target Corretion: dp = %f,\t rms_dp = %f", exttargcorr_dp / fNRawData, TMath::Sqrt(rms_exttargcorr_dp / fNRawData));

//     // make sure kCalcTh, kCalcPh is filled, although not necessary
//     SumSquareDTh();
//     SumSquareDPhi();
//     SumSquareDp();

//     DEBUG_INFO("PrepareDp", "Done!");
// }

// TCanvas * LOpticsOpt::CheckDp()
// {
//     // Visualize 1D hitogram of dp_kin

//     DEBUG_INFO("CheckDp", "Entry Point");

//     // calculate Data[kCalcDpKin] for all events
//     SumSquareDp(kTRUE);

//     const Double_t DpRange = .04;
//     const UInt_t NDpRange = 500;

//     TH1D * hDpKinCalib[NKine];
//     TH1D * hDpKinAll[NKine];
//     TH1D * hRealReactZ[NKine];
//     TH1D * hElossTgBefore[NKine];
//     TH2D * hElossTg_ReactZ_Before[NKine];
//     TH2D * hElossTg_ReactZ_After[NKine];
//     TH2D * hTravelLength_ReactZ_Before[NKine];
//     TH2D * hTravelLength_ReactZ_After[NKine];
//     Double_t RealDpKin[NKine] = {0};
//     Double_t AverCalcDpKin[NKine] = {0};
//     UInt_t NEvntDpKin[NKine] = {0};
//     Double_t RealDpKinAllExcit[NExcitationStates][NKine] = {
//         {0}
//     };
//     Double_t NewArbitaryDpKinShift[NKine];
//     Double_t AveRealDpKinMatrix[NKine]={0};

//     for (UInt_t KineID = 0; KineID < NKine; KineID++) {
//       hDpKinCalib[KineID] = new TH1D(Form("hDpKinCalib%d", KineID), Form("Dp_Kin for Delta Scan Kine. %d%%", 2*(KineID-2)), NDpRange, -2*DpRange, 2*DpRange);
//       hDpKinAll[KineID] = new TH1D(Form("hDpKinAll%d", KineID), Form("Dp_Kin for Delta Scan Kine. #%d", KineID), NDpRange, -DpRange, DpRange);

//       hRealReactZ[KineID] = new TH1D(Form("hRealReactZ%d", KineID), Form("ReacZ for Delta Scan Kine. %d%%", 2*(KineID-2)), 400, -0.15, 0.15);
//       hElossTgBefore[KineID] = new TH1D(Form("hElossTg%d", KineID), Form("ElossTg for Delta Scan Kine. %d%%", 2*(KineID-2)), 400, 0, 0.0055);
//       hElossTg_ReactZ_Before[KineID] = new TH2D(Form("hElossTg_ReacZ%d", KineID), Form("ElossTg for Delta Scan Kine. %d%%", 2*(KineID-2)), 400,-0.15,0.15,400, 0, 7);
//       hElossTg_ReactZ_After[KineID] = new TH2D(Form("hElossTg_ReacZ%d", KineID), Form("ElossTg for Delta Scan Kine. %d%%", 2*(KineID-2)), 400,-0.15,0.15,400, 0, 3.5);
//       hTravelLength_ReactZ_Before[KineID] = new TH2D(Form("hTravelLength_ReacZ%d", KineID), Form("TravelLength vs ReacZ Before for Delta Scan Kine. %d%%", 2*(KineID-2)), 400,-0.15,0.15,800, 0, 0.15);
//       hTravelLength_ReactZ_After[KineID] = new TH2D(Form("hElossTg_ReacZ%d", KineID), Form("TravelLength vs ReacZ After for Delta Scan Kine. %d%%", 2*(KineID-2)), 400,-0.15,0.15,800, 0, 0.07);

//         assert(hDpKinCalib[KineID]); //pointer check
//         assert(hDpKinAll[KineID]); //pointer check
//     }

//     for (UInt_t idx = 0; idx < fNRawData; idx++) {
//         const EventData &eventdata = fRawData[idx];
//         const UInt_t ExtraDataFlag = (UInt_t) (eventdata.Data[kExtraDataFlag]);
//         assert(ExtraDataFlag == 0 || ExtraDataFlag == 1); //flag definition consistency check
//         const UInt_t KineID = (UInt_t) (eventdata.Data[kKineID]);

//         assert(KineID < NKine);
//         if (!ExtraDataFlag) {
// 	  AveRealDpKinMatrix[KineID] += eventdata.Data[kRealDpKinMatrix];
// 	  hDpKinCalib[KineID]->Fill((eventdata.Data[kCalcDpKinMatrix]));

// 	  hRealReactZ[KineID] ->Fill( eventdata.Data[kRealReactZ]);
// 	  hElossTgBefore[KineID]->Fill( eventdata.Data[kElossTgBefore]);
// 	  hElossTg_ReactZ_Before[KineID]->Fill( eventdata.Data[kRealReactZ], eventdata.Data[kElossTgBefore]*1.e3);
// 	  hElossTg_ReactZ_After[KineID]->Fill( eventdata.Data[kRealReactZ], eventdata.Data[kElossTgAfter]*1.e3);
// 	  hTravelLength_ReactZ_Before[KineID]->Fill( eventdata.Data[kRealReactZ], eventdata.Data[kTravelLengthBefore]);
// 	  hTravelLength_ReactZ_After[KineID]->Fill( eventdata.Data[kRealReactZ], eventdata.Data[kTravelLengthAfter]);

// 	  AverCalcDpKin[KineID] += eventdata.Data[kCalcDpKinMatrix];
// 	  NEvntDpKin[KineID]++;
//         }
// 	//        hDpKinAll[KineID]->Fill(eventdata.Data[kCalcDpKin] + eventdata.Data[kElossDp]);
// 	//        hDpKinAll[KineID]->Fill(eventdata.Data[kCalcDpKin]);

// 	//        RealDpKin[KineID] = eventdata.Data[kRealDpKin] + eventdata.Data[kElossDp];
//         RealDpKin[KineID] = eventdata.Data[kRealDpKinMatrix];
//         for (UInt_t ExcitID = 0; ExcitID < NExcitationStates; ExcitID++) {
//             assert(kRealDpKinExcitations + ExcitID < kRealTh); //index check
//             RealDpKinAllExcit[ExcitID][KineID] = eventdata.Data[kRealDpKinExcitations + ExcitID];
//        }
//     }
//     for (UInt_t KineID = 0; KineID < NKine; KineID++) {
//       AveRealDpKinMatrix[KineID] =  AveRealDpKinMatrix[KineID]/ NEvntDpKin[KineID];
//     }


//     TCanvas * c1 = new TCanvas("CheckDp", "Check Dp Kin Reconstruction", 1800, 900);
//     c1->Divide(3, 2);

//     for (UInt_t KineID = 0; KineID < NKine; KineID++) {
//         c1->cd(KineID+1);

//         AverCalcDpKin[KineID] /= NEvntDpKin[KineID];
//         DEBUG_MASSINFO("CheckDp", "AverCalcDpKin[%d] = %f", KineID, AverCalcDpKin[KineID]);

//         // Histograms
// 	hDpKinCalib[KineID]->SetXTitle("Dp_Kin");
// 	hDpKinCalib[KineID]->GetXaxis()->SetLabelSize(0.028);
// 	hDpKinCalib[KineID]->GetXaxis()->SetRangeUser(AverCalcDpKin[KineID]-0.01
// ,  AverCalcDpKin[KineID] +0.01 );
// 	hDpKinCalib[KineID]->Draw();
	
//         // expectation lines
//         Double_t MaxPlot = 20000;
// 	c1->cd(KineID+1)->Update();
// 	MaxPlot = c1->cd(KineID+1)->GetUymax();

// 	TLine *l = new TLine(AveRealDpKinMatrix[KineID], 0,AveRealDpKinMatrix[KineID] , MaxPlot);
// 	//	TLine *l = new TLine(RealDpKin[KineID], 0, RealDpKin[KineID], MaxPlot);
// 	cout<<"KineID:"<<KineID<<endl;
// 	cout<<"AveRealDpKinMatrix[KineID]"<<AveRealDpKinMatrix[KineID]<<endl;
//         l->SetLineColor(6);
//         l->SetLineWidth(2);
//         l->Draw("");

//         // Fits
//         const Double_t DefResolution = 1e-4;
//         const Double_t FitRangeMultiply = 20;
			
//         TString FitFunc = Form("DpPeak%d", KineID);
//         TF1 *f = new TF1(FitFunc, "gaus", AverCalcDpKin[KineID] - DefResolution*FitRangeMultiply, AverCalcDpKin[KineID] + DefResolution * FitRangeMultiply);
//         f->SetParameter(1, AverCalcDpKin[KineID]);
// 	//        f->SetParameter(2, DefResolution);
// 	//        hDpKinAll[KineID] -> Fit(FitFunc, "RN0");
//         hDpKinCalib[KineID] -> Fit(FitFunc, "R");
//         f->SetLineColor(2);
// 	f->Draw("SAME");
// 	l->Draw();

//         TLatex *t1 = new TLatex(f->GetParameter(1) + 2*f->GetParameter(2), f->GetParameter(0), Form("\\Delta = %2.1f \\times 10^{-4}", 10000 * (f->GetParameter(1) - AveRealDpKinMatrix[KineID]) ));
//         TLatex *t2 = new TLatex(f->GetParameter(1) + 2*f->GetParameter(2), 0.9*f->GetParameter(0), Form("\\sigma =  %2.1f  \\times 10^{-4}", 10000 * f->GetParameter(2)));
//         t1->SetTextSize(0.055);
//         t1->SetTextAlign(12);
//         t1->SetTextColor(2);
//         t1->Draw("same");
//         t2->SetTextSize(0.055);
//         t2->SetTextAlign(12);
//         t2->SetTextColor(2);
//         t2->Draw("same");

//         NewArbitaryDpKinShift[KineID] = f->GetParameter(1) - AveRealDpKinMatrix[KineID] + fArbitaryDpKinShift[KineID];
			
//     }

//     c1->cd(6);
//     TPaveText *t0 = new TPaveText(0.03,0.05,0.97,0.95,"NDC");
//     t0->SetShadowColor(0);
//     t0->AddText("DpKin_{Real} = #frac{P_{#theta_{HRS}} - P_{Central}}{P_{Central}}");
//     t0->AddText("DpKin = dp - #frac{(P_{#theta} - P_{Loss}) - P_{#theta_{HRS}} }{ P_{Central} }");
//     t0->AddText("DpKin - ti_{Real} = dp - #frac{(P_{#theta} - P_{Loss}) - P_{Central}}{P_{Central}}");
//     t0->Draw();

//     Info("CheckDp", "New set of arbitary dp shifts:");
//     for (UInt_t KineID = 0; KineID < NKine; KineID++)
//         printf("opt->fArbitaryDpKinShift[%d] = %e;\n", KineID, NewArbitaryDpKinShift[KineID]);

       
//     TCanvas * c2 = new TCanvas("Check reconstructed Z", "Check reconstructed Z", 1800, 900);
//     c2->Divide(3, 2);

//     for (UInt_t KineID = 0; KineID < NKine; KineID++) {
//         c2->cd(KineID+1);
// 	hRealReactZ[KineID]->Draw();
//     }
          
//     TCanvas * c3 = new TCanvas("Check ElossTg", "Check ElossTg Before ", 1800, 900);
//     c3->Divide(3, 2);

//     for (UInt_t KineID = 0; KineID < NKine; KineID++) {
//         c3->cd(KineID+1);
// 	hElossTgBefore[KineID]->Draw();
//     }
       
//     TCanvas * c4 = new TCanvas("c4", "Check ElossTg vs ReactZ Before", 1800, 900);
//     /* c4->Divide(3, 2);

//     for (UInt_t KineID = 0; KineID < NKine; KineID++) {
//         c4->cd(KineID+1);
// 	hElossTg_ReactZ_Before[KineID]->Draw();
// 	}*/
//     hElossTg_ReactZ_Before[2]->Draw();
//     TCanvas * c5 = new TCanvas("c5", "Check ElossTg vs ReactZ After", 1800, 900);
//     /*    c5->Divide(3, 2);

//     for (UInt_t KineID = 0; KineID < NKine; KineID++) {
//         c5->cd(KineID+1);
// 	hElossTg_ReactZ_After[KineID]->Draw();
// 	}*/
//     hElossTg_ReactZ_After[2]->Draw();
    
//     TCanvas * c6 = new TCanvas("c6", "Check ElossTg vs ReactZ Before", 1800, 900);
//     c6->Divide(3, 2);

//     for (UInt_t KineID = 0; KineID < NKine; KineID++) {
//         c6->cd(KineID+1);
// 	hTravelLength_ReactZ_Before[KineID]->Draw();
//     }
//     TCanvas * c7 = new TCanvas("c7", "Check ElossTg vs ReactZ After", 1800, 900);
//     c7->Divide(3, 2);

//     for (UInt_t KineID = 0; KineID < NKine; KineID++) {
//         c7->cd(KineID+1);
// 	hTravelLength_ReactZ_After[KineID]->Draw();
//     }
        
//     return c1;
// }
// 

// TCanvas * LOpticsOpt::CheckDpGlobal()
// {
//     // Visualize 1D hitogram of dp_kin

//     DEBUG_INFO("CheckDp", "Entry Point");

//     SumSquareDp(kTRUE);

//     const Double_t DpRange = .05;
//     const UInt_t NDpRange = 5000;

//     TH1D * hDpKinCalib[NKine];
//     TH1D * hDpKinAll[NKine];
//     Double_t RealDpKin[NKine];
//     Double_t AverCalcDpKin[NKine] = {0};
//     UInt_t NEvntDpKin[NKine] = {0};
//     Double_t RealDpKinAllExcit[NExcitationStates][NKine] = {
//         {0}
//     };

//     for (UInt_t KineID = 0; KineID < NKine; KineID++) {
//         hDpKinCalib[KineID] = new TH1D(Form("hDpKinCalib%d", KineID), "Dp_Kin for Delta Scan", NDpRange, -DpRange, DpRange);
//         hDpKinAll[KineID] = new TH1D(Form("hDpKinAll%d", KineID), Form("Dp_Kin for Delta Scan Kine. #%d", KineID), NDpRange, -DpRange, DpRange);

//         assert(hDpKinCalib[KineID]); //pointer check
//         assert(hDpKinAll[KineID]); //pointer check
//     }

//     for (UInt_t idx = 0; idx < fNRawData; idx++) {
//         const EventData &eventdata = fRawData[idx];
//         const UInt_t ExtraDataFlag = (UInt_t) (eventdata.Data[kExtraDataFlag]);
//         assert(ExtraDataFlag == 0 || ExtraDataFlag == 1); //flag definition consistency check
//         const UInt_t KineID = (UInt_t) (eventdata.Data[kKineID]);
//         assert(KineID < NKine);

//         if (!ExtraDataFlag) {
//             hDpKinCalib[KineID]->Fill(eventdata.Data[kCalcDpKin] + eventdata.Data[kElossDp]);
//             AverCalcDpKin[KineID] += eventdata.Data[kCalcDpKin] + eventdata.Data[kElossDp];
//             NEvntDpKin[KineID]++;
//         }
//         hDpKinAll[KineID]->Fill(eventdata.Data[kCalcDpKin] + eventdata.Data[kElossDp]);

// 	//        RealDpKin[KineID] = eventdata.Data[kRealDpKin] + eventdata.Data[kElossDp];
//         RealDpKin[KineID] = eventdata.Data[kRealDpKin];

//         for (UInt_t ExcitID = 0; ExcitID < NExcitationStates; ExcitID++) {
//             assert(kRealDpKinExcitations + ExcitID < kRealTh); //index check
//             RealDpKinAllExcit[ExcitID][KineID] = eventdata.Data[kRealDpKinExcitations + ExcitID];
//         }
//     }

//     TCanvas * c1 = new TCanvas("CheckDp", "Check Dp Kin Reconstruction", 1800, 900);

//     for (UInt_t KineID = 0; KineID < NKine; KineID++) {
//         c1->cd();
//         // gPad -> SetLogy();

//         AverCalcDpKin[KineID] /= NEvntDpKin[KineID];
//         DEBUG_MASSINFO("CheckDp", "AverCalcDpKin[%d] = %f", KineID, AverCalcDpKin[KineID]);

//         // Histograms
//         hDpKinCalib[KineID]->SetLineColor(4);
//         hDpKinCalib[KineID]->SetFillColor(4);
//         hDpKinCalib[KineID]->SetFillStyle(3008);

//         hDpKinAll[KineID]->SetLineColor(1);
//         // hDpKinAll[KineID]->SetFillColor(1);

//         if (KineID == 0) {
//             const Double_t dpRange = 0.04;
//             hDpKinCalib[KineID]->SetAxisRange(-dpRange, +dpRange);
//             hDpKinAll[KineID]->SetAxisRange(RealDpKin[KineID] - dpRange, RealDpKin[KineID] + dpRange);

//             hDpKinCalib[KineID]->SetXTitle("dp_kin (energy loss corrected)");

//             hDpKinCalib[KineID]->Draw();
//         } else {
//             hDpKinCalib[KineID]->Draw("same");
//         }

//         // expectation lines
//         const Double_t MaxPlot = 20000;
//         for (UInt_t ExcitID = 0; ExcitID < NExcitationStates; ExcitID++) {
//             const Double_t x = RealDpKinAllExcit[ExcitID][KineID];
//             TLine *l = new TLine(x, 0, x, +MaxPlot);
//             l->SetLineColor(3);
//             l->SetLineWidth(2);
//             l->Draw();
//         }
//         TLine *l = new TLine(RealDpKin[KineID], 0, RealDpKin[KineID], +MaxPlot);
//         l->SetLineColor(6);
//         l->SetLineWidth(2);
//         l->Draw();
//     }

//     return c1;
// }

Double_t LOpticsOpt::SumSquareDp(Bool_t IncludeExtraData)
{
    // return square sum of diff between calculated dp_kin and expected dp_kin

    Double_t d_dp = 0; //Difference
    Double_t rms_dp = 0; //mean square

    static UInt_t NCall = 0;
    NCall++;

    UInt_t NCalibData = 0;

    if (IncludeExtraData) {
        Warning("SumSquareDp", "Data Beyond selected excitation state is included in this calculation");
    }
    for (UInt_t idx = 0; idx < fNRawData; idx++) {
        Double_t dp, dp_kin;

        EventData &eventdata = fRawData[idx];

        //jump through data beyond selected excitation states
        if (eventdata.Data[kExtraDataFlag] > 0 && !IncludeExtraData) continue;
        NCalibData++;

        Double_t x_fp = eventdata.Data[kX];
        const Double_t(*powers)[5] = eventdata.powers;

        // calculate the matrices we need
        CalcMatrix(x_fp, fDMatrixElems);
        // CalcMatrix(x_fp, fTMatrixElems);
        // CalcMatrix(x_fp, fYMatrixElems);
        // CalcMatrix(x_fp, fYTAMatrixElems);
        // CalcMatrix(x_fp, fPMatrixElems);
        // CalcMatrix(x_fp, fPTAMatrixElems);

        // calculate momentum
        dp = CalcTargetVar(fDMatrixElems, powers);
	dp_kin = dp - eventdata.Data[kDpKinOffsets];
	cout<<eventdata.Data[kDpKinOffsets]<<endl;
        const UInt_t KineID = (UInt_t) (eventdata.Data[kKineID]);
        assert(KineID < NKine); //check array index size
        const Double_t ArbitaryDpKinShift = fArbitaryDpKinShift[KineID];

        d_dp += dp_kin - eventdata.Data[kRealDpKinMatrix] + ArbitaryDpKinShift;
	rms_dp += (dp_kin - eventdata.Data[kRealDpKinMatrix] + ArbitaryDpKinShift)*(dp_kin - eventdata.Data[kRealDpKinMatrix] + ArbitaryDpKinShift);

        DEBUG_MASSINFO("SumSquareDp", "d_dp = %f = \t%f - \t%f", dp_kin - eventdata.Data[kRealDpKinMatrix], dp_kin, eventdata.Data[kRealDpKinMatrix]);

        // save the results
        eventdata.Data[kCalcDpKinMatrix] = dp_kin;
        eventdata.Data[kCalcDpKin] = dp_kin + eventdata.Data[kRealTgX] / ExtTarCor_DeltaCorr;
    }

    if (!IncludeExtraData)
        assert(fNCalibData == NCalibData); // check number of event for calibration

    // DEBUG_INFO("SumSquareDp", "#%d : d_dp = %f, rms_dp = %f", NCall, d_dp / NCalibData, TMath::Sqrt(rms_dp / NCalibData));
    printf("SumSquareDp: #%d : d_dp = %f, rms_dp = %f\n", NCall, d_dp / NCalibData, 100*TMath::Sqrt(rms_dp / NCalibData));

    return rms_dp;
}

///////////////////////////////////////////////////////////////////////////////
// declarations for target vertex reconstruction
///////////////////////////////////////////////////////////////////////////////

void LOpticsOpt::CalcMatrix(const Double_t x, vector<THaMatrixElement>& matrix)
{
    // calculates the values of the matrix elements for a given location
    // by evaluating a polynomial in x of order it->order with
    // coefficients given by it->poly

    for (vector<THaMatrixElement>::iterator it = matrix.begin();
            it != matrix.end(); it++) {
        it->v = 0.0;

        if (it->order > 0) {
            for (int i = it->order - 1; i >= 1; i--)
                it->v = x * (it->v + it->poly[i]);
            it->v += it->poly[0];
        }
    }
}

Double_t LOpticsOpt::CalcTargetVar(const vector<THaMatrixElement>& matrix, const Double_t powers[][5])
{
    // calculates the value of a variable at the target
    // the x-dependence is already in the matrix, so only 1-3 (or np) used
    Double_t retval = 0.0;
    Double_t v = 0;
    for (vector<THaMatrixElement>::const_iterator it = matrix.begin();
            it != matrix.end(); it++)
        if (it->v != 0.0) {
            v = it->v;
            unsigned int np = it->pw.size(); // generalize for extra matrix elems.
            for (unsigned int i = 0; i < np; i++)
                v *= powers[it->pw[i]][i + 1];
            retval += v;
            //      retval += it->v * powers[it->pw[0]][1]
            //	              * powers[it->pw[1]][2]
            //	              * powers[it->pw[2]][3];
        }

    return retval;
}




///////////////////////////////////////////////////////////////////////////////
// class for storing matrix element data
///////////////////////////////////////////////////////////////////////////////

bool THaMatrixElement::match(const THaMatrixElement& rhs) const
{
    // Compare coefficients of this matrix element to another

    if (pw.size() != rhs.pw.size())
        return false;
    for (vector<int>::size_type i = 0; i < pw.size(); i++) {
        if (pw[i] != rhs.pw[i])
            return false;
    }
    return true;
}

void THaMatrixElement::SkimPoly()
{
    // reduce order to highest non-zero poly

    if (iszero) return;

    while (!poly[order - 1] && order > 0) {
        poly.pop_back();
        order = order - 1;
    }

    if (order == 0) iszero = kTRUE;
}

ClassImp(LOpticsOpt);
