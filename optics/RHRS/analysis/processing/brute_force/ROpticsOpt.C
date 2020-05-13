///////////////////////////////////////////////////////////////////////////////
//
// ROpticsOpt
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
#include "TChain.h"
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

#include "THaGlobals.h"
#include "THaEvData.h"
#include "THaDetMap.h"
#include "THaTrack.h"
#include "THaScintillator.h"
#include "THaSpectrometer.h"
#include "VarDef.h"

#include "ROpticsOpt.h"

#ifdef WITH_DEBUG
#include <iostream>
#endif

#define tgy_hole_info false

using namespace std;
using THaString::Split;

///////////////////////////////////////////////////////////////////////////////
// Input Sections
///////////////////////////////////////////////////////////////////////////////

#include "InputR.h"

///////////////////////////////////////////////////////////////////////////////
// Constructors
///////////////////////////////////////////////////////////////////////////////

ROpticsOpt::ROpticsOpt(const char* name, const char* description, THaApparatus* apparatus) :
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


    DEBUG_INFO("ROpticsOpt", "Read in configuration " + InputID);
    DEBUG_INFO("ROpticsOpt", "HRS @ %f Degree, PointingOffset = (%f,%f,%f), SievePos = (%f,%f,%f)", HRSAngle / TMath::Pi()*180, fPointingOffset.X(), fPointingOffset.Y(), fPointingOffset.Z(), SieveOffX, SieveOffY, ZPos);

    fNRawData = 0;

    for (UInt_t i = 0; i < 100; i++)
        fArbitaryDpKinShift[i] = fArbitaryVertexShift[i] = 0;
}

ROpticsOpt::~ROpticsOpt()
{
    // Destructor.
}



Int_t ROpticsOpt::LoadDataBase(TString DataBaseName,TString xfp_range)
{
  DataBaseName = DataBaseName + "_" + xfp_range;
  
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


    std::vector<THaMatrixElement> fTMatrixElems;
    std::vector<THaMatrixElement> fDMatrixElems;
    std::vector<THaMatrixElement> fPMatrixElems;
    std::vector<THaMatrixElement> fPTAMatrixElems;
    std::vector<THaMatrixElement> fYMatrixElems;
    std::vector<THaMatrixElement> fYTAMatrixElems;
    std::vector<THaMatrixElement> fFPMatrixElems;
    std::vector<THaMatrixElement> fRMatrixElems;

    if(xfp_range == "-50_-30"){
      fTMatrixElems = fTMatrixElems_50_30;
      fDMatrixElems = fDMatrixElems_50_30;
      fPMatrixElems = fPMatrixElems_50_30;
      fPTAMatrixElems = fPTAMatrixElems_50_30;
      fYMatrixElems = fYMatrixElems_50_30;
      fYTAMatrixElems = fYTAMatrixElems_50_30;
      fRMatrixElems = fRMatrixElems_50_30;
      fFPMatrixElems = fFPMatrixElems_50_30;
    }

    if(xfp_range == "-30_-10"){
      fTMatrixElems = fTMatrixElems_30_10;
      fDMatrixElems = fDMatrixElems_30_10;
      fPMatrixElems = fPMatrixElems_30_10;
      fPTAMatrixElems = fPTAMatrixElems_30_10;
      fYMatrixElems = fYMatrixElems_30_10;
      fYTAMatrixElems = fYTAMatrixElems_30_10;
      fRMatrixElems = fRMatrixElems_30_10;
      fFPMatrixElems = fFPMatrixElems_30_10;
    }

    if(xfp_range == "-10_10"){
      fTMatrixElems = fTMatrixElems_10_10;
      fDMatrixElems = fDMatrixElems_10_10;
      fPMatrixElems = fPMatrixElems_10_10;
      fPTAMatrixElems = fPTAMatrixElems_10_10;
      fYMatrixElems = fYMatrixElems_10_10;
      fYTAMatrixElems = fYTAMatrixElems_10_10;
      fRMatrixElems = fRMatrixElems_10_10;
      fFPMatrixElems = fFPMatrixElems_10_10;
    }

    if(xfp_range == "10_30"){
      fTMatrixElems = fTMatrixElems_10_30;
      fDMatrixElems = fDMatrixElems_10_30;
      fPMatrixElems = fPMatrixElems_10_30;
      fPTAMatrixElems = fPTAMatrixElems_10_30;
      fYMatrixElems = fYMatrixElems_10_30;
      fYTAMatrixElems = fYTAMatrixElems_10_30;
      fRMatrixElems = fRMatrixElems_10_30;
      fFPMatrixElems = fFPMatrixElems_10_30;
    }

    if(xfp_range == "30_50"){
      fTMatrixElems = fTMatrixElems_30_50;
      fDMatrixElems = fDMatrixElems_30_50;
      fPMatrixElems = fPMatrixElems_30_50;
      fPTAMatrixElems = fPTAMatrixElems_30_50;
      fYMatrixElems = fYMatrixElems_30_50;
      fYTAMatrixElems = fYTAMatrixElems_30_50;
      fRMatrixElems = fRMatrixElems_30_50;
      fFPMatrixElems = fFPMatrixElems_30_50;
    }

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



    if(xfp_range == "-50_-30"){
      fTMatrixElems_50_30 = fTMatrixElems;
      fDMatrixElems_50_30 = fDMatrixElems;
      fPMatrixElems_50_30 = fPMatrixElems;
      fPTAMatrixElems_50_30 = fPTAMatrixElems;
      fYMatrixElems_50_30 = fYMatrixElems;
      fYTAMatrixElems_50_30 = fYTAMatrixElems;
      fRMatrixElems_50_30 = fRMatrixElems;
      fFPMatrixElems_50_30 = fFPMatrixElems;
    }

     if(xfp_range == "-30_-10"){
      fTMatrixElems_30_10 = fTMatrixElems;
      fDMatrixElems_30_10 = fDMatrixElems;
      fPMatrixElems_30_10 = fPMatrixElems;
      fPTAMatrixElems_30_10 = fPTAMatrixElems;
      fYMatrixElems_30_10 = fYMatrixElems;
      fYTAMatrixElems_30_10 = fYTAMatrixElems;
      fRMatrixElems_30_10 = fRMatrixElems;
      fFPMatrixElems_30_10 = fFPMatrixElems;
    }

     if(xfp_range == "-10_10"){
      fTMatrixElems_10_10 = fTMatrixElems;
      fDMatrixElems_10_10 = fDMatrixElems;
      fPMatrixElems_10_10 = fPMatrixElems;
      fPTAMatrixElems_10_10 = fPTAMatrixElems;
      fYMatrixElems_10_10 = fYMatrixElems;
      fYTAMatrixElems_10_10 = fYTAMatrixElems;
      fRMatrixElems_10_10 = fRMatrixElems;
      fFPMatrixElems_10_10 = fFPMatrixElems;
     }

     if(xfp_range == "10_30"){
      fTMatrixElems_10_30 = fTMatrixElems;
      fDMatrixElems_10_30 = fDMatrixElems;
      fPMatrixElems_10_30 = fPMatrixElems;
      fPTAMatrixElems_10_30 = fPTAMatrixElems;
      fYMatrixElems_10_30 = fYMatrixElems;
      fYTAMatrixElems_10_30 = fYTAMatrixElems;
      fRMatrixElems_10_30 = fRMatrixElems;
      fFPMatrixElems_10_30 = fFPMatrixElems;
    }

     if(xfp_range == "30_50"){
      fTMatrixElems_30_50 = fTMatrixElems;
      fDMatrixElems_30_50 = fDMatrixElems;
      fPMatrixElems_30_50 = fPMatrixElems;
      fPTAMatrixElems_30_50 = fPTAMatrixElems;
      fYMatrixElems_30_50 = fYMatrixElems;
      fYTAMatrixElems_30_50 = fYTAMatrixElems;
      fRMatrixElems_30_50 = fRMatrixElems;
      fFPMatrixElems_30_50 = fFPMatrixElems;
    }


   
    fIsInit = true;
    fclose(file);
    return kOK;
}




UInt_t ROpticsOpt::LoadRawData(TChain* t)
{

  int entries = t->GetEntries();

  double R_tr_x_rot[100];
  double R_tr_y_rot[100];
  double R_tr_th_rot[100];
  double R_tr_ph_rot[100];
  double R_x[100];
  double R_y[100];
  double Ru_x[100];
  double Ru_y[100];
  
  t->SetBranchStatus("*",0);
 
  t->SetBranchStatus("R.tr.r_x",1);
  t->SetBranchStatus("R.tr.r_y",1);
  t->SetBranchStatus("R.tr.r_th",1);
  t->SetBranchStatus("R.tr.r_ph",1);
  t->SetBranchStatus("Rrb.x",1);
  t->SetBranchStatus("Rrb.y",1);
  t->SetBranchStatus("Rurb.x",1);
  t->SetBranchStatus("Rurb.y",1);

  t->SetBranchAddress("R.tr.r_x",R_tr_x_rot);
  t->SetBranchAddress("R.tr.r_y",R_tr_y_rot);
  t->SetBranchAddress("R.tr.r_th",R_tr_th_rot);
  t->SetBranchAddress("R.tr.r_ph",R_tr_ph_rot);
  t->SetBranchAddress("Rrb.x",R_x);
  t->SetBranchAddress("Rrb.y",R_y);
  t->SetBranchAddress("Rurb.x",Ru_x);
  t->SetBranchAddress("Rurb.y",Ru_y);

  for(UInt_t NRead = 0; NRead<entries; NRead++){
    Double_t * eventdata = fRawData[NRead].Data;

    t->GetEntry(NRead);

        
    eventdata[0] = R_tr_x_rot[0];
    eventdata[1] = R_tr_th_rot[0];
    eventdata[2] = R_tr_y_rot[0];
    eventdata[3] = R_tr_ph_rot[0];
    eventdata[4] = R_x[0];
    eventdata[5] = R_y[0];
  

 
    Double_t(*powers)[5] = fRawData[NRead].powers;
    Double_t x_fp = eventdata[kX];
    Double_t th_fp = eventdata[kTh];
    Double_t y_fp = eventdata[kY];
    Double_t ph_fp = eventdata[kPhi];

    
    // calculate the powers we need
    for (int i = 0; i < kNUM_PRECOMP_POW; i++) {
      powers[i][0] = pow(x_fp, i);
      powers[i][1] = pow(th_fp, i);
      powers[i][2] = pow(y_fp, i);
      powers[i][3] = pow(ph_fp, i);
      powers[i][4] = pow(TMath::Abs(th_fp), i);
    }
    
    
  }
  
    return entries;
}

///////////////////////////////////////////////////////////////////////////////
// declarations for target vertex reconstruction
///////////////////////////////////////////////////////////////////////////////

void ROpticsOpt::CalcMatrix(const Double_t x, vector<THaMatrixElement>& matrix)
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


 

Double_t ROpticsOpt::CalcTargetVar(const vector<THaMatrixElement>& matrix, const Double_t powers[][5])
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

double ROpticsOpt::calc_tgy(int event){

  std::vector<THaMatrixElement> fTMatrixElems;
  std::vector<THaMatrixElement> fDMatrixElems;
  std::vector<THaMatrixElement> fPMatrixElems;
  std::vector<THaMatrixElement> fPTAMatrixElems;
  std::vector<THaMatrixElement> fYMatrixElems;
  std::vector<THaMatrixElement> fYTAMatrixElems;
  std::vector<THaMatrixElement> fFPMatrixElems;
  std::vector<THaMatrixElement> fRMatrixElems;
  

  double tg_y;
  
  EventData &eventdata = fRawData[event];
  
  Double_t x_fp = eventdata.Data[kX];
  const Double_t(*powers)[5] = eventdata.powers;

  if(x_fp > -0.50 && x_fp < -0.30){
    fTMatrixElems = fTMatrixElems_50_30;
    fDMatrixElems = fDMatrixElems_50_30;
    fPMatrixElems = fPMatrixElems_50_30;
    fPTAMatrixElems = fPTAMatrixElems_50_30;
    fYMatrixElems = fYMatrixElems_50_30;
    fYTAMatrixElems = fYTAMatrixElems_50_30;
    fRMatrixElems = fRMatrixElems_50_30;
    fFPMatrixElems = fFPMatrixElems_50_30;
  }
  else if(x_fp > -0.30 && x_fp < -0.10){
    fTMatrixElems = fTMatrixElems_30_10;
    fDMatrixElems = fDMatrixElems_30_10;
    fPMatrixElems = fPMatrixElems_30_10;
    fPTAMatrixElems = fPTAMatrixElems_30_10;
    fYMatrixElems = fYMatrixElems_30_10;
    fYTAMatrixElems = fYTAMatrixElems_30_10;
    fRMatrixElems = fRMatrixElems_30_10;
    fFPMatrixElems = fFPMatrixElems_30_10;
  }
  else if(x_fp > -0.10 && x_fp < 0.10){
    fTMatrixElems = fTMatrixElems_10_10;
    fDMatrixElems = fDMatrixElems_10_10;
    fPMatrixElems = fPMatrixElems_10_10;
    fPTAMatrixElems = fPTAMatrixElems_10_10;
    fYMatrixElems = fYMatrixElems_10_10;
    fYTAMatrixElems = fYTAMatrixElems_10_10;
    fRMatrixElems = fRMatrixElems_10_10;
    fFPMatrixElems = fFPMatrixElems_10_10;
  }
  else if(x_fp > 0.10 && x_fp < 0.30){
    fTMatrixElems = fTMatrixElems_10_30;
    fDMatrixElems = fDMatrixElems_10_30;
    fPMatrixElems = fPMatrixElems_10_30;
    fPTAMatrixElems = fPTAMatrixElems_10_30;
    fYMatrixElems = fYMatrixElems_10_30;
    fYTAMatrixElems = fYTAMatrixElems_10_30;
    fRMatrixElems = fRMatrixElems_10_30;
    fFPMatrixElems = fFPMatrixElems_10_30;
  }
  else {
      fTMatrixElems = fTMatrixElems_30_50;
      fDMatrixElems = fDMatrixElems_30_50;
      fPMatrixElems = fPMatrixElems_30_50;
      fPTAMatrixElems = fPTAMatrixElems_30_50;
      fYMatrixElems = fYMatrixElems_30_50;
      fYTAMatrixElems = fYTAMatrixElems_30_50;
      fRMatrixElems = fRMatrixElems_30_50;
      fFPMatrixElems = fFPMatrixElems_30_50;
    }
  
  // calculate the matrices we need
  // CalcMatrix(x_fp, fDMatrixElems);
  // CalcMatrix(x_fp, fTMatrixElems);
  CalcMatrix(x_fp, fYMatrixElems);
  CalcMatrix(x_fp, fYTAMatrixElems);
  // CalcMatrix(x_fp, fPMatrixElems);
  // CalcMatrix(x_fp, fPTAMatrixElems);
  
  // calculate the coordinates at the target
  tg_y = CalcTargetVar(fYMatrixElems, powers) + CalcTargetVar(fYTAMatrixElems, powers);

  //if(tg_y > 0 && tg_y < 0.0002) cout<<tg_y*1000<<" "<<x_fp*1000<<endl;
  
  return tg_y;
  
}

double ROpticsOpt::calc_tgth(int event){

  std::vector<THaMatrixElement> fTMatrixElems;
  std::vector<THaMatrixElement> fDMatrixElems;
  std::vector<THaMatrixElement> fPMatrixElems;
  std::vector<THaMatrixElement> fPTAMatrixElems;
  std::vector<THaMatrixElement> fYMatrixElems;
  std::vector<THaMatrixElement> fYTAMatrixElems;
  std::vector<THaMatrixElement> fFPMatrixElems;
  std::vector<THaMatrixElement> fRMatrixElems;
 

  double theta;
  
  
  EventData &eventdata = fRawData[event];
  
  Double_t x_fp = eventdata.Data[kX];
  const Double_t(*powers)[5] = eventdata.powers;

  if(x_fp > -0.50 && x_fp < -0.30){
    fTMatrixElems = fTMatrixElems_50_30;
    fDMatrixElems = fDMatrixElems_50_30;
    fPMatrixElems = fPMatrixElems_50_30;
    fPTAMatrixElems = fPTAMatrixElems_50_30;
    fYMatrixElems = fYMatrixElems_50_30;
    fYTAMatrixElems = fYTAMatrixElems_50_30;
    fRMatrixElems = fRMatrixElems_50_30;
    fFPMatrixElems = fFPMatrixElems_50_30;
  }
  else if(x_fp > -0.30 && x_fp < -0.10){
    fTMatrixElems = fTMatrixElems_30_10;
    fDMatrixElems = fDMatrixElems_30_10;
    fPMatrixElems = fPMatrixElems_30_10;
    fPTAMatrixElems = fPTAMatrixElems_30_10;
    fYMatrixElems = fYMatrixElems_30_10;
    fYTAMatrixElems = fYTAMatrixElems_30_10;
    fRMatrixElems = fRMatrixElems_30_10;
    fFPMatrixElems = fFPMatrixElems_30_10;
  }
  else if(x_fp > -0.10 && x_fp < 0.10){
    fTMatrixElems = fTMatrixElems_10_10;
    fDMatrixElems = fDMatrixElems_10_10;
    fPMatrixElems = fPMatrixElems_10_10;
    fPTAMatrixElems = fPTAMatrixElems_10_10;
    fYMatrixElems = fYMatrixElems_10_10;
    fYTAMatrixElems = fYTAMatrixElems_10_10;
    fRMatrixElems = fRMatrixElems_10_10;
    fFPMatrixElems = fFPMatrixElems_10_10;
  }
  else if(x_fp > 0.10 && x_fp < 0.30){
    fTMatrixElems = fTMatrixElems_10_30;
    fDMatrixElems = fDMatrixElems_10_30;
    fPMatrixElems = fPMatrixElems_10_30;
    fPTAMatrixElems = fPTAMatrixElems_10_30;
    fYMatrixElems = fYMatrixElems_10_30;
    fYTAMatrixElems = fYTAMatrixElems_10_30;
    fRMatrixElems = fRMatrixElems_10_30;
    fFPMatrixElems = fFPMatrixElems_10_30;
  }
  else{
    fTMatrixElems = fTMatrixElems_30_50;
    fDMatrixElems = fDMatrixElems_30_50;
    fPMatrixElems = fPMatrixElems_30_50;
    fPTAMatrixElems = fPTAMatrixElems_30_50;
    fYMatrixElems = fYMatrixElems_30_50;
    fYTAMatrixElems = fYTAMatrixElems_30_50;
    fRMatrixElems = fRMatrixElems_30_50;
    fFPMatrixElems = fFPMatrixElems_30_50;
  }

  
  // calculate the matrices we need
  // CalcMatrix(x_fp, fDMatrixElems);
  CalcMatrix(x_fp, fTMatrixElems);
  // CalcMatrix(x_fp, fYMatrixElems);
  // CalcMatrix(x_fp, fYTAMatrixElems);
  // CalcMatrix(x_fp, fPMatrixElems);
  // CalcMatrix(x_fp, fPTAMatrixElems);
  
  // calculate the coordinates at the target
  theta = CalcTargetVar(fTMatrixElems, powers);
  
  return theta; 
}

double ROpticsOpt::calc_tgph(int event){

  std::vector<THaMatrixElement> fTMatrixElems;
  std::vector<THaMatrixElement> fDMatrixElems;
  std::vector<THaMatrixElement> fPMatrixElems;
  std::vector<THaMatrixElement> fPTAMatrixElems;
  std::vector<THaMatrixElement> fYMatrixElems;
  std::vector<THaMatrixElement> fYTAMatrixElems;
  std::vector<THaMatrixElement> fFPMatrixElems;
  std::vector<THaMatrixElement> fRMatrixElems;
  

  double phi;
  
  
  EventData &eventdata = fRawData[event];
  
  Double_t x_fp = eventdata.Data[kX];
  const Double_t(*powers)[5] = eventdata.powers;

  if(x_fp > -0.50 && x_fp < -0.30){
    fTMatrixElems = fTMatrixElems_50_30;
    fDMatrixElems = fDMatrixElems_50_30;
    fPMatrixElems = fPMatrixElems_50_30;
    fPTAMatrixElems = fPTAMatrixElems_50_30;
    fYMatrixElems = fYMatrixElems_50_30;
    fYTAMatrixElems = fYTAMatrixElems_50_30;
    fRMatrixElems = fRMatrixElems_50_30;
    fFPMatrixElems = fFPMatrixElems_50_30;
  }
  else if(x_fp > -0.30 && x_fp < -0.10){
    fTMatrixElems = fTMatrixElems_30_10;
    fDMatrixElems = fDMatrixElems_30_10;
    fPMatrixElems = fPMatrixElems_30_10;
    fPTAMatrixElems = fPTAMatrixElems_30_10;
    fYMatrixElems = fYMatrixElems_30_10;
    fYTAMatrixElems = fYTAMatrixElems_30_10;
    fRMatrixElems = fRMatrixElems_30_10;
    fFPMatrixElems = fFPMatrixElems_30_10;
  }

  if(x_fp > -0.10 && x_fp < 0.10){
    fTMatrixElems = fTMatrixElems_10_10;
    fDMatrixElems = fDMatrixElems_10_10;
    fPMatrixElems = fPMatrixElems_10_10;
    fPTAMatrixElems = fPTAMatrixElems_10_10;
    fYMatrixElems = fYMatrixElems_10_10;
    fYTAMatrixElems = fYTAMatrixElems_10_10;
    fRMatrixElems = fRMatrixElems_10_10;
    fFPMatrixElems = fFPMatrixElems_10_10;
  }
  else if(x_fp > 0.10 && x_fp < 0.30){
    fTMatrixElems = fTMatrixElems_10_30;
    fDMatrixElems = fDMatrixElems_10_30;
    fPMatrixElems = fPMatrixElems_10_30;
    fPTAMatrixElems = fPTAMatrixElems_10_30;
    fYMatrixElems = fYMatrixElems_10_30;
    fYTAMatrixElems = fYTAMatrixElems_10_30;
    fRMatrixElems = fRMatrixElems_10_30;
    fFPMatrixElems = fFPMatrixElems_10_30;
  }
  else {
    fTMatrixElems = fTMatrixElems_30_50;
    fDMatrixElems = fDMatrixElems_30_50;
    fPMatrixElems = fPMatrixElems_30_50;
    fPTAMatrixElems = fPTAMatrixElems_30_50;
    fYMatrixElems = fYMatrixElems_30_50;
    fYTAMatrixElems = fYTAMatrixElems_30_50;
    fRMatrixElems = fRMatrixElems_30_50;
    fFPMatrixElems = fFPMatrixElems_30_50;
  }

  // calculate the matrices we need
  // CalcMatrix(x_fp, fDMatrixElems);
  //CalcMatrix(x_fp, fTMatrixElems);
  // CalcMatrix(x_fp, fYMatrixElems);
  // CalcMatrix(x_fp, fYTAMatrixElems);
  CalcMatrix(x_fp, fPMatrixElems);
  CalcMatrix(x_fp, fPTAMatrixElems);
    
  // calculate the coordinates at the target
  phi = CalcTargetVar(fPMatrixElems, powers) + CalcTargetVar(fPTAMatrixElems, powers);
  
  
  return phi; 

}


double ROpticsOpt::calc_tgdp(int event){

  std::vector<THaMatrixElement> fTMatrixElems;
  std::vector<THaMatrixElement> fDMatrixElems;
  std::vector<THaMatrixElement> fPMatrixElems;
  std::vector<THaMatrixElement> fPTAMatrixElems;
  std::vector<THaMatrixElement> fYMatrixElems;
  std::vector<THaMatrixElement> fYTAMatrixElems;
  std::vector<THaMatrixElement> fFPMatrixElems;
  std::vector<THaMatrixElement> fRMatrixElems;
  
  double dp;
  
  
  EventData &eventdata = fRawData[event];
  
  Double_t x_fp = eventdata.Data[kX];
  const Double_t(*powers)[5] = eventdata.powers;

  if(x_fp > -0.50 && x_fp < -0.30){
    fTMatrixElems = fTMatrixElems_50_30;
    fDMatrixElems = fDMatrixElems_50_30;
    fPMatrixElems = fPMatrixElems_50_30;
    fPTAMatrixElems = fPTAMatrixElems_50_30;
    fYMatrixElems = fYMatrixElems_50_30;
    fYTAMatrixElems = fYTAMatrixElems_50_30;
    fRMatrixElems = fRMatrixElems_50_30;
    fFPMatrixElems = fFPMatrixElems_50_30;
  }
  else if(x_fp > -0.30 && x_fp < -0.10){
    fTMatrixElems = fTMatrixElems_30_10;
    fDMatrixElems = fDMatrixElems_30_10;
    fPMatrixElems = fPMatrixElems_30_10;
    fPTAMatrixElems = fPTAMatrixElems_30_10;
    fYMatrixElems = fYMatrixElems_30_10;
    fYTAMatrixElems = fYTAMatrixElems_30_10;
    fRMatrixElems = fRMatrixElems_30_10;
    fFPMatrixElems = fFPMatrixElems_30_10;
  }
  else if(x_fp > -0.10 && x_fp < 0.10){
    fTMatrixElems = fTMatrixElems_10_10;
    fDMatrixElems = fDMatrixElems_10_10;
    fPMatrixElems = fPMatrixElems_10_10;
    fPTAMatrixElems = fPTAMatrixElems_10_10;
    fYMatrixElems = fYMatrixElems_10_10;
    fYTAMatrixElems = fYTAMatrixElems_10_10;
    fRMatrixElems = fRMatrixElems_10_10;
    fFPMatrixElems = fFPMatrixElems_10_10;
  }
  else if(x_fp > 0.10 && x_fp < 0.30){
    fTMatrixElems = fTMatrixElems_10_30;
    fDMatrixElems = fDMatrixElems_10_30;
    fPMatrixElems = fPMatrixElems_10_30;
    fPTAMatrixElems = fPTAMatrixElems_10_30;
    fYMatrixElems = fYMatrixElems_10_30;
    fYTAMatrixElems = fYTAMatrixElems_10_30;
    fRMatrixElems = fRMatrixElems_10_30;
    fFPMatrixElems = fFPMatrixElems_10_30;
  }
  else {
    fTMatrixElems = fTMatrixElems_30_50;
    fDMatrixElems = fDMatrixElems_30_50;
    fPMatrixElems = fPMatrixElems_30_50;
    fPTAMatrixElems = fPTAMatrixElems_30_50;
    fYMatrixElems = fYMatrixElems_30_50;
    fYTAMatrixElems = fYTAMatrixElems_30_50;
    fRMatrixElems = fRMatrixElems_30_50;
    fFPMatrixElems = fFPMatrixElems_30_50;
  }
  
  // calculate the matrices we need
  CalcMatrix(x_fp, fDMatrixElems);
  //CalcMatrix(x_fp, fTMatrixElems);
  // CalcMatrix(x_fp, fYMatrixElems);
  // CalcMatrix(x_fp, fYTAMatrixElems);
  //CalcMatrix(x_fp, fPMatrixElems);
  //CalcMatrix(x_fp, fPTAMatrixElems);
  
  // calculate the coordinates at the target
  dp = CalcTargetVar(fDMatrixElems, powers);
  
  return dp; 

}


double ROpticsOpt::calc_vz(int event, double y, double ph){

  EventData &eventdata = fRawData[event];
  TVector3 BeamSpotHCS(eventdata.Data[kBeamX], eventdata.Data[kBeamY], 0);
  
  const Int_t a = (HRSAngle > 0) ? 1 : -1;
  double CalcReacZ = - ( y -a*MissPointZ)*TMath::Cos(TMath::ATan(ph))/TMath::Sin(HRSAngle + TMath::ATan(ph)) + BeamSpotHCS.X()*TMath::Cos(HRSAngle+TMath::ATan(ph))/TMath::Sin(HRSAngle+TMath::ATan(ph));
  
   
  return CalcReacZ; 

}



double ROpticsOpt::sieve_x(int event){

  std::vector<THaMatrixElement> fTMatrixElems;
  std::vector<THaMatrixElement> fDMatrixElems;
  std::vector<THaMatrixElement> fPMatrixElems;
  std::vector<THaMatrixElement> fPTAMatrixElems;
  std::vector<THaMatrixElement> fYMatrixElems;
  std::vector<THaMatrixElement> fYTAMatrixElems;
  std::vector<THaMatrixElement> fFPMatrixElems;
  std::vector<THaMatrixElement> fRMatrixElems;

  EventData &eventdata = fRawData[event];


  UInt_t FoilID = 1; //starting 0!

  

  Double_t x_fp = eventdata.Data[kX];
  const Double_t(*powers)[5] = eventdata.powers;


if(x_fp > -0.50 && x_fp < -0.30){
    fTMatrixElems = fTMatrixElems_50_30;
    fDMatrixElems = fDMatrixElems_50_30;
    fPMatrixElems = fPMatrixElems_50_30;
    fPTAMatrixElems = fPTAMatrixElems_50_30;
    fYMatrixElems = fYMatrixElems_50_30;
    fYTAMatrixElems = fYTAMatrixElems_50_30;
    fRMatrixElems = fRMatrixElems_50_30;
    fFPMatrixElems = fFPMatrixElems_50_30;
  }

  if(x_fp > -0.30 && x_fp < -0.10){
    fTMatrixElems = fTMatrixElems_30_10;
    fDMatrixElems = fDMatrixElems_30_10;
    fPMatrixElems = fPMatrixElems_30_10;
    fPTAMatrixElems = fPTAMatrixElems_30_10;
    fYMatrixElems = fYMatrixElems_30_10;
    fYTAMatrixElems = fYTAMatrixElems_30_10;
    fRMatrixElems = fRMatrixElems_30_10;
    fFPMatrixElems = fFPMatrixElems_30_10;
  }

  if(x_fp > -0.10 && x_fp < 0.10){
    fTMatrixElems = fTMatrixElems_10_10;
    fDMatrixElems = fDMatrixElems_10_10;
    fPMatrixElems = fPMatrixElems_10_10;
    fPTAMatrixElems = fPTAMatrixElems_10_10;
    fYMatrixElems = fYMatrixElems_10_10;
    fYTAMatrixElems = fYTAMatrixElems_10_10;
    fRMatrixElems = fRMatrixElems_10_10;
    fFPMatrixElems = fFPMatrixElems_10_10;
  }

  if(x_fp > 0.10 && x_fp < 0.30){
    fTMatrixElems = fTMatrixElems_10_30;
    fDMatrixElems = fDMatrixElems_10_30;
    fPMatrixElems = fPMatrixElems_10_30;
    fPTAMatrixElems = fPTAMatrixElems_10_30;
    fYMatrixElems = fYMatrixElems_10_30;
    fYTAMatrixElems = fYTAMatrixElems_10_30;
    fRMatrixElems = fRMatrixElems_10_30;
    fFPMatrixElems = fFPMatrixElems_10_30;
  }

  if(x_fp > 0.30 && x_fp < 0.50){
    fTMatrixElems = fTMatrixElems_30_50;
    fDMatrixElems = fDMatrixElems_30_50;
    fPMatrixElems = fPMatrixElems_30_50;
    fPTAMatrixElems = fPTAMatrixElems_30_50;
    fYMatrixElems = fYMatrixElems_30_50;
    fYTAMatrixElems = fYTAMatrixElems_30_50;
    fRMatrixElems = fRMatrixElems_30_50;
    fFPMatrixElems = fFPMatrixElems_30_50;
  }
  

  CalcMatrix(x_fp, fTMatrixElems);

  eventdata.Data[kCalcTh] = CalcTargetVar(fTMatrixElems, powers);

  const TVector3 BeamSpotHCS(eventdata.Data[kBeamX], eventdata.Data[kBeamY], targetfoils[FoilID]);
  const TVector3 BeamSpotTCS = fTCSInHCS.Inverse()*(BeamSpotHCS - fPointingOffset);
  const Double_t x_tg = BeamSpotTCS.X() - BeamSpotTCS.Z() * tan(eventdata.Data[kCalcTh]);
  
  double ProjectionX = x_tg + (tan(eventdata.Data[kCalcTh]) + x_tg * ExtTarCor_ThetaCorr) * (ZPos);


  return ProjectionX;
}

double ROpticsOpt::sieve_y(int event){

  std::vector<THaMatrixElement> fTMatrixElems;
  std::vector<THaMatrixElement> fDMatrixElems;
  std::vector<THaMatrixElement> fPMatrixElems;
  std::vector<THaMatrixElement> fPTAMatrixElems;
  std::vector<THaMatrixElement> fYMatrixElems;
  std::vector<THaMatrixElement> fYTAMatrixElems;
  std::vector<THaMatrixElement> fFPMatrixElems;
  std::vector<THaMatrixElement> fRMatrixElems;
    
  EventData &eventdata = fRawData[event];

  UInt_t FoilID = 1; //starting 0!
    
  Double_t x_fp = eventdata.Data[kX];
  const Double_t(*powers)[5] = eventdata.powers;


if(x_fp > -0.50 && x_fp < -0.30){
    fTMatrixElems = fTMatrixElems_50_30;
    fDMatrixElems = fDMatrixElems_50_30;
    fPMatrixElems = fPMatrixElems_50_30;
    fPTAMatrixElems = fPTAMatrixElems_50_30;
    fYMatrixElems = fYMatrixElems_50_30;
    fYTAMatrixElems = fYTAMatrixElems_50_30;
    fRMatrixElems = fRMatrixElems_50_30;
    fFPMatrixElems = fFPMatrixElems_50_30;
  }

  if(x_fp > -0.30 && x_fp < -0.10){
    fTMatrixElems = fTMatrixElems_30_10;
    fDMatrixElems = fDMatrixElems_30_10;
    fPMatrixElems = fPMatrixElems_30_10;
    fPTAMatrixElems = fPTAMatrixElems_30_10;
    fYMatrixElems = fYMatrixElems_30_10;
    fYTAMatrixElems = fYTAMatrixElems_30_10;
    fRMatrixElems = fRMatrixElems_30_10;
    fFPMatrixElems = fFPMatrixElems_30_10;
  }

  if(x_fp > -0.10 && x_fp < 0.10){
    fTMatrixElems = fTMatrixElems_10_10;
    fDMatrixElems = fDMatrixElems_10_10;
    fPMatrixElems = fPMatrixElems_10_10;
    fPTAMatrixElems = fPTAMatrixElems_10_10;
    fYMatrixElems = fYMatrixElems_10_10;
    fYTAMatrixElems = fYTAMatrixElems_10_10;
    fRMatrixElems = fRMatrixElems_10_10;
    fFPMatrixElems = fFPMatrixElems_10_10;
  }

  if(x_fp > 0.10 && x_fp < 0.30){
    fTMatrixElems = fTMatrixElems_10_30;
    fDMatrixElems = fDMatrixElems_10_30;
    fPMatrixElems = fPMatrixElems_10_30;
    fPTAMatrixElems = fPTAMatrixElems_10_30;
    fYMatrixElems = fYMatrixElems_10_30;
    fYTAMatrixElems = fYTAMatrixElems_10_30;
    fRMatrixElems = fRMatrixElems_10_30;
    fFPMatrixElems = fFPMatrixElems_10_30;
  }

  if(x_fp > 0.30 && x_fp < 0.50){
    fTMatrixElems = fTMatrixElems_30_50;
    fDMatrixElems = fDMatrixElems_30_50;
    fPMatrixElems = fPMatrixElems_30_50;
    fPTAMatrixElems = fPTAMatrixElems_30_50;
    fYMatrixElems = fYMatrixElems_30_50;
    fYTAMatrixElems = fYTAMatrixElems_30_50;
    fRMatrixElems = fRMatrixElems_30_50;
    fFPMatrixElems = fFPMatrixElems_30_50;
  }


  CalcMatrix(x_fp, fPMatrixElems);
  CalcMatrix(x_fp, fPTAMatrixElems);
  
  // calculate the coordinates at the target
  eventdata.Data[kCalcPh] = CalcTargetVar(fPMatrixElems, powers) + CalcTargetVar(fPTAMatrixElems, powers);

  const TVector3 BeamSpotHCS(eventdata.Data[kBeamX], eventdata.Data[kBeamY], targetfoils[FoilID]);
  const TVector3 BeamSpotTCS = fTCSInHCS.Inverse()*(BeamSpotHCS - fPointingOffset);
  const Double_t y_tg = BeamSpotTCS.Y() - BeamSpotTCS.Z() * tan(eventdata.Data[kCalcPh]);

  
  double ProjectionY = y_tg + tan(eventdata.Data[kCalcPh]) * (ZPos);

  

  return ProjectionY;
  
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

ClassImp(ROpticsOpt);
