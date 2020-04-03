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

#ifndef ROOT_ROpticsOpt
#define ROOT_ROpticsOpt

#include <vector>

#include "TRotation.h"
#include "TMath.h"

#include "THaTrackingDetector.h"
#include "THaString.h"

///////////////////////////////////////////////////////////////////////////////
// Debug Definitions
// place this section below any other head files
///////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG_LEVEL
#undef DEBUG_LEVEL
#endif
//     DEBUG_LEVEL;
//     =0      or not define: no debug, full speed
//     >=1     enable debug extra warning (suggested setting)
//     >=2     above + enable debug assert
//     >=3     above + enable debug extra info
//     >=4     above + massive info (in a for or while)
//     >=5     Decode dump
#define DEBUG_LEVEL 3
#include "DebugDef.h"
#include "InputR.h"

///////////////////////////////////////////////////////////////////////////////
class TCanvas;
class THaTrack;
class TClonesArray;
class TTree;
class TVector3;

class THaMatrixElement;

class ROpticsOpt : public THaTrackingDetector {
public:
    ROpticsOpt(const char* name = "Optimizer", const char* description = "Optimizer for HRS Optics", THaApparatus* a = NULL);
    virtual ~ROpticsOpt();

    ///////////////////////////////////////////////////////////////////////////
    // Database input/output
    ///////////////////////////////////////////////////////////////////////////
    TString OldComments;
    Int_t LoadDataBase(TString DataBaseName, TString xfp_range); // Database file -> Memory
    Int_t SaveDataBase(TString DataBaseName); // Memory -> Database file

       
    enum {
        MaxNEventData = 50, MaxNRawData = 2000000, kNUM_PRECOMP_POW = 10, kMaxDataGroup = 180 * 5 * 5
    };
    
    
    
    //UInt_t LoadRawData(TString DataFileName, UInt_t NLoad = MaxNRawData, UInt_t MaxDataPerGroup = (UInt_t) - 1); // load data to Rawdata[]
    UInt_t LoadRawData(TTree *t); // load data to Rawdata[]
    
    //typedef struct {
    struct EventData {
        Double_t Data[MaxNEventData]; // [CommonIdx]
        Double_t powers[kNUM_PRECOMP_POW][5]; // {(x), th, y, ph, abs(th) }
    };
    //} EventData;
    EventData fRawData[MaxNRawData]; // [fNRawData]
    UInt_t fNRawData;
    UInt_t fNCalibData; // for dp calib only

    enum CommonIdx {
        kX = 0, // R.tr.r_x
        kTh = 1, // R.tr.r_th
        kY = 2, // R.tr.r_y
        kPhi = 3, // R.tr.r_ph
        kBeamX = 5, // urb.x or rb.x
        kBeamY = 6, // urb.y or rb.y
	kBeamVZ = 7, //R.tr.vz
	kCutID = 8, // cut ID in order of tree2ascii cut file
    };

    enum ExtraSieveIdx {
        kRealTh = 30, // real target th from survey
        kRealPhi, // real target ph from survey
        kRealTgX, // real target x from survey, beam
        kRealTgY, // real target y from survey, beam
        kRealThMatrix, // expected target th before extended target corrections
        kCalcTh, // calculated th from matrix
        kCalcPh, // calculated ph from matrix
        kSieveX,
        kSieveY,
        kSieveZ,
        kBeamZ
    };

    enum ExtraVertexIdx {
        kRealReactZ = 10, //expected ReactZ
        kCalcTgY, //calculated Tg_y
        kCalcReactZ //calculated ReactZ
    };

    enum ExtraDpIdx {
        kExtraDataFlag = 15, //Whether this event is for optimization; 0=used for optimization, 1=for plotting only
        kKineID, //Delta Scan Kinematics ID
        kCentralp, //Central Momentum
        kElossDp, //Radiation Loss for this event in unit of dp
        kScatterAngle, //Scattering Angle
        kDpKinOffsets, //=dp-dp_kin, based on sieve hole survey
        kRealDpKin, //expected dp_kin, before radiation correction
        kRealDpKinMatrix, //real dp kin before extended target corrections
        kCalcDpKinMatrix, //calculated dp kin before extended target corrections
        kCalcDpKin, //calculated dp_kin, before radiation correction
        kRealDpKinExcitations,/*Do not append more index*/ //first index of expected dp_kins for all excitation states
	kElossTgBefore,
	kElossTgAfter,
	kTravelLengthBefore,
	kTravelLengthAfter
    };

    ///////////////////////////////////////////////////////////////////////////
    // Optimization related Commands
    ///////////////////////////////////////////////////////////////////////////
  
    Double_t fArbitaryVertexShift[100]; // compensate bias due to event selections, array of [FoilID]
    
    Double_t fArbitaryDpKinShift[100]; // compensate bias due to dp event selections, array of [KineID]
    
    TRotation fTCSInHCS; // transformations vector from TCS to HCS
    TVector3 fPointingOffset; // Optical point in lab coordinate system

    inline Double_t ScatMom(Double_t DM, Double_t Ma, Double_t P0, Double_t Theta)
    {
        // Calc Scattered Electron Momentum
        // Assuming Electron with P0, Scattering on a fix target of Mass Ma Assuming
        // recoil system is a resonance of DM and scattering angle is Theta.
        return (-DM * DM - 2 * Ma * DM + 2 * Ma * P0) / (2 * (Ma + P0 - P0 * TMath::Cos(Theta)));
    }

    ///////////////////////////////////////////////////////////////////////////
    // declarations for target vertex reconstruction
    ///////////////////////////////////////////////////////////////////////////

    enum ECoordTypes {
        kTransport, kRotatingTransport
    };

    enum EFPMatrixElemTags {
        T000 = 0, Y000, P000
    };

    enum {
        kPORDER = 7
    };

    friend class THaMatrixElement;
    std::vector<THaMatrixElement> * fCurrentMatrixElems;
    // initial matrix elements
    std::vector<THaMatrixElement> fTMatrixElems_50_30;
    std::vector<THaMatrixElement> fDMatrixElems_50_30;
    std::vector<THaMatrixElement> fPMatrixElems_50_30;
    std::vector<THaMatrixElement> fPTAMatrixElems_50_30; // involves abs(theta_fp)
    std::vector<THaMatrixElement> fYMatrixElems_50_30;
    std::vector<THaMatrixElement> fYTAMatrixElems_50_30; // involves abs(theta_fp)
    std::vector<THaMatrixElement> fFPMatrixElems_50_30; // matrix elements used in
    std::vector<THaMatrixElement> fRMatrixElems_50_30; // Path-length corrections (meters)

    
    std::vector<THaMatrixElement> * fCurrenteMatrixElems_30_10;
    std::vector<THaMatrixElement> fTMatrixElems_30_10;
    std::vector<THaMatrixElement> fDMatrixElems_30_10;
    std::vector<THaMatrixElement> fPMatrixElems_30_10;
    std::vector<THaMatrixElement> fPTAMatrixElems_30_10; // involves abs(theta_fp)
    std::vector<THaMatrixElement> fYMatrixElems_30_10;
    std::vector<THaMatrixElement> fYTAMatrixElems_30_10; // involves abs(theta_fp)
    std::vector<THaMatrixElement> fFPMatrixElems_30_10; // matrix elements used in
    std::vector<THaMatrixElement> fRMatrixElems_30_10; // Pathlength corrections (meters)

    std::vector<THaMatrixElement> * fCurrentMatrixElems_10_10;
    // initial matrix elements
    std::vector<THaMatrixElement> fTMatrixElems_10_10;
    std::vector<THaMatrixElement> fDMatrixElems_10_10;
    std::vector<THaMatrixElement> fPMatrixElems_10_10;
    std::vector<THaMatrixElement> fPTAMatrixElems_10_10; // involves abs(theta_fp)
    std::vector<THaMatrixElement> fYMatrixElems_10_10;
    std::vector<THaMatrixElement> fYTAMatrixElems_10_10; // involves abs(theta_fp)
    std::vector<THaMatrixElement> fFPMatrixElems_10_10; // matrix elements used in
    std::vector<THaMatrixElement> fRMatrixElems_10_10; // Path-length corrections (meters)

    std::vector<THaMatrixElement> * fCurrentMatrixElems_10_30;
    // initial matrix elements
    std::vector<THaMatrixElement> fTMatrixElems_10_30;
    std::vector<THaMatrixElement> fDMatrixElems_10_30;
    std::vector<THaMatrixElement> fPMatrixElems_10_30;
    std::vector<THaMatrixElement> fPTAMatrixElems_10_30; // involves abs(theta_fp)
    std::vector<THaMatrixElement> fYMatrixElems_10_30;
    std::vector<THaMatrixElement> fYTAMatrixElems_10_30; // involves abs(theta_fp)
    std::vector<THaMatrixElement> fFPMatrixElems_10_30; // matrix elements used in

    std::vector<THaMatrixElement> fRMatrixElems_10_30; // Path-length corrections (meters)

    std::vector<THaMatrixElement> * fCurrentMatrixElems_30_50;
    // initial matrix elements
    std::vector<THaMatrixElement> fTMatrixElems_30_50;
    std::vector<THaMatrixElement> fDMatrixElems_30_50;
    std::vector<THaMatrixElement> fPMatrixElems_30_50;
    std::vector<THaMatrixElement> fPTAMatrixElems_30_50; // involves abs(theta_fp)
    std::vector<THaMatrixElement> fYMatrixElems_30_50;
    std::vector<THaMatrixElement> fYTAMatrixElems_30_50; // involves abs(theta_fp)
    std::vector<THaMatrixElement> fFPMatrixElems_30_50; // matrix elements used in

    std::vector<THaMatrixElement> fRMatrixElems_30_50; // Path-length corrections (meters)
    
    void CalcMatrix(const double x, std::vector<THaMatrixElement> &matrix);
    // Double_t DoPoly(const int n, const std::vector<double> &a, const double x);
    // Double_t PolyInv(const double x1, const double x2, const double xacc, const double y, const int norder, const std::vector<double> &a);
    Double_t CalcTargetVar(const std::vector<THaMatrixElement> &matrix, const double powers[][5]);
    double calc_tgy(int event);
    double calc_tgth(int event);
    double calc_tgph(int event);
    double calc_tgdp(int event);
    double sieve_x(int event);
    double sieve_y(int event);

    ///////////////////////////////////////////////////////////////////////////
    // Inherited from THaTrackingDetector
    ///////////////////////////////////////////////////////////////////////////

    virtual Int_t Decode(const THaEvData&)
    {
        return 0;
    }

    virtual Int_t CoarseTrack(TClonesArray&)
    {
        return 0;
    }

    virtual Int_t FineTrack(TClonesArray&)
    {
        return 0;
    }

    virtual EStatus Init(const TDatime&)
    {
        return fStatus = kOK;
    };

private:
    ClassDef(ROpticsOpt, 0) // HRS Optics Optimizer
};

///////////////////////////////////////////////////////////////////////////////
// class for storing matrix element data
///////////////////////////////////////////////////////////////////////////////

class THaMatrixElement {
public:

    THaMatrixElement() : iszero(true), pw(3), order(0), v(0), poly(ROpticsOpt::kPORDER), OptOrder(0)
    {
    }
    bool match(const THaMatrixElement& rhs) const;

    bool iszero; // whether the element is zero
    std::vector<int> pw; // exponents of matrix element
    //   e.g. D100 = { 1, 0, 0 }
    int order;
    double v; // its computed value
    std::vector<double> poly; // the associated polynomial

    void SkimPoly(); //reduce order to highest non-zero poly

    UInt_t OptOrder; //order optimize to
};
#endif
