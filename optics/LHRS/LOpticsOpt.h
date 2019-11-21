///////////////////////////////////////////////////////////////////////////////
//
// LOpticsOpt
//
// http://www.jlab.org/~jinhuang/HRSOptics/
//
// HRS optics matrix optimization class
// Based on THaVDC
//
//
// Units used:
//        For X, Y, and Z coordinates of track    -  meters
//        For Theta and Phi angles of track       -  tan(angle)
//
// Author: Jin Huang <jinhuang@jlab.org>
// Modification:
//			Jun 25, 2010 Updated for APEX optics calibration
//
///////////////////////////////////////////////////////////////////////////////



#ifndef ROOT_LOpticsOpt
#define ROOT_LOpticsOpt


#include <THaTrackingDetector.h>
#include "TRotation.h"
#include "TMath.h"
#include <vector>
#include "InputAPEXL.h"

//------------------------------------------------------//
//
//      Debug Definitions
//      place this section below any other head files
//
//------------------------------------------------------//
#ifdef DEBUG_LEVEL
#   undef DEBUG_LEVEL
#endif

//      DEBUG_LEVEL;
//      =0      or not define: no debug, full speed
//      >=1     enable debug extra warning (suggested setting)
//      >=2     above + enable debug assert
//      >=3     above + enable debug extra info
//      >=4     above + massive info (in a for or while)
//  >=5 Decode dump
#define DEBUG_LEVEL  3
#include    "DebugDef.h"
//------------------------------------------------------//


#include "THaString.h"
class TCanvas;
class THaTrack;
class TClonesArray;
class TTree;
class TVector3;

class THaMatrixElement;

class LOpticsOpt : public THaTrackingDetector {

	public:
		/////////////////////////////////////////////////////////////////////////////////
		// Database input/output

		TString OldComments;
		Int_t LoadDataBase(TString DataBaseName);//Database file -> Memory
		Int_t SaveDataBase(TString DataBaseName);//Memory -> Database file

		virtual void Print(const Option_t* opt) const;

		UInt_t Matrix2Array(Double_t Array[],Bool_t FreeParaFlag[]=NULL)  //fCurrentMatrixElems -> Array
		{assert(fCurrentMatrixElems);return Matrix2Array(Array,(*fCurrentMatrixElems),FreeParaFlag);}
		UInt_t Matrix2Array(Double_t Array[], const std::vector<THaMatrixElement> &Matrix,Bool_t FreeParaFlag[]=NULL);

		UInt_t Array2Matrix(const Double_t Array[])  // Array -> fCurrentMatrixElems
		{assert(fCurrentMatrixElems);return Array2Matrix(Array,(*fCurrentMatrixElems));}
		UInt_t Array2Matrix(const Double_t Array[], std::vector<THaMatrixElement> &Matrix);

		/////////////////////////////////////////////////////////////////////////////////
		// Data storage

		UInt_t LoadRawData(TString DataFileName,UInt_t NLoad = MaxNRawData, UInt_t MaxDataPerGroup = (UInt_t)-1); //load data to Rawdata[]
		UInt_t SaveDataBuffer(TTree * T); //save Rawdata[] to T Tree. return N event written
		UInt_t SaveDataBuffer(TString fname, TString tree = "T"); //save Rawdata[] to file

		enum {MaxNEventData = 100,MaxNRawData = 2000000, kNUM_PRECOMP_POW= 10, kMaxDataGroup=180*5*5};
		typedef struct
		{
			Double_t Data[MaxNEventData];	//[CommonIdx]
			Double_t powers[kNUM_PRECOMP_POW][5];  // {(x), th, y, ph, abs(th) }
		} EventData;
		EventData fRawData[MaxNRawData];	//[fNRawData]
		UInt_t fNRawData;
		UInt_t fNCalibData;	// for dp calib only


		// JW changed this from having one kCutID to 3 parameters: kFoilID, kColID & kRowID
		enum CommonIdx{
			kFoilID = 0,
			kColID = 1,
			kRowID = 2,
	   		kX = 3,		//L.tr.r_x
	  		kTh = 4,	//L.tr.r_th
   			kY = 5,		//L.tr.r_y
   			kPhi = 6,		//L.tr.r_ph
   			kBeamX = 7,	//urb.x or rb.x
   			kBeamY = 8,	//urb.y or rb.y
		   	kL_tr_tg_th = 9,	//L.tr.tg_th
   			kL_tr_tg_ph = 10,	//L.tr.tg_ph
			kL_tr_tg_y = 11	        //L.tr.tg_y
		};
		enum ExtraSieveIdx{
			kRealTh = 40,	//real target th from survey
	   		kRealPhi,	//real target ph from survey
	   		kRealTgX,	//real target x from survey, beam
	   		kRealThMatrix,	//expected target th before extended target corrections
	   		kCalcTh,		//calculated th from matrix
	   		kCalcPh		//calculated ph from matrix
		};
		enum ExtraVertexIdx{
		  kBeamDirX = 12,//urb.dir.y or rb.dir.y
		  kBeamDirY,//urb.dir.y or rb.dir.y
		  kBeamDirZ,//urb.dir.y or rb.dir.y
		  kRealTgY=60,	//Expected Tg_y from Survey and
		  kRealReactZ,	//expected ReactZ
		  kCalcTgY,	//calculated Tg_y
		  kCalcReactZ	//calculated ReactZ
		};
		enum ExtraDpIdx{
		  kL_tr_tg_dp=71,	//L.tr.tg_dp
		  kL_tr_p,	//L.tr.p
		  kurb_e,	//Beam energy
		  kRunNum,	//Run number
		  kExtraDataFlag,	//Whether this event is for optimization; 0=used for optimization, 1=for plotting only
		  kKineID,		//Delta Scan Kinematics ID
		  kCentralp,	//Central Momentum
		  kRadiLossDp,	//Radiation Loss for this event in unit of dp
		  kScatterAngle,	//Scattering Angle
		  kDpKinOffsets,	//=dp-dp_kin, based on sieve hole survey
		  kRealDpKin,		//expected dp_kin, before radiation correction
		  kRealDpKinMatrix,	//real dp kin before extended target corrections
		  kCalcDpKinMatrix,	//calculated dp kin before extended target corrections
		  kCalcDpKin,	//calculated dp_kin, before radiation correction
		  kRealDpKinExcitations/*Do not append more index*/ //first index of expected dp_kins for all excitation states
		};

		/////////////////////////////////////////////////////////////////////////////////
		//Optimization related Commands
		/* const TVector3 GetSieveHoleTCS(UInt_t Col, UInt_t Row) /\*const*\/; */

		/* const TVector3 GetSieveHoleTCS(UInt_t Hole) /\*const*\/; */

		/* std::vector<int> Get_Col_Row(Int_t Hole); */


		/* Int_t Get_Hole(Int_t Col, Int_t Row); */


		void PrepareSieve(void);
		Double_t VerifyMatrix_Sieve(void);
		// JW: temporarily commented out CheckSieve
		TCanvas * CheckSieve(UInt_t SpecialNLoad= NFiols);
		Double_t SumSquareDTh(Bool_t PrintEachHole = kFALSE);
		Double_t SumSquareDPhi(Bool_t PrintEachHole = kFALSE);

		Double_t fArbitaryVertexShift[100];	//compensate bias due to event selections, array of [FoilID]
		void PrepareVertex(void);
		Double_t VerifyMatrix_Vertex(void);

		//		typedef TCanvas* (canv_array)[2];
		//		TCanvas** CheckVertex(void);
		TCanvas * CheckVertex(void);
		TCanvas * CheckVertexRes(void);
		TCanvas * CheckVertexDiff(void);


		Double_t SumSquareDTgY();
		Double_t SumSquareDTgYAverFoils();

		Double_t fArbitaryDpKinShift[100];	//compensate bias due to dp event selections, array of [KineID]
		void PrepareDp(void);
		Double_t VerifyMatrix_Dp(void);
		TCanvas* CheckDp(void);
		TCanvas* CheckDpGlobal(void);
		TCanvas* CheckDpVSAngle(void);
		TCanvas* CheckDpVSCutID(void);
		Double_t SumSquareDp(Bool_t IncludeExtraData=kFALSE);


		TRotation fTCSInHCS;		//transformations vector from TCS to HCS
		TVector3 fPointingOffset;	//Optical point in lab coordinate system


		inline Double_t ScatMom(Double_t DM, Double_t Ma, Double_t P0, Double_t Theta)
		{
			//Calc Scattered Electron Momentum
			// Assuming Electron with P0, Scattering on a fix target of Mass Ma Assuming
			// recoil system is a resonance of DM and scattering angle is Theta.
			return (-DM*DM - 2*Ma*DM + 2*Ma*P0)/(2*(Ma + P0 - P0*TMath::Cos(Theta)));
		}

		/////////////////////////////////////////////////////////////////////////////////
		// Heritaged from THaVDC

		LOpticsOpt( const char* name="Optimizer", const char* description = "Optimizer for HRS Optics",
					THaApparatus* a = NULL );

		virtual ~LOpticsOpt();



		/////////////////////////////////////////////////////////////////////////////////

  // Bits & and bit masks for THaTrack
		enum {
			kStageMask     = BIT(14) | BIT(15),  // Track processing stage bits
								 kInvalid       = 0x0000,  // Not processed
		 kCoarse        = BIT(14), // Coarse track
							  kFine          = BIT(15), // Fine track
									  kReassigned    = BIT(16), // Track is a new track in Fine stage
											  kMultiTrack    = BIT(17), // Track was generated in the multitrack analysis
													  kBadTrack      = BIT(18)  // Track prematurely exits the spectrometer or similar
		};

  // Bits and bit masks for this object
		enum {
			kOnlyFastest    = BIT(14), // Only use fastest hit for each wire (highest TDC)
								  kTDCbits        = BIT(15) | BIT(16),  // Mask for TDC mode bits
										  kHardTDCcut     = BIT(15), // Use hard TDC cuts (fMinTime, fMaxTime)
												  kSoftTDCcut     = BIT(16), // Use soft TDC cut (reasonable estimated drifts)
														  kIgnoreNegDrift = BIT(17), // Completely ignore negative drift times

																  kCoarseOnly     = BIT(23) // Do only coarse tracking
		};


  // declarations for target vertex reconstruction
		enum ECoordTypes { kTransport, kRotatingTransport };
		enum EFPMatrixElemTags { T000 = 0, Y000, P000 };
		enum { kPORDER = 7 };


		friend class THaMatrixElement;
		std::vector<THaMatrixElement> * fCurrentMatrixElems;
  // initial matrix elements
		std::vector<THaMatrixElement> fTMatrixElems;
		std::vector<THaMatrixElement> fDMatrixElems;
		std::vector<THaMatrixElement> fPMatrixElems;
		std::vector<THaMatrixElement> fPTAMatrixElems; // involves abs(theta_fp)
		std::vector<THaMatrixElement> fYMatrixElems;
		std::vector<THaMatrixElement> fYTAMatrixElems; // involves abs(theta_fp)
		std::vector<THaMatrixElement> fFPMatrixElems;  // matrix elements used in
                                            // focal plane transformations
                                            // { T, Y, P }

		std::vector<THaMatrixElement> fLMatrixElems;   // Path-length corrections (meters)

		/////////////////////////////////////////////////////////////////////////////////
		void CalcTargetCoords(THaTrack *the_track, const ECoordTypes mode);
		void CalcMatrix(const double x, std::vector<THaMatrixElement> &matrix);
// 		  Double_t DoPoly(const int n, const std::vector<double> &a, const double x);
// 		  Double_t PolyInv(const double x1, const double x2, const double xacc,
// 						   const double y, const int norder,
// 		 const std::vector<double> &a);
		Double_t CalcTargetVar(const std::vector<THaMatrixElement> &matrix,
							   const double powers[][5]);
		Double_t CalcTarget2FPLen(const std::vector<THaMatrixElement>& matrix,
								  const Double_t powers[][5]);


		/////////////////////////////////////////////////////////////////////////////////
		virtual Int_t ConstructTracks( TClonesArray* tracks = NULL, Int_t flag = 0 );


		virtual void  Clear( Option_t* opt="" );
		virtual Int_t Decode( const THaEvData& );
		virtual Int_t CoarseTrack( TClonesArray& tracks );
		virtual Int_t FineTrack( TClonesArray& tracks );
		virtual Int_t FindVertices( TClonesArray& tracks );
		virtual EStatus Init( const TDatime& date );


		ClassDef(LOpticsOpt,0)             // HRS Optics Optimizer
};

////////////////////////////////////////////////////////////////////////////////

// class for storing matrix element data
class THaMatrixElement {
	  public:
		  THaMatrixElement() : iszero(true), pw(3), order(0), v(0), poly(LOpticsOpt::kPORDER),OptOrder(0) {}
		  bool match( const THaMatrixElement& rhs ) const;

		  bool iszero;             // whether the element is zero
		  std::vector<int> pw;     // exponents of matrix element
							 //   e.g. D100 = { 1, 0, 0 }
		  int  order;
		  double v;                // its computed value
		  std::vector<double> poly;// the associated polynomial

		  void SkimPoly(); //reduce order to highest non-zero poly

		  UInt_t OptOrder;		//order optimize to
  };
#endif
