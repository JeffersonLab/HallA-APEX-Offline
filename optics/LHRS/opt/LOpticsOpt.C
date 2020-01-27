///////////////////////////////////////////////////////////////////////////////
//
// LOpticsOpt
//
// http://www.jlab.org/~jinhuang/HRSOptics/
//
// HRS optics matrix optimization class
// Based on THaVDC
//
// Units used:
//        For X, Y, and Z coordinates of track    -  meters
//        For Theta and Phi angles of track       -  tan(angle)
//        For Momentums, Masses                   -  GeV, GeV/c^2
//
// Author: Jin Huang <jinhuang@jlab.org>
// Modification:
//			Jun 25, 2010 Updated for APEX optics calibration
//
//////////////////////////////////////////////////////////////////////////////

#include "THaGlobals.h"
#include "THaEvData.h"
#include "THaDetMap.h"
#include "THaTrack.h"
#include "THaScintillator.h"
#include "THaSpectrometer.h"
#include "TGraphErrors.h"
#include "TClonesArray.h"
#include "TList.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH2.h"
#include "TH3.h"
#include "TH1.h"
#include "TF1.h"
#include "TLatex.h"
#include "TVector3.h"
#include "TLine.h"
#include "TArrow.h"
#include "VarDef.h"
#include "TEllipse.h"
#include "TPaveStats.h"
#include "TStyle.h"


//#include <algorithm>
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TDatime.h"
#include <map>
#include <cstdio>
#include <cstdlib>


#include "LOpticsOpt.h"
#include "APEX_Sieve.h"



#ifdef WITH_DEBUG
#include <iostream>
#endif

using namespace std;
using THaString::Split;

/////////////////////////////////////////////////////////////////////////
// Input Sections
///////////////////////////////////////////////////////////////////////
// #include "InputE06010.h"
#include "InputAPEXL.h"
//#include "InputPREXR.h"
//
//#include "InputAPEXL.h"


//JW: commented here temporarily:

// const UInt_t NHoles = 225;
// const UInt_t NoinEachRow[] = {15, 12, 15, 11, 15, 11, 15, 11, 15, 11, 15, 11, 15, 11, 15, 12, 15};


 


//_____________________________________________________________________________
Double_t LOpticsOpt::SumSquareDp(Bool_t IncludeExtraData)
{
	//return square sum of diff between calculated dp_kin and expected dp_kin

	Double_t d_dp = 0/*, dphi=0*/;	//Difference
	Double_t rms_dp = 0/*, rmsphi=0*/; //mean square

	static UInt_t NCall = 0;
	NCall++;

	UInt_t NCalibData=0;

	if (IncludeExtraData)
	{
		Warning("SumSquareDp","Data Beyond selected excitation state is included in this calculation");
	}
	for (UInt_t idx = 0; idx<fNRawData; idx++)
	{

		Double_t  dp, dp_kin;

		EventData &eventdata = fRawData[idx];

		//jump through data beyond selected excitation states
		if (eventdata.Data[kExtraDataFlag]>0 && !IncludeExtraData) continue;
		NCalibData++;

		Double_t x_fp = eventdata.Data[kX];
		const Double_t (*powers)[5] = eventdata.powers;

  // calculate the matrices we need
		CalcMatrix(x_fp, fDMatrixElems);
// 		CalcMatrix(x_fp, fTMatrixElems);
// 		CalcMatrix(x_fp, fYMatrixElems);
// 		CalcMatrix(x_fp, fYTAMatrixElems);
// 		CalcMatrix(x_fp, fPMatrixElems);
// 		CalcMatrix(x_fp, fPTAMatrixElems);

  // calculate the coordinates at the target
// 		theta = CalcTargetVar(fTMatrixElems, powers);
// 		phi = CalcTargetVar(fPMatrixElems, powers)+CalcTargetVar(fPTAMatrixElems,powers);
// 		y = CalcTargetVar(fYMatrixElems, powers)+CalcTargetVar(fYTAMatrixElems,powers);

  // calculate momentum
		dp = CalcTargetVar(fDMatrixElems, powers);
		dp_kin = dp - eventdata.Data[kDpKinOffsets];

		const UInt_t KineID = (UInt_t)(eventdata.Data[kKineID]);
		assert(KineID<NKine);//check array index size
		const Double_t ArbitaryDpKinShift = fArbitaryDpKinShift[KineID];

		d_dp += dp_kin - eventdata.Data[kRealDpKinMatrix]+ArbitaryDpKinShift;
		rms_dp += (dp_kin - eventdata.Data[kRealDpKinMatrix]+ArbitaryDpKinShift)
				*(dp_kin - eventdata.Data[kRealDpKinMatrix]+ArbitaryDpKinShift);

		DEBUG_MASSINFO("SumSquareDp","d_dp = %f = \t%f - \t%f",
					   dp_kin - eventdata.Data[kRealDpKinMatrix], dp_kin , eventdata.Data[kRealDpKinMatrix]  );

		//save the results
		eventdata.Data[kCalcDpKinMatrix] = dp_kin;
		eventdata.Data[kCalcDpKin] = dp_kin+eventdata.Data[kRealTgX]/ExtTarCor_DeltaCorr;
	}

	if (!IncludeExtraData)
		assert(fNCalibData==NCalibData); // check number of event for calibration

	DEBUG_INFO("SumSquareDp","#%d : d_dp = %f,rms_dp=%f",NCall,
			   d_dp/NCalibData,TMath::Sqrt(rms_dp/NCalibData));

	return rms_dp;
}

//_____________________________________________________________________________
Double_t LOpticsOpt::SumSquareDTh(Bool_t PrintEachHole)
{
	//return square sum of diff between calculated tg_th and expected tg_th

	Double_t dth = 0/*, dphi=0*/;	//Difference
	Double_t rmsth = 0/*, rmsphi=0*/; //mean square

	// same parameters but for sieve x
	Double_t dx = 0; //difference
	Double_t rmsdx = 0; // mean square


	// parameters for each individual hole 

	Double_t dth_hole[NHoles] = {0}; 
	Double_t rmsth_hole[NHoles] = {0};
	Double_t hole_stat[NHoles] = {0};
	Double_t rmsth_hole_av = 0;

	Double_t no_hole_used = 0;

	Double_t dx_hole[NHoles] = {0}; 
	Double_t rmsx_hole[NHoles] = {0}; 
	Double_t rmsx_hole_av = 0;

	
	// parameters for foil optimisation

	Double_t dth_foil[NFiols] = {0}; 
	Double_t rmsth_foil[NFiols] = {0};
	Double_t foil_stat[NFiols] = {0};
	Double_t rmsth_foil_av = 0;

	Double_t no_foil_used = 0;

	Double_t dx_foil[NFiols] = {0}; 
	Double_t rmsx_foil[NFiols] = {0}; 
	Double_t rmsx_foil_av = 0;


	static UInt_t NCall = 0;
	NCall++;

	Double_t theta/*, phi, dp, p, pathl*/;

	// stat. for each hole
	Double_t * stat_dev = NULL;
	Double_t * stat_rms = NULL;
	Int_t * stat_cnt = NULL;
	if (PrintEachHole)
	{
		stat_dev = new Double_t[kMaxDataGroup];
		stat_rms = new Double_t[kMaxDataGroup];
		stat_cnt = new Int_t[kMaxDataGroup];
		

		for(int i=0; i<kMaxDataGroup; i++) stat_dev[i] = stat_rms[i] = stat_cnt[i] = 0;
	}
	

	for (UInt_t idx = 0; idx<fNRawData; idx++)
	{
	  EventData &eventdata = fRawData[idx];



	  Double_t x_fp = eventdata.Data[kX];
	  const Double_t (*powers)[5] = eventdata.powers;

  // calculate the matrices we need
// 		CalcMatrix(x_fp, fDMatrixElems);
		CalcMatrix(x_fp, fTMatrixElems);
// 		CalcMatrix(x_fp, fYMatrixElems);
// 		CalcMatrix(x_fp, fYTAMatrixElems);
// 		CalcMatrix(x_fp, fPMatrixElems);
// 		CalcMatrix(x_fp, fPTAMatrixElems);

  // calculate the coordinates at the target
		theta = CalcTargetVar(fTMatrixElems, powers);
// 		phi = CalcTargetVar(fPMatrixElems, powers)+CalcTargetVar(fPTAMatrixElems,powers);
// 		y = CalcTargetVar(fYMatrixElems, powers)+CalcTargetVar(fYTAMatrixElems,powers);

  // calculate momentum
// 		dp = CalcTargetVar(fDMatrixElems, powers);



		// retrieve FoilID of event (necessary for per-foil optimisation)
		const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID];

		// necessary to get hole id of event
		const UInt_t Col = (UInt_t)eventdata.Data[kColID]; 

		const UInt_t Row = (UInt_t)eventdata.Data[kRowID]; 
		
		const UInt_t HoleID = Get_Hole(Col,Row);
	    
		

		//save the results
		eventdata.Data[kCalcTh] = theta;
// 		eventdata.Data[kCalcPh] = eventdata.Data[kL_tr_tg_ph];


		dth += theta - eventdata.Data[kRealThMatrix];
		rmsth += (theta - eventdata.Data[kRealThMatrix])*(theta - eventdata.Data[kRealThMatrix]);



		// save per-hole results

		dth_hole[HoleID] += theta - eventdata.Data[kRealThMatrix];

		rmsth_hole[HoleID] += (theta - eventdata.Data[kRealThMatrix])*(theta - eventdata.Data[kRealThMatrix]);

		hole_stat[HoleID]++;


		
		// save per-foil results

		dth_foil[FoilID] += theta - eventdata.Data[kRealThMatrix];

		rmsth_foil[FoilID] += (theta - eventdata.Data[kRealThMatrix])*(theta - eventdata.Data[kRealThMatrix]);

		foil_stat[FoilID]++;
		

	       
		
		double ProjectionX = eventdata.Data[kRealTgX] + (TMath::Tan(eventdata.Data[kCalcTh]) + eventdata.Data[kRealTgX] * ExtTarCor_ThetaCorr) * (ZPos);
		

		


		

		// below gives default vales for row and column if only foils have been cut (and thus sieve hole row and column unknown)
		// this makes the information about kRealTh and kRealPh meaningless


		TVector3 SieveHoleTCS;
		



		if(Col > 0 && Col < NSieveCol && Row > 0 && Row < NSieveRow){
		  
		  SieveHoleTCS = GetSieveHoleTCS(Col,Row);
		}
		else{
		  SieveHoleTCS = GetSieveHoleTCS(1,1);
		}


		const Double_t posx = SieveHoleTCS.X();
	

	     		//		const Double_t posx = GetSieveHoleTCS(Col,Row).X();


	
		dx += ProjectionX - posx;
		rmsdx += (ProjectionX - posx)*(ProjectionX - posx);

		
		// per-hole results
		dx_hole[HoleID] += ProjectionX - posx;
		rmsx_hole[HoleID] += (ProjectionX - posx)*(ProjectionX - posx);


		// per-foil results
		dx_foil[FoilID] += ProjectionX - posx;
		rmsx_foil[FoilID] += (ProjectionX - posx)*(ProjectionX - posx);


	
		// Calculate difference in calculated x_sieve and x_sieve (from survey)

		// const TVector3 SieveHoleTCS = GetSieveHoleTCS(Col,Row);


		// TVector3 BeamSpotHCS(eventdata.Data[kBeamX],eventdata.Data[kBeamY],targetfoils[FoilID]);

		// TRotation fTCSInHCS;
		// TVector3 TCSX(0,-1,0);
		// TVector3 TCSZ(TMath::Sin(HRSAngle),0,TMath::Cos(HRSAngle));
		// TVector3 TCSY = TCSZ.Cross(TCSX);
		// fTCSInHCS.RotateAxes(TCSX,TCSY,TCSZ);
		// TVector3 BeamSpotTCS=fTCSInHCS.Inverse()*(BeamSpotHCS-fPointingOffset);

		// double ProjectionX = x_tg + (tan(eventdata.Data[kCalcTh]) + x_tg * ExtTarCor_ThetaCorr) * (ZPos);



		



		if (PrintEachHole)
		  {
		    // JW: FIX this -> replaced kCutID with kFoilID		  
		    UInt_t HoleID = (UInt_t)eventdata.Data[kFoilID];
		    assert(HoleID<kMaxDataGroup);
		    
		    stat_dev[HoleID] += theta - eventdata.Data[kRealThMatrix];
		    stat_rms[HoleID] += (theta - eventdata.Data[kRealThMatrix])*(theta - eventdata.Data[kRealThMatrix]);
		    stat_cnt[HoleID] += 1;
		  
		
		// 		dphi += phi - eventdata.Data[kRealPhi];
		// 		rmsphi += (phi - eventdata.Data[kRealPhi])*(phi - eventdata.Data[kRealPhi]);
		
		DEBUG_MASSINFO("SumSquareDTh","D_Th = %f = \t%f - \t%f",
			       theta - eventdata.Data[kRealThMatrix], theta , eventdata.Data[kRealThMatrix]  );
		DEBUG_MASSINFO("SumSquareDTh","%d : %f, %f, %f, %f, %f"
			       ,kRealTh
			       , eventdata.Data[kRealTh-2]
			       , eventdata.Data[kRealTh-1]
			       , eventdata.Data[kRealTh]
			       , eventdata.Data[kRealTh+1]
			       , eventdata.Data[kRealTh+2]
			       );
		  }
#ifdef DISABLED_SIEVE_HOLE
		// JW: altered following lines
		// UInt_t res = (UInt_t)eventdata.Data[kCutID];
		// const UInt_t FoilID = res/(NSieveRow*NSieveCol); //starting 0!
		//		const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID]; //starting 0!
		// res = res%(NSieveRow*NSieveCol);
		// const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
		// const UInt_t Row = (UInt_t)eventdata.Data[kRowID]; //starting 0!

		if (DISABLED_SIEVE_HOLE) continue;
#endif



	}


	
	// loop to go through individual hole data and calculate average rms of all holes

	for(Int_t j = 0; j<NHoles; j++){

	  //cout << "For hole " << j << " dth_hole = " << dth_hole[j] << endl;
	  //	  cout << "          rmsth_hole = " << rmsth_hole[j] << endl;


	  if( hole_stat[j] > 0){

	    dth_hole[j] = dth_hole[j]/hole_stat[j];     
	    rmsth_hole[j] = rmsth_hole[j]/hole_stat[j];

	    dx_hole[j] = dx_hole[j]/hole_stat[j];
	    rmsx_hole[j] = rmsx_hole[j]/hole_stat[j];
	    no_hole_used++;
	  }
	  else{
	    dth_hole[j] = 0;     
	    rmsth_hole[j] = 0;
	    
	    dx_hole[j] = 0;
	    rmsx_hole[j] = 0;

	  }

	  rmsth_hole_av += rmsth_hole[j];
	  rmsx_hole_av += rmsx_hole[j];
	  


	}



	
	rmsth_hole_av = rmsth_hole_av/(Double_t)no_hole_used;
	
	rmsx_hole_av = rmsx_hole_av/no_hole_used;






	// loop to go through per-foil data and calculate average rms of foils (effectively weighs all foils equally)

	for(Int_t j = 0; j<NFiols; j++){
	  
	  //cout << "For hole " << j << " dth_hole = " << dth_hole[j] << endl;
	  //	  cout << "          rmsth_hole = " << rmsth_hole[j] << endl;


	  if( foil_stat[j] > 0){

	    dth_foil[j] = dth_foil[j]/foil_stat[j];     
	    rmsth_foil[j] = rmsth_foil[j]/foil_stat[j];

	    dx_foil[j] = dx_foil[j]/foil_stat[j];
	    rmsx_foil[j] = rmsx_foil[j]/foil_stat[j];
	    no_foil_used++;
	  }
	  else{
	    dth_foil[j] = 0;     
	    rmsth_foil[j] = 0;
	    
	    dx_foil[j] = 0;
	    rmsx_foil[j] = 0;

	  }

	  rmsth_foil_av += rmsth_foil[j];
	  rmsx_foil_av += rmsx_foil[j];
	  


	}

	
	rmsth_foil_av = rmsth_foil_av/(Double_t)no_foil_used;

	rmsx_foil_av = rmsx_foil_av/no_foil_used;





	//	cout << "rmsth_hole_av = " << rmsth_hole_av << endl;



	DEBUG_INFO("SumSquareDTh","#%d : dth = %f,rmsth=%f",NCall,dth/fNRawData,TMath::Sqrt(rmsth/fNRawData));
	DEBUG_INFO("SumSquareDTh average hole deviation","#%d : rmsth=%f",NCall,TMath::Sqrt(rmsth_hole_av));
	DEBUG_INFO("SumSquareDTh average foil deviation","#%d : rmsth=%f",NCall,TMath::Sqrt(rmsth_foil_av));

	
	DEBUG_INFO("SumSquareDX","#%d : dx = %f,rmsx=%f",NCall,dx/fNRawData,TMath::Sqrt(rmsdx/fNRawData));
	DEBUG_INFO("SumSquareDX average hole deviation","#%d : rmsx=%f",NCall,TMath::Sqrt(rmsx_hole_av));
	DEBUG_INFO("SumSquareDX average foil deviation","#%d : rmsth=%f \n",NCall,TMath::Sqrt(rmsx_foil_av));

// 	DEBUG_INFO("VerifyMatrix_Sieve","dphi = %f, rmsphi=%f",
// 			   dphi/fNRawData,TMath::Sqrt(rmsphi/fNRawData));



	return rmsth;
//	return rmsdx;

//	return rmsth_hole_av;


//	return rmsth_foil_av;
	
}

//_____________________________________________________________________________
Double_t LOpticsOpt::SumSquareDPhi(Bool_t PrintEachHole)
{
	//return square sum of diff between calculated tg_ph and expected tg_ph

	Double_t /*dth = 0, */dphi=0;	//Difference
	Double_t /*rmsth = 0, */rmsphi=0; //mean square

	// same parameters but for sieve y
	Double_t dy = 0; //difference
	Double_t rmsdy = 0; // mean square



	// parameters for each individual hole 

	Double_t dph_hole[NHoles] = {0}; 
	Double_t rmsph_hole[NHoles] = {0};
	Double_t hole_stat[NHoles] = {0};
	Double_t rmsph_hole_av = 0;

	Double_t no_hole_used = 0;

	Double_t dy_hole[NHoles] = {0}; 
	Double_t rmsy_hole[NHoles] = {0}; 
	Double_t rmsy_hole_av = 0;


	// parameters for foil optimisation

	Double_t dph_foil[NFiols] = {0}; 
	Double_t rmsph_foil[NFiols] = {0};
	Double_t foil_stat[NFiols] = {0};
	Double_t rmsph_foil_av = 0;

	Double_t no_foil_used = 0;

	Double_t dy_foil[NFiols] = {0}; 
	Double_t rmsy_foil[NFiols] = {0}; 
	Double_t rmsy_foil_av = 0;



	static UInt_t NCall = 0;
	NCall++;

	Double_t /*theta, */phi/*, dp, p, pathl*/;

	// stat. for each hole
	Double_t * stat_dev = NULL;
	Double_t * stat_rms = NULL;
	Int_t * stat_cnt = NULL;
	if (PrintEachHole)
	{
		stat_dev = new Double_t[kMaxDataGroup];
		stat_rms = new Double_t[kMaxDataGroup];
		stat_cnt = new Int_t[kMaxDataGroup];

		for(int i=0; i<kMaxDataGroup; i++) stat_dev[i] = stat_rms[i] = stat_cnt[i] = 0;
	}

	for (UInt_t idx = 0; idx<fNRawData; idx++)
	{
		EventData &eventdata = fRawData[idx];

		Double_t x_fp = eventdata.Data[kX];
		const Double_t (*powers)[5] = eventdata.powers;

  // calculate the matrices we need
// 		CalcMatrix(x_fp, fDMatrixElems);
// 		CalcMatrix(x_fp, fTMatrixElems);
// 		CalcMatrix(x_fp, fYMatrixElems);
// 		CalcMatrix(x_fp, fYTAMatrixElems);
		CalcMatrix(x_fp, fPMatrixElems);
// 		CalcMatrix(x_fp, fPTAMatrixElems);

  // calculate the coordinates at the target
// 		theta = CalcTargetVar(fTMatrixElems, powers);
		phi = CalcTargetVar(fPMatrixElems, powers)/*+CalcTargetVar(fPTAMatrixElems,powers)*/;
// 		y = CalcTargetVar(fYMatrixElems, powers)+CalcTargetVar(fYTAMatrixElems,powers);

  // calculate momentum
// 		dp = CalcTargetVar(fDMatrixElems, powers);

// 		dth += theta - eventdata.Data[kRealTh];
// 		rmsth += (theta - eventdata.Data[kRealTh])*(theta - eventdata.Data[kRealTh]);


		if (PrintEachHole)
		{
		        // JW: FIX this: just replaced kCutID with kFoilID
			UInt_t HoleID = (UInt_t)eventdata.Data[kFoilID];
			assert(HoleID<kMaxDataGroup);

			stat_dev[HoleID] += phi - eventdata.Data[kRealPhi];
			stat_rms[HoleID] += (phi - eventdata.Data[kRealPhi])*(phi - eventdata.Data[kRealPhi]);
			stat_cnt[HoleID] += 1;
		}

		//save the results
		eventdata.Data[kCalcPh] = phi;

#ifdef DISABLED_SIEVE_HOLE
		// JW: altered following lines
		// UInt_t res = (UInt_t)eventdata.Data[kCutID];
		// const UInt_t FoilID = res/(NSieveRow*NSieveCol); //starting 0!
		const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID]; //starting 0!
		// res = res%(NSieveRow*NSieveCol);
		const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
		const UInt_t Row = (UInt_t)eventdata.Data[kRowID]; //starting 0!

		if (DISABLED_SIEVE_HOLE) continue;

#endif

		const UInt_t HoleID = Get_Hole(Col,Row);
		


		dphi += phi - eventdata.Data[kRealPhi];
		rmsphi += (phi - eventdata.Data[kRealPhi])*(phi - eventdata.Data[kRealPhi]);


		// save per-hole results

		dph_hole[HoleID] += phi - eventdata.Data[kRealPhi];
		rmsph_hole[HoleID] += (phi - eventdata.Data[kRealPhi])*(phi - eventdata.Data[kRealPhi]);

		hole_stat[HoleID]++;


		// save per-foil results
		dph_foil[FoilID] += phi - eventdata.Data[kRealPhi];

		rmsph_foil[FoilID] += (phi - eventdata.Data[kRealPhi])*(phi - eventdata.Data[kRealPhi]);

		foil_stat[FoilID]++;






		double ProjectionY = eventdata.Data[kRealTgY] + tan(eventdata.Data[kCalcPh]) * (ZPos);

		// const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
		// const UInt_t Row = (UInt_t)eventdata.Data[kRowID]; //starting 0!

		TVector3 SieveHoleTCS;
		
		if(Col > 0 && Col < NSieveCol && Row > 0 && Row < NSieveRow){

		  SieveHoleTCS = GetSieveHoleTCS(Col,Row);
		}
		else{
		  SieveHoleTCS = GetSieveHoleTCS(1,1);
		}


		const Double_t posy = SieveHoleTCS.Y();

		
		//		const Double_t posy = GetSieveHoleTCS(Col,Row).Y();
	
		dy += ProjectionY - posy;
		rmsdy += (ProjectionY - posy)*(ProjectionY - posy);


		// per-hole results
		dy_hole[HoleID] += ProjectionY - posy;
		rmsy_hole[HoleID] += (ProjectionY - posy)*(ProjectionY - posy);
		

		// per-foil results
		dy_foil[FoilID] += ProjectionY - posy;
		rmsy_foil[FoilID] += (ProjectionY - posy)*(ProjectionY - posy);



	}

// 	DEBUG_INFO("SumSquareDTh","#%d : dth = %f,rmsth=%f",NCall,
// 			   dth/fNRawData,TMath::Sqrt(rmsth/fNRawData));
	DEBUG_INFO("SumSquareDPhi","#%d : dphi = %f, rmsphi=%f",NCall,dphi/fNRawData,TMath::Sqrt(rmsphi/fNRawData));
	DEBUG_INFO("SumSquareDY","#%d : dy = %f, rmsy=%f",NCall,dy/fNRawData,TMath::Sqrt(rmsdy/fNRawData));

	if (PrintEachHole)
	{
		const Double_t PrintLimit = 3e-4; //rad
		DEBUG_INFO("SumSquareDPhi","Print Deviation and RMS of each hole (|Dev| > %f):",PrintLimit);
		cout<<"Config\tCol\tRow\tStat.\tDeviation\tRMS\n";

		Double_t sum_dev=0, sum_absdev=0, sum_rms=0, sum_cnt=0, weight=0, weighted_dev=0;
		Double_t sum_devabs[NFiols] = {0}, sum_nhole[NFiols]={0}, sum_devabserr[NFiols]={0};
		for(UInt_t HoleID=0; HoleID<kMaxDataGroup; HoleID++)
		{
			if (stat_cnt[HoleID] <= 0) continue;

			UInt_t res = HoleID;
			const UInt_t FoilID = res/(NSieveRow*NSieveCol); //starting 0!
			assert(FoilID<NFiols);
			res = res%(NSieveRow*NSieveCol);
			const UInt_t Col = res/(NSieveRow); //starting 0!
			const UInt_t Row = res%(NSieveRow); //starting 0!

//			Col, Row;

			Double_t dev = stat_dev[HoleID]/stat_cnt[HoleID];
			Double_t rms  = TMath::Sqrt(stat_rms[HoleID]/stat_cnt[HoleID]);

//			if (TMath::Abs(dev)-rms/TMath::Sqrt(stat_cnt[HoleID]) > PrintLimit)
//			{
//				cout<<FoilID<<"\t"<<Col<<"\t"<<Row<<"\t"
//						<<stat_cnt[HoleID]<<"\t"
//						<<dev <<" +/- "<<rms/TMath::Sqrt(stat_cnt[HoleID])<<"\t"
//						<<rms <<"\n";
//			}

			sum_cnt ++;
			sum_dev += dev;
			sum_rms += rms;
			sum_absdev += TMath::Abs(dev);

			weight += stat_cnt[HoleID] / rms /rms;
			weighted_dev += TMath::Abs(dev) * stat_cnt[HoleID] / rms /rms ;

			sum_devabs[FoilID] += TMath::Abs(dev);
			sum_nhole[FoilID] += 1;
			sum_devabserr[FoilID] += rms*rms/ stat_cnt[HoleID];
		}

		cout<<"All(Same wt each hole):\t"<<(int)sum_cnt<<"\t"
				<<sum_dev/sum_cnt<<"\t"<<sum_rms/sum_cnt<<"\n";


//		DEBUG_INFO("SumSquareDPhi","Average |deviation| = %f +/- %f (stat. w.)"
//				,weighted_dev/weight,TMath::Sqrt(1/weight));
		DEBUG_INFO("SumSquareDPhi","Average |deviation| = %f "
				,sum_absdev/sum_cnt);

		TString tmp = "Average |deviation| =\n ";
		for (UInt_t i = 0; i< NFiols; i++)
			tmp += Form("#%d: %f +/- %f \n ",i ,sum_devabs[i]/sum_nhole[i],TMath::Sqrt(sum_devabserr[i])/sum_nhole[i]);
		DEBUG_INFO("SumSquareDPhi",tmp.Data());

//		DEBUG_INFO("SumSquareDPhi","Average |deviation| =\n %f +/- %f \n %f +/- %f \n %f +/- %f \n %f +/- %f \n %f +/- %f \n %f +/- %f \n %f +/- %f ",
//				sum_devabs[0]/sum_nhole[0],TMath::Sqrt(sum_devabserr[0])/sum_nhole[0],
//				sum_devabs[1]/sum_nhole[1],TMath::Sqrt(sum_devabserr[1])/sum_nhole[1],
//				sum_devabs[2]/ssum_nhole[2],TMath::Sqrt(sum_devabserr[2])/sum_nhole[2],
//				sum_devabs[3]/sum_nhole[3],TMath::Sqrt(sum_devabserr[3])/sum_nhole[3],
//				sum_devabs[4]/sum_nhole[4],TMath::Sqrt(sum_devabserr[4])/sum_nhole[4],
//				sum_devabs[5]/sum_nhole[5],TMath::Sqrt(sum_devabserr[5])/sum_nhole[5],
//				sum_devabs[6]/sum_nhole[6],TMath::Sqrt(sum_devabserr[6])/sum_nhole[6]
//				                         );

		delete [] stat_dev;
		delete [] stat_rms;
		delete [] stat_cnt;
	}





	fstream col_test;
	col_test.open("Foil5_col_test_1_2_all_foil", ios_base::out);
	  
	
	// loop to go through individual hole data

	for(Int_t j = 0; j<NHoles; j++){

	  //cout << "For hole " << j << " dth_hole = " << dth_hole[j] << endl;
	  //	  cout << "          rmsth_hole = " << rmsth_hole[j] << endl;
	  
	  

	  if( hole_stat[j] > 0){
	    
	    dph_hole[j] = dph_hole[j]/hole_stat[j];     
	    rmsph_hole[j] = rmsph_hole[j]/hole_stat[j];
	    
	    dy_hole[j] = dy_hole[j]/hole_stat[j];
	    rmsy_hole[j] = rmsy_hole[j]/hole_stat[j];
	    no_hole_used++;

	    
	    col_test <<  "Hole " << j << " : " <<  "rmsph = " << rmsph_hole[j] << ", rmsy = " <<  rmsy_hole[j] << endl;

	    
	  }
	  else{
	    dph_hole[j] = 0;     
	    rmsph_hole[j] = 0;
	    
	    dy_hole[j] = 0;
	    rmsy_hole[j] = 0;

	  }
	  
	  rmsph_hole_av += rmsph_hole[j];
	  rmsy_hole_av += rmsy_hole[j];
	  

	  
	}

	

	// loop to go through per-foil data and calculate average rms of foils (effectively weighs all foils equally)


	for(Int_t j = 0; j<NFiols; j++){
	  
	  //cout << "For hole " << j << " dph_hole = " << dph_hole[j] << endl;
	  //	  cout << "          rmsph_hole = " << rmsph_hole[j] << endl;

	  
	  if( foil_stat[j] > 0){

	    dph_foil[j] = dph_foil[j]/foil_stat[j];     
	    rmsph_foil[j] = rmsph_foil[j]/foil_stat[j];

	    dy_foil[j] = dy_foil[j]/foil_stat[j];
	    rmsy_foil[j] = rmsy_foil[j]/foil_stat[j];
	    no_foil_used++;
	  }
	  else{
	    dph_foil[j] = 0;     
	    rmsph_foil[j] = 0;
	    
	    dy_foil[j] = 0;
	    rmsy_foil[j] = 0;

	  }

	  rmsph_foil_av += rmsph_foil[j];
	  rmsy_foil_av += rmsy_foil[j];
	  


	}

	
	rmsph_foil_av = rmsph_foil_av/(Double_t)no_foil_used;

	rmsy_foil_av = rmsy_foil_av/no_foil_used;



	DEBUG_INFO("SumSquareDPhi average foil deviation","#%d : rmsphi=%f",NCall,TMath::Sqrt(rmsph_foil_av));
	DEBUG_INFO("SumSquareDY average foil deviation","#%d :  rmsy=%f",NCall,TMath::Sqrt(rmsy_foil_av));


	rmsph_hole_av = rmsph_hole_av/(Double_t)no_hole_used;

	rmsy_hole_av = rmsy_hole_av/no_hole_used;


	DEBUG_INFO("SumSquareDPhi average hole deviation","#%d : rmsphi=%f",NCall,TMath::Sqrt(rmsph_hole_av));
	DEBUG_INFO("SumSquareDY average hole deviation","#%d :  rmsy=%f \n",NCall,TMath::Sqrt(rmsy_hole_av));



	
	col_test << "Overall rmsphi = " << rmsphi << endl;
	col_test << "Overall rmsdy = " << rmsdy << endl;


	return rmsphi;
	//	return  rmsdy;
	//	return rmsph_hole_av;
	// return rmsph_foil_av;

}

//_____________________________________________________________________________

Double_t LOpticsOpt::SumSquareDTgY(void)
{

	//return square sum of diff between calculated tg_y and expected tg_y

	Double_t /*dth = 0, */dtg_y=0;	//Difference
	Double_t /*rmsth = 0, */dtg_y_rms=0; //mean square


	// add differene between calculated and 'real' reactz


	Double_t /*dth = 0, */dz=0;	//Difference
	Double_t /*rmsth = 0, */dz_rms=0; //mean square

	Double_t /*dth = 0, */dz_calc=0;	//Difference
	Double_t /*rmsth = 0, */dz_calc_rms=0; //mean square


	// parameters for equal-foil-weighted optimisation

	
	Double_t dz_foil[NFiols] = {0}; 
	Double_t rmsz_foil[NFiols] = {0};
	Double_t foil_stat[NFiols] = {0};
	Double_t rmsz_foil_av = 0;

	Double_t no_foil_used = 0;

	static UInt_t NCall = 0;
	NCall++;

	Double_t /*theta, */tg_y/*, dp, p, pathl*/;

	for (UInt_t idx = 0; idx<fNRawData; idx++)
	{
		EventData &eventdata = fRawData[idx];

		Double_t x_fp = eventdata.Data[kX];
		const Double_t (*powers)[5] = eventdata.powers;

  // calculate the matrices we need
// 		CalcMatrix(x_fp, fDMatrixElems);
// 		CalcMatrix(x_fp, fTMatrixElems);
		CalcMatrix(x_fp, fYMatrixElems);
// 		CalcMatrix(x_fp, fYTAMatrixElems);
// 		CalcMatrix(x_fp, fPMatrixElems);
// 		CalcMatrix(x_fp, fPTAMatrixElems);

  // calculate the coordinates at the target
// 		theta = CalcTargetVar(fTMatrixElems, powers);
// 		phi = CalcTargetVar(fPMatrixElems, powers)/*+CalcTargetVar(fPTAMatrixElems,powers)*/;
		tg_y = CalcTargetVar(fYMatrixElems, powers)
				/*+CalcTargetVar(fYTAMatrixElems,powers)*/;

  // calculate momentum
// 		dp = CalcTargetVar(fDMatrixElems, powers);
		
		const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID];

		

		//const Double_t ArbitaryVertexShift = fArbitaryVertexShift[FoilID] * (-TMath::Sin(HRSAngle));

		const Double_t ArbitaryVertexShift = 0;



		// JW: added debug statements for dtg_y
		
		// if(NCall==5 && eventdata.Data[kRealTgY]*eventdata.Data[kRealTgY]>1 ){
		//   DEBUG_INFO("SumSquareDTgY","#%d : tg_y = %f, Real-TgY=%f, ArbitaryVertexShift=%f",NCall,
		// 	     tg_y,eventdata.Data[kRealTgY],ArbitaryVertexShift);
		// }



		// calculate reactz
		
		const Int_t a = (HRSAngle > 0) ? 1 : -1;

		


		


		dtg_y += tg_y - eventdata.Data[kRealTgY]+ArbitaryVertexShift;
		dtg_y_rms += (tg_y - eventdata.Data[kRealTgY]+ArbitaryVertexShift)*(tg_y - eventdata.Data[kRealTgY]+ArbitaryVertexShift);


		// if(dtg_y > 0.5){
		//   DEBUG_INFO("SumSquareDTgY","#%d : tg_y = %f, Real-TgY=%f, ArbitaryVertexShift=%f & dtg_y=%f",NCall,tg_y,eventdata.Data[kRealTgY],ArbitaryVertexShift,dtg_y);
		// }

		//save the results
		eventdata.Data[kCalcTgY] = tg_y;

		
		TVector3 BeamSpotHCS(eventdata.Data[kBeamX], eventdata.Data[kBeamY], targetfoils[FoilID]);


		//		Double_t reactz = - ( eventdata.Data[kCalcTgY] -a*MissPointZ)*TMath::Cos(TMath::ATan(eventdata.Data[kRealPhi]))/TMath::Sin(HRSAngle + TMath::ATan(eventdata.Data[kRealPhi]))    +          BeamSpotHCS.X()*TMath::Cos(HRSAngle+TMath::ATan(eventdata.Data[kRealPhi]))/TMath::Sin(HRSAngle + TMath::ATan(eventdata.Data[kRealPhi]));

		Double_t reactz_calc = - ( eventdata.Data[kCalcTgY] -a*MissPointZ)*TMath::Cos(TMath::ATan(eventdata.Data[kCalcPh]))/TMath::Sin(HRSAngle + TMath::ATan(eventdata.Data[kCalcPh]))    +          BeamSpotHCS.X()*TMath::Cos(HRSAngle+TMath::ATan(eventdata.Data[kCalcPh]))/TMath::Sin(HRSAngle + TMath::ATan(eventdata.Data[kCalcPh]));


		Double_t Real_Tg_Phi = eventdata.Data[kRealPhi];

		// def found in Sean's code and adapted
		Double_t reactz = - ( eventdata.Data[kCalcTgY] -a*MissPointZ)*TMath::Cos(TMath::ATan(Real_Tg_Phi))/TMath::Sin(HRSAngle + TMath::ATan(Real_Tg_Phi)) + BeamSpotHCS.X()*TMath::Cos(HRSAngle+TMath::ATan(Real_Tg_Phi))/TMath::Sin(HRSAngle+TMath::ATan(Real_Tg_Phi));




		// re-write reactz

		//		reactz = 
		  // (  - ( eventdata.Data[kCalcTgY] + MissPointY)   +   BeamSpotHCS.X()*TMath::Cos(HRSAngle+TMath::ATan(eventdata.Data[kRealPhi]))  )
		  // 		/TMath::Sin(HRSAngle + TMath::ATan(eventdata.Data[kRealPhi]));

		dz += reactz - eventdata.Data[kRealReactZ];

		dz_rms += (reactz - eventdata.Data[kRealReactZ])*(reactz - eventdata.Data[kRealReactZ]);


		dz_calc += reactz_calc - eventdata.Data[kRealReactZ];

		dz_calc_rms += (reactz_calc - eventdata.Data[kRealReactZ])*(reactz_calc - eventdata.Data[kRealReactZ]);



		// equal foil-weighted calc

		dz_foil[FoilID] += reactz - eventdata.Data[kRealReactZ];
		rmsz_foil[FoilID]  += (reactz - eventdata.Data[kRealReactZ])*(reactz - eventdata.Data[kRealReactZ]);

		foil_stat[FoilID]++;



	}



	
	// loop to go through per-foil data and calculate average rms of foils (effectively weighs all foils equally)


	for(Int_t j = 0; j<NFiols; j++){
	  
	  	  
	  if( foil_stat[j] > 0){

	    dz_foil[j] = dz_foil[j]/foil_stat[j];     
	    rmsz_foil[j] = rmsz_foil[j]/foil_stat[j];

	    no_foil_used++;
	  }
	  else{
	    dz_foil[j] = 0;     
	    rmsz_foil[j] = 0;	   

	  }

	  rmsz_foil_av += rmsz_foil[j];
	}

	
	rmsz_foil_av = rmsz_foil_av/(Double_t)no_foil_used;





// 	DEBUG_INFO("SumSquareDTh","#%d : dth = %f,rmsth=%f",NCall,
// 			   dth/fNRawData,TMath::Sqrt(rmsth/fNRawData));
	DEBUG_INFO("SumSquareDTgY","#%d : dtg_y = %f, dtg_y_rms=%f",NCall,dtg_y/fNRawData,TMath::Sqrt(dtg_y_rms/fNRawData));

	DEBUG_INFO("SumSquareDTgY","#%d : dz = %f, dz_rms=%f",NCall,dz/fNRawData,TMath::Sqrt(dz_rms/fNRawData));

	DEBUG_INFO("SumSquareDTgY","#%d : dz_calc = %f, dz_calc_rms=%f",NCall,dz_calc/fNRawData,TMath::Sqrt(dz_calc_rms/fNRawData));

	DEBUG_INFO("SumSquareDTgY","#%d  dz_foil_weighted_rms=%f",NCall,TMath::Sqrt(rmsz_foil_av/fNRawData));

	//	return dtg_y_rms;

	return dz_rms;

	// return dz_calc_rms; 


	//	return rmsz_foil_av;
}

//_____________________________________________________________________________
Double_t LOpticsOpt::SumSquareDTgYAverFoils(void)
{

	//return square sum of diff between calculated tg_y and expected tg_y
	//Statistical Weight on each foil is same

	Double_t /*dth = 0, */dtg_y=0;	//Difference
	Double_t /*rmsth = 0, */dtg_y_rms=0; //mean square

	static UInt_t NCall = 0;
	NCall++;

	Double_t /*theta, */tg_y/*, dp, p, pathl*/;

	Double_t dtg_y_foil[NFiols] ={0};
	Double_t rmstg_y_foil[NFiols] ={0};
	UInt_t ndata_foil[NFiols] ={0};


	for (UInt_t idx = 0; idx<fNRawData; idx++)
	{
		EventData &eventdata = fRawData[idx];

		Double_t x_fp = eventdata.Data[kX];
		const Double_t (*powers)[5] = eventdata.powers;

  // calculate the matrices we need
// 		CalcMatrix(x_fp, fDMatrixElems);
// 		CalcMatrix(x_fp, fTMatrixElems);
		CalcMatrix(x_fp, fYMatrixElems);
// 		CalcMatrix(x_fp, fYTAMatrixElems);
// 		CalcMatrix(x_fp, fPMatrixElems);
// 		CalcMatrix(x_fp, fPTAMatrixElems);

  // calculate the coordinates at the target
// 		theta = CalcTargetVar(fTMatrixElems, powers);
// 		phi = CalcTargetVar(fPMatrixElems, powers)/*+CalcTargetVar(fPTAMatrixElems,powers)*/;
		tg_y = CalcTargetVar(fYMatrixElems, powers)/*+CalcTargetVar(fYTAMatrixElems,powers)*/;

  // calculate momentum
// 		dp = CalcTargetVar(fDMatrixElems, powers);

		dtg_y += tg_y - eventdata.Data[kRealTgY];
		dtg_y_rms += (tg_y - eventdata.Data[kRealTgY])*(tg_y - eventdata.Data[kRealTgY]);

		//save the results
		eventdata.Data[kCalcTgY] = tg_y;


		//statics by foils
		const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID];
		const Double_t ArbitaryVertexShift =		//arbitary shifts
				fArbitaryVertexShift[FoilID] * (-TMath::Sin(HRSAngle));
		ndata_foil[FoilID]++;
		dtg_y_foil[FoilID] += eventdata.Data[kCalcTgY] - eventdata.Data[kRealTgY]+ArbitaryVertexShift;
		rmstg_y_foil[FoilID] +=
				(eventdata.Data[kCalcTgY] - eventdata.Data[kRealTgY]+ArbitaryVertexShift)*
				(eventdata.Data[kCalcTgY] - eventdata.Data[kRealTgY]+ArbitaryVertexShift);
	}

	Double_t /*dth = 0, */dtg_y_foilaver=0;	//Difference
	Double_t /*rmsth = 0, */dtg_y_rms_foilaver=0; //mean square
	for(UInt_t FoilID = 0; FoilID<NFiols; FoilID++)
	{
		assert(ndata_foil[FoilID]);// make sure there are at least one event on this foil
		const UInt_t N = ndata_foil[FoilID];

		dtg_y_foilaver += dtg_y_foil[FoilID]/N;
		dtg_y_rms_foilaver += rmstg_y_foil[FoilID]/N;;
	}
	dtg_y_foilaver = dtg_y_foilaver/NFiols*fNRawData;
	dtg_y_rms_foilaver = dtg_y_rms_foilaver/NFiols*fNRawData;

	DEBUG_INFO("SumSquareDTgY","#%d Foil Ave: dtg_y = %f, dtg_y_rms=%f",NCall,
			   dtg_y_foilaver/fNRawData,TMath::Sqrt(dtg_y_rms_foilaver/fNRawData));
	DEBUG_INFO("SumSquareDTgY","#%d : dtg_y = %f, dtg_y_rms=%f",NCall,
			   dtg_y/fNRawData,TMath::Sqrt(dtg_y_rms/fNRawData));

	return dtg_y_rms_foilaver;
}

//_____________________________________________________________________________

void LOpticsOpt::PrepareSieve(void)
{
	//calculate kRealTh, kRealPhi

// 	DEBUG_INFO("PrepareSieve","Entry Point");

	Double_t dth = 0, dphi=0;
	Double_t exttargcorr_th=0, rms_exttargcorr_th=0;




	// JW: attempt at plot of realth and realph (to check if they are reasonable)


	TH2D* Real_angles[NFiols];


	for(UInt_t idx = 0; idx<NFiols; idx++){
		  
	  Real_angles[idx]  = new TH2D(Form("Real_angles_%d",idx),Form("Real_angles_%d",idx),500,-0.04,0.04,500,-0.08,0.08);
	}



	for (UInt_t idx = 0; idx<fNRawData; idx++)
	{
		EventData &eventdata = fRawData[idx];

		// JW: adjusted
		const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID]; //starting 0!
		// res = res%(NSieveRow*NSieveCol);
		const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
		const UInt_t Row = (UInt_t)eventdata.Data[kRowID]; //starting 0!

		assert(FoilID<NFiols);//check array index size

		//		cout << "Col = " << Col << endl;


		

		TVector3 SieveHoleTCS;
		
		
		// below gives default vales for row and column if only foils have been cut (and thus sieve hole row and column unknown)
		// this makes the information about kRealTh and kRealPh meaningless
		if(Col > 0 && Col < NSieveCol && Row > 0 && Row < NSieveRow){
		  //		  cout << "entered default loop!" << endl;
		  //		  cout << "Col = " << Col << " and Row = " << Row << endl;
		  SieveHoleTCS = GetSieveHoleTCS(Col,Row);
		}
		else{
		  SieveHoleTCS = GetSieveHoleTCS(1,1);

		}

		//		const TVector3 SieveHoleTCS = GetSieveHoleTCS(Col,Row);

		const TVector3 BeamSpotHCS(eventdata.Data[kBeamX],eventdata.Data[kBeamY],targetfoils[FoilID]);

		const TVector3 BeamSpotTCS=fTCSInHCS.Inverse()*(BeamSpotHCS-fPointingOffset);

		


		const TVector3 MomDirectionTCS = SieveHoleTCS - BeamSpotTCS;

		eventdata.Data[kRealTh] = MomDirectionTCS.X()/MomDirectionTCS.Z();
		eventdata.Data[kRealPhi] = MomDirectionTCS.Y()/MomDirectionTCS.Z();


		// if(FoilID != 1){
		//   Real_angles[FoilID]->Fill(eventdata.Data[kRealPhi],eventdata.Data[kRealTh]);
		// }


		const Double_t x_tg = BeamSpotTCS.X()-BeamSpotTCS.Z()*tan(eventdata.Data[kRealTh]); // default optimiser calc for x_tg


		// new x_tg definition taken from sagdh tech note

		//		const Double_t x_tg = - eventdata.Data[kRealTh]*(targetfoils[FoilID] * TMath::Cos(HRSAngle) )/( TMath::Cos( TMath::ATan(eventdata.Data[kRealPhi]) ) ) -  eventdata.Data[kBeamY];   

		

		const Double_t y_tg = BeamSpotTCS.Y() - BeamSpotTCS.Z() * tan(eventdata.Data[kRealPhi]);
		

		// additional calculations necessary for x_tg used by E97-110 experiment with septum detailed in technote (can be found as tech-note 8 here: https://hallaweb.jlab.org/experiment/E97-110/tech/tech_notes.html)

		// TVector3 Tg_YSpotTCS(0,eventdata.Data[kRealTgY],0);
		// TVector3 MomDirectionTCS_Z(0,eventdata.Data[kL_tr_tg_ph],1);

		// TVector3 Tg_YSpotHCS=fTCSInHCS*Tg_YSpotTCS+fPointingOffset;
		// TVector3 MomDirectionHCS=fTCSInHCS*MomDirectionTCS_Z;
		
		// Double_t reactz = (Tg_YSpotHCS.X() - eventdata.Data[kBeamX] - (MomDirectionHCS.X()/MomDirectionHCS.Z())*Tg_YSpotHCS.Z() )/( eventdata.Data[kBeamDirX]/eventdata.Data[kBeamDirZ] - (MomDirectionHCS.X()/MomDirectionHCS.Z()) );


		// Double_t x_tg = (-eventdata.Data[kRealTh] * reactz * TMath::Cos(HRSAngle) )/ ( TMath::Cos( TMath::ATan(eventdata.Data[kRealPhi])) );



		
		eventdata.Data[kRealTgX] = x_tg;
		eventdata.Data[kRealTgY] = y_tg;


		//Expected th ph before ext. target correction
// 		fDeltaTh = fThetaCorr * x_tg;
// 		Double_t theta = trkifo->GetTheta() + fDeltaTh;
		eventdata.Data[kRealThMatrix]=eventdata.Data[kRealTh] - x_tg*ExtTarCor_ThetaCorr;


//		eventdata.Data[kRealThMatrix]=eventdata.Data[kRealTh];

		Real_angles[FoilID]->Fill(eventdata.Data[kRealPhi],eventdata.Data[kRealThMatrix]);
		


			exttargcorr_th += x_tg*ExtTarCor_ThetaCorr;
		rms_exttargcorr_th += x_tg*ExtTarCor_ThetaCorr*x_tg*ExtTarCor_ThetaCorr;


		DEBUG_MASSINFO("PrepareSieve","%d,%d,%d: D_Th = %f,\t D_Phi = %f",FoilID,Col,Row,
					   eventdata.Data[kRealThMatrix]-eventdata.Data[kL_tr_tg_th],
		eventdata.Data[kRealPhi]-eventdata.Data[kL_tr_tg_ph]
					  );
		DEBUG_MASSINFO("PrepareSieve","%f,\t%f",
					   eventdata.Data[kRealThMatrix],eventdata.Data[kL_tr_tg_th]
					  );

		dth+=eventdata.Data[kRealThMatrix]-eventdata.Data[kL_tr_tg_th];
		dphi+=eventdata.Data[kRealPhi]-eventdata.Data[kL_tr_tg_ph];

	}

	DEBUG_INFO("PrepareSieve","Average : D_Th = %f,\t D_Phi = %f",dth/fNRawData,dphi/fNRawData);
	DEBUG_INFO("PrepareSieve","Average Extended Target Corretion: th = %f,\t rms_th = %f"
			,exttargcorr_th/fNRawData,TMath::Sqrt(rms_exttargcorr_th/fNRawData));

	TCanvas* c1 = new TCanvas("c1","c1");

	
	c1->Divide(3);

	for(UInt_t idx = 0; idx<NFiols; idx++){
	  
	  c1->cd(idx+1);
	  Real_angles[idx]->Draw("colz");

	}

	//make sure kCalcTh, kCalcPh is filled
	SumSquareDTh();
	SumSquareDPhi();
}
//_____________________________________________________________________________

void LOpticsOpt::PrepareVertex(void)
{
	//calculate kRealTgY, kRealReactZ

	//set fYMatrixElems as current matrix to optimize

	static UInt_t NCall = 0;
	NCall++;



	fCurrentMatrixElems = & fYMatrixElems;

	Double_t dtg_y = 0,dtg_y_rms = 0;

	for (UInt_t idx = 0; idx<fNRawData; idx++)
	{
		EventData &eventdata = fRawData[idx];

		// UInt_t res = (UInt_t)eventdata.Data[kCutID];
		// const UInt_t FoilID = res;
		
		const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID];

		assert(FoilID<NFiols);//check array index size


		//		TVector3 BeamSpotHCS(eventdata.Data[kBeamX] + targetfoils[FoilID]/eventdata.Data[kBeamDirZ]*eventdata.Data[kBeamDirX],eventdata.Data[kBeamY] + targetfoils[FoilID]/eventdata.Data[kBeamDirZ]*eventdata.Data[kBeamDirY],targetfoils[FoilID]);

		//TVector3 BeamSpotHCS(eventdata.Data[kBeamX] + targetfoils[FoilID]*(eventdata.Data[kBeamDirX]/eventdata.Data[kBeamDirX]), eventdata.Data[kBeamY] + targetfoils[FoilID]*(eventdata.Data[kBeamDirY]/eventdata.Data[kBeamDirZ]),targetfoils[FoilID]);

		


		TVector3 BeamSpotHCS(eventdata.Data[kBeamX], eventdata.Data[kBeamY], targetfoils[FoilID]);


		TVector3 BeamSpotTCS= fTCSInHCS.Inverse() * (BeamSpotHCS-fPointingOffset);



		const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
		const UInt_t Row = (UInt_t)eventdata.Data[kRowID];

		TVector3 SieveHoleTCS = GetSieveHoleTCS(Col,Row);
		

		TVector3 MomDirectionTCS = SieveHoleTCS - BeamSpotTCS;

		Double_t Real_Tg_Phi = MomDirectionTCS.Y() / MomDirectionTCS.Z();
				
		
		// original tg y calculation
		//		Double_t Real_Tg_Y = BeamSpotTCS.Y() + eventdata.Data[kL_tr_tg_ph] * (0-BeamSpotTCS.Z());

		// copied from Sean's version 


		//		Double_t Real_Tg_Y = BeamSpotTCS.Y() + Real_Tg_Phi * (0 - BeamSpotTCS.Z());

		// JW: test
		
		Double_t Real_Tg_Y = BeamSpotTCS.Y() + TMath::Tan(Real_Tg_Phi) * (0 - BeamSpotTCS.Z());

		// cout << "Real_Tg_Y calculation: " << endl;

		// cout << "BeamSpotTCS.Y() (" << BeamSpotTCS.Y() << " + eventdata.Data[kL_tr_tg_ph] (" << eventdata.Data[kL_tr_tg_ph] <<  ") * (0-BeamSpotTCS.Z() ( " << BeamSpotTCS.Z()  <<  " )  )  = " << Real_Tg_Y << endl;
		

		eventdata.Data[kRealTgY] = Real_Tg_Y;

		if(NCall <5 && eventdata.Data[kRealTgY]*eventdata.Data[kRealTgY]>1 ){
		  DEBUG_INFO("PrepareVertex","#%d : Real_tg_y = %f, BeamSpotTCS.Y() = %f, eventdata.Data[kL_tr_tg_ph] = %f, 0-BeamSpotTCS.Z() = %f",NCall,eventdata.Data[kRealTgY], BeamSpotTCS.Y(), eventdata.Data[kL_tr_tg_ph],  (0-BeamSpotTCS.Z()));
		}

		eventdata.Data[kRealReactZ] = targetfoils[FoilID];

		dtg_y += (eventdata.Data[kL_tr_tg_y] - Real_Tg_Y);
		dtg_y_rms += (eventdata.Data[kL_tr_tg_y] - Real_Tg_Y)
				*(eventdata.Data[kL_tr_tg_y] - Real_Tg_Y);

		////////////////////////////////////////////////////////////
		//redundant checks
		////////////////////////////////////////////////////////////

		TVector3 Tg_YSpotTCS(0,eventdata.Data[kRealTgY],0);
		//		TVector3 MomDirectionTCS(0,eventdata.Data[kL_tr_tg_ph],1);

		TVector3 Tg_YSpotHCS=fTCSInHCS*Tg_YSpotTCS+fPointingOffset;
		TVector3 MomDirectionHCS=fTCSInHCS*MomDirectionTCS;

		//		assert(Tg_YSpotHCS.Y()==MissPointY);	// check coordinates conversions
		//		assert(MomDirectionHCS.Y()==0);			// check coordinates conversions

		//		Double_t reactz = Tg_YSpotHCS.Z()
		//				- (Tg_YSpotHCS.X()-eventdata.Data[kBeamX])
		//				/MomDirectionHCS.X()*MomDirectionHCS.Z();
		
		//		Double_t reactz = (Tg_YSpotHCS.X() - eventdata.Data[kBeamX] - (MomDirectionHCS.X()/MomDirectionHCS.Z())*Tg_YSpotHCS.Z() )/( eventdata.Data[kBeamDirX]/eventdata.Data[kBeamDirZ] - (MomDirectionHCS.X()/MomDirectionHCS.Z()) );

		// definition adopoted from Sean's code


		
		const Int_t a = (HRSAngle > 0) ? 1 : -1;

		//		Double_t reactz = - ( eventdata.Data[kCalcTgY] -a*MissPointZ)*TMath::Cos(TMath::ATan(Real_Tg_Phi))/TMath::Sin(HRSAngle + TMath::ATan(Real_Tg_Phi)) + BeamSpotHCS.X()*TMath::Cos(HRSAngle+TMath::ATan(Real_Tg_Phi))/TMath::Sin(HRSAngle+TMath::ATan(Real_Tg_Phi));
		
		// DEBUG_MASSINFO("PrepareVertex","reactz =%f, eventdata.Data[kRealReactZ]=%f",reactz,eventdata.Data[kRealReactZ]);
		// DEBUG_MASSINFO("PrepareVertex","Real_Tg_Y =%f, eventdata.Data[kL_tr_tg_ph]=%f, targetfoils[FoilID]=%f",
		// 	       Real_Tg_Y,eventdata.Data[kL_tr_tg_ph],targetfoils[FoilID]);
		//JW removed assertion
		// assert(TMath::Abs(reactz - eventdata.Data[kRealReactZ])<1e-4); //check internal calculation consistency
		
		////////////////////////////////////////////////////////////
		
		
	}
	
	DEBUG_INFO("PrepareVertex","Average : dtg_y = %f, dtg_y_rms=%f"
		   ,dtg_y/fNRawData,dtg_y_rms/fNRawData);
	
	//make sure kCalcTh, kCalcPh is filled
	SumSquareDTgY();
}
//_____________________________________________________________________________

void LOpticsOpt::PrepareDp(void)
{
	//calate expected dp_kin, dp_kin offsets ....
	//Fill up fRawData[].Data[] kKineID thr kRealDpKin

	//print Central Momentums
	printf("HRSCentralMom[%d] (GeV) = {",NKine);
	for(UInt_t KineID = 0; KineID<NKine; KineID++)
		printf("%f  ",HRSCentralMom[KineID]*1000);
	printf("}\n");

	//print radiation loss numbers
	printf("RadiationLossByFoil[%d] (MeV) = {",NFiols);
	for(UInt_t FoilID = 0; FoilID<NFiols; FoilID++)
		printf("%f  ",RadiationLossByFoil[FoilID]*1000);
	printf("}\n");


	//set fDMatrixElems as current matrix to optimize
	fCurrentMatrixElems = & fDMatrixElems;

	Double_t dth = 0, dphi=0, scatang=0;
	Double_t dpkinoff = 0, dpkinoff_rms=0;
	fNCalibData = 0;
	Double_t exttargcorr_dp=0, rms_exttargcorr_dp=0;

	for (UInt_t idx = 0; idx<fNRawData; idx++)
	{
		DEBUG_MASSINFO("PrepareDp","=========== Event %d ===========",idx);

		EventData &eventdata = fRawData[idx];

		//decoding kCutID
		// JW: FIX just swapped kFoilID for kCutID here
		UInt_t res = (UInt_t)eventdata.Data[kFoilID];
		//		UInt_t res = (UInt_t)eventdata.Data[kCutID];
		


		const UInt_t ExtraDataFlag = res/(NSieveRow*NSieveCol*NFiols*NKine);
		res = res%(NSieveRow*NSieveCol*NFiols*NKine);
		const UInt_t KineID = res/(NSieveRow*NSieveCol*NFiols); //starting 0!
		res = res%(NSieveRow*NSieveCol*NFiols);
		const UInt_t FoilID = res/(NSieveRow*NSieveCol); //starting 0!
		res = res%(NSieveRow*NSieveCol);
		const UInt_t Col = res/(NSieveRow); //starting 0!
		const UInt_t Row = res%(NSieveRow); //starting 0!
		assert(ExtraDataFlag<2);		//Check flag range. beyond 2 is not used
		assert(KineID<NKine);//check array index size
		assert(FoilID<NFiols);//check array index size

		// DEBUG_MASSINFO("PrepareDp","%d => KineID=%d,\tFoilID=%d,\tCol=%d,\tRow=%d",
		// 			   (UInt_t)eventdata.Data[kCutID],KineID,FoilID,Col,Row);

		//write some variables
 		eventdata.Data[kExtraDataFlag] = ExtraDataFlag;
		if (!ExtraDataFlag) fNCalibData++;
		eventdata.Data[kKineID] = KineID;
		eventdata.Data[kCentralp] = HRSCentralMom[KineID];
		const Double_t EPICS_P = eventdata.Data[kL_tr_p]/(1.+eventdata.Data[kL_tr_tg_dp]);
		DEBUG_MASSINFO("PrepareDp","Central_P/GeV: EPICS=%f (%f, %f), NMR=%f",
					   EPICS_P, eventdata.Data[kL_tr_p],eventdata.Data[kL_tr_tg_dp],eventdata.Data[kCentralp]);
		assert(TMath::Abs(EPICS_P - eventdata.Data[kCentralp])<5e-2);  // WARNING: Change to 1e-3 for normal replays

		const TVector3 SieveHoleTCS = GetSieveHoleTCS(Col,Row);

		TVector3 BeamSpotHCS(
							 eventdata.Data[kBeamX],
		eventdata.Data[kBeamY],
  targetfoils[FoilID]
							);

		// Calculate Real Scattering Angles by Sieve Holes cuts
		TVector3 BeamSpotTCS=fTCSInHCS.Inverse()*(BeamSpotHCS-fPointingOffset);

		TVector3 MomDirectionTCS = SieveHoleTCS - BeamSpotTCS;

		eventdata.Data[kRealTh] = MomDirectionTCS.X()/MomDirectionTCS.Z();
		eventdata.Data[kRealPhi] = MomDirectionTCS.Y()/MomDirectionTCS.Z();

		const Double_t x_tg = BeamSpotTCS.X()-BeamSpotTCS.Z()*eventdata.Data[kRealTh];
		eventdata.Data[kRealTgX] = x_tg;

		DEBUG_MASSINFO("PrepareDp","%d,%d,%d: D_Th = %f,\t D_Phi = %f",FoilID,Col,Row,
					   eventdata.Data[kRealThMatrix]-eventdata.Data[kL_tr_tg_th],
						eventdata.Data[kRealPhi]-eventdata.Data[kL_tr_tg_ph]);
		DEBUG_MASSINFO("PrepareDp","RealTh=%f,\tL_tr_tg_th=%f",
					   eventdata.Data[kRealThMatrix],eventdata.Data[kL_tr_tg_th]);
		DEBUG_MASSINFO("PrepareDp","RealPh=%f,\tkL_tr_tg_ph=%f",
					   eventdata.Data[kRealPhi],eventdata.Data[kL_tr_tg_ph]);
		DEBUG_MASSINFO("PrepareDp","SieveHoleY=%f,\tMom.Y=%f,\tMom.Z=%f",
					   SieveHoleTCS.y(),MomDirectionTCS.Y(),MomDirectionTCS.Z());

		dth+=eventdata.Data[kRealThMatrix]-eventdata.Data[kL_tr_tg_th];
		dphi+=eventdata.Data[kRealPhi]-eventdata.Data[kL_tr_tg_ph];

		TVector3 MomDirectionHCS = fTCSInHCS*MomDirectionTCS;
		TVector3 BeamDirection(0,0,1);
		const Double_t ScatteringAngle = BeamDirection.Angle(MomDirectionHCS);
		eventdata.Data[kScatterAngle] = ScatteringAngle;
		scatang+=ScatteringAngle;

		// calculate difference between dp_kin and dp
		// dp_kin + kDpKinOffsets = dp
		const Double_t DM = ExcitationEnergy[KineID];
		const Double_t Ma = GroundNuclearMass;
		const Double_t P0 = eventdata.Data[kurb_e];
		const Double_t DpKinOffsets =
				(ScatMom(DM, Ma, P0, ScatteringAngle)
				-ScatMom(DM, Ma, P0, TMath::Abs(HRSAngle)))/eventdata.Data[kCentralp];
		eventdata.Data[kDpKinOffsets] = DpKinOffsets;

		dpkinoff += DpKinOffsets;
		dpkinoff_rms += DpKinOffsets*DpKinOffsets;

		// calculate kRealDpKin, should be same for same kine settings
		eventdata.Data[kRadiLossDp]
				= RadiationLossByFoil[FoilID] / eventdata.Data[kCentralp];
		eventdata.Data[kRealDpKin] =
				ScatMom(DM, Ma, P0, TMath::Abs(HRSAngle))/eventdata.Data[kCentralp]- 1 -eventdata.Data[kRadiLossDp] ;

		DEBUG_MASSINFO("PrepareDp","ScatterMom=%f,\t Centralp=%f,radloss = %f, ebeam=%f",
				ScatMom(DM, Ma, P0, TMath::Abs(HRSAngle)), eventdata.Data[kCentralp], eventdata.Data[kRadiLossDp], P0);

		//Expected th ph before ext. target correction
// 		fDeltaTh = fThetaCorr * x_tg;
// 		fDeltaDp = x_tg / fDeltaCorr;
// 		Double_t theta = trkifo->GetTheta() + fDeltaTh;
// 		Double_t dp = trkifo->GetDp() + fDeltaDp;
		eventdata.Data[kRealDpKinMatrix]=eventdata.Data[kRealDpKin] - x_tg/ExtTarCor_DeltaCorr;

		exttargcorr_dp+=x_tg/ExtTarCor_DeltaCorr;
		rms_exttargcorr_dp+=(x_tg/ExtTarCor_DeltaCorr)*(x_tg/ExtTarCor_DeltaCorr);

		//calcalculate expected dp_kin for all other exciation states
		for (UInt_t ExcitID = 0; ExcitID<NExcitationStates; ExcitID++)
		{
			assert(kRealDpKinExcitations+ExcitID<kRealTh);//check array index size
			eventdata.Data[kRealDpKinExcitations+ExcitID] =
					ScatMom(ExcitationEnergyList[ExcitID], Ma, P0, TMath::Abs(HRSAngle))/
					eventdata.Data[kCentralp]-1;
		}

		DEBUG_MASSINFO("PrepareDp","ScatterAngle=%f,\tDpKinOffsets=%f",
					   ScatteringAngle/TMath::Pi()*180,DpKinOffsets);

		// data self consistency check
		if (idx>0)
		{
			const EventData &lasteventdata = fRawData[idx-1];
			if (eventdata.Data[kRunNum] == lasteventdata.Data[kRunNum])
			{
				assert(eventdata.Data[kKineID] == lasteventdata.Data[kKineID]);//check data continuity; check cut definition consistency
				assert(TMath::Abs(eventdata.Data[kCentralp]-lasteventdata.Data[kCentralp])<1e-5);//check data continuity; check cut definition consistency
				assert(TMath::Abs(eventdata.Data[kRealDpKin]-lasteventdata.Data[kRealDpKin])<4e-3);//check data continuity; check cut definition consistency
			}
			else
			{//new run
				DEBUG_INFO("PrepareDp","Run %4.0f : Kinematics #%1.0f, Central p = %fGeV, Excit. State Selected=%f MeV, Dp Kin=%f%%"
						,eventdata.Data[kRunNum]
								,eventdata.Data[kKineID]
										,eventdata.Data[kCentralp]
												,1000*ExcitationEnergy[KineID]
														,100*eventdata.Data[kRealDpKin]);
			}
		}
		else
		{
			//first run
			DEBUG_INFO("PrepareDp","Run %4.0f : Kinematics #%1.0f, Central p = %fGeV, Excit. State Selected=%f MeV, Dp Kin=%f%%"
					,eventdata.Data[kRunNum]
							,eventdata.Data[kKineID]
									,eventdata.Data[kCentralp]
											,1000*ExcitationEnergy[KineID]
													,100*eventdata.Data[kRealDpKin]);
		}
	}

	DEBUG_INFO("PrepareDp","%d out of %d data is for calibration"
			,fNCalibData,fNRawData);
	DEBUG_INFO("PrepareDp","Average : D_Th = %f,\t D_Phi = %f",dth/fNRawData,dphi/fNRawData);
	DEBUG_INFO("PrepareDp","Average : ScatteringAngle = %f",scatang/fNRawData/TMath::Pi()*180);
	DEBUG_INFO("PrepareDp","Average DpKinOffsets = %f, RMS DpKinOffsets = %f"
			,dpkinoff/fNRawData,TMath::Sqrt(dpkinoff_rms/fNRawData));
	DEBUG_INFO("PrepareDp","Average Extended Target Corretion: dp = %f,\t rms_dp = %f"
			,exttargcorr_dp/fNRawData,TMath::Sqrt(rms_exttargcorr_dp/fNRawData));

	//make sure kCalcTh, kCalcPh is filled, although not necessary
	SumSquareDTh();
	SumSquareDPhi();
	SumSquareDp();
}
//_____________________________________________________________________________



// const TVector3 LOpticsOpt::GetSieveHoleTCS(UInt_t Col, UInt_t Row) /*const*/
// {

//   // Double_t new_Z = ZPos;
//   // Double_t X_SH = SieveXbyRow[Row];
//   // cout << "Col = " << Col << endl;

//   // //  Double_t Y_SH = SieveYbyCol[Col];
//   // cout << "X_SH = " << X_SH << " and Y_SH = " << 999 << " and new_Z = " << ZPos << endl;


//   // Double_t Y_SH = SieveYbyCol[Col];
//   // X_SH = SieveXbyRow[Row];
//   // Y_SH = SieveYbyCol[Col];

//   // cout << "X_SH = " << X_SH << " and Y_SH = " << Y_SH << endl;

//   TVector3 SieveHoleTCS(SieveXbyRow[Row]+SieveOffX, SieveYbyCol[Col]+SieveOffY, ZPos);

//   return SieveHoleTCS;



// }

// //_____________________________________________________________________________


// const TVector3 LOpticsOpt::GetSieveHoleTCS(UInt_t Hole) /*const*/
// {

//   std::vector<int> x_y = {};
//   x_y = Get_Col_Row( Hole);

//   UInt_t Col = x_y[0];
//   UInt_t Row = x_y[1];

  
//   TVector3 SieveHoleTCS = {};

//   SieveHoleTCS =  GetSieveHoleTCS(Col, Row);

//   return SieveHoleTCS;



// }

// //_____________________________________________________________________________


// std::vector<int> LOpticsOpt::Get_Col_Row(Int_t Hole){


//   Int_t row_comp = 0;
//   Int_t no_col = 0;
//   Int_t col = 0;
//   Int_t row = 0;


//   for(Int_t i = 0; i<NSieveRow; i++){

//     row_comp += NoinEachRow[i];

//     if( (row_comp-1) >= Hole){
//       row = i;
//       no_col = Hole - ( row_comp - NoinEachRow[i]);
      
//       break;
//     }

//   }



//   if(row%2 == 0){
//     if(no_col==13){
//       col = 25;
//     }
//     else if (no_col==14){
//       col = 26;
//     }
//     else{

//     col = no_col *2;
//     }
//   }

//   if(row%2 == 1){

//     if(row > 1 && row < 15){
//       if(no_col >5){
// 	col = (no_col*2)+3;
//       }
//       else{
// 	col = (no_col*2) +1;
//       }
//     }
//     else{
//       col = (no_col*2) +1;
//     }

    
//   }


//   cout << "For Hole " << Hole << ": row = " << row << " and col = " << col << endl;
//     //  hole_no += Col;
//   std::vector<int> rowcol{col, row};

//   return rowcol;
// }




// //____________________________________________________________________________



// Int_t LOpticsOpt::Get_Hole(Int_t Col, Int_t Row){


//   Int_t hole_no = 0;

//   for(Int_t i = 0; i<Row; i++){

//     hole_no += NoinEachRow[i];

//   }

//   // else if conditions here deal with columns at right edge of sieve slit (area odd rows do not have holes)

//   if(Row%2 == 0){
//     if(Col < 25){
//       hole_no += (Col/2);
//     }
//     else if(Col == 25){
//       hole_no += 13;      
//     }
//     else if(Col == 26){
//       hole_no += 14;
//     }      
//   }


//   if(Row%2 == 1){
//     if(Row > 1 && Row < 15){
//       if(Col < 13){
// 	Col = ((Col+1)/2) - 1;
//       }
//       else if (Col >= 13){
// 	Col = ((Col-3)/2);
//       }      
//     }
//     else{
//       Col = ((Col+1)/2) -1;
//     }    
//     hole_no += Col;
//   }
  

//   //  hole_no += Col;



//   return hole_no;
// }





// //____________________________________________________________________________



void LOpticsOpt::CheckSieve(UInt_t SpecialNLoad, Int_t print){


	Bool_t special = kFALSE;

	if (SpecialNLoad<NFiols)
	{
	  //		DEBUG_INFO("CheckSieve","Special Plot only Setting 0-%d ",SpecialNLoad-1);
		special = kTRUE;
	}

	const UInt_t nplot = SpecialNLoad<NFiols?SpecialNLoad:NFiols;

	//Visualize Sieve Plane

	//	DEBUG_INFO("CheckSieve","Entry Point");

	TH2D * HSievePlane[NFiols] = {0};
	Double_t x_lim[NFiols] = {0};
	Double_t y_lim[NFiols] = {0};

	//  theta vs phi plot addition 

	TH2D* HSieveAngle[NFiols] = {0};
	//	= new TH2F("h2", "theta_target vs. phi_target", 600, -0.04, 0.01, 600, 0.0, 0.08);

	for(UInt_t idx = 0; idx<NFiols; idx++)
	{
		x_lim[idx] = 1.3*TMath::Max(TMath::Abs(SieveYbyCol[0]),TMath::Abs(SieveYbyCol[NSieveCol-1]));
		y_lim[idx] = 1.3*TMath::Max(TMath::Abs(SieveXbyRow[0]),TMath::Abs(SieveXbyRow[NSieveRow-1]));

		if (!special)
		  {
		    HSievePlane[idx] = new TH2D(Form("Sieve_Foil%d",idx),Form("Sieve Plane Proj. (tg_X vs tg_Y) for Data set #%d",idx), 500,-x_lim[idx],x_lim[idx], 500,-y_lim[idx],y_lim[idx]);

		    HSieveAngle[idx] = new TH2D(Form("Angle_Sieve_Foil%d",idx),Form("theta_target vs. phi_target for Data set #%d",idx), 300, -0.04, 0.04, 300, -0.08, 0.08);
		}
		else
		  {
		    
		    HSievePlane[idx] = new TH2D(Form("Sieve_Foil%d",idx),Form("Sieve Plane Proj. (tg_X vs tg_Y) for Data set #%d",idx),200,-x_lim[idx],x_lim[idx],200,-y_lim[idx],y_lim[idx]);

		    HSieveAngle[idx] = new TH2D(Form("Angle_Sieve_Foil%d",idx),Form("theta_target vs. phi_target for Data set #%d",idx), 600, -0.04, 0.04, 600, -0.08, 0.08);

		}

		HSievePlane[idx]->SetXTitle("Sieve H [m]");
		HSievePlane[idx]->SetYTitle("Sieve V [m]");

		HSieveAngle[idx]->SetXTitle("Phi");
		HSieveAngle[idx]->SetYTitle("Theta");

		assert(HSievePlane[idx]);//assure memory allocation
	}



	Double_t dX = 0, dY=0;

	enum {kEventID, kRealSieveX, kRealSieveY, kCalcSieveX, kCalcSieveY};
	//	Double_t SieveEventID[NFiols][NSieveCol][NSieveRow][5]={{{{0}}}};
	Double_t SieveEventID[NFiols][NHoles][5]={{{0}}}; 




	for (UInt_t idx = 0; idx<fNRawData; idx++)
	  {
	    EventData &eventdata = fRawData[idx];
	    
	    // UInt_t res = (UInt_t)eventdata.Data[kCutID];
	    // const UInt_t FoilID = res/(NSieveRow*NSieveCol); //starting 0!
	    // res = res%(NSieveRow*NSieveCol);
	    // const UInt_t Col = res/(NSieveRow); //starting 0!
	    // const UInt_t Row = res%(NSieveRow); //starting 0!

	    
	    const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID]; //starting 0!
	    const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
	    const UInt_t Row = (UInt_t)eventdata.Data[kRowID]; //starting 0!
	  
	    const UInt_t HoleID = Get_Hole(Col,Row);
	    
	    assert(FoilID<NFiols);//array index check
	    

	    const TVector3 SieveHoleTCS = GetSieveHoleTCS(Col,Row);


	    TVector3 BeamSpotHCS(eventdata.Data[kBeamX],eventdata.Data[kBeamY],targetfoils[FoilID]);

	    // TRotation fTCSInHCS;
	    // TVector3 TCSX(0,-1,0);
	    // TVector3 TCSZ(TMath::Sin(HRSAngle),0,TMath::Cos(HRSAngle));
	    // TVector3 TCSY = TCSZ.Cross(TCSX);
	    // fTCSInHCS.RotateAxes(TCSX,TCSY,TCSZ);
	    TVector3 BeamSpotTCS=fTCSInHCS.Inverse()*(BeamSpotHCS-fPointingOffset);
	    
	    TVector3 MomDirectionTCS = SieveHoleTCS - BeamSpotTCS;
	    
	    eventdata.Data[kRealTh] = MomDirectionTCS.X()/MomDirectionTCS.Z();
	    eventdata.Data[kRealPhi] = MomDirectionTCS.Y()/MomDirectionTCS.Z();
	    


	    
	    const Double_t x_tg = BeamSpotTCS.X() - BeamSpotTCS.Z() * tan(eventdata.Data[kCalcTh]);
	    const Double_t y_tg = BeamSpotTCS.Y() - BeamSpotTCS.Z() * tan(eventdata.Data[kCalcPh]);


	    
	    if (special)
	      {
		// fill the Tg_x
		cout << "  SPECIAAAAAL" <<endl;
		assert(eventdata.Data[kRealTgX] == 0);
		
		eventdata.Data[kRealTgX] = BeamSpotTCS.X();
	      }
	    

	    
	    // JW: temp

	    //	    eventdata.Data[kRealTgX] = 0;

	    //normal filling
	    //	    Double_t ProjectionX = BeamSpotTCS.X() + (eventdata.Data[kCalcTh]+eventdata.Data[kRealTgX]*ExtTarCor_ThetaCorr) * (SieveHoleTCS.Z() - BeamSpotTCS.Z());

	    // newer version from Sean's optics

	    Double_t ProjectionX =  x_tg + (TMath::Tan(eventdata.Data[kCalcTh]) + x_tg*ExtTarCor_ThetaCorr) * (ZPos);


	    // JW: temp
	    //	    Double_t ProjectionX = BeamSpotTCS.X() + (eventdata.Data[kCalcTh]) * (SieveHoleTCS.Z() - BeamSpotTCS.Z());

	    //	    Double_t ProjectionY = BeamSpotTCS.Y() + eventdata.Data[kCalcPh] * (SieveHoleTCS.Z() - BeamSpotTCS.Z());

	    Double_t ProjectionY = y_tg + TMath::Tan(eventdata.Data[kCalcPh]) * (ZPos);
	    
	    //		//special check filling
	    //		const Double_t pi = TMath::Pi();
	    ////		ExTgtCor_L_urbb.th*(0.8054) - urbb.y
	    //		Double_t ProjectionX = (eventdata.Data[kCalcTh]+eventdata.Data[kRealTgX]*ExtTarCor_ThetaCorr)*(0.8054) - eventdata.Data[kBeamY];
	    ////		ExTgtCor_l_urbb.ph*(0.8054) + urbb.x*cos(pi/36) - (-0.00444)*sin(pi/36)
	    //		Double_t ProjectionY = eventdata.Data[kCalcPh] *(0.8054) + eventdata.Data[kBeamX] - (-0.00444)*sin(pi/36);


	    // JW: new angular picture of sieve slit
	    HSieveAngle[FoilID]->Fill(eventdata.Data[kCalcPh],eventdata.Data[kCalcTh]);
	    //	    HSieveAngle[FoilID]->Fill(eventdata.Data[kCalcTh],eventdata.Data[kCalcPh]);

	    HSievePlane[FoilID]->Fill(ProjectionY,ProjectionX);

	    
	    
	    dX+= ProjectionX - SieveHoleTCS.X();
	    dY+= ProjectionY - SieveHoleTCS.Y();


	    SieveEventID[FoilID][HoleID][kEventID] = idx;
	    SieveEventID[FoilID][HoleID][kRealSieveX] = SieveHoleTCS.X();
	    SieveEventID[FoilID][HoleID][kRealSieveY] = SieveHoleTCS.Y();
	    SieveEventID[FoilID][HoleID][kCalcSieveX] = ProjectionX;
	    SieveEventID[FoilID][HoleID][kCalcSieveY] = ProjectionY;
	}
	



	//	DEBUG_INFO("CheckSieve","Average : D_X = %f,\t D_Y = %f",dX/fNRawData,dY/fNRawData);


	
	// stat. for each hole
	// TH1D * hDeltaTheta[kMaxDataGroup] = {NULL};
	// TH1D * hDeltaPhi[kMaxDataGroup] = {NULL};
	// TF1 * fDeltaTheta[kMaxDataGroup] = {NULL};
	// TF1 * fDeltaPhi[kMaxDataGroup] = {NULL};

	TH1D * hDeltaTheta[NFiols*NHoles] = {NULL};
	TH1D * hDeltaPhi[NFiols*NHoles] = {NULL};
	TF1 * fDeltaTheta[NFiols*NHoles] = {NULL};
	TF1 * fDeltaPhi[NFiols*NHoles] = {NULL};
	

	const Double_t PhiPrintLimit = 3e-4; //rad
	
	if (!special)
	{
		for (UInt_t i=0; i <NFiols; i++ )
		{
		  for(UInt_t j=0; j <NHoles; j++ )
		    {
		      UInt_t FoilID = i;
		      
		      
		      // const UInt_t Col = res/(NSieveRow); //starting 0!
		      // const UInt_t Row = res%(NSieveRow); //starting 0!
		      
		      std::vector<int> x_y = {};
		      x_y = Get_Col_Row(j);
		      
		      UInt_t Col = x_y[0];
		      UInt_t Row = x_y[1];

		      TString idstring = Form("_%d_%d_%d",FoilID,Col,Row);

		      Int_t id = (FoilID*NHoles) + j;
		      
		      hDeltaTheta[id] = new TH1D("hDeltaTheta"+idstring,"hDeltaTheta"+idstring,200,-.02,.02);
		      hDeltaPhi[id] = new TH1D("hDeltaPhi"+idstring,"hDeltaPhi"+idstring,200,-.01,.01);
		    }
		}

		for (UInt_t idx = 0; idx<fNRawData; idx++)
		{
			EventData &eventdata = fRawData[idx];

			UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID];
			UInt_t ColID = (UInt_t)eventdata.Data[kColID];
			UInt_t RowID = (UInt_t)eventdata.Data[kRowID];

			UInt_t HoleID =  Get_Hole( ColID, RowID);
			
			
			UInt_t id = (FoilID*NHoles) + HoleID;

			assert(hDeltaTheta[id]);
			assert(hDeltaPhi[id]);

			hDeltaTheta[id] -> Fill(eventdata.Data[kCalcTh] - eventdata.Data[kRealTh]);
			hDeltaPhi[id] -> Fill(eventdata.Data[kCalcPh] - eventdata.Data[kRealPhi]);
		}


		for (UInt_t i=0; i <NFiols; i++ )
		  {
		    for(UInt_t j=0; j <NHoles; j++ )
		      {
		// for (UInt_t i=0; i <=NSieveRow*NSieveCol*NFiols; i++ )
		//   {


			const UInt_t FoilID = i; //starting 0!

			std::vector<int> x_y = {};
			x_y = Get_Col_Row(j);
			
			UInt_t Col = x_y[0];
			UInt_t Row = x_y[1];
			
			Int_t id = (FoilID*NHoles) + j;


			if (hDeltaTheta[id]->GetSum()>10)
			  {
			    TString idstring = Form("_%d_%d_%d",FoilID,Col,Row);

			    fDeltaTheta[id] = new TF1("gausth"+idstring,"gaus",-.02,.02);
			    fDeltaPhi[id] = new TF1("gausph"+idstring,"gaus",-.01,.01);
			    
			    fDeltaTheta[id]->SetParameter(1,hDeltaTheta[id]->GetMean());
			    fDeltaTheta[id]->SetParameter(2,hDeltaTheta[id]->GetRMS());
			    
			    fDeltaPhi[id]->SetParameter(1,hDeltaPhi[id]->GetMean());
			    fDeltaPhi[id]->SetParameter(2,hDeltaPhi[id]->GetRMS());
			    
			    hDeltaTheta[id] -> Fit(fDeltaTheta[id], "RQ0+");
			    hDeltaPhi[id] -> Fit(fDeltaPhi[id], "RQ0+");
			    
			    
			    if (TMath::Abs(fDeltaPhi[id]->GetParameter(1))-fDeltaPhi[id]->GetParError(1) > PhiPrintLimit)
			      {
				//					cout<<FoilID<<"\t"<<Col<<"\t"<<Row<<"\t"
				//							<<hDeltaTheta[i]->GetSum()<<"\t"
				//							<<fDeltaPhi[i]->GetParameter(1) <<" +/- "<<fDeltaPhi[i]->GetParError(1)<<"\t"
				//							<<fDeltaPhi[i]->GetParameter(2) <<"\n";
			      }
			    
			  }
			
			
		      } // here
		  }
		
	} 



	TCanvas * c1 = NULL; //= new TCanvas("SieveCheck","SieveCheck",1800,1100);
	if (nplot<=1)
	{
		c1 = new TCanvas("SieveCheck","SieveCheck",800,1100);
		c1->Divide(1,1);
	}
	else if (nplot<=3)
	{
		c1 = new TCanvas("SieveCheck","SieveCheck",1800,1100);
		c1->Divide(3,1);
	}
	else if (nplot<=6)
	{
		c1 = new TCanvas("SieveCheck","SieveCheck",1800,1100);
		c1->Divide(3,2);
	}
	else if (nplot<=8)
	{
		c1 = new TCanvas("SieveCheck","SieveCheck",1800,1100);
		c1->Divide(4,2);
	}
	else if (nplot<=12)
	{
		c1 = new TCanvas("SieveCheck","SieveCheck",1800,1100);
		c1->Divide(4,3);
	}

	assert (c1);

	//JW: additional canvas to be used for angular sieve slit picture

	TCanvas * c2 = NULL; //= new TCanvas("SieveCheck_angles","SieveCheck_angles",1800,1100);
	if (nplot<=1)
	{
		c2 = new TCanvas("SieveCheck_angles","SieveCheck_angles",800,1100);
		c2->Divide(1,1);
	}
	else if (nplot<=3)
	{
		c2 = new TCanvas("SieveCheck_angles","SieveCheck_angles",1800,1100);
		c2->Divide(3,1);
	}
	else if (nplot<=6)
	{
		c2 = new TCanvas("SieveCheck_angles","SieveCheck_angles",1800,1100);
		c2->Divide(3,2);
	}
	else if (nplot<=8)
	{
		c2 = new TCanvas("SieveCheck_angles","SieveCheck_angles",1800,1100);
		c2->Divide(4,2);
	}
	else if (nplot<=12)
	{
		c2 = new TCanvas("SieveCheck_angles","SieveCheck_angles",1800,1100);
		c2->Divide(4,3);
	}

	assert (c2);


	for(UInt_t idx = 0; idx<nplot; idx++)
	{
		UInt_t FoilID = idx;

		c1->cd(idx+1);
		assert(HSievePlane[idx]);//pointer check


		if (!special)
		{
			HSievePlane[idx]->Draw("COLZ");
		}
		else
		{
			HSievePlane[idx]->Draw("COLZ");
			HSievePlane[idx]->Draw("contzsame");
		}



		// Double_t MaxPlot = .05;

		// for(UInt_t Row = 0; Row < NSieveRow; Row++)
		// {
		// 	Double_t SievePos =  SieveOffX + SieveXbyRow[Row];

		// 	TLine *l = new TLine(-MaxPlot,SievePos,+MaxPlot,SievePos);

		// 	if (Row % 2 ==0 )
		// 		l->SetLineColor(kRed);
		// 	else
		// 		l->SetLineColor(kBlue);

		// 	l->Draw();
		// }

		// MaxPlot = .1;
		// for(UInt_t Col = 0; Col < NSieveCol; Col++)
		// {
		// 	Double_t SievePos =  SieveOffY + SieveYbyCol[Col];
		// 	TLine *l = new TLine(SievePos,-MaxPlot,SievePos,MaxPlot);

		// 	if (Col==0 or Col==2 or Col==5 or Col==8 or Col==11)
		// 		l->SetLineColor(kRed);
		// 	else
		// 		l->SetLineColor(kBlue);

		// 	l->Draw();
		// }

		const Double_t plotwidth = 0.0015;
		for(UInt_t Hole = 0; Hole < NHoles; Hole++){


				// Bool_t redr = Row % 2 ==0;
				// Bool_t redc = Col==0 or Col==2 or Col==5 or Col==8 or Col==11;

				Color_t color = kBlack;
				Int_t width = 1;

				if (Hole == 112){
				  color = kRed;
				  width = 2;
				    }
				if( Hole == 160){
				  color = kRed;
				  width = 2;
				}

				// if (redr and redc) color=kRed;
				// else if (!redr and !redc) color=kBlue;
				// else continue;


				// std::vector<int> x_y = {};
				// x_y = Get_Col_Row(j);
		      
				// UInt_t Col = x_y[0];
				// UInt_t Row = x_y[1];

				

				// const Double_t posx = SieveOffY + SieveYbyCol[Col];
				// const Double_t posy = SieveOffX + SieveXbyRow[Row];

				TVector3 Hole_pos = GetSieveHoleTCS(Hole);

				Double_t posy = Hole_pos.X();
				Double_t posx = Hole_pos.Y();

				//				cout << "For Hole: " << Hole << ": posy = " << posy << " and posx = " << posx << endl;


				// + type crosses
				// TLine *lh = new TLine(posx-plotwidth,posy,posx+plotwidth,posy);
				// TLine *lv = new TLine(posx,posy-plotwidth,posx,posy+plotwidth);
				
				
				// Saltire-type crosses
				TLine *lc1 = new TLine(posx-plotwidth,posy-plotwidth,posx+plotwidth,posy+plotwidth);
				lc1->SetLineWidth(width);
				TLine *lc2 = new TLine(posx-plotwidth,posy+plotwidth,posx+plotwidth,posy-plotwidth);
				lc2->SetLineWidth(width);




				// lh -> SetLineColor(color);
				// lv -> SetLineColor(color);

				// lh -> Draw();
				// lv -> Draw();

				lc1->SetLineColor(color);
				lc2->SetLineColor(color);

				lc1->Draw();
				lc2->Draw();





				// mark bad spots

				int i = (FoilID*NHoles) + Hole;

				if (fDeltaPhi[i]!=NULL)
				{
					if (TMath::Abs(fDeltaPhi[i]->GetParameter(1))-fDeltaPhi[i]->GetParError(1) > PhiPrintLimit)
					{
						TEllipse * e = new TEllipse(
								posx+fDeltaPhi[i]->GetParameter(1)*ZPos,
								posy+fDeltaTheta[i]->GetParameter(1)*ZPos,
								fDeltaPhi[i]->GetParameter(2)*ZPos,
								fDeltaTheta[i]->GetParameter(2)*ZPos
								);

						e->SetLineColor(kMagenta);
						e->SetLineWidth(2);
						e->SetFillStyle(0);
						e->Draw();

						cout<<FoilID<<"\t"<<Hole<<"\t"
						    <<hDeltaTheta[i]->GetSum()<<"\t"
								<<fDeltaPhi[i]->GetParameter(1) <<" +/- "<<fDeltaPhi[i]->GetParError(1)<<"\t"
								<<fDeltaPhi[i]->GetParameter(2) <<"\n";
					}

				}


		}
		
		
	 

		
		if (!special)
		  {
		    //draw arrows
		    for(UInt_t Hole = 0; Hole < NHoles; Hole++){
		      
		      if (SieveEventID[FoilID][Hole][kEventID]>0)
			{
			  
			  
			  std::vector<int> x_y = {};
			  x_y = Get_Col_Row(Hole);
			  
			  UInt_t Col = x_y[0];
			  UInt_t Row = x_y[1];
			  
			  assert(SieveEventID[FoilID][Hole][kEventID] < fNRawData);//array index bondary check
			  TArrow * ar2 = new TArrow(SieveEventID[FoilID][Hole][kCalcSieveY]
						    ,SieveEventID[FoilID][Hole][kCalcSieveX]
						    ,SieveEventID[FoilID][Hole][kRealSieveY]
						    ,SieveEventID[FoilID][Hole][kRealSieveX]
						    ,0.008,"|>");
			  ar2->SetAngle(40);
			  ar2->SetLineColor(kMagenta);
			  ar2->SetFillColor(kMagenta);
			  
			  const Double_t ignorelimit = 0.006;
			  // if ((ar2->GetX1()-ar2->GetX2())*(ar2->GetX1()-ar2->GetX2())
			  //     +(ar2->GetY1()-ar2->GetY2())*(ar2->GetY1()-ar2->GetY2())
			  //     >ignorelimit*ignorelimit ) // JW change sign here
			  if (1==1){
			    ar2->Draw();
			  }
			  
			}
		    }
		  }
	}
	//	return c1;


	for(UInt_t idx = 0; idx<nplot; idx++)
	{
		UInt_t FoilID = idx;

		c2->cd(idx+1);


		if (!special)
		{
			HSieveAngle[idx]->Draw("COLZ");
		}
		else
		{
			HSieveAngle[idx]->Draw("COLZ");
			HSieveAngle[idx]->Draw("contzsame");
		}


		const Double_t plotwidth = 0.00075;
		for(UInt_t Hole = 0; Hole < NHoles; Hole++){


				// Bool_t redr = Row % 2 ==0;
				// Bool_t redc = Col==0 or Col==2 or Col==5 or Col==8 or Col==11;

				Color_t color = kBlack;
				Int_t width = 1;

				if (Hole == 112){
				  color = kRed;
				  width = 2;
				    }
				if( Hole == 160){
				  color = kRed;
				  width = 2;
				}

				// if (redr and redc) color=kRed;
				// else if (!redr and !redc) color=kBlue;
				// else continue;


				// std::vector<int> x_y = {};
				// x_y = Get_Col_Row(j);
		      
				// UInt_t Col = x_y[0];
				// UInt_t Row = x_y[1];

				

				// const Double_t posx = SieveOffY + SieveYbyCol[Col];
				// const Double_t posy = SieveOffX + SieveXbyRow[Row];

				TVector3 Hole_pos = GetSieveHoleTCS(Hole);


				Double_t BeamX_average = 0;
				Double_t BeamY_average = 0;
				
				if (FoilID == 0){
				  BeamX_average = BeamX_average_V1;
				  BeamY_average = BeamY_average_V1;
				}
				if (FoilID == 1){
				  BeamX_average = BeamX_average_V2;
				  BeamY_average = BeamY_average_V2;
				}
				if (FoilID == 2){
				  BeamX_average = BeamX_average_V3;
				  BeamY_average = BeamY_average_V3;
				}

				
				// BeamX_average = -0.0004606;
				// BeamY_average = 0.002448;

				    
				BeamX_average = -0.0006361;
				BeamY_average = 0.002419;    


				TVector3 BeamSpotHCS_average(BeamX_average,BeamY_average,targetfoils[FoilID]);


				TVector3 BeamSpotTCS_average = fTCSInHCS.Inverse()*(BeamSpotHCS_average-fPointingOffset);

				TVector3 MomDirectionTCS_hole = Hole_pos - BeamSpotTCS_average;

				Double_t theta_hole = MomDirectionTCS_hole.X()/MomDirectionTCS_hole.Z();
				Double_t phi_hole = MomDirectionTCS_hole.Y()/MomDirectionTCS_hole.Z();



				// + type crosses
				// TLine *lh = new TLine(posx-plotwidth,posy,posx+plotwidth,posy);
				// TLine *lv = new TLine(posx,posy-plotwidth,posx,posy+plotwidth);
				
				
				// Saltire-type crosses
				TLine *lc1 = new TLine(phi_hole-plotwidth,theta_hole-plotwidth,phi_hole+plotwidth,theta_hole+plotwidth);
				lc1->SetLineWidth(width);
				TLine *lc2 = new TLine(phi_hole-plotwidth,theta_hole+plotwidth,phi_hole+plotwidth,theta_hole-plotwidth);
				lc2->SetLineWidth(width);




				// lh -> SetLineColor(color);
				// lv -> SetLineColor(color);

				// lh -> Draw();
				// lv -> Draw();

				lc1 -> SetLineColor(color);
				lc2-> SetLineColor(color);

				lc1 -> Draw();
				lc2 -> Draw();
	        



		}


	}

	cout << "Overall number of entries = " << fNRawData <<endl;

	cout << "Number of entries in Sieve plot (x vs y) = " << HSievePlane[0]->GetEntries() << endl;

	cout << "Number of entries in Sieve plot (theta vs phi) = " << HSieveAngle[0]->GetEntries() << endl;


	//	return c1;


	if(print != 0){

	  c1->Print(fDestDataBaseName + ".Sieve_angles.Opt.png", "png");
	  c2->Print(fDestDataBaseName + ".Sieve_XY.Opt.png", "png");


	}

		
}






void  LOpticsOpt::CheckVertex(Int_t print)
{
	//Visualize ReactZ spectrum

	DEBUG_INFO("CheckVertex","Entry Point");

	const Double_t VertexRange=.5;
	TH2D * hYVSReactZ = new TH2D("hYVSReactZ","Rot. Y VS ReactZ",1000,-VertexRange,VertexRange, 1000,-.1,.1);
	TH2D * hPhVSReactZ = new TH2D("hPhVSReactZ","Rot. Ph VS ReactZ",1000,-VertexRange,VertexRange, 1000,-.05,.05);
	TH2D * hrbxVSReactZ = new TH2D("hrbxVSReactZ","Beam X VS ReactZ",1000,-VertexRange,VertexRange, 1000,-.01,.01);
	TH1D * hReactZ = new TH1D("hReactZ","2 LHRS ReactZ",400,-VertexRange,VertexRange);

	Double_t dtg_vz = 0, dtg_vz_rms = 0, reactz=0;
	Double_t dtg_vz_foil[NFiols] ={0};
	Double_t dtg_y_foil[NFiols] ={0};
	Double_t rmstg_vz_foil[NFiols] ={0};
	Double_t rmstg_y_foil[NFiols] ={0};
	UInt_t ndata_foil[NFiols] ={0};




	// add plots for individual holes dtg_vz and dtg_y

	Double_t dtg_vz_rms_hole[NHoles][NFiols] = {0};
	Double_t dtg_vz_hole[NHoles][NFiols] = {0};


	Double_t dtg_y_rms_hole[NHoles][NFiols] = {0};
	Double_t dtg_y_hole[NHoles][NFiols] = {0};

	// keep track of number of events for each hole
	Int_t hole_events[NHoles][NFiols] = {0};


	// 3d plots of difference between 'real' and epxted target y and reactz

	TH3F* h_z_diff[NFiols];

	TH3F* h_ytg_diff[NFiols];

	
	// 2D plots for averages of tgy and reectz rms
	TH2D* h_dtgy_rms[NFiols];

	TH2D* h_dvz_rms[NFiols];


	for(UInt_t FoilID = 0; FoilID<NFiols; FoilID++){
	  
	  h_dtgy_rms[FoilID]  = new TH2D(Form("dTgY av rms all holes Foil_%d",FoilID),Form("dTgY av rms all holes Foil_%d",FoilID),27,0,26,17,0,16);

	  h_dtgy_rms[FoilID]->GetXaxis()->SetTitle("Column #");
	  h_dtgy_rms[FoilID]->GetYaxis()->SetTitle("Row #");


	  h_dvz_rms[FoilID] = new TH2D(Form("dVz rms all holes Foil_%d",FoilID),Form("dVz av rms all holes Foil_%d",FoilID),27,0,26,17,0,16);
	  h_dvz_rms[FoilID]->GetXaxis()->SetTitle("Column #");
	  h_dvz_rms[FoilID]->GetYaxis()->SetTitle("Row #");


	  h_z_diff[FoilID]  = new TH3F(Form("dVZ for all holes Foil_%d",FoilID),Form("dVZ for all holes Foil_%d",FoilID),27,0,26,17,0,16,100,-0.1,0.1);


	  //	  h_ytg_diff[NFiols];
	}


	//	TH2D* h_dtgy = new TH2D("dTgY all holes","dTgY all holes",27,0,26,17,0,16);

	//	TH2D* h_dvz = new TH2D("dVz all holes","dVz all holes",27,0,26,17,0,16);





	// add plots for tg_y optimisation

	TH2D * hYVSrealTgY[NFiols];
	TH2D * hYVScalcTgY[NFiols];


	for(UInt_t FoilID = 0; FoilID<NFiols; FoilID++){
	  hYVSrealTgY[FoilID] = new TH2D(Form("hYVSrealTgY_%d",FoilID),Form("Rot. Y VS 'real' tg_y (Foil_%d)",FoilID),1000,-0.1,0.1, 1000,-0.04,0.04);

	  hYVScalcTgY[FoilID] = new TH2D(Form("hYVScalcTgY_%d",FoilID),Form("Rot. Y VS 'calc' tg_y (Foil_%d)",FoilID),1000,-0.1,0.1, 1000,-0.04,0.04);
	
	}


	// tg_y spectrum (calc)
	TH1D * hTgY_calc = new TH1D("hTgY_calc","LHRS target_y",400,-0.04,0.04);
	TH1D * hTgY_real = new TH1D("hTgY_real","LHRS target_y",400,-0.04,0.04);

	TH1D * hTgY_diff = new TH1D("hTgY_diff","LHRS tg_y: calc - real",400,-0.01,0.01);


	// tf1s to fit target y difference

	TF1*  tgy_fits;
	TF1* tgy_fits_2;
	

	
	//  plots for real and calc Z-spectrum
	

	TH1D * hVz_real = new TH1D("hVZ_real","2 LHRS ReactZ",400,-VertexRange,VertexRange);

	


	//calculate kCalcReactZ
	for (UInt_t idx = 0; idx<fNRawData; idx++)
	{
	  EventData &eventdata = fRawData[idx];
	  
	  // 		TVector3 BeamSpotHCS(
	  // 							 eventdata.Data[kBeamX],
	  // 		eventdata.Data[kBeamY],
	  //   eventdata.Data[kRealReactZ]
	  // 							);
	  
	  TVector3 Tg_YSpotTCS(0,eventdata.Data[kCalcTgY],0);
	  //	  TVector3 MomDirectionTCS(0,eventdata.Data[kL_tr_tg_ph],1);

	  const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID];
	  

	  TVector3 BeamSpotHCS(eventdata.Data[kBeamX],eventdata.Data[kBeamY],targetfoils[FoilID]);
	  
	  TVector3 BeamSpotTCS=fTCSInHCS.Inverse()*(BeamSpotHCS-fPointingOffset);
	  
	  

	  const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
	  const UInt_t Row = (UInt_t)eventdata.Data[kRowID];

	  const UInt_t Hole = Get_Hole(Col,Row);

	  TVector3 SieveHoleTCS = GetSieveHoleTCS(Col,Row);
      
	  const TVector3 MomDirectionTCS = SieveHoleTCS - BeamSpotTCS;

	  
	  
	  //TVector3 Tg_YSpotHCS=fTCSInHCS*Tg_YSpotTCS+fPointingOffset;
	  TVector3 MomDirectionHCS=fTCSInHCS*MomDirectionTCS;
	  
	  //	  assert(Tg_YSpotHCS.Y()==MissPointY);//internal consistency check
	  //	  assert(MomDirectionHCS.Y()==0);//internal consistency check
	  
	  //		reactz = Tg_YSpotHCS.Z()
	  //				- (Tg_YSpotHCS.X()-eventdata.Data[kBeamX])
	  //				/MomDirectionHCS.X()*MomDirectionHCS.Z();
	  

	  
	  //	  reactz = (Tg_YSpotHCS.X() - eventdata.Data[kBeamX] - (MomDirectionHCS.X()/MomDirectionHCS.Z())*Tg_YSpotHCS.Z() )/( eventdata.Data[kBeamDirX]/eventdata.Data[kBeamDirZ] - (MomDirectionHCS.X()/MomDirectionHCS.Z()) );


	  // second reactz definition

	  const Int_t a = (HRSAngle > 0) ? 1 : -1;
      
	  Double_t Real_Tg_Phi = MomDirectionTCS.Y() / MomDirectionTCS.Z();
      
	  reactz = - ( eventdata.Data[kCalcTgY] -a*MissPointZ)*TMath::Cos(TMath::ATan(Real_Tg_Phi))/TMath::Sin(HRSAngle + TMath::ATan(Real_Tg_Phi)) + BeamSpotHCS.X()*TMath::Cos(HRSAngle+TMath::ATan(Real_Tg_Phi))/TMath::Sin(HRSAngle+TMath::ATan(Real_Tg_Phi));

	  
	  //	for (UInt_t idx = 0; idx<fNRawData; idx++)
	  // JW: Checking react-z calcultaion

	  



	  eventdata.Data[kCalcReactZ] = reactz;
	  
	  dtg_vz += reactz - eventdata.Data[kRealReactZ];
	  dtg_vz_rms += (reactz - eventdata.Data[kRealReactZ])*(reactz - eventdata.Data[kRealReactZ]);


	  h_z_diff[FoilID]->Fill(Col,Row,reactz - eventdata.Data[kRealReactZ]);


	  // fill arrays which give rms values for each hole (in zreact and tg_y)


	  dtg_vz_hole[Hole][FoilID] += reactz - eventdata.Data[kRealReactZ];
	  dtg_vz_rms_hole[Hole][FoilID] += (reactz - eventdata.Data[kRealReactZ])*(reactz - eventdata.Data[kRealReactZ]);


	  dtg_y_hole[Hole][FoilID] += eventdata.Data[kCalcTgY] - eventdata.Data[kRealTgY];
	  dtg_y_rms_hole[Hole][FoilID] += (eventdata.Data[kCalcTgY] - eventdata.Data[kRealTgY])*(eventdata.Data[kCalcTgY] - eventdata.Data[kRealTgY]);


	  hole_events[Hole][FoilID]++;

	  
	  
	  
	  hYVSReactZ->Fill(reactz,eventdata.Data[kY]);
	  hPhVSReactZ->Fill(reactz,eventdata.Data[kPhi]);
	  hrbxVSReactZ->Fill(reactz,eventdata.Data[kBeamX]);
	  hReactZ->Fill(reactz);

	  hVz_real->Fill(targetfoils[FoilID]);
	  
	  //statics by foils
	  //	  const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID];
	  ndata_foil[FoilID]++;
	  
	  dtg_vz_foil[FoilID] += eventdata.Data[kCalcReactZ] - eventdata.Data[kRealReactZ];
	  rmstg_vz_foil[FoilID] += (eventdata.Data[kCalcReactZ] - eventdata.Data[kRealReactZ])*(eventdata.Data[kCalcReactZ] - eventdata.Data[kRealReactZ]);
	  
	  dtg_y_foil[FoilID] += eventdata.Data[kCalcTgY] - eventdata.Data[kRealTgY];
	  rmstg_y_foil[FoilID] += (eventdata.Data[kCalcTgY] - eventdata.Data[kRealTgY])*(eventdata.Data[kCalcTgY] - eventdata.Data[kRealTgY]);
	 
	  // fill tg_y plots

	  hYVSrealTgY[FoilID]->Fill(eventdata.Data[kRealTgY],eventdata.Data[kY]);
	  hYVScalcTgY[FoilID]->Fill(eventdata.Data[kCalcTgY],eventdata.Data[kY]);

	  hTgY_calc->Fill(eventdata.Data[kCalcTgY]);
	  hTgY_real->Fill(eventdata.Data[kRealTgY]);	 	  	  
	  hTgY_diff->Fill(eventdata.Data[kCalcTgY] -eventdata.Data[kRealTgY]);
	}
	
	DEBUG_INFO("CheckVertex","dtg_vz = %f,\t dtg_vz_rms = %f",dtg_vz/fNRawData,dtg_vz_rms/fNRawData);
	
	TCanvas * c_vertcheck = new TCanvas("CheckVertex","SieveCheck",1800,900);
	c_vertcheck->Divide(3,1);
	UInt_t idx=1;
	
	c_vertcheck->cd(idx++);
	hYVSReactZ->Draw("COLZ");
	c_vertcheck->cd(idx++);
	hPhVSReactZ->Draw("COLZ");
	c_vertcheck->cd(idx++);
	hrbxVSReactZ->Draw("COLZ");
	c_vertcheck->cd(idx++);
	hReactZ->Draw();
	
	for(idx = 1; idx<=4; idx++)
	  {
	    c_vertcheck->cd(idx);
	    
	    for(UInt_t FoilID = 0; FoilID<NFiols; FoilID++)
	      {
		const Double_t MaxPlot = 2000;
		
		TLine *l = new TLine(targetfoils[FoilID],-MaxPlot,targetfoils[FoilID],+MaxPlot);
		l->SetLineColor(6);
		l->Draw();
	      }
	  }
	
	cout<<"Statistic on each foils:\n";
	cout<<"Foil ID\td_vz\trms_vz\td_tgy\trms_tgy\tNev\n";
	for(UInt_t FoilID = 0; FoilID<NFiols; FoilID++)
	  {
	    // 		assert(ndata_foil[FoilID]);//at least 1 event on each foil
	    UInt_t N = ndata_foil[FoilID];
	    if (N<=0) N=1;
	    
	    printf("%d\t%f\t%f\t%f\t%f\t%d\n",FoilID
		   ,dtg_vz_foil[FoilID]/N
		   ,rmstg_vz_foil[FoilID]/N
		   ,dtg_y_foil[FoilID]/N
		   ,rmstg_y_foil[FoilID]/N
		   ,N);
	  }
	

	// second canvas for tg_y 

	
	TCanvas * c_tg_y = new TCanvas("CheckTgY","SieveCheckRes TgY",1800,900);

	c_tg_y->Divide(3,3);

	c_tg_y->cd(1);
	hYVSrealTgY[0]->Draw("col");
	c_tg_y->cd(2);
	hYVScalcTgY[0]->Draw("col");

	c_tg_y->cd(4);
	hYVSrealTgY[1]->Draw("col");
	c_tg_y->cd(5);
	hYVScalcTgY[1]->Draw("col");


	c_tg_y->cd(7);
	hYVSrealTgY[2]->Draw("col");
	c_tg_y->cd(8);
	hYVScalcTgY[2]->Draw("col");


	c_tg_y->cd(3);

	hTgY_calc->SetLineColor(kRed);
	hTgY_calc->Draw();

	hTgY_real->SetLineColor(kBlue);
	hTgY_real->Draw("same");

	c_tg_y->cd(6);
	hTgY_diff->Draw();

	
	tgy_fits = new TF1("TgY fit 1","gaus",-0.01,0.01);

	hTgY_diff->Fit("TgY fit 1","R0");

	tgy_fits_2 = new TF1("TgY fit 2","gaus",tgy_fits->GetParameter(1)-2.*tgy_fits->GetParameter(2),tgy_fits->GetParameter(1)+2.*tgy_fits->GetParameter(2));
	hTgY_diff->Fit("TgY fit 2","R");
	
	tgy_fits_2->SetLineColor(kRed);
	tgy_fits_2->Draw("same");
	

	c_tg_y->cd(9);

	hReactZ->SetLineColor(kRed);
	hReactZ->Draw();

	hVz_real->SetLineColor(kBlue);
	hVz_real->Draw("same");


	


	TCanvas* c_vert_holes = new TCanvas("Hole check","Hole check",1800,900);

	c_vert_holes->Divide(3,2);


	for(Int_t foil_no = 0; foil_no < NFiols; foil_no++){
	  for(Int_t hole_no = 0; hole_no < NHoles; hole_no++){

	    
	    std::vector<int> hole_cr = Get_Col_Row(hole_no);
	    
	    Int_t hole_col = hole_cr[0];
	    
	    Int_t hole_row = hole_cr[1];
	    
	    
	    if(hole_events[hole_no][foil_no]>0){
	      	      

	      h_dtgy_rms[foil_no]->Fill(hole_col,hole_row,dtg_y_rms_hole[hole_no][foil_no]/hole_events[hole_no][foil_no]);

	      h_dvz_rms[foil_no]->Fill(hole_col,hole_row,dtg_vz_rms_hole[hole_no][foil_no]/hole_events[hole_no][foil_no]);

	      //	      h_z_diff[foil_no]->Fill(hole_col,hole_row,dtg_vz_hole[hole_no][foil_no]);
	      




	    }
	  //	  h_dtgy->Fill(dtg_y_hole
	    

	  }
	}


	c_vert_holes->cd(1);
	h_dtgy_rms[0]->Draw("lego2");

	c_vert_holes->cd(4);
	h_dvz_rms[0]->Draw("lego2");

	c_vert_holes->cd(2);
	h_dtgy_rms[1]->Draw("lego2");

	c_vert_holes->cd(5);
	h_dvz_rms[1]->Draw("lego2");


	c_vert_holes->cd(3);
	h_dtgy_rms[2]->Draw("lego2");

	c_vert_holes->cd(6);
	h_dvz_rms[2]->Draw("lego2");



	TCanvas* c_holes_3D = new TCanvas("3D Hole check","3D Hole check",1800,900);

	c_holes_3D->Divide(3,2);

	c_holes_3D->cd(1);
	h_z_diff[0]->Draw();

	c_holes_3D->cd(2);
	h_z_diff[1]->Draw();

	c_holes_3D->cd(3);
	h_z_diff[2]->Draw();



	// print canvases
	if(print != 0){


	  c_vertcheck->Print(fDestDataBaseName + ".Vertex.Opt.png", "png"); // 

	  c_tg_y->Print(fDestDataBaseName + ".Vertex_TgY_check.Opt.png", "png");   // "check TgY"

	    // c_vert_holes // 2D plots
	    // c_holes_3D // 3D plots
	}





	
}

void LOpticsOpt::CheckVertexRes(Int_t print)
{
	//Visualize ReactZ spectrum


  gStyle->SetOptFit(1111);
  
  DEBUG_INFO("CheckVertex","Entry Point");

	const Double_t VertexRange=.5;
	TH2D * hYVSReactZ = new TH2D("hYVSReactZ","Rot. Y VS ReactZ",1000,-VertexRange,VertexRange, 1000,-.1,.1);
	TH2D * hPhVSReactZ = new TH2D("hPhVSReactZ","Rot. Ph VS ReactZ",1000,-VertexRange,VertexRange, 1000,-.05,.05);
	TH2D * hrbxVSReactZ = new TH2D("hrbxVSReactZ","Beam X VS ReactZ",1000,-VertexRange,VertexRange, 1000,-.01,.01);
	TH1D * hReactZ = new TH1D("hReactZ","LHRS ReactZ",400,-VertexRange,VertexRange);

	//JW: add array of TH1Ds to plot hReactZ for each foil




	TH1D * hReacts_foils[NFiols];

	for(UInt_t FoilID = 0; FoilID<NFiols; FoilID++){
	  hReacts_foils[FoilID] = new TH1D(Form("hReactZ_%d",FoilID),Form("LHRS_ReactZ_%d",FoilID),400,-VertexRange,VertexRange);
	}

	  

	Double_t dtg_vz = 0, dtg_vz_rms = 0, reactz=0;
	Double_t dtg_vz_foil[NFiols] ={0};
	Double_t dtg_y_foil[NFiols] ={0};
	Double_t rmstg_vz_foil[NFiols] ={0};
	Double_t rmstg_y_foil[NFiols] ={0};
	UInt_t ndata_foil[NFiols] ={0};

	//calculate kCalcReactZ
	for (UInt_t idx = 0; idx<fNRawData; idx++)
	{
	  EventData &eventdata = fRawData[idx];
	  
	  const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID];

	  // 		TVector3 BeamSpotHCS(
	  // 							 eventdata.Data[kBeamX],
	  // 		eventdata.Data[kBeamY],
	  //   eventdata.Data[kRealReactZ]
	  // 							);
	  
	  TVector3 Tg_YSpotTCS(0,eventdata.Data[kCalcTgY],0);

	  

	  // original version of MomDirectionTCS
	  //	  TVector3  MomDirectionTCS(0,eventdata.Data[kL_tr_tg_ph],1);

	   // TRotation fTCSInHCS;
	   // TVector3 TCSX(0,-1,0);
	   // TVector3 TCSZ(TMath::Sin(HRSAngle),0,TMath::Cos(HRSAngle));
	   // TVector3 TCSY = TCSZ.Cross(TCSX);
	   // fTCSInHCS.RotateAxes(TCSX,TCSY,TCSZ);
	  


	  TVector3 BeamSpotHCS(eventdata.Data[kBeamX],eventdata.Data[kBeamY],targetfoils[FoilID]);


	  TVector3 BeamSpotTCS=fTCSInHCS.Inverse()*(BeamSpotHCS-fPointingOffset);

	  const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
	  const UInt_t Row = (UInt_t)eventdata.Data[kRowID];
	  TVector3 SieveHoleTCS = GetSieveHoleTCS(Col,Row);
	  
	  const TVector3 MomDirectionTCS = SieveHoleTCS - BeamSpotTCS;

	  
	  TVector3 Tg_YSpotHCS=fTCSInHCS*Tg_YSpotTCS+fPointingOffset;
	  TVector3 MomDirectionHCS=fTCSInHCS*MomDirectionTCS;
	  
	  //	  assert(Tg_YSpotHCS.Y()==MissPointY);//internal consistency check
	  //	  assert(MomDirectionHCS.Y()==0);//internal consistency check
	  
	  //		reactz = Tg_YSpotHCS.Z()
	  //				- (Tg_YSpotHCS.X()-eventdata.Data[kBeamX])
	  //				/MomDirectionHCS.X()*MomDirectionHCS.Z();
	  
	  
	  //	  reactz= (Tg_YSpotHCS.X() - eventdata.Data[kBeamX] - (MomDirectionHCS.X()/MomDirectionHCS.Z())*Tg_YSpotHCS.Z() )/( eventdata.Data[kBeamDirX]/eventdata.Data[kBeamDirZ] - (MomDirectionHCS.X()/MomDirectionHCS.Z()) );

	  // Alternative definition of z-vertex in Sean's optics

	  const Int_t a = (HRSAngle > 0) ? 1 : -1;
	  //	  CalcReacZ = - ( eventdata.Data[kCalcTgY] -a*MissPointZ)*TMath::Cos(TMath::ATan(Real_Tg_Phi))/TMath::Sin(HRSAngle + TMath::ATan(Real_Tg_Phi)) + BeamSpotHCS.X()*TMath::Cos(HRSAngle+TMath::ATan(Real_Tg_Phi))/TMath::Sin(HRSAngle+TMath::ATan(Real_Tg_Phi));

	  Double_t Real_Tg_Phi = MomDirectionTCS.Y() / MomDirectionTCS.Z();


	 
	  
	  reactz = - ( eventdata.Data[kCalcTgY] -a*MissPointZ)*TMath::Cos(TMath::ATan(Real_Tg_Phi))/TMath::Sin(HRSAngle + TMath::ATan(Real_Tg_Phi)) + BeamSpotHCS.X()*TMath::Cos(HRSAngle+TMath::ATan(Real_Tg_Phi))/TMath::Sin(HRSAngle+TMath::ATan(Real_Tg_Phi));





	  
	  //	for (UInt_t idx = 0; idx<fNRawData; idx++)
	  // JW: Checking react-z 	  

	  eventdata.Data[kCalcReactZ] = reactz;
	  
	  dtg_vz += reactz - eventdata.Data[kRealReactZ];
	  dtg_vz_rms += (reactz - eventdata.Data[kRealReactZ])*(reactz - eventdata.Data[kRealReactZ]);
	  
	  hYVSReactZ->Fill(reactz,eventdata.Data[kY]);
	  hPhVSReactZ->Fill(reactz,eventdata.Data[kPhi]);
	  hrbxVSReactZ->Fill(reactz,eventdata.Data[kBeamX]);
	  hReactZ->Fill(reactz);


	  
	  //statics by foils
	  //	  const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID];
	  ndata_foil[FoilID]++;

	  hReacts_foils[FoilID]->Fill(reactz);
	  
	  dtg_vz_foil[FoilID] += eventdata.Data[kCalcReactZ] - eventdata.Data[kRealReactZ];
	  rmstg_vz_foil[FoilID] += (eventdata.Data[kCalcReactZ] - eventdata.Data[kRealReactZ])*(eventdata.Data[kCalcReactZ] - eventdata.Data[kRealReactZ]);
	  
	  dtg_y_foil[FoilID] += eventdata.Data[kCalcTgY] - eventdata.Data[kRealTgY];
	  rmstg_y_foil[FoilID] += (eventdata.Data[kCalcTgY] - eventdata.Data[kRealTgY])*(eventdata.Data[kCalcTgY] - eventdata.Data[kRealTgY]);
	  
	}
	
	DEBUG_INFO("CheckVertex","dtg_vz = %f,\t dtg_vz_rms = %f",dtg_vz/fNRawData,dtg_vz_rms/fNRawData);
	
	TCanvas * c_vertcheck_res = new TCanvas("CheckVertexRes","SieveCheckRes",1800,900);
	c_vertcheck_res->Divide(3,1);





	// UInt_t idx=1;
	
	// c_vertcheck_res->cd(idx++);
	// hYVSReactZ->Draw("COLZ");
	// c_vertcheck_res->cd(idx++);
	// hPhVSReactZ->Draw("COLZ");
	// c_vertcheck_res->cd(idx++);
	// hrbxVSReactZ->Draw("COLZ");
	// c_vertcheck_res->cd(idx++);
	// hReactZ->Draw();

	TF1* foil_fits[NFiols];
	TF1* foil_fits_2[NFiols];
	


	TPaveStats *st[NFiols];

	
	//  JW: very temporary change to looping through 2 foils
	for(UInt_t FoilID = 0; FoilID<NFiols; FoilID++){


	  cout << " For foil " << FoilID << endl;
	  c_vertcheck_res->cd(FoilID+1);

	  //	  hReacts_foils[FoilID]->Fit(
	  foil_fits[FoilID] = new TF1(Form("Foil fit %i",FoilID),"gaus",targetfoils[FoilID]-0.15,targetfoils[FoilID]+0.15);
	  hReacts_foils[FoilID]->Fit(Form("Foil fit %i",FoilID),"R0");

	  foil_fits_2[FoilID] = new TF1(Form("Foil fit 2 %i",FoilID),"gaus",foil_fits[FoilID]->GetParameter(1)-2.*foil_fits[FoilID]->GetParameter(2),foil_fits[FoilID]->GetParameter(1)+2.*foil_fits[FoilID]->GetParameter(2));
	  hReacts_foils[FoilID]->Fit(Form("Foil fit 2 %i",FoilID),"R");


	  

	  

	  std::cout << "Test of fit results: " << foil_fits[FoilID]->IsValid() << std::endl;


	  hReacts_foils[FoilID]->Draw();
	  if(hReacts_foils[FoilID]->GetEntries()!=0){
	    cout << "entered function drawing " << endl;
	    foil_fits_2[FoilID]->SetLineColor(kRed);
	    foil_fits_2[FoilID]->Draw("same");
	  }


	  // adjust stats box to show fit
	  //	  st[FoilID] = *(hReacts_foils[FoilID]->FindObject("stats"));
	  //	  st[FoilID] = *(hReacts_foils[FoilID]->FindObject(Form("Foil fit%f",FoilID)));

	  gPad->Modified(); gPad->Update(); // make sure it’s (re)drawn
	  
	  //	  st[FoilID] = *(static_cast<TPaveStats*>((hReacts_foils[FoilID]->FindObject(Form("Foil fit%f",FoilID)))));
	  st[FoilID] = (static_cast<TPaveStats*>((hReacts_foils[FoilID]->FindObject("stats"))));

	  TIter next(hReacts_foils[FoilID]->GetListOfFunctions());

	  //	  GetListOfFunctions()
	  TObject* func;
	    
	  while(func = next()){

	    cout << "For foil " << FoilID << " function name is " << func->GetName() << endl;
	  }


	  cout << "Exited while loop" << endl;

	  // st[FoilID]->SetY1NDC(0.7);
	  // st[FoilID]->SetY2NDC(0.90);

	  // if(FoilID <2){
	  //   st[FoilID]->SetX1NDC(0.6);
	  //   st[FoilID]->SetX2NDC(0.9);
	  // }
	  // if(FoilID == 2){
	  //   st[FoilID]->SetX1NDC(0.1);
	  //   st[FoilID]->SetX2NDC(0.4);
	  // }


	  // st[FoilID]->SetName(Form("Foil%d stats",FoilID));
	  // st[FoilID]->SetOptFit(1111);
	  // hReacts_foils[FoilID]->SetStats(0);
	  

	  gPad->Modified(); gPad->Update(); // make sure it’s (re)drawn

	  cout << "Reached after modified" << endl;
	  //	  st[FoilID] = dynamic_cast<TPaveStats>(*(hReacts_foils[FoilID]->FindObject("stats")));
	  //static_cast<THaSpectrometer*>
	 
	  

	}


	

// 	    TF1* f2a=new TF1("f2a","gaus",bin_gaus_min,bin_gaus_max);
// 	  f2a->SetLineColorAlpha(3,1.0);
// 	  TF1* f31=new TF1("f31","pol1",bin_lin_min,bin_lin_max);
// 	  f31->SetLineColorAlpha(2,1.0);
	  
// 	  TF1* total=new TF1("total","gaus(0)+pol1(3)",bin_fit_min,bin_fit_max);
	  
// 	  total->SetLineColor(1);
// 	  double par[5];
// //f2a = ht1->GetFunction("gaus");
// 	  ht1->Fit(f2a,"R");
// 	  ht1->Fit(f31,"R+");
// 	  f2a->GetParameters(&par[0]);


	
	// for(idx = 1; idx<=4; idx++)
	//   {
	//     c_vertcheck_res->cd(idx);
	    
	//     for(UInt_t FoilID = 0; FoilID<NFiols; FoilID++)
	//       {
	// 	const Double_t MaxPlot = 2000;
		
	// 	TLine *l = new TLine(targetfoils[FoilID],-MaxPlot,targetfoils[FoilID],+MaxPlot);
	// 	l->SetLineColor(6);
	// 	l->Draw();
	//       }
	//   }
	
	cout<<"Statistic on each foils:\n";
	cout<<"Foil ID\td_vz\trms_vz\td_tgy\trms_tgy\tNev\n";
	for(UInt_t FoilID = 0; FoilID<NFiols; FoilID++)
	  {
	    // 		assert(ndata_foil[FoilID]);//at least 1 event on each foil
	    UInt_t N = ndata_foil[FoilID];
	    if (N<=0) N=1;
	    
	    printf("%d\t%f\t%f\t%f\t%f\t%d\n",FoilID
		   ,dtg_vz_foil[FoilID]/N
		   ,rmstg_vz_foil[FoilID]/N
		   ,dtg_y_foil[FoilID]/N
		   ,rmstg_y_foil[FoilID]/N
		   ,N);
	  }
	
	//	return c_vertcheck_res;

	// print plots
	if(print != 0){

	  c_vertcheck_res->Print(fDestDataBaseName + ".VertexRes.Opt.png", "png");

	}



}


void LOpticsOpt::CheckVertexDiff(Int_t print)
{
	//Visualize CalcZ - ReactZ spectrum


  gStyle->SetOptFit(1111);


  const Double_t VertexRange=.5;

  TH2D * hrbxVSReactZ = new TH2D("hrbxVSReactZ","Beam X VS ReactZ",1000,-VertexRange,VertexRange, 1000,-.007,.003);
  hrbxVSReactZ->GetYaxis()->SetTitle("Beam x [m]");
  hrbxVSReactZ->GetXaxis()->SetTitle("React Z [m]");


  TH2D * hrbyVSReactZ = new TH2D("hrbyVSReactZ","Beam Y VS ReactZ",1000,-VertexRange,VertexRange, 1000,.0015,.005);
  hrbyVSReactZ->GetYaxis()->SetTitle("Beam y [m]");
  hrbyVSReactZ->GetXaxis()->SetTitle("React Z [m]");



  TH1D * hReactZ_diff = new TH1D("hReactZ_diff","Calculated minus real Zreact",400,-0.2,0.2);
  hReactZ_diff->GetXaxis()->SetTitle("#Deltaz-coord [m]");

  
  
  
  Double_t dtg_vz = 0, dtg_vz_rms = 0, reactz=0;
  Double_t dtg_vz_foil[NFiols] ={0};
  Double_t dtg_y_foil[NFiols] ={0};
  Double_t rmstg_vz_foil[NFiols] ={0};
  Double_t rmstg_y_foil[NFiols] ={0};
  UInt_t ndata_foil[NFiols] ={0};
  


  //calculate kCalcReactZ
  for (UInt_t idx = 0; idx<fNRawData; idx++)
    {
      
      EventData &eventdata = fRawData[idx];
	  
      const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID];
      
	  
      TVector3 Tg_YSpotTCS(0,eventdata.Data[kCalcTgY],0);
      
	  

      TVector3 BeamSpotHCS(eventdata.Data[kBeamX],eventdata.Data[kBeamY],targetfoils[FoilID]);
      
      
      TVector3 BeamSpotTCS=fTCSInHCS.Inverse()*(BeamSpotHCS-fPointingOffset);
      
      const UInt_t Col = (UInt_t)eventdata.Data[kColID]; //starting 0!
      const UInt_t Row = (UInt_t)eventdata.Data[kRowID];
      TVector3 SieveHoleTCS = GetSieveHoleTCS(Col,Row);
      
      const TVector3 MomDirectionTCS = SieveHoleTCS - BeamSpotTCS;

      
      TVector3 Tg_YSpotHCS=fTCSInHCS*Tg_YSpotTCS+fPointingOffset;
      TVector3 MomDirectionHCS=fTCSInHCS*MomDirectionTCS;
      
      //      assert(Tg_YSpotHCS.Y()==MissPointY);//internal consistency check
      
      const Int_t a = (HRSAngle > 0) ? 1 : -1;
      
      Double_t Real_Tg_Phi = MomDirectionTCS.Y() / MomDirectionTCS.Z();
      
      reactz = - ( eventdata.Data[kCalcTgY] -a*MissPointZ)*TMath::Cos(TMath::ATan(Real_Tg_Phi))/TMath::Sin(HRSAngle + TMath::ATan(Real_Tg_Phi)) + BeamSpotHCS.X()*TMath::Cos(HRSAngle+TMath::ATan(Real_Tg_Phi))/TMath::Sin(HRSAngle+TMath::ATan(Real_Tg_Phi));
      
      

      //      reactz = 
	// (  - ( eventdata.Data[kCalcTgY] + MissPointY)   +   BeamSpotHCS.X()*TMath::Cos(HRSAngle+TMath::ATan(eventdata.Data[kRealPhi]))  )/TMath::Sin(HRSAngle + TMath::ATan(eventdata.Data[kRealPhi]));




      eventdata.Data[kCalcReactZ] = reactz;
      
      dtg_vz +=  reactz - eventdata.Data[kRealReactZ];
      dtg_vz_rms += (reactz - eventdata.Data[kRealReactZ])*(reactz - eventdata.Data[kRealReactZ]);
      



	  
      //      hYVSReactZ->Fill(reactz,eventdata.Data[kY]);
      //      hPhVSReactZ->Fill(reactz,eventdata.Data[kPhi]);
      hrbxVSReactZ->Fill(reactz,eventdata.Data[kBeamX]);
      hrbyVSReactZ->Fill(reactz,eventdata.Data[kBeamY]);
      //      hReactZ_diff->Fill(dtg_vz);
      hReactZ_diff->Fill(reactz - eventdata.Data[kRealReactZ]);
      
      
      
      //statics by foils
      //	  const UInt_t FoilID = (UInt_t)eventdata.Data[kFoilID];
      // ndata_foil[FoilID]++;
      
      // hReacts_foils[FoilID]->Fill(reactz);
      
      // dtg_vz_foil[FoilID] += eventdata.Data[kCalcReactZ] - eventdata.Data[kRealReactZ];
      // rmstg_vz_foil[FoilID] += (eventdata.Data[kCalcReactZ] - eventdata.Data[kRealReactZ])*(eventdata.Data[kCalcReactZ] - eventdata.Data[kRealReactZ]);
      
      // dtg_y_foil[FoilID] += eventdata.Data[kCalcTgY] - eventdata.Data[kRealTgY];
      // 	  rmstg_y_foil[FoilID] += (eventdata.Data[kCalcTgY] - eventdata.Data[kRealTgY])*(eventdata.Data[kCalcTgY] - eventdata.Data[kRealTgY]);
	  
	}
	
  //	DEBUG_INFO("CheckVertex","dtg_vz = %f,\t dtg_vz_rms = %f",dtg_vz/fNRawData,dtg_vz_rms/fNRawData);
	
	TCanvas * c_vert_diff = new TCanvas("CheckVertexDiff","SieveCheckDiff",1800,900);
	c_vert_diff->Divide(1,3);





	// UInt_t idx=1;
	
	// c_vertcheck_res->cd(idx++);
	// hYVSReactZ->Draw("COLZ");
	// c_vertcheck_res->cd(idx++);
	// hPhVSReactZ->Draw("COLZ");
	// c_vertcheck_res->cd(idx++);
	// hrbxVSReactZ->Draw("COLZ");
	// c_vertcheck_res->cd(idx++);
	// hReactZ->Draw();

	TF1* foil_fits;
	TF1* foil_fits_2;
	


	//	TPaveStats *st[NFiols];

	
	//  JW: very temporary change to looping through 2 foils
	
	c_vert_diff->cd(1);
	

	// create stats object for reactz diff plot
	TPaveStats *st[NFiols];

	  //	  hReacts_foils[FoilID]->Fit(
	foil_fits = new TF1("Foil fit 1","gaus",-0.2,0.2);
	hReactZ_diff->Fit("Foil fit 1","R0");

	foil_fits_2 = new TF1("Foil fit 2","gaus",foil_fits->GetParameter(1)-2.*foil_fits->GetParameter(2),foil_fits->GetParameter(1)+2.*foil_fits->GetParameter(2));
	hReactZ_diff->Fit("Foil fit 2","R");	  	
			   
			     
			   
	hReactZ_diff->Draw();
	if(hReactZ_diff->GetEntries()!=0){
	  cout << "entered function drawing " << endl;
	  foil_fits_2->SetLineColor(kRed);
	  foil_fits_2->Draw("same");
	}
	

	//	st[FoilID] = (static_cast<TPaveStats*>((hReacts_foils[FoilID]->FindObject("stats"))));

	



	gPad->Modified();
	gPad->Update(); // make sure it’s (re)drawn
	


	
	c_vert_diff->cd(2);
	hrbxVSReactZ->Draw("col");

	c_vert_diff->cd(3);
	hrbyVSReactZ->Draw("col");



	
	
	
	gPad->Modified(); gPad->Update(); // make sure it’s (re)drawn

	
	//	  st[FoilID] = dynamic_cast<TPaveStats>(*(hReacts_foils[FoilID]->FindObject("stats")));
	//static_cast<THaSpectrometer*>
	
	
	
	


	

// 	    TF1* f2a=new TF1("f2a","gaus",bin_gaus_min,bin_gaus_max);
// 	  f2a->SetLineColorAlpha(3,1.0);
// 	  TF1* f31=new TF1("f31","pol1",bin_lin_min,bin_lin_max);
// 	  f31->SetLineColorAlpha(2,1.0);
	  
// 	  TF1* total=new TF1("total","gaus(0)+pol1(3)",bin_fit_min,bin_fit_max);
	  
// 	  total->SetLineColor(1);
// 	  double par[5];
// //f2a = ht1->GetFunction("gaus");
// 	  ht1->Fit(f2a,"R");
// 	  ht1->Fit(f31,"R+");
// 	  f2a->GetParameters(&par[0]);


	
	// for(idx = 1; idx<=4; idx++)
	//   {
	//     c_vertcheck_res->cd(idx);
	    
	//     for(UInt_t FoilID = 0; FoilID<NFiols; FoilID++)
	//       {
	// 	const Double_t MaxPlot = 2000;
		
	// 	TLine *l = new TLine(targetfoils[FoilID],-MaxPlot,targetfoils[FoilID],+MaxPlot);
	// 	l->SetLineColor(6);
	// 	l->Draw();
	//       }
	//   }
	
	cout<<"Statistic on each foils:\n";
	cout<<"Foil ID\td_vz\trms_vz\td_tgy\trms_tgy\tNev\n";
	for(UInt_t FoilID = 0; FoilID<NFiols; FoilID++)
	  {
	    // 		assert(ndata_foil[FoilID]);//at least 1 event on each foil
	    UInt_t N = ndata_foil[FoilID];
	    if (N<=0) N=1;
	    
	    printf("%d\t%f\t%f\t%f\t%f\t%d\n",FoilID
		   ,dtg_vz_foil[FoilID]/N
		   ,rmstg_vz_foil[FoilID]/N
		   ,dtg_y_foil[FoilID]/N
		   ,rmstg_y_foil[FoilID]/N
		   ,N);
	  }
	
	//	return c_vert_diff;


	// print canvasses
		
	if(print != 0){

	  c_vert_diff->Print(fDestDataBaseName + ".VertexDiff.Opt.png", "png"); 

	}



}

TCanvas * LOpticsOpt::CheckDp()
{
	//Visualize 1D hitogram of dp_kin

	DEBUG_INFO("CheckDp","Entry Point");

	//calculate Data[kCalcDpKin] for all events
	SumSquareDp(kTRUE);

	const Double_t DpRange=.05;
	const UInt_t NDpRange=5000;

	TH1D * hDpKinCalib[NKine];
	TH1D * hDpKinAll[NKine];
	Double_t RealDpKin[NKine];
	Double_t AverCalcDpKin[NKine]={0};
	UInt_t NEvntDpKin[NKine]={0};
	Double_t RealDpKinAllExcit[NExcitationStates][NKine];
	Double_t NewArbitaryDpKinShift[NKine];

	for(UInt_t KineID=0; KineID<NKine; KineID++)
	{
		hDpKinCalib[KineID]
				= new TH1D(Form("hDpKinCalib%d",KineID)
				,Form("Dp_Kin for Delta Scan Kine. #%d (Selected Exct. State)",KineID)
				,NDpRange,-DpRange,DpRange);
		hDpKinAll[KineID]
				= new TH1D(Form("hDpKinAll%d",KineID)
				,Form("Dp_Kin for Delta Scan Kine. #%d (All Data)",KineID)
				,NDpRange,-DpRange,DpRange);

		assert(hDpKinCalib[KineID]);//pointer check
		assert(hDpKinAll[KineID]);//pointer check
	}

	for (UInt_t idx = 0; idx<fNRawData; idx++)
	{
		const EventData &eventdata = fRawData[idx];
		const UInt_t ExtraDataFlag = (UInt_t)(eventdata.Data[kExtraDataFlag]);
		assert(ExtraDataFlag==0 || ExtraDataFlag==1); //flag definition consistency check
		const UInt_t KineID = (UInt_t)(eventdata.Data[kKineID]);
		assert(KineID<NKine);

		if (!ExtraDataFlag)
		{
			hDpKinCalib[KineID]->Fill(eventdata.Data[kCalcDpKin]+eventdata.Data[kRadiLossDp]);
			AverCalcDpKin[KineID]+=
					eventdata.Data[kCalcDpKin]+eventdata.Data[kRadiLossDp];
			NEvntDpKin[KineID]++;
		}
		hDpKinAll[KineID]->Fill(eventdata.Data[kCalcDpKin]+eventdata.Data[kRadiLossDp]);

		RealDpKin[KineID] = eventdata.Data[kRealDpKin]+eventdata.Data[kRadiLossDp];

		for (UInt_t ExcitID = 0; ExcitID<NExcitationStates; ExcitID++)
		{
			assert(kRealDpKinExcitations+ExcitID<kRealTh);//index check
			RealDpKinAllExcit[ExcitID][KineID]
					=eventdata.Data[kRealDpKinExcitations+ExcitID];
		}
	}

	TCanvas * c1 = new TCanvas("CheckDp","Check Dp Kin Reconstruction",1800,900);

	if (NKine<=3)
		c1->Divide(3,1);
	else if (NKine<=6)
		c1->Divide(3,2);
	else
		c1->Divide(3,3);

	UInt_t idx=1;

	for(UInt_t KineID=0; KineID<NKine; KineID++)
	{
		c1->cd(idx++);
		gPad -> SetLogy();

		AverCalcDpKin[KineID]/=NEvntDpKin[KineID];
		DEBUG_MASSINFO("CheckDp","AverCalcDpKin[%d] = %f"
				,KineID,AverCalcDpKin[KineID]);

		// Histograms
		hDpKinCalib[KineID]->SetLineColor(4);
		hDpKinCalib[KineID]->SetFillColor(4);
		hDpKinCalib[KineID]->SetFillStyle(3008);

		hDpKinAll[KineID]->SetLineColor(1);
// 		hDpKinAll[KineID]->SetFillColor(1);

		const Double_t dpRange=0.01;
		hDpKinCalib[KineID]->
				SetAxisRange(RealDpKin[KineID]-dpRange,RealDpKin[KineID]+dpRange);
		hDpKinAll[KineID]->
				SetAxisRange(RealDpKin[KineID]-dpRange,RealDpKin[KineID]+dpRange);

		hDpKinCalib[KineID]
				->SetXTitle("radiation corrected dp_kin (angular independant dp)");
		hDpKinAll[KineID]
				->SetXTitle("radiation corrected dp_kin (angular independant dp)");

		hDpKinAll[KineID]->Draw();
		hDpKinCalib[KineID]->Draw("SAME");

		// expectation lines
		const Double_t MaxPlot = 20000;
		for (UInt_t ExcitID = 0; ExcitID<NExcitationStates; ExcitID++)
		{
			const Double_t x = RealDpKinAllExcit[ExcitID][KineID];
			TLine *l = new TLine(x,0,x,+MaxPlot);
			l->SetLineColor(3);
			l->SetLineWidth(2);
			l->Draw();
		}
		TLine *l = new TLine(RealDpKin[KineID],0,RealDpKin[KineID],+MaxPlot);
		l->SetLineColor(6);
		l->SetLineWidth(2);
		l->Draw();


		// Fits
		const Double_t DefResolution = 1e-4;
		const Double_t FitRangeMultiply = 5;

		TString FitFunc = Form("DpPeak%d",KineID);
		TF1 *f = new TF1(FitFunc,"gaus+[3]+[4]*x"
				,AverCalcDpKin[KineID]-DefResolution*FitRangeMultiply
						,AverCalcDpKin[KineID]+DefResolution*FitRangeMultiply);
		f->SetParameter(1,AverCalcDpKin[KineID]);
		f->SetParameter(2,DefResolution);
		hDpKinAll[KineID] -> Fit(FitFunc,"RN0");
// 		Info("CheckDp","Fit for delta scan #%d peak:",KineID);
// 		f->Print();
		f->SetLineColor(2);
		f->Draw("SAME");

		TLatex *t = new TLatex(f->GetParameter(1)+DefResolution
				,f->GetParameter(0)+f->GetParameter(3)+f->GetParameter(4)*f->GetParameter(1)
				,Form("\\Delta \\pm \\sigma = %2.1f \\pm %2.1f \\times 10^{-4}"
						,10000*(f->GetParameter(1)-RealDpKin[KineID])
								,10000*f->GetParameter(2)));
		t->SetTextSize(0.05);
		t->SetTextAlign(12);
		t->SetTextColor(2);
		t->Draw();

		NewArbitaryDpKinShift[KineID] =
				f->GetParameter(1)-RealDpKin[KineID]+fArbitaryDpKinShift[KineID];

	}

	Info("CheckDp","New set of arbitary dp shifts:");
	for(UInt_t KineID=0; KineID<NKine; KineID++)
		printf("opt->fArbitaryDpKinShift[%d] = %e;\n"
				,KineID,NewArbitaryDpKinShift[KineID]);

	return c1;

}

TCanvas * LOpticsOpt::CheckDpGlobal()
{
	//Visualize 1D hitogram of dp_kin

	DEBUG_INFO("CheckDp","Entry Point");

	//calculate Data[kCalcDpKin] for all events
	SumSquareDp(kTRUE);

	const Double_t DpRange=.05;
	const UInt_t NDpRange=5000;

	TH1D * hDpKinCalib[NKine];
	TH1D * hDpKinAll[NKine];
	Double_t RealDpKin[NKine];
	Double_t AverCalcDpKin[NKine]={0};
	UInt_t NEvntDpKin[NKine]={0};
	Double_t RealDpKinAllExcit[NExcitationStates][NKine];
//	Double_t NewArbitaryDpKinShift[NKine];

	TH1D * hSumDpKin = new TH1D("hSumDpKin"
			,"Overlay of dp_{kin} data from elastic peak scans"
			,NDpRange,-DpRange,DpRange);
	TGraphErrors * geSum = new TGraphErrors(NKine);

	for(UInt_t KineID=0; KineID<NKine; KineID++)
	{
		hDpKinCalib[KineID]
				= new TH1D(Form("hDpKinCalib%d",KineID)
				,Form("Dp_Kin for Delta Scan Kine. #%d (Selected Exct. State)",KineID)
				,NDpRange,-DpRange,DpRange);
		hDpKinAll[KineID]
				= new TH1D(Form("hDpKinAll%d",KineID)
				,Form("Dp_Kin for Delta Scan Kine. #%d (All Data)",KineID)
				,NDpRange,-DpRange,DpRange);

		assert(hDpKinCalib[KineID]);//pointer check
		assert(hDpKinAll[KineID]);//pointer check
	}

	for (UInt_t idx = 0; idx<fNRawData; idx++)
	{
		const EventData &eventdata = fRawData[idx];
		const UInt_t ExtraDataFlag = (UInt_t)(eventdata.Data[kExtraDataFlag]);
		assert(ExtraDataFlag==0 || ExtraDataFlag==1); //flag definition consistency check
		const UInt_t KineID = (UInt_t)(eventdata.Data[kKineID]);
		assert(KineID<NKine);

		if (!ExtraDataFlag)
		{
			hDpKinCalib[KineID]->Fill(eventdata.Data[kCalcDpKin]+eventdata.Data[kRadiLossDp]);
			AverCalcDpKin[KineID]+=
					eventdata.Data[kCalcDpKin]+eventdata.Data[kRadiLossDp];
			NEvntDpKin[KineID]++;
		}
		hDpKinAll[KineID]->Fill(eventdata.Data[kCalcDpKin]+eventdata.Data[kRadiLossDp]);

		RealDpKin[KineID] = eventdata.Data[kRealDpKin]+eventdata.Data[kRadiLossDp];

		for (UInt_t ExcitID = 0; ExcitID<NExcitationStates; ExcitID++)
		{
			assert(kRealDpKinExcitations+ExcitID<kRealTh);//index check
			RealDpKinAllExcit[ExcitID][KineID]
					=eventdata.Data[kRealDpKinExcitations+ExcitID];
		}
	}


	//fit & sum
	for(UInt_t KineID=0; KineID<NKine; KineID++)
	{
		AverCalcDpKin[KineID]/=NEvntDpKin[KineID];
		DEBUG_MASSINFO("CheckDp","AverCalcDpKin[%d] = %f"
				,KineID,AverCalcDpKin[KineID]);



		// Fits
		const Double_t DefResolution = 1e-4;
		const Double_t FitRangeMultiply = 5;

		TString FitFunc = Form("DpPeak%d",KineID);
		TF1 *f = new TF1(FitFunc,"gaus+[3]+[4]*x"
				,AverCalcDpKin[KineID]-DefResolution*FitRangeMultiply
						,AverCalcDpKin[KineID]+DefResolution*FitRangeMultiply);
		f->SetParameter(1,AverCalcDpKin[KineID]);
		f->SetParameter(2,DefResolution);
		hDpKinAll[KineID] -> Fit(FitFunc,"RN0");
// 		Info("CheckDp","Fit for delta scan #%d peak:",KineID);
// 		f->Print();
		f->SetLineColor(2);

		hSumDpKin -> Add(hDpKinAll[KineID]);
//		hSumDpKin -> Add(hDpKinCalib[KineID]);

		geSum -> GetX()[KineID] = RealDpKin[KineID];
		geSum -> GetY()[KineID] = (f->GetParameter(1)-RealDpKin[KineID]);
		geSum -> GetEY()[KineID] = f->GetParameter(2);

	}

	TCanvas * c1 = new TCanvas("CheckDpGlobal","Check Dp Kin Reconstruction",1600,900);
	c1->Divide(1,2);
	UInt_t idx=1;

	const Double_t xlim1=-.005, xlim2=.045, ddplim=.0005;

	c1->cd(idx++);
	hSumDpKin->GetXaxis()->SetRangeUser(xlim1,xlim2);
	hSumDpKin->GetXaxis()->SetTitle("dp_{kin} (Angular indep. mom. dev.)");
	hSumDpKin->Draw();

	for(UInt_t KineID=0; KineID<NKine; KineID++)
	{
		TLine *l = new TLine(RealDpKin[KineID],hSumDpKin->GetMinimum(),RealDpKin[KineID],hSumDpKin->GetMaximum()*1.1);
		l->SetLineColor(kBlue);
		l->Draw();
	}


	c1->cd(idx++);
	TH1 * axis = TVirtualPad::Pad()->DrawFrame(xlim1,-ddplim,xlim2,ddplim,"");
//	axis->SetTitle("Summary of reconstructed elastic peaks");
	axis->GetXaxis()->SetTitle("dp_{kin} (Angular indep. mom. dev.)");
	axis->GetYaxis()->SetTitle("offset and width of elastic peaks");

	geSum->SetMarkerStyle(3);
	geSum->Draw("pe");

	TLine *l = new TLine(xlim1,0,xlim2,0);
	l->SetLineColor(kBlack);
	l->Draw();

	return c1;

}

TCanvas * LOpticsOpt::CheckDpVSAngle()
{

	//Visualize 2D hitogram of Scattering Angle VS dp_kin

	DEBUG_INFO("CheckDpVSAngle","Entry Point");

	//calculate Data[kCalcDpKin] for all events
	SumSquareDp(kTRUE);

	const Double_t DpRange=.05;
	const UInt_t NDpRange=800*5;
	const Double_t AngleRange=2;
	const UInt_t NAngleRange=100*5;

	TH2D * hDpKinCalib[NKine];
	TH2D * hDpKinAll[NKine];
	Double_t RealDpKin[NKine];
// 	Double_t AverCalcDpKin[NKine]={0};
// 	UInt_t NEvntDpKin[NKine]={0};
	Double_t RealDpKinAllExcit[NExcitationStates][NKine];
// 	Double_t NewArbitaryDpKinShift[NKine];

	for(UInt_t KineID=0; KineID<NKine; KineID++)
	{
		hDpKinCalib[KineID]
				= new TH2D(Form("hDpKinCalibVSAngle%d",KineID)
				,Form("Dp_Kin for Delta Scan Kine. #%d (Selected Exct. State)",KineID)
				,NDpRange,-DpRange,DpRange
						,NAngleRange,TMath::Abs(HRSAngle)/TMath::Pi()*180-AngleRange,TMath::Abs(HRSAngle)/TMath::Pi()*180+AngleRange);
		hDpKinAll[KineID]
				= new TH2D(Form("hDpKinAllVSAngle%d",KineID)
				,Form("Dp_Kin for Delta Scan Kine. #%d (All Data)",KineID)
				,NDpRange,-DpRange,DpRange
						,NAngleRange,TMath::Abs(HRSAngle)/TMath::Pi()*180-AngleRange,TMath::Abs(HRSAngle)/TMath::Pi()*180+AngleRange);

		assert(hDpKinCalib[KineID]);//pointer check
		assert(hDpKinAll[KineID]);//pointer check
	}

	for (UInt_t idx = 0; idx<fNRawData; idx++)
	{
		const EventData &eventdata = fRawData[idx];
		const UInt_t ExtraDataFlag = (UInt_t)(eventdata.Data[kExtraDataFlag]);
		assert(ExtraDataFlag==0 || ExtraDataFlag==1); //flag definition consistency check
		const UInt_t KineID = (UInt_t)(eventdata.Data[kKineID]);
		assert(KineID<NKine);

		if (!ExtraDataFlag)
		{
			hDpKinCalib[KineID]
					->Fill(eventdata.Data[kCalcDpKin]+eventdata.Data[kRadiLossDp]
					,+eventdata.Data[kScatterAngle]/TMath::Pi()*180);
// 			AverCalcDpKin[KineID]+=
// 					eventdata.Data[kCalcDpKin]+eventdata.Data[kRadiLossDp];
// 			NEvntDpKin[KineID]++;
		}
		hDpKinAll[KineID]
				->Fill(eventdata.Data[kCalcDpKin]+eventdata.Data[kRadiLossDp]
				,+eventdata.Data[kScatterAngle]/TMath::Pi()*180);

		RealDpKin[KineID] = eventdata.Data[kRealDpKin]+eventdata.Data[kRadiLossDp];

		for (UInt_t ExcitID = 0; ExcitID<NExcitationStates; ExcitID++)
		{
			assert(kRealDpKinExcitations+ExcitID<kRealTh);//index check
			RealDpKinAllExcit[ExcitID][KineID]
					=eventdata.Data[kRealDpKinExcitations+ExcitID];
		}
	}

	TCanvas * c1 = new TCanvas("CheckDpVSAngle","Check Dp Kin Reconstruction VS Scattering Angle",1800,900);

	if (NKine<=3)
		c1->Divide(3,1);
	else if (NKine<=6)
		c1->Divide(3,2);
	else
		c1->Divide(3,3);

	UInt_t idx=1;

	for(UInt_t KineID=0; KineID<NKine; KineID++)
	{
		c1->cd(idx++);
// 		gPad -> SetLogy();

// 		AverCalcDpKin[KineID]/=NEvntDpKin[KineID];
// 		DEBUG_MASSINFO("CheckDp","AverCalcDpKin[%d] = %f"
// 				,KineID,AverCalcDpKin[KineID]);

		// Histograms
// 		hDpKinCalib[KineID]->SetLineColor(4);
// 		hDpKinCalib[KineID]->SetFillColor(4);
// 		hDpKinCalib[KineID]->SetFillStyle(3008);

// 		hDpKinAll[KineID]->SetLineColor(1);
// 		hDpKinAll[KineID]->SetFillColor(1);

		const Double_t dpRange=0.01;
		hDpKinCalib[KineID]->
				SetAxisRange(RealDpKin[KineID]-dpRange,RealDpKin[KineID]+dpRange);
		hDpKinAll[KineID]->
				SetAxisRange(RealDpKin[KineID]-dpRange,RealDpKin[KineID]+dpRange);

		hDpKinCalib[KineID]
				->SetXTitle("radiation corrected dp_kin (angular independant dp)");
		hDpKinCalib[KineID]
				->SetYTitle("Scattering Angle (Degree)");
		hDpKinAll[KineID]
				->SetXTitle("radiation corrected dp_kin (angular independant dp)");
		hDpKinAll[KineID]
				->SetYTitle("Scattering Angle (Degree)");

		hDpKinAll[KineID]->Draw("COLZ");
		TVirtualPad::Pad()->SetLogz();
// 		hDpKinCalib[KineID]->Draw("SAME");

		// expectation lines
		const Double_t MinPlot = TMath::Abs(HRSAngle)/TMath::Pi()*180-AngleRange;
		const Double_t MaxPlot = TMath::Abs(HRSAngle)/TMath::Pi()*180+AngleRange;
		for (UInt_t ExcitID = 0; ExcitID<NExcitationStates; ExcitID++)
		{
			const Double_t x = RealDpKinAllExcit[ExcitID][KineID];
			TLine *l = new TLine(x,MinPlot,x,+MaxPlot);
			l->SetLineColor(3);
			l->SetLineWidth(2);
			l->Draw();
		}
		TLine *l = new TLine(RealDpKin[KineID],MinPlot,RealDpKin[KineID],+MaxPlot);
		l->SetLineColor(6);
		l->SetLineWidth(2);
		l->Draw();
	}


	return c1;

}

// TCanvas * LOpticsOpt::CheckDpVSCutID()
// {
// 	//Visualize 2D hitogram of Sieve Hole+Foil ID VS dp_kin

// 	DEBUG_INFO("CheckDpVSCutID","Entry Point");

// 	//calculate Data[kCalcDpKin] for all events
// 	SumSquareDp(kTRUE);

// 	const Double_t DpRange=.05;
// 	const UInt_t NDpRange=800*5;

// 	const UInt_t NCutID = NSieveRow*NSieveCol*NFiols;
// 	const Double_t CutIDRangeLow = -.5;
// 	const Double_t CutIDRangeHigh= -.5+NCutID;

// 	TH2D * hDpKinCalib[NKine];
// 	TH2D * hDpKinAll[NKine];
// 	Double_t RealDpKin[NKine];
// // 	Double_t AverCalcDpKin[NKine]={0};
// // 	UInt_t NEvntDpKin[NKine]={0};
// 	Double_t RealDpKinAllExcit[NExcitationStates][NKine];
// // 	Double_t NewArbitaryDpKinShift[NKine];

// 	for(UInt_t KineID=0; KineID<NKine; KineID++)
// 	{
// 		hDpKinCalib[KineID]
// 				= new TH2D(Form("hDpKinCalibVSCutID%d",KineID)
// 				,Form("Dp_Kin for Delta Scan Kine. #%d (Selected Exct. State)",KineID)
// 				,NDpRange,-DpRange,DpRange
// 						,NCutID,CutIDRangeLow,CutIDRangeHigh);
// 		hDpKinAll[KineID]
// 				= new TH2D(Form("hDpKinAllVSCutID%d",KineID)
// 				,Form("Dp_Kin for Delta Scan Kine. #%d (All Data)",KineID)
// 				,NDpRange,-DpRange,DpRange
// 						,NCutID,CutIDRangeLow,CutIDRangeHigh);

// 		assert(hDpKinCalib[KineID]);//pointer check
// 		assert(hDpKinAll[KineID]);//pointer check
// 	}

// 	for (UInt_t idx = 0; idx<fNRawData; idx++)
// 	{
// 		const EventData &eventdata = fRawData[idx];
// 		const UInt_t ExtraDataFlag = (UInt_t)(eventdata.Data[kExtraDataFlag]);
// 		assert(ExtraDataFlag==0 || ExtraDataFlag==1); //flag definition consistency check
// 		const UInt_t KineID = (UInt_t)(eventdata.Data[kKineID]);
// 		assert(KineID<NKine);

// 		if (!ExtraDataFlag)
// 		{
// // 			assert((UInt_t)eventdata.Data[kCutID]<NCutID); // cut definition check
// 			hDpKinCalib[KineID]
// 					->Fill(eventdata.Data[kCalcDpKin]+eventdata.Data[kRadiLossDp]
// 					,((UInt_t)eventdata.Data[kCutID])%NCutID);
// // 			AverCalcDpKin[KineID]+=
// // 					eventdata.Data[kCalcDpKin]+eventdata.Data[kRadiLossDp];
// // 			NEvntDpKin[KineID]++;
// 		}
// 		hDpKinAll[KineID]
// 				->Fill(eventdata.Data[kCalcDpKin]+eventdata.Data[kRadiLossDp]
// 				,((UInt_t)eventdata.Data[kCutID])%NCutID);

// 		RealDpKin[KineID] = eventdata.Data[kRealDpKin]+eventdata.Data[kRadiLossDp];

// 		for (UInt_t ExcitID = 0; ExcitID<NExcitationStates; ExcitID++)
// 		{
// 			assert(kRealDpKinExcitations+ExcitID<kRealTh);//index check
// 			RealDpKinAllExcit[ExcitID][KineID]
// 					=eventdata.Data[kRealDpKinExcitations+ExcitID];
// 		}
// 	}

// 	TCanvas * c1 = new TCanvas("CheckDpVSCutID","Check Dp Kin Reconstruction VS Scattering Angle",1800,1100);

// 	if (NKine<=3)
// 		c1->Divide(3,1);
// 	else if (NKine<=6)
// 		c1->Divide(3,2);
// 	else
// 		c1->Divide(3,3);


// 	c1->Update();
// 	UInt_t idx=1;

// 	for(UInt_t KineID=0; KineID<NKine; KineID++)
// 	{
// 		c1->cd(idx++);
// // 		gPad -> SetLogy();

// // 		AverCalcDpKin[KineID]/=NEvntDpKin[KineID];
// // 		DEBUG_MASSINFO("CheckDp","AverCalcDpKin[%d] = %f"
// // 				,KineID,AverCalcDpKin[KineID]);

// 		// Histograms
// // 		hDpKinCalib[KineID]->SetLineColor(4);
// // 		hDpKinCalib[KineID]->SetFillColor(4);
// // 		hDpKinCalib[KineID]->SetFillStyle(3008);

// // 		hDpKinAll[KineID]->SetLineColor(1);
// // 		hDpKinAll[KineID]->SetFillColor(1);

// 		const Double_t dpRange=0.01;
// 		hDpKinCalib[KineID]->
// 				SetAxisRange(RealDpKin[KineID]-dpRange,RealDpKin[KineID]+dpRange);
// 		hDpKinAll[KineID]->
// 				SetAxisRange(RealDpKin[KineID]-dpRange,RealDpKin[KineID]+dpRange);

// 		hDpKinCalib[KineID]
// 				->SetXTitle("radiation corrected dp_kin (angular independant dp)");
// 		hDpKinCalib[KineID]
// 				->SetYTitle("combination of foil and sieve hole ID");
// 		hDpKinAll[KineID]
// 				->SetXTitle("radiation corrected dp_kin (angular independant dp)");
// 		hDpKinAll[KineID]
// 				->SetYTitle("combination of foil and sieve hole ID");

// 		hDpKinAll[KineID]->Draw("COLZ");
// 		TVirtualPad::Pad()->SetLogz();
// // 		hDpKinCalib[KineID]->Draw("SAME");

// 		// expectation lines
// 		const Double_t MinPlot = CutIDRangeLow;
// 		const Double_t MaxPlot = CutIDRangeHigh;
// 		for (UInt_t ExcitID = 0; ExcitID<NExcitationStates; ExcitID++)
// 		{//all exciation states
// 			const Double_t x = RealDpKinAllExcit[ExcitID][KineID];
// 			TLine *l = new TLine(x,MinPlot,x,+MaxPlot);
// 			l->SetLineColor(3);
// 			l->SetLineWidth(2);
// 			l->Draw();
// 		}

// //		for(UInt_t FoilID = 0; FoilID<NFiols; FoilID++)
// //		{//Foils Divide
// //			const Double_t y = FoilID*NSieveRow*NSieveCol-.5;
// //			TLine *l = new TLine(RealDpKin[KineID]-dpRange,y,RealDpKin[KineID]+dpRange,y);
// //			l->SetLineColor(5);
// //// 			l->SetLineStyle(2);
// //			l->SetLineWidth(1);
// //			l->Draw();
// //		}
// 		for(UInt_t Col = 0; Col<NSieveCol; Col++)
// 		{//Foils Divide
// 			const Double_t y = Col*NSieveRow-.5;
// 			TLine *l = new TLine(RealDpKin[KineID]-dpRange,y,RealDpKin[KineID]+dpRange,y);
// 			l->SetLineColor(5);
//  			l->SetLineStyle(2);
// 			l->SetLineWidth(1);
// 			l->Draw();
// 		}

// 		TLine *l = new TLine(RealDpKin[KineID],MinPlot,RealDpKin[KineID],+MaxPlot);
// 		l->SetLineColor(6);
// 		l->SetLineWidth(2);
// 		l->Draw();
// 	}


// 	return c1;

// }

//_____________________________________________________________________________
Double_t LOpticsOpt::VerifyMatrix_Sieve(void)
{
	//static summarize difference between tg_th, th_ph caculated from current database and those in root file

	Double_t dth = 0, dphi=0;	//Difference
	Double_t rmsth = 0, rmsphi=0; //mean square

	Double_t theta, phi/*, dp, p, pathl*/;

	for (UInt_t idx = 0; idx<fNRawData; idx++)
	{
		const EventData eventdata = fRawData[idx];

		const Double_t x_fp = eventdata.Data[kX];
		const Double_t (*powers)[5] = eventdata.powers;

  // calculate the matrices we need
		CalcMatrix(x_fp, fDMatrixElems);
		CalcMatrix(x_fp, fTMatrixElems);
		CalcMatrix(x_fp, fYMatrixElems);
		CalcMatrix(x_fp, fYTAMatrixElems);
		CalcMatrix(x_fp, fPMatrixElems);
		CalcMatrix(x_fp, fPTAMatrixElems);

  // calculate the coordinates at the target
		theta = CalcTargetVar(fTMatrixElems, powers);
		phi = CalcTargetVar(fPMatrixElems, powers)+CalcTargetVar(fPTAMatrixElems,powers);
// 		y = CalcTargetVar(fYMatrixElems, powers)+CalcTargetVar(fYTAMatrixElems,powers);

  // calculate momentum
// 		dp = CalcTargetVar(fDMatrixElems, powers);

		dth += theta - eventdata.Data[kL_tr_tg_th];
		rmsth += (theta - eventdata.Data[kL_tr_tg_th])*(theta - eventdata.Data[kL_tr_tg_th]);

		dphi += phi - eventdata.Data[kL_tr_tg_ph];
		rmsphi += (phi - eventdata.Data[kL_tr_tg_ph])*(phi - eventdata.Data[kL_tr_tg_ph]);
	}

	DEBUG_INFO("VerifyMatrix_Sieve","dth = %f,rmsth=%f",
			   dth/fNRawData,TMath::Sqrt(rmsth/fNRawData));
	DEBUG_INFO("VerifyMatrix_Sieve","dphi = %f, rmsphi=%f",
			   dphi/fNRawData,TMath::Sqrt(rmsphi/fNRawData));

	return TMath::Sqrt(rmsth/fNRawData+rmsphi/fNRawData);
}

//_____________________________________________________________________________
Double_t LOpticsOpt::VerifyMatrix_Vertex(void)
{
	//static summarize difference between tg_y caculated from current database and those in root file

	Double_t /*dth = 0, */dtg_y=0;	//Difference
	Double_t /*rmsth = 0, */dtg_y_rms=0; //mean square

	Double_t y/*, dp, p*//*, pathl*/;

	for (UInt_t idx = 0; idx<fNRawData; idx++)
	{
		const EventData eventdata = fRawData[idx];

		Double_t x_fp = eventdata.Data[kX];
		const Double_t (*powers)[5] = eventdata.powers;

  // calculate the matrices we need
// 		CalcMatrix(x_fp, fDMatrixElems);
// 		CalcMatrix(x_fp, fTMatrixElems);
		CalcMatrix(x_fp, fYMatrixElems);
// 		CalcMatrix(x_fp, fYTAMatrixElems);
// 		CalcMatrix(x_fp, fPMatrixElems);
// 		CalcMatrix(x_fp, fPTAMatrixElems);

  // calculate the coordinates at the target
// 		theta = CalcTargetVar(fTMatrixElems, powers);
// 		phi = CalcTargetVar(fPMatrixElems, powers)+CalcTargetVar(fPTAMatrixElems,powers);
		y = CalcTargetVar(fYMatrixElems, powers)/*+CalcTargetVar(fYTAMatrixElems,powers)*/;

  // calculate momentum
// 		dp = CalcTargetVar(fDMatrixElems, powers);

		dtg_y += y - eventdata.Data[kL_tr_tg_y];
		dtg_y_rms += (y - eventdata.Data[kL_tr_tg_y])*(y - eventdata.Data[kL_tr_tg_y]);

		DEBUG_MASSINFO("VerifyMatrix_Vertex","y = %f; eventdata.Data[kL_tr_tg_y] = %f",
					   y , eventdata.Data[kL_tr_tg_y]);
	}

	DEBUG_INFO("VerifyMatrix_Vertex","dtg_y = %f,dtg_y_rms=%f",
			   dtg_y/fNRawData,TMath::Sqrt(dtg_y_rms/fNRawData));

	return TMath::Sqrt(dtg_y_rms/fNRawData);
}

//_____________________________________________________________________________
Double_t LOpticsOpt::VerifyMatrix_Dp(void)
{

	//static summarize difference between tg_dp caculated from current database and those in root file

	Double_t d_dp = 0/*, dphi=0*/;	//Difference
	Double_t rms_dp = 0/*, rmsphi=0*/; //mean square

	static UInt_t NCall = 0;
	NCall++;

	for (UInt_t idx = 0; idx<fNRawData; idx++)
	{

		Double_t  dp, dp_kin;

		EventData &eventdata = fRawData[idx];

		Double_t x_fp = eventdata.Data[kX];
		const Double_t (*powers)[5] = eventdata.powers;

  // calculate the matrices we need
		CalcMatrix(x_fp, fDMatrixElems);
// 		CalcMatrix(x_fp, fTMatrixElems);
// 		CalcMatrix(x_fp, fYMatrixElems);
// 		CalcMatrix(x_fp, fYTAMatrixElems);
// 		CalcMatrix(x_fp, fPMatrixElems);
// 		CalcMatrix(x_fp, fPTAMatrixElems);

  // calculate the coordinates at the target
// 		theta = CalcTargetVar(fTMatrixElems, powers);
// 		phi = CalcTargetVar(fPMatrixElems, powers)+CalcTargetVar(fPTAMatrixElems,powers);
// 		y = CalcTargetVar(fYMatrixElems, powers)+CalcTargetVar(fYTAMatrixElems,powers);

  // calculate momentum
		dp = CalcTargetVar(fDMatrixElems, powers);
		dp_kin = dp - eventdata.Data[kDpKinOffsets];

		d_dp += dp - eventdata.Data[kL_tr_tg_dp];
		rms_dp += (dp - eventdata.Data[kL_tr_tg_dp])*(dp - eventdata.Data[kL_tr_tg_dp]);

		DEBUG_MASSINFO("SumSquareDp","d_dp = %f = \t%f - \t%f",
					   dp_kin - eventdata.Data[kRealDpKin], dp_kin , eventdata.Data[kRealDpKin]  );

	}

	DEBUG_INFO("VerifyMatrix_Dp","(dp_Calc - dp in root tree) = d_dp = %f, rms_d_dp=%f",
			   d_dp/fNRawData,TMath::Sqrt(rms_dp/fNRawData));

	return TMath::Sqrt(rms_dp/fNRawData);
}

//_____________________________________________________________________________
LOpticsOpt::LOpticsOpt( const char* name, const char* description,
		THaApparatus* apparatus ) :
		THaTrackingDetector(name,description,apparatus)
{
  // Constructor

	fPrefix = new char [1000];
	sprintf(fPrefix, "%s", Prefix);

	fCurrentMatrixElems = NULL;

	TVector3 TCSX(0,-1,0);
	TVector3 TCSZ(TMath::Sin(HRSAngle),0,TMath::Cos(HRSAngle));
	TVector3 TCSY = TCSZ.Cross(TCSX);
	fTCSInHCS.RotateAxes(TCSX,TCSY,TCSZ);


	// original definition with this code
	// 	fPointingOffset.SetXYZ(-MissPointZ*TMath::Sin(HRSAngle)*TMath::Cos(HRSAngle),(Double_t)MissPointY,MissPointZ*TMath::Sin(HRSAngle)*TMath::Sin(HRSAngle));

	// seperate definition found in Sean's code
 	fPointingOffset.SetXYZ(-MissPointZ*TMath::Cos(HRSAngle),(Double_t)MissPointY,-MissPointZ*TMath::Sin(HRSAngle));

	//	fPointingOffset.SetXYZ(0,0,0);
	

	DEBUG_INFO("LOpticsOpt","Read in configuration "+InputID);

#ifdef DISABLED_SIEVE_HOLE
	DEBUG_WARNING("LOpticsOpt","Some of the sieve data are excluded from fitting");
#endif

	DEBUG_INFO("LOpticsOpt","HRS @ %f Degree, PointingOffset = (%f,%f,%f), SievePos = (%f,%f,%f)",
			   HRSAngle/TMath::Pi()*180,
			   fPointingOffset.X(),fPointingOffset.Y(),fPointingOffset.Z(),
			   SieveOffX, SieveOffY, ZPos);

	fNRawData=0;

	for(UInt_t i=0; i<100; i++)
		fArbitaryDpKinShift[i] = fArbitaryVertexShift[i] = 0;
}

//_____________________________________________________________________________
Int_t LOpticsOpt::LoadDataBase( TString DataBaseName )
{
	OldComments = "";

  // Read VDC database

	FILE* file = fopen( DataBaseName,"r" );
	if( !file ) {
		Error("LoadDataBase","%s can not be opened", DataBaseName.Data());
		assert(0);//
		 return kFileError;
	}
	else DEBUG_INFO("LoadDataBase","Parsing Database %s", DataBaseName.Data());

  // load global VDC parameters
	static const char* const here = "LoadDataBase";
	const int LEN = 200;
	char buff[LEN];

  // Look for the section [<prefix>.global] in the file, e.g. [ R.global ]
// 	TString tag(fPrefix);
// 	Ssiz_t pos = tag.Index(".");
// 	if( pos != kNPOS )
// 		tag = tag(0,pos+1);
// 	else
// 		tag.Append(".");
// 	tag.Prepend("[");
// 	tag.Append("global]");
// 	TString line, tag2(tag);
// 	tag.ToLower();

    TString tag(fPrefix);
    Ssiz_t tagpos = tag.Index(".");
    if (tagpos != kNPOS)
        tag = tag(0, tagpos + 1);
    else
        tag.Append(".");
    //    tag.Prepend("[");
    tag.Append("vdc.matrixelem=");
 	TString line, tag2(tag);
	//    TString tag2(tag);
    tag.ToLower();

	bool found = false;
	while (!found && fgets(buff, LEN, file) != NULL) {

		//read in comments
		TString tmpline = buff;
		if ( tmpline.BeginsWith("#"))
		{
			OldComments += tmpline;
//			OldComments += "\n";
		}

		char* buf = ::Compress(buff);  //strip blanks
		line = buf;
		delete [] buf;
		if( line.EndsWith("\n") ) line.Chop();


		line.ToLower();
		if ( tag == line )
			found = true;
	}
	if( !found ) {
		Error(Here(here), "Database entry %s not found!", tag2.Data() );
		fclose(file);
		assert(0);//
		return kInitError;
	}

  // We found the section, now read the data

  // read in some basic constants first
  //  fscanf(file, "%lf", &fSpacing);
  // fSpacing is calculated from the actual z-positions in Init()
	fgets(buff, LEN, file); // Skip rest of line	
	fgets(buff, LEN, file); // Skip comment line

	fTMatrixElems.clear();
	fDMatrixElems.clear();
	fPMatrixElems.clear();
	fPTAMatrixElems.clear();
	fYMatrixElems.clear();
	fYTAMatrixElems.clear();
	fLMatrixElems.clear();

	fFPMatrixElems.clear();
	fFPMatrixElems.resize(3);

	typedef vector<string>::size_type vsiz_t;
	map<string,vsiz_t> power;
	power["t"] = 3;  // transport to focal-plane tensors
	power["y"] = 3;
	power["p"] = 3;
	power["D"] = 3;  // focal-plane to target tensors
	power["T"] = 3;
	power["Y"] = 3; // (JW: Defined above references)
	power["YTA"] = 4;
	power["P"] = 3;
	power["PTA"] = 4;
	power["L"] = 4;  // pathlength from z=0 (target) to focal plane (meters) (JW: Defined above references)
	power["XF"] = 5; // forward: target to focal-plane (I think)
	power["TF"] = 5;
	power["PF"] = 5;
	power["YF"] = 5;




	// JW added to input order being optimised to
	// where order = i + j + k + l for optics elements
	map<string,Int_t> power_opt;
	power_opt["t"] = 4;  // transport to focal-plane tensors
	power_opt["y"] = 4;
	power_opt["p"] = 4;
	power_opt["D"] = 5;  // focal-plane to target tensors
	power_opt["T"] = 5;
	power_opt["Y"] = 5; 

	power_opt["YTA"] = 4;
	power_opt["P"] = 5;
	power_opt["PTA"] = 4;
	power_opt["L"] = 4;  

	power_opt["XF"] = 5; 
	power_opt["TF"] = 5;
	power_opt["PF"] = 5;
	power_opt["YF"] = 5;

	


	map<string,vector<THaMatrixElement>*> matrix_map;
	matrix_map["t"] = &fFPMatrixElems;
	matrix_map["y"] = &fFPMatrixElems;
	matrix_map["p"] = &fFPMatrixElems;
	matrix_map["D"] = &fDMatrixElems;
	matrix_map["T"] = &fTMatrixElems;
	matrix_map["Y"] = &fYMatrixElems;
	matrix_map["YTA"] = &fYTAMatrixElems;
	matrix_map["P"] = &fPMatrixElems;
	matrix_map["PTA"] = &fPTAMatrixElems;
	matrix_map["L"] = &fLMatrixElems;

	map <string,int> fp_map;
	fp_map["t"] = 0;
	fp_map["y"] = 1;
	fp_map["p"] = 2;

  // Read in as many of the matrix elements as there are.
  // Read in line-by-line, so as to be able to handle tensors of
  // different orders.
	while( fgets(buff, LEN, file) ) {
		string line(buff);
    // Erase trailing newline
		if( line.size() > 0 && line[line.size()-1] == '\n' ) {
			buff[line.size()-1] = 0;
			line.erase(line.size()-1,1);
		}
    // Split the line into whitespace-separated fields
		vector<string> line_spl = Split(line);

    // Stop if the line does not start with a string referring to
    // a known type of matrix element. In particular, this will
    // stop on a subsequent timestamp or configuration tag starting with "["
		if(line_spl.empty())
			continue; //ignore empty lines
		const char* w = line_spl[0].c_str();
		vsiz_t npow = power[w]; // JW line picks out 'L' or 'Y' and gets power as defined above (JW: defined above)
		if( npow == 0 )
			break;


#if DEBUG_LEVEL>=4
		cout<<"Matrix Line = ";
		for (Int_t pos=1; (UInt_t)pos<(UInt_t)line_spl.size(); pos++) {
			cout<< pos <<"("<<line_spl[pos].c_str()<<"), ";
		}
		cout<<endl;
#endif

    // Looks like a good line, go parse it.
		THaMatrixElement ME;
		ME.pw.resize(npow);
		ME.iszero = true;  ME.order = 0;
		vsiz_t pos;
		for (pos=1; pos<=npow && pos<line_spl.size(); pos++) {
			ME.pw[pos-1] = atoi(line_spl[pos].c_str());
		}
		vsiz_t p_cnt;
		for ( p_cnt=0; pos<line_spl.size() && p_cnt<kPORDER && pos<=npow+kPORDER;
					pos++,p_cnt++ )
		{
			ME.poly[p_cnt] = atof(line_spl[pos].c_str());
			if (ME.poly[p_cnt] != 0.0) {
				ME.iszero = false;
				ME.order = p_cnt+1;
			}
		}
		if (p_cnt < 1) {
			Error(Here(here), "Could not read in Matrix Element %s%d%d%d!",
				  w, ME.pw[0], ME.pw[1], ME.pw[2]);
			Error(Here(here), "Line looks like: %s",line.c_str());
			fclose(file);
			return kInitError;
		}


		// JW: Altered from old reading of OptOrder to setting above in this function
		
		// Olden way
		//order optimize to
		// ~~~~~~~~
		//		ME.OptOrder = atoi(line_spl[line_spl.size()-1].c_str());
		// ~~~~~~~~


		Int_t opt_order = power_opt[w];
		
		Int_t sum_powers = 0;

		for( Int_t pc = 0; pc < npow; pc++){ //pc acronym for power count (goes through i, j and k to get size of these)
		sum_powers  += ME.pw[pc];
		}

		ME.OptOrder = opt_order - sum_powers; 
		cout << "Type: " << w << " :  opt_order = " << opt_order << ", sum_powers = " << sum_powers << ", ME.OptOrder = " << ME.OptOrder << endl;




    // Don't bother with all-zero matrix elements
		if( ME.iszero )
			continue;

    // Add this matrix element to the appropriate array
		vector<THaMatrixElement> *mat = matrix_map[w];
		if (mat) {
      // Special checks for focal plane matrix elements
			if( mat == &fFPMatrixElems ) {
				if( ME.pw[0] == 0 && ME.pw[1] == 0 && ME.pw[2] == 0 ) {
					THaMatrixElement& m = (*mat)[fp_map[w]];
					if( m.order > 0 ) {
						Warning(Here(here), "Duplicate definition of focal plane "
								"matrix element: %s. Using first definition.", buff);
					} else
						m = ME;
				} else
					Warning(Here(here), "Bad coefficients of focal plane matrix "
							"element %s", buff);
			}
			else {
	// All other matrix elements are just appended to the respective array
	// but ensure that they are defined only once!
				bool match = false;
				for( vector<THaMatrixElement>::iterator it = mat->begin();
								 it != mat->end() && !(match = it->match(ME)); it++ ) {}
				if( match ) {
					Warning(Here(here), "Duplicate definition of "
							"matrix element: %s. Using first definition.", buff);
				} else
					mat->push_back(ME);
			}
		}
		else if ( fDebug > 0 )
			Warning(Here(here), "Not storing matrix for: %s !", w);

	} //while(fgets)

//   // Compute derived quantities and set some hardcoded parameters
//   const Double_t degrad = TMath::Pi()/180.0;
//   fTan_vdc  = fFPMatrixElems[T000].poly[0];
//   fVDCAngle = TMath::ATan(fTan_vdc);
//   fSin_vdc  = TMath::Sin(fVDCAngle);
//   fCos_vdc  = TMath::Cos(fVDCAngle);
	//
//   // Define the VDC coordinate axes in the "detector system". By definition,
//   // the detector system is identical to the VDC origin in the Hall A HRS.
//   DefineAxes(0.0*degrad);
	//
//   fNumIter = 1;      // Number of iterations for FineTrack()
//   fErrorCutoff = 1e100;
	//
//   // figure out the track length from the origin to the s1 plane
	//
//   // since we take the VDC to be the origin of the coordinate
//   // space, this is actually pretty simple
//   const THaDetector* s1 = GetApparatus()->GetDetector("s1");
//   if(s1 == NULL)
//     fCentralDist = 0;
//   else
//     fCentralDist = s1->GetOrigin().Z();

	CalcMatrix(1.,fLMatrixElems); // tensor without explicit polynomial in x_fp


	fIsInit = true;
	fclose(file);
	return kOK;
}

//_____________________________________________________________________________
LOpticsOpt::~LOpticsOpt()
{
  // Destructor.

}

//_____________________________________________________________________________
void LOpticsOpt::CalcTargetCoords(THaTrack *track, const ECoordTypes mode)
{
  // calculates target coordinates from focal plane coordinates

// 	const Int_t kNUM_PRECOMP_POW = 10;

	Double_t x_fp, y_fp, th_fp, ph_fp;
	Double_t powers[kNUM_PRECOMP_POW][5];  // {(x), th, y, ph, abs(th) }
	Double_t x, y, theta, phi, dp, p, pathl;

  // first select the coords to use
	if(mode == kTransport) {
		x_fp = track->GetX();
		y_fp = track->GetY();
		th_fp = track->GetTheta();
		ph_fp = track->GetPhi();
	} else {//if(mode == kRotatingTransport) {
		x_fp = track->GetRX();
		y_fp = track->GetRY();
		th_fp = track->GetRTheta();
		ph_fp = track->GetRPhi();
	}

  // calculate the powers we need
	for(int i=0; i<kNUM_PRECOMP_POW; i++) {
		powers[i][0] = pow(x_fp, i);
		powers[i][1] = pow(th_fp, i);
		powers[i][2] = pow(y_fp, i);
		powers[i][3] = pow(ph_fp, i);
		powers[i][4] = pow(TMath::Abs(th_fp),i);
	}

  // calculate the matrices we need
	CalcMatrix(x_fp, fDMatrixElems);
	CalcMatrix(x_fp, fTMatrixElems);
	CalcMatrix(x_fp, fYMatrixElems);
	CalcMatrix(x_fp, fYTAMatrixElems);
	CalcMatrix(x_fp, fPMatrixElems);
	CalcMatrix(x_fp, fPTAMatrixElems);

  // calculate the coordinates at the target
	theta = CalcTargetVar(fTMatrixElems, powers);
	phi = CalcTargetVar(fPMatrixElems, powers)+CalcTargetVar(fPTAMatrixElems,powers);
	y = CalcTargetVar(fYMatrixElems, powers)+CalcTargetVar(fYTAMatrixElems,powers);

	THaSpectrometer *app = static_cast<THaSpectrometer*>(GetApparatus());
  // calculate momentum
	dp = CalcTargetVar(fDMatrixElems, powers);
	p  = app->GetPcentral() * (1.0+dp);

  // pathlength matrix is for the Transport coord plane
	pathl = CalcTarget2FPLen(fLMatrixElems, powers);

  //FIXME: estimate x ??
	x = 0.0;

  // Save the target quantities with the tracks
	track->SetTarget(x, y, theta, phi);
	track->SetDp(dp);
	track->SetMomentum(p);
	track->SetPathLen(pathl);

	app->TransportToLab( p, theta, phi, track->GetPvect() );

}

//_____________________________________________________________________________
void LOpticsOpt::CalcMatrix( const Double_t x, vector<THaMatrixElement>& matrix )
{
  // calculates the values of the matrix elements for a given location
  // by evaluating a polynomial in x of order it->order with
  // coefficients given by it->poly

	for( vector<THaMatrixElement>::iterator it=matrix.begin();
			it!=matrix.end(); it++ ) {
				it->v = 0.0;

				if(it->order > 0) {
					for(int i=it->order-1; i>=1; i--)
						it->v = x * (it->v + it->poly[i]);
					it->v += it->poly[0];
				}
			}
}

//_____________________________________________________________________________
Double_t LOpticsOpt::CalcTargetVar(const vector<THaMatrixElement>& matrix,
								   const Double_t powers[][5])
{
  // calculates the value of a variable at the target
  // the x-dependence is already in the matrix, so only 1-3 (or np) used
	Double_t retval=0.0;
	Double_t v=0;
	for( vector<THaMatrixElement>::const_iterator it=matrix.begin();
			it!=matrix.end(); it++ )
		if(it->v != 0.0) {
		v = it->v;
		unsigned int np = it->pw.size(); // generalize for extra matrix elems.
		for (unsigned int i=0; i<np; i++)
			v *= powers[it->pw[i]][i+1];
		retval += v;
  //      retval += it->v * powers[it->pw[0]][1]
  //	              * powers[it->pw[1]][2]
  //	              * powers[it->pw[2]][3];
		}

		return retval;
}

//_____________________________________________________________________________
Double_t LOpticsOpt::CalcTarget2FPLen(const vector<THaMatrixElement>& matrix,
									  const Double_t powers[][5])
{
  // calculates distance from the nominal target position (z=0)
  // to the transport plane

	Double_t retval=0.0;
	for( vector<THaMatrixElement>::const_iterator it=matrix.begin();
			it!=matrix.end(); it++ )
		if(it->v != 0.0)
			retval += it->v * powers[it->pw[0]][0]
					* powers[it->pw[1]][1]
					* powers[it->pw[2]][2]
					* powers[it->pw[3]][3];

	return retval;
}

//_____________________________________________________________________________
UInt_t LOpticsOpt::LoadRawData(TString DataFileName, UInt_t NLoad, UInt_t MaxDataPerGroup)
{
	//load "f51" ascii data file to Rawdata[]

	DEBUG_INFO("LoadRawData","Loading %s",DataFileName.Data());

	if (BeamShiftX!=0)
		DEBUG_WARNING("LoadRawData","Shift Beam X = %f",BeamShiftX);

	UInt_t datagrpcnt[kMaxDataGroup]={0};

	FILE* file = fopen( DataFileName,"r" );
	if( !file ) return kFileError;

	UInt_t NRead = 0;
	// JW: changed line below
	const int LEN = 2000;
	char buff[LEN];

	Double_t NDataRead = 0;
	int NLineRead = 0;

	while( fgets(buff, LEN, file) ) {

		NLineRead++;

		if (NLineRead%100000==0) DEBUG_INFO("LoadRawData","%d/%d Entries Loaded",NRead,NLineRead);

		assert(NRead<MaxNRawData);//too much data if fails

		if (NRead>=NLoad) break;

		Double_t * eventdata = fRawData[NRead].Data;

		string line(buff);
    // Erase trailing newline
		if( line.size() > 0 && line[line.size()-1] == '\n' ) {
			buff[line.size()-1] = 0;
			line.erase(line.size()-1,1);
		}
    // Split the line into whitespace-separated fields
		vector<string> line_spl = Split(line);

		assert(line_spl.size()<=MaxNEventData);//array size check
		for(UInt_t idx = 0; idx<line_spl.size(); idx++)
			eventdata[idx] = atof(line_spl[idx].c_str());

		//WARNING : shift beam x
		if (BeamShiftX!=0)
			eventdata[kBeamX] += BeamShiftX;

		//determine whether to save this data
		// JW: FIX this (just replaced kCutID with KFoilID
		UInt_t cutid = (UInt_t)eventdata[kFoilID];



		// FoilID true/false to exlucde foils not wanted in optimisation

		
		
		

		if( !foil_on[cutid]){

		  cout << "Excluding event from foil " << cutid << endl;
		  continue;
		}



		assert(cutid<kMaxDataGroup); // to many cuts
		UInt_t & grpcnt = datagrpcnt[cutid];
		grpcnt++;
		if (grpcnt>MaxDataPerGroup) {
			DEBUG_MASSINFO("LoadRawData","ignore data %d from cutid %d (%d ev total)",NLineRead,cutid,grpcnt);
			continue;
		}



		UInt_t foilid = cutid;
		UInt_t colid  = eventdata[kColID];
		UInt_t rowid  = eventdata[kRowID];
		UInt_t holeid = Get_Hole(colid,rowid);


		hole_opt_select();
		// ignore prescibed holes 
		if(!hole_select[holeid]){		
		  continue;
		}


		NDataRead+=line_spl.size();

		Double_t (*powers)[5] = fRawData[NRead].powers;
		Double_t x_fp = 	eventdata[kX];
		Double_t th_fp =	eventdata[kTh];
		Double_t y_fp =		eventdata[kY];
		Double_t ph_fp = 	eventdata[kPhi];

  // calculate the powers we need
		for(int i=0; i<kNUM_PRECOMP_POW; i++) {
			powers[i][0] = pow(x_fp, i);
			powers[i][1] = pow(th_fp, i);
			powers[i][2] = pow(y_fp, i);
			powers[i][3] = pow(ph_fp, i);
			powers[i][4] = pow(TMath::Abs(th_fp),i);
		}

		NRead++;
	}

	fclose(file);
	fNRawData = NRead;
	fNCalibData = NRead; //fNCalibData shall be updated later if only part of data read in is for calibration use

	UInt_t goodstatcut=0, actcutcnt=0;
	for(int i=0;i<kMaxDataGroup; i++)
	{
		if (datagrpcnt[i]>0) actcutcnt++;

		if (datagrpcnt[i]>MaxDataPerGroup) goodstatcut++;
	}

	DEBUG_INFO("LoadRawData","Event Limit/Cut = %d, %d / %d ev read, %d / %d cut have enough ev",
			MaxDataPerGroup,NRead,NLineRead,goodstatcut,actcutcnt);
	DEBUG_INFO("LoadRawData","%d events x %f record/event read from %s",
			   fNRawData,NDataRead/fNRawData,DataFileName.Data());
	return NRead;
}

UInt_t LOpticsOpt::SaveDataBuffer(TTree * T)
{
	//save Rawdata[] to T Tree. return N event written

	assert(T);

	DEBUG_INFO("SaveDataBuffer","Saving Data buffer to Tree %s ... ", T->GetName());

	EventData eventdata;

	// JW: Added new cuts definitions/IDs to tree
	//	T->Branch("CutID",&(eventdata.Data[kCutID]),"CutID/D");
	T->Branch("FoilID",&(eventdata.Data[kFoilID]),"FoilID/D");
	T->Branch("ColID",&(eventdata.Data[kColID]),"ColID/D");
	T->Branch("RowID",&(eventdata.Data[kRowID]),"RowID/D");
	T->Branch("X",&(eventdata.Data[kX]),"X/D");
	T->Branch("Th",&(eventdata.Data[kTh]),"Th/D");
	T->Branch("Y",&(eventdata.Data[kY]),"Y/D");
	T->Branch("Phi",&(eventdata.Data[kPhi]),"Phi/D");
	T->Branch("BeamX",&(eventdata.Data[kBeamX]),"BeamX/D");
	T->Branch("BeamY",&(eventdata.Data[kBeamY]),"BeamY/D");
	T->Branch("L_tr_tg_th",&(eventdata.Data[kL_tr_tg_th]),"L_tr_tg_th/D");
	T->Branch("L_tr_tg_ph",&(eventdata.Data[kL_tr_tg_ph]),"L_tr_tg_ph/D");

	T->Branch("RealTh",&(eventdata.Data[kRealTh]),"RealTh/D");
	T->Branch("RealPhi",&(eventdata.Data[kRealPhi]),"RealPhi/D");
	T->Branch("RealTgX",&(eventdata.Data[kRealTgX]),"RealTgX/D");
	T->Branch("RealThMatrix",&(eventdata.Data[kRealThMatrix]),"RealThMatrix/D");
	T->Branch("CalcTh",&(eventdata.Data[kCalcTh]),"CalcTh/D");
	T->Branch("CalcPh",&(eventdata.Data[kCalcPh]),"CalcPh/D");

	T->Branch("L_tr_tg_y",&(eventdata.Data[kL_tr_tg_y]),"L_tr_tg_y/D");
	T->Branch("RealTgY",&(eventdata.Data[kRealTgY]),"RealTgY/D");
	T->Branch("RealReactZ",&(eventdata.Data[kRealReactZ]),"RealReactZ/D");
	T->Branch("CalcTgY",&(eventdata.Data[kCalcTgY]),"CalcTgY/D");
	T->Branch("CalcReactZ",&(eventdata.Data[kCalcReactZ]),"CalcReactZ/D");

	T->Branch("L_tr_tg_dp",&(eventdata.Data[kL_tr_tg_dp]),"L_tr_tg_dp/D");
	T->Branch("L_tr_p",&(eventdata.Data[kL_tr_p]),"L_tr_p/D");
	T->Branch("urb_e",&(eventdata.Data[kurb_e]),"urb_e/D");
	T->Branch("RunNum",&(eventdata.Data[kRunNum]),"RunNum/D");
	T->Branch("ExtraDataFlag",&(eventdata.Data[kExtraDataFlag]),"ExtraDataFlag/D");
	T->Branch("KineID",&(eventdata.Data[kKineID]),"KineID/D");
	T->Branch("Centralp",&(eventdata.Data[kCentralp]),"Centralp/D");
	T->Branch("RadiLossDp",&(eventdata.Data[kRadiLossDp]),"RadiLossDp/D");
	T->Branch("ScatterAngle",&(eventdata.Data[kScatterAngle]),"ScatterAngle/D");
	T->Branch("DpKinOffsets",&(eventdata.Data[kDpKinOffsets]),"DpKinOffsets/D");
	T->Branch("RealDpKin",&(eventdata.Data[kRealDpKin]),"RealDpKin/D");
	T->Branch("RealDpKinMatrix",&(eventdata.Data[kRealDpKinMatrix]),"RealDpKinMatrix/D");
	T->Branch("CalcDpKinMatrix",&(eventdata.Data[kCalcDpKinMatrix]),"CalcDpKinMatrix/D");
	T->Branch("CalcDpKin",&(eventdata.Data[kCalcDpKin]),"CalcDpKin/D");
	T->Branch("RealDpKinExcitations",&(eventdata.Data[kRealDpKinExcitations]),"RealDpKinExcitations/D");

	UInt_t idx=0;
	for ( idx = 0; idx<fNRawData; idx++)
	{
		eventdata = fRawData[idx];
		T->Fill();
	}

	T -> Write();

	DEBUG_INFO("SaveDataBuffer","Done : Saving Data buffer to Tree");
	return idx;
}

UInt_t LOpticsOpt::SaveDataBuffer(TString fname, TString tree)
{
	DEBUG_INFO("SaveDataBuffer","Saving Data buffer to File %s",fname.Data());

	TFile *f = new TFile(fname,"recreate");
	assert(f);
	f->cd();

	TTree *T = new TTree(tree.Data(),"Data buffer of HRS Optics Optimizer");
	assert(T);
	UInt_t cnt = SaveDataBuffer( T );

	f->cd();
	f->Write();

	return cnt;
}

//_____________________________________________________________________________
UInt_t LOpticsOpt::Matrix2Array(Double_t Array[], const std::vector<THaMatrixElement> &Matrix, Bool_t FreeParaFlag[])
{
	  //Matrix -> Array
// 	DEBUG_INFO("Matrix2Array","Entry Point");
	typedef vector<THaMatrixElement>::size_type vsiz_t;

	UInt_t idx =0;

	
	//JW: have altered OptOrder to OptOrder + 1 below


	for (vsiz_t i=0; i<Matrix.size(); i++) {
		const THaMatrixElement& m = Matrix[i];
		UInt_t j;
		//		for (j=0; (int)j<m.order; j++) {
		for (j=0; (int)j<m.OptOrder+1; j++) {
		  if (FreeParaFlag) FreeParaFlag[idx] = j<(m.OptOrder+1)?kTRUE:kFALSE;
			Array[idx++]=m.poly[j];
		}
		for (; j<kPORDER; j++) {
		  if (FreeParaFlag) FreeParaFlag[idx] = j<(m.OptOrder+1)?kTRUE:kFALSE;
			Array[idx++]=0;
		}

		cout << "for matrix element " << i << " Optorder = " << m.OptOrder << " and kPORDER = " << kPORDER << endl;
	}

	DEBUG_INFO("Matrix2Array","Fill Size = %d",idx);

	return idx;
}

//_____________________________________________________________________________
UInt_t LOpticsOpt::Array2Matrix(const Double_t Array[], std::vector<THaMatrixElement> &Matrix)
{
	//Array -> fCurrentMatrixElems
// 	DEBUG_INFO("Array2Matrix","Entry Point");
	typedef vector<THaMatrixElement>::size_type vsiz_t;

	UInt_t idx =0;
	Double_t maxele = 0, sumele=0;
	vsiz_t 	max_i=-1;
	int 	max_j=-1;


	//	cout << "Matrix size is " << Matrix.size() << " and max_i = " << Int_t(max_i) << endl;
	for (vsiz_t i=0; i<Matrix.size(); i++) {
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

     
		// for (j=0; j<m.order; j++) {
		// 	int totalpow = j + m.pw[0] + m.pw[1] + m.pw[2];			
		// 	if (totalpow <1) totalpow = 1;
		// 	const Double_t effsize = TMath::Power(TMath::Abs(Array[idx]), 1./totalpow);

		// 	if (effsize>maxele) {max_i=i; max_j=j;maxele=effsize;}
			
		// 	if(max_j<0){

		// 	  cout << " max_j < 0 for " << endl;
		// 	  for (vsiz_t j=0; j<m.pw.size(); j++) {
		// 	    cout << Form("%d ",m.pw[j]);
		// 	  }

		// 	  cout << " where effsize = " << effsize << " and Array[idx] = " << Array[idx] << " and totalpow = " << totalpow << endl;
		// 	}
		// 	else{
		// 	  cout << " max_j > 0 for " << endl;			  
		// 	  			  for (vsiz_t j=0; j<m.pw.size(); j++) {
		// 	    cout << Form("%d ",m.pw[j]);
		// 	  }

		// 	  cout << " where effsize = " << effsize << " and Array[idx] = " << Array[idx] << " and totalpow = " << totalpow << endl;
		// 	}




		// 	//			cout << "i = " << i << " & max_i = " << Int_t(max_i) << " & effsize = " << effsize << endl;
		// 	sumele += effsize;// TMath::Abs(Array[idx]);
		// 	m.poly[j]=Array[idx];
		// 	idx++;
		// }


		m.SkimPoly();
	}

	// std::cout << "max_i =" << int(max_i) << ", & Matrix.size() = " << int(Matrix.size()) << std::endl;
	// cout << "int(max_j) = " << max_j << endl;

	
	assert(int(max_i)<int(Matrix.size()));

	assert( max_j >=0);
	THaMatrixElement& m = Matrix[max_i];
	//	DEBUG_INFO("Array2Matrix","Load Size = %d, max ele = (%d %d %d %d) = %f (%f), Eff average = %f"
	//			,idx, max_j, m.pw[0], m.pw[1], m.pw[2], m.poly[max_j], maxele, sumele/idx);

	return idx;
}

//_____________________________________________________________________________
void LOpticsOpt::Print(const Option_t* opt) const
{
	//Print current matrix

	THaTrackingDetector::Print(opt);
	typedef vector<THaMatrixElement>::size_type vsiz_t;

  // Print out the optics matrices, to verify they make sense
// 	printf("Matrix FP (t000, y000, p000)\n");
// 	for (vsiz_t i=0; i<fFPMatrixElems.size(); i++) {
// 		const THaMatrixElement& m = fFPMatrixElems[i];
// 		for (vsiz_t j=0; j<m.pw.size(); j++) {
// 			printf("  %2d",m.pw[j]);
// 		}
// 		for (int j=0; j<m.order; j++) {
// 			printf("  %g",m.poly[j]);
// 		}
// 		printf(" : Opt -> %d",m.OptOrder);
// 		printf("\n");
// 	}

	printf("LOpticsOpt::Print: Transport Matrix:  D-terms\n");
	for (vsiz_t i=0; i<fDMatrixElems.size(); i++) {
		const THaMatrixElement& m = fDMatrixElems[i];
		for (vsiz_t j=0; j<m.pw.size(); j++) {
			printf("  %2d",m.pw[j]);
		}
		for (int j=0; j<m.order; j++) {
			printf("\t%g",m.poly[j]);
		}
		printf(" : Opt -> %d",m.OptOrder);
		if ((UInt_t)m.order!=m.OptOrder) printf(" != Matrix Order !!");
		printf("\n");
	}

	printf("LOpticsOpt::Print: Transport Matrix:  T-terms\n");
	for (vsiz_t i=0; i<fTMatrixElems.size(); i++) {
		const THaMatrixElement& m = fTMatrixElems[i];
		for (vsiz_t j=0; j<m.pw.size(); j++) {
			printf("  %2d",m.pw[j]);
		}
		for (int j=0; j<m.order; j++) {
			printf("\t%g",m.poly[j]);
		}
		printf(" : Opt -> %d",m.OptOrder);
		if ((UInt_t)m.order!=m.OptOrder) printf(" != Matrix Order !!");
		printf("\n");
	}

	printf("LOpticsOpt::Print: Transport Matrix:  Y-terms\n");
	for (vsiz_t i=0; i<fYMatrixElems.size(); i++) {
		const THaMatrixElement& m = fYMatrixElems[i];
		for (vsiz_t j=0; j<m.pw.size(); j++) {
			printf("  %2d",m.pw[j]);
		}
		for (int j=0; j<m.order; j++) {
			printf("\t%g",m.poly[j]);
		}
		printf(" : Opt -> %d",m.OptOrder);
		if ((UInt_t)m.order!=m.OptOrder) printf(" != Matrix Order !!");
		printf("\n");
	}

// 	printf("Transport Matrix:  YTA-terms (abs(theta))\n");
// 	for (vsiz_t i=0; i<fYTAMatrixElems.size(); i++) {
// 		const THaMatrixElement& m = fYTAMatrixElems[i];
// 		for (vsiz_t j=0; j<m.pw.size(); j++) {
// 			printf("  %2d",m.pw[j]);
// 		}
// 		for (int j=0; j<m.order; j++) {
// 			printf("\t%g",m.poly[j]);
// 		}
// 		printf(" : Opt -> %d",m.OptOrder);
// 		printf("\n");
// 	}

	printf("LOpticsOpt::Print: Transport Matrix:  P-terms\n");
	for (vsiz_t i=0; i<fPMatrixElems.size(); i++) {
		const THaMatrixElement& m = fPMatrixElems[i];
		for (vsiz_t j=0; j<m.pw.size(); j++) {
			printf("  %2d",m.pw[j]);
		}
		for (int j=0; j<m.order; j++) {
			printf("\t%g",m.poly[j]);
		}
		printf(" : Opt -> %d",m.OptOrder);
		if ((UInt_t)m.order!=m.OptOrder) printf(" != Matrix Order !!");
		printf("\n");
	}

// 	printf("Transport Matrix:  PTA-terms\n");
// 	for (vsiz_t i=0; i<fPTAMatrixElems.size(); i++) {
// 		const THaMatrixElement& m = fPTAMatrixElems[i];
// 		for (vsiz_t j=0; j<m.pw.size(); j++) {
// 			printf("  %2d",m.pw[j]);
// 		}
// 		for (int j=0; j<m.order; j++) {
// 			printf("\t%g",m.poly[j]);
// 		}
// 		printf(" : Opt -> %d",m.OptOrder);
// 		printf("\n");
// 	}

// 	printf("Matrix L\n");
// 	for (vsiz_t i=0; i<fLMatrixElems.size(); i++) {
// 		const THaMatrixElement& m = fLMatrixElems[i];
// 		for (vsiz_t j=0; j<m.pw.size(); j++) {
// 			printf("  %2d",m.pw[j]);
// 		}
// 		for (int j=0; j<m.order; j++) {
// 			printf("\t%g",m.poly[j]);
// 		}
// // 		printf(" : Opt -> %d",m.OptOrder);
// 		printf("\n");
// 	}

	printf("fArbitaryVertexShift[%d] = {",NFiols);
	for(UInt_t FoilID = 0; FoilID<NFiols; FoilID++)
		printf("%f  ",fArbitaryVertexShift[FoilID]);
	printf("}\n");
	printf("fArbitaryDpKinShift[%d] = {",NKine);
	for(UInt_t KineID=0; KineID<NKine; KineID++)
		printf("%f  ",fArbitaryDpKinShift[KineID]);
	printf("}\n");

	return;
}

//_____________________________________________________________________________
Int_t LOpticsOpt::SaveDataBase(TString DataBaseName)
{

	//output database in memory to new database file
	//WARNING: Hard coded text included

  
        // for saving canvasses later 
        fDestDataBaseName = DataBaseName;

	DEBUG_INFO("SaveDataBase","Saving to %s",DataBaseName.Data());

	typedef vector<THaMatrixElement>::size_type vsiz_t;

	FILE* file = fopen( DataBaseName,"w" );
	if( !file ) {
		Info("SaveDataBase","Error Openin %s",DataBaseName.Data());
		return kFileError;
	}
	TDatime dt;
	fprintf(file,"# -------------------------------------------------------------");	fprintf(file,"\n");
	fprintf(file,"# Optimized by John Williamson @ %s",dt.AsString());fprintf(file,"\n");
	fprintf(file,"# Saved to %s",DataBaseName.Data());fprintf(file,"\n");
	fprintf(file,OldComments);

	//Header Part
	// 	[ L.global ]
	// 	0.3327 1 0.0 270.2 0.0 -1.6e-03        VDC Angle, Plane Spacing, Gamma Coefficents
	// 	matrix elements
	// 	t 0 0 0  -1.001135e+00 -3.313373e-01 -4.290819e-02  4.470852e-03  0.000000e+00  0.000000e+00  0.000000e+00  0
	// 			y 0 0 0  -8.060915e-03  1.071977e-03  9.019102e-04 -3.239615e-04  0.000000e+00  0.000000e+00  0.000000e+00  0
	// 			p 0 0 0  -2.861912e-03 -2.469069e-03  8.427172e-03  2.274635e-03  0.000000e+00  0.000000e+00  0.000000e+00  0
//	fprintf(file,"[ L.global ]");fprintf(file,"\n");
//	fprintf(file,"0.3327 1 0.0 270.2 0.0 -1.6e-03        VDC Angle, Plane Spacing, Gamma Coefficents");fprintf(file,"\n");
//	fprintf(file,"matrix elements");fprintf(file,"\n");
//	fprintf(file,"t 0 0 0  -1.001135e+00 -3.313373e-01 -4.290819e-02  4.470852e-03  0.000000e+00  0.000000e+00  0.000000e+00  0");fprintf(file,"\n");
//	fprintf(file,"y 0 0 0  -8.060915e-03  1.071977e-03  9.019102e-04 -3.239615e-04  0.000000e+00  0.000000e+00  0.000000e+00  0");fprintf(file,"\n");
//	fprintf(file,"p 0 0 0  -2.861912e-03 -2.469069e-03  8.427172e-03  2.274635e-03  0.000000e+00  0.000000e+00  0.000000e+00  0");fprintf(file,"\n");

	fprintf(file,DatabaseHeader);

	DEBUG_INFO("SaveDataBase","Transport Matrix:  D-terms");
	for (vsiz_t i=0; i<fDMatrixElems.size(); i++) {
		fprintf(file,"D ");
		const THaMatrixElement& m = fDMatrixElems[i];
		for (vsiz_t j=0; j<m.pw.size(); j++) {
			fprintf(file,"%d ",m.pw[j]);
		}
		int j;
		for (j=0; j<m.order; j++) {
			fprintf(file," %13.6e",m.poly[j]);
		}
		for (; j<kPORDER; j++) {
			fprintf(file," %13.6e",0.0);
		}
		fprintf(file,"  %d",m.OptOrder);
		fprintf(file,"\n");
	}

	DEBUG_INFO("SaveDataBase","Transport Matrix:  T-terms");
	for (vsiz_t i=0; i<fTMatrixElems.size(); i++) {
		fprintf(file,"T ");
		const THaMatrixElement& m = fTMatrixElems[i];
		for (vsiz_t j=0; j<m.pw.size(); j++) {
			fprintf(file,"%d ",m.pw[j]);
		}
		int j;
		for (j=0; j<m.order; j++) {
			fprintf(file," %13.6e",m.poly[j]);
		}
		for (; j<kPORDER; j++) {
			fprintf(file," %13.6e",0.0);
		}
		fprintf(file,"  %d",m.OptOrder);
		fprintf(file,"\n");
	}

	DEBUG_INFO("SaveDataBase","Transport Matrix:  P-terms");
	for (vsiz_t i=0; i<fPMatrixElems.size(); i++) {
		fprintf(file,"P ");
		const THaMatrixElement& m = fPMatrixElems[i];
		for (vsiz_t j=0; j<m.pw.size(); j++) {
			fprintf(file,"%d ",m.pw[j]);
		}
		int j;
		//JW changed
		// j<m.order to j<(m.OptOrder+1)
		for (j=0; j<(m.OptOrder+1); j++) {
			fprintf(file," %13.6e",m.poly[j]);
		}
		for (; j<kPORDER; j++) {
			fprintf(file," %13.6e",0.0);
		}
		fprintf(file,"  %d",m.OptOrder);
		fprintf(file,"\n");
	}

	DEBUG_INFO("SaveDataBase","Transport Matrix:  Y-terms");
	for (vsiz_t i=0; i<fYMatrixElems.size(); i++) {
		fprintf(file,"Y ");
		const THaMatrixElement& m = fYMatrixElems[i];
		for (vsiz_t j=0; j<m.pw.size(); j++) {
			fprintf(file,"%d ",m.pw[j]);
		}
		int j;
		for (j=0; j<m.order; j++) {
			fprintf(file," %13.6e",m.poly[j]);
		}
		for (; j<kPORDER; j++) {
			fprintf(file," %13.6e",0.0);
		}
		fprintf(file,"  %d",m.OptOrder);
		fprintf(file,"\n");
	}

	// L and XF Matrix
//	fprintf(file,"L 0 0 0 0  25.713");fprintf(file,"\n");
//	fprintf(file,"L 1 0 0 0  0.1650 ");fprintf(file,"\n");
//	fprintf(file,"L 2 0 0 0 -0.05");fprintf(file,"\n");
//	fprintf(file,"L 0 1 0 0 -11.6554");fprintf(file,"\n");
//	fprintf(file,"L 0 2 0 0 -9.4951");fprintf(file,"\n");
//	fprintf(file,"L 0 0 1 0  0.0");fprintf(file,"\n");
//	fprintf(file,"L 0 0 2 0  0.0");fprintf(file,"\n");
//	fprintf(file,"L 0 0 0 1  0.0");fprintf(file,"\n");
//	fprintf(file,"L 0 0 0 2  0.0");fprintf(file,"\n");
//	fprintf(file,"XF 1 0 0 0 0 -2.181E+00");fprintf(file,"\n");
//	fprintf(file,"XF 0 1 0 0 0 -1.980E-01");fprintf(file,"\n");
//	fprintf(file,"XF 0 0 0 0 1  1.191E+01");fprintf(file,"\n");
//	fprintf(file,"TF 1 0 0 0 0 -1.000E-01");fprintf(file,"\n");
//	fprintf(file,"TF 0 1 0 0 0 -4.690E-01");fprintf(file,"\n");
//	fprintf(file,"TF 0 0 0 0 1  1.967E+00");fprintf(file,"\n");
//	fprintf(file,"PF 0 0 1 0 0  3.630E-01");fprintf(file,"\n");
//	fprintf(file,"PF 0 0 0 1 0 -0.902E+00");fprintf(file,"\n");
//	fprintf(file,"YF 0 0 1 0 0 -5.950E-01");fprintf(file,"\n");
//	fprintf(file,"YF 0 0 0 1 0 -1.274E+00");fprintf(file,"\n");

	fprintf(file,DatabaseFooter);

	fclose(file);

	return kOK;
}

//_____________________________________________________________________________
bool THaMatrixElement::match(const THaMatrixElement& rhs) const
{
  // Compare coefficients of this matrix element to another

	if( pw.size() != rhs.pw.size() )
		return false;
	for( vector<int>::size_type i=0; i<pw.size(); i++ ) {
		if( pw[i] != rhs.pw[i] )
			return false;
	}
	return true;
}
//_____________________________________________________________________________
void THaMatrixElement::SkimPoly()
{
  //reduce order to highest non-zero poly

	if (iszero) return;

	while(!poly[order-1] && order >0)
	{
		poly.pop_back();
		order = order-1;
	}

	if (order==0) iszero = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
//_____________________________________________________________________________
THaAnalysisObject::EStatus LOpticsOpt::Init( const TDatime& /*date*/ )
{
  // Initialize VDC. Calls standard Init(), then initializes subdetectors.


	return fStatus = kOK;
}

//_____________________________________________________________________________
Int_t LOpticsOpt::ConstructTracks( TClonesArray* /*tracks*/, Int_t /*mode*/ )
{

	return 0;
}

//_____________________________________________________________________________
void LOpticsOpt::Clear( Option_t* /*opt*/ )
{
}

//_____________________________________________________________________________
Int_t LOpticsOpt::Decode( const THaEvData& /*evdata*/ )
{
	return 0;
}

//_____________________________________________________________________________
Int_t LOpticsOpt::CoarseTrack( TClonesArray& /*tracks*/ )
{
	return 0;
}

//_____________________________________________________________________________
Int_t LOpticsOpt::FineTrack( TClonesArray& /*tracks*/ )
{

	return 0;
}

//_____________________________________________________________________________
Int_t LOpticsOpt::FindVertices( TClonesArray& tracks )
{
  // Calculate the target location and momentum at the target.
  // Assumes that CoarseTrack() and FineTrack() have both been called.

	Int_t n_exist = tracks.GetLast()+1;
	for( Int_t t = 0; t < n_exist; t++ ) {
		THaTrack* theTrack = static_cast<THaTrack*>( tracks.At(t) );
		CalcTargetCoords(theTrack, kRotatingTransport);
	}


	return 0;
}

ClassImp(LOpticsOpt);
