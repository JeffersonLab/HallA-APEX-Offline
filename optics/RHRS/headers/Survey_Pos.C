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


#include "Survey_Pos.h"
#include "InputR.h"



void Initialize(){


  TVector3 TCSX(0, -1, 0);
  TVector3 TCSZ(TMath::Sin(HRSAngle), 0, TMath::Cos(HRSAngle));
  TVector3 TCSY = TCSZ.Cross(TCSX);
  fTCSInHCS.RotateAxes(TCSX, TCSY, TCSZ);
  
  const Int_t a = (HRSAngle > 0) ? 1 : -1;
  //fPointingOffset, HCS
  fPointingOffset.SetXYZ(-a*MissPointZ*TMath::Cos(HRSAngle), MissPointY, -MissPointZ * TMath::Sin(HRSAngle)); 

}


void Hole_exp(int n_foil,int &row_min,int &row_max,int &col_min,int &col_max){

  if(n_foil == 1){
      row_min = 3;
      row_max = 13;
      col_min = 4;
      col_max = 19;
    }
    
    if(n_foil == 0){
      row_min = 2;
      row_max = 14;
      col_min = 6;
      col_max = 24;
    }
    
    if(n_foil == 2){
      row_min = 4;
      row_max = 12;
      col_min = 2;
      col_max = 11;
    }

    if(n_foil == 7){
      row_min = 5;
      row_max = 13;
      col_min = 8;
      col_max = 22;
    }
    
    if(n_foil == 8){
      row_min = 3;
      row_max = 13;
      col_min = 6;
      col_max = 20;
    }
    
    if(n_foil == 9){
      row_min = 4;
      row_max = 12;
      col_min = 3;
      col_max = 11;
    }
    
    if(n_foil == 10){
      row_min = 5;
      row_max = 12;
      col_min = 0;
      col_max = 3;
    }

    if(n_foil == 3){
      row_min = 3;
      row_max = 10;
      col_min = 10;
      col_max = 23;
    }
    
    if(n_foil == 4){
      row_min = 3;
      row_max = 13;
      col_min = 6;
      col_max = 22;
    }
    
    if(n_foil == 5){
      row_min = 4;
      row_max = 13;
      col_min = 3;
      col_max = 16;
    }

    if(n_foil == 6){
      row_min = 5;
      row_max = 12;
      col_min = 1;
      col_max = 10;
    }




}




TVector3 BeamSpotHCS_Correction(UInt_t FoilID, double beam_y, double beam_z){
  
  TVector3 foil(targetfoilsX[FoilID], beam_y, beam_z);

  foil.RotateY(target_yaw);

  return foil;
}

const TVector3 GetSieveHoleTCS(UInt_t Col, UInt_t Row)
{
    /// Calculate position with survey info ////

    TVector3 SieveHoleTCS(SieveXbyRow[Row], SieveYbyCol[Col], 0);
    
    SieveHoleTCS.RotateX(-(yaw - HRSAngle));
    SieveHoleTCS.RotateY(pitch - TMath::Pi()/2);

    SieveHoleTCS.SetXYZ( SieveHoleTCS.X() + SieveOffX,  SieveHoleTCS.Y() + SieveOffY, SieveHoleTCS.Z() + ZPos + SieveOffZ);
   
 
    return SieveHoleTCS;
}


const TVector3 GetSieveHoleCorrectionTCS(UInt_t nfoil, UInt_t Col, UInt_t Row)
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

    //const TVector3 BeamSpotHCS_average(BeamX_average[nfoil], BeamY_average[nfoil], targetfoils[nfoil]);
    const TVector3 BeamSpotHCS_average = BeamSpotHCS_Correction(nfoil,BeamY_average[nfoil], targetfoils[nfoil]);
    const TVector3 BeamSpotTCS_average = fTCSInHCS.Inverse()*(BeamSpotHCS_average - fPointingOffset);        
    BeamXHCS= BeamSpotTCS_average.X();
    BeamYHCS= BeamSpotTCS_average.Y();

    const TVector3 SieveHoleTCS = GetSieveHoleTCS(Col, Row);

    Z_distance = SieveHoleTCS.Z() - BeamSpotTCS_average.Z();
    Y_p = SieveHoleTCS.Y() +  SieveRadius/2.* 25.4e-3; // .157/2. * 25.4e-3 is the sieve hole radius
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




void Sieve_hole_pos(int FoilID, int Col, int Row, double sieve_ph_th[], double  sieve_yx[]){


  //Initialize();
  

  const TVector3 SieveHoleCorrectionTCS = GetSieveHoleCorrectionTCS(FoilID, Col, Row);

  const TVector3 BeamSpotHCS(BeamX_average[FoilID], BeamY_average[FoilID], targetfoils[FoilID]);
  
  //const TVector3 BeamSpotHCS = BeamSpotHCS_Correction(FoilID, BeamY_average[FoilID], targetfoils[FoilID]);
  
  const TVector3 BeamSpotTCS = fTCSInHCS.Inverse()*(BeamSpotHCS - fPointingOffset);        
  
  const TVector3 MomDirectionTCS = SieveHoleCorrectionTCS - BeamSpotTCS;
  
  //cout<<Col<<endl;
  
  
  sieve_ph_th[1] = MomDirectionTCS.X() / MomDirectionTCS.Z();
  sieve_ph_th[0] = MomDirectionTCS.Y() / MomDirectionTCS.Z();
  
  sieve_yx[1] = BeamSpotTCS.X() - BeamSpotTCS.Z() * sieve_ph_th[1];
  sieve_yx[0] = BeamSpotTCS.Y() - BeamSpotTCS.Z() * sieve_ph_th[0];

   
}
