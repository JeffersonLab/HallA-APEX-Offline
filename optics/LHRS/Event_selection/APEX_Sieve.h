/////////////////////////////////////////////
/////////////////////////////////////////////
//
//
// APEX_Sieve
//
// small class used to get Sieve information:
// translates between hole number to column and row and vice vers
// (slighlty more complex arangement in APEX sieve)
// 
// Also has function to return TCS (Target Co-ordinate System) coords
//
/////////////////////////////////////////////



#ifndef ROOT_APEX_Sieve
#define ROOT_APEX_Sieve


#include "file_def.h"
#include "../opt_new/InputAPEXL.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TCut.h"


std::vector<int> Get_Col_Row(Int_t Hole);




const TVector3 GetSieveHoleTCS(Int_t Col, Int_t Row) /*const*/
{


  if(Col < 0 || Col > NSieveCol){
    
    Col = 1;
  }


  if(Row < 0 || Row > NSieveRow){
    
    Row = 1;
  }



  // following lines get Sieve offset in HCS (Hall Co-ordinate System) and transfer to TCS (Target Co-ordinate System)

  TVector3 TCSX(0, -1, 0);
  TVector3 TCSZ(TMath::Sin(HRSAngle), 0, TMath::Cos(HRSAngle));
  TVector3 TCSY = TCSZ.Cross(TCSX);

  TRotation fTCSInHCS;
  fTCSInHCS.RotateAxes(TCSX, TCSY, TCSZ);

  TVector3 Sieve_offset_HCS(SieveOffX_HCS,SieveOffY_HCS,ZPos_HCS);

  TVector3 Sieve_offset_TCS = fTCSInHCS.Inverse()*(Sieve_offset_HCS);
    

  //  TVector3 SieveHoleTCS(SieveXbyRow[Row]+SieveOffX, SieveYbyCol[Col]+SieveOffY, ZPos);


  TVector3 SieveHoleTCS( SieveXbyRow[Row] + Sieve_offset_TCS.X(), Sieve_offset_TCS.Y() + SieveYbyCol[Col], Sieve_offset_TCS.Z());


  return SieveHoleTCS;



}

// //_____________________________________________________________________________



const TVector3 GetSieveHoleTCS(UInt_t Hole) /*const*/
{

  std::vector<int> x_y = {};
  x_y = Get_Col_Row( Hole);

  UInt_t Col = x_y[0];
  UInt_t Row = x_y[1];

  
  TVector3 SieveHoleTCS = {};

  SieveHoleTCS =  GetSieveHoleTCS(Col, Row);

  return SieveHoleTCS;



}


// corrected sieve hole functions


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

    //const TVector3 BeamSpotHCS_average(BeamX_average, BeamY_average, targetfoils[nfoil]);
    //    const TVector3 BeamSpotHCS_average(BeamX_average[nfoil], BeamY_average, targetfoils[nfoil]);


    TVector3 TCSX(0, -1, 0);
    TVector3 TCSZ(TMath::Sin(HRSAngle), 0, TMath::Cos(HRSAngle));
    TVector3 TCSY = TCSZ.Cross(TCSX);

    TRotation fTCSInHCS;
    fTCSInHCS.RotateAxes(TCSX, TCSY, TCSZ);



    const TVector3 BeamSpotHCS_average(BeamX_average, BeamY_average, targetfoils[nfoil]);

    const Int_t a = (HRSAngle > 0) ? 1 : -1;
    //    fPointingOffset.SetXYZ(a*-MissPointZ*TMath::Sin(HRSAngle)*TMath::Cos(HRSAngle),(Double_t)MissPointY,MissPointZ*TMath::Sin(HRSAngle)*TMath::Sin(HRSAngle));
    fPointingOffset.SetXYZ(-a*MissPointZ*TMath::Cos(HRSAngle), MissPointY, -MissPointZ * TMath::Sin(HRSAngle));


    const TVector3 BeamSpotTCS_average = fTCSInHCS.Inverse()*(BeamSpotHCS_average - fPointingOffset);        
    BeamXHCS= BeamSpotTCS_average.X();
    BeamYHCS= BeamSpotTCS_average.Y();

    const TVector3 SieveHoleTCS = GetSieveHoleTCS(Col, Row);

    Z_distance = SieveHoleTCS.Z() - BeamSpotTCS_average.Z();
    Y_p = SieveHoleTCS.Y() +  .157/2.* 25.4e-3; //  // .157/2. * 25.4e-3 is the sieve hole radius
    Y_m = SieveHoleTCS.Y() - .157/2. * 25.4e-3;
    Yback_p = Z_distance /(Z_distance + 25.4e-3) * (Y_p-BeamYHCS) + BeamYHCS;
    Yback_m = Z_distance /(Z_distance + 25.4e-3) * (Y_m-BeamYHCS) + BeamYHCS;
    Yreal_p = (Y_p >= Yback_p) ? Yback_p : Y_p;
    Yreal_m = (Y_m >= Yback_m) ? Y_m : Yback_m;
    X_p = SieveHoleTCS.X() +  .157/2. * 25.4e-3; 
    X_m = SieveHoleTCS.X() -  .157/2. * 25.4e-3;
    Xback_p = Z_distance /(Z_distance + 25.4e-3) * (X_p-BeamXHCS) + BeamXHCS;
    Xback_m = Z_distance /(Z_distance + 25.4e-3) * (X_m-BeamXHCS) + BeamXHCS;
    Xreal_p = (X_p >= Xback_p) ? Xback_p : X_p;
    Xreal_m = (X_m >= Xback_m) ? X_m : Xback_m;

    const TVector3 SieveHoleCorrectionTCS((Xreal_p + Xreal_m)/2, (Yreal_p + Yreal_m)/2 ,SieveHoleTCS.Z());
    //const TVector3 SieveHoleCorrectionTCS( SieveHoleTCS.X(), SieveHoleTCS.Y() ,SieveHoleTCS.Z() );
    return SieveHoleCorrectionTCS;

   
}



const TVector3 GetSieveHoleCorrectionTCS(UInt_t nfoil, UInt_t Hole) /*const*/
{

  std::vector<int> x_y = {};
  x_y = Get_Col_Row( Hole);

  UInt_t Col = x_y[0];
  UInt_t Row = x_y[1];

  
  TVector3 SieveHoleTCS = {};

  SieveHoleTCS =  GetSieveHoleCorrectionTCS(nfoil, Col, Row);

  return SieveHoleTCS;



}



// //_____________________________________________________________________________


std::vector<int> Get_Col_Row(Int_t Hole){


  Int_t row_comp = 0;
  Int_t no_col = 0;
  Int_t col = 0;
  Int_t row = 0;


  for(Int_t i = 0; i<NSieveRow; i++){

    row_comp += NoinEachRow[i];

    if( (row_comp-1) >= Hole){
      row = i;
      no_col = Hole - ( row_comp - NoinEachRow[i]);
      
      break;
    }

  }



  if(row%2 == 0){
    if(no_col==13){
      col = 25;
    }
    else if (no_col==14){
      col = 26;
    }
    else{

    col = no_col *2;
    }
  }

  if(row%2 == 1){

    if(row > 1 && row < 15){
      if(no_col >5){
	col = (no_col*2)+3;
      }
      else{
	col = (no_col*2) +1;
      }
    }
    else{
      col = (no_col*2) +1;
    }

    
  }

  std::vector<int> rowcol{col, row};

  return rowcol;
}




// //____________________________________________________________________________



Int_t Get_Hole(Int_t Col, Int_t Row){


  Int_t hole_no = 0;

  for(Int_t i = 0; i<Row; i++){

    hole_no += NoinEachRow[i];

  }

  // else if conditions here deal with columns at right edge of sieve slit (area odd rows do not have holes)

  if(Row%2 == 0){
    if(Col < 25){
      hole_no += (Col/2);
    }
    else if(Col == 25){
      hole_no += 13;      
    }
    else if(Col == 26){
      hole_no += 14;
    }      
  }


  if(Row%2 == 1){
    if(Row > 1 && Row < 15){
      if(Col < 13){
	Col = ((Col+1)/2) - 1;
      }
      else if (Col >= 13){
	Col = ((Col-3)/2);
      }      
    }
    else{
      Col = ((Col+1)/2) -1;
    }    
    hole_no += Col;
  }
  

  //  hole_no += Col;



  return hole_no;
}








#endif

