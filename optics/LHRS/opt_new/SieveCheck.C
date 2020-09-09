//#include "file_def.h"
#include "InputAPEXL.h"
#include <algorithm>
#include "TVector3.h"
//#include "ROpticsOpt.h"
#include  "APEX_Sieve.h"


// const TVector3 GetSieveHoleTCS(UInt_t Col, UInt_t Row);

// const TVector3 GetSieveHoleCorrectionTCS(UInt_t nfoil, UInt_t Col, UInt_t Row);


TRotation fTCSInHCS;
//TVector3 fPointingOffset;

void SieveCheck(Int_t FoilID){

  // const Double_t BeamX_average = -0.0006391;
  // const Double_t BeamY_average = 0.002405; 
 
  
  //  TRotation fTCSInHCS;
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


  
  
  
  // for(vector<Double_t>::size_type i = 0; i != Row_vals.size(); i++){
  //   cout << Row_vals[i] << " ";
  // }


  
}







// const TVector3 GetSieveHoleTCS(UInt_t Col, UInt_t Row)
// {
//     assert(Col < NSieveCol);
//     assert(Row < NSieveRow);

//     Double_t XbyRow =0;
//     /*
//     if((Col%2) == 0){
//       XbyRow = SieveXbyRowEven[Row];
//     }
//     if((Col%2) == 1){
//       XbyRow = SieveXbyRowOdd[Row];
//     }
//     */

//     XbyRow = SieveXbyRow[Row];

//     TVector3 SieveHoleTCS(SieveOffX + XbyRow, SieveOffY + SieveYbyCol[Col], ZPos);
//     /*
//     cout<<"Col%2:"<<Col%2<<endl;
//     cout<<"Col:"<<Col<<endl;
//     cout<<"XbyRow:"<<XbyRow<<endl;
//     cout<<"Row:"<<Row<<endl;
//     cout<<"YbyCol:"<<SieveYbyCol[Col]<<endl;
//     cout<<"*******************"<<endl;
//     */
//     return SieveHoleTCS;
// }



// const TVector3 GetSieveHoleCorrectionTCS(UInt_t nfoil, UInt_t Col, UInt_t Row){
//     assert(nfoil <NFoils);
//     assert(Col < NSieveCol);
//     assert(Row < NSieveRow);
     
//     //consider the difference of real distribution and hole center
//     Double_t Y_p=0, Y_m=0, Yback_p=0, Yback_m=0; // Y the sieve hole Y position
//     Double_t X_p=0, X_m=0, Xback_p=0, Xback_m=0; // X the sieve hole X position
//     Double_t Yreal_p=0, Yreal_m=0, Xreal_p=0, Xreal_m=0; // real distrituion limit
//     Double_t Z_distance=0;
//     Double_t BeamYHCS=0,BeamXHCS=0;
//     //    Double_t SieveY_Correction[NFoils][NSieveCol][NSieveRow] ={{{0}}};
//     //    Double_t SieveX_Correction[NFoils][NSieveCol][NSieveRow] ={{{0}}};

//     //const TVector3 BeamSpotHCS_average(BeamX_average, BeamY_average, targetfoils[nfoil]);
//     //    const TVector3 BeamSpotHCS_average(BeamX_average[nfoil], BeamY_average, targetfoils[nfoil]);
//     const TVector3 BeamSpotHCS_average(BeamX_average[nfoil], BeamY_average[nfoil], targetfoils[nfoil]);
//     const TVector3 BeamSpotTCS_average = fTCSInHCS.Inverse()*(BeamSpotHCS_average - fPointingOffset);        
//     BeamXHCS= BeamSpotTCS_average.X();
//     BeamYHCS= BeamSpotTCS_average.Y();

//     const TVector3 SieveHoleTCS = GetSieveHoleTCS(Col, Row);

//     Z_distance = SieveHoleTCS.Z() - BeamSpotTCS_average.Z();
//     Y_p = SieveHoleTCS.Y() +  .157/2.* 25.4e-3; //  // .157/2. * 25.4e-3 is the sieve hole radius
//     Y_m = SieveHoleTCS.Y() - .157/2. * 25.4e-3;
//     Yback_p = Z_distance /(Z_distance + 25.4e-3) * (Y_p-BeamYHCS) + BeamYHCS;
//     Yback_m = Z_distance /(Z_distance + 25.4e-3) * (Y_m-BeamYHCS) + BeamYHCS;
//     Yreal_p = (Y_p >= Yback_p) ? Yback_p : Y_p;
//     Yreal_m = (Y_m >= Yback_m) ? Y_m : Yback_m;
//     X_p = SieveHoleTCS.X() +  .157/2. * 25.4e-3; 
//     X_m = SieveHoleTCS.X() -  .157/2. * 25.4e-3;
//     Xback_p = Z_distance /(Z_distance + 25.4e-3) * (X_p-BeamXHCS) + BeamXHCS;
//     Xback_m = Z_distance /(Z_distance + 25.4e-3) * (X_m-BeamXHCS) + BeamXHCS;
//     Xreal_p = (X_p >= Xback_p) ? Xback_p : X_p;
//     Xreal_m = (X_m >= Xback_m) ? X_m : Xback_m;

//     const TVector3 SieveHoleCorrectionTCS((Xreal_p + Xreal_m)/2, (Yreal_p + Yreal_m)/2 ,SieveHoleTCS.Z());
//     //const TVector3 SieveHoleCorrectionTCS( SieveHoleTCS.X(), SieveHoleTCS.Y() ,SieveHoleTCS.Z() );
//     return SieveHoleCorrectionTCS;

   
// }
