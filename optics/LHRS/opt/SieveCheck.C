#include "file_def.h"
#include "InputAPEXL.h"
#include  "APEX_Sieve.h"





void SieveCheck(Int_t FoilID){

  const Double_t BeamX_average = -0.0006391;
  const Double_t BeamY_average = 0.002405; 
 
  
  TRotation fTCSInHCS;
  TVector3 TCSX(0,-1,0);
  TVector3 TCSZ(TMath::Sin(HRSAngle),0,TMath::Cos(HRSAngle));
  TVector3 TCSY = TCSZ.Cross(TCSX);
  fTCSInHCS.RotateAxes(TCSX,TCSY,TCSZ);

  TVector3 fPointingOffset;
  fPointingOffset.SetXYZ(-MissPointZ*TMath::Sin(HRSAngle)*TMath::Cos(HRSAngle),(Double_t)MissPointY,MissPointZ*TMath::Sin(HRSAngle)*TMath::Sin(HRSAngle));


  TVector3 BeamSpotHCS_average(BeamX_average,BeamY_average,targetfoils[FoilID]);
  TVector3 BeamSpotTCS_average = fTCSInHCS.Inverse()*(BeamSpotHCS_average-fPointingOffset);

  for( Int_t row = 0; row<NSieveRow; row++){
    for( Int_t col = 0; col<NSieveCol; col++){


      Int_t Hole = Get_Hole(col, row);
      
      TVector3 Hole_pos = GetSieveHoleTCS(Hole);

       TVector3 Hole_pos_nocorrect = GetSieveHoleTCS(col,row);
      
      TVector3 MomDirectionTCS_hole = Hole_pos - BeamSpotTCS_average;
      
      Double_t theta_hole = MomDirectionTCS_hole.X()/MomDirectionTCS_hole.Z();
      Double_t phi_hole = MomDirectionTCS_hole.Y()/MomDirectionTCS_hole.Z();
      //      
      //      cout << "For row " << row << ", column " << col << " theta = " << theta_hole << ", phi = " << phi_hole << ", Hole_pos.Y() = " << Hole_pos.Y() << ", BeamSpotTCS_average = " << BeamSpotTCS_average.Y() << endl;

      cout << "For row " << row << ", column " << col << " theta = " << theta_hole << ", phi = " << phi_hole << ", Hole_pos_nocorrect.Y() = " << Hole_pos_nocorrect.Y() << ", Hole_pos.Y() = " << Hole_pos.Y() << ", BeamSpotTCS_average = " << BeamSpotTCS_average.Y() << endl;
      
      
    }
  }
}
