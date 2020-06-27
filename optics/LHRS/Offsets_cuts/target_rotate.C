#include "InputAPEXL.h"



void target_rotate(){

  // define Vertical target foils (set y to 0)

  TVector3 Vert_foils[3];

  Int_t offset = 8; // offset for vertical schem where first 8 foils are Optics Foils


  cout << "All values in mm" << endl;
  
  for(Int_t i = 0; i<3; i++){
    
    Vert_foils[i].SetXYZ(orig_targetfoilsX[i+offset]*1000,0*1000,targetfoils[i+offset]*1000);
   
    cout  << Single_foil[i] << ": Vert_foils[" << i << "] = [" << Vert_foils[i].X() << "," << Vert_foils[i].Y() << "," << Vert_foils[i].Z() << "]" << endl << endl;;


    // adding yaw

    Double_t new_x = Vert_foils[i].X()*TMath::Cos(target_yaw) + Vert_foils[i].Z()*TMath::Sin(target_yaw);
    Double_t new_z = -Vert_foils[i].X()*TMath::Sin(target_yaw) + Vert_foils[i].Z()*TMath::Cos(target_yaw);
    
    Vert_foils[i].RotateY(target_yaw);
    
    cout << "Adding target_yaw..." << endl;

    cout  << "Manual target_yaw " << Single_foil[i] << ": Vert_foils[" << i << "] = [" << new_x << "," << 0  << "," << new_z << "]" << endl;
    
    cout  << "Artificial target_yaw " << Single_foil[i] << ": Vert_foils[" << i << "] = [" << Vert_foils[i].X() << "," << Vert_foils[i].Y() << "," << Vert_foils[i].Z() << "]" << endl << endl << endl;
  }


  
  
  


}
