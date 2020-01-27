///////////////////////////////////////////////////////////////////////////////
//
// Allows chosen hole not to be used in optimisation
//
//

// Author: John Williamson <jwilli@jlab.org>
// Modification:
//			Jan 13th 2020
//
///////////////////////////////////////////////////////////////////////////////


#include "hole_opt_select.h"



void hole_opt_select(){


  Int_t foil_list[] = {}; // foils to be ignored from optimisation

  Int_t col_list[] = {}; // columns to be missed from optimisation

  Int_t hole_list[] = {};
  //  Int_t hole_list[] = {79,105,131};



  //  Int_t row_list = [0,1]; // rows to be missed from optimisation

  
  // setting foils to be ignored

  for(Int_t i = 0; i<NFoils; i++){
    
    foil_select[i] = true;

    for(const Int_t foil_ign : foil_list){
      
      if( i == foil_ign){
    	foil_select[i] = false;
      }
    }
    
  }



  // next loop sets holes and columsn to be ignored during optimisation

  // for(Int_t i = 0; i<NHoles; i++){
    
  //   hole_select[i] = true;

  //  for(const Int_t hole_ign : hole_list){
      
  //     if( i == hole_ign){
  //   	hole_select[i] = false;
  //     }
  //   }
    

   /// / setting columns to be ignored
    // std::vector<int> x_y = {};
    // x_y = Get_Col_Row(i);
    
    // UInt_t Col = x_y[0];
    // UInt_t Row = x_y[1];
    


    
    // for(const Int_t col_ign : col_list){
      
    //   if( Col == col_ign){
    // 	hole_select[i] = false;
    //   }
    // }
   
  //  }




  for(Int_t i = 0; i<NSieveCol; i++){
    
    col_select[i] = true;

    for(const Int_t col_ign : col_list){
      
      if( i == col_ign){
    	col_select[i] = false;
      }
    }
    
  }




}


