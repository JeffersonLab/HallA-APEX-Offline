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


  //  Int_t foil_list[] = {8,9}; // foils to be ignored from optimisation
  Int_t foil_list[] = {}; // foils to be ignored from optimisation

  Int_t col_list[] = {}; // columns to be missed from optimisation

  
 
  
  
  //  holes_fcr hole_list[] = {{10,0,4},{10,0,12},{10,1,5},{10,2,4},{10,2,6},{10,2,12},{10,3,5},{10,4,6},{10,4,12},{10,5,7},{10,6,8},{10,6,10},{10,6,12}};
  //  holes_fcr hole_list[] = {{10,1,5},{10,1,7},{10,1,9},{10,1,11},{10,2,4},{10,2,6},{10,2,8},{10,2,10},{10,2,12},{10,3,5},{10,3,7},{10,3,9},{10,3,11},{10,4,6},{10,4,8},{10,4,10},{10,4,12},{10,5,7},{10,5,9},{10,5,11},{10,6,8},{10,6,10},{10,6,12}};
  //  holes_fcr hole_list[] = {{10,2,4},{10,2,6},{10,2,8},{10,2,10},{10,2,12},{10,3,5},{10,3,7},{10,3,9},{10,3,11},{10,4,6},{10,4,8},{10,4,10},{10,4,12},{10,5,7},{10,5,9},{10,5,11},{10,6,8},{10,6,10},{10,6,12}};
  //  holes_fcr hole_list[] = {};
  
  holes_fcr hole_list[] = {{9,11,5},{9,11,11},{9,12,4},{9,12,6},{9,12,10},{9,14,4},{9,14,6},{9,14,10},{10,3,11},{10,4,10}};


  /// columns defined by foil and column to be ignored
  cols_fc col_foil_list[] = {{8,6},{8,7},{8,8},{8,9},{8,18},{8,19},{8,20},{9,2},{9,15},{9,16},{9,17},{9,18},{10,5},{10,6}};


  /// rows defined by foil and row to be ignored
  rows_fr row_foil_list[] = {{8,3},{8,12},{9,3},{9,12},{9,13},{10,3},{10,4},{10,5},{10,12},{10,13}};
  

  Int_t row_list[] = {13,14,15}; // rows to be missed from optimisation



  // set holes to be ignored based foil, column and hole

  for(Int_t i = 0; i<11; i++){
    for(Int_t j = 0; j<NSieveCol; j++){
      for(Int_t k = 0; k<NSieveRow; k++){      
	
	hole_select[i][j][k] = true;
	
	
	for(const holes_fcr hole_ign : hole_list){

	  if(i == hole_ign.Foil && j == hole_ign.Col && k == hole_ign.Row){
	    hole_select[i][j][k] = false;
	    // condition checks if hole (from its foil, column and row number) is on list of holes to be exluded in optimisation (and sets bool to fales if it is to be excluded)
	  }	  
	}

	
      }
    }
  }
  
  
  // setting foils to be ignored

  for(Int_t i = 0; i<11; i++){
    
    foil_select[i] = true;

    for(const Int_t foil_ign : foil_list){
      
      if( i == foil_ign){
    	foil_select[i] = false;
      }
    }
    
  }



  // next loop sets holes and columns to be ignored during optimisation

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

 

  for(Int_t i = 0; i<11; i++){
    for(Int_t j = 0; j<NSieveCol; j++){
    
      col_select[i][j] = true;

      // this loop sets columsn to be ignored based on column list (will ignore for all foils)
      for(const Int_t col_ign : col_list){
	
	if( j == col_ign){
	  col_select[i][j] = false;
	}      	
      }

      // this loop sets columsn to be ignored based on column list (will ignore for all foils)
      for(const cols_fc col_ign : col_foil_list){
	
	if( i==col_ign.Foil && j == col_ign.Col){
	  col_select[i][j] = false;
	}      	
      }      


      
    }

  }



  

  for(Int_t i = 0; i<11; i++){
    for(Int_t j = 0; j<NSieveRow; j++){
    
      row_select[i][j] = true;

      // this loop sets rows to be ignored based on row list (will ignore for all foils)
      for(const Int_t row_ign : row_list){
	
	if( j == row_ign){
	  row_select[i][j] = false;
	}      	
      }

      // this loop sets rows to be ignored based on row-foil list (will ignore row for particular foil)
      for(const rows_fr row_ign : row_foil_list){
	
	if( i==row_ign.Foil && j == row_ign.Row){
	  row_select[i][j] = false;
	}      	
      }      


      
    }

  }


}


