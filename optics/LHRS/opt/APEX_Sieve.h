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
#include "InputAPEXL.h"


std::vector<int> Get_Col_Row(Int_t Hole);



const TVector3 GetSieveHoleTCS(Int_t Col, Int_t Row) /*const*/
{

  // Double_t new_Z = ZPos;
  // Double_t X_SH = SieveXbyRow[Row];
  // cout << "Col = " << Col << endl;

  // //  Double_t Y_SH = SieveYbyCol[Col];
  // cout << "X_SH = " << X_SH << " and Y_SH = " << 999 << " and new_Z = " << ZPos << endl;


  // Double_t Y_SH = SieveYbyCol[Col];
  // X_SH = SieveXbyRow[Row];
  // Y_SH = SieveYbyCol[Col];

  // cout << "X_SH = " << X_SH << " and Y_SH = " << Y_SH << endl;


  //  Double_t Sieve_Z = ZPos - 


  if(Col < 0 || Col > NSieveCol){
    
    Col = 1;
  }


  if(Row < 0 || Row > NSieveRow){
    
    Row = 1;
  }


  TVector3 SieveHoleTCS(SieveXbyRow[Row]+SieveOffX, SieveYbyCol[Col]+SieveOffY, ZPos);

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


//  cout << "For Hole " << Hole << ": row = " << row << " and col = " << col << endl;
    //  hole_no += Col;
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

