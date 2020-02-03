////////////////////////////////////////////////////
//          LoadDatabase
//
//   Copied from OpticOpt Script:
//   - designed to Load VDC DataBase matrix 
//    coeffecients into matrices
//    
//   
//////////////////////////////////////////////////


#ifndef ROOT_LoadDB
#define ROOT_LoadDB
//#include "THaTrackingDetector.h"
#include "TRotation.h"
#include "TMath.h"
#include <vector>
#include "DebugDef.h"
#include "THaString.h"
#include <iostream>
#include <fstream>



void LoadDataBase(TString DataBaseName);

enum { kPORDER = 7 };


// class for storing matrix element data
class THaMatrixElement{
 public:
 THaMatrixElement() : iszero(true), pw(3), order(0), v(0), poly(kPORDER),OptOrder(0) {}
  bool match( const THaMatrixElement& rhs ) const;
  
  bool iszero;             // whether the element is zero
  std::vector<int> pw;     // exponents of matrix element
  //   e.g. D100 = { 1, 0, 0 }
  int  order; // order of x?
  double v;                // its computed value
  std::vector<double> poly;// the associated polynomial
  
  void SkimPoly(); //reduce order to highest non-zero poly
  
  UInt_t OptOrder;		//order optimize to

  /* std::vector<double> values; //  */
  //  double values[5][5][5][5] = {{{{0}}}};
  

  // define true/false matrix for existence of elements in DB ie does (x^i)(theta^j)(y^k)(phi^l) exist

  //  Int_t Exist[5][5][5][5] = {{{{0}}}};

  


};

Double_t powers[4][5] = {0};


void CalcMatrix(const double x, vector<THaMatrixElement>& matrix);












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


// structs to hold matrix elems with other info (all values for all of one kind of matrix elements (all T elements)

/* double values[5][5][5][5] = {{{{0}}}}; */
/* Int_t Exist[5][5][5][5] = {{{{0}}}}; */


struct DB_entry{
  // stores database coeffecient with indices of FP variables with which its multiplied

  Int_t x_p = 0;
  Int_t th_p = 0;
  Int_t y_p = 0;
  Int_t ph_p = 0;


  Double_t DB_coeff = 0.0;
};



struct MatrixElems_Vals{
  //  std::vector<THaMatrixElement> fTMatrixEl;
  //  double values[5][5][5][5] = {{{{0}}}};

  // vector with values of all matrix elems (for one event)
  vector< double> values = {0};

  // 4D array with 0/1 represnting if ME is used in optimisation (non-zero)
  Int_t Exist[5][5][5][5] = {{{{0}}}};

  // total number of parameters that are non-zero
  Int_t no_elements = 0;

  
  Double_t tg_v = 0;


  vector< DB_entry> DB_entries;
};


// holders for matrix elemtns of T,Y,P and D

MatrixElems_Vals TMatrix_vals;
MatrixElems_Vals YMatrix_vals;
MatrixElems_Vals PMatrix_vals;
MatrixElems_Vals DMatrix_vals;




//void LoadDataBase(TString DataBaseName);

TString OldComments;

auto fPrefix = new char [1000];


#endif


