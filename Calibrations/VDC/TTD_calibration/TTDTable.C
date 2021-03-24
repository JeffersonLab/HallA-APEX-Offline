#include "TTDTable.h"



// performs lookup table ttd conversion without angular correction
double TTDTable::Convert(double dtime){

  Double_t bin_res = 0.5e-9;

  dtime -= Low;
  Int_t bin_no = dtime/(bin_res);
  Double_t dist = LTable[bin_no];
  
  return dist;  
}




// performs lookup table ttd conversion with angular correction
// obselete version
/*
double TTDTable::ConvertAngleCorr(double dtime, double tanTheta){

  double dist = Convert(dtime);
  
 
  Double_t a1 = 0.0, a2 = 0.0;

  const Double_t* fA1tdcCor = &(par[0]);
  const Double_t* fA2tdcCor = &(par[4]);

  
  if( fabs(tanTheta) > 1e-7 ){
    //   tanTheta = 1.0 / tanTheta;
    tanTheta = 1.0 / tanTheta - invTanTheta0;
  } else {
    //    cerr << "TTDform: Invalid tanTheta = " << tanTheta << endl;
    return 0.0;
  }

  for (Int_t i = 3; i >= 1; i--) {
    a1 = tanTheta * (a1 + fA1tdcCor[i]);
    a2 = tanTheta * (a2 + fA2tdcCor[i]);
  }
  a1 += fA1tdcCor[0];
  a2 += fA2tdcCor[0];


  if (dist < 0) {
    //    cerr << "ttdForm: invalid dist = " << dist << endl;
    return 0;
    // } else if (a2<0 || a1<0 || dist < 0) {

  //   return 1e32;
  } else if (dist < a1 ) { 
    dist *= ( 1.0 + a2 / a1);
  }  else {
    dist +=  a2;
  }
      
  return dist;

}
*/


// Angular correction for TTD lookup value
// here tanTheta argument is given as slope of cluster 
double TTDTable::ConvertAngleCorr(double dtime, double tanTheta){


  double dist = Convert(dtime); // get lookup table value for central angle
  
   
  const double* fR = &(par[0]); // radius of electric field of wire
  //  const double* fslope0 = &(par[1]); // central angle
  const double* fslope0 = &(par[1]); // central angle
  
  
  const double fTheta0 = TMath::ATan(*fslope0);
  //  const double fTheta0 = 1/(*fslope0);
  /* const double fTheta0 = (*fslope0); */

  tanTheta = TMath::ATan(tanTheta);
  //  tanTheta = 1/tanTheta;
  
  if (dist >= *fR) {
    dist += (*fR)*( (1/TMath::Cos(tanTheta) -  (1/TMath::Cos(fTheta0))));

  } else if (dist < *fR ) { 
    
    dist *= (1/TMath::Cos(tanTheta) /  (1/TMath::Cos(fTheta0)) );
  }
		 
  return dist;
}


// Print angular correction parameters
void TTDTable::PrintParams(){
  
  cout << "R = " << par[0] << endl;
  cout << "Theta0 = " << par[1] << endl;
}


