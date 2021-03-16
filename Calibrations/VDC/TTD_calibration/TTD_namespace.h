#ifndef TTD_NAME_H
#define TTD_NAME_H

// Free function for saving TTD calibraton

namespace TTD_func

{



  /* perform analytic TTD conversion based on difference from slope to central slope */
  double TTDform( double dtime, double tanTheta, double fDriftVel, const double *par, double invTanTheta0 ){
  double a1 = 0.0, a2 = 0.0;

  const double* fA1tdcCor = &(par[0]);
  const double* fA2tdcCor = &(par[4]);

  if( fabs(tanTheta) > 1e-7 ){
    tanTheta = 1.0 / tanTheta - invTanTheta0;
  } else {
    cerr << "TTDform: Invalid tanTheta = " << tanTheta << endl;
    return 0.0;
  }

  for (Int_t i = 3; i >= 1; i--) {
    a1 = tanTheta * (a1 + fA1tdcCor[i]);
    a2 = tanTheta * (a2 + fA2tdcCor[i]);
  }
  a1 += fA1tdcCor[0];
  a2 += fA2tdcCor[0];

  double dist = fDriftVel * dtime;

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


    /* perform analytic TTD conversion based on difference from slope to central slope */

  Double_t TTDAna( Double_t dtime, Double_t tanTheta, Double_t fDriftVel, const Double_t *par){
  
  Double_t a1 = 0.0, a2 = 0.0;

  const Double_t* fA1tdcCor = &(par[0]);
  const Double_t* fA2tdcCor = &(par[4]);

  if( fabs(tanTheta) > 1e-7 ){
    tanTheta = 1.0 / tanTheta;
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

  Double_t dist = fDriftVel * dtime;

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



  
  template <typename T>
  int SaveNewTTDData(vector<T> table, Double_t NBins, TString arm, TString plane, Int_t runnumber, TString name){

    
    TString filename = Form("DB/lookup_tables/db_%s_%s_lookup_TTD_%s.vdc.%d.dat",arm.Data(), plane.Data(), name.Data(), runnumber);
  
    std::ofstream* outp = new std::ofstream;

    outp->open(filename.Data());

    *outp << Form("%s.vdc.%s.ttd_table.nbins = ",arm.Data(),plane.Data()) << endl <<NBins << endl << endl;


    *outp << Form("%s.vdc.%s.ttd_table.table =",arm.Data(),plane.Data()) << endl;

    for(Int_t j=0; j<NBins; j++){
      if (j%10 == 0 && j>0){
	*outp << endl;
      }
      *outp << table[j] << " ";
    }
    *outp << endl;

  
    return 1;
  }


  
  bool passtrg(Int_t evttype, Int_t trg){
    //  cout<<evttype<<"   "<<trg<<endl;
    return evttype&(1<<trg);
  }



  std::vector<Double_t> ReadLookupTable(std::ifstream& DBfile, Int_t line_no, Int_t NBins){

    std::string line;
    
    Int_t counter = 0;
    
    std::vector<Double_t> table;
    
    while(std::getline(DBfile, line) && counter < NBins){
      
      std::istringstream ss(line);
      
      Double_t x;
      while(ss>>x && counter < NBins){	
	table.push_back(x);
	counter++;
      }            
    }
    
    return table;
  }



  template <typename T>
T ReadSingleVal(std::string line){


  std::istringstream iss(line);

  T SingeVal = 0;
  
  iss >> SingeVal;

  return SingeVal;


  }


  Int_t ReadNBins(std::string line){
    
    
  Int_t NBins = 0.0;
  
  NBins = ReadSingleVal<Int_t>(line);
  
  return NBins;
  
  }
  







/*

  perform analytic TTD conversion based on difference from slope to central slope.
  Different from standard Hall A method 

  Based on two MIT theses:
  - V. Jordan, MIT, 1994
  - W. Schmitt, MIT, 1993



  */

double TTD_Corr( double dist, double tanTheta, const double *par){
  
  
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



}

#endif
