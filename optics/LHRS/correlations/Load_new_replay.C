////////////////////////////////////////////////////
//          Load_new_replay
//
//   Script designed to load new replay, 
//   created by apex_ME_calc.C script which
//   stores matrix variables etc
//
//
//  John Williamson
//  25/10/2019
///////////////////////////////////////////////////


// #ifndef Load_new_replay
// #define Load_new_replay

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>
#include "file_def.h"
#include <iostream>
//#include "Load_new_replay.h"



TChain* Load_new_replay(TString DB_name, int runnum_1=4179, int runnum_2= -1){

 
  TChain *T = new TChain("T");
  

  
  TString filenamebase;
  TString filename;

  Long_t split = 0;

  Int_t No_of_runs = 0;

  if( runnum_2 == -1){
    No_of_runs = 0;
  }
  else{
    No_of_runs = runnum_2-runnum_1;
  }



  for(Int_t i = 0; i< No_of_runs+1; i++){
    

    filenamebase = "rootfiles/" + DB_name + Form("/%d_replay.root",runnum_1+i);
    filename = filenamebase;
    
    cout << "Current runnumber = " << runnum_1+i << endl;
    
    filenamebase.Remove(filenamebase.Last('.'),5);
        
    split = 0;
    while ( !gSystem->AccessPathName(filename.Data()) ) {
      cout << "Added root file: " << filename << " to Tree" << endl;
      T->Add(filename);
      split++;
      filename = filenamebase + "_" + split + ".root";
    }
    
  }
  
  if( runnum_2 == -1){
    
    cout << "Opened Runs " << runnum_1 << " with "  << T->GetEntries() << " events"<<  endl;
  }
  else{
    
    cout << "Opened Runs " << runnum_1 << " to " << runnum_2 <<  " with "  << T->GetEntries() << " events"<<  endl;
  }
  

  return T;
  


}


//#endif
