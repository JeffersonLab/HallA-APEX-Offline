#include "TString.h"
#include "../Load_more_rootfiles.C"
#include "../InputAPEXL.h"

void Vertical_foil_plots(){

  // ----------------------------------------------------------------------
  //

  Int_t run_nos[] = {4766,4768,4769};
  
  TChain* T;

  TChain* T_ind[3];

  
  Int_t run_po = 0;
  for(auto run_no : run_nos){

    TChain* T_add = Load_more_rootfiles(run_no);

       
    // *T_ind[run_po] = *T_add;
    
    if(run_no == run_nos[0]){
      T = T_add;
    }
    else{
      T->Add(T_add);
    }

    cout << "T->GetEntries() = " << T->GetEntries() << endl;

    run_po++;
  }

  

  
  TCanvas *c1 = new TCanvas("c1"); 

  TH2F* h2 = new TH2F("h2","FP Vertex plot", 400, -0.05, 0.05, 200,-0.05, 0.04);


  h2->GetXaxis()->SetTitle("y_{FP} [m]");
  h2->GetYaxis()->SetTitle("#phi_{FP} [rad]");

  TCut GenrealCut = GeneralSieveCut + PID_cuts + FP_cuts;
  
  T->Draw("L.tr.r_ph:L.tr.r_y>>h2",GenrealCut, "COLZ");


  TCanvas *c2 = new TCanvas("c2");

  TH2F* h_ind[3];

  c2->Divide(3,1);
  
  for (int i=0; i<3; ++i){

    h_ind[i] = new TH2F(Form("h2_%d",i),"FP Vertex plot", 400, -0.05, 0.05, 200,-0.05, 0.04);

    c2->cd(i+1);

    TCut beam_cut(get_Beamcut(run_nos[i]));
    T->Draw(Form("L.tr.r_ph:L.tr.r_y>>h2_%d",i),beam_cut + GenrealCut, "COLZ");

  }

 
  
  

  
}
    
  


