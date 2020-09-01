#include "TString.h"
#include "../Load_more_rootfiles.C"
#include "../InputAPEXL.h"

void Vertical_foil_plots(){

  // ----------------------------------------------------------------------
  //

  Int_t run_nos[] = {4766,4768,4769};
  
  TChain* T;

  TChain* T_ind[3];

  TFile* f1[3];

  f1[0] = new TFile("/w/work3/home/johnw/Rootfiles/apex_4766_15_8_2020_V1.root.FullCut.root","READ");
  f1[1] = new TFile("/w/work3/home/johnw/Rootfiles/apex_4768_15_8_2020_V2.root.FullCut.root","READ");
  f1[2] = new TFile("/w/work3/home/johnw/Rootfiles/apex_4769_15_8_2020_V3.root.FullCut.root","READ");

  TCutG* foil_tg_gcuts[3];
  TCutG* foil_FP_gcuts[3];
  
  Int_t run_po = 0;
  for(auto run_no : run_nos){

    foil_tg_gcuts[run_po] = (TCutG*) f1[run_po]->GetObjectChecked(Form("fcut_L_%d", run_po + 8), "TCutG"); //looking for foil cut definition

    foil_FP_gcuts[run_po] = (TCutG*) f1[run_po]->GetObjectChecked(Form("fcut_L_FP_%d", run_po + 8), "TCutG"); //looking for foil cut definition

    
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

  for(Int_t i = 0; i<3; i++){
    foil_FP_gcuts[i]->SetLineColor(kRed);
    foil_FP_gcuts[i]->Draw("PL same");
  }

  c1->Print("vertex_results/Vert_FP_foils_all.pdf");


  TCanvas *c2 = new TCanvas("c2");

  TH2F* h_ind[3];

  c2->Divide(3,1);
  
  for (int i=0; i<3; ++i){

    h_ind[i] = new TH2F(Form("h2_%d",i),"FP Vertex plot", 400, -0.05, 0.05, 200,-0.05, 0.04);

    c2->cd(i+1);

    TCut beam_cut(get_Beamcut(run_nos[i]));
    T->Draw(Form("L.tr.r_ph:L.tr.r_y>>h2_%d",i),beam_cut + GenrealCut, "COLZ");

    foil_FP_gcuts[i]->SetLineColor(kRed);
    foil_FP_gcuts[i]->Draw("PL same");

  }

  c2->Print("vertex_results/Vert_FP_foils_ind.pdf");

  TCanvas *c3 = new TCanvas("c3");

  TH2F* h_f_tg = new TH2F("h_f_tg", "Vertex Tg", 400, -0.05, 0.05, 200,-0.5, 0.5);
  h_f_tg->GetXaxis()->SetTitle("\\phi_{tg} [rad]");
  h_f_tg->GetYaxis()->SetTitle("Reactz [m]");

  T->Draw("reactz:ph_tgt>>h_f_tg", GenrealCut, "COLZ");

  for(Int_t i = 0; i<3; i++){
    foil_tg_gcuts[i]->SetLineColor(kRed);
    foil_tg_gcuts[i]->Draw("PL same");
  }

  c3->Print("vertex_results/Vert_tg_foils_all.pdf");
  
  
  TCanvas *c4 = new TCanvas("c4");
  c4->Divide(3,1);
  
  TH2F* h_f_tg_ind[3];

  for (int i=0; i<3; ++i){

    c4->cd(i+1);
    
    h_f_tg_ind[i] = new TH2F(Form("h_f_tg_ind_%d",i+8), Form("theta_target vs. phi_target, Foil %d", i+8),  400, -0.05, 0.05, 200,-0.5, 0.5);

    h_f_tg_ind[i]->GetXaxis()->SetTitle("#phi_{tg} [rad]");
    h_f_tg_ind[i]->GetYaxis()->SetTitle("Reactz [m]");

    TCut beam_cut(get_Beamcut(run_nos[i]));
    T->Draw(Form("reactz:ph_tgt>>h_f_tg_ind_%d",i+8),beam_cut + GenrealCut, "COLZ");
    
    foil_tg_gcuts[i]->SetLineColor(kRed);
    foil_tg_gcuts[i]->Draw("PL same");


  }

  c4->Print("vertex_results/Vert_tg_foils_ind.pdf");

  TCanvas *c5 = new TCanvas("c5");
  c5->Divide(3,1);

  TH2F* h_tg_sieve[3];

  for (int i=0; i<3; ++i){

    c5->cd(i+1);
    Int_t FoilID = i + 8;
    
    h_tg_sieve[i] = new TH2F(Form("h_tg_sieve_%d",i+8), Form("theta_target vs. phi_target, Foil %d", FoilID),  400,-30,30,400,-65,65);

    h_tg_sieve[i]->GetXaxis()->SetTitle("#phi_{tg} [mrad]");
    h_tg_sieve[i]->GetYaxis()->SetTitle("#theta_{tg} [mrad]");

    TCut beam_cut(get_Beamcut(run_nos[i]));
    T->Draw(Form("(1000*th_tgt):(1000*ph_tgt)>>h_tg_sieve_%d",i+8),beam_cut + GenrealCut + TCut(Form("fcut_L_%d", FoilID)) + TCut(Form("fcut_L_FP_%d", FoilID)), "COLZ");
    
    // foil_tg_gcuts[i]->SetLineColor(kRed);
    // foil_tg_gcuts[i]->Draw("PL same");


  }

  c5->Print("vertex_results/Vert_tg_sieve_ind.pdf");
  
  
}
    
  


