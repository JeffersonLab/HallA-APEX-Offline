////////////////////////////////////////////////////
//          corr_plot
//
//   Script designed to plot, for a given tree,
//   various correlation plots between target 
//   and focal-plane variables
//
//
//  John Williamson
//  25/10/2019
///////////////////////////////////////////////////


//#include "Load_new_replay.C"
#include "InputAPEXL.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "file_def.h"
#include "TSystem.h"
#include <iostream>
#include "plot_parameters.h"

//void corr_plot(TString DB_name, TChain* T, Int_t runnumber = 4179, Int_t col = -1, Int_t row = -1){

//void corr_plot(DB_info DB_run, event_info event = {4179,1,-1,-1}){

void corr_plot(TString DB_name, TChain* T, Int_t runnumber = 4179, Int_t col = -1, Int_t row = -1){

  gROOT->SetBatch(kTRUE);

  gStyle->SetOptStat(0);

  /* TChain* T = Load_new_replay(4179); */

  /* Int_t runnumber = 4179; */



  // define histogram limits 
  // (could be altered to vary depending on run)

  
  Double_t xfp_h = 0;
  Double_t xfp_l = 0;

  Double_t yfp_h = 0;
  Double_t yfp_l = 0;

  Double_t thfp_h = 0;
  Double_t thfp_l = 0;

  Double_t phfp_h = 0;
  Double_t phfp_l = 0;


  Double_t ytg_h = 0;
  Double_t ytg_l = 0;

  Double_t thtg_h = 0;
  Double_t thtg_l = 0;

  Double_t phtg_h = 0;
  Double_t phtg_l = 0;

  Double_t dptg_h = 0;
  Double_t dptg_l = 0;

  Double_t zreact_h = 0.5;
  Double_t zreact_l = -0.5;


  // if(runnumber == 4179){

    xfp_h = 0.8;
    xfp_l = -0.8;

    yfp_h = 0.05;
    yfp_l = -0.05;
    
    thfp_h = 30;
    thfp_l = -35;

    phfp_h = 50;
    phfp_l = -50;

    if(runnumber == 4766){
      // V1
      ytg_h = 0.04;
      ytg_l = 0.00;
    }
    else if(runnumber == 4768){
      // V2
      ytg_h = 0.02;
      ytg_l = -0.02;
    }
    else if(runnumber == 4769){
      // V3
      ytg_h = 0.00;
      ytg_l = -0.04;
    }
    else{
      ytg_h = 0.05;
      ytg_l = -0.05;
    }

    thtg_h = 70;
    thtg_l = -50;
    
    phtg_h = 40;
    phtg_l = -50;

    dptg_h = 0.05;
    dptg_l = -0.05;
    

  // }



  // plot theta_target against focal plane parameters

  
  TCanvas *c1 = new TCanvas("c1","LHRS theta_target vs FP variables",1000,800); 
  c1->Divide(2,2); 


  c1->cd(1);
  TH2D* Thtgt_v_ThFp = new TH2D("Thtgt_v_ThFp","\\theta_{tg} v \\theta_{FP}",300,thfp_l,thfp_h,300,thtg_l,thtg_h);
  
  Thtgt_v_ThFp->GetYaxis()->SetTitleOffset(1.0);
  Thtgt_v_ThFp->GetXaxis()->SetTitleSize(0.05);
  Thtgt_v_ThFp->GetYaxis()->SetTitleSize(0.05);
  Thtgt_v_ThFp->GetXaxis()->SetTitle("\\theta_{FP}");
  Thtgt_v_ThFp->GetYaxis()->SetTitle("\\theta_{tg}");

  
  T->Draw("L.tr.tg_th*1000:L.tr.r_th*1000>>Thtgt_v_ThFp",PID_cuts + GeneralSieveCut,"colz");

  
  c1->cd(2);
  TH2D* Thtgt_v_PhFp = new TH2D("Thtgt_v_PhFp","\\theta_{tg} v \\phi_{FP}",300,phfp_l,phfp_h,300,thtg_l,thtg_h);
  
  Thtgt_v_PhFp->GetYaxis()->SetTitleOffset(1.0);
  Thtgt_v_PhFp->GetXaxis()->SetTitleSize(0.05);
  Thtgt_v_PhFp->GetYaxis()->SetTitleSize(0.05);
  Thtgt_v_PhFp->GetXaxis()->SetTitle("\\phi_{FP}");
  Thtgt_v_PhFp->GetYaxis()->SetTitle("\\theta_{tg}");

  
  T->Draw("L.tr.tg_th*1000:L.tr.r_ph*1000>>Thtgt_v_PhFp",PID_cuts + GeneralSieveCut,"colz");


  c1->cd(3);
  TH2D* Thtgt_v_xFp = new TH2D("Thtgt_v_xFp","\\theta_{tg} v x_{FP}",300,xfp_l,xfp_h,300,thtg_l,thtg_h);
  
  Thtgt_v_xFp->GetYaxis()->SetTitleOffset(1.0);
  Thtgt_v_xFp->GetXaxis()->SetTitleSize(0.05);
  Thtgt_v_xFp->GetYaxis()->SetTitleSize(0.05);
  Thtgt_v_xFp->GetXaxis()->SetTitle("x_{FP}");
  Thtgt_v_xFp->GetYaxis()->SetTitle("\\theta_{tg}");

  
  T->Draw("L.tr.tg_th*1000:L.tr.r_x>>Thtgt_v_xFp",PID_cuts + GeneralSieveCut,"colz");


  c1->cd(4);
  TH2D* Thtgt_v_yFp = new TH2D("Thtgt_v_yFp","\\theta_{tg} v y_{FP}",300,yfp_l,yfp_h,300,thtg_l,thtg_h);
  
  Thtgt_v_yFp->GetYaxis()->SetTitleOffset(1.0);
  Thtgt_v_yFp->GetXaxis()->SetTitleSize(0.05);
  Thtgt_v_yFp->GetYaxis()->SetTitleSize(0.05);
  Thtgt_v_yFp->GetXaxis()->SetTitle("y_{FP}");
  Thtgt_v_yFp->GetYaxis()->SetTitle("\\theta_{tg}");

  
  T->Draw("L.tr.tg_th*1000:L.tr.r_y>>Thtgt_v_yFp",PID_cuts + GeneralSieveCut,"colz");


  // plot phi_target against focal plane parameters

  TCanvas *c2 = new TCanvas("c2","LHRS phi_target vs FP variables",1000,800); 
  c2->Divide(2,2); 


  c2->cd(1);
  TH2D* Phtgt_v_ThFp = new TH2D("Phtgt_v_ThFp","\\phi_{tg} v \\theta_{FP}",300,thfp_l,thfp_h,300,phtg_l,phtg_h);
  
  Phtgt_v_ThFp->GetYaxis()->SetTitleOffset(1.0);
  Phtgt_v_ThFp->GetXaxis()->SetTitleSize(0.05);
  Phtgt_v_ThFp->GetYaxis()->SetTitleSize(0.05);
  Phtgt_v_ThFp->GetXaxis()->SetTitle("\\theta_{FP}");
  Phtgt_v_ThFp->GetYaxis()->SetTitle("\\phi_{tg}");

  
  T->Draw("L.tr.tg_ph*1000:L.tr.r_th*1000>>Phtgt_v_ThFp",PID_cuts + GeneralSieveCut,"colz");

  
  c2->cd(2);
  TH2D* Phtgt_v_PhFp = new TH2D("Phtgt_v_PhFp","\\phi_{tg} v \\phi_{FP}",300,phfp_l,phfp_h,300,phtg_l,phtg_h);
  
  Phtgt_v_PhFp->GetYaxis()->SetTitleOffset(1.0);
  Phtgt_v_PhFp->GetXaxis()->SetTitleSize(0.05);
  Phtgt_v_PhFp->GetYaxis()->SetTitleSize(0.05);
  Phtgt_v_PhFp->GetXaxis()->SetTitle("\\phi_{FP}");
  Phtgt_v_PhFp->GetYaxis()->SetTitle("\\phi_{tg}");

  
  T->Draw("L.tr.tg_ph*1000:L.tr.r_ph*1000>>Phtgt_v_PhFp",PID_cuts + GeneralSieveCut,"colz");


  c2->cd(3);
  TH2D* Phtgt_v_xFp = new TH2D("Phtgt_v_xFp","\\phi_{tg} v x_{FP}",300,xfp_l,xfp_h,300,phtg_l,phtg_h);
  
  Phtgt_v_xFp->GetYaxis()->SetTitleOffset(1.0);
  Phtgt_v_xFp->GetXaxis()->SetTitleSize(0.05);
  Phtgt_v_xFp->GetYaxis()->SetTitleSize(0.05);
  Phtgt_v_xFp->GetXaxis()->SetTitle("x_{FP}");
  Phtgt_v_xFp->GetYaxis()->SetTitle("\\phi_{tg}");

  
  T->Draw("L.tr.tg_ph*1000:L.tr.r_x>>Phtgt_v_xFp",PID_cuts + GeneralSieveCut,"colz");


  c2->cd(4);
  TH2D* Phtgt_v_yFp = new TH2D("Phtgt_v_yFp","\\phi_{tg} v y_{FP}",300,yfp_l,yfp_h,300,phtg_l,phtg_h);
  
  Phtgt_v_yFp->GetYaxis()->SetTitleOffset(1.0);
  Phtgt_v_yFp->GetXaxis()->SetTitleSize(0.05);
  Phtgt_v_yFp->GetYaxis()->SetTitleSize(0.05);
  Phtgt_v_yFp->GetXaxis()->SetTitle("y_{FP}");
  Phtgt_v_yFp->GetYaxis()->SetTitle("\\phi_{tg}");

  
  T->Draw("L.tr.tg_ph*1000:L.tr.r_y>>Phtgt_v_yFp",PID_cuts + GeneralSieveCut,"colz");


  // plot y_target against focal plane parameters

  TCanvas *c3 = new TCanvas("c3","LHRS y_target vs FP variables",1000,800); 
  c3->Divide(2,2); 


  c3->cd(1);
  TH2D* ytgt_v_ThFp = new TH2D("ytgt_v_ThFp","y_{tg} v \\theta_{FP}",300,thfp_l,thfp_h,300,ytg_l,ytg_h);
  
  ytgt_v_ThFp->GetYaxis()->SetTitleOffset(1.0);
  ytgt_v_ThFp->GetXaxis()->SetTitleSize(0.05);
  ytgt_v_ThFp->GetYaxis()->SetTitleSize(0.05);
  ytgt_v_ThFp->GetXaxis()->SetTitle("\\theta_{FP}");
  ytgt_v_ThFp->GetYaxis()->SetTitle("y_{tg}");

  
  T->Draw("L.tr.tg_y:L.tr.r_th*1000>>ytgt_v_ThFp",PID_cuts + GeneralSieveCut,"colz");

  
  c3->cd(2);
  TH2D* ytgt_v_PhFp = new TH2D("ytgt_v_PhFp","y_{tg} v \\phi_{FP}",300,phfp_l,phfp_h,300,ytg_l,ytg_h);
  
  ytgt_v_PhFp->GetYaxis()->SetTitleOffset(1.0);
  ytgt_v_PhFp->GetXaxis()->SetTitleSize(0.05);
  ytgt_v_PhFp->GetYaxis()->SetTitleSize(0.05);
  ytgt_v_PhFp->GetXaxis()->SetTitle("\\phi_{FP}");
  ytgt_v_PhFp->GetYaxis()->SetTitle("y_{tg}");

  
  T->Draw("L.tr.tg_y:L.tr.r_ph*1000>>ytgt_v_PhFp",PID_cuts + GeneralSieveCut,"colz");


  c3->cd(3);
  TH2D* ytgt_v_xFp = new TH2D("ytgt_v_xFp","y_{tg} v x_{FP}",300,xfp_l,xfp_h,300,ytg_l,ytg_h);
  
  ytgt_v_xFp->GetYaxis()->SetTitleOffset(1.0);
  ytgt_v_xFp->GetXaxis()->SetTitleSize(0.05);
  ytgt_v_xFp->GetYaxis()->SetTitleSize(0.05);
  ytgt_v_xFp->GetXaxis()->SetTitle("x_{FP}");
  ytgt_v_xFp->GetYaxis()->SetTitle("y_{tg}");

  
  T->Draw("L.tr.tg_y:L.tr.r_x>>ytgt_v_xFp",PID_cuts + GeneralSieveCut,"colz");


  c3->cd(4);
  TH2D* ytgt_v_yFp = new TH2D("ytgt_v_yFp","y_{tg} v y_{FP}",300,yfp_l,yfp_h,300,ytg_l,ytg_h);
  
  ytgt_v_yFp->GetYaxis()->SetTitleOffset(1.0);
  ytgt_v_yFp->GetXaxis()->SetTitleSize(0.05);
  ytgt_v_yFp->GetYaxis()->SetTitleSize(0.05);
  ytgt_v_yFp->GetXaxis()->SetTitle("y_{FP}");
  ytgt_v_yFp->GetYaxis()->SetTitle("y_{tg}");

  
  T->Draw("L.tr.tg_y:L.tr.r_y>>ytgt_v_yFp",PID_cuts + GeneralSieveCut,"colz");
  

  

  // plot dp_target against FP variables


  TCanvas *c4 = new TCanvas("c4","LHRS dp_target vs FP variables",1000,800); 
  c4->Divide(2,2); 


  c4->cd(1);
  TH2D* dptgt_v_ThFp = new TH2D("dptgt_v_ThFp","\\deltap_{tg} v \\theta_{FP}",300,thfp_l,thfp_h,300,-0.05,0.05);
  
  dptgt_v_ThFp->GetYaxis()->SetTitleOffset(1.0);
  dptgt_v_ThFp->GetXaxis()->SetTitleSize(0.05);
  dptgt_v_ThFp->GetYaxis()->SetTitleSize(0.05);
  dptgt_v_ThFp->GetXaxis()->SetTitle("\\theta_{FP}");
  dptgt_v_ThFp->GetYaxis()->SetTitle("\\deltap_{tg}");

  
  T->Draw("L.tr.tg_dp:L.tr.r_th*1000>>dptgt_v_ThFp",PID_cuts + GeneralSieveCut,"colz");

  
  c4->cd(2);
  TH2D* dptgt_v_PhFp = new TH2D("dptgt_v_PhFp","\\deltap_{tg} v \\phi_{FP}",300,phfp_l,phfp_h,300,-0.05,0.05);
  
  dptgt_v_PhFp->GetYaxis()->SetTitleOffset(1.0);
  dptgt_v_PhFp->GetXaxis()->SetTitleSize(0.05);
  dptgt_v_PhFp->GetYaxis()->SetTitleSize(0.05);
  dptgt_v_PhFp->GetXaxis()->SetTitle("\\phi_{FP}");
  dptgt_v_PhFp->GetYaxis()->SetTitle("\\deltap_{tg}");

  
  T->Draw("L.tr.tg_dp:L.tr.r_ph*1000>>dptgt_v_PhFp",PID_cuts + GeneralSieveCut,"colz");


  c4->cd(3);
  TH2D* dptgt_v_xFp = new TH2D("dptgt_v_xFp","\\deltap_{tg} v x_{FP}",300,xfp_l,xfp_h,300,-0.05,0.05);
  
  dptgt_v_xFp->GetYaxis()->SetTitleOffset(1.0);
  dptgt_v_xFp->GetXaxis()->SetTitleSize(0.05);
  dptgt_v_xFp->GetYaxis()->SetTitleSize(0.05);
  dptgt_v_xFp->GetXaxis()->SetTitle("x_{FP}");
  dptgt_v_xFp->GetYaxis()->SetTitle("\\deltap_{tg}");

  
  T->Draw("L.tr.tg_dp:L.tr.r_x>>dptgt_v_xFp",PID_cuts + GeneralSieveCut,"colz");


  c4->cd(4);
  TH2D* dptgt_v_yFp = new TH2D("dptgt_v_yFp","\\deltap_{tg} v y_{FP}",300,-0.05,0.05,300,-0.05,0.05);
  
  dptgt_v_yFp->GetYaxis()->SetTitleOffset(1.0);
  dptgt_v_yFp->GetXaxis()->SetTitleSize(0.05);
  dptgt_v_yFp->GetYaxis()->SetTitleSize(0.05);
  dptgt_v_yFp->GetXaxis()->SetTitle("y_{FP}");
  dptgt_v_yFp->GetYaxis()->SetTitle("\\deltap_{tg}");

  
  T->Draw("L.tr.tg_dp:L.tr.r_y>>dptgt_v_yFp",PID_cuts + GeneralSieveCut,"colz");




  // plot x-target against FP variables 

  /* TCanvas *c5 = new TCanvas("c5","LHRS ",1000,800);  */
  /* c5->Divide(2,2);  */


  /* c5->cd(1); */
  /* TH2D* xtgt_v_ThFp = new TH2D("xtgt_v_ThFp","xtgt_v_\\theta_{FP}",100,-50,70,100,-0.05,0.05); */
  
  /* xtgt_v_ThFp->GetYaxis()->SetTitleOffset(1.0); */
  /* xtgt_v_ThFp->GetXaxis()->SetTitleSize(0.05); */
  /* xtgt_v_ThFp->GetYaxis()->SetTitleSize(0.05); */
  /* xtgt_v_ThFp->GetXaxis()->SetTitle("theta (FP)"); */
  /* xtgt_v_ThFp->GetYaxis()->SetTitle("\\deltap_{tg}"); */

  
  /* T->Draw("L.tr.tg_dp:L.tr.r_th*1000>>xtgt_v_ThFp",PID_cuts + GeneralSieveCut,"colz"); */

  
  /* c5->cd(2); */
  /* TH2D* xtgt_v_PhFp = new TH2D("xtgt_v_PhFp","xtgt_v_\\phi_{FP}",100,phfp_l,phfp_h,100,-0.05,0.05); */
  
  /* xtgt_v_PhFp->GetYaxis()->SetTitleOffset(1.0); */
  /* xtgt_v_PhFp->GetXaxis()->SetTitleSize(0.05); */
  /* xtgt_v_PhFp->GetYaxis()->SetTitleSize(0.05); */
  /* xtgt_v_PhFp->GetXaxis()->SetTitle("\\phi_{FP}"); */
  /* xtgt_v_PhFp->GetYaxis()->SetTitle("\\deltap_{tg}"); */

  
  /* T->Draw("L.tr.tg_dp:L.tr.r_ph*1000>>xtgt_v_PhFp",PID_cuts + GeneralSieveCut,"colz"); */


  /* c5->cd(3); */
  /* TH2D* xtgt_v_xFp = new TH2D("xtgt_v_xFp","xtgt_v_x_{FP}",100,xfp_l,xfp_h,100,-0.05,0.05); */
  
  /* xtgt_v_xFp->GetYaxis()->SetTitleOffset(1.0); */
  /* xtgt_v_xFp->GetXaxis()->SetTitleSize(0.05); */
  /* xtgt_v_xFp->GetYaxis()->SetTitleSize(0.05); */
  /* xtgt_v_xFp->GetXaxis()->SetTitle("x_{FP}"); */
  /* xtgt_v_xFp->GetYaxis()->SetTitle("\\deltap_{tg}"); */

  
  /* T->Draw("L.tr.tg_dp:L.tr.r_x>>xtgt_v_xFp",PID_cuts + GeneralSieveCut,"colz"); */


  /* c5->cd(4); */
  /* TH2D* xtgt_v_yFp = new TH2D("xtgt_v_yFp","xtgt_v_y_{FP}",100,-0.05,0.05,100,-0.05,0.05); */
  
  /* xtgt_v_yFp->GetYaxis()->SetTitleOffset(1.0); */
  /* xtgt_v_yFp->GetXaxis()->SetTitleSize(0.05); */
  /* xtgt_v_yFp->GetYaxis()->SetTitleSize(0.05); */
  /* xtgt_v_yFp->GetXaxis()->SetTitle("y_{FP}"); */
  /* xtgt_v_yFp->GetYaxis()->SetTitle("\\deltap_{tg}"); */

  
  /* T->Draw("L.tr.tg_dp:L.tr.r_y>>xtgt_v_yFp",PID_cuts + GeneralSieveCut,"colz"); */






  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Add target variable canvases (plotting target variables against one another)
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  TCanvas *c6 = new TCanvas("c6","LHRS target variables",1000,800); 
  c6->Divide(2,2); 
  

  c6->cd(1);
  TH1D* y_tgt = new TH1D("y_tgt","y_{tg}",300,-0.06,0.06);

  y_tgt->GetYaxis()->SetTitleOffset(1.0);
  y_tgt->GetXaxis()->SetTitleSize(0.05);
  y_tgt->GetYaxis()->SetTitleSize(0.05);
  y_tgt->GetXaxis()->SetTitle("y_{tg}");
  y_tgt->GetYaxis()->SetTitle("Entries");

  T->Draw("L.tr.tg_y>>y_tgt");

  c6->cd(2);
  TH1D* theta_tgt = new TH1D("theta_tgt","\\theta_{tg}",300,-70,90);

  theta_tgt->GetYaxis()->SetTitleOffset(1.0);
  theta_tgt->GetXaxis()->SetTitleSize(0.05);
  theta_tgt->GetYaxis()->SetTitleSize(0.05);
  theta_tgt->GetXaxis()->SetTitle("\\theta_{tg} [mrad]");
  theta_tgt->GetYaxis()->SetTitle("Entries");

  T->Draw("L.tr.tg_th*3000>>theta_tgt");


  c6->cd(3);
  TH1D* phi_tgt = new TH1D("phi_tgt","\\phi_{tg}",300,-70,60);

  phi_tgt->GetYaxis()->SetTitleOffset(1.0);
  phi_tgt->GetXaxis()->SetTitleSize(0.05);
  phi_tgt->GetYaxis()->SetTitleSize(0.05);
  phi_tgt->GetXaxis()->SetTitle("\\phi_{tg} [mrad]");
  phi_tgt->GetYaxis()->SetTitle("Entries");

  T->Draw("L.tr.tg_ph*1000>>phi_tgt");


  c6->cd(4);
  TH1D* dp_tgt = new TH1D("dp_tgt","\\deltap_{tg}",300,-0.07,0.07);

  dp_tgt->GetYaxis()->SetTitleOffset(1.0);
  dp_tgt->GetXaxis()->SetTitleSize(0.05);
  dp_tgt->GetYaxis()->SetTitleSize(0.05);
  dp_tgt->GetXaxis()->SetTitle("\\deltap_{tg}");
  dp_tgt->GetYaxis()->SetTitle("Entries");

  T->Draw("L.tr.tg_dp>>dp_tgt");
  
  
  
  // add target variables plotted against one another
  
  TCanvas *c7 = new TCanvas("c7","LHRS target variables",1000,800); 
  c7->Divide(2,3); 


  c7->cd(1);
  TH2D* ytgt_v_thtgt = new TH2D("ytgt_v_thtgt","y_{tg} v \\theta_{tg}",300,thtg_l,thtg_h,300,ytg_l,ytg_h);
  
  ytgt_v_thtgt->GetYaxis()->SetTitleOffset(1.0);
  ytgt_v_thtgt->GetXaxis()->SetTitleSize(0.05);
  ytgt_v_thtgt->GetYaxis()->SetTitleSize(0.05);
  ytgt_v_thtgt->GetXaxis()->SetTitle("\\theta_{tg} [mrad]");
  ytgt_v_thtgt->GetYaxis()->SetTitle("y_{tg}");

  
  T->Draw("L.tr.tg_y:L.tr.tg_th*1000>>ytgt_v_thtgt",PID_cuts + GeneralSieveCut,"colz");



  c7->cd(2);
  TH2D* ytgt_v_phtgt = new TH2D("ytgt_v_phtgt","y_{tg} v \\phi_{tg}",300,phtg_l,phtg_h,300,ytg_l,ytg_h);
  
  ytgt_v_phtgt->GetYaxis()->SetTitleOffset(1.0);
  ytgt_v_phtgt->GetXaxis()->SetTitleSize(0.05);
  ytgt_v_phtgt->GetYaxis()->SetTitleSize(0.05);
  ytgt_v_phtgt->GetXaxis()->SetTitle("\\phi_{tg} [mrad]");
  ytgt_v_phtgt->GetYaxis()->SetTitle("y_{tg}");

  
  T->Draw("L.tr.tg_y:L.tr.tg_ph*1000>>ytgt_v_phtgt",PID_cuts + GeneralSieveCut,"colz");


  c7->cd(3);
  TH2D* ytgt_v_dptgt = new TH2D("ytgt_v_dptgt","y_{tg} v \\deltap_{tg}",300,dptg_l,dptg_h,300,ytg_l,ytg_h);
  
  ytgt_v_dptgt->GetYaxis()->SetTitleOffset(1.0);
  ytgt_v_dptgt->GetXaxis()->SetTitleSize(0.05);
  ytgt_v_dptgt->GetYaxis()->SetTitleSize(0.05);
  ytgt_v_dptgt->GetXaxis()->SetTitle("\\deltap_{tg}");
  ytgt_v_dptgt->GetYaxis()->SetTitle("y_{tg}");

  
  T->Draw("L.tr.tg_y:L.tr.tg_dp>>ytgt_v_dptgt",PID_cuts + GeneralSieveCut,"colz");


  c7->cd(4);
  TH2D* thtgt_v_phtgt = new TH2D("thtgt_v_phtgt","\\theta_{tg} v \\phi_{tg}",300,phtg_l,phtg_h,300,thtg_l,thtg_h);


  thtgt_v_phtgt->GetYaxis()->SetTitleOffset(1.0);
  thtgt_v_phtgt->GetXaxis()->SetTitleSize(0.05);
  thtgt_v_phtgt->GetYaxis()->SetTitleSize(0.05);
  thtgt_v_phtgt->GetXaxis()->SetTitle("ph (trgt) [mrad]");
  thtgt_v_phtgt->GetYaxis()->SetTitle("th (trgt) [mrad]");

  
  T->Draw("L.tr.tg_th*1000:L.tr.tg_ph*1000>>thtgt_v_phtgt",PID_cuts + GeneralSieveCut,"colz");




  c7->cd(5);
  

  TH2D* thtgt_v_dptgt = new TH2D("thtgt_v_dptgt","\\theta_{tg} v \\deltap_{tg}",300,dptg_l,dptg_h,300,thtg_l,thtg_h);


  thtgt_v_dptgt->GetYaxis()->SetTitleOffset(1.0);
  thtgt_v_dptgt->GetXaxis()->SetTitleSize(0.05);
  thtgt_v_dptgt->GetYaxis()->SetTitleSize(0.05);
  thtgt_v_dptgt->GetXaxis()->SetTitle("\\deltap_{tg}");
  thtgt_v_dptgt->GetYaxis()->SetTitle("th (trgt) [mrad]");

  
  T->Draw("L.tr.tg_th*1000:L.tr.tg_dp>>thtgt_v_dptgt",PID_cuts + GeneralSieveCut,"colz");




  c7->cd(6);
  

  TH2D* phtgt_v_dptgt = new TH2D("phtgt_v_dptgt","\\phi_{tg} v \\deltap_{tg}",300,dptg_l,dptg_h,300,phtg_l,phtg_h);


  phtgt_v_dptgt->GetYaxis()->SetTitleOffset(1.0);
  phtgt_v_dptgt->GetXaxis()->SetTitleSize(0.05);
  phtgt_v_dptgt->GetYaxis()->SetTitleSize(0.05);
  phtgt_v_dptgt->GetXaxis()->SetTitle("\\deltap_{tg}");
  phtgt_v_dptgt->GetYaxis()->SetTitle("ph (trgt) [mrad]");

  
  T->Draw("L.tr.tg_ph*1000:L.tr.tg_dp>>phtgt_v_dptgt",PID_cuts + GeneralSieveCut,"colz");



  





  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Add focal-plane canvases (plotting focal plane variables against one another)
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  
  TCanvas *c8 = new TCanvas("c8","LHRS FP variables",1000,800); 
  c8->Divide(2,2); 

  c8->cd(1);
  TH1D* x_fp = new TH1D("x_fp","x_{FP}",300,xfp_l,xfp_h);
  
  x_fp->GetYaxis()->SetTitleOffset(1.0);
  x_fp->GetXaxis()->SetTitleSize(0.05);
  x_fp->GetYaxis()->SetTitleSize(0.05);
  x_fp->GetXaxis()->SetTitle("x_{FP}");
  x_fp->GetYaxis()->SetTitle("Entries");

  T->Draw("L.tr.r_x>>x_fp");


  c8->cd(2);
  TH1D* y_fp = new TH1D("y_fp","y_{FP}",300,yfp_l,yfp_h);
  
  y_fp->GetYaxis()->SetTitleOffset(1.0);
  y_fp->GetXaxis()->SetTitleSize(0.05);
  y_fp->GetYaxis()->SetTitleSize(0.05);
  y_fp->GetXaxis()->SetTitle("y_{FP}");
  y_fp->GetYaxis()->SetTitle("Entries");

  T->Draw("L.tr.r_y>>y_fp");


  c8->cd(3);
  TH1D* th_fp = new TH1D("th_fp","\\theta_{FP}",300,thfp_l,thfp_h);
  
  th_fp->GetYaxis()->SetTitleOffset(1.0);
  th_fp->GetXaxis()->SetTitleSize(0.05);
  th_fp->GetYaxis()->SetTitleSize(0.05);
  th_fp->GetXaxis()->SetTitle("\\theta_{FP} [mrad]");
  th_fp->GetYaxis()->SetTitle("Entries");

  T->Draw("L.tr.r_th*1000>>th_fp");


  c8->cd(4);
  TH1D* ph_fp = new TH1D("ph_fp","\\phi_{FP}",300,phfp_l,phfp_h);
  
  ph_fp->GetYaxis()->SetTitleOffset(1.0);
  ph_fp->GetXaxis()->SetTitleSize(0.05);
  ph_fp->GetYaxis()->SetTitleSize(0.05);
  ph_fp->GetXaxis()->SetTitle("\\phi_{FP} [mrad]");
  ph_fp->GetYaxis()->SetTitle("Entries");

  T->Draw("L.tr.r_ph*1000>>ph_fp");
  


  // add FP variables plotted against one another
  
  TCanvas *c9 = new TCanvas("c9","LHRS FP variables",1000,800); 
  c9->Divide(2,3); 


  c9->cd(1);
  TH2D* yfp_v_xfp  = new TH2D("yfp_v_xfp","y_{FP} v x_{FP}",300,xfp_l,xfp_h,300,yfp_l,yfp_h);

  yfp_v_xfp->GetYaxis()->SetTitleOffset(1.0);
  yfp_v_xfp->GetXaxis()->SetTitleSize(0.05);
  yfp_v_xfp->GetYaxis()->SetTitleSize(0.05);
  yfp_v_xfp->GetXaxis()->SetTitle("x_{FP}");
  yfp_v_xfp->GetYaxis()->SetTitle("y_{FP}");

  T->Draw("L.tr.r_y:L.tr.r_x>>yfp_v_xfp",PID_cuts + GeneralSieveCut,"colz");



  c9->cd(2);
  TH2D* yfp_v_thfp  = new TH2D("thfp_v_yfp","\\theta_{FP} v y_{FP}",300,yfp_l,yfp_h,300,thfp_l,thfp_h);

  yfp_v_thfp->GetYaxis()->SetTitleOffset(1.0);
  yfp_v_thfp->GetXaxis()->SetTitleSize(0.05);
  yfp_v_thfp->GetYaxis()->SetTitleSize(0.05);
  yfp_v_thfp->GetXaxis()->SetTitle("y_{FP}");
  yfp_v_thfp->GetYaxis()->SetTitle("\\theta_{FP} [mrad]");

  T->Draw("L.tr.r_th*1000:L.tr.r_y>>thfp_v_yfp",PID_cuts + GeneralSieveCut,"colz");


  c9->cd(3);
  TH2D* yfp_v_phfp  = new TH2D("phfp_v_yfp","\\phi_{FP} v y_{FP}",300,yfp_l,yfp_h,300,phfp_l,phfp_h);

  yfp_v_phfp->GetYaxis()->SetTitleOffset(1.0);
  yfp_v_phfp->GetXaxis()->SetTitleSize(0.05);
  yfp_v_phfp->GetYaxis()->SetTitleSize(0.05);
  yfp_v_phfp->GetXaxis()->SetTitle("y_{FP}");
  yfp_v_phfp->GetYaxis()->SetTitle("\\phi_{FP} [mrad]");

  T->Draw("L.tr.r_ph*1000:L.tr.r_y>>phfp_v_yfp",PID_cuts + GeneralSieveCut,"colz");


  c9->cd(4);
  TH2D* xfp_v_thfp  = new TH2D("xfp_v_thfp","x_{FP} v \\theta_{FP}",300,thfp_l,thfp_h,300,xfp_l,xfp_h);

  xfp_v_thfp->GetYaxis()->SetTitleOffset(1.0);
  xfp_v_thfp->GetXaxis()->SetTitleSize(0.05);
  xfp_v_thfp->GetYaxis()->SetTitleSize(0.05);
  xfp_v_thfp->GetXaxis()->SetTitle("\\theta_{FP} [mrad]");
  xfp_v_thfp->GetYaxis()->SetTitle("x_{FP}");

  T->Draw("L.tr.r_x:L.tr.r_th*1000>>xfp_v_thfp",PID_cuts + GeneralSieveCut,"colz");


  c9->cd(5);
  TH2D* xfp_v_phfp  = new TH2D("xfp_v_phfp","x_{FP} v \\phi_{FP}",300,phfp_l,phfp_h,300,xfp_l,xfp_h);

  xfp_v_phfp->GetYaxis()->SetTitleOffset(1.0);
  xfp_v_phfp->GetXaxis()->SetTitleSize(0.05);
  xfp_v_phfp->GetYaxis()->SetTitleSize(0.05);
  xfp_v_phfp->GetXaxis()->SetTitle("\\phi_{FP} [mrad]");
  xfp_v_phfp->GetYaxis()->SetTitle("x_{FP}");

  T->Draw("L.tr.r_x:L.tr.r_ph*1000>>xfp_v_phfp",PID_cuts + GeneralSieveCut,"colz");


  c9->cd(6);
  TH2D* thfp_v_phfp  = new TH2D("thfp_v_phfp","\\theta_{FP} v \\phi_{FP}",300,phfp_l,phfp_h,300,thfp_l,thfp_h);

  thfp_v_phfp->GetYaxis()->SetTitleOffset(1.0);
  thfp_v_phfp->GetXaxis()->SetTitleSize(0.05);
  thfp_v_phfp->GetYaxis()->SetTitleSize(0.05);
  thfp_v_phfp->GetXaxis()->SetTitle("\\phi_{FP} [mrad]");
  thfp_v_phfp->GetYaxis()->SetTitle("\\theta_{FP} [mrad]");

  T->Draw("L.tr.r_th*1000:L.tr.r_ph*1000>>thfp_v_phfp",PID_cuts + GeneralSieveCut,"colz");



  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Add ReactZ plotted against target variables
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  TCanvas *c10 = new TCanvas("c10","LHRS ZReact vs target variables",1000,800); 
  c10->Divide(2,2); 

  c10->cd(1);

  TH2D* ReactZ_vs_Thtgt = new TH2D("ReactZ_vs_Thtgt","ReactZ v \\theta_{tg}",300,thtg_l,thtg_h,300,zreact_l,zreact_h);

  ReactZ_vs_Thtgt->GetYaxis()->SetTitleOffset(1.0);
  ReactZ_vs_Thtgt->GetXaxis()->SetTitleSize(0.05);
  ReactZ_vs_Thtgt->GetYaxis()->SetTitleSize(0.05);
  ReactZ_vs_Thtgt->GetXaxis()->SetTitle("\\theta_{tg}");
  ReactZ_vs_Thtgt->GetYaxis()->SetTitle("ReactZ");

  
  T->Draw("reactz:L.tr.tg_th*1000>>ReactZ_vs_Thtgt",PID_cuts + GeneralSieveCut,"colz");

  c10->cd(2);
  //  TH2D* ReactZ_vs_Phtgt = new TH2D("ReactZ_vs_Phtgt","ReactZ v \\phi_{tg}",300,zreact_l,zreact_h,300,phtg_l,phtg_h);
  TH2D* ReactZ_vs_Phtgt = new TH2D("ReactZ_vs_Phtgt","ReactZ v \\phi_{tg}",300,phtg_l,phtg_h,300,zreact_l,zreact_h);

  ReactZ_vs_Phtgt->GetYaxis()->SetTitleOffset(1.0);
  ReactZ_vs_Phtgt->GetXaxis()->SetTitleSize(0.05);
  ReactZ_vs_Phtgt->GetYaxis()->SetTitleSize(0.05);
  ReactZ_vs_Phtgt->GetXaxis()->SetTitle("\\phi_{tg}");
  ReactZ_vs_Phtgt->GetYaxis()->SetTitle("ReactZ");

  
  T->Draw("reactz:L.tr.tg_ph*1000>>ReactZ_vs_Phtgt",PID_cuts + GeneralSieveCut,"colz");


  c10->cd(3);
  //  TH2D* ReactZ_vs_ytgt = new TH2D("ReactZ_vs_ytgt","ReactZ v y_{tg}",300,zreact_l,zreact_h,300,ytg_l,ytg_h);
  TH2D* ReactZ_vs_ytgt = new TH2D("ReactZ_vs_ytgt","ReactZ v y_{tg}",300,ytg_l,ytg_h,300,zreact_l,zreact_h);

  ReactZ_vs_ytgt->GetYaxis()->SetTitleOffset(1.0);
  ReactZ_vs_ytgt->GetXaxis()->SetTitleSize(0.05);
  ReactZ_vs_ytgt->GetYaxis()->SetTitleSize(0.05);
  ReactZ_vs_ytgt->GetXaxis()->SetTitle("y_{tg}");
  ReactZ_vs_ytgt->GetYaxis()->SetTitle("ReactZ");

  
  T->Draw("reactz:L.tr.tg_y>>ReactZ_vs_ytgt",PID_cuts + GeneralSieveCut,"colz");



  c10->cd(4);
  //  TH2D* ReactZ_vs_dptgt = new TH2D("ReactZ_vs_dptgt","ReactZ v \\deltap_{tg}",300,zreact_l,zreact_h,300,dptg_l,dptg_h);
  TH2D* ReactZ_vs_dptgt = new TH2D("ReactZ_vs_dptgt","ReactZ v \\deltap_{tg}",300,dptg_l,dptg_h,300,zreact_l,zreact_h);

  ReactZ_vs_dptgt->GetYaxis()->SetTitleOffset(1.0);
  ReactZ_vs_dptgt->GetXaxis()->SetTitleSize(0.05);
  ReactZ_vs_dptgt->GetYaxis()->SetTitleSize(0.05);
  ReactZ_vs_dptgt->GetXaxis()->SetTitle("\\deltap_{tg}");
  ReactZ_vs_dptgt->GetYaxis()->SetTitle("ReactZ");

  
  T->Draw("reactz:L.tr.tg_dp>>ReactZ_vs_dptgt",PID_cuts + GeneralSieveCut,"colz");

  





  
  // create suffix for pdf file name based on cuts used

  TString pdf_suffix;



  if( row == -1 && col == -1){
    // case for all holes/ data being plotted
    pdf_suffix = "_all_holes";

  }
  else if(row >= 0 && col < 0){
    // case for row being selected 
    pdf_suffix = Form("_row_%i",row);

  }
  else if(row < 0 && col >= 0){
    // case for column being selected 
    pdf_suffix = Form("_col_%i",col);
  }
  else if(row == 1 && col == 1){
    // case for hole being selected 
    pdf_suffix = "_all_events";
  }
  else if(row >= 0 && col >= 0){
    // case for hole being selected 
    pdf_suffix = Form("_hole_row_%i_col_%i",row,col);
  }

  // Print canvas

  cout << "pdf_suffix = " << pdf_suffix << endl;
  



  c1->Print("rootfiles/" + DB_name + Form("/Correlation_plots_%i",runnumber) + pdf_suffix + ".pdf(");
  c2->Print("rootfiles/" + DB_name + Form("/Correlation_plots_%i",runnumber) + pdf_suffix + ".pdf");
  c3->Print("rootfiles/" + DB_name + Form("/Correlation_plots_%i",runnumber) + pdf_suffix + ".pdf");
  c4->Print("rootfiles/" + DB_name + Form("/Correlation_plots_%i",runnumber) + pdf_suffix + ".pdf");
  c6->Print("rootfiles/" + DB_name + Form("/Correlation_plots_%i",runnumber) + pdf_suffix + ".pdf");
  c7->Print("rootfiles/" + DB_name + Form("/Correlation_plots_%i",runnumber) + pdf_suffix + ".pdf");
  c8->Print("rootfiles/" + DB_name + Form("/Correlation_plots_%i",runnumber) + pdf_suffix + ".pdf");
  c9->Print("rootfiles/" + DB_name + Form("/Correlation_plots_%i",runnumber) + pdf_suffix + ".pdf");
  c10->Print("rootfiles/" + DB_name + Form("/Correlation_plots_%i",runnumber) + pdf_suffix + ".pdf)");


  gSystem->Exec("cp " + string(gSystem->pwd()) + "/rootfiles/" + DB_name + Form("/Correlation_plots_%i" + pdf_suffix + ".pdf",runnumber) + " /home/johnw/public_html/correlation_plots/latest.pdf");


 // gSystem->Exec("ln -fs " + string(gSystem->pwd()) + "/rootfiles/" + file_date + Form("/Correlation_plots_%i" + pdf_suffix + ".pdf",runnumber) + " /home/johnw/public_html/correlation_plots/latest.pdf");




  // copy pdf to specific location dependent on cuts used


  gSystem->Exec("mkdir /home/johnw/public_html/correlation_plots/" + DB_name + "/" + Form("%d",runnumber));



  if( row == -1 && col == -1){
    // case for all holes/ data being plotted


    gSystem->Exec("mkdir /home/johnw/public_html/correlation_plots/" + DB_name + "/" + Form("%d",runnumber) + "/all_hole_plots");

    gSystem->Exec("cp " + string(gSystem->pwd()) + "/rootfiles/" + DB_name + Form("/Correlation_plots_%i" + pdf_suffix + ".pdf",runnumber) + " /home/johnw/public_html/correlation_plots/" + DB_name + "/" + Form("%d",runnumber) + "/all_hole_plots" + Form("/Correlation_plots_%i" + pdf_suffix + ".pdf",runnumber));
    
  }
  else if(row >= 0 && col < 0){
    // case for row being selected 

    gSystem->Exec("mkdir /home/johnw/public_html/correlation_plots/" + DB_name + "/" + Form("%d",runnumber) + "/row_plots");

    gSystem->Exec("cp " + string(gSystem->pwd()) + "/rootfiles/" + DB_name + Form("/Correlation_plots_%i" + pdf_suffix + ".pdf",runnumber) + " /home/johnw/public_html/correlation_plots/" + DB_name + "/" + Form("%d",runnumber) +  "/row_plots" + Form("/Correlation_plots_%i" + pdf_suffix + ".pdf",runnumber));
    
  }
  else if(row < 0 && col >= 0){
    // case for column being selected 
    gSystem->Exec("mkdir /home/johnw/public_html/correlation_plots/" + DB_name + "/" + Form("%d",runnumber) + "/column_plots");

    gSystem->Exec("cp " + string(gSystem->pwd()) + "/rootfiles/" + DB_name + Form("/Correlation_plots_%i" + pdf_suffix + ".pdf",runnumber) + " /home/johnw/public_html/correlation_plots/" + DB_name + "/" + Form("%d",runnumber) + "/column_plots" + Form("/Correlation_plots_%i" + pdf_suffix + ".pdf",runnumber));

    cout << "Passed column plots " << endl;

  }
  else if(row == 1 && col == 1){

    gSystem->Exec("mkdir /home/johnw/public_html/correlation_plots/" + DB_name + "/" + Form("%d",runnumber) + "/all_events_plots");
    
    
    gSystem->Exec("cp " + string(gSystem->pwd()) + "/rootfiles/" + DB_name + Form("/Correlation_plots_%i" + pdf_suffix + ".pdf",runnumber) + " /home/johnw/public_html/correlation_plots/" + DB_name + "/" + Form("%d",runnumber) + "/all_events_plots" + Form("/Correlation_plots_%i" + pdf_suffix + ".pdf",runnumber));

  }
  else if(row >= 0 && col >= 0){
    // case for hole being selected 
    gSystem->Exec("mkdir /home/johnw/public_html/correlation_plots/" + DB_name + "/" + Form("%d",runnumber) + "/hole_plots");

    gSystem->Exec("cp " + string(gSystem->pwd()) + "/rootfiles/" + DB_name + Form("/Correlation_plots_%i" + pdf_suffix + ".pdf",runnumber) + " /home/johnw/public_html/correlation_plots/" + DB_name + "/" + Form("%d",runnumber) + "/hole_plots" + Form("/Correlation_plots_%i" + pdf_suffix + ".pdf",runnumber));
    

  }


  
  



}



void corr_plot(Run_spec run_info, TChain* T){



  TString DB_name = run_info.DB_info[DB_NAME];
  TString Cut_name = run_info.DB_info[CUT_NAME];
  
  Int_t runnumber = run_info.event_info[RUN_NO];
  Int_t foil_no = run_info.event_info[FOIL_NO];
  Int_t col = run_info.event_info[COL_NO];
  Int_t row = run_info.event_info[ROW_NO];



  corr_plot(DB_name, T, runnumber, col, row);

}
