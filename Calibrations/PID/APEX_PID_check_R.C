// John Williamson 6/10/2019
// APEX PID script designed to find correct level of cut in Cherenkov and calorimeters for maximissing ratio of electrons to mesons
// 


#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <stdexcept>
#include <cassert>

#include "TTree.h"
//#include "TFile.h"
#include "TString.h"
#include "TStyle.h"
#include "TLine.h"
#include "TBox.h"
#include "Load_more_rootfiles.C"



using namespace std;


Double_t peak(Double_t *x, Double_t *par) {
  return  par[0]*TMath::Exp(-(((x[0]-par[1])*(x[0]-par[1]))/(2*par[2]*par[2])));
}

Double_t bg(Double_t *x, Double_t *par) {
  return par[0]*TMath::Exp(par[1]*x[0]);
}

// sum of peaks and backgorund
Double_t overall(Double_t *x, Double_t *par) {
  return peak(x,par) + peak(x,&par[3]) + bg(x,&par[6]);
}

void APEX_PID_check_R()
{


  



  Int_t run_number;    

  cout << "enter run_number: ";
  cin >> run_number;

  
  TChain* T_2 = Load_more_rootfiles(run_number);
  // TChain* T = Load_more_rootfiles(run_number);
  //
  TTree* T = T_2->CloneTree(1e4);
  





  cout << T->GetEntries() << endl;


  



 
  // Right-arm PID plots


  // right cut 

  TCut rcut = "R.tr.n==1";
  gStyle->SetOptStat(0);







  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~1
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~1

  // Initial 'good' electron, pion and other backgorund cuts
  // where e_cuts are electron cuts, bg_cuts are background cuts (with signal in cherenkov)
  // pi_cuts are pion cuts and zer_cuts are cuts with zero signal in both detectors
  // EP cuts are calorimeter cuts based on total energy in both layers
  // _h are higher bound cuts, _l are lower bound cuts

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~1
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~1



  TCanvas *c1 = new TCanvas("c1","Track variables and relation to electrons, mesons, bg",1200,1200);

  c1->Divide(1,2);


  c1->cd(1);


  
  // Double_t e_cut_cal_h = 1.05;
  // Double_t e_cut_cal_l = 0.75;
  
  // Double_t e_cut_cer_h  = 8500;
  // Double_t e_cut_cer_l  = 2000;


  // Double_t pi_cut_cal_h = 0.33;
  // Double_t pi_cut_cal_l = 0.165;
  
  // Double_t pi_cut_cer_h  = 200;
  // Double_t pi_cut_cer_l  = 0.0;


  // Double_t zer_cut_cal_h = 0.13;
  // Double_t zer_cut_cal_l = -0.03;
  
  // Double_t zer_cut_cer_h  = 200;
  // Double_t zer_cut_cer_l  = 0.0;


  // Double_t bg_cut_cal_h = 0.15;
  // Double_t bg_cut_cal_l = -0.03;
  
  // Double_t bg_cut_cer_h  = 6200;
  // Double_t bg_cut_cer_l  = 1750;


  Double_t e_cut_cal_h = 1.12;
  Double_t e_cut_cal_l = 0.87;
  
  Double_t e_cut_cer_h  = 7000;
  Double_t e_cut_cer_l  = 2500;


  Double_t pi_cut_cal_h = 0.30;
  Double_t pi_cut_cal_l = 0.165;
  
  Double_t pi_cut_cer_h  = 200;
  Double_t pi_cut_cer_l  = 0.0;


  Double_t zer_cut_cal_h = 0.13;
  Double_t zer_cut_cal_l = -0.03;
  
  Double_t zer_cut_cer_h  = 200;
  Double_t zer_cut_cer_l  = 0.0;


  Double_t bg_cut_cal_h = 0.13;
  Double_t bg_cut_cal_l = -0.03;
  
  Double_t bg_cut_cer_h  = 6200;
  Double_t bg_cut_cer_l  = 1750;


  

  // Cherenkov sum vs Calorimeter energy plot

  TH2D *ChVsEP = new TH2D("ChVsEP","RHRS PID Plot",1000,-0.1,1.4,1000,-10,10000);

  

  ChVsEP->GetXaxis()->SetTitle("E_{tot}/P");
  ChVsEP->GetXaxis()->CenterTitle();
  ChVsEP->GetXaxis()->SetTitleSize(0.048);
  ChVsEP->GetXaxis()->SetTitleOffset(0.98);
  
  
  
  ChVsEP->GetYaxis()->SetTitle("Cherenkov sum");
  ChVsEP->GetYaxis()->CenterTitle();
  ChVsEP->GetYaxis()->SetTitleSize(0.048);
  ChVsEP->GetYaxis()->SetTitleOffset(0.85);


  
  T->Draw("R.cer.asum_c:((R.ps.e+R.sh.e)/(R.gold.p*1000))>>ChVsEP",rcut,"col");


  // box defining electrons

  TBox* e_cut_box = new TBox(e_cut_cal_l,e_cut_cer_l,e_cut_cal_h,e_cut_cer_h);

  e_cut_box->SetLineColor(kRed);
  e_cut_box->SetFillStyle(0);
  e_cut_box->SetLineWidth(3);
  

  e_cut_box->Draw("same");



  // box defining pions

  TBox* pi_cut_box = new TBox(pi_cut_cal_l,pi_cut_cer_l,pi_cut_cal_h,pi_cut_cer_h);

  pi_cut_box->SetLineColor(kBlue);
  pi_cut_box->SetFillStyle(0);
  pi_cut_box->SetLineWidth(2);

  pi_cut_box->DrawClone("same");


  // box defining zero hit events

  TBox* zer_cut_box = new TBox(zer_cut_cal_l,zer_cut_cer_l,zer_cut_cal_h,zer_cut_cer_h);

  zer_cut_box->SetLineColor(kBlack);
  zer_cut_box->SetFillStyle(0);
  zer_cut_box->SetLineWidth(2);

  zer_cut_box->DrawClone("same");


  // box defining cherenkov background

  TBox* bg_cut_box = new TBox(bg_cut_cal_l,bg_cut_cer_l,bg_cut_cal_h,bg_cut_cer_h);

  bg_cut_box->SetLineColor(kGreen);
  bg_cut_box->SetFillStyle(0);
  bg_cut_box->SetLineWidth(3);
  bg_cut_box->Draw("same");



  // reset line widths for purpose of legend
  pi_cut_box->SetLineWidth(3);
  zer_cut_box->SetLineWidth(3);

  
  // legend for cherenkov sum vs Calorimeter energy plot

  TLegend* leg = new TLegend(.1,.65,.37,.9,"Key");
  leg->SetFillColor(0);
  leg->AddEntry(e_cut_box,"Electrons","l");
  leg->AddEntry(pi_cut_box,"Pions","l");
  leg->AddEntry(zer_cut_box,"No hit","l");
  leg->AddEntry(bg_cut_box,"Background","l");

  leg->Draw("same");




  // create plot 'zoomed' in to focus on position of smaller pion and no-hit boxes on Cherenkov sum vs Calorimeter enrgy plot


  c1->cd(2);

  TH2D* ChVsEP_zoom = (TH2D*)ChVsEP->Clone("ChVsEP_zoom");

  ChVsEP_zoom->GetXaxis()->SetRangeUser(-0.1,0.40);
  ChVsEP_zoom->GetYaxis()->SetRangeUser(-10,250);

  ChVsEP_zoom->Draw("col");  

  
  pi_cut_box->SetLineWidth(6);
  pi_cut_box->DrawClone("same");

  zer_cut_box->SetLineWidth(6);
  zer_cut_box->DrawClone("same");



  
  // reset line widths for purpose of legend
  pi_cut_box->SetLineWidth(3);
  zer_cut_box->SetLineWidth(3);



  // method of saving aimed to produce better resolution images

  c1->SaveAs("plots/L_event_selection.pdf");
  
  
  gSystem->Exec("convert -density 700 -trim plots/L_event_selection.pdf plots/L_event_selection.png");
  



  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~2
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~2

  // Plot different types of events selected in first canvas against various track properties (test variation of these properties for different species)
  // plot here against target variables: theta,phi,y and dp

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~2
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~2


  // set-up R.gold.y (y_tg of golden track) plots for all 4 categories of events


  TCanvas* c2 = new TCanvas("c2","Track and Reconstruction variables",1200,1200);

  c2->Divide(2,2);
  
  // set-up R.gold.th (theta angle of golden track) plots for all 4 categories of events

  c2->cd(1);
  gPad->SetLogy();

  TH1D* Theta_plot_e =  new TH1D("Theta_plot_e","",100,-0.2,0.2);

  Theta_plot_e->GetXaxis()->SetTitle("\\theta_{tg} [rad]");
  Theta_plot_e->GetXaxis()->CenterTitle();
  Theta_plot_e->GetYaxis()->SetTitle("Entries");
  Theta_plot_e->GetYaxis()->CenterTitle();
  
  
  Theta_plot_e->SetLineColor(kRed);
  Theta_plot_e->SetLineWidth(2);


  // electron cut
  TCut e_cut = Form("R.cer.asum_c > %f &&  R.cer.asum_c < %f && ((R.ps.e+R.sh.e)/(R.gold.p*1000)) > %f && ((R.ps.e+R.sh.e)/(R.gold.p*1000)) < %f",e_cut_cer_l,e_cut_cer_h,e_cut_cal_l,e_cut_cal_h);
  
  
  T->Draw("R.gold.th>>Theta_plot_e",e_cut);


  // draw theta cuts (used for 'good' tracks)


  Double_t theta_cut_l = -0.05;
  Double_t theta_cut_h = 0.06;

  TLine* theta_cut_line_l = new TLine(theta_cut_l,0,theta_cut_l,Theta_plot_e->GetMaximum());

  TLine* theta_cut_line_h = new TLine(theta_cut_h,0,theta_cut_h,Theta_plot_e->GetMaximum());


  theta_cut_line_l->SetLineColor(kMagenta);  
  theta_cut_line_l->SetLineStyle(9);
  theta_cut_line_l->SetLineWidth(2);
  theta_cut_line_l->Draw("same");

  theta_cut_line_h->SetLineColor(kMagenta);  
  theta_cut_line_h->SetLineStyle(9);
  theta_cut_line_h->SetLineWidth(2);
  theta_cut_line_h->Draw("same");



  //pion theta plot

  TH1D* Theta_plot_pi = (TH1D*)Theta_plot_e->Clone("Theta_plot_pi");
  Theta_plot_pi->SetLineColor(kBlue);

  TCut pi_cut = Form("R.cer.asum_c > %f &&  R.cer.asum_c < %f && ((R.ps.e+R.sh.e)/(R.gold.p*1000)) > %f && ((R.ps.e+R.sh.e)/(R.gold.p*1000)) < %f",pi_cut_cer_l,pi_cut_cer_h,pi_cut_cal_l,pi_cut_cal_h);
  
  T->Draw("R.gold.th>>Theta_plot_pi",pi_cut,"same");



  // no hits theta plot

  TH1D* Theta_plot_zer = (TH1D*)Theta_plot_e->Clone("Theta_plot_zer");

  TCut zer_cut = Form("R.cer.asum_c > %f &&  R.cer.asum_c < %f && ((R.ps.e+R.sh.e)/(R.gold.p*1000)) > %f && ((R.ps.e+R.sh.e)/(R.gold.p*1000)) < %f",zer_cut_cer_l,zer_cut_cer_h,zer_cut_cal_l,zer_cut_cal_h);

  Theta_plot_zer->SetLineColor(kBlack);
  T->Draw("R.gold.th>>Theta_plot_zer",zer_cut,"same");


  // background theta plot

  TH1D* Theta_plot_bg = (TH1D*)Theta_plot_e->Clone("Theta_plot_bg");


  TCut bg_cut = Form("R.cer.asum_c > %f &&  R.cer.asum_c < %f && ((R.ps.e+R.sh.e)/(R.gold.p*1000)) > %f && ((R.ps.e+R.sh.e)/(R.gold.p*1000)) < %f",bg_cut_cer_l,bg_cut_cer_h,bg_cut_cal_l,bg_cut_cal_h);

  
  Theta_plot_bg->SetLineColor(kGreen);
  T->Draw("R.gold.th>>Theta_plot_bg",bg_cut,"same");


  cout << "Number of Entries in \\theta plots: " << endl;
  cout << "Electrons " << Theta_plot_e->GetEntries() << endl;
  cout << "Pions " << Theta_plot_pi->GetEntries() << endl;
  cout << "Zeros " << Theta_plot_zer->GetEntries() << endl;
  cout << "Bg " << Theta_plot_bg->GetEntries() << endl;
  



  // set-up R.gold.ph (phi angle of golden track) plots for all 4 categories of events

  c2->cd(2);
  gPad->SetLogy();

  TH1D* Phi_plot_e =  new TH1D("Phi_plot_e","",100,-0.2,0.2);

  Phi_plot_e->GetXaxis()->SetTitle("\\phi_{tg} [rad]");
  Phi_plot_e->GetXaxis()->CenterTitle();
  Phi_plot_e->GetYaxis()->SetTitle("Entries");
  Phi_plot_e->GetYaxis()->CenterTitle();
  
  
  Phi_plot_e->SetLineColor(kRed);
  Phi_plot_e->SetLineWidth(2);;

  T->Draw("R.gold.ph>>Phi_plot_e",e_cut);




  // phi cut (used later on to select 'good' electrons)

  Double_t phi_cut_l = -0.045;
  Double_t phi_cut_h = 0.055;

  TLine* phi_cut_line_l = new TLine(phi_cut_l,0,phi_cut_l,Phi_plot_e->GetMaximum());

  TLine* phi_cut_line_h = new TLine(phi_cut_h,0,phi_cut_h,Phi_plot_e->GetMaximum());


  phi_cut_line_l->SetLineColor(kMagenta);  
  phi_cut_line_l->SetLineStyle(9);
  phi_cut_line_l->SetLineWidth(2);
  phi_cut_line_l->Draw("same");

  phi_cut_line_h->SetLineColor(kMagenta);  
  phi_cut_line_h->SetLineStyle(9);
  phi_cut_line_h->SetLineWidth(2);
  phi_cut_line_h->Draw("same");



  //pion phi plot

  TH1D* Phi_plot_pi = (TH1D*)Phi_plot_e->Clone("Phi_plot_pi");
  Phi_plot_pi->SetLineColor(kBlue);
  T->Draw("R.gold.ph>>Phi_plot_pi",pi_cut,"same");



  // no hits phi plot

  TH1D* Phi_plot_zer = (TH1D*)Phi_plot_e->Clone("Phi_plot_zer");


  Phi_plot_zer->SetLineColor(kBlack);
  T->Draw("R.gold.ph>>Phi_plot_zer",zer_cut,"same");


  // background phi plot

  TH1D* Phi_plot_bg = (TH1D*)Phi_plot_e->Clone("Phi_plot_bg");

  Phi_plot_bg->SetLineColor(kGreen);
  T->Draw("R.gold.ph>>Phi_plot_bg",bg_cut,"same");




  // set-up R.gold.dp (deviation from central momentum of golden track) plots for all 4 categories of events

  c2->cd(3);
  gPad->SetLogy();

  TH1D* Dp_plot_e =  new TH1D("Dp_plot_e","",100,-0.2,0.2);

  Dp_plot_e->GetXaxis()->SetTitle("\\delta p");
  Dp_plot_e->GetXaxis()->CenterTitle();
  Dp_plot_e->GetYaxis()->SetTitle("Entries");
  Dp_plot_e->GetYaxis()->CenterTitle();
  
  
  Dp_plot_e->SetLineColor(kRed);
  Dp_plot_e->SetLineWidth(2);
  
  T->Draw("R.gold.dp>>Dp_plot_e",e_cut);

  //  leg_track->Draw("same");



  // dp cut (used later on to select 'good' electrons

  Double_t dp_cut_l = -0.048;
  Double_t dp_cut_h = 0.052;

  TLine* dp_cut_line_l = new TLine(dp_cut_l,0,dp_cut_l,Dp_plot_e->GetMaximum());

  TLine* dp_cut_line_h = new TLine(dp_cut_h,0,dp_cut_h,Dp_plot_e->GetMaximum());


  dp_cut_line_l->SetLineColor(kMagenta);  
  dp_cut_line_l->SetLineStyle(9);
  dp_cut_line_l->SetLineWidth(2);
  dp_cut_line_l->Draw("same");

  dp_cut_line_h->SetLineColor(kMagenta);  
  dp_cut_line_h->SetLineStyle(9);
  dp_cut_line_h->SetLineWidth(2);
  dp_cut_line_h->Draw("same");




  //pion dp plot

  TH1D* Dp_plot_pi = (TH1D*)Dp_plot_e->Clone("Dp_plot_pi");
  Dp_plot_pi->SetLineColor(kBlue);
  T->Draw("R.gold.dp>>Dp_plot_pi",pi_cut,"same");



  // no hits dp plot

  TH1D* Dp_plot_zer = (TH1D*)Dp_plot_e->Clone("Dp_plot_zer");


  Dp_plot_zer->SetLineColor(kBlack);
  T->Draw("R.gold.dp>>Dp_plot_zer",zer_cut,"same");


  // background dp plot

  TH1D* Dp_plot_bg = (TH1D*)Dp_plot_e->Clone("Dp_plot_bg");

  Dp_plot_bg->SetLineColor(kGreen);
  T->Draw("R.gold.dp>>Dp_plot_bg",bg_cut,"same");




  c2->cd(4);
  gPad->SetLogy();

  TH1D* Y_plot_e =  new TH1D("Y_plot_e","",100,-0.2,0.2);

  Y_plot_e->GetXaxis()->SetTitle("y_{tg} [rad]");
  Y_plot_e->GetXaxis()->CenterTitle();
  Y_plot_e->GetYaxis()->SetTitle("Entries");
  Y_plot_e->GetYaxis()->CenterTitle();
  
  
  Y_plot_e->SetLineColor(kRed);
  Y_plot_e->SetLineWidth(2);

  T->Draw("R.gold.y>>Y_plot_e",e_cut);


  // draw theta cuts (used for 'good' tracks)


  Double_t y_cut_l = -0.05;
  Double_t y_cut_h = 0.05;

  TLine* y_cut_line_l = new TLine(y_cut_l,0,y_cut_l,Y_plot_e->GetMaximum());

  TLine* y_cut_line_h = new TLine(y_cut_h,0,y_cut_h,Y_plot_e->GetMaximum());


  y_cut_line_l->SetLineColor(kMagenta);  
  y_cut_line_l->SetLineStyle(9);
  y_cut_line_l->SetLineWidth(2);
  y_cut_line_l->Draw("same");

  y_cut_line_h->SetLineColor(kMagenta);  
  y_cut_line_h->SetLineStyle(9);
  y_cut_line_h->SetLineWidth(2);
  y_cut_line_h->Draw("same");



  //pion y plot

  TH1D* Y_plot_pi = (TH1D*)Y_plot_e->Clone("Y_plot_pi");
  Y_plot_pi->SetLineColor(kBlue);
  T->Draw("R.gold.y>>Y_plot_pi",pi_cut,"same");



  // no hits y plot

  TH1D* Y_plot_zer = (TH1D*)Y_plot_e->Clone("Y_plot_zer");


  Y_plot_zer->SetLineColor(kBlack);
  T->Draw("R.gold.y>>Y_plot_zer",zer_cut,"same");


  // background y plot

  TH1D* Y_plot_bg = (TH1D*)Y_plot_e->Clone("Y_plot_bg");

  Y_plot_bg->SetLineColor(kGreen);
  T->Draw("R.gold.y>>Y_plot_bg",bg_cut,"same");


  cout << "Number of Entries in \\y plots: " << endl;
  cout << "Electrons " << Y_plot_e->GetEntries() << endl;
  cout << "Pions " << Y_plot_pi->GetEntries() << endl;
  cout << "Zeros " << Y_plot_zer->GetEntries() << endl;
  cout << "Bg " << Y_plot_bg->GetEntries() << endl;
  


  




  

  // save canvas

  c2->SaveAs("plots/L_track_var_plots.pdf");
   
  
  gSystem->Exec("convert -density 700 -trim plots/L_track_var_plots.pdf plots/L_track_var_plots.png");
  




  // for remianing plots change default legend and title sizes 

  gStyle->SetTitleFontSize(0.058);
  gStyle->SetLegendTextSize(0.038);




  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~2_b
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~2_b

  // Plot track properties: Beta and chi^2

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~2_b

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~2_b


  

  TCanvas* c2_b = new TCanvas("c2_b","Track and Reconstruction variables",1200,1200);

  c2_b->Divide(2);
  c2_b->cd(1);
  gPad->SetLogy();


  // set-up Beta plots for 4 categories of events described above: electrons, pions, no_signal and cherenkov-only signal

  TH1D *Beta_plot_e = new TH1D("Beta_plot_e","",100,-3,6);
  Beta_plot_e->GetXaxis()->SetTitle("\\beta");
  Beta_plot_e->GetXaxis()->CenterTitle();
  Beta_plot_e->GetYaxis()->SetTitle("Entries");
  Beta_plot_e->GetYaxis()->CenterTitle();

 
  Beta_plot_e->SetLineColor(kRed);
  Beta_plot_e->SetLineWidth(2);
  
  
  T->Draw("R.gold.beta>>Beta_plot_e",e_cut);


 

  // beta cuts (for later plots) defining 'good electron'

  Double_t beta_cut_l = 0.4;

  TLine* beta_cut_line = new TLine(beta_cut_l,0,beta_cut_l,Beta_plot_e->GetMaximum());


  beta_cut_line->SetLineColor(kMagenta);  
  beta_cut_line->SetLineStyle(9);
  beta_cut_line->SetLineWidth(2);
  beta_cut_line->Draw("same");

 


  TLegend* leg_track = (TLegend*)leg->Clone("leg_track");
  
  leg_track->AddEntry(beta_cut_line,"Cut for 'good' electrons","l");

 


  //pion beta plot

  TH1D* Beta_plot_pi = (TH1D*)Beta_plot_e->Clone("Beta_plot_pi");


  Beta_plot_pi->SetLineColor(kBlue);
  T->Draw("R.gold.beta>>Beta_plot_pi",pi_cut,"same");


  // no hits beta plot

  TH1D* Beta_plot_zer = (TH1D*)Beta_plot_e->Clone("Beta_plot_zer");


  Beta_plot_zer->SetLineColor(kBlack);
  T->Draw("R.gold.beta>>Beta_plot_zer",zer_cut,"same");


  // background beta plot

  TH1D* Beta_plot_bg = (TH1D*)Beta_plot_e->Clone("Beta_plot_bg");


  Beta_plot_bg->SetLineColor(kGreen);
  T->Draw("R.gold.beta>>Beta_plot_bg",bg_cut,"same");




  c2_b->cd(2);
  gPad->SetLogy();


  // set-up Chi^2 plots for 4 categories of events described above: electrons, pions, no_signal and cherenkov-only signal

  TH1D *Chi_plot_e = new TH1D("Chi_plot_e","",100,0,0.01);
  Chi_plot_e->GetXaxis()->SetTitle("\\chi^2");
  Chi_plot_e->GetXaxis()->CenterTitle();
  Chi_plot_e->GetYaxis()->SetTitle("Entries");
  Chi_plot_e->GetYaxis()->CenterTitle();

 
  Chi_plot_e->SetLineColor(kRed);
  Chi_plot_e->SetLineWidth(2);
  
  
  T->Draw("R.tr.chi2>>Chi_plot_e",e_cut);




  
  // chi^2 cuts (for later plots) defining 'good electron'

  Double_t chi2_cut_l = 0.0035;

  TLine* chi2_cut_line = new TLine(chi2_cut_l,0,chi2_cut_l,Chi_plot_e->GetMaximum());


  chi2_cut_line->SetLineColor(kMagenta);  
  chi2_cut_line->SetLineStyle(9);
  chi2_cut_line->SetLineWidth(2);
  chi2_cut_line->Draw("same");

  


  //pion chi2 plot

  TH1D* Chi_plot_pi = (TH1D*)Chi_plot_e->Clone("Chi_plot_pi");


  Chi_plot_pi->SetLineColor(kBlue);
  T->Draw("R.tr.chi2>>Chi_plot_pi",pi_cut,"same");


  // no hits chi2 plot

  TH1D* Chi_plot_zer = (TH1D*)Chi_plot_e->Clone("Chi_plot_zer");


  Chi_plot_zer->SetLineColor(kBlack);
  T->Draw("R.tr.chi2>>Chi_plot_zer",zer_cut,"same");


  // background chi2 plot

  TH1D* Chi_plot_bg = (TH1D*)Chi_plot_e->Clone("Chi_plot_bg");


  Chi_plot_bg->SetLineColor(kGreen);
  T->Draw("R.tr.chi2>>Chi_plot_bg",bg_cut,"same");




  //~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // From above plots on 2nd canvas create an over 'good' electron/general cut


  TCut gen_cut = Form("R.tr.n==1 && R.gold.beta>%f && R.tr.chi2<%f && R.gold.th>%f && R.gold.th <%f && R.gold.ph>%f && R.gold.ph <%f && R.gold.dp>%f && R.gold.dp<%f", beta_cut_l, chi2_cut_l, theta_cut_l, theta_cut_h, phi_cut_l, phi_cut_h, dp_cut_l, dp_cut_h);








  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~3
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~3

  // Plot track projection distribution in s0

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~3

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~3
  
  TCanvas* c3 = new TCanvas("c3","S0 variables",1200,1200);
    
  c3->Divide(2,2);
  c3->cd(1);
  // gPad->SetLogy();


  // set-up s0 y-dim for al categories of events

  Double_t s0_ymin = -0.4;
  Double_t s0_ymax = 0.4;
  
  Double_t s0_xmin = -1.3;
  Double_t s0_xmax = 1.3;


  
  TH2D *s0_plot_e = new TH2D("s0_plot_e","Electrons",200,s0_ymin,s0_ymax,200,s0_xmin,s0_xmax);
  
  s0_plot_e->GetXaxis()->SetTitle("s0 Y track projection [m]");
  s0_plot_e->GetXaxis()->CenterTitle();
  s0_plot_e->GetYaxis()->SetTitle("s0 X track projection [m]");
  s0_plot_e->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.s0.trx:R.s0.try>>s0_plot_e",e_cut,"colz");

  // draw box to represent s0 dimensions

  Double_t  s0_x_height = 1.7;
  Double_t  s0_y_width = 0.25;

  TBox* s0_box = new TBox(-s0_y_width/2,-s0_x_height/2,+s0_y_width/2,+s0_x_height/2);

  s0_box->SetFillStyle(0);
  s0_box->SetLineColor(kRed);
  // TLine* s0_line_1 = new TLine(-s0_y_width/2,-s0_x_width/2,-s0_y_width/2,+s0_x_width/2);
  // TLine* s0_line_2 = new TLine(-s0_y_width/2,+s0_x_width/2,+s0_y_width/2,+s0_x_width/2);
  // TLine* s0_line_3 = new TLine(+s0_y_width/2,-s0_x_width/2,+s0_y_width/2,+s0_x_width/2);
  // TLine* s0_line_4 = new TLine(-s0_y_width/2,-s0_x_width/2,+s0_y_width/2,-s0_x_width/2);

  s0_box->Draw("same");


  c3->Update();

  c3->cd(2);
  // gPad->SetLogy();

  TH2D *s0_plot_pi = new TH2D("s0_plot_pi","Pions",200,s0_ymin,s0_ymax,200,s0_xmin,s0_xmax);
  
  s0_plot_pi->GetXaxis()->SetTitle("s0 Y track projection [m]");
  s0_plot_pi->GetXaxis()->CenterTitle();
  s0_plot_pi->GetYaxis()->SetTitle("s0 X track projection [m]");
  s0_plot_pi->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.s0.trx:R.s0.try>>s0_plot_pi",pi_cut,"colz");

  s0_box->Draw("same");


  c3->cd(3);
  // gPad->SetLogy();

  TH2D *s0_plot_zer = new TH2D("s0_plot_zer","Zero hits",200,s0_ymin,s0_ymax,200,s0_xmin,s0_xmax);
  
  s0_plot_zer->GetXaxis()->SetTitle("s0 Y track projection [m]");
  s0_plot_zer->GetXaxis()->CenterTitle();
  s0_plot_zer->GetYaxis()->SetTitle("s0 X track projection [m]");
  s0_plot_zer->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.s0.trx:R.s0.try>>s0_plot_zer",zer_cut,"colz");

  s0_box->Draw("same");




  c3->cd(4);
  // gPad->SetLogy();

  TH2D *s0_plot_bg = new TH2D("s0_plot_bg","Bg hits",200,s0_ymin,s0_ymax,200,s0_xmin,s0_xmax);
  
  s0_plot_bg->GetXaxis()->SetTitle("s0 Y track projection [m]");
  s0_plot_bg->GetXaxis()->CenterTitle();
  s0_plot_bg->GetYaxis()->SetTitle("s0 X track projection [m]");
  s0_plot_bg->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.s0.trx:R.s0.try>>s0_plot_bg",bg_cut,"colz");

  s0_box->Draw("same");


  // gPad->SetLogy();



  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~4
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~4

  // Plot track projection distribution in s2

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~4

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~4
  
  TCanvas* c4 = new TCanvas("c4","S2 variables",1200,1200);
    
  c4->Divide(2,2);
  c4->cd(1);
  // gPad->SetLogy();


  // set-up s2 y-dim for al categories of events

  Double_t s2_ymin = -0.4;
  Double_t s2_ymax = 0.4;
  
  Double_t s2_xmin = -1.3;
  Double_t s2_xmax = 1.3;


  
  TH2D *s2_plot_e = new TH2D("s2_plot_e","Electrons",200,s2_ymin,s2_ymax,200,s2_xmin,s2_xmax);
  
  s2_plot_e->GetXaxis()->SetTitle("s2 Y track projection [m]");
  s2_plot_e->GetXaxis()->CenterTitle();
  s2_plot_e->GetYaxis()->SetTitle("s2 X track projection [m]");
  s2_plot_e->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.s2.trx:R.s2.try>>s2_plot_e",e_cut,"colz");

  // draw box to represent s2 dimensions

  Double_t  s2_x_height = 2.2353;
  Double_t  s2_y_width = 0.4318;

  TBox* s2_box = new TBox(-s2_y_width/2,-s2_x_height/2,+s2_y_width/2,+s2_x_height/2);

  s2_box->SetFillStyle(0);
  s2_box->SetLineColor(kRed);
  // TLine* s2_line_1 = new TLine(-s2_y_width/2,-s2_x_width/2,-s2_y_width/2,+s2_x_width/2);
  // TLine* s2_line_2 = new TLine(-s2_y_width/2,+s2_x_width/2,+s2_y_width/2,+s2_x_width/2);
  // TLine* s2_line_3 = new TLine(+s2_y_width/2,-s2_x_width/2,+s2_y_width/2,+s2_x_width/2);
  // TLine* s2_line_4 = new TLine(-s2_y_width/2,-s2_x_width/2,+s2_y_width/2,-s2_x_width/2);

  s2_box->Draw("same");


  c4->Update();

  c4->cd(2);
  // gPad->SetLogy();

  TH2D *s2_plot_pi = new TH2D("s2_plot_pi","Pions",200,s2_ymin,s2_ymax,200,s2_xmin,s2_xmax);
  
  s2_plot_pi->GetXaxis()->SetTitle("s2 Y track projection [m]");
  s2_plot_pi->GetXaxis()->CenterTitle();
  s2_plot_pi->GetYaxis()->SetTitle("s2 X track projection [m]");
  s2_plot_pi->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.s2.trx:R.s2.try>>s2_plot_pi",pi_cut,"colz");

  s2_box->Draw("same");


  c4->cd(3);
  // gPad->SetLogy();

  TH2D *s2_plot_zer = new TH2D("s2_plot_zer","Zero hits",200,s2_ymin,s2_ymax,200,s2_xmin,s2_xmax);
  
  s2_plot_zer->GetXaxis()->SetTitle("s2 Y track projection [m]");
  s2_plot_zer->GetXaxis()->CenterTitle();
  s2_plot_zer->GetYaxis()->SetTitle("s2 X track projection [m]");
  s2_plot_zer->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.s2.trx:R.s2.try>>s2_plot_zer",zer_cut,"colz");

  s2_box->Draw("same");




  c4->cd(4);
  // gPad->SetLogy();

  TH2D *s2_plot_bg = new TH2D("s2_plot_bg","Bg hits",200,s2_ymin,s2_ymax,200,s2_xmin,s2_xmax);
  
  s2_plot_bg->GetXaxis()->SetTitle("s2 Y track projection [m]");
  s2_plot_bg->GetXaxis()->CenterTitle();
  s2_plot_bg->GetYaxis()->SetTitle("s2 X track projection [m]");
  s2_plot_bg->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.s2.trx:R.s2.try>>s2_plot_bg",bg_cut,"colz");

  s2_box->Draw("same");


  // gPad->SetLogy();




  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~5
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~5

  // Plot track projection distribution in PS

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~5

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~5

  
  TCanvas* c5 = new TCanvas("c5","PS variables",1200,1200);
    
  c5->Divide(2,2);
  c5->cd(1);
  // gPad->SetLogy();


  // set-up ps y-dim for al categories of events

  Double_t ps_ymin = -0.4;
  Double_t ps_ymax = 0.4;
  
  Double_t ps_xmin = -1.5;
  Double_t ps_xmax = 1.5;


  
  TH2D *ps_plot_e = new TH2D("ps_plot_e","Electrons",200,ps_ymin,ps_ymax,200,ps_xmin,ps_xmax);
  
  ps_plot_e->GetXaxis()->SetTitle("ps Y track projection [m]");
  ps_plot_e->GetXaxis()->CenterTitle();
  ps_plot_e->GetYaxis()->SetTitle("ps X track projection [m]");
  ps_plot_e->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.ps.trx:R.ps.try>>ps_plot_e",e_cut,"colz");

  // draw box to represent ps dimensions

  Double_t  ps_x_height = 2.4;
  Double_t  ps_y_width = 0.7;

  TBox* ps_box = new TBox(-ps_y_width/2,-ps_x_height/2,+ps_y_width/2,+ps_x_height/2);

  ps_box->SetFillStyle(0);
  ps_box->SetLineColor(kRed);


  ps_box->Draw("same");


  c5->Update();

  c5->cd(2);
  // gPad->SetLogy();

  TH2D *ps_plot_pi = new TH2D("ps_plot_pi","Pions",200,ps_ymin,ps_ymax,200,ps_xmin,ps_xmax);
  
  ps_plot_pi->GetXaxis()->SetTitle("ps Y track projection [m]");
  ps_plot_pi->GetXaxis()->CenterTitle();
  ps_plot_pi->GetYaxis()->SetTitle("ps X track projection [m]");
  ps_plot_pi->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.ps.trx:R.ps.try>>ps_plot_pi",pi_cut,"colz");

  ps_box->Draw("same");


  c5->cd(3);
  // gPad->SetLogy();

  TH2D *ps_plot_zer = new TH2D("ps_plot_zer","Zero hits",200,ps_ymin,ps_ymax,200,ps_xmin,ps_xmax);
  
  ps_plot_zer->GetXaxis()->SetTitle("ps Y track projection [m]");
  ps_plot_zer->GetXaxis()->CenterTitle();
  ps_plot_zer->GetYaxis()->SetTitle("ps X track projection [m]");
  ps_plot_zer->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.ps.trx:R.ps.try>>ps_plot_zer",zer_cut,"colz");

  ps_box->Draw("same");




  c5->cd(4);
  // gPad->SetLogy();

  TH2D *ps_plot_bg = new TH2D("ps_plot_bg","Bg hits",200,ps_ymin,ps_ymax,200,ps_xmin,ps_xmax);
  
  ps_plot_bg->GetXaxis()->SetTitle("ps Y track projection [m]");
  ps_plot_bg->GetXaxis()->CenterTitle();
  ps_plot_bg->GetYaxis()->SetTitle("ps X track projection [m]");
  ps_plot_bg->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.ps.trx:R.ps.try>>ps_plot_bg",bg_cut,"colz");

  ps_box->Draw("same");

  

  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~6
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~6

  // Plot track projection distribution in SH

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~6

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~6

  
  TCanvas* c6 = new TCanvas("c6","SH variables",1200,1200);
    
  c6->Divide(2,2);
  c6->cd(1);
  // gPad->SetLogy();


  // set-up sh y-dim for al categories of events

  Double_t sh_ymin = -0.4;
  Double_t sh_ymax = 0.4;
  
  Double_t sh_xmin = -1.5;
  Double_t sh_xmax = 1.5;


  
  TH2D *sh_plot_e = new TH2D("sh_plot_e","Electrons",200,sh_ymin,sh_ymax,200,sh_xmin,sh_xmax);
  
  sh_plot_e->GetXaxis()->SetTitle("sh Y track projection [m]");
  sh_plot_e->GetXaxis()->CenterTitle();
  sh_plot_e->GetYaxis()->SetTitle("sh X track projection [m]");
  sh_plot_e->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.sh.trx:R.sh.try>>sh_plot_e",e_cut,"colz");

  // draw box to represent sh dimensions

  Double_t  sh_x_height = 2.25;
  Double_t  sh_y_width = 0.75;

  TBox* sh_box = new TBox(-sh_y_width/2,-sh_x_height/2,+sh_y_width/2,+sh_x_height/2);

  sh_box->SetFillStyle(0);
  sh_box->SetLineColor(kRed);


  sh_box->Draw("same");


  c6->Update();

  c6->cd(2);
  // gPad->SetLogy();

  TH2D *sh_plot_pi = new TH2D("sh_plot_pi","Pions",200,sh_ymin,sh_ymax,200,sh_xmin,sh_xmax);
  
  sh_plot_pi->GetXaxis()->SetTitle("sh Y track projection [m]");
  sh_plot_pi->GetXaxis()->CenterTitle();
  sh_plot_pi->GetYaxis()->SetTitle("sh X track projection [m]");
  sh_plot_pi->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.sh.trx:R.sh.try>>sh_plot_pi",pi_cut,"colz");

  sh_box->Draw("same");


  c6->cd(3);
  // gPad->SetLogy();

  TH2D *sh_plot_zer = new TH2D("sh_plot_zer","Zero hits",200,sh_ymin,sh_ymax,200,sh_xmin,sh_xmax);
  
  sh_plot_zer->GetXaxis()->SetTitle("sh Y track projection [m]");
  sh_plot_zer->GetXaxis()->CenterTitle();
  sh_plot_zer->GetYaxis()->SetTitle("sh X track projection [m]");
  sh_plot_zer->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.sh.trx:R.sh.try>>sh_plot_zer",zer_cut,"colz");

  sh_box->Draw("same");




  c6->cd(4);
  // gPad->SetLogy();

  TH2D *sh_plot_bg = new TH2D("sh_plot_bg","Bg hits",200,sh_ymin,sh_ymax,200,sh_xmin,sh_xmax);
  
  sh_plot_bg->GetXaxis()->SetTitle("sh Y track projection [m]");
  sh_plot_bg->GetXaxis()->CenterTitle();
  sh_plot_bg->GetYaxis()->SetTitle("sh X track projection [m]");
  sh_plot_bg->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.sh.trx:R.sh.try>>sh_plot_bg",bg_cut,"colz");

  sh_box->Draw("same");




  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~7
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~7

  // Plot track projection vs Cal1 (PS) Energy x distribution 

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~7

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~7
  
  TCanvas* c7 = new TCanvas("c7","PS E_x vs tr_x variables",1200,1200);
    
  c7->Divide(2,2);
  c7->cd(1);
  // gPad->SetLogy();


  
  TH2D *ps_EvTx_plot_e = new TH2D("ps_EvTx_plot_e","Electrons",200,ps_xmin,ps_xmax,200,ps_xmin,ps_xmax);
  
  ps_EvTx_plot_e->GetXaxis()->SetTitle("VDC X projection [m]");
  ps_EvTx_plot_e->GetXaxis()->CenterTitle();
  ps_EvTx_plot_e->GetYaxis()->SetTitle("PS Energy-weighted x [m]");
  ps_EvTx_plot_e->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.ps.x:R.ps.trx>>ps_EvTx_plot_e",e_cut,"colz");



  c7->Update();

  c7->cd(2);
  
  TH2D *ps_EvTx_plot_pi = new TH2D("ps_EvTx_plot_pi","Pions",200,ps_xmin,ps_xmax,200,ps_xmin,ps_xmax);
  
  ps_EvTx_plot_pi->GetXaxis()->SetTitle("VDC X projection [m]");
  ps_EvTx_plot_pi->GetXaxis()->CenterTitle();
  ps_EvTx_plot_pi->GetYaxis()->SetTitle("PS Energy-weighted x [m]");
  ps_EvTx_plot_pi->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.ps.x:R.ps.trx>>ps_EvTx_plot_pi",pi_cut,"colz");

  

  c7->cd(3);


  TH2D *ps_EvTx_plot_zer = new TH2D("ps_EvTx_plot_zer","Zero hits",200,ps_xmin,ps_xmax,200,ps_xmin,ps_xmax);
  
  ps_EvTx_plot_zer->GetXaxis()->SetTitle("VDC X projection [m]");
  ps_EvTx_plot_zer->GetXaxis()->CenterTitle();
  ps_EvTx_plot_zer->GetYaxis()->SetTitle("PS Energy-weighted x [m]");
  ps_EvTx_plot_zer->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.ps.x:R.ps.trx>>ps_EvTx_plot_zer",zer_cut,"colz");


  

  c7->cd(4);


  TH2D *ps_EvTx_plot_bg = new TH2D("ps_EvTx_plot_bg","Bg hits",200,ps_xmin,ps_xmax,200,ps_xmin,ps_xmax);
  
  ps_EvTx_plot_bg->GetXaxis()->SetTitle("VDC X projection [m]");
  ps_EvTx_plot_bg->GetXaxis()->CenterTitle();
  ps_EvTx_plot_bg->GetYaxis()->SetTitle("PS Energy-weighted x [m]");
  ps_EvTx_plot_bg->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.ps.x:R.ps.trx>>ps_EvTx_plot_bg",bg_cut,"colz");





  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~8
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~8

  // Plot track projection vs Cal1 (PS) Energy y distribution 

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~8

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~8
  
  TCanvas* c8 = new TCanvas("c8","PS E_y vs tr_y variables",1200,1200);
  
  c8->Divide(2,2);
  c8->cd(1);
  // gPad->SetLogy();


  
  TH2D *ps_EvTy_plot_e = new TH2D("ps_EvTy_plot_e","Electrons",200,ps_ymin,ps_ymax,200,ps_ymin,ps_ymax);
  
  ps_EvTy_plot_e->GetXaxis()->SetTitle("VDC Y projection [m]");
  ps_EvTy_plot_e->GetXaxis()->CenterTitle();
  ps_EvTy_plot_e->GetYaxis()->SetTitle("PS Energy-weighted y [m]");
  ps_EvTy_plot_e->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.ps.y:R.ps.try>>ps_EvTy_plot_e",e_cut,"colz");



  c8->Update();

  c8->cd(2);
  
  TH2D *ps_EvTy_plot_pi = new TH2D("ps_EvTy_plot_pi","Pions",200,ps_ymin,ps_ymax,200,ps_ymin,ps_ymax);
  
  ps_EvTy_plot_pi->GetXaxis()->SetTitle("VDC Y projection [m]");
  ps_EvTy_plot_pi->GetXaxis()->CenterTitle();
  ps_EvTy_plot_pi->GetYaxis()->SetTitle("PS Energy-weighted y [m]");
  ps_EvTy_plot_pi->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.ps.y:R.ps.try>>ps_EvTy_plot_pi",pi_cut,"colz");

  

  c8->cd(3);


  TH2D *ps_EvTy_plot_zer = new TH2D("ps_EvTy_plot_zer","Zero hits",200,ps_ymin,ps_ymax,200,ps_ymin,ps_ymax);
  
  ps_EvTy_plot_zer->GetXaxis()->SetTitle("VDC Y projection [m]");
  ps_EvTy_plot_zer->GetXaxis()->CenterTitle();
  ps_EvTy_plot_zer->GetYaxis()->SetTitle("PS Energy-weighted y [m]");
  ps_EvTy_plot_zer->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.ps.y:R.ps.try>>ps_EvTy_plot_zer",zer_cut,"colz");


  

  c8->cd(4);


  TH2D *ps_EvTy_plot_bg = new TH2D("ps_EvTy_plot_bg","Bg hits",200,ps_ymin,ps_ymax,200,ps_ymin,ps_ymax);
  
  ps_EvTy_plot_bg->GetXaxis()->SetTitle("VDC Y projection [m]");
  ps_EvTy_plot_bg->GetXaxis()->CenterTitle();
  ps_EvTy_plot_bg->GetYaxis()->SetTitle("PS Energy-weighted y [m]");
  ps_EvTy_plot_bg->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.ps.y:R.ps.try>>ps_EvTy_plot_bg",bg_cut,"colz");





  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~9
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~9

  // Plot track projection vs Cal1 (SH) Energy x distribution 

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~9

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~9
  
  TCanvas* c9 = new TCanvas("c9","SH E_x vs tr_x variables",1200,1200);
    
  c9->Divide(2,2);
  c9->cd(1);
  // gPad->SetLogy();


  
  TH2D *sh_EvTx_plot_e = new TH2D("sh_EvTx_plot_e","Electrons",200,sh_xmin,sh_xmax,200,sh_xmin,sh_xmax);
  
  sh_EvTx_plot_e->GetXaxis()->SetTitle("VDC X projection [m]");
  sh_EvTx_plot_e->GetXaxis()->CenterTitle();
  sh_EvTx_plot_e->GetYaxis()->SetTitle("SH Energy-weighted x [m]");
  sh_EvTx_plot_e->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.sh.x:R.sh.trx>>sh_EvTx_plot_e",e_cut,"colz");



  c9->Update();

  c9->cd(2);
  
  TH2D *sh_EvTx_plot_pi = new TH2D("sh_EvTx_plot_pi","Pions",200,sh_xmin,sh_xmax,200,sh_xmin,sh_xmax);
  
  sh_EvTx_plot_pi->GetXaxis()->SetTitle("VDC X projection [m]");
  sh_EvTx_plot_pi->GetXaxis()->CenterTitle();
  sh_EvTx_plot_pi->GetYaxis()->SetTitle("SH Energy-weighted x [m]");
  sh_EvTx_plot_pi->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.sh.x:R.sh.trx>>sh_EvTx_plot_pi",pi_cut,"colz");

  

  c9->cd(3);


  TH2D *sh_EvTx_plot_zer = new TH2D("sh_EvTx_plot_zer","Zero hits",200,sh_xmin,sh_xmax,200,sh_xmin,sh_xmax);
  
  sh_EvTx_plot_zer->GetXaxis()->SetTitle("VDC X projection [m]");
  sh_EvTx_plot_zer->GetXaxis()->CenterTitle();
  sh_EvTx_plot_zer->GetYaxis()->SetTitle("SH Energy-weighted x [m]");
  sh_EvTx_plot_zer->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.sh.x:R.sh.trx>>sh_EvTx_plot_zer",zer_cut,"colz");


  

  c9->cd(4);


  TH2D *sh_EvTx_plot_bg = new TH2D("sh_EvTx_plot_bg","Bg hits",200,sh_xmin,sh_xmax,200,sh_xmin,sh_xmax);
  
  sh_EvTx_plot_bg->GetXaxis()->SetTitle("VDC X projection [m]");
  sh_EvTx_plot_bg->GetXaxis()->CenterTitle();
  sh_EvTx_plot_bg->GetYaxis()->SetTitle("SH Energy-weighted x [m]");
  sh_EvTx_plot_bg->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.sh.x:R.sh.trx>>sh_EvTx_plot_bg",bg_cut,"colz");





  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~10
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~10

  // Plot track projection vs Cal1 (SH) Energy y distribution 

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~10

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~10
  
  TCanvas* c10 = new TCanvas("c10","SH E_y vs tr_y variables",1200,1200);
    
  c10->Divide(2,2);
  c10->cd(1);
  // gPad->SetLogy();


  
  TH2D *sh_EvTy_plot_e = new TH2D("sh_EvTy_plot_e","Electrons",200,sh_ymin,sh_ymax,200,sh_ymin,sh_ymax);
  
  sh_EvTy_plot_e->GetXaxis()->SetTitle("VDC Y projection [m]");
  sh_EvTy_plot_e->GetXaxis()->CenterTitle();
  sh_EvTy_plot_e->GetYaxis()->SetTitle("SH Energy-weighted y [m]");
  sh_EvTy_plot_e->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.sh.y:R.sh.try>>sh_EvTy_plot_e",e_cut,"colz");



  c10->Update();

  c10->cd(2);
  
  TH2D *sh_EvTy_plot_pi = new TH2D("sh_EvTy_plot_pi","Pions",200,sh_ymin,sh_ymax,200,sh_ymin,sh_ymax);
  
  sh_EvTy_plot_pi->GetXaxis()->SetTitle("VDC Y projection [m]");
  sh_EvTy_plot_pi->GetXaxis()->CenterTitle();
  sh_EvTy_plot_pi->GetYaxis()->SetTitle("SH Energy-weighted y [m]");
  sh_EvTy_plot_pi->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.sh.y:R.sh.try>>sh_EvTy_plot_pi",pi_cut,"colz");

  

  c10->cd(3);


  TH2D *sh_EvTy_plot_zer = new TH2D("sh_EvTy_plot_zer","Zero hits",200,sh_ymin,sh_ymax,200,sh_ymin,sh_ymax);
  
  sh_EvTy_plot_zer->GetXaxis()->SetTitle("VDC Y projection [m]");
  sh_EvTy_plot_zer->GetXaxis()->CenterTitle();
  sh_EvTy_plot_zer->GetYaxis()->SetTitle("SH Energy-weighted y [m]");
  sh_EvTy_plot_zer->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.sh.y:R.sh.try>>sh_EvTy_plot_zer",zer_cut,"colz");


  

  c10->cd(4);


  TH2D *sh_EvTy_plot_bg = new TH2D("sh_EvTy_plot_bg","Bg hits",200,sh_ymin,sh_ymax,200,sh_ymin,sh_ymax);
  
  sh_EvTy_plot_bg->GetXaxis()->SetTitle("VDC Y projection [m]");
  sh_EvTy_plot_bg->GetXaxis()->CenterTitle();
  sh_EvTy_plot_bg->GetYaxis()->SetTitle("SH Energy-weighted y [m]");
  sh_EvTy_plot_bg->GetYaxis()->CenterTitle(); 

    
  T->Draw("R.sh.y:R.sh.try>>sh_EvTy_plot_bg",bg_cut,"colz");




  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~11
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~11

  // Plot 1D difference between track projection and Cal x/y

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~11

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~11


  
  
  TCanvas* c11 = new TCanvas("c11","Track - Calorimeter x/y",1200,1200);

  c11->Divide(2);

  c11->cd(1);

  TH1D *ps_E_min_Tx_plot_e = new TH1D("ps_E_min_Tx_plot_e","PS X",200,-0.1,0.1);

  ps_E_min_Tx_plot_e->GetXaxis()->SetTitle("PS - VDC X [m]");
  ps_E_min_Tx_plot_e->GetXaxis()->CenterTitle();


  ps_E_min_Tx_plot_e->SetLineColor(kRed);
  ps_E_min_Tx_plot_e->SetLineWidth(2);
  
  T->Draw("R.ps.x-R.ps.trx>>ps_E_min_Tx_plot_e",e_cut);



  // draw cuts for ps_E_x - ps_track_x


  Double_t ps_E_min_Tx_l = -0.055;
  Double_t ps_E_min_Tx_h = 0.050;
  
  TLine* ps_E_min_Tx_cut_line_l = new TLine(ps_E_min_Tx_l,0,ps_E_min_Tx_l,ps_E_min_Tx_plot_e->GetMaximum());


  ps_E_min_Tx_cut_line_l->SetLineColor(kMagenta);  
  ps_E_min_Tx_cut_line_l->SetLineStyle(9);
  ps_E_min_Tx_cut_line_l->SetLineWidth(2);
  ps_E_min_Tx_cut_line_l->Draw("same");


  TLine* ps_E_min_Tx_cut_line_h = new TLine(ps_E_min_Tx_h,0,ps_E_min_Tx_h,ps_E_min_Tx_plot_e->GetMaximum());


  ps_E_min_Tx_cut_line_h->SetLineColor(kMagenta);  
  ps_E_min_Tx_cut_line_h->SetLineStyle(9);
  ps_E_min_Tx_cut_line_h->SetLineWidth(2);
  ps_E_min_Tx_cut_line_h->Draw("same");

  
  
  TH1D* ps_E_min_Tx_plot_pi = (TH1D*)ps_E_min_Tx_plot_e->Clone("ps_E_min_Tx_plot_pi");
  ps_E_min_Tx_plot_pi->SetLineColor(kBlue);

  T->Draw("R.ps.x-R.ps.trx>>ps_E_min_Tx_plot_pi",pi_cut,"same");

  
  TH1D* ps_E_min_Tx_plot_zer = (TH1D*)ps_E_min_Tx_plot_e->Clone("ps_E_min_Tx_plot_zer");
  ps_E_min_Tx_plot_zer->SetLineColor(kBlack);

  T->Draw("R.ps.x-R.ps.trx>>ps_E_min_Tx_plot_zer",zer_cut,"same");

  
  TH1D* ps_E_min_Tx_plot_bg = (TH1D*)ps_E_min_Tx_plot_e->Clone("ps_E_min_Tx_plot_bg");
  ps_E_min_Tx_plot_bg->SetLineColor(kGreen);


  T->Draw("R.ps.x-R.ps.trx>>ps_E_min_Tx_plot_bg",bg_cut,"same");

  
  
  c11->cd(2);


  TH1D *sh_E_min_Tx_plot_e = new TH1D("sh_E_min_Tx_plot_e","SH X",200,-0.1,0.1);

  sh_E_min_Tx_plot_e->GetXaxis()->SetTitle("SH - VDC X [m]");
  sh_E_min_Tx_plot_e->GetXaxis()->CenterTitle();


  sh_E_min_Tx_plot_e->SetLineColor(kRed);
  sh_E_min_Tx_plot_e->SetLineWidth(2);
  
  T->Draw("R.sh.x-R.sh.trx>>sh_E_min_Tx_plot_e",e_cut);


  // draw cuts for sh_E_x - sh_track_x


  Double_t sh_E_min_Tx_l = -0.065;
  Double_t sh_E_min_Tx_h = 0.055;
  
  TLine* sh_E_min_Tx_cut_line_l = new TLine(sh_E_min_Tx_l,0,sh_E_min_Tx_l,sh_E_min_Tx_plot_e->GetMaximum());


  sh_E_min_Tx_cut_line_l->SetLineColor(kMagenta);  
  sh_E_min_Tx_cut_line_l->SetLineStyle(9);
  sh_E_min_Tx_cut_line_l->SetLineWidth(2);
  sh_E_min_Tx_cut_line_l->Draw("same");


  TLine* sh_E_min_Tx_cut_line_h = new TLine(sh_E_min_Tx_h,0,sh_E_min_Tx_h,sh_E_min_Tx_plot_e->GetMaximum());


  sh_E_min_Tx_cut_line_h->SetLineColor(kMagenta);  
  sh_E_min_Tx_cut_line_h->SetLineStyle(9);
  sh_E_min_Tx_cut_line_h->SetLineWidth(2);
  sh_E_min_Tx_cut_line_h->Draw("same");




  
  TH1D* sh_E_min_Tx_plot_pi = (TH1D*)sh_E_min_Tx_plot_e->Clone("sh_E_min_Tx_plot_pi");
  sh_E_min_Tx_plot_pi->SetLineColor(kBlue);

  T->Draw("R.sh.x-R.sh.trx>>sh_E_min_Tx_plot_pi",pi_cut,"same");

  
  TH1D* sh_E_min_Tx_plot_zer = (TH1D*)sh_E_min_Tx_plot_e->Clone("sh_E_min_Tx_plot_zer");
  sh_E_min_Tx_plot_zer->SetLineColor(kBlack);

  T->Draw("R.sh.x-R.sh.trx>>sh_E_min_Tx_plot_zer",zer_cut,"same");

  
  TH1D* sh_E_min_Tx_plot_bg = (TH1D*)sh_E_min_Tx_plot_e->Clone("sh_E_min_Tx_plot_bg");
  sh_E_min_Tx_plot_bg->SetLineColor(kGreen);


  T->Draw("R.sh.x-R.sh.trx>>sh_E_min_Tx_plot_bg",bg_cut,"same");

 

  // currently commented out y diff between prl and track
  // only 2 columns so this give little info
  
  // c11->cd(3);

  // TH1D *ps_E_min_Ty_plot_e = new TH1D("ps_E_min_Ty_plot_e","PS Y",200,-0.1,0.1);

  // ps_E_min_Ty_plot_e->GetXaxis()->SetTitle("PS - VDC Y [m]");
  // ps_E_min_Ty_plot_e->GetXaxis()->CenterTitle();


  // ps_E_min_Ty_plot_e->SetLineColor(kRed);
  // ps_E_min_Ty_plot_e->SetLineWidth(2);
  
  // T->Draw("R.ps.y-R.ps.try>>ps_E_min_Ty_plot_e",e_cut);

  
  // TH1D* ps_E_min_Ty_plot_pi = (TH1D*)ps_E_min_Ty_plot_e->Clone("ps_E_min_Ty_plot_pi");
  // ps_E_min_Ty_plot_pi->SetLineColor(kBlue);

  // T->Draw("R.ps.y-R.ps.try>>ps_E_min_Ty_plot_pi",pi_cut,"same");

  
  // TH1D* ps_E_min_Ty_plot_zer = (TH1D*)ps_E_min_Ty_plot_e->Clone("ps_E_min_Ty_plot_zer");
  // ps_E_min_Ty_plot_zer->SetLineColor(kBlack);

  // T->Draw("R.ps.y-R.ps.try>>ps_E_min_Ty_plot_zer",zer_cut,"same");

  
  // TH1D* ps_E_min_Ty_plot_bg = (TH1D*)ps_E_min_Ty_plot_e->Clone("ps_E_min_Ty_plot_bg");
  // ps_E_min_Ty_plot_bg->SetLineColor(kGreen);


  // T->Draw("R.ps.y-R.ps.try>>ps_E_min_Ty_plot_bg",bg_cut,"same");



  // c11->cd(4);

  // TH1D *sh_E_min_Ty_plot_e = new TH1D("sh_E_min_Ty_plot_e","SH Y",200,-0.1,0.1);

  // sh_E_min_Ty_plot_e->GetXaxis()->SetTitle("SH - VDC Y [m]");
  // sh_E_min_Ty_plot_e->GetXaxis()->CenterTitle();


  // sh_E_min_Ty_plot_e->SetLineColor(kRed);
  // sh_E_min_Ty_plot_e->SetLineWidth(2);
  
  // T->Draw("R.sh.y-R.sh.try>>sh_E_min_Ty_plot_e",e_cut);

  
  // TH1D* sh_E_min_Ty_plot_pi = (TH1D*)sh_E_min_Ty_plot_e->Clone("sh_E_min_Ty_plot_pi");
  // sh_E_min_Ty_plot_pi->SetLineColor(kBlue);

  // T->Draw("R.sh.y-R.sh.try>>sh_E_min_Ty_plot_pi",pi_cut,"same");

  
  // TH1D* sh_E_min_Ty_plot_zer = (TH1D*)sh_E_min_Ty_plot_e->Clone("sh_E_min_Ty_plot_zer");
  // sh_E_min_Ty_plot_zer->SetLineColor(kBlack);

  // T->Draw("R.sh.y-R.sh.try>>sh_E_min_Ty_plot_zer",zer_cut,"same");

  
  // TH1D* sh_E_min_Ty_plot_bg = (TH1D*)sh_E_min_Ty_plot_e->Clone("sh_E_min_Ty_plot_bg");
  // sh_E_min_Ty_plot_bg->SetLineColor(kGreen);


  // T->Draw("R.sh.y-R.sh.try>>sh_E_min_Ty_plot_bg",bg_cut,"same");



 

  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~50
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // Cherenkov Sum Scan
  // aim to find level of Cherenkov sum cut that optimises electron to pion ratio

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~50

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~50



  TCanvas* c50 = new TCanvas("c50","Cherenkov Cut Scan ",1200,1200);

  c50->Divide(3,2);
  


  // firstly display cuts used to select clean sample of 'good' electrons from Calorimeters for determining cherenkov PID effenciecy


  // First plot sum of Calorimeter energies (displaying cuts for electron + pion) (energy normalised to track momentum)

  c50->cd(1);

  gPad->SetLogy();
  
  TH1D* Cal_EP_plot =  new TH1D("Cal_EP_plot","PS + SH Energy",100,-.05,1.5);


  Cal_EP_plot->SetLineColor(kBlack);
  

  Cal_EP_plot->GetXaxis()->SetTitle("(PS.e + SH.e)/R.gold.p");
  Cal_EP_plot->GetXaxis()->CenterTitle();


  T->Draw("(R.sh.e + R.ps.e)/(R.gold.p*1000)>>Cal_EP_plot",gen_cut,"");



  // electron cut on PRL sum
    

  Double_t e_ps_2_min = 0.85;
  Double_t e_ps_2_max = 1.15;


  TLine* e_prl_com_l = new TLine(e_ps_2_min,0,e_ps_2_min,Cal_EP_plot->GetMaximum());

  e_prl_com_l->SetLineWidth(2);
  e_prl_com_l->SetLineColor(kRed);
  
  e_prl_com_l->Draw("same");


  TLine* e_prl_com_h = new TLine(e_ps_2_max,0,e_ps_2_max,Cal_EP_plot->GetMaximum());

  e_prl_com_h->SetLineWidth(2);
  e_prl_com_h->SetLineColor(kRed);
  
  e_prl_com_h->Draw("same");

  
  
  // add legend for Calorimeter event selection

  TLegend* leg_Cher_scan = new TLegend(.1,.7,.3,.9,"Key");
  leg_Cher_scan->SetFillColor(0);
  leg_Cher_scan->AddEntry(e_cut_box,"Electron cut","l");
  leg_Cher_scan->AddEntry(pi_cut_box,"Pion cut","l");
  leg_Cher_scan->Draw("same");






  // plot second calorimeter energy (energy normalised to track momentum) with cuts for electrons and pions
  
  
  c50->cd(2);

  
  gPad->SetLogy();

  TH1D* Cal2_EP_plot =  new TH1D("Cal2_EP_plot","SH Energy",100,-0.05,1.2);
  
  Cal2_EP_plot->GetXaxis()->SetTitle("(SH.e)/R.gold.p");
  Cal2_EP_plot->GetXaxis()->CenterTitle();
  
  Cal2_EP_plot->SetLineColor(kBlack);


  T->Draw("(R.sh.e)/(R.gold.p*1000)>>Cal2_EP_plot",gen_cut,"");
  


  // electron cuts based on SH


 Double_t e_sh_min = 0.15;
 Double_t e_sh_max = 0.80;

 
 TLine* e_sh_s_l= new TLine(e_sh_min,0,e_sh_min,Cal2_EP_plot->GetMaximum());
 TLine* e_sh_s_h= new TLine(e_sh_max,0,e_sh_max,Cal2_EP_plot->GetMaximum());
  
 
 e_sh_s_l->SetLineWidth(2);
 e_sh_s_l->SetLineColor(kRed);
 
 e_sh_s_l->Draw("same");
 
 e_sh_s_h->SetLineWidth(2);
 e_sh_s_h->SetLineColor(kRed);
 
 e_sh_s_h->Draw("same");
  


 // pion cuts based on SH

 Double_t pi_sh_min = 0.09;
 Double_t pi_sh_max = 0.15;

 TLine* pi_sh_s_l= new TLine(pi_sh_min,0,pi_sh_min,Cal2_EP_plot->GetMaximum());
 TLine* pi_sh_s_h= new TLine(pi_sh_max,0,pi_sh_max,Cal2_EP_plot->GetMaximum());


 pi_sh_s_l->SetLineWidth(2);
 pi_sh_s_l->SetLineColor(kBlue);
 
 pi_sh_s_l->Draw("same");
 
 pi_sh_s_h->SetLineWidth(2);
 pi_sh_s_h->SetLineColor(kBlue);
 
 pi_sh_s_h->Draw("same");




  // plot first calorimeter energy (energy normalised to track momentum) with cuts for electrons and pions

 
 c50->cd(3);

 
 gPad->SetLogy();
 
 TH1D* Cal1_EP_plot =  new TH1D("Cal1_EP_plot","PS Energy",100,-0.05,1.2);
 
 
 Cal1_EP_plot->GetXaxis()->SetTitle("(PS.e)/R.gold.p");
 Cal1_EP_plot->GetXaxis()->CenterTitle();
  
 Cal1_EP_plot->SetLineColor(kBlack);
  
 T->Draw("(R.ps.e)/(R.gold.p*1000)>>Cal1_EP_plot",gen_cut,"");
    
 


  
 // electron cuts based on PS
  

  Double_t e_ps_min_1 = 0.25;
  Double_t e_ps_max_1 = 0.8;
  
 
  TLine* e_ps_s_l= new TLine(e_ps_min_1,0,e_ps_min_1,Cal2_EP_plot->GetMaximum());
  TLine* e_ps_s_h= new TLine(e_ps_max_1,0,e_ps_max_1,Cal2_EP_plot->GetMaximum());
  
  
  e_ps_s_l->SetLineWidth(2);
  e_ps_s_l->SetLineColor(kRed);
  
  e_ps_s_l->Draw("same");
  
  e_ps_s_h->SetLineWidth(2);
  e_ps_s_h->SetLineColor(kRed);
  
  e_ps_s_h->Draw("same");
 

  
  // pion cuts for PS


  Double_t pi_ps_min = 0.075;
  Double_t pi_ps_max = 0.16;
  
  TLine* pi_ps_s_l= new TLine(pi_ps_min,0,pi_ps_min,Cal2_EP_plot->GetMaximum());
  TLine* pi_ps_s_h= new TLine(pi_ps_max,0,pi_ps_max,Cal2_EP_plot->GetMaximum());
  
  
  pi_ps_s_l->SetLineWidth(2);
  pi_ps_s_l->SetLineColor(kBlue);
  
  pi_ps_s_l->Draw("same");
 
  pi_ps_s_h->SetLineWidth(2);
  pi_ps_s_h->SetLineColor(kBlue);
  
  pi_ps_s_h->Draw("same");





  // plot second vs first calorimeter energy (normalised to track momentum) with cuts for electrons and pions



  c50->cd(4);
  

  TH2D *Cal_2Vs1 = new TH2D("Cal_2Vs1","SH vs PS Energy",1000,0.0,1.2,1000,0.0,1.2);

  Cal_2Vs1->GetXaxis()->SetTitle("PS.e/R.gold.p");
  Cal_2Vs1->GetXaxis()->CenterTitle();
  Cal_2Vs1->GetYaxis()->SetTitle("SH.e/R.gold.p");
  Cal_2Vs1->GetYaxis()->CenterTitle();
  
  
  T->Draw("(R.sh.e/(R.gold.p*1000)):(R.ps.e/(R.gold.p*1000))>>Cal_2Vs1",gen_cut,"col");
  

  // Draw horizontal line for electron based on sh energy depostition


  TLine* e_sh_l= new TLine(0,e_sh_min,1.2,e_sh_min);
  TLine* e_sh_h= new TLine(0,e_sh_max,1.2,e_sh_max);
  
  e_sh_l->SetLineWidth(2);
  e_sh_l->SetLineColor(kRed);
  
  e_sh_l->Draw("same");

  e_sh_h->SetLineWidth(2);
  e_sh_h->SetLineColor(kRed);

  e_sh_h->Draw("same");
  
  

  // draw vertical cuts

  TLine* e_ps_l= new TLine(e_ps_min_1,0,e_ps_min_1,1.2);
  TLine* e_ps_h= new TLine(e_ps_max_1,0,e_ps_max_1,1.2);


  e_ps_l->SetLineWidth(2);
  e_ps_l->SetLineColor(kRed);
  
  e_ps_l->Draw("same");

  e_ps_h->SetLineWidth(2);
  e_ps_h->SetLineColor(kRed);

  e_ps_h->Draw("same");



  // draw slopes in Calorimeter 2 vs Calorimeter 1 (based on cuts on Cal2_energy + Cal1_energy >(<) x

  Double_t e_sh_1_min   =  e_ps_2_min;
  Double_t e_sh_1_max   =  e_ps_2_max;

  Double_t e_ps_min = e_ps_2_min;
  Double_t e_ps_max = e_ps_2_max;


  
  TLine* e_sh_1_l = new TLine(0,e_sh_1_min,e_ps_min,0);
  TLine* e_sh_1_h = new TLine(0,e_sh_1_max,e_ps_max,0);
  
  
  e_sh_1_l->SetLineWidth(2);
  e_sh_1_l->SetLineColor(kRed);
  
  e_sh_1_l->Draw("same");

  e_sh_1_h->SetLineWidth(2);
  e_sh_1_h->SetLineColor(kRed);

  e_sh_1_h->Draw("same");



  // draw pion cuts

  
  TLine* pi_sh_h = new TLine(0, pi_sh_max,1.2,pi_sh_max);

  pi_sh_h->SetLineWidth(2);
  pi_sh_h->SetLineColor(kBlue);
  pi_sh_h->Draw("same");


  TLine* pi_sh_l = new TLine(0, pi_sh_min,1.2,pi_sh_min);

  pi_sh_l->SetLineWidth(2);
  pi_sh_l->SetLineColor(kBlue);
  pi_sh_l->Draw("same");


  
  TLine* pi_ps_h = new TLine( pi_ps_max,0,pi_ps_max,1.2);

  pi_ps_h->SetLineWidth(2);
  pi_ps_h->SetLineColor(kBlue);
  pi_ps_h->Draw("same");
  
  
  TLine* pi_ps_l = new TLine(pi_ps_min,0,pi_ps_min,1.2);

  pi_ps_l->SetLineWidth(2);
  pi_ps_l->SetLineColor(kBlue);
  pi_ps_l->Draw("same");




  

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Plot Cherenkov Sum spectrum for chosen electrons
  // and pions (and total)

  c50->cd(5);
  gPad->SetLogy();

  TH1D *cer_scan = new TH1D("cer_scan","Cherenkov Sum",100,0,10000);
  


  // plot all events (that pass general cut) for sanity

  cer_scan->SetLineColor(kBlack);
  

  cer_scan->GetXaxis()->SetTitle("Cherenkov Sum");
  cer_scan->GetXaxis()->CenterTitle();
  cer_scan->GetYaxis()->SetTitle("Entries");
  cer_scan->GetYaxis()->CenterTitle();

  T->Draw("R.cer.asum_c>>cer_scan",gen_cut,"");



  // plot electrons general cut as depicted in the 2nd canvas and Calorimeter based cuts in this canvas


  TCut e_c_scan_cut = gen_cut + Form("(R.ps.e/(R.gold.p*1000)) > %f && (R.ps.e/(R.gold.p*1000)) < %f && (R.sh.e/(R.gold.p*1000)) > %f && (R.sh.e/(R.gold.p*1000)) < %f && ((R.sh.e + R.ps.e)/(R.gold.p*1000)) > %f && ((R.sh.e + R.ps.e)/(R.gold.p*1000)) < %f",e_ps_min_1,e_ps_max_1,e_sh_min,e_sh_max, e_ps_2_min, e_ps_2_max);

  cout << "e_c_scan_cut = " << endl;
  cout << e_c_scan_cut << endl << endl;
  
 

   TH1D* cer_scan_e = (TH1D*)cer_scan->Clone("cer_scan_e");


   cer_scan_e->Reset();
   

   cer_scan_e->SetLineColor(kRed);
   
   T->Draw("R.cer.asum_c>>cer_scan_e",e_c_scan_cut,"same");


   

  // plot pions with general cut as displayed in canvas 2 and pions cuts from the calorimeter as shown in this canvas
  

  TCut pi_c_scan_cut = gen_cut + Form("(R.ps.e/(R.gold.p*1000)) > %f && (R.ps.e/(R.gold.p*1000)) < %f && (R.sh.e/(R.gold.p*1000)) > %f && (R.sh.e/(R.gold.p*1000)) < %f",pi_ps_min,pi_ps_max,pi_sh_min,pi_sh_max);


  cout << "pi_c_scan_cut = " << endl;
  cout << pi_c_scan_cut << endl << endl;
  

  TH1D* cer_scan_pi = (TH1D*)cer_scan->Clone("cer_scan_pi");

  cer_scan_pi->Reset();

  cer_scan_pi->SetLineColor(kBlue);
  
  T->Draw("R.cer.asum_c>>cer_scan_pi",pi_c_scan_cut,"same");



   
  TLegend* leg_Cher_sum = new TLegend(.1,.7,.3,.9,"Key");
  leg_Cher_sum->SetFillColor(0);
  leg_Cher_sum->AddEntry(cer_scan,"Total","l");
  leg_Cher_sum->AddEntry(cer_scan_e,"Electrons","l");
  leg_Cher_sum->AddEntry(cer_scan_pi,"Pions","l");






  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Calculate Cherenkov electron effeciency and pion rejection effeciency and optimises ratio e/(1-pi) where e is the electron effeciency, pi the pion rejection effenciency and plots these values for different values of Cherenkov cut



  c50->cd(6);
  
  

  // number of points to sample effeciencies over
  Int_t No_points = 1000;


  // variables to hold electron, pion effeciencies, ratip of electrons to pions and cherenkov cut level 

  Double_t e_eff[No_points];
  Double_t pi_eff[No_points];
  Double_t eff_ratio[No_points];
  Double_t cer_cut_level[No_points];




  // some pions are only cut at very great levels of cherenkov sum - these are almost certainly not real pions but some other particle that has happened to pass the god electron cuts and pion targeted PRL cuts
  // solution: when increasing cut level if level of pions remains constant for X )cer_eff_lim) or more cuts then take this as final level of pions and final effeciency
  Double_t cer_eff_freeze = 0;
  Int_t cer_eff_count = 0;
  Int_t cer_eff_lim = 10;
  Bool_t freeze = false;



  // create histograms to store electron effeciency, pion rejection effeciency and ratio 

  TH1D* e_eff_h = new TH1D("e_eff_h","Effeciency ratios",1000,0,10000);
  TH1D* pi_eff_h = new TH1D("pi_eff_h","",1000,0,10000);
  TH1D* eff_ratio_h = new TH1D("eff_ratio_h","",1000,0,10000);
  




  for(Int_t j = 0; j<No_points; j++){

    e_eff[j] = 0.0;
    pi_eff[j] = 0.0;
    eff_ratio[j] = 0.0;

  }


  // varaible to hold level of cut

  Double_t sum_cut = 0;


  // define variables to hold the largest electron to pion ratio and the cherenkov cut that produces this

  Double_t ch_largest_ratio = 0;
  Double_t ch_lr_cut = 0;
  

  // loop through all values of cut and calculate effeciencies

  for(Int_t isum = 0; isum<No_points; isum++){
    
    sum_cut = (isum+1) * (10000 / No_points);
    
    cer_cut_level[isum] = sum_cut;
    

    // calculation of electron effecicency
    e_eff[isum] = cer_scan_e->Integral(cer_scan_e->FindBin(sum_cut),cer_scan_e->FindBin(10000))/cer_scan_e->Integral(cer_scan_e->FindBin(0),cer_scan_e->FindBin(10000));

    // calculation of pion rejection effecicency
    pi_eff[isum] =  1.0 - cer_scan_pi->Integral(cer_scan_pi->FindBin(sum_cut),cer_scan_pi->FindBin(10000))/cer_scan_pi->Integral(cer_scan_pi->FindBin(0),cer_scan_pi->FindBin(10000));
    
    

    
    // Loops to check if effeciency is static:
    // if true 'freeze' effeciency.
    // This is designed so that small number of 
    // 'pions' in sample that survive larger cuts are 
    // not included in pion rejection effeciency 
    // (as they are not real pions).
    
    if(pi_eff[isum] == pi_eff[isum-1] && !freeze){
      cer_eff_count++;
    }
    else{
      cer_eff_count = 0;
    }


    if(cer_eff_count >= cer_eff_lim && !freeze){
      cer_eff_freeze = pi_eff[isum];
      freeze = true;
      cout << endl << "Frozen at cut = " << isum << " with pi_Eff = " << pi_eff[isum] << endl << endl;
    }
    

    if(freeze){
      pi_eff[isum] = cer_eff_freeze;
    }
    


    // calculate ratio of electrons to pions
    eff_ratio[isum] = Double_t(e_eff[isum]/(1.0 -  pi_eff[isum]));


 
    // fill effeciency histograms
    e_eff_h->Fill(sum_cut,e_eff[isum]);
    pi_eff_h->Fill(sum_cut,pi_eff[isum]);
    eff_ratio_h->Fill(sum_cut,eff_ratio[isum]);



    // test if new ratio is largest 
    if(eff_ratio[isum]>ch_largest_ratio){

      ch_largest_ratio = eff_ratio[isum];
      ch_lr_cut = sum_cut;
      

    }

    
  }


  // set properties of histograms

  e_eff_h->SetLineColor(kRed);
  e_eff_h->SetMarkerColor(kRed);
  e_eff_h->SetMarkerStyle(kFullCircle);
  

  e_eff_h->GetXaxis()->SetTitle("Cherenkov Sum Cut");
  e_eff_h->GetXaxis()->CenterTitle();
  e_eff_h->GetYaxis()->SetTitle("Effeciency");
  e_eff_h->GetYaxis()->CenterTitle();


  e_eff_h->Draw("hist ");

  pi_eff_h->SetMarkerColor(kBlue);
  pi_eff_h->SetMarkerStyle(kFullTriangleUp);
  pi_eff_h->SetLineColor(kBlue);
  pi_eff_h->Draw("hist  same");




  // draw optimal cut

  TLine* ch_eff_line = new TLine(ch_lr_cut,.0,ch_lr_cut,gPad->GetUymax()*1.05);
  ch_eff_line->SetLineColor(kMagenta);
  ch_eff_line->SetLineWidth(2);

  ch_eff_line->Draw("same");


  
  // plot ratio on same canvas but with different y-scale (effeciencies are from 0-1 but ratio of electrons to pions can be much (orders of magnitude) greater)

  
  Double_t right_max = 1.1*eff_ratio_h->GetMaximum();
  Double_t ratio_scale = gPad->GetUymax()/right_max;
  
  eff_ratio_h->SetMarkerColor(kGreen+3);
  eff_ratio_h->SetMarkerStyle(kFullSquare);
  eff_ratio_h->Scale(ratio_scale);
  eff_ratio_h->SetLineColor(kGreen+3);
  eff_ratio_h->Draw("same hist ");

  
  // draw seperate axis for ratio

  TGaxis*axis = new TGaxis(10000,gPad->GetUymin(),10000,gPad->GetUymax(),0,right_max,510,"+L");
  
  cout << "  gPad->GetUxmax() = " << gPad->GetUxmax() << endl;

  axis->SetLineColor(kGreen+3);
  axis->SetLabelColor(kGreen+3);

  axis->SetTitle("Ratio (e:pi)");
  axis->CenterTitle();

  axis->SetTitleColor(kGreen+3);
  axis->Draw();

  



  TLegend* leg_Cher_eff = new TLegend(.20,.2,.50,.4,"Key");

  
  leg_Cher_eff->SetFillColor(0);
  leg_Cher_eff->AddEntry(e_eff_h,"\\epsilon_{e}","l");
  leg_Cher_eff->AddEntry(pi_eff_h,"\\epsilon_{\\pi}","l");
  leg_Cher_eff->AddEntry(eff_ratio_h,"\\epsilon_{e} /(1 - \\epsilon_{\\pi})","l");
  leg_Cher_eff->AddEntry(ch_eff_line,"Optimal cut","l");

  leg_Cher_eff->Draw("same");




  //~~~~~~~~~~~~~~~
  // Draw cut line for Cherenkov sum line

  c50->cd(5);

  TLine* ch_eff_line_cp = new TLine(ch_lr_cut,.0,ch_lr_cut,1.2*(cer_scan->GetMaximum()));


  ch_eff_line_cp->SetLineColor(kMagenta);
  ch_eff_line_cp->SetLineWidth(2);
  ch_eff_line_cp->Draw("same");

  leg_Cher_sum->AddEntry(ch_eff_line,"Optimal cut","l");
  leg_Cher_sum->Draw("same"); // draw legend here so as to appear in front of cut



  // save canvas


  c50->SaveAs("plots/L_Cherenkov_scan.pdf");
   
  
  gSystem->Exec("convert -density 700 -trim plots/L_Cherenkov_scan.pdf plots/L_Cherenkov_scan.png");







  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~60
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~60

  // Plots for Calorimeter cut scan 
  // 
  // by selecting 'clean' samples of pions and electrons from the Cherenkov, hope to optimise the level of cuts for the Calorimeters such that the level of electrons to pions is optimised

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~60
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~60



  TCanvas* c60 = new TCanvas("c60","Calorimeter Cut Scan",1200,1200);

  c60->Divide(3,2);



  // first plot cherenkov and relevant cuts usd for selecting 'good' electrons and pions for the Calorimeter cuts

  c60->cd(1);
  gPad->SetLogy();

  T->Draw("R.cer.asum_c>>cer_scan",gen_cut,"");

  
  // levels of electron cut (sampling from cherenkov)
  
  Double_t e_cer_l = 3000;
  Double_t e_cer_h = 9000;

  TLine* e_cer_lline = new TLine(e_cer_l, cer_scan->GetMinimum(), e_cer_l,cer_scan->GetMaximum());
  
  e_cer_lline->SetLineColor(kRed);
  
  e_cer_lline->Draw("same");
  
  
  TLine* e_cer_hline = new TLine(e_cer_h, cer_scan->GetMinimum(), e_cer_h,cer_scan->GetMaximum());

  e_cer_hline->SetLineColor(kRed);

  e_cer_hline->Draw("same");


  // levels of pion cut (sampling from cherenkov)
 
  Double_t pi_cer_l = 5;
  Double_t pi_cer_h = 400;
  
  TLine* pi_cer_lline = new TLine(pi_cer_l, cer_scan->GetMinimum(), pi_cer_l,cer_scan->GetMaximum());
  
  pi_cer_lline->SetLineColor(kBlue);

  pi_cer_lline->Draw("same");

  
  TLine* pi_cer_hline = new TLine(pi_cer_h, cer_scan->GetMinimum(), pi_cer_h,cer_scan->GetMaximum());

  pi_cer_hline->SetLineColor(kBlue);

  pi_cer_hline->Draw("same");


  leg_Cher_scan->Draw("same");
  
  



  // Plot of PS with electron sample and pion sample
  // later add determined cut


  c60->cd(2);

  gPad->SetLogy();
  

  Cal1_EP_plot->Draw();

  TH1D* Cal1_EP_plot_e = (TH1D*)Cal1_EP_plot->Clone("Cal1_EP_plot_e");



  Cal1_EP_plot_e->SetLineColor(kRed);
  
  Cal1_EP_plot_e->SetTitle("PS Energy");
  


  // set cherenkov-determined electron cut and plot electrons 

  TCut e_cal_scan_cut = gen_cut + Form("R.cer.asum_c > %f && R.cer.asum_c < %f", e_cer_l, e_cer_h);

  T->Draw("(R.ps.e)/(R.gold.p*1000)>>Cal1_EP_plot_e",e_cal_scan_cut,"same");


  
  

  // set cherenkov-determined pion cut and plot pions

  
  TH1D* Cal1_EP_plot_pi = (TH1D*)Cal1_EP_plot->Clone("Cal1_EP_plot_pi");


  Cal1_EP_plot_pi->SetLineColor(kBlue);

  TCut pi_cal_scan_cut = gen_cut + Form("R.cer.asum_c > %f && R.cer.asum_c < %f", pi_cer_l, pi_cer_h);

  T->Draw("(R.ps.e)/(R.gold.p*1000)>>Cal1_EP_plot_pi",pi_cal_scan_cut,"same");

  //attempt to fit gaussians to very small peak around zero and 'true' muon peak

  // TF1* PS_f_zero = new TF1("PS_f_zero",peak,-0.05,0.05,3);

  // PS_f_zero->SetParameter(0,500);
  // PS_f_zero->SetParameter(1,0.01);
  // PS_f_zero->SetParLimits(1,-0.0001,0.06);
  
  // PS_f_zero->SetParameter(2,0.01);
  
  // Cal1_EP_plot_pi->Fit(PS_f_zero,"R+");

  TF1* PS_fit = new TF1("PS_fit",overall,0,0.4,8);

  // PS_fit->SetParameter(0,500);
  PS_fit->SetParameter(1,0.01);
  // PS_fit->SetParLimits(1,-0.0001,0.06);
   PS_fit->SetParameter(2,0.01);


  // PS_fit->SetParameter(3,1000);
  PS_fit->SetParameter(4,0.10);
  // PS_fit->SetParLimits(4,0.09,0.12);
  PS_fit->SetParameter(5,0.03);

  PS_fit->SetParameter(6,200);
  //  PS_fit->SetParameter(7,1000);

  Cal1_EP_plot_pi->Fit(PS_fit,"R0+");

  // retrieve Pion peak and background from overall fit
  
  TF1* PS_Pi_fit = new TF1("PS_Pi_fit",peak,0,0.4,3);

  PS_Pi_fit->SetParameter(0,PS_fit->GetParameter(3));
  PS_Pi_fit->SetParameter(1,PS_fit->GetParameter(4));
  PS_Pi_fit->SetParameter(2,PS_fit->GetParameter(5));

  PS_Pi_fit->SetLineColor(kGreen+2);
  PS_Pi_fit->SetLineStyle(9);
  
  PS_Pi_fit->Draw("same");

  
  TF1* PS_bg_fit = new TF1("PS_bg_fit",bg,0,0.4,2);

  PS_bg_fit->SetParameter(0,PS_fit->GetParameter(6));
  PS_bg_fit->SetParameter(1,PS_fit->GetParameter(7));

  PS_bg_fit->SetLineColor(kOrange+2);
  PS_bg_fit->SetLineStyle(9);
  PS_bg_fit->Draw("same");
  
  // Int_t nopar = PS_f_zero->GetNumberFreeParameters();
  // Double_t* pi_pars = PS_f_zero->GetParameters();

  // for(Int_t i = 0; i<nopar; i++){
  //   cout << "parameter " << i << " =  " << pi_pars[i] << endl;
  // }
  
  //  PS_f_zero->Draw("same");

  // Effeciency plot for PS energy cut


  c60->cd(3);

  // set-up variables

  Double_t PS_cut = 0;
  Double_t PS_cut_level[No_points];

  // define variables to hold the largest electron to pion ratio and the cherenkov cut that produces this
  // set effecieny variables back to zero


  for(Int_t j = 0; j<No_points; j++){
    e_eff[j] = 0.0;
    pi_eff[j] = 0.0;
    eff_ratio[j] = 0.0;     
    PS_cut_level[No_points];
  }



  // create histograms to store electron effeciency, pion rejection effecieny and ratio 
  
  // upper limit of PS cut
  Double_t PS_lim = 1;


  TH1D* PS_e_eff_h = new TH1D("PS_e_eff_h","Effeciency ratios",1000,0,PS_lim);
  TH1D* PS_pi_eff_h = new TH1D("PS_pi_eff_h","",1000,0,PS_lim);
  TH1D* PS_eff_ratio_h = new TH1D("PS_eff_ratio_h","",1000,0,PS_lim);



  // define variables to hold the largest electron to pion ratio and the PS cut that produces this

  Double_t PS_largest_ratio = 0;
  Double_t PS_lr_cut = 0;
  


  // loop through values PS cut and calculate effeciencies
  for(Int_t isum = 0; isum<No_points; isum++){
    

    // define PRL cut as going from to 1
    sum_cut = (isum+1) * (PS_lim / No_points);
    
    cer_cut_level[isum] = sum_cut;


    // electron effeciency calculation
    e_eff[isum] =  Cal1_EP_plot_e->Integral(Cal1_EP_plot_e->FindBin(sum_cut),Cal1_EP_plot_e->FindBin(PS_lim))/Cal1_EP_plot_e->Integral(Cal1_EP_plot_e->FindBin(0),Cal1_EP_plot_e->FindBin(PS_lim));
    
    PS_e_eff_h->Fill(sum_cut,e_eff[isum]);


    // pion rejection effeciency calculation

    pi_eff[isum] =  1.0 - Cal1_EP_plot_pi->Integral(Cal1_EP_plot_pi->FindBin(sum_cut),Cal1_EP_plot_pi->FindBin(PS_lim))/Cal1_EP_plot_pi->Integral(Cal1_EP_plot_pi->FindBin(0),Cal1_EP_plot_pi->FindBin(PS_lim));

    PS_pi_eff_h->Fill(sum_cut,pi_eff[isum]);


    // electron pion ratio calculation

    
    eff_ratio[isum] = Double_t(e_eff[isum]/(1.0 -  pi_eff[isum]));


    PS_eff_ratio_h->Fill(sum_cut,eff_ratio[isum]);

    
    // test if new ratio is largest 
    if(eff_ratio[isum]>PS_largest_ratio){
      
      PS_largest_ratio = eff_ratio[isum];
      PS_lr_cut = sum_cut;
    }





  }

  

  // draw histograms

  PS_e_eff_h->SetLineColor(kRed);
  
  PS_e_eff_h->GetXaxis()->SetTitle("PS cut level");
  PS_e_eff_h->GetXaxis()->CenterTitle();
  PS_e_eff_h->GetYaxis()->SetTitle("Effeciency");
  PS_e_eff_h->GetYaxis()->CenterTitle();
  PS_e_eff_h->Draw("hist");

  
  PS_pi_eff_h->SetLineColor(kBlue);
  PS_pi_eff_h->Draw("hist  same");


  // plot ratio on same canvas but with different y-scale (effeciencies are from 0-1 but ratio of electrons to pions can be much (orders of magnitude) greater)

  right_max = 1.1*PS_eff_ratio_h->GetMaximum();
  ratio_scale = gPad->GetUymax()/right_max;
  
  PS_eff_ratio_h->SetMarkerColor(kGreen+3);
  PS_eff_ratio_h->SetMarkerStyle(kFullSquare);
  PS_eff_ratio_h->Scale(ratio_scale);
  PS_eff_ratio_h->SetLineColor(kGreen+3);
  PS_eff_ratio_h->Draw("same hist ");


  // draw seperate axis for ratio
  
  cout << "right_max for PS = " << right_max << endl;

  TGaxis*PSaxis = new TGaxis(PS_lim,gPad->GetUymin(),PS_lim,gPad->GetUymax(),0,right_max,510,"+L");
  
  cout << "  gPad->GetUxmax() = " << gPad->GetUxmax() << endl;

  PSaxis->SetLineColor(kGreen+3);
  PSaxis->SetLabelColor(kGreen+3);

  PSaxis->SetTitle("Ratio (e:pi)");
  PSaxis->CenterTitle();
  ///  PSaxis->Rotate();
  PSaxis->SetTitleColor(kGreen+3);
  PSaxis->Draw();


  
  // later add determined cut




  

  // draw PS + SH energy (total, electrons sampled from cherenkov and pions sampled from cherenkov)


  c60->cd(4);

  gPad->SetLogy();
  
  // draw total

  Cal_EP_plot->Draw();



  // plot electrons 

  TH1D* Cal1_2_EP_plot_e = (TH1D*)Cal_EP_plot->Clone("Cal1_2_EP_plot_e");


  Cal1_2_EP_plot_e->SetLineColor(kRed);
  
  Cal1_2_EP_plot_e->SetTitle("PS + SH Energy");

  T->Draw("(R.ps.e+R.sh.e)/(R.gold.p*1000)>>Cal1_2_EP_plot_e",e_cal_scan_cut,"same");



  //  plot pions

  
  TH1D* Cal1_2_EP_plot_pi = (TH1D*)Cal_EP_plot->Clone("Cal1_2_EP_plot_pi");


  Cal1_2_EP_plot_pi->SetLineColor(kBlue);


  T->Draw("(R.ps.e+R.sh.e)/(R.gold.p*1000)>>Cal1_2_EP_plot_pi",pi_cal_scan_cut,"same");
  




  // Effeciency plot for PS + SH energy cut


  c60->cd(5);

  // set-up variables

  Double_t PS_2_cut = 0;
  Double_t PS_2_cut_level[No_points];

  // define variables to hold the largest electron to pion ratio and the cherenkov cut that produces this
  // set effecieny variables back to zero


  for(Int_t j = 0; j<No_points; j++){
    e_eff[j] = 0.0;
    pi_eff[j] = 0.0;
    eff_ratio[j] = 0.0;     
    PS_2_cut_level[No_points];
  }



  
  // create histograms to store electron effeciency, pion rejection effecieny and ratio 
  
  // upper limit of PS_2 cut
  Double_t PS_2_lim = 1;


  TH1D* PS_2_e_eff_h = new TH1D("PS_2_e_eff_h","Effeciency ratios",1000,0,PS_2_lim);
  TH1D* PS_2_pi_eff_h = new TH1D("PS_2_pi_eff_h","",1000,0,PS_2_lim);
  TH1D* PS_2_eff_ratio_h = new TH1D("PS_2_eff_ratio_h","",1000,0,PS_2_lim);



  
  // define variables to hold the largest electron to pion ratio and the PS_2 cut that produces this

  Double_t PS_2_largest_ratio = 0;
  Double_t PS_2_lr_cut = 0;
  

  // loop through values PS + SH cut and calculate effeciencies
  for(Int_t isum = 0; isum<No_points; isum++){
    

    // define PRL cut as going from to 1
    sum_cut = (isum+1) * (PS_2_lim / No_points);
    
    cer_cut_level[isum] = sum_cut;


    // electron effeciency calculation
    e_eff[isum] =  Cal1_2_EP_plot_e->Integral(Cal1_2_EP_plot_e->FindBin(sum_cut),Cal1_2_EP_plot_e->FindBin(PS_2_lim))/Cal1_2_EP_plot_e->Integral(Cal1_2_EP_plot_e->FindBin(0),Cal1_2_EP_plot_e->FindBin(PS_2_lim));
    
    PS_2_e_eff_h->Fill(sum_cut,e_eff[isum]);


    // pion rejection effeciency calculation

    pi_eff[isum] =  1.0 - Cal1_2_EP_plot_pi->Integral(Cal1_2_EP_plot_pi->FindBin(sum_cut),Cal1_2_EP_plot_pi->FindBin(PS_2_lim))/Cal1_2_EP_plot_pi->Integral(Cal1_2_EP_plot_pi->FindBin(0),Cal1_2_EP_plot_pi->FindBin(PS_2_lim));

    PS_2_pi_eff_h->Fill(sum_cut,pi_eff[isum]);


    // electron pion ratio calculation

    
    eff_ratio[isum] = Double_t(e_eff[isum]/(1.0 -  pi_eff[isum]));


    PS_2_eff_ratio_h->Fill(sum_cut,eff_ratio[isum]);



    
    // test if new ratio is largest 
    if(eff_ratio[isum]>PS_2_largest_ratio){
      
      PS_2_largest_ratio = eff_ratio[isum];
      PS_2_lr_cut = sum_cut;
    }





  }
  

  // draw histograms

  PS_2_e_eff_h->SetLineColor(kRed);
  
  PS_2_e_eff_h->GetXaxis()->SetTitle("PS + SH cut level");
  PS_2_e_eff_h->GetXaxis()->CenterTitle();
  PS_2_e_eff_h->GetYaxis()->SetTitle("Effeciency");
  PS_2_e_eff_h->GetYaxis()->CenterTitle();
  PS_2_e_eff_h->Draw("hist");

  
  PS_2_pi_eff_h->SetLineColor(kBlue);
  PS_2_pi_eff_h->Draw("hist  same");


  // plot ratio on same canvas but with different y-scale (effeciencies are from 0-1 but ratio of electrons to pions can be much (orders of magnitude) greater)

  right_max = 1.1*PS_2_eff_ratio_h->GetMaximum();
  ratio_scale = gPad->GetUymax()/right_max;
  
  PS_2_eff_ratio_h->SetMarkerColor(kGreen+3);
  PS_2_eff_ratio_h->SetMarkerStyle(kFullSquare);
  PS_2_eff_ratio_h->Scale(ratio_scale);
  PS_2_eff_ratio_h->SetLineColor(kGreen+3);
  PS_2_eff_ratio_h->Draw("same hist ");


  // draw seperate axis for ratio
  
  cout << "right_max for PS_2 = " << right_max << endl;

  TGaxis*PS_2axis = new TGaxis(PS_2_lim,gPad->GetUymin(),PS_2_lim,gPad->GetUymax(),0,right_max,510,"+L");
  
  cout << "  gPad->GetUxmax() = " << gPad->GetUxmax() << endl;

  PS_2axis->SetLineColor(kGreen+3);
  PS_2axis->SetLabelColor(kGreen+3);

  PS_2axis->SetTitle("Ratio (e:pi)");
  PS_2axis->CenterTitle();
  ///  PS_2axis->Rotate();
  PS_2axis->SetTitleColor(kGreen+3);
  PS_2axis->Draw();


  leg_Cher_eff->Draw("same");

  // later add determined cut




  // Plot of PS vs SH with cuts determined previously in canvas


  c60->cd(6);
  
  TH2D* Cal1_Vs_Cal2_EP_plot = new TH2D("Cal1_Vs_Cal2_EP_plot","",1000,-0.05,1.2,1000,-0.05,1.2);



  Cal1_Vs_Cal2_EP_plot->GetXaxis()->SetTitle("(PS.e)/R.gold.p");
  Cal1_Vs_Cal2_EP_plot->GetXaxis()->CenterTitle();
  
  Cal1_Vs_Cal2_EP_plot->GetYaxis()->SetTitle("(SH.e)/R.gold.p");
  Cal1_Vs_Cal2_EP_plot->GetYaxis()->CenterTitle();
 
  T->Draw("(R.sh.e)/(R.gold.p*1000):(R.ps.e)/(R.gold.p*1000)>>Cal1_Vs_Cal2_EP_plot",gen_cut,"col");
    



  // write out level of PS cut(s)
 
  Double_t PS_cut_l = 0.2;
 
   
  TLine* PS_fin_cut = new TLine(PS_cut_l,0,PS_cut_l,Cal1_EP_plot_e->GetMaximum());
  PS_fin_cut->SetLineColor(kMagenta);
  PS_fin_cut->SetLineWidth(2);
 

  c60->cd(2);
  PS_fin_cut->Draw("same");

  // draw legend
  leg_Cher_sum->Draw("same");



  c60->cd(3);
  
  TLine* PS_fin_cut_eff = new TLine(PS_cut_l,0,PS_cut_l,1);
  PS_fin_cut_eff->SetLineColor(kMagenta);
  PS_fin_cut_eff->SetLineWidth(2);
  

  PS_fin_cut_eff->Draw("same");
  
  leg_Cher_eff->Draw("same");





  // write out level of PS + SH cut

  Double_t PS_2_cut_l = 0.51;
  Double_t PS_2_cut_h = 1.25;
  

  
  TLine* PS_2_fin_cut = new TLine(PS_2_cut_l,0,PS_2_cut_l,Cal1_2_EP_plot_e->GetMaximum());
  PS_2_fin_cut->SetLineColor(kMagenta);
  PS_2_fin_cut->SetLineWidth(2);
  
  c60->cd(4);

  PS_2_fin_cut->Draw("same");

  // draw legend
  leg_Cher_sum->Draw("same");

  c60->cd(5);


  TLine* PS_2_fin_cut_cp = new TLine(PS_2_cut_l,0,PS_2_cut_l,1);
  PS_2_fin_cut_cp->SetLineColor(kMagenta);
  PS_2_fin_cut_cp->SetLineWidth(2);
  
  PS_2_fin_cut_cp->Draw("same");
  
  leg_Cher_eff->Draw("same");


  


  c60->cd(6);

  // draw diagonal cut on SH vs PS plot

  TLine* PS_2_fin_cut_diag_l = new TLine(0,PS_2_cut_l,PS_2_cut_l,0);

  
  PS_2_fin_cut_diag_l->SetLineColor(kMagenta);
  PS_2_fin_cut_diag_l->SetLineWidth(2);
  PS_2_fin_cut_diag_l->Draw("same");

  
  TLine* PS_2_fin_cut_diag_h = new TLine(0,PS_2_cut_h,PS_2_cut_h,0);

  
  PS_2_fin_cut_diag_h->SetLineColor(kMagenta);
  PS_2_fin_cut_diag_h->SetLineWidth(2);
  PS_2_fin_cut_diag_h->Draw("same");
  
  
  PS_fin_cut_eff->Draw("same");


  c60->SaveAs("plots/L_Calorimeter_scan.pdf");
     
  gSystem->Exec("convert -density 700 -trim plots/L_Calorimeter_scan.pdf plots/L_Calorimeter_scan.png");


  gSystem->Exec("pdfunite plots/L_event_selection.pdf plots/L_track_var_plots.pdf plots/L_Cherenkov_scan.pdf plots/L_Calorimeter_scan.pdf plots/L_PID_plots_combined.pdf");
  


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Finally print values of cuts to terminal


  cout << "RHRS final cuts are: " << endl;
  cout << "gen_cut: " << gen_cut << endl << endl;
  cout << "PID cuts " << Form("(R.ps.e/(R.gold.p*1000)) > %f  && ((R.sh.e + R.ps.e)/(R.gold.p*1000)) > %f && ((R.sh.e + R.ps.e)/(R.gold.p*1000)) < %f && R.cer.asum_c > %f",PS_cut_l,PS_2_cut_l,PS_2_cut_h,ch_lr_cut) << endl << endl;
  cout << "Track-detector cuts: " << Form("R.s0.trx>%f && R.s0.trx<%f && R.s0.try>%f && R.s0.try<%f && R.s2.trx>%f && R.s2.trx<%f && R.s2.try>%f && R.s2.try<%f && R.ps.trx>%f && R.ps.trx<%f && R.ps.try>%f && R.ps.try<%f && R.sh.trx>%f && R.sh.trx<%f && R.sh.try>%f && R.sh.try<%f",-0.5*s0_x_height,0.5*s0_x_height,-0.5*s0_y_width,0.5*s0_y_width,-0.5*s2_x_height,0.5*s2_x_height,-0.5*s2_y_width,0.5*s2_y_width,-0.5*ps_x_height,0.5*ps_x_height,-0.5*ps_y_width,0.5*ps_y_width,-0.5*sh_x_height,0.5*sh_x_height,-0.5*sh_y_width,0.5*sh_y_width) << endl << endl;
  cout << "Track-Calorimeter agreement cuts" <<  Form("(R.ps.x-R.ps.trx)>%f && (R.ps.x-R.ps.trx)<%f && (R.sh.x-R.sh.trx)>%f && (R.sh.x-R.sh.trx)<%f",ps_E_min_Tx_l,ps_E_min_Tx_h,sh_E_min_Tx_l,sh_E_min_Tx_h) << endl << endl;
  
  


}



  



