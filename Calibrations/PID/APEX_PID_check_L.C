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

void APEX_PID_check_L()
{


  



  Int_t run_number;    

  cout << "enter run_number: ";
  cin >> run_number;

  
  // TChain* T_2 = Load_more_rootfiles(run_number);
  TChain* T = Load_more_rootfiles(run_number);
  //
  // TTree* T = T_2->CloneTree(1e4);
  





  cout << T->GetEntries() << endl;


  



 
  // Left-arm PID plots


  // left cut 

  TCut lcut = "L.tr.n==1";
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


  
  Double_t e_cut_cal_h = 1.05;
  Double_t e_cut_cal_l = 0.75;
  
  Double_t e_cut_cer_h  = 8500;
  Double_t e_cut_cer_l  = 2000;


  Double_t pi_cut_cal_h = 0.33;
  Double_t pi_cut_cal_l = 0.165;
  
  Double_t pi_cut_cer_h  = 200;
  Double_t pi_cut_cer_l  = 0.0;


  Double_t zer_cut_cal_h = 0.13;
  Double_t zer_cut_cal_l = -0.03;
  
  Double_t zer_cut_cer_h  = 200;
  Double_t zer_cut_cer_l  = 0.0;


  Double_t bg_cut_cal_h = 0.15;
  Double_t bg_cut_cal_l = -0.03;
  
  Double_t bg_cut_cer_h  = 6200;
  Double_t bg_cut_cer_l  = 1750;



  // Cherenkov sum vs Calorimeter energy plot

  TH2D *ChVsEP = new TH2D("ChVsEP","LHRS PID Plot",1000,-0.1,1.4,1000,-10,10000);

  

  ChVsEP->GetXaxis()->SetTitle("E_{tot}/P");
  ChVsEP->GetXaxis()->CenterTitle();
  ChVsEP->GetXaxis()->SetTitleSize(0.048);
  ChVsEP->GetXaxis()->SetTitleOffset(0.98);
  
  
  
  ChVsEP->GetYaxis()->SetTitle("Cherenkov sum");
  ChVsEP->GetYaxis()->CenterTitle();
  ChVsEP->GetYaxis()->SetTitleSize(0.048);
  ChVsEP->GetYaxis()->SetTitleOffset(0.85);


  
  T->Draw("L.cer.asum_c:((L.prl1.e+L.prl2.e)/(L.gold.p*1000))>>ChVsEP",lcut,"col");


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


  // set-up L.gold.y (y_tg of golden track) plots for all 4 categories of events


  TCanvas* c2 = new TCanvas("c2","Track and Reconstruction variables",1200,1200);

  c2->Divide(2,2);
  
  // set-up L.gold.th (theta angle of golden track) plots for all 4 categories of events

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
  TCut e_cut = Form("L.cer.asum_c > %f &&  L.cer.asum_c < %f && ((L.prl1.e+L.prl2.e)/(L.gold.p*1000)) > %f && ((L.prl1.e+L.prl2.e)/(L.gold.p*1000)) < %f",e_cut_cer_l,e_cut_cer_h,e_cut_cal_l,e_cut_cal_h);
  
  
  T->Draw("L.gold.th>>Theta_plot_e",e_cut);


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

  TCut pi_cut = Form("L.cer.asum_c > %f &&  L.cer.asum_c < %f && ((L.prl1.e+L.prl2.e)/(L.gold.p*1000)) > %f && ((L.prl1.e+L.prl2.e)/(L.gold.p*1000)) < %f",pi_cut_cer_l,pi_cut_cer_h,pi_cut_cal_l,pi_cut_cal_h);
  
  T->Draw("L.gold.th>>Theta_plot_pi",pi_cut,"same");



  // no hits theta plot

  TH1D* Theta_plot_zer = (TH1D*)Theta_plot_e->Clone("Theta_plot_zer");

  TCut zer_cut = Form("L.cer.asum_c > %f &&  L.cer.asum_c < %f && ((L.prl1.e+L.prl2.e)/(L.gold.p*1000)) > %f && ((L.prl1.e+L.prl2.e)/(L.gold.p*1000)) < %f",zer_cut_cer_l,zer_cut_cer_h,zer_cut_cal_l,zer_cut_cal_h);

  Theta_plot_zer->SetLineColor(kBlack);
  T->Draw("L.gold.th>>Theta_plot_zer",zer_cut,"same");


  // background theta plot

  TH1D* Theta_plot_bg = (TH1D*)Theta_plot_e->Clone("Theta_plot_bg");


  TCut bg_cut = Form("L.cer.asum_c > %f &&  L.cer.asum_c < %f && ((L.prl1.e+L.prl2.e)/(L.gold.p*1000)) > %f && ((L.prl1.e+L.prl2.e)/(L.gold.p*1000)) < %f",bg_cut_cer_l,bg_cut_cer_h,bg_cut_cal_l,bg_cut_cal_h);

  
  Theta_plot_bg->SetLineColor(kGreen);
  T->Draw("L.gold.th>>Theta_plot_bg",bg_cut,"same");


  cout << "Number of Entries in \\theta plots: " << endl;
  cout << "Electrons " << Theta_plot_e->GetEntries() << endl;
  cout << "Pions " << Theta_plot_pi->GetEntries() << endl;
  cout << "Zeros " << Theta_plot_zer->GetEntries() << endl;
  cout << "Bg " << Theta_plot_bg->GetEntries() << endl;
  



  // set-up L.gold.ph (phi angle of golden track) plots for all 4 categories of events

  c2->cd(2);
  gPad->SetLogy();

  TH1D* Phi_plot_e =  new TH1D("Phi_plot_e","",100,-0.2,0.2);

  Phi_plot_e->GetXaxis()->SetTitle("\\phi_{tg} [rad]");
  Phi_plot_e->GetXaxis()->CenterTitle();
  Phi_plot_e->GetYaxis()->SetTitle("Entries");
  Phi_plot_e->GetYaxis()->CenterTitle();
  
  
  Phi_plot_e->SetLineColor(kRed);
  Phi_plot_e->SetLineWidth(2);;

  T->Draw("L.gold.ph>>Phi_plot_e",e_cut);




  // phi cut (used later on to select 'good' electrons)

  Double_t phi_cut_l = -0.05;
  Double_t phi_cut_h = 0.05;

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
  T->Draw("L.gold.ph>>Phi_plot_pi",pi_cut,"same");



  // no hits phi plot

  TH1D* Phi_plot_zer = (TH1D*)Phi_plot_e->Clone("Phi_plot_zer");


  Phi_plot_zer->SetLineColor(kBlack);
  T->Draw("L.gold.ph>>Phi_plot_zer",zer_cut,"same");


  // background phi plot

  TH1D* Phi_plot_bg = (TH1D*)Phi_plot_e->Clone("Phi_plot_bg");

  Phi_plot_bg->SetLineColor(kGreen);
  T->Draw("L.gold.ph>>Phi_plot_bg",bg_cut,"same");




  // set-up L.gold.dp (deviation from central momentum of golden track) plots for all 4 categories of events

  c2->cd(3);
  gPad->SetLogy();

  TH1D* Dp_plot_e =  new TH1D("Dp_plot_e","",100,-0.2,0.2);

  Dp_plot_e->GetXaxis()->SetTitle("\\delta p");
  Dp_plot_e->GetXaxis()->CenterTitle();
  Dp_plot_e->GetYaxis()->SetTitle("Entries");
  Dp_plot_e->GetYaxis()->CenterTitle();
  
  
  Dp_plot_e->SetLineColor(kRed);
  Dp_plot_e->SetLineWidth(2);
  
  T->Draw("L.gold.dp>>Dp_plot_e",e_cut);

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
  T->Draw("L.gold.dp>>Dp_plot_pi",pi_cut,"same");



  // no hits dp plot

  TH1D* Dp_plot_zer = (TH1D*)Dp_plot_e->Clone("Dp_plot_zer");


  Dp_plot_zer->SetLineColor(kBlack);
  T->Draw("L.gold.dp>>Dp_plot_zer",zer_cut,"same");


  // background dp plot

  TH1D* Dp_plot_bg = (TH1D*)Dp_plot_e->Clone("Dp_plot_bg");

  Dp_plot_bg->SetLineColor(kGreen);
  T->Draw("L.gold.dp>>Dp_plot_bg",bg_cut,"same");




  c2->cd(4);
  gPad->SetLogy();

  TH1D* Y_plot_e =  new TH1D("Y_plot_e","",100,-0.2,0.2);

  Y_plot_e->GetXaxis()->SetTitle("y_{tg} [rad]");
  Y_plot_e->GetXaxis()->CenterTitle();
  Y_plot_e->GetYaxis()->SetTitle("Entries");
  Y_plot_e->GetYaxis()->CenterTitle();
  
  
  Y_plot_e->SetLineColor(kRed);
  Y_plot_e->SetLineWidth(2);

  T->Draw("L.gold.y>>Y_plot_e",e_cut);


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
  T->Draw("L.gold.y>>Y_plot_pi",pi_cut,"same");



  // no hits y plot

  TH1D* Y_plot_zer = (TH1D*)Y_plot_e->Clone("Y_plot_zer");


  Y_plot_zer->SetLineColor(kBlack);
  T->Draw("L.gold.y>>Y_plot_zer",zer_cut,"same");


  // background y plot

  TH1D* Y_plot_bg = (TH1D*)Y_plot_e->Clone("Y_plot_bg");

  Y_plot_bg->SetLineColor(kGreen);
  T->Draw("L.gold.y>>Y_plot_bg",bg_cut,"same");


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
  
  
  T->Draw("L.gold.beta>>Beta_plot_e",e_cut);


 

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
  T->Draw("L.gold.beta>>Beta_plot_pi",pi_cut,"same");


  // no hits beta plot

  TH1D* Beta_plot_zer = (TH1D*)Beta_plot_e->Clone("Beta_plot_zer");


  Beta_plot_zer->SetLineColor(kBlack);
  T->Draw("L.gold.beta>>Beta_plot_zer",zer_cut,"same");


  // background beta plot

  TH1D* Beta_plot_bg = (TH1D*)Beta_plot_e->Clone("Beta_plot_bg");


  Beta_plot_bg->SetLineColor(kGreen);
  T->Draw("L.gold.beta>>Beta_plot_bg",bg_cut,"same");




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
  
  
  T->Draw("L.tr.chi2>>Chi_plot_e",e_cut);




  
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
  T->Draw("L.tr.chi2>>Chi_plot_pi",pi_cut,"same");


  // no hits chi2 plot

  TH1D* Chi_plot_zer = (TH1D*)Chi_plot_e->Clone("Chi_plot_zer");


  Chi_plot_zer->SetLineColor(kBlack);
  T->Draw("L.tr.chi2>>Chi_plot_zer",zer_cut,"same");


  // background chi2 plot

  TH1D* Chi_plot_bg = (TH1D*)Chi_plot_e->Clone("Chi_plot_bg");


  Chi_plot_bg->SetLineColor(kGreen);
  T->Draw("L.tr.chi2>>Chi_plot_bg",bg_cut,"same");




  //~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // From above plots on 2nd canvas create an over 'good' electron/general cut


  TCut gen_cut = Form("L.tr.n==1 && L.gold.beta>%f && L.tr.chi2<%f && L.gold.th>%f && L.gold.th <%f && L.gold.ph>%f && L.gold.ph <%f && L.gold.dp>%f && L.gold.dp<%f", beta_cut_l, chi2_cut_l, theta_cut_l, theta_cut_h, phi_cut_l, phi_cut_h, dp_cut_l, dp_cut_h);








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

    
  T->Draw("L.s0.trx:L.s0.try>>s0_plot_e",e_cut,"colz");

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
  
  s0_plot_pi->GetXaxis()->SetTitle("s0 X track projection [m]");
  s0_plot_pi->GetXaxis()->CenterTitle();
  s0_plot_pi->GetYaxis()->SetTitle("s0 Y track projection [m]");
  s0_plot_pi->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.s0.trx:L.s0.try>>s0_plot_pi",pi_cut,"colz");

  s0_box->Draw("same");


  c3->cd(3);
  // gPad->SetLogy();

  TH2D *s0_plot_zer = new TH2D("s0_plot_zer","Zero hits",200,s0_ymin,s0_ymax,200,s0_xmin,s0_xmax);
  
  s0_plot_zer->GetXaxis()->SetTitle("s0 X track projection [m]");
  s0_plot_zer->GetXaxis()->CenterTitle();
  s0_plot_zer->GetYaxis()->SetTitle("s0 Y track projection [m]");
  s0_plot_zer->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.s0.trx:L.s0.try>>s0_plot_zer",zer_cut,"colz");

  s0_box->Draw("same");




  c3->cd(4);
  // gPad->SetLogy();

  TH2D *s0_plot_bg = new TH2D("s0_plot_bg","Bg hits",200,s0_ymin,s0_ymax,200,s0_xmin,s0_xmax);
  
  s0_plot_bg->GetXaxis()->SetTitle("s0 X track projection [m]");
  s0_plot_bg->GetXaxis()->CenterTitle();
  s0_plot_bg->GetYaxis()->SetTitle("s0 Y track projection [m]");
  s0_plot_bg->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.s0.trx:L.s0.try>>s0_plot_bg",bg_cut,"colz");

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
  
  s2_plot_e->GetXaxis()->SetTitle("s2 X track projection [m]");
  s2_plot_e->GetXaxis()->CenterTitle();
  s2_plot_e->GetYaxis()->SetTitle("s2 Y track projection [m]");
  s2_plot_e->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.s2.trx:L.s2.try>>s2_plot_e",e_cut,"colz");

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
  
  s2_plot_pi->GetXaxis()->SetTitle("s2 X track projection [m]");
  s2_plot_pi->GetXaxis()->CenterTitle();
  s2_plot_pi->GetYaxis()->SetTitle("s2 Y track projection [m]");
  s2_plot_pi->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.s2.trx:L.s2.try>>s2_plot_pi",pi_cut,"colz");

  s2_box->Draw("same");


  c4->cd(3);
  // gPad->SetLogy();

  TH2D *s2_plot_zer = new TH2D("s2_plot_zer","Zero hits",200,s2_ymin,s2_ymax,200,s2_xmin,s2_xmax);
  
  s2_plot_zer->GetXaxis()->SetTitle("s2 X track projection [m]");
  s2_plot_zer->GetXaxis()->CenterTitle();
  s2_plot_zer->GetYaxis()->SetTitle("s2 Y track projection [m]");
  s2_plot_zer->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.s2.trx:L.s2.try>>s2_plot_zer",zer_cut,"colz");

  s2_box->Draw("same");




  c4->cd(4);
  // gPad->SetLogy();

  TH2D *s2_plot_bg = new TH2D("s2_plot_bg","Bg hits",200,s2_ymin,s2_ymax,200,s2_xmin,s2_xmax);
  
  s2_plot_bg->GetXaxis()->SetTitle("s2 X track projection [m]");
  s2_plot_bg->GetXaxis()->CenterTitle();
  s2_plot_bg->GetYaxis()->SetTitle("s2 Y track projection [m]");
  s2_plot_bg->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.s2.trx:L.s2.try>>s2_plot_bg",bg_cut,"colz");

  s2_box->Draw("same");


  // gPad->SetLogy();




  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~5
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~5

  // Plot track projection distribution in PRL1

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~5

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~5

  
  TCanvas* c5 = new TCanvas("c5","PRL1 variables",1200,1200);
    
  c5->Divide(2,2);
  c5->cd(1);
  // gPad->SetLogy();


  // set-up prl1 y-dim for al categories of events

  Double_t prl1_ymin = -0.4;
  Double_t prl1_ymax = 0.4;
  
  Double_t prl1_xmin = -1.5;
  Double_t prl1_xmax = 1.5;


  
  TH2D *prl1_plot_e = new TH2D("prl1_plot_e","Electrons",200,prl1_ymin,prl1_ymax,200,prl1_xmin,prl1_xmax);
  
  prl1_plot_e->GetXaxis()->SetTitle("prl1 X track projection [m]");
  prl1_plot_e->GetXaxis()->CenterTitle();
  prl1_plot_e->GetYaxis()->SetTitle("prl1 Y track projection [m]");
  prl1_plot_e->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl1.trx:L.prl1.try>>prl1_plot_e",e_cut,"colz");

  // draw box to represent prl1 dimensions

  Double_t  prl1_x_height = 2.6;
  Double_t  prl1_y_width = 0.7;

  TBox* prl1_box = new TBox(-prl1_y_width/2,-prl1_x_height/2,+prl1_y_width/2,+prl1_x_height/2);

  prl1_box->SetFillStyle(0);
  prl1_box->SetLineColor(kRed);


  prl1_box->Draw("same");


  c5->Update();

  c5->cd(2);
  // gPad->SetLogy();

  TH2D *prl1_plot_pi = new TH2D("prl1_plot_pi","Pions",200,prl1_ymin,prl1_ymax,200,prl1_xmin,prl1_xmax);
  
  prl1_plot_pi->GetXaxis()->SetTitle("prl1 X track projection [m]");
  prl1_plot_pi->GetXaxis()->CenterTitle();
  prl1_plot_pi->GetYaxis()->SetTitle("prl1 Y track projection [m]");
  prl1_plot_pi->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl1.trx:L.prl1.try>>prl1_plot_pi",pi_cut,"colz");

  prl1_box->Draw("same");


  c5->cd(3);
  // gPad->SetLogy();

  TH2D *prl1_plot_zer = new TH2D("prl1_plot_zer","Zero hits",200,prl1_ymin,prl1_ymax,200,prl1_xmin,prl1_xmax);
  
  prl1_plot_zer->GetXaxis()->SetTitle("prl1 X track projection [m]");
  prl1_plot_zer->GetXaxis()->CenterTitle();
  prl1_plot_zer->GetYaxis()->SetTitle("prl1 Y track projection [m]");
  prl1_plot_zer->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl1.trx:L.prl1.try>>prl1_plot_zer",zer_cut,"colz");

  prl1_box->Draw("same");




  c5->cd(4);
  // gPad->SetLogy();

  TH2D *prl1_plot_bg = new TH2D("prl1_plot_bg","Bg hits",200,prl1_ymin,prl1_ymax,200,prl1_xmin,prl1_xmax);
  
  prl1_plot_bg->GetXaxis()->SetTitle("prl1 X track projection [m]");
  prl1_plot_bg->GetXaxis()->CenterTitle();
  prl1_plot_bg->GetYaxis()->SetTitle("prl1 Y track projection [m]");
  prl1_plot_bg->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl1.trx:L.prl1.try>>prl1_plot_bg",bg_cut,"colz");

  prl1_box->Draw("same");

  

  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~6
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~6

  // Plot track projection distribution in PRL2

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~6

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~6

  
  TCanvas* c6 = new TCanvas("c6","PRL2 variables",1200,1200);
    
  c6->Divide(2,2);
  c6->cd(1);
  // gPad->SetLogy();


  // set-up prl2 y-dim for al categories of events

  Double_t prl2_ymin = -0.4;
  Double_t prl2_ymax = 0.4;
  
  Double_t prl2_xmin = -1.5;
  Double_t prl2_xmax = 1.5;


  
  TH2D *prl2_plot_e = new TH2D("prl2_plot_e","Electrons",200,prl2_ymin,prl2_ymax,200,prl2_xmin,prl2_xmax);
  
  prl2_plot_e->GetXaxis()->SetTitle("prl2 X track projection [m]");
  prl2_plot_e->GetXaxis()->CenterTitle();
  prl2_plot_e->GetYaxis()->SetTitle("prl2 Y track projection [m]");
  prl2_plot_e->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl2.trx:L.prl2.try>>prl2_plot_e",e_cut,"colz");

  // draw box to represent prl2 dimensions

  Double_t  prl2_x_height = 2.6;
  Double_t  prl2_y_width = 0.7;

  TBox* prl2_box = new TBox(-prl2_y_width/2,-prl2_x_height/2,+prl2_y_width/2,+prl2_x_height/2);

  prl2_box->SetFillStyle(0);
  prl2_box->SetLineColor(kRed);


  prl2_box->Draw("same");


  c6->Update();

  c6->cd(2);
  // gPad->SetLogy();

  TH2D *prl2_plot_pi = new TH2D("prl2_plot_pi","Pions",200,prl2_ymin,prl2_ymax,200,prl2_xmin,prl2_xmax);
  
  prl2_plot_pi->GetXaxis()->SetTitle("prl2 X track projection [m]");
  prl2_plot_pi->GetXaxis()->CenterTitle();
  prl2_plot_pi->GetYaxis()->SetTitle("prl2 Y track projection [m]");
  prl2_plot_pi->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl2.trx:L.prl2.try>>prl2_plot_pi",pi_cut,"colz");

  prl2_box->Draw("same");


  c6->cd(3);
  // gPad->SetLogy();

  TH2D *prl2_plot_zer = new TH2D("prl2_plot_zer","Zero hits",200,prl2_ymin,prl2_ymax,200,prl2_xmin,prl2_xmax);
  
  prl2_plot_zer->GetXaxis()->SetTitle("prl2 X track projection [m]");
  prl2_plot_zer->GetXaxis()->CenterTitle();
  prl2_plot_zer->GetYaxis()->SetTitle("prl2 Y track projection [m]");
  prl2_plot_zer->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl2.trx:L.prl2.try>>prl2_plot_zer",zer_cut,"colz");

  prl2_box->Draw("same");




  c6->cd(4);
  // gPad->SetLogy();

  TH2D *prl2_plot_bg = new TH2D("prl2_plot_bg","Bg hits",200,prl2_ymin,prl2_ymax,200,prl2_xmin,prl2_xmax);
  
  prl2_plot_bg->GetXaxis()->SetTitle("prl2 X track projection [m]");
  prl2_plot_bg->GetXaxis()->CenterTitle();
  prl2_plot_bg->GetYaxis()->SetTitle("prl2 Y track projection [m]");
  prl2_plot_bg->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl2.trx:L.prl2.try>>prl2_plot_bg",bg_cut,"colz");

  prl2_box->Draw("same");




  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~7
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~7

  // Plot track projection vs Cal1 (PRL1) Energy x distribution 

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~7

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~7
  
  TCanvas* c7 = new TCanvas("c7","PRL1 E_x vs tr_x variables",1200,1200);
    
  c7->Divide(2,2);
  c7->cd(1);
  // gPad->SetLogy();


  
  TH2D *prl1_EvTx_plot_e = new TH2D("prl1_EvTx_plot_e","Electrons",200,prl1_xmin,prl1_xmax,200,prl1_xmin,prl1_xmax);
  
  prl1_EvTx_plot_e->GetXaxis()->SetTitle("VDC X projection [m]");
  prl1_EvTx_plot_e->GetXaxis()->CenterTitle();
  prl1_EvTx_plot_e->GetYaxis()->SetTitle("PRL1 Energy-weighted x [m]");
  prl1_EvTx_plot_e->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl1.x:L.prl1.trx>>prl1_EvTx_plot_e",e_cut,"colz");



  c7->Update();

  c7->cd(2);
  
  TH2D *prl1_EvTx_plot_pi = new TH2D("prl1_EvTx_plot_pi","Pions",200,prl1_xmin,prl1_xmax,200,prl1_xmin,prl1_xmax);
  
  prl1_EvTx_plot_pi->GetXaxis()->SetTitle("VDC X projection [m]");
  prl1_EvTx_plot_pi->GetXaxis()->CenterTitle();
  prl1_EvTx_plot_pi->GetYaxis()->SetTitle("PRL1 Energy-weighted x [m]");
  prl1_EvTx_plot_pi->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl1.x:L.prl1.trx>>prl1_EvTx_plot_pi",pi_cut,"colz");

  

  c7->cd(3);


  TH2D *prl1_EvTx_plot_zer = new TH2D("prl1_EvTx_plot_zer","Zero hits",200,prl1_xmin,prl1_xmax,200,prl1_xmin,prl1_xmax);
  
  prl1_EvTx_plot_zer->GetXaxis()->SetTitle("VDC X projection [m]");
  prl1_EvTx_plot_zer->GetXaxis()->CenterTitle();
  prl1_EvTx_plot_zer->GetYaxis()->SetTitle("PRL1 Energy-weighted x [m]");
  prl1_EvTx_plot_zer->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl1.x:L.prl1.trx>>prl1_EvTx_plot_zer",zer_cut,"colz");


  

  c7->cd(4);


  TH2D *prl1_EvTx_plot_bg = new TH2D("prl1_EvTx_plot_bg","Bg hits",200,prl1_xmin,prl1_xmax,200,prl1_xmin,prl1_xmax);
  
  prl1_EvTx_plot_bg->GetXaxis()->SetTitle("VDC X projection [m]");
  prl1_EvTx_plot_bg->GetXaxis()->CenterTitle();
  prl1_EvTx_plot_bg->GetYaxis()->SetTitle("PRL1 Energy-weighted x [m]");
  prl1_EvTx_plot_bg->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl1.x:L.prl1.trx>>prl1_EvTx_plot_bg",bg_cut,"colz");





  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~8
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~8

  // Plot track projection vs Cal1 (PRL1) Energy y distribution 

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~8

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~8
  
  TCanvas* c8 = new TCanvas("c8","PRL1 E_y vs tr_y variables",1200,1200);
  
  c8->Divide(2,2);
  c8->cd(1);
  // gPad->SetLogy();


  
  TH2D *prl1_EvTy_plot_e = new TH2D("prl1_EvTy_plot_e","Electrons",200,prl1_ymin,prl1_ymax,200,prl1_ymin,prl1_ymax);
  
  prl1_EvTy_plot_e->GetXaxis()->SetTitle("VDC Y projection [m]");
  prl1_EvTy_plot_e->GetXaxis()->CenterTitle();
  prl1_EvTy_plot_e->GetYaxis()->SetTitle("PRL1 Energy-weighted y [m]");
  prl1_EvTy_plot_e->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl1.y:L.prl1.try>>prl1_EvTy_plot_e",e_cut,"colz");



  c8->Update();

  c8->cd(2);
  
  TH2D *prl1_EvTy_plot_pi = new TH2D("prl1_EvTy_plot_pi","Pions",200,prl1_ymin,prl1_ymax,200,prl1_ymin,prl1_ymax);
  
  prl1_EvTy_plot_pi->GetXaxis()->SetTitle("VDC Y projection [m]");
  prl1_EvTy_plot_pi->GetXaxis()->CenterTitle();
  prl1_EvTy_plot_pi->GetYaxis()->SetTitle("PRL1 Energy-weighted y [m]");
  prl1_EvTy_plot_pi->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl1.y:L.prl1.try>>prl1_EvTy_plot_pi",pi_cut,"colz");

  

  c8->cd(3);


  TH2D *prl1_EvTy_plot_zer = new TH2D("prl1_EvTy_plot_zer","Zero hits",200,prl1_ymin,prl1_ymax,200,prl1_ymin,prl1_ymax);
  
  prl1_EvTy_plot_zer->GetXaxis()->SetTitle("VDC Y projection [m]");
  prl1_EvTy_plot_zer->GetXaxis()->CenterTitle();
  prl1_EvTy_plot_zer->GetYaxis()->SetTitle("PRL1 Energy-weighted y [m]");
  prl1_EvTy_plot_zer->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl1.y:L.prl1.try>>prl1_EvTy_plot_zer",zer_cut,"colz");


  

  c8->cd(4);


  TH2D *prl1_EvTy_plot_bg = new TH2D("prl1_EvTy_plot_bg","Bg hits",200,prl1_ymin,prl1_ymax,200,prl1_ymin,prl1_ymax);
  
  prl1_EvTy_plot_bg->GetXaxis()->SetTitle("VDC Y projection [m]");
  prl1_EvTy_plot_bg->GetXaxis()->CenterTitle();
  prl1_EvTy_plot_bg->GetYaxis()->SetTitle("PRL1 Energy-weighted y [m]");
  prl1_EvTy_plot_bg->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl1.y:L.prl1.try>>prl1_EvTy_plot_bg",bg_cut,"colz");





  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~9
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~9

  // Plot track projection vs Cal1 (PRL2) Energy x distribution 

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~9

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~9
  
  TCanvas* c9 = new TCanvas("c9","PRL2 E_x vs tr_x variables",1200,1200);
    
  c9->Divide(2,2);
  c9->cd(1);
  // gPad->SetLogy();


  
  TH2D *prl2_EvTx_plot_e = new TH2D("prl2_EvTx_plot_e","Electrons",200,prl2_xmin,prl2_xmax,200,prl2_xmin,prl2_xmax);
  
  prl2_EvTx_plot_e->GetXaxis()->SetTitle("VDC X projection [m]");
  prl2_EvTx_plot_e->GetXaxis()->CenterTitle();
  prl2_EvTx_plot_e->GetYaxis()->SetTitle("PRL2 Energy-weighted x [m]");
  prl2_EvTx_plot_e->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl2.x:L.prl2.trx>>prl2_EvTx_plot_e",e_cut,"colz");



  c9->Update();

  c9->cd(2);
  
  TH2D *prl2_EvTx_plot_pi = new TH2D("prl2_EvTx_plot_pi","Pions",200,prl2_xmin,prl2_xmax,200,prl2_xmin,prl2_xmax);
  
  prl2_EvTx_plot_pi->GetXaxis()->SetTitle("VDC X projection [m]");
  prl2_EvTx_plot_pi->GetXaxis()->CenterTitle();
  prl2_EvTx_plot_pi->GetYaxis()->SetTitle("PRL2 Energy-weighted x [m]");
  prl2_EvTx_plot_pi->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl2.x:L.prl2.trx>>prl2_EvTx_plot_pi",pi_cut,"colz");

  

  c9->cd(3);


  TH2D *prl2_EvTx_plot_zer = new TH2D("prl2_EvTx_plot_zer","Zero hits",200,prl2_xmin,prl2_xmax,200,prl2_xmin,prl2_xmax);
  
  prl2_EvTx_plot_zer->GetXaxis()->SetTitle("VDC X projection [m]");
  prl2_EvTx_plot_zer->GetXaxis()->CenterTitle();
  prl2_EvTx_plot_zer->GetYaxis()->SetTitle("PRL2 Energy-weighted x [m]");
  prl2_EvTx_plot_zer->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl2.x:L.prl2.trx>>prl2_EvTx_plot_zer",zer_cut,"colz");


  

  c9->cd(4);


  TH2D *prl2_EvTx_plot_bg = new TH2D("prl2_EvTx_plot_bg","Bg hits",200,prl2_xmin,prl2_xmax,200,prl2_xmin,prl2_xmax);
  
  prl2_EvTx_plot_bg->GetXaxis()->SetTitle("VDC X projection [m]");
  prl2_EvTx_plot_bg->GetXaxis()->CenterTitle();
  prl2_EvTx_plot_bg->GetYaxis()->SetTitle("PRL2 Energy-weighted x [m]");
  prl2_EvTx_plot_bg->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl2.x:L.prl2.trx>>prl2_EvTx_plot_bg",bg_cut,"colz");





  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~10
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~10

  // Plot track projection vs Cal1 (PRL2) Energy y distribution 

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~10

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~10
  
  TCanvas* c10 = new TCanvas("c10","PRL2 E_y vs tr_y variables",1200,1200);
    
  c10->Divide(2,2);
  c10->cd(1);
  // gPad->SetLogy();


  
  TH2D *prl2_EvTy_plot_e = new TH2D("prl2_EvTy_plot_e","Electrons",200,prl2_ymin,prl2_ymax,200,prl2_ymin,prl2_ymax);
  
  prl2_EvTy_plot_e->GetXaxis()->SetTitle("VDC Y projection [m]");
  prl2_EvTy_plot_e->GetXaxis()->CenterTitle();
  prl2_EvTy_plot_e->GetYaxis()->SetTitle("PRL2 Energy-weighted y [m]");
  prl2_EvTy_plot_e->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl2.y:L.prl2.try>>prl2_EvTy_plot_e",e_cut,"colz");



  c10->Update();

  c10->cd(2);
  
  TH2D *prl2_EvTy_plot_pi = new TH2D("prl2_EvTy_plot_pi","Pions",200,prl2_ymin,prl2_ymax,200,prl2_ymin,prl2_ymax);
  
  prl2_EvTy_plot_pi->GetXaxis()->SetTitle("VDC Y projection [m]");
  prl2_EvTy_plot_pi->GetXaxis()->CenterTitle();
  prl2_EvTy_plot_pi->GetYaxis()->SetTitle("PRL2 Energy-weighted y [m]");
  prl2_EvTy_plot_pi->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl2.y:L.prl2.try>>prl2_EvTy_plot_pi",pi_cut,"colz");

  

  c10->cd(3);


  TH2D *prl2_EvTy_plot_zer = new TH2D("prl2_EvTy_plot_zer","Zero hits",200,prl2_ymin,prl2_ymax,200,prl2_ymin,prl2_ymax);
  
  prl2_EvTy_plot_zer->GetXaxis()->SetTitle("VDC Y projection [m]");
  prl2_EvTy_plot_zer->GetXaxis()->CenterTitle();
  prl2_EvTy_plot_zer->GetYaxis()->SetTitle("PRL2 Energy-weighted y [m]");
  prl2_EvTy_plot_zer->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl2.y:L.prl2.try>>prl2_EvTy_plot_zer",zer_cut,"colz");


  

  c10->cd(4);


  TH2D *prl2_EvTy_plot_bg = new TH2D("prl2_EvTy_plot_bg","Bg hits",200,prl2_ymin,prl2_ymax,200,prl2_ymin,prl2_ymax);
  
  prl2_EvTy_plot_bg->GetXaxis()->SetTitle("VDC Y projection [m]");
  prl2_EvTy_plot_bg->GetXaxis()->CenterTitle();
  prl2_EvTy_plot_bg->GetYaxis()->SetTitle("PRL2 Energy-weighted y [m]");
  prl2_EvTy_plot_bg->GetYaxis()->CenterTitle(); 

    
  T->Draw("L.prl2.y:L.prl2.try>>prl2_EvTy_plot_bg",bg_cut,"colz");




  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~11
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~11

  // Plot 1D difference between track projection and Cal x/y

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~11

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~11


  
  
  TCanvas* c11 = new TCanvas("c11","Track - Calorimeter x/y",1200,1200);

  c11->Divide(2);

  c11->cd(1);

  TH1D *prl1_E_min_Tx_plot_e = new TH1D("prl1_E_min_Tx_plot_e","PRL1 X",200,-0.1,0.1);

  prl1_E_min_Tx_plot_e->GetXaxis()->SetTitle("PRL1 - VDC X [m]");
  prl1_E_min_Tx_plot_e->GetXaxis()->CenterTitle();


  prl1_E_min_Tx_plot_e->SetLineColor(kRed);
  prl1_E_min_Tx_plot_e->SetLineWidth(2);
  
  T->Draw("L.prl1.x-L.prl1.trx>>prl1_E_min_Tx_plot_e",e_cut);



  // draw cuts for prl1_E_x - prl1_track_x


  Double_t prl1_E_min_Tx_l = -0.065;
  Double_t prl1_E_min_Tx_h = 0.065;
  
  TLine* prl1_E_min_Tx_cut_line_l = new TLine(prl1_E_min_Tx_l,0,prl1_E_min_Tx_l,prl1_E_min_Tx_plot_e->GetMaximum());


  prl1_E_min_Tx_cut_line_l->SetLineColor(kMagenta);  
  prl1_E_min_Tx_cut_line_l->SetLineStyle(9);
  prl1_E_min_Tx_cut_line_l->SetLineWidth(2);
  prl1_E_min_Tx_cut_line_l->Draw("same");


  TLine* prl1_E_min_Tx_cut_line_h = new TLine(prl1_E_min_Tx_h,0,prl1_E_min_Tx_h,prl1_E_min_Tx_plot_e->GetMaximum());


  prl1_E_min_Tx_cut_line_h->SetLineColor(kMagenta);  
  prl1_E_min_Tx_cut_line_h->SetLineStyle(9);
  prl1_E_min_Tx_cut_line_h->SetLineWidth(2);
  prl1_E_min_Tx_cut_line_h->Draw("same");

  
  
  TH1D* prl1_E_min_Tx_plot_pi = (TH1D*)prl1_E_min_Tx_plot_e->Clone("prl1_E_min_Tx_plot_pi");
  prl1_E_min_Tx_plot_pi->SetLineColor(kBlue);

  T->Draw("L.prl1.x-L.prl1.trx>>prl1_E_min_Tx_plot_pi",pi_cut,"same");

  
  TH1D* prl1_E_min_Tx_plot_zer = (TH1D*)prl1_E_min_Tx_plot_e->Clone("prl1_E_min_Tx_plot_zer");
  prl1_E_min_Tx_plot_zer->SetLineColor(kBlack);

  T->Draw("L.prl1.x-L.prl1.trx>>prl1_E_min_Tx_plot_zer",zer_cut,"same");

  
  TH1D* prl1_E_min_Tx_plot_bg = (TH1D*)prl1_E_min_Tx_plot_e->Clone("prl1_E_min_Tx_plot_bg");
  prl1_E_min_Tx_plot_bg->SetLineColor(kGreen);


  T->Draw("L.prl1.x-L.prl1.trx>>prl1_E_min_Tx_plot_bg",bg_cut,"same");

  
  
  c11->cd(2);


  TH1D *prl2_E_min_Tx_plot_e = new TH1D("prl2_E_min_Tx_plot_e","PRL2 X",200,-0.1,0.1);

  prl2_E_min_Tx_plot_e->GetXaxis()->SetTitle("PRL2 - VDC X [m]");
  prl2_E_min_Tx_plot_e->GetXaxis()->CenterTitle();


  prl2_E_min_Tx_plot_e->SetLineColor(kRed);
  prl2_E_min_Tx_plot_e->SetLineWidth(2);
  
  T->Draw("L.prl2.x-L.prl2.trx>>prl2_E_min_Tx_plot_e",e_cut);


  // draw cuts for prl2_E_x - prl2_track_x


  Double_t prl2_E_min_Tx_l = -0.045;
  Double_t prl2_E_min_Tx_h = 0.065;
  
  TLine* prl2_E_min_Tx_cut_line_l = new TLine(prl2_E_min_Tx_l,0,prl2_E_min_Tx_l,prl2_E_min_Tx_plot_e->GetMaximum());


  prl2_E_min_Tx_cut_line_l->SetLineColor(kMagenta);  
  prl2_E_min_Tx_cut_line_l->SetLineStyle(9);
  prl2_E_min_Tx_cut_line_l->SetLineWidth(2);
  prl2_E_min_Tx_cut_line_l->Draw("same");


  TLine* prl2_E_min_Tx_cut_line_h = new TLine(prl2_E_min_Tx_h,0,prl2_E_min_Tx_h,prl2_E_min_Tx_plot_e->GetMaximum());


  prl2_E_min_Tx_cut_line_h->SetLineColor(kMagenta);  
  prl2_E_min_Tx_cut_line_h->SetLineStyle(9);
  prl2_E_min_Tx_cut_line_h->SetLineWidth(2);
  prl2_E_min_Tx_cut_line_h->Draw("same");




  
  TH1D* prl2_E_min_Tx_plot_pi = (TH1D*)prl2_E_min_Tx_plot_e->Clone("prl2_E_min_Tx_plot_pi");
  prl2_E_min_Tx_plot_pi->SetLineColor(kBlue);

  T->Draw("L.prl2.x-L.prl2.trx>>prl2_E_min_Tx_plot_pi",pi_cut,"same");

  
  TH1D* prl2_E_min_Tx_plot_zer = (TH1D*)prl2_E_min_Tx_plot_e->Clone("prl2_E_min_Tx_plot_zer");
  prl2_E_min_Tx_plot_zer->SetLineColor(kBlack);

  T->Draw("L.prl2.x-L.prl2.trx>>prl2_E_min_Tx_plot_zer",zer_cut,"same");

  
  TH1D* prl2_E_min_Tx_plot_bg = (TH1D*)prl2_E_min_Tx_plot_e->Clone("prl2_E_min_Tx_plot_bg");
  prl2_E_min_Tx_plot_bg->SetLineColor(kGreen);


  T->Draw("L.prl2.x-L.prl2.trx>>prl2_E_min_Tx_plot_bg",bg_cut,"same");

 

  // currently commented out y diff between prl and track
  // only 2 columns so this give little info
  
  // c11->cd(3);

  // TH1D *prl1_E_min_Ty_plot_e = new TH1D("prl1_E_min_Ty_plot_e","PRL1 Y",200,-0.1,0.1);

  // prl1_E_min_Ty_plot_e->GetXaxis()->SetTitle("PRL1 - VDC Y [m]");
  // prl1_E_min_Ty_plot_e->GetXaxis()->CenterTitle();


  // prl1_E_min_Ty_plot_e->SetLineColor(kRed);
  // prl1_E_min_Ty_plot_e->SetLineWidth(2);
  
  // T->Draw("L.prl1.y-L.prl1.try>>prl1_E_min_Ty_plot_e",e_cut);

  
  // TH1D* prl1_E_min_Ty_plot_pi = (TH1D*)prl1_E_min_Ty_plot_e->Clone("prl1_E_min_Ty_plot_pi");
  // prl1_E_min_Ty_plot_pi->SetLineColor(kBlue);

  // T->Draw("L.prl1.y-L.prl1.try>>prl1_E_min_Ty_plot_pi",pi_cut,"same");

  
  // TH1D* prl1_E_min_Ty_plot_zer = (TH1D*)prl1_E_min_Ty_plot_e->Clone("prl1_E_min_Ty_plot_zer");
  // prl1_E_min_Ty_plot_zer->SetLineColor(kBlack);

  // T->Draw("L.prl1.y-L.prl1.try>>prl1_E_min_Ty_plot_zer",zer_cut,"same");

  
  // TH1D* prl1_E_min_Ty_plot_bg = (TH1D*)prl1_E_min_Ty_plot_e->Clone("prl1_E_min_Ty_plot_bg");
  // prl1_E_min_Ty_plot_bg->SetLineColor(kGreen);


  // T->Draw("L.prl1.y-L.prl1.try>>prl1_E_min_Ty_plot_bg",bg_cut,"same");



  // c11->cd(4);

  // TH1D *prl2_E_min_Ty_plot_e = new TH1D("prl2_E_min_Ty_plot_e","PRL2 Y",200,-0.1,0.1);

  // prl2_E_min_Ty_plot_e->GetXaxis()->SetTitle("PRL2 - VDC Y [m]");
  // prl2_E_min_Ty_plot_e->GetXaxis()->CenterTitle();


  // prl2_E_min_Ty_plot_e->SetLineColor(kRed);
  // prl2_E_min_Ty_plot_e->SetLineWidth(2);
  
  // T->Draw("L.prl2.y-L.prl2.try>>prl2_E_min_Ty_plot_e",e_cut);

  
  // TH1D* prl2_E_min_Ty_plot_pi = (TH1D*)prl2_E_min_Ty_plot_e->Clone("prl2_E_min_Ty_plot_pi");
  // prl2_E_min_Ty_plot_pi->SetLineColor(kBlue);

  // T->Draw("L.prl2.y-L.prl2.try>>prl2_E_min_Ty_plot_pi",pi_cut,"same");

  
  // TH1D* prl2_E_min_Ty_plot_zer = (TH1D*)prl2_E_min_Ty_plot_e->Clone("prl2_E_min_Ty_plot_zer");
  // prl2_E_min_Ty_plot_zer->SetLineColor(kBlack);

  // T->Draw("L.prl2.y-L.prl2.try>>prl2_E_min_Ty_plot_zer",zer_cut,"same");

  
  // TH1D* prl2_E_min_Ty_plot_bg = (TH1D*)prl2_E_min_Ty_plot_e->Clone("prl2_E_min_Ty_plot_bg");
  // prl2_E_min_Ty_plot_bg->SetLineColor(kGreen);


  // T->Draw("L.prl2.y-L.prl2.try>>prl2_E_min_Ty_plot_bg",bg_cut,"same");



 

  
  
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
  
  TH1D* Cal_EP_plot =  new TH1D("Cal_EP_plot","PRL1 + PRL2 Energy",100,-.05,1.5);


  Cal_EP_plot->SetLineColor(kBlack);
  

  Cal_EP_plot->GetXaxis()->SetTitle("(PRL1.e + PRL2.e)/L.gold.p");
  Cal_EP_plot->GetXaxis()->CenterTitle();


  T->Draw("(L.prl2.e + L.prl1.e)/(L.gold.p*1000)>>Cal_EP_plot",gen_cut,"");



  // electron cut on PRL sum
    

  Double_t e_prl1_2_min = 0.8;
  Double_t e_prl1_2_max = 1.01;


  TLine* e_prl_com_l = new TLine(e_prl1_2_min,0,e_prl1_2_min,Cal_EP_plot->GetMaximum());

  e_prl_com_l->SetLineWidth(2);
  e_prl_com_l->SetLineColor(kRed);
  
  e_prl_com_l->Draw("same");


  TLine* e_prl_com_h = new TLine(e_prl1_2_max,0,e_prl1_2_max,Cal_EP_plot->GetMaximum());

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

  TH1D* Cal2_EP_plot =  new TH1D("Cal2_EP_plot","PRL2 Energy",100,-0.05,1.2);
  
  Cal2_EP_plot->GetXaxis()->SetTitle("(PRL2.e)/L.gold.p");
  Cal2_EP_plot->GetXaxis()->CenterTitle();
  
  Cal2_EP_plot->SetLineColor(kBlack);


  T->Draw("(L.prl2.e)/(L.gold.p*1000)>>Cal2_EP_plot",gen_cut,"");
  


  // electron cuts based on PRL2


 Double_t e_prl2_min = 0.14;
 Double_t e_prl2_max = 0.33;

 
 TLine* e_prl2_s_l= new TLine(e_prl2_min,0,e_prl2_min,Cal2_EP_plot->GetMaximum());
 TLine* e_prl2_s_h= new TLine(e_prl2_max,0,e_prl2_max,Cal2_EP_plot->GetMaximum());
  
 
 e_prl2_s_l->SetLineWidth(2);
 e_prl2_s_l->SetLineColor(kRed);
 
 e_prl2_s_l->Draw("same");
 
 e_prl2_s_h->SetLineWidth(2);
 e_prl2_s_h->SetLineColor(kRed);
 
 e_prl2_s_h->Draw("same");
  


 // pion cuts based on PRL2

 Double_t pi_prl2_min = 0.09;
 Double_t pi_prl2_max = 0.15;

 TLine* pi_prl2_s_l= new TLine(pi_prl2_min,0,pi_prl2_min,Cal2_EP_plot->GetMaximum());
 TLine* pi_prl2_s_h= new TLine(pi_prl2_max,0,pi_prl2_max,Cal2_EP_plot->GetMaximum());


 pi_prl2_s_l->SetLineWidth(2);
 pi_prl2_s_l->SetLineColor(kBlue);
 
 pi_prl2_s_l->Draw("same");
 
 pi_prl2_s_h->SetLineWidth(2);
 pi_prl2_s_h->SetLineColor(kBlue);
 
 pi_prl2_s_h->Draw("same");




  // plot first calorimeter energy (energy normalised to track momentum) with cuts for electrons and pions

 
 c50->cd(3);

 
 gPad->SetLogy();
 
 TH1D* Cal1_EP_plot =  new TH1D("Cal1_EP_plot","PRL1 Energy",100,-0.05,1.2);
 
 
 Cal1_EP_plot->GetXaxis()->SetTitle("(PRL1.e)/L.gold.p");
 Cal1_EP_plot->GetXaxis()->CenterTitle();
  
 Cal1_EP_plot->SetLineColor(kBlack);
  
 T->Draw("(L.prl1.e)/(L.gold.p*1000)>>Cal1_EP_plot",gen_cut,"");
    
 


  
 // electron cuts based on PRL1
  

  Double_t e_prl1_min_1 = 0.55;
  Double_t e_prl1_max_1 = 0.8;
  
 
  TLine* e_prl1_s_l= new TLine(e_prl1_min_1,0,e_prl1_min_1,Cal2_EP_plot->GetMaximum());
  TLine* e_prl1_s_h= new TLine(e_prl1_max_1,0,e_prl1_max_1,Cal2_EP_plot->GetMaximum());
  
  
  e_prl1_s_l->SetLineWidth(2);
  e_prl1_s_l->SetLineColor(kRed);
  
  e_prl1_s_l->Draw("same");
  
  e_prl1_s_h->SetLineWidth(2);
  e_prl1_s_h->SetLineColor(kRed);
  
  e_prl1_s_h->Draw("same");
 

  
  // pion cuts for PRL1


  Double_t pi_prl1_min = 0.075;
  Double_t pi_prl1_max = 0.16;
  
  TLine* pi_prl1_s_l= new TLine(pi_prl1_min,0,pi_prl1_min,Cal2_EP_plot->GetMaximum());
  TLine* pi_prl1_s_h= new TLine(pi_prl1_max,0,pi_prl1_max,Cal2_EP_plot->GetMaximum());
  
  
  pi_prl1_s_l->SetLineWidth(2);
  pi_prl1_s_l->SetLineColor(kBlue);
  
  pi_prl1_s_l->Draw("same");
 
  pi_prl1_s_h->SetLineWidth(2);
  pi_prl1_s_h->SetLineColor(kBlue);
  
  pi_prl1_s_h->Draw("same");





  // plot second vs first calorimeter energy (normalised to track momentum) with cuts for electrons and pions



  c50->cd(4);
  

  TH2D *Cal_2Vs1 = new TH2D("Cal_2Vs1","PRL2 vs PRL1 Energy",1000,0.0,1.2,1000,0.0,1.2);

  Cal_2Vs1->GetXaxis()->SetTitle("PRL1.e/L.gold.p");
  Cal_2Vs1->GetXaxis()->CenterTitle();
  Cal_2Vs1->GetYaxis()->SetTitle("PRL2.e/L.gold.p");
  Cal_2Vs1->GetYaxis()->CenterTitle();
  
  
  T->Draw("(L.prl2.e/(L.gold.p*1000)):(L.prl1.e/(L.gold.p*1000))>>Cal_2Vs1",gen_cut,"col");
  

  // Draw horizontal line for electron based on prl2 energy depostition


  TLine* e_prl2_l= new TLine(0,e_prl2_min,1.2,e_prl2_min);
  TLine* e_prl2_h= new TLine(0,e_prl2_max,1.2,e_prl2_max);
  
  e_prl2_l->SetLineWidth(2);
  e_prl2_l->SetLineColor(kRed);
  
  e_prl2_l->Draw("same");

  e_prl2_h->SetLineWidth(2);
  e_prl2_h->SetLineColor(kRed);

  e_prl2_h->Draw("same");
  
  

  // draw vertical cuts

  TLine* e_prl1_l= new TLine(e_prl1_min_1,0,e_prl1_min_1,1.2);
  TLine* e_prl1_h= new TLine(e_prl1_max_1,0,e_prl1_max_1,1.2);


  e_prl1_l->SetLineWidth(2);
  e_prl1_l->SetLineColor(kRed);
  
  e_prl1_l->Draw("same");

  e_prl1_h->SetLineWidth(2);
  e_prl1_h->SetLineColor(kRed);

  e_prl1_h->Draw("same");



  // draw slopes in Calorimeter 2 vs Calorimeter 1 (based on cuts on Cal2_energy + Cal1_energy >(<) x

  Double_t e_prl2_1_min   =  e_prl1_2_min;
  Double_t e_prl2_1_max   =  e_prl1_2_max;

  Double_t e_prl1_min = e_prl1_2_min;
  Double_t e_prl1_max = e_prl1_2_max;


  
  TLine* e_prl2_1_l = new TLine(0,e_prl2_1_min,e_prl1_min,0);
  TLine* e_prl2_1_h = new TLine(0,e_prl2_1_max,e_prl1_max,0);
  
  
  e_prl2_1_l->SetLineWidth(2);
  e_prl2_1_l->SetLineColor(kRed);
  
  e_prl2_1_l->Draw("same");

  e_prl2_1_h->SetLineWidth(2);
  e_prl2_1_h->SetLineColor(kRed);

  e_prl2_1_h->Draw("same");



  // draw pion cuts

  
  TLine* pi_prl2_h = new TLine(0, pi_prl2_max,1.2,pi_prl2_max);

  pi_prl2_h->SetLineWidth(2);
  pi_prl2_h->SetLineColor(kBlue);
  pi_prl2_h->Draw("same");


  TLine* pi_prl2_l = new TLine(0, pi_prl2_min,1.2,pi_prl2_min);

  pi_prl2_l->SetLineWidth(2);
  pi_prl2_l->SetLineColor(kBlue);
  pi_prl2_l->Draw("same");


  
  TLine* pi_prl1_h = new TLine( pi_prl1_max,0,pi_prl1_max,1.2);

  pi_prl1_h->SetLineWidth(2);
  pi_prl1_h->SetLineColor(kBlue);
  pi_prl1_h->Draw("same");
  
  
  TLine* pi_prl1_l = new TLine(pi_prl1_min,0,pi_prl1_min,1.2);

  pi_prl1_l->SetLineWidth(2);
  pi_prl1_l->SetLineColor(kBlue);
  pi_prl1_l->Draw("same");




  

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

  T->Draw("L.cer.asum_c>>cer_scan",gen_cut,"");



  // plot electrons general cut as depicted in the 2nd canvas and Calorimeter based cuts in this canvas


  TCut e_c_scan_cut = gen_cut + Form("(L.prl1.e/(L.gold.p*1000)) > %f && (L.prl1.e/(L.gold.p*1000)) < %f && (L.prl2.e/(L.gold.p*1000)) > %f && (L.prl2.e/(L.gold.p*1000)) < %f && ((L.prl2.e + L.prl1.e)/(L.gold.p*1000)) > %f && ((L.prl2.e + L.prl1.e)/(L.gold.p*1000)) < %f",e_prl1_min_1,e_prl1_max_1,e_prl2_min,e_prl2_max, e_prl1_2_min, e_prl1_2_max);

  cout << "e_c_scan_cut = " << endl;
  cout << e_c_scan_cut << endl << endl;
  
 

   TH1D* cer_scan_e = (TH1D*)cer_scan->Clone("cer_scan_e");


   cer_scan_e->Reset();
   

   cer_scan_e->SetLineColor(kRed);
   
   T->Draw("L.cer.asum_c>>cer_scan_e",e_c_scan_cut,"same");


   

  // plot pions with general cut as displayed in canvas 2 and pions cuts from the calorimeter as shown in this canvas
  

  TCut pi_c_scan_cut = gen_cut + Form("(L.prl1.e/(L.gold.p*1000)) > %f && (L.prl1.e/(L.gold.p*1000)) < %f && (L.prl2.e/(L.gold.p*1000)) > %f && (L.prl2.e/(L.gold.p*1000)) < %f",pi_prl1_min,pi_prl1_max,pi_prl2_min,pi_prl2_max);


  cout << "pi_c_scan_cut = " << endl;
  cout << pi_c_scan_cut << endl << endl;
  

  TH1D* cer_scan_pi = (TH1D*)cer_scan->Clone("cer_scan_pi");

  cer_scan_pi->Reset();

  cer_scan_pi->SetLineColor(kBlue);
  
  T->Draw("L.cer.asum_c>>cer_scan_pi",pi_c_scan_cut,"same");



   
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

  T->Draw("L.cer.asum_c>>cer_scan",gen_cut,"");

  
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
  
  



  // Plot of PRL1 with electron sample and pion sample
  // later add determined cut


  c60->cd(2);

  gPad->SetLogy();
  

  Cal1_EP_plot->Draw();

  TH1D* Cal1_EP_plot_e = (TH1D*)Cal1_EP_plot->Clone("Cal1_EP_plot_e");



  Cal1_EP_plot_e->SetLineColor(kRed);
  
  Cal1_EP_plot_e->SetTitle("PRL1 Energy");
  


  // set cherenkov-determined electron cut and plot electrons 

  TCut e_cal_scan_cut = gen_cut + Form("L.cer.asum_c > %f && L.cer.asum_c < %f", e_cer_l, e_cer_h);

  T->Draw("(L.prl1.e)/(L.gold.p*1000)>>Cal1_EP_plot_e",e_cal_scan_cut,"same");


  
  

  // set cherenkov-determined pion cut and plot pions

  
  TH1D* Cal1_EP_plot_pi = (TH1D*)Cal1_EP_plot->Clone("Cal1_EP_plot_pi");


  Cal1_EP_plot_pi->SetLineColor(kBlue);

  TCut pi_cal_scan_cut = gen_cut + Form("L.cer.asum_c > %f && L.cer.asum_c < %f", pi_cer_l, pi_cer_h);

  T->Draw("(L.prl1.e)/(L.gold.p*1000)>>Cal1_EP_plot_pi",pi_cal_scan_cut,"same");

  //attempt to fit gaussians to very small peak around zero and 'true' muon peak

  // TF1* PRL1_f_zero = new TF1("PRL1_f_zero",peak,-0.05,0.05,3);

  // PRL1_f_zero->SetParameter(0,500);
  // PRL1_f_zero->SetParameter(1,0.01);
  // PRL1_f_zero->SetParLimits(1,-0.0001,0.06);
  
  // PRL1_f_zero->SetParameter(2,0.01);
  
  // Cal1_EP_plot_pi->Fit(PRL1_f_zero,"R+");

  TF1* PRL1_fit = new TF1("PRL1_fit",overall,0,0.4,8);

  // PRL1_fit->SetParameter(0,500);
  PRL1_fit->SetParameter(1,0.01);
  // PRL1_fit->SetParLimits(1,-0.0001,0.06);
   PRL1_fit->SetParameter(2,0.01);


  // PRL1_fit->SetParameter(3,1000);
  PRL1_fit->SetParameter(4,0.10);
  // PRL1_fit->SetParLimits(4,0.09,0.12);
  PRL1_fit->SetParameter(5,0.03);

  PRL1_fit->SetParameter(6,200);
  //  PRL1_fit->SetParameter(7,1000);

  Cal1_EP_plot_pi->Fit(PRL1_fit,"R0+");

  // retrieve Pion peak and background from overall fit
  
  TF1* PRL1_Pi_fit = new TF1("PRL1_Pi_fit",peak,0,0.4,3);

  PRL1_Pi_fit->SetParameter(0,PRL1_fit->GetParameter(3));
  PRL1_Pi_fit->SetParameter(1,PRL1_fit->GetParameter(4));
  PRL1_Pi_fit->SetParameter(2,PRL1_fit->GetParameter(5));

  PRL1_Pi_fit->SetLineColor(kGreen+2);
  PRL1_Pi_fit->SetLineStyle(9);
  
  PRL1_Pi_fit->Draw("same");

  
  TF1* PRL1_bg_fit = new TF1("PRL1_bg_fit",bg,0,0.4,2);

  PRL1_bg_fit->SetParameter(0,PRL1_fit->GetParameter(6));
  PRL1_bg_fit->SetParameter(1,PRL1_fit->GetParameter(7));

  PRL1_bg_fit->SetLineColor(kOrange+2);
  PRL1_bg_fit->SetLineStyle(9);
  PRL1_bg_fit->Draw("same");
  
  // Int_t nopar = PRL1_f_zero->GetNumberFreeParameters();
  // Double_t* pi_pars = PRL1_f_zero->GetParameters();

  // for(Int_t i = 0; i<nopar; i++){
  //   cout << "parameter " << i << " =  " << pi_pars[i] << endl;
  // }
  
  //  PRL1_f_zero->Draw("same");

  // Effeciency plot for PRL1 energy cut


  c60->cd(3);

  // set-up variables

  Double_t PRL1_cut = 0;
  Double_t PRL1_cut_level[No_points];

  // define variables to hold the largest electron to pion ratio and the cherenkov cut that produces this
  // set effecieny variables back to zero


  for(Int_t j = 0; j<No_points; j++){
    e_eff[j] = 0.0;
    pi_eff[j] = 0.0;
    eff_ratio[j] = 0.0;     
    PRL1_cut_level[No_points];
  }



  // create histograms to store electron effeciency, pion rejection effecieny and ratio 
  
  // upper limit of PRL1 cut
  Double_t PRL1_lim = 1;


  TH1D* PRL1_e_eff_h = new TH1D("PRL1_e_eff_h","Effeciency ratios",1000,0,PRL1_lim);
  TH1D* PRL1_pi_eff_h = new TH1D("PRL1_pi_eff_h","",1000,0,PRL1_lim);
  TH1D* PRL1_eff_ratio_h = new TH1D("PRL1_eff_ratio_h","",1000,0,PRL1_lim);



  // define variables to hold the largest electron to pion ratio and the PRL1 cut that produces this

  Double_t PRL1_largest_ratio = 0;
  Double_t PRL1_lr_cut = 0;
  


  // loop through values PRL1 cut and calculate effeciencies
  for(Int_t isum = 0; isum<No_points; isum++){
    

    // define PRL cut as going from to 1
    sum_cut = (isum+1) * (PRL1_lim / No_points);
    
    cer_cut_level[isum] = sum_cut;


    // electron effeciency calculation
    e_eff[isum] =  Cal1_EP_plot_e->Integral(Cal1_EP_plot_e->FindBin(sum_cut),Cal1_EP_plot_e->FindBin(PRL1_lim))/Cal1_EP_plot_e->Integral(Cal1_EP_plot_e->FindBin(0),Cal1_EP_plot_e->FindBin(PRL1_lim));
    
    PRL1_e_eff_h->Fill(sum_cut,e_eff[isum]);


    // pion rejection effeciency calculation

    pi_eff[isum] =  1.0 - Cal1_EP_plot_pi->Integral(Cal1_EP_plot_pi->FindBin(sum_cut),Cal1_EP_plot_pi->FindBin(PRL1_lim))/Cal1_EP_plot_pi->Integral(Cal1_EP_plot_pi->FindBin(0),Cal1_EP_plot_pi->FindBin(PRL1_lim));

    PRL1_pi_eff_h->Fill(sum_cut,pi_eff[isum]);


    // electron pion ratio calculation

    
    eff_ratio[isum] = Double_t(e_eff[isum]/(1.0 -  pi_eff[isum]));


    PRL1_eff_ratio_h->Fill(sum_cut,eff_ratio[isum]);

    
    // test if new ratio is largest 
    if(eff_ratio[isum]>PRL1_largest_ratio){
      
      PRL1_largest_ratio = eff_ratio[isum];
      PRL1_lr_cut = sum_cut;
    }





  }

  

  // draw histograms

  PRL1_e_eff_h->SetLineColor(kRed);
  
  PRL1_e_eff_h->GetXaxis()->SetTitle("PRL1 cut level");
  PRL1_e_eff_h->GetXaxis()->CenterTitle();
  PRL1_e_eff_h->GetYaxis()->SetTitle("Effeciency");
  PRL1_e_eff_h->GetYaxis()->CenterTitle();
  PRL1_e_eff_h->Draw("hist");

  
  PRL1_pi_eff_h->SetLineColor(kBlue);
  PRL1_pi_eff_h->Draw("hist  same");


  // plot ratio on same canvas but with different y-scale (effeciencies are from 0-1 but ratio of electrons to pions can be much (orders of magnitude) greater)

  right_max = 1.1*PRL1_eff_ratio_h->GetMaximum();
  ratio_scale = gPad->GetUymax()/right_max;
  
  PRL1_eff_ratio_h->SetMarkerColor(kGreen+3);
  PRL1_eff_ratio_h->SetMarkerStyle(kFullSquare);
  PRL1_eff_ratio_h->Scale(ratio_scale);
  PRL1_eff_ratio_h->SetLineColor(kGreen+3);
  PRL1_eff_ratio_h->Draw("same hist ");


  // draw seperate axis for ratio
  
  cout << "right_max for PRL1 = " << right_max << endl;

  TGaxis*PRL1axis = new TGaxis(PRL1_lim,gPad->GetUymin(),PRL1_lim,gPad->GetUymax(),0,right_max,510,"+L");
  
  cout << "  gPad->GetUxmax() = " << gPad->GetUxmax() << endl;

  PRL1axis->SetLineColor(kGreen+3);
  PRL1axis->SetLabelColor(kGreen+3);

  PRL1axis->SetTitle("Ratio (e:pi)");
  PRL1axis->CenterTitle();
  ///  PRL1axis->Rotate();
  PRL1axis->SetTitleColor(kGreen+3);
  PRL1axis->Draw();


  
  // later add determined cut




  

  // draw PRL1 + PRL2 energy (total, electrons sampled from cherenkov and pions sampled from cherenkov)


  c60->cd(4);

  gPad->SetLogy();
  
  // draw total

  Cal_EP_plot->Draw();



  // plot electrons 

  TH1D* Cal1_2_EP_plot_e = (TH1D*)Cal_EP_plot->Clone("Cal1_2_EP_plot_e");


  Cal1_2_EP_plot_e->SetLineColor(kRed);
  
  Cal1_2_EP_plot_e->SetTitle("PRL1 + PRL2 Energy");

  T->Draw("(L.prl1.e+L.prl2.e)/(L.gold.p*1000)>>Cal1_2_EP_plot_e",e_cal_scan_cut,"same");



  //  plot pions

  
  TH1D* Cal1_2_EP_plot_pi = (TH1D*)Cal_EP_plot->Clone("Cal1_2_EP_plot_pi");


  Cal1_2_EP_plot_pi->SetLineColor(kBlue);


  T->Draw("(L.prl1.e+L.prl2.e)/(L.gold.p*1000)>>Cal1_2_EP_plot_pi",pi_cal_scan_cut,"same");
  




  // Effeciency plot for PRL1 + PRL2 energy cut


  c60->cd(5);

  // set-up variables

  Double_t PRL1_2_cut = 0;
  Double_t PRL1_2_cut_level[No_points];

  // define variables to hold the largest electron to pion ratio and the cherenkov cut that produces this
  // set effecieny variables back to zero


  for(Int_t j = 0; j<No_points; j++){
    e_eff[j] = 0.0;
    pi_eff[j] = 0.0;
    eff_ratio[j] = 0.0;     
    PRL1_2_cut_level[No_points];
  }



  
  // create histograms to store electron effeciency, pion rejection effecieny and ratio 
  
  // upper limit of PRL1_2 cut
  Double_t PRL1_2_lim = 1;


  TH1D* PRL1_2_e_eff_h = new TH1D("PRL1_2_e_eff_h","Effeciency ratios",1000,0,PRL1_2_lim);
  TH1D* PRL1_2_pi_eff_h = new TH1D("PRL1_2_pi_eff_h","",1000,0,PRL1_2_lim);
  TH1D* PRL1_2_eff_ratio_h = new TH1D("PRL1_2_eff_ratio_h","",1000,0,PRL1_2_lim);



  
  // define variables to hold the largest electron to pion ratio and the PRL1_2 cut that produces this

  Double_t PRL1_2_largest_ratio = 0;
  Double_t PRL1_2_lr_cut = 0;
  

  // loop through values PRL1 + PRL2 cut and calculate effeciencies
  for(Int_t isum = 0; isum<No_points; isum++){
    

    // define PRL cut as going from to 1
    sum_cut = (isum+1) * (PRL1_2_lim / No_points);
    
    cer_cut_level[isum] = sum_cut;


    // electron effeciency calculation
    e_eff[isum] =  Cal1_2_EP_plot_e->Integral(Cal1_2_EP_plot_e->FindBin(sum_cut),Cal1_2_EP_plot_e->FindBin(PRL1_2_lim))/Cal1_2_EP_plot_e->Integral(Cal1_2_EP_plot_e->FindBin(0),Cal1_2_EP_plot_e->FindBin(PRL1_2_lim));
    
    PRL1_2_e_eff_h->Fill(sum_cut,e_eff[isum]);


    // pion rejection effeciency calculation

    pi_eff[isum] =  1.0 - Cal1_2_EP_plot_pi->Integral(Cal1_2_EP_plot_pi->FindBin(sum_cut),Cal1_2_EP_plot_pi->FindBin(PRL1_2_lim))/Cal1_2_EP_plot_pi->Integral(Cal1_2_EP_plot_pi->FindBin(0),Cal1_2_EP_plot_pi->FindBin(PRL1_2_lim));

    PRL1_2_pi_eff_h->Fill(sum_cut,pi_eff[isum]);


    // electron pion ratio calculation

    
    eff_ratio[isum] = Double_t(e_eff[isum]/(1.0 -  pi_eff[isum]));


    PRL1_2_eff_ratio_h->Fill(sum_cut,eff_ratio[isum]);



    
    // test if new ratio is largest 
    if(eff_ratio[isum]>PRL1_2_largest_ratio){
      
      PRL1_2_largest_ratio = eff_ratio[isum];
      PRL1_2_lr_cut = sum_cut;
    }





  }
  

  // draw histograms

  PRL1_2_e_eff_h->SetLineColor(kRed);
  
  PRL1_2_e_eff_h->GetXaxis()->SetTitle("PRL1 + PRL2 cut level");
  PRL1_2_e_eff_h->GetXaxis()->CenterTitle();
  PRL1_2_e_eff_h->GetYaxis()->SetTitle("Effeciency");
  PRL1_2_e_eff_h->GetYaxis()->CenterTitle();
  PRL1_2_e_eff_h->Draw("hist");

  
  PRL1_2_pi_eff_h->SetLineColor(kBlue);
  PRL1_2_pi_eff_h->Draw("hist  same");


  // plot ratio on same canvas but with different y-scale (effeciencies are from 0-1 but ratio of electrons to pions can be much (orders of magnitude) greater)

  right_max = 1.1*PRL1_2_eff_ratio_h->GetMaximum();
  ratio_scale = gPad->GetUymax()/right_max;
  
  PRL1_2_eff_ratio_h->SetMarkerColor(kGreen+3);
  PRL1_2_eff_ratio_h->SetMarkerStyle(kFullSquare);
  PRL1_2_eff_ratio_h->Scale(ratio_scale);
  PRL1_2_eff_ratio_h->SetLineColor(kGreen+3);
  PRL1_2_eff_ratio_h->Draw("same hist ");


  // draw seperate axis for ratio
  
  cout << "right_max for PRL1_2 = " << right_max << endl;

  TGaxis*PRL1_2axis = new TGaxis(PRL1_2_lim,gPad->GetUymin(),PRL1_2_lim,gPad->GetUymax(),0,right_max,510,"+L");
  
  cout << "  gPad->GetUxmax() = " << gPad->GetUxmax() << endl;

  PRL1_2axis->SetLineColor(kGreen+3);
  PRL1_2axis->SetLabelColor(kGreen+3);

  PRL1_2axis->SetTitle("Ratio (e:pi)");
  PRL1_2axis->CenterTitle();
  ///  PRL1_2axis->Rotate();
  PRL1_2axis->SetTitleColor(kGreen+3);
  PRL1_2axis->Draw();


  leg_Cher_eff->Draw("same");

  // later add determined cut




  // Plot of PRL1 vs PRL2 with cuts determined previously in canvas


  c60->cd(6);
  
  TH2D* Cal1_Vs_Cal2_EP_plot = new TH2D("Cal1_Vs_Cal2_EP_plot","",1000,-0.05,1.2,1000,-0.05,1.2);



  Cal1_Vs_Cal2_EP_plot->GetXaxis()->SetTitle("(PRL1.e)/L.gold.p");
  Cal1_Vs_Cal2_EP_plot->GetXaxis()->CenterTitle();
  
  Cal1_Vs_Cal2_EP_plot->GetYaxis()->SetTitle("(PRL2.e)/L.gold.p");
  Cal1_Vs_Cal2_EP_plot->GetYaxis()->CenterTitle();
 
  T->Draw("(L.prl2.e)/(L.gold.p*1000):(L.prl1.e)/(L.gold.p*1000)>>Cal1_Vs_Cal2_EP_plot",gen_cut,"col");
    



  // write out level of PRL1 cut(s)
 
  Double_t PRL1_cut_l = 0.2;
 
   
  TLine* PRL1_fin_cut = new TLine(PRL1_cut_l,0,PRL1_cut_l,Cal1_EP_plot_e->GetMaximum());
  PRL1_fin_cut->SetLineColor(kMagenta);
  PRL1_fin_cut->SetLineWidth(2);
 

  c60->cd(2);
  PRL1_fin_cut->Draw("same");

  // draw legend
  leg_Cher_sum->Draw("same");



  c60->cd(3);
  
  TLine* PRL1_fin_cut_eff = new TLine(PRL1_cut_l,0,PRL1_cut_l,1);
  PRL1_fin_cut_eff->SetLineColor(kMagenta);
  PRL1_fin_cut_eff->SetLineWidth(2);
  

  PRL1_fin_cut_eff->Draw("same");
  
  leg_Cher_eff->Draw("same");





  // write out level of PRL1 + PRL2 cut

  Double_t PRL1_2_cut_l = 0.51;
  Double_t PRL1_2_cut_h = 1.25;
  

  
  TLine* PRL1_2_fin_cut = new TLine(PRL1_2_cut_l,0,PRL1_2_cut_l,Cal1_2_EP_plot_e->GetMaximum());
  PRL1_2_fin_cut->SetLineColor(kMagenta);
  PRL1_2_fin_cut->SetLineWidth(2);
  
  c60->cd(4);

  PRL1_2_fin_cut->Draw("same");

  // draw legend
  leg_Cher_sum->Draw("same");

  c60->cd(5);


  TLine* PRL1_2_fin_cut_cp = new TLine(PRL1_2_cut_l,0,PRL1_2_cut_l,1);
  PRL1_2_fin_cut_cp->SetLineColor(kMagenta);
  PRL1_2_fin_cut_cp->SetLineWidth(2);
  
  PRL1_2_fin_cut_cp->Draw("same");
  
  leg_Cher_eff->Draw("same");


  


  c60->cd(6);

  // draw diagonal cut on PRL2 vs PRL1 plot

  TLine* PRL1_2_fin_cut_diag_l = new TLine(0,PRL1_2_cut_l,PRL1_2_cut_l,0);

  
  PRL1_2_fin_cut_diag_l->SetLineColor(kMagenta);
  PRL1_2_fin_cut_diag_l->SetLineWidth(2);
  PRL1_2_fin_cut_diag_l->Draw("same");

  
  TLine* PRL1_2_fin_cut_diag_h = new TLine(0,PRL1_2_cut_h,PRL1_2_cut_h,0);

  
  PRL1_2_fin_cut_diag_h->SetLineColor(kMagenta);
  PRL1_2_fin_cut_diag_h->SetLineWidth(2);
  PRL1_2_fin_cut_diag_h->Draw("same");
  
  
  PRL1_fin_cut_eff->Draw("same");


  c60->SaveAs("plots/L_Calorimeter_scan.pdf");
     
  gSystem->Exec("convert -density 700 -trim plots/L_Calorimeter_scan.pdf plots/L_Calorimeter_scan.png");


  gSystem->Exec("pdfunite plots/L_event_selection.pdf plots/L_track_var_plots.pdf plots/L_Cherenkov_scan.pdf plots/L_Calorimeter_scan.pdf plots/L_PID_plots_combined.pdf");
  


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Finally print values of cuts to terminal


  cout << "Final cuts are: " << endl;
  cout << "gen_cut: " << gen_cut << endl << endl;
  cout << "PID cuts " << Form("(L.prl1.e/(L.gold.p*1000)) > %f  && ((L.prl2.e + L.prl1.e)/(L.gold.p*1000)) > %f && ((L.prl2.e + L.prl1.e)/(L.gold.p*1000)) < %f && L.cer.asum_c > %f",PRL1_cut_l,PRL1_2_cut_l,PRL1_2_cut_h,ch_lr_cut) << endl << endl;
  cout << "Track-detector cuts: " << Form("L.s0.trx>%f && L.s0.trx<%f && L.s0.try>%f && L.s0.try<%f && L.s2.trx>%f && L.s2.trx<%f && L.s2.try>%f && L.s2.try<%f && L.prl1.trx>%f && L.prl1.trx<%f && L.prl1.try>%f && L.prl1.try<%f && L.prl2.trx>%f && L.prl2.trx<%f && L.prl2.try>%f && L.prl2.try<%f",-0.5*s0_x_height,0.5*s0_x_height,-0.5*s0_y_width,0.5*s0_y_width,-0.5*s2_x_height,0.5*s2_x_height,-0.5*s2_y_width,0.5*s2_y_width,-0.5*prl1_x_height,0.5*prl1_x_height,-0.5*prl1_y_width,0.5*prl1_y_width,-0.5*prl2_x_height,0.5*prl2_x_height,-0.5*prl2_y_width,0.5*prl2_y_width) << endl << endl;
  cout << "Track-Calorimeter agreement cuts" <<  Form("(L.prl1.x-L.prl1.trx)>%f && (L.prl1.x-L.prl1.trx)<%f && (L.prl2.x-L.prl2.trx)>%f && (L.prl2.x-L.prl2.trx)<%f",prl1_E_min_Tx_l,prl1_E_min_Tx_h,prl2_E_min_Tx_l,prl2_E_min_Tx_h) << endl << endl;
  
  


}



  



