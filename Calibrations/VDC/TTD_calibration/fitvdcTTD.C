/**************************************************************
  fitvdcTTD.C
  Seamus Riordan
  spr4y@virginia.edu
  June 16, 2010

  This script fits the VDC analytic drift time-to-distance 
  form in the THaVDCAnalyticTTDConv class.

  It's best run compiled

  analyzer
  .L fitvdcTTD.C+
  fitvdcTTD("L",runno)

  or

  fitvdcTTD("R",runno)

  for the left and right arms
  
*************************************************************/

#include <TMinuit.h>
#include <TChain.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include "THaEvent.h"
#include "TSystem.h"
#include "TLine.h"
#include "TF1.h"
#include "TProfile.h"
#include "TLegend.h"

#include "Load_more_rootfiles.C"
#include "file_def.h"

#include "TTD_namespace.h"
#include "NoiseCut.h"

#define MAX_ENTRIES 1000000
#define MAX_HIT 1000
#define NPLANE 4

#define DEG_TO_RAD 0.017453278

// Approx center of 1/tanTheta distribution
static const Double_t TT0 = 1.4;
// Old known good calibration, for comparison
static const Double_t test_vel[] = { 49.28e3, 49.73e3, 49.13e3, 49.19e3 };
static const Double_t test_apars[] = { 2.12e-3, 0, 0, 0,
				       -4.2e-4, 1.3e-3, 1.06e-4, 0 };
using namespace std;

Int_t    this_plane;
Int_t    nent[NPLANE];


std::vector<Double_t> wtime[NPLANE];
std::vector<Double_t> tanTh[NPLANE];
std::vector<Double_t> trdist[NPLANE];


bool passtrg(Int_t a,Int_t b);



Double_t TTDform( Double_t dtime, Double_t tanTheta, Double_t fDriftVel, const Double_t *par, Double_t invTanTheta0 ){
  Double_t a1 = 0.0, a2 = 0.0;

  const Double_t* fA1tdcCor = &(par[0]);
  const Double_t* fA2tdcCor = &(par[4]);

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




Double_t TTD_eval(Double_t *x, Double_t* par){

  Double_t dtime = x[0]/(1e9);
  
  Double_t* TTD_pars = &(par[0]);
  Double_t fDriftVel = par[8];
  Double_t tanTheta = 1./par[9];
  Double_t invTanTheta0 = par[10];
  
  Double_t Distance  = TTD_func::TTDform(dtime, tanTheta, fDriftVel, TTD_pars, invTanTheta0);
  
  return Distance;
  //  cout << "returning form TTD_func " << Distance << endl;  
}

Int_t Nbin = 100;

// Function to minimize
void fcn(Int_t& /*npar*/, Double_t* /*gin*/, Double_t& f, Double_t* par, Int_t /*iflag*/) {
  Int_t i;

  Double_t chisq = 0.0;
  Double_t delta;

  for (i=0; i<nent[this_plane]; i++) {
    delta =  TTD_func::TTDform(wtime[this_plane][i] - par[9], 
		     tanTh[this_plane][i],
		     par[0], &(par[1]), TT0 )
      - trdist[this_plane][i];

    //		printf("%g %g %g -> delta %f\n", wtime[this_plane][i] - par[9],  tanTh[this_plane][i], trdist[this_plane][i],  delta );

    chisq += delta*delta;
  }

  //	printf("%f\n", chisq );

  f = chisq;

  return;
}



void fitvdcTTD( const char *arm = "L", Int_t runnumber = -1 ){


  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

 
  TChain* T = new TChain("T");

  if(!strcmp(arm,"L")){
    T = Load_more_rootfiles(runnumber);
  }
  else if (!strcmp(arm,"R")){  
    T = Load_more_rootfiles(runnumber);
  }
  else{
    cout << "arm must be L or R, " << arm << " not acceptable" << endl;
    return;
  }


  // PID cuts (different for Left and Right arms)
  Double_t Cer_cut = 0.0;
  Double_t Ps_cut = 0.0;
  Double_t Ps_Sh_cut_l = 0.0;
  Double_t Ps_Sh_cut_h = 0.0;
  

  Int_t trg = 0;
  // here is where trigger is defined
  if(!strcmp(arm,"L")){
    trg = 1;
    Cer_cut = 1500.0;
    Ps_cut = 0.3;
    Ps_Sh_cut_l = 0.625;
    Ps_Sh_cut_h = 1.1;
    cout << "left trigger" << endl;
  }
  else if (!strcmp(arm,"R")){
    cout << "right trigger" << endl;
    trg = 5;
    Cer_cut = 650.0;
    Ps_cut = 0.2;
    Ps_Sh_cut_l = 0.51;
    Ps_Sh_cut_h = 1.11;
    // trg = 6;
  }
  else{
    cout << "arm must be L or R, " << arm << " not acceptable" << endl;
    return;

  }

  

  const char plane[NPLANE][8] = {"u1", "u2", "v1", "v2"};
  Double_t ang[NPLANE] = {-45.0, -45.0, 45.0, 45.0};

  Double_t nhit[NPLANE], ntr;
  Double_t hittime[NPLANE][MAX_HIT], hittrknum[NPLANE][MAX_HIT],
    hittrdist[NPLANE][MAX_HIT];
  Double_t d_th[MAX_HIT], d_ph[MAX_HIT];
  Double_t cer_sum, ps_e, sh_e;
  Double_t tr_p[100];
  THaEvent* evt = 0;

  Int_t i, j, hit;

  Double_t evttype;

  // Set up branches

  T->SetBranchStatus("Event_Branch*", kTRUE);
  T->SetBranchAddress("Event_Branch", &evt);

  T->SetBranchStatus("DR.evtypebits", kTRUE);
  T->SetBranchAddress("DR.evtypebits", &evttype);

  T->SetBranchStatus(Form("%s.tr.n", arm), kTRUE);
  T->SetBranchAddress(Form("%s.tr.n", arm), &ntr);

  T->SetBranchStatus(Form("%s.tr.d_th", arm), kTRUE);
  T->SetBranchAddress(Form("%s.tr.d_th", arm), &d_th);
  T->SetBranchStatus(Form("%s.tr.d_ph", arm), kTRUE);
  T->SetBranchAddress(Form("%s.tr.d_ph", arm), &d_ph);

  T->SetBranchStatus(Form("%s.tr.p",arm),kTRUE);
  T->SetBranchAddress(Form("%s.tr.p",arm),tr_p);
  
  T->SetBranchStatus(Form("%s.cer.asum_c",arm),kTRUE);
  T->SetBranchAddress(Form("%s.cer.asum_c",arm),&cer_sum);
  
  if(!strcmp(arm,"L")){
    T->SetBranchStatus("L.prl1.e",kTRUE);
    T->SetBranchAddress("L.prl1.e",&ps_e);
    T->SetBranchStatus("L.prl2.e",kTRUE);
    T->SetBranchAddress("L.prl2.e",&sh_e);
  }
  else if(!strcmp(arm,"R")){
    T->SetBranchStatus("R.ps.e",kTRUE);
    T->SetBranchAddress("R.ps.e",&ps_e);
    T->SetBranchStatus("R.ps.e",kTRUE);
    T->SetBranchAddress("R.sh.e",&sh_e);
  }
  

  for( i = 0; i < NPLANE; i++ ){
    T->SetBranchStatus(Form("%s.vdc.%s.nhit", arm, plane[i]), kTRUE);
    T->SetBranchAddress(Form("%s.vdc.%s.nhit", arm, plane[i]), &nhit[i]);

    T->SetBranchStatus(Form("%s.vdc.%s.time", arm, plane[i]), kTRUE);
    T->SetBranchAddress(Form("%s.vdc.%s.time", arm, plane[i]), hittime[i]);

    T->SetBranchStatus(Form("%s.vdc.%s.trknum", arm, plane[i]), kTRUE);
    T->SetBranchAddress(Form("%s.vdc.%s.trknum", arm, plane[i]), hittrknum[i]);

//     T->SetBranchStatus(Form("%s.vdc.%s.ltrdist", arm, plane[i]), kTRUE);
//     T->SetBranchAddress(Form("%s.vdc.%s.ltrdist", arm, plane[i]), hittrdist[i]);
    T->SetBranchStatus(Form("%s.vdc.%s.trdist", arm, plane[i]), kTRUE);
    T->SetBranchAddress(Form("%s.vdc.%s.trdist", arm, plane[i]), hittrdist[i]);

    nent[i] = 0;
  }

  Double_t this_slope;

  // Load up drift times and tangents
  for( i = 0; i < T->GetEntries(); i++ ){
    //  for( i = 0; i < 10; i++ ){
    T->GetEntry(i);
    if( (i%5000)==0 ) { cout << "Entry " << i << endl; }
    //    if( ntr == 1 && TTD_func::passtrg(Int_t(evttype), trg) && cer_sum > Cer_cut && ps_e/(1e3*tr_p[0]) > Ps_cut && (ps_e+sh_e)/(1e3*tr_p[0]) > Ps_Sh_cut_l &&  (ps_e+sh_e)/(1e3*tr_p[0]) < Ps_Sh_cut_h ){
    if( ntr == 1 && TTD_func::passtrg(Int_t(evttype), trg)  && cer_sum > Cer_cut){
      for( j = 0; j < NPLANE; j++ ){
	this_slope = d_th[0]*cos(ang[j]*DEG_TO_RAD) 
	  + d_ph[0]*sin(ang[j]*DEG_TO_RAD);
	for( hit = 0; hit < nhit[j] && nent[j] < MAX_ENTRIES; hit++ ){
	  if( 0 < hittime[j][hit] 
	      &&  hittime[j][hit]  < 260.0e-9
	      && hittrdist[j][hit]< 0.015
	      && hittrknum[j][hit] == 1 ){
	    wtime[j].push_back(hittime[j][hit]);
	    tanTh[j].push_back(this_slope);
	    trdist[j].push_back(hittrdist[j][hit]);
	    
	    nent[j]++;
	  }
	}
      }
    }
  }

  printf("Start fitting\n");
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  // Drift velocity + 6 constants + offset
  TMinuit *gMinuit;

  Double_t velstart = 50.0e3;
  Double_t velstep  =  1.0e1;

  Double_t conststep = 5e-7;

  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;

  TH2F *hfit[NPLANE];
  TLegend *leg_fit[NPLANE];
  
  TH2F *htime_dist[NPLANE];
  TLegend *leg_tdist[NPLANE];
  
  // plot drift time spectra for all planes
  TH1F *htime[NPLANE];

  // plot 'real' drift distance spectra for all planes
  TH1F *hdist[NPLANE];

  // plot 'real' drift distance versus time spectra for all planes
  TH2D *htime_dist_real[NPLANE];
  TH1D *hprof_time_dist_real[NPLANE];


  // track distance corrected for time versus time
  TH2D *htime_dist_Corr[NPLANE];
  // track distance corrected for time
  TH1F *hdist_Corr[NPLANE];
  // track distance corrected for time versus angle
  TH2F *hdist_Corr_ang[NPLANE];

  // slice distributions
  TH1D *h_means[NPLANE];
  TH1D *h_sigma[NPLANE];
  TH1D *h_const[NPLANE];

  // plot sigma limits above and below means of distribution
  TH1D *h_fit_up[NPLANE];
  TH1D *h_fit_low[NPLANE];
  

  
  // plot spectrum of angles/ slopes for all events
  TH1F *hslope = new TH1F("hslope_%s", "slope distribution", 100,0.5,1.0);
  hslope->GetXaxis()->SetTitle("Slope");
  hslope->GetXaxis()->CenterTitle();
  TH1F *hslope_inv = new TH1F("hslope_inv_%s", "1/slope distribution", 100,1,1.9);
  hslope_inv->GetXaxis()->SetTitle("1/Slope");
  hslope_inv->GetXaxis()->CenterTitle();
  
  TCanvas *c1[NPLANE];


  TCanvas *c2[NPLANE];

  // canvas for fitting real distance vs time distrib 
  TCanvas *c_slices[NPLANE];

  std::vector<Double_t> means[NPLANE];
  std::vector<Double_t> sigmas[NPLANE];
  
  TLine *limit[NPLANE];

  TLine *limit_time[NPLANE];

  TF1 *TTD_line[NPLANE];

  cout << "Declare lin fucntions " << endl;
  
  // function that display perfect linear correlation between 'real' distance and distance caluclated from time (y = x)
  //  TF1 *Lin_cor = new TF1("Lin_cor","pol1(0)", 0.0, 0.015);
  TF1 *Lin_cor = new TF1("Lin_cor","pol1", 0.0, 0.015);
  Lin_cor->SetParameter(0,0.0);
  Lin_cor->SetParameter(1,1.0);
  Lin_cor->SetLineStyle(4); // dashed line


  Double_t XSigma_up = 3.0; // how many sigma to cut away from mean of 'real' distance distrib

  Double_t XSigma_down = 2.0; // how many sigma to cut away from mean of 'real' distance distrib
    
    
  Double_t a1_a2[NPLANE] = {0.0};
  
  Double_t finalvel[NPLANE], apars[NPLANE][8], finaloff[NPLANE], dummy;


  NoiseCut* Real_cut[NPLANE];
  
  
  for( i = 0; i < NPLANE; i++ ){

    hfit[i] = new TH2F(Form("hfit_%s", plane[i]), Form("%s TTD", plane[i]), 200, 0.0, 0.015, 200, -0.002, 0.015 );
    hfit[i]->GetXaxis()->SetTitle("Track Dist (m)");
    hfit[i]->GetXaxis()->CenterTitle();
    hfit[i]->GetYaxis()->SetTitle("Dist from Time (m)");
    hfit[i]->GetYaxis()->CenterTitle();

    htime_dist[i] = new TH2F(Form("htime_dist_%s", plane[i]), Form("%s TTD", plane[i]), 290, -10,280, 200, 0.0, 0.015);
    htime_dist[i]->GetXaxis()->SetTitle("Drift Time(ns)");
    htime_dist[i]->GetXaxis()->CenterTitle();
    htime_dist[i]->GetYaxis()->SetTitle("Dist from Time (m)");
    htime_dist[i]->GetYaxis()->CenterTitle();


    htime[i] = new TH1F(Form("htime_%s", plane[i]), Form("%s TTD", plane[i]), 290, -10,280);
    htime[i]->GetXaxis()->SetTitle("Drift Time(ns)");
    htime[i]->GetXaxis()->CenterTitle();


    hdist[i] = new TH1F(Form("hdist_%s", plane[i]), Form("%s TTD", plane[i]), 200, 0.0, 0.015);
    hdist[i]->GetXaxis()->SetTitle("Track Dist (m)");
    hdist[i]->GetXaxis()->CenterTitle();


    htime_dist_real[i] = new TH2D(Form("htime_dist_real_%s", plane[i]), Form("%s TTD", plane[i]), 290, -10,280, 200, 0.0, 0.015);
    htime_dist_real[i]->GetXaxis()->SetTitle("Drift Time(ns)");
    htime_dist_real[i]->GetXaxis()->CenterTitle();
    htime_dist_real[i]->GetYaxis()->SetTitle("Track Dist (m)");
    htime_dist_real[i]->GetYaxis()->CenterTitle();


    htime_dist_Corr[i] = new TH2D(Form("htime_dist_Corr_%s", plane[i]), Form("%s TTD", plane[i]), 290, -10,280, 200, -0.005, 0.005);
    htime_dist_Corr[i]->GetXaxis()->SetTitle("Drift Time(ns)");
    htime_dist_Corr[i]->GetXaxis()->CenterTitle();
    htime_dist_Corr[i]->GetYaxis()->SetTitle("Track Dist (m) (Corrected for time)");
    htime_dist_Corr[i]->GetYaxis()->CenterTitle();

    hdist_Corr[i] = new TH1F(Form("hdist_Corr_%s", plane[i]), Form("%s TTD", plane[i]), 200, -0.005, 0.005);
    hdist_Corr[i]->GetXaxis()->SetTitle("Track Dist (m) (corrected for time)");
    hdist_Corr[i]->GetXaxis()->CenterTitle();


    hdist_Corr_ang[i] = new TH2F(Form("hdist_Corr_ang_%s", plane[i]), Form("%s TTD", plane[i]), 100, 1.0, 1.9, 200, -0.005, 0.005);
    hdist_Corr_ang[i]->GetXaxis()->SetTitle("1/Slope");
    hdist_Corr_ang[i]->GetXaxis()->CenterTitle();
    hdist_Corr_ang[i]->GetYaxis()->SetTitle("Track Dist (m) (corrected for time)");
    hdist_Corr_ang[i]->GetYaxis()->CenterTitle();
    
    
    for( j = 0; j < nent[i]; j++ ){

      htime[i]->Fill(1e9*(wtime[i][j] - finaloff[i]));
      hdist[i]->Fill(trdist[i][j]);
      htime_dist_real[i]->Fill(1e9*(wtime[i][j] - finaloff[i]),trdist[i][j]);
      hslope->Fill(tanTh[i][j]);
      hslope_inv->Fill(1./tanTh[i][j]);            
    }
    
   
    c2[i] = new TCanvas(Form("c2_%s",plane[i]),Form("Real Distributions %s",plane[i]), 640, 480);
    c2[i]->SetGridx();
    c2[i]->SetGridy();
    c2[i]->Divide(3,1);
    c2[i]->cd(1);
    htime[i]->Draw();

    c2[i]->cd(2);
    hdist[i]->Draw();

    c2[i]->cd(3);
    htime_dist_real[i]->Draw("colz");


    Real_cut[i] =  new NoiseCut(htime_dist_real[i]);
    
    c_slices[i] = new TCanvas(Form("c_slices_%s",plane[i]), Form("Slices distributions %s",plane[i]),640,480);

    c_slices[i]->Divide(2,2);

    c_slices[i]->cd(1);
    htime_dist_real[i]->Draw("colz");

    c_slices[i]->cd(2);
    h_const[i] = (TH1D*)gDirectory->Get(Form("htime_dist_real_%s_0",plane[i]));
    h_const[i]->Draw();

    c_slices[i]->cd(3);
    h_means[i] = (TH1D*)gDirectory->Get(Form("htime_dist_real_%s_1",plane[i]));
    h_means[i]->Draw();

    c_slices[i]->cd(4);
    h_sigma[i] = (TH1D*)gDirectory->Get(Form("htime_dist_real_%s_2",plane[i]));
    h_sigma[i]->Draw();

    h_fit_up[i] = (TH1D*) h_means[i]->Clone();
    h_fit_up[i]->Add(h_sigma[i],XSigma_up);
    h_fit_up[i]->SetLineColor(kRed);

    h_fit_low[i] = (TH1D*) h_means[i]->Clone();
    h_fit_low[i]->Add(h_sigma[i],-XSigma_down);
    h_fit_low[i]->SetLineColor(kRed);


    c_slices[i]->cd(1);
    h_means[i]->SetLineColor(kBlack);
    h_means[i]->Draw("same");
    h_fit_up[i]->Draw("hist  L same");
    h_fit_low[i]->Draw("hist L same");
    

  }



  // perform noise cuts
  for( i = 0; i < NPLANE; i++ ){
    
    Real_cut[i]->PassNoiseCut(wtime[i],trdist[i]);
    nent[i] = wtime[i].size();
    
  }

 
  // fit line to 'real' distance versus time spectra to extract drift-velocity, v_d
  // then use as input paramter to optimisation



  TF1* f_vd_mean[NPLANE];

  Double_t f_vd_start = 75.0;
  Double_t f_vd_stop = 200.0;

  Double_t vd_fit[NPLANE] = {0};

  for( i = 0; i < NPLANE; i++ ){



    f_vd_mean[i] = new TF1(Form("f_vd_mean_%i",i),"pol1",f_vd_start,f_vd_stop);
    h_means[i]->Fit(Form("f_vd_mean_%i",i),"QR","");

    cout << "Drift (mean calculated) velocity for " << plane[i] << " = " << 1e9*(f_vd_mean[i]->GetParameter(1)) << endl;

    vd_fit[i] = 1e9*(f_vd_mean[i]->GetParameter(1));
    
  }
  
  
  for( i = 0; i < NPLANE; i++ ){

    this_plane = i;
    cout << "Minimisation for " << plane[i] << endl;

    gMinuit	= new TMinuit(10); 


    gMinuit->Command("SET PRINT -1");
    gMinuit->SetFCN(fcn);

    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

    
    // gMinuit->mnparm(0, "vel", velstart, velstep, 4e4, 5.001e4,ierflg);
    gMinuit->mnparm(0, "vel", vd_fit[i], velstep, 4e4, 6.000e4,ierflg);
    //    gMinuit->FixParameter(0);

    //		gMinuit->mnparm(1, "a1_0", 0.0, conststep, 0,0,ierflg);
    gMinuit->mnparm(1, "a1_0", 0.002, 1e-6, 0.001,0.003,ierflg);
    //    gMinuit->mnparm(2, "a1_1", 0.0, conststep, -2e-3,2e-3,ierflg);
    // gMinuit->mnparm(3, "a1_2", 0.0, conststep, -1e-3,1e-3,ierflg);
    gMinuit->mnparm(2, "a1_1", 0.0, conststep,0,0,ierflg);
    gMinuit->FixParameter(2);
    gMinuit->mnparm(3, "a1_2", 0.0, conststep,0,0,ierflg);
    gMinuit->FixParameter(3);
    //    gMinuit->mnparm(4, "a1_3", 0.0, conststep, -1e-3,1e-3,ierflg);
    //    gMinuit->mnparm(4, "a1_3", -0.000358453, conststep, -1e-3,1e-3,ierflg);
    gMinuit->mnparm(4, "a1_3", 0.0, conststep, 0,0,ierflg);
    gMinuit->FixParameter(4);
    /*
    gMinuit->mnparm(2, "a1_1", 0.0, 0.0, 0,0,ierflg);
    gMinuit->mnparm(3, "a1_2", 0.0, 0.0, 0,0,ierflg);
    gMinuit->mnparm(4, "a1_3", 0.0, 0.0, 0,0,ierflg);
    */

    //		gMinuit->mnparm(5, "a2_0", 0.0, conststep, 0,0,ierflg);
    //    gMinuit->mnparm(5, "a2_0", 0.0, 1e-6, -2e-3, 2e-3,ierflg);
    gMinuit->mnparm(5, "a2_0", -4.20e-04, 1e-6, -2e-3, 2e-3,ierflg);
    gMinuit->mnparm(6, "a2_1", 0.0, conststep, -2e-3, 2e-3,ierflg);
    gMinuit->FixParameter(6);
    gMinuit->mnparm(7, "a2_2", 0.0, conststep, -1e-3, 1e-3,ierflg);
    gMinuit->FixParameter(7);
    //    gMinuit->mnparm(8, "a2_3", 0.0, conststep, -1e-3, 1e-3,ierflg);
    gMinuit->mnparm(8, "a2_3", 0.0, conststep,0,0,ierflg);
    gMinuit->FixParameter(8);
    /*
    gMinuit->mnparm(6, "a2_1", 0.0, 0.0, 0,0,ierflg);
    gMinuit->mnparm(7, "a2_2", 0.0, 0.0, 0,0,ierflg);
    gMinuit->mnparm(8, "a2_3", 0.0, 0.0, 0,0,ierflg);
    */

    //gMinuit->mnparm(9, "offset", 0.0e-7, 0.5e-9, -25e-9, 25e-9,ierflg);
    gMinuit->mnparm(9, "offset", 0.0e-7, 0.0e-9, 0,0,ierflg);

    arglist[0] = 20000;
    arglist[1] = 1.;
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

    gMinuit->GetParameter(0, finalvel[i], dummy );
    gMinuit->GetParameter(9, finaloff[i], dummy );

    for( j = 0; j < 8; j++ ){
      gMinuit->GetParameter(j+1, apars[i][j], dummy );
    }

    
    TTD_line[i] = new TF1(Form("TTD_eval_%s",plane[i]),TTD_eval,-10,280, 11);
    TTD_line[i]->SetParameters(apars[i]);
    TTD_line[i]->SetParameter(8,finalvel[i]); 
    TTD_line[i]->SetParameter(9,TT0);
    TTD_line[i]->SetParameter(10,TT0);
    TTD_line[i]->SetLineColor(kBlack);
    TTD_line[i]->SetLineStyle(4); // dashed line
    

  }



  for( i = 0; i < NPLANE; i++ ){
        
    for( j = 0; j < nent[i]; j++ ){      
      
      Double_t time = wtime[i][j];
      Double_t dist = trdist[i][j];
      
      hfit[i]->Fill( trdist[i][j], TTD_func::TTDform( wtime[i][j] - finaloff[i], tanTh[i][j], finalvel[i], apars[i], TT0 ) );
      htime_dist[i]->Fill(1e9*(wtime[i][j] - finaloff[i]), TTD_func::TTDform( wtime[i][j] - finaloff[i], tanTh[i][j], finalvel[i], apars[i], TT0 ));
      
    }
  }

  for( i = 0; i < NPLANE; i++ ){
    c1[i] = new TCanvas(Form("c1_%s",plane[i]), Form("TTD Fit %s",plane[i]), 640, 480);
    c1[i]->SetGridx();
    c1[i]->SetGridy();
    c1[i]->Divide(2,1);
    c1[i]->cd(1);
    hfit[i]->Draw("COLZ");
    
    c1[i]->cd(2);
    htime_dist[i]->Draw("COLZ");
  }

  
  
  // retrieve parameters in terms of tan(theta), as opposed to tan(theta)-tan(theta_0) where theta_0 is the central angle

  
  double apars_ex[NPLANE][8];

  for( i = 0; i < NPLANE; i++ ){
    printf("\n\n");

    printf("Plane %s\n", plane[i]);
    printf("vel    = %g\n", finalvel[i] );
    printf("offset = %5.1f ns\n", finaloff[i]*1e9 );

    apars_ex[i][0] = apars[i][0] - apars[i][1]*TT0 + apars[i][2]*pow(TT0,2) - apars[i][3]*pow(TT0, 3);
    apars_ex[i][1] = apars[i][1] - 2.0*apars[i][2]*pow(TT0,1) + 3.0*apars[i][3]*pow(TT0, 2);
    apars_ex[i][2] = apars[i][2] - 3.0*apars[i][3]*pow(TT0, 1);
    apars_ex[i][3] = apars[i][3];

    apars_ex[i][4] = apars[i][4] - apars[i][5]*TT0 + apars[i][6]*pow(TT0,2) - apars[i][7]*pow(TT0, 3);
    apars_ex[i][5] = apars[i][5] - 2.0*apars[i][6]*pow(TT0,1) + 3.0*apars[i][7]*pow(TT0, 2);
    apars_ex[i][6] = apars[i][6] - 3.0*apars[i][7]*pow(TT0, 1);
    apars_ex[i][7] = apars[i][7];

    for( j = 0; j < 4; j++ ){
      printf("fA1tdcCor[%d] = %g;\n", j, apars_ex[i][j] );
    }
    for( j = 4; j < 8; j++ ){
      printf("fA2tdcCor[%d] = %g;\n", j-4, apars_ex[i][j] );
    }

    // for each plane
    // draw line for distance = a1 + a2 for central angle
    // this is where for non-linear effect has been taken into account (for central angle/slope)
    
    
    Double_t a1 = 0.0;
    Double_t a2 = 0.0;
    //    Double_t a1_a2 = 0.0;
    

    for (Int_t j = 3; j >= 1; j--) {
      a1 = TT0 * (a1 + apars_ex[i][j]);
      a2 = TT0 * (a2 + apars_ex[i][j+4]);
    }
    a1 += apars_ex[i][0];
    a2 += apars_ex[i][4];

    a1_a2[i] = a1 + a2;
    cout << "for " << plane[i] << ": " << endl;
    cout << "a1 + a2 = " << a1_a2[i] << endl;
    cout << "a1 = " << a1 << endl;
    cout << "a2 = " << a2 << endl;


    limit[i] = new TLine(0.000, a1_a2[i], 0.015, a1_a2[i]);
    limit[i]->SetLineColor(kRed);

    Double_t time_lim = 1e9*(a1/finalvel[i]);
    limit_time[i] = new TLine(time_lim,0.0,time_lim,0.015);
    limit_time[i]->SetLineColor(kRed);

    
    c1[i]->cd(1);
    limit[i]->Draw();
    Lin_cor->Draw("same");
    
    leg_fit[i]= new TLegend(.1,.65,.37,.9,"Key");
    leg_fit[i]->SetFillColor(0);    
    leg_fit[i]->AddEntry(Lin_cor,"Exact linear correlation","l");
    leg_fit[i]->AddEntry(limit[i],"a1+a2 for central slope","l");
    leg_fit[i]->Draw("same");

    c1[i]->cd(2);
    limit_time[i]->Draw();
    TTD_line[i]->Draw("same");
    leg_tdist[i]= new TLegend(.1,.65,.37,.9,"Key");
    leg_tdist[i]->SetFillColor(0);
    leg_tdist[i]->AddEntry(TTD_line[i],"TTD function for central slope","l");
    leg_tdist[i]->AddEntry(limit_time[i],"#a_1/v_d for central slope","l");
    leg_tdist[i]->Draw("same");

    
    c1[i]->Update();    

    
  }


  
  // draw slopes for all planes

  TCanvas *c3 = new TCanvas(Form("Slopes %s",plane[i]),Form("Slopes %s",plane[i]), 640, 480);
  c3->Divide(2,1);
  c3->cd(1);
  hslope->Draw();
  c3->cd(2);
  hslope_inv->Draw();



 
  // // Save results to DB file

  std::ofstream* outp = new std::ofstream;

  outp->open(Form("DB/analytic_TTD/db_%s_TTD.vdc.%d.dat", arm, runnumber) );


  for( i = 0; i < NPLANE; i++ ){
    *outp<<arm<<".vdc."<<plane[i]<<".driftvel ="<<endl;
    *outp<<finalvel[i]<<endl<<endl;
  }


  *outp<<endl;
  
  for( i = 0; i < NPLANE; i++ ){
    *outp<<arm<<".vdc."<<plane[i]<<".ttd.param ="<<endl;

    // A1 coeffecients
    for( j = 0; j < 4; j++ ){
      *outp<<apars_ex[i][j]<<" ";
    }

    *outp<<endl;
    
    // A2 coeffecients
    for( j = 0; j < 4; j++ ){
      *outp<<apars_ex[i][j+4]<<" ";
    }
    
    *outp<<endl<<endl;
    
  }

  outp->close();


  // Save results in terms of corrections from central angle

  std::ofstream* outp_cen = new std::ofstream;

  outp_cen->open(Form("DB/analytic_TTD/db_%s_TTD_cen.vdc.%d.dat", arm, runnumber) );


  for( i = 0; i < NPLANE; i++ ){
    *outp_cen<<arm<<".vdc."<<plane[i]<<".driftvel ="<<endl;
    *outp_cen<<finalvel[i]<<endl<<endl;
  }


  *outp_cen<<endl;
  
  for( i = 0; i < NPLANE; i++ ){
    *outp_cen<<arm<<".vdc."<<plane[i]<<".ttd.param ="<<endl;

    // A1 coeffecients
    for( j = 0; j < 4; j++ ){
      *outp_cen<<apars[i][j]<<" ";
    }

    *outp_cen<<endl;
    
    // A2 coeffecients
    for( j = 0; j < 4; j++ ){
      *outp_cen<<apars[i][j+4]<<" ";
    }
    
    *outp_cen<<endl<<endl;
    
  }

  outp_cen->close();

  
}


bool passtrg(Int_t evttype, Int_t trg){
  //  cout<<evttype<<"   "<<trg<<endl;
  return evttype&(1<<trg);
}

