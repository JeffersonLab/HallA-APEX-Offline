/*
*************************************************************
14/7/21 John Williamson
Script that tests Coincidence peak and pl corrections from S2 analyzer class


*************************************************************
*/

#include "file_def.h"
#include "Load_more_rootfiles.C"
#include "CsvParser.h"

#include <iostream>

void s2TW(Int_t runno=-1);


void Coinc_test(Int_t runno){


  gStyle->SetOptFit(0011);

  gStyle->SetOptStat(0);

  const Int_t NS2Pad = 16;
  const Int_t MaxS2Hit = 15;
  
 
  // calculate width of coincidence peak
  
  TChain* T = Load_more_rootfiles(runno);
    
  
  Double_t tmin = 1690, tmax = 1705;  // the limits for 'peak'
  
  TCut TimingCut  = Form("DR.rrawt2<%.4f&&DR.rrawt2>%.4f", tmax, tmin);



  const Double_t fTdc2T = 0.5e-9;      // seconds/channel


  // variables used for cutting
  Double_t L_tr_n,L_cer_asum_c,L_ps_e,L_sh_e;
  Double_t L_tr_p[100],L_s0_trx[100],L_s2_try[100];

  Double_t R_tr_n,R_cer_asum_c,R_ps_e,R_sh_e;
  Double_t R_tr_p[100],R_s0_trx[100],R_s2_try[100];

  // pl variables
  Double_t L_r_x[100], L_r_y[100], L_tg_dp[100], L_r_th[100], L_r_ph[100];
  Double_t R_r_x[100], R_r_y[100], R_tg_dp[100], R_r_th[100], R_r_ph[100];
  
  
  Double_t L_s2_lt[NS2Pad],L_s2_rt[NS2Pad];
  Double_t L_s2_lt_c[NS2Pad],L_s2_rt_c[NS2Pad];
  Double_t L_s2_time[NS2Pad];
  Double_t L_s2_plcorr[NS2Pad];
  Double_t L_s2_nthit;
  Double_t L_s2_t_pads[NS2Pad];
  Double_t L_s2_la_p[NS2Pad];
  Double_t L_s2_la_c[NS2Pad];
  Double_t L_s2_ra_p[NS2Pad];
  Double_t L_s2_ra_c[NS2Pad];
  Double_t L_s2_matches[NS2Pad];

  Double_t R_s2_lt[NS2Pad],R_s2_rt[NS2Pad];
  Double_t R_s2_lt_c[NS2Pad],R_s2_rt_c[NS2Pad];
  Double_t R_s2_time[NS2Pad];
  Double_t R_s2_plcorr[NS2Pad];
  Double_t R_s2_nthit;
  Double_t R_s2_t_pads[NS2Pad];
  Double_t R_s2_la_p[NS2Pad];
  Double_t R_s2_la_c[NS2Pad];
  Double_t R_s2_ra_p[NS2Pad];
  Double_t R_s2_ra_c[NS2Pad];
  Double_t R_s2_matches[NS2Pad];

  
  //trigger
  Double_t Trig_type;
  
  Double_t Trig_time_T2[MaxS2Hit];
  Double_t Trig_No_T2;
  Double_t Trig_time_T5[MaxS2Hit];
  Double_t Trig_No_T5;

  
  
  T->SetBranchStatus("*",0);

  T->SetBranchStatus("DL.evtype",1);
  
  T->SetBranchStatus("L.tr.n",1);
  T->SetBranchStatus("L.tr.p",1);
  T->SetBranchStatus("L.cer.asum_c",1);
  T->SetBranchStatus("L.prl1.e",1);
  T->SetBranchStatus("L.prl2.e",1);

  T->SetBranchStatus("L.tr.r_x",1);
  T->SetBranchStatus("L.tr.r_y",1);
  T->SetBranchStatus("L.tr.r_th",1);
  T->SetBranchStatus("L.tr.r_ph",1);
  T->SetBranchStatus("L.tr.tg_dp",1);
  
  T->SetBranchAddress("L.tr.r_x",L_r_x);
  T->SetBranchAddress("L.tr.r_y",L_r_y);
  T->SetBranchAddress("L.tr.r_th",L_r_th);
  T->SetBranchAddress("L.tr.r_ph",L_r_ph);
  T->SetBranchAddress("L.tr.tg_dp",L_tg_dp);
  
  T->SetBranchStatus("R.tr.n",1);
  T->SetBranchStatus("R.tr.p",1);
  T->SetBranchStatus("R.cer.asum_c",1);
  T->SetBranchStatus("R.ps.e",1);
  T->SetBranchStatus("R.sh.e",1);


  T->SetBranchStatus("R.tr.r_x",1);
  T->SetBranchStatus("R.tr.r_y",1);
  T->SetBranchStatus("R.tr.r_th",1);
  T->SetBranchStatus("R.tr.r_ph",1);
  T->SetBranchStatus("R.tr.tg_dp",1);
  
  T->SetBranchAddress("R.tr.r_x",R_r_x);
  T->SetBranchAddress("R.tr.r_y",R_r_y);
  T->SetBranchAddress("R.tr.r_th",R_r_th);
  T->SetBranchAddress("R.tr.r_ph",R_r_ph);
  T->SetBranchAddress("R.tr.tg_dp",R_tg_dp);
  
  
  T->SetBranchStatus("L.s2.nthit",1);
  T->SetBranchStatus("L.s2.t_pads",1);
  T->SetBranchStatus("L.s2.lt",1);
  T->SetBranchStatus("L.s2.rt",1);
  T->SetBranchStatus("L.s2.lt_c",1);
  T->SetBranchStatus("L.s2.rt_c",1);
  T->SetBranchStatus("L.s2.time",1);
  T->SetBranchStatus("L.s2.pl_corr",1);
  T->SetBranchStatus("L.s2.la_p",1);
  T->SetBranchStatus("L.s2.ra_p",1);
  T->SetBranchStatus("L.s2.la_c",1);
  T->SetBranchStatus("L.s2.ra_c",1);
  T->SetBranchStatus("L.s2.tr_Matches",1); 
  T->SetBranchStatus("R.s2.nthit",1);
  T->SetBranchStatus("R.s2.t_pads",1);
  T->SetBranchStatus("R.s2.lt",1);
  T->SetBranchStatus("R.s2.rt",1);
  T->SetBranchStatus("R.s2.lt_c",1);
  T->SetBranchStatus("R.s2.rt_c",1);
  T->SetBranchStatus("R.s2.time",1);
  T->SetBranchStatus("R.s2.pl_corr",1);
  T->SetBranchStatus("R.s2.la_p",1);
  T->SetBranchStatus("R.s2.ra_p",1);
  T->SetBranchStatus("R.s2.la_c",1);
  T->SetBranchStatus("R.s2.ra_c",1);
  T->SetBranchStatus("R.s2.tr_Matches",1); 

  T->SetBranchStatus("DR.rrawt2",1); 
  T->SetBranchStatus("Ndata.DR.rrawt2",1);
  T->SetBranchStatus("DR.rrawt5",1);
  T->SetBranchStatus("Ndata.DR.rrawt5",1);
  
  
  T->SetBranchAddress("DL.evtype",&Trig_type);
  
  T->SetBranchAddress("L.tr.n",&L_tr_n);
  T->SetBranchAddress("L.tr.p",L_tr_p);
  T->SetBranchAddress("L.cer.asum_c",&L_cer_asum_c);
  T->SetBranchAddress("L.prl1.e",&L_ps_e);
  T->SetBranchAddress("L.prl2.e",&L_sh_e);

  T->SetBranchAddress("R.tr.n",&R_tr_n);
  T->SetBranchAddress("R.tr.p",R_tr_p);
  T->SetBranchAddress("R.cer.asum_c",&R_cer_asum_c);
  T->SetBranchAddress("R.ps.e",&R_ps_e);
  T->SetBranchAddress("R.sh.e",&R_sh_e);
  
  T->SetBranchAddress("L.s2.lt",L_s2_lt);
  T->SetBranchAddress("L.s2.rt",L_s2_rt);
  T->SetBranchAddress("L.s2.lt_c",L_s2_lt_c);
  T->SetBranchAddress("L.s2.rt_c",L_s2_rt_c);
  T->SetBranchAddress("L.s2.time",L_s2_time);
  T->SetBranchAddress("L.s2.pl_corr",L_s2_plcorr);
  T->SetBranchAddress("L.s2.nthit",&L_s2_nthit);
  T->SetBranchAddress("L.s2.t_pads",L_s2_t_pads);
  T->SetBranchAddress("L.s2.la_p",L_s2_la_p);
  T->SetBranchAddress("L.s2.ra_p",L_s2_ra_p);
  T->SetBranchAddress("L.s2.la_c",L_s2_la_c);
  T->SetBranchAddress("L.s2.ra_c",L_s2_ra_c);
  T->SetBranchAddress("L.s2.tr_Matches",L_s2_matches); 
  
  T->SetBranchAddress("R.s2.lt",R_s2_lt);
  T->SetBranchAddress("R.s2.rt",R_s2_rt);
  T->SetBranchAddress("R.s2.lt_c",R_s2_lt_c);
  T->SetBranchAddress("R.s2.rt_c",R_s2_rt_c);
  T->SetBranchAddress("R.s2.time",R_s2_time);
  T->SetBranchAddress("R.s2.pl_corr",R_s2_plcorr);
  T->SetBranchAddress("R.s2.nthit",&R_s2_nthit);
  T->SetBranchAddress("R.s2.t_pads",R_s2_t_pads);
  T->SetBranchAddress("R.s2.la_p",R_s2_la_p);
  T->SetBranchAddress("R.s2.ra_p",R_s2_ra_p);
  T->SetBranchAddress("R.s2.la_c",R_s2_la_c);
  T->SetBranchAddress("R.s2.ra_c",R_s2_ra_c);
  T->SetBranchAddress("R.s2.tr_Matches",R_s2_matches); 
  
  T->SetBranchAddress("DR.rrawt2", Trig_time_T2); 
  T->SetBranchAddress("Ndata.DR.rrawt2", &Trig_No_T2);
  T->SetBranchAddress("DR.rrawt5",Trig_time_T5);
  T->SetBranchAddress("Ndata.DR.rrawt5",&Trig_No_T5);
  

  Int_t nentries = T->GetEntries();
  //  nentries = 10e4;
  //  Int_t nentries = 100;


  // pl variable limits

  Double_t L_x_low = -0.5;
  Double_t L_x_up = 0.4;
  
  Double_t L_th_low = -0.022;
  Double_t L_th_up = 0.020;

  Double_t L_ph_low = -0.03;
  Double_t L_ph_up = 0.025;

  Double_t L_dp_low = -0.04;
  Double_t L_dp_up = 0.045;

  
  Double_t R_x_low = -0.62;
  Double_t R_x_up = 0.5;
    
  Double_t R_th_low = -0.02;
  Double_t R_th_up = 0.02;

  Double_t R_ph_low = -0.015;
  Double_t R_ph_up = 0.02;

  Double_t R_dp_low = -0.04;
  Double_t R_dp_up = 0.045;


  
  
  // coincidence histogram widths
  //  Int_t coinc_bins = 140;
  Int_t coinc_bins = 400;
  Double_t coinc_start = 0;
  Double_t coinc_end = 2e-7;


  // smaller coincidence window

  Int_t scoinc_bins = 50;
  Double_t scoinc_start = 1.60e-7; //uncorrected
  Double_t scoinc_end = 1.70e-7; // uncorrected

  
  // intitialise timing coincidence histogram
  //  TH1F *h1 = new TH1F("h1"," S2-Time difference (corrected)",coinc_bins,coinc_start-coin_off,coinc_end-coin_off);
  TH1F *h1 = new TH1F("h1"," S2-Time difference (pl corrected)",coinc_bins,coinc_start,coinc_end);    

  // intitialise timing coincidence histogram
  TH1F *h2 = new TH1F("h2"," S2-Time difference (uncorrected)",coinc_bins,coinc_start,coinc_end);

    
  // itnitialise paddle versus coincidence histograms
  TH2F* h_L_lvT = new TH2F("h_L_lvT","S2: Coincidence time versus l-paddle (pl corrected)",NS2Pad,1,NS2Pad,coinc_bins,coinc_start,coinc_end);
  h_L_lvT->GetXaxis()->SetTitle("LHRS S2 l-paddle #");
  h_L_lvT->GetYaxis()->SetTitle("Coincidence time (s)");
  

  TH2F* h_R_lvT = new TH2F("h_R_lvT","S2: Coincidence time versus l-paddle (pl corrected)",NS2Pad,1,NS2Pad,coinc_bins,coinc_start,coinc_end);
  h_R_lvT->GetXaxis()->SetTitle("RHRS S2 l-paddle #");
  h_R_lvT->GetYaxis()->SetTitle("Coincidence time (s)");
      

  TH2F* h_L_lvT_un = new TH2F("h_L_lvT_un","S2: Coincidence time versus l-paddle (uncorrected)",NS2Pad,1,NS2Pad,coinc_bins,coinc_start,coinc_end);
  h_L_lvT_un->GetXaxis()->SetTitle("LHRS S2 l-paddle #");
  h_L_lvT_un->GetYaxis()->SetTitle("Coincidence time (s)");

  TH2F* h_R_lvT_un = new TH2F("h_R_lvT_un","S2: Coincidence time versus l-paddle (uncorrected)",NS2Pad,1,NS2Pad,coinc_bins,coinc_start,coinc_end);
  h_R_lvT_un->GetXaxis()->SetTitle("RHRS S2 l-paddle #");
  double t = 0.0;
  h_R_lvT_un->GetYaxis()->SetTitle("Coincidence time (s)");

  
  
  // plots of coincidence time vs 1D variables

  TH2F* h_L_th_vs = new TH2F("h_L_th_vs","Coinc time vs Theta",50,L_th_low,L_th_up,scoinc_bins,scoinc_start,scoinc_end);
  h_L_th_vs->GetXaxis()->SetTitle("#theta (rad)");
  h_L_th_vs->GetYaxis()->SetTitle("Coinc time (s)");
    
  TH2F* h_L_ph_vs = new TH2F("h_L_ph_vs","Coinc time vs Phi",50,L_ph_low,L_ph_up,scoinc_bins,scoinc_start,scoinc_end);
  h_L_ph_vs->GetXaxis()->SetTitle("#phi (rad)");
  h_L_ph_vs->GetYaxis()->SetTitle("Coinc time (s)");

  TH2F* h_L_x_vs = new TH2F("h_L_x_vs","Coinc time vs x",50,L_x_low,L_x_up,scoinc_bins,scoinc_start,scoinc_end);
  h_L_x_vs->GetXaxis()->SetTitle("x (m)");
  h_L_x_vs->GetYaxis()->SetTitle("Coinc time (s)");

  TH2F* h_L_dp_vs = new TH2F("h_L_dp_vs","Coinc time vs dp",50,L_dp_low,L_dp_up,scoinc_bins,scoinc_start,scoinc_end);
  h_L_dp_vs->GetXaxis()->SetTitle("#delta p");
  h_L_dp_vs->GetYaxis()->SetTitle("Coinc time (s)");

  TH2F* h_L_pad_vs = new TH2F("h_L_pad_vs","Coinc time vs pad",16,-0.5,15.5,scoinc_bins,scoinc_start,scoinc_end);
  h_L_pad_vs->GetXaxis()->SetTitle("pad #");
  h_L_pad_vs->GetYaxis()->SetTitle("Coinc time (s)");

  TH2F* h_R_th_vs = new TH2F("h_R_th_vs","Coinc time vs Theta",50,R_th_low,R_th_up,scoinc_bins,scoinc_start,scoinc_end);
  h_R_th_vs->GetXaxis()->SetTitle("#theta (rad)");
  h_R_th_vs->GetYaxis()->SetTitle("Coinc time (s)");
    
  TH2F* h_R_ph_vs = new TH2F("h_R_ph_vs","Coinc time vs Phi",50,R_ph_low,R_ph_up,scoinc_bins,scoinc_start,scoinc_end);
  h_R_ph_vs->GetXaxis()->SetTitle("#phi (rad)");
  h_R_ph_vs->GetYaxis()->SetTitle("Coinc time (s)");

  TH2F* h_R_x_vs = new TH2F("h_R_x_vs","Coinc time vs x",50,R_x_low,R_x_up,scoinc_bins,scoinc_start,scoinc_end);
  h_R_x_vs->GetXaxis()->SetTitle("x (m)");
  h_R_x_vs->GetYaxis()->SetTitle("Coinc time (s)");

  TH2F* h_R_dp_vs = new TH2F("h_R_dp_vs","Coinc time vs dp",50,R_dp_low,R_dp_up,scoinc_bins,scoinc_start,scoinc_end);
  h_R_dp_vs->GetXaxis()->SetTitle("#delta p");
  h_R_dp_vs->GetYaxis()->SetTitle("Coinc time (s)");

  TH2F* h_R_pad_vs = new TH2F("h_R_pad_vs","Coinc time vs pad",16,-0.5,15.5,scoinc_bins,scoinc_start,scoinc_end);
  h_R_pad_vs->GetXaxis()->SetTitle("pad #");
  h_R_pad_vs->GetYaxis()->SetTitle("Coinc time (s)");
  
  
  Double_t LTime = 0.0;
  Double_t RTime = 0.0;

  Double_t LTime_un = 0.0;
  Double_t RTime_un = 0.0;


  // record which paddle is hit for an entry
  Int_t LHRS_pad = -1;
  Int_t RHRS_pad = -1;
  

  // parameters for S2 and trigger timing cuts

  Double_t T2LowCut = 1683;
  Double_t T2HighCut = 1710;

  Double_t T5LowCut = 1625;
  Double_t T5HighCut = 1634;

  Bool_t T2Pass = false;
  Bool_t T5Pass = false;
  

  
  for(Int_t i=0;i<nentries;i++){
    T->GetEntry(i);


    //    if(L_tr_n==1 && L_cer_asum_c>1500 && (L_ps_e+L_sh_e)/(1000.*L_tr_p[0])>0.8 && L_s0_nthit==1){
    LTime = 0.0;
    RTime = 0.0;
    
    LTime_un = 0.0;
    RTime_un = 0.0;
    
    LHRS_pad = -1;
    RHRS_pad = -1;

    T2Pass = false;
    T5Pass = false;
    
    for(Int_t j = 0; j < Trig_No_T2; j++){      
      if (Trig_time_T2[j] > T2LowCut && Trig_time_T2[j] < T2HighCut){	
	T2Pass = true;
      }
    }

    
    for(Int_t j = 0; j < Trig_No_T5; j++){      
      if (Trig_time_T5[j] > T5LowCut && Trig_time_T5[j] < T5HighCut){	
	T5Pass = true;
      }
    }
    

    // if(L_s2_t_pads[0] == 4 || L_s2_t_pads[0] == 4){
    //   continue;
    // }

    
    if(L_s2_nthit ==1 && R_s2_nthit==1 && L_tr_n==1 && L_cer_asum_c>1500 && (L_ps_e+L_sh_e)/(1000.*L_tr_p[0])>0.8 && R_tr_n==1 && R_cer_asum_c>1500 && (R_ps_e+R_sh_e)/(1000.*R_tr_p[0])>0.8&& Trig_type==6){

      
      for(int j=0;j<NS2Pad;j++){

          // if(j==4 || j ==5){
          //  continue;
          // }
	
	if (L_s2_t_pads[0]==j && (L_s2_lt_c[j]+L_s2_rt_c[j]) != 0 && L_s2_matches[j] == 1){ //&& L_s2_la_c[j]>50 && L_s2_ra_c[j]>50){

	  
	  
	  LTime = L_s2_time[j];
	  LTime_un = L_s2_time[j] + L_s2_plcorr[j];

	  LHRS_pad = j+1;
	
	}

      }


      for(int j=0;j<NS2Pad;j++){
	
	if (R_s2_t_pads[0]==j && (R_s2_lt_c[j]+R_s2_rt_c[j]) != 0 && R_s2_la_c[j]>50 && R_s2_ra_c[j]>50 &&  R_s2_matches[j] == 1){
	  //	  	if (R_s2_t_pads[0]==j && (R_s2_lt_c[j]+R_s2_rt_c[j]) != 0){// && R_s2_la_c[j]>50 && R_s2_ra_c[j]>50){

	  RTime = R_s2_time[j];
	  RTime_un = R_s2_time[j] + R_s2_plcorr[j];

	  RHRS_pad = j+1;

	  if(i%10000 == 0){
	    cout << "Event " << i << endl;
	    cout << endl;
	  }

	}

      }


      if( (LTime-RTime) != 0 ){
      
	h1->Fill(LTime-RTime);
	h2->Fill(LTime_un-RTime_un);
      
	h_L_lvT->Fill(LHRS_pad,LTime-RTime);
	h_R_lvT->Fill(RHRS_pad,LTime-RTime);

	h_L_lvT_un->Fill(LHRS_pad,LTime_un-RTime_un);
	h_R_lvT_un->Fill(RHRS_pad,LTime_un-RTime_un);


	h_L_th_vs->Fill(L_r_th[0],LTime-RTime);
	h_L_ph_vs->Fill(L_r_ph[0],LTime-RTime);
	h_L_x_vs->Fill(L_r_x[0],LTime-RTime);
	h_L_dp_vs->Fill(L_tg_dp[0],LTime-RTime);
	h_L_pad_vs->Fill(LHRS_pad*1.0,LTime-RTime);


	h_R_th_vs->Fill(R_r_th[0],LTime-RTime);
	h_R_ph_vs->Fill(R_r_ph[0],LTime-RTime);
	h_R_x_vs->Fill(R_r_x[0],LTime-RTime);
	h_R_dp_vs->Fill(R_tg_dp[0],LTime-RTime);
	h_R_pad_vs->Fill(RHRS_pad*1.0,LTime-RTime);
	
      }
      
    }
    
      
  }
  


  // fit peak of corrected coincidence timing

  TCanvas* c1 = new TCanvas("c1","c1",1000,800);
  
  h1->Draw();
  cout << "h1 entries = " << h1->GetEntries() << endl;;
  
  Double_t max= h1->GetBinCenter(h1->GetMaximumBin());
  cout << "h1 max = " << max << endl;

  
  TF1* f1 = new TF1("f1","pol1",max-(6e-8),max-(1e-8));
  h1->Fit(f1,"FRQ");
  
  Double_t par[5];   
  f1->GetParameters(&par[0]);
  
  h1->Fit(f1,"FRQ");
  
  f1->GetParameters(&par[0]);
  
  TF1* f2 = new TF1("f2","gaus",max-(2.5e-9),max+(2.5e-9));
 
  h1->Fit(f2,"FRQ");
 
  
  f2->GetParameters(&par[2]);
 
 

  TF1* f3 = new TF1("f3","pol1(0)+gaus(2)",max-(6e-8),max+(1e-8));
 
  f3->SetParameters(par);

  
  h1->Fit(f3,"FRQ");

  f3->GetParameters(&par[0]);

  Double_t corr_width = par[4];
  Double_t corr_mean = par[3];

  cout << "corr_width = " << corr_width << ", corr_mean = " << corr_mean << endl;

  gPad->Update();
  
  TPaveStats *st_1 = (TPaveStats*)h1->FindObject("stats");
  st_1->SetX1NDC(0.15); 
  st_1->SetX2NDC(0.45);
  st_1->SetY1NDC(0.45); 
  st_1->SetY2NDC(0.90);
  st_1->Draw();
    
  c1->Update();
  

  // fit peak of ucorrected coincidence timing

  TCanvas* c2 = new TCanvas("c2","c2",1000,800);
  
  h2->Draw();
  
  Double_t max2= h2->GetBinCenter(h2->GetMaximumBin());


  
  TF1* f1_un = new TF1("f1_un","pol1",max2-(6e-8),max2-(1e-8));
  h2->Fit(f1_un,"FRQ");
  
  Double_t par2[5];   
  f1_un->GetParameters(&par2[0]);
  
  h2->Fit(f1_un,"FRQ");
  
  f1_un->GetParameters(&par2[0]);
  
  TF1* f2_un = new TF1("f2_un","gaus",max2-(2.5e-9),max2+(2.5e-9));
 
  h2->Fit(f2_un,"FRQ");
 
  
  f2_un->GetParameters(&par2[2]);
 
 

  TF1* f3_un = new TF1("f3_un","pol1(0)+gaus(2)",max2-(6e-8),max2+(1e-8));
 
  f3_un->SetParameters(par2);
  h2->Fit(f3_un,"FRQ");

  f3_un->GetParameters(&par2[0]);
  
  Double_t uncorr_width = par2[4];
  Double_t uncorr_mean = par2[3];


  gPad->Update();

  
  TPaveStats *st_2 = (TPaveStats*)h2->FindObject("stats");
  st_2->SetX1NDC(0.15); 
  st_2->SetX2NDC(0.45); 
  st_2->SetY1NDC(0.45); 
  st_2->SetY2NDC(0.90);
  st_2->Draw();
  
  c2->Update();

  
  

  // draw coincidence time vs paddle number plots

  TCanvas* c5 = new TCanvas("c5","c5",1000,800);

  h_L_lvT->Draw("colz");


  TCanvas* c6 = new TCanvas("c6","c6",1000,800);

  h_R_lvT->Draw("colz");

  
  TCanvas* c7 = new TCanvas("c7","c7",1000,800);

  h_L_lvT_un->Draw("colz");


  TCanvas* c8 = new TCanvas("c8","c8",1000,800);

  h_R_lvT_un->Draw("colz");


  TCanvas* c9 = new TCanvas("c9","c9",1000,800);
  c9->Divide(5,2);

  c9->cd(1);
  h_L_th_vs->Draw("colz");

  c9->cd(2);
  h_L_ph_vs->Draw("colz");

  c9->cd(3);
  h_L_x_vs->Draw("colz");

  c9->cd(4);
  h_L_dp_vs->Draw("colz");

  c9->cd(5);
  h_L_pad_vs->Draw("colz");

  c9->cd(6);
  h_L_th_vs->FitSlicesY(0,0,-1,10,"QNR");

  TH1D* h_L_th_vs_mean = (TH1D*)gDirectory->Get("h_L_th_vs_1");
  h_L_th_vs_mean->GetYaxis()->SetRangeUser(scoinc_start,scoinc_end);
  h_L_th_vs_mean->Draw();

  c9->cd(7);
  h_L_ph_vs->FitSlicesY(0,0,-1,10,"QNR");

  TH1D* h_L_ph_vs_mean = (TH1D*)gDirectory->Get("h_L_ph_vs_1");
  h_L_ph_vs_mean->GetYaxis()->SetRangeUser(scoinc_start,scoinc_end);
  h_L_ph_vs_mean->Draw();

  c9->cd(8);
  h_L_x_vs->FitSlicesY(0,0,-1,10,"QNR");

  TH1D* h_L_x_vs_mean = (TH1D*)gDirectory->Get("h_L_x_vs_1");
  h_L_x_vs_mean->GetYaxis()->SetRangeUser(scoinc_start,scoinc_end);
  h_L_x_vs_mean->Draw();

  c9->cd(9);
  h_L_dp_vs->FitSlicesY(0,0,-1,10,"QNR");

  TH1D* h_L_dp_vs_mean = (TH1D*)gDirectory->Get("h_L_dp_vs_1");
  h_L_dp_vs_mean->GetYaxis()->SetRangeUser(scoinc_start,scoinc_end);
  h_L_dp_vs_mean->Draw();

  c9->cd(10);
  h_L_pad_vs->FitSlicesY(0,0,-1,10,"QNR");

  TH1D* h_L_pad_vs_mean = (TH1D*)gDirectory->Get("h_L_pad_vs_1");
  h_L_pad_vs_mean->GetYaxis()->SetRangeUser(scoinc_start,scoinc_end);
  h_L_pad_vs_mean->Draw();


  TCanvas* c10 = new TCanvas("c10","c10",1000,800);
  c10->Divide(5,2);

  c10->cd(1);
  h_R_th_vs->Draw("colz");

  c10->cd(2);
  h_R_ph_vs->Draw("colz");

  c10->cd(3);
  h_R_x_vs->Draw("colz");

  c10->cd(4);
  h_R_dp_vs->Draw("colz");

  c10->cd(5);
  h_R_pad_vs->Draw("colz");

  c10->cd(6);
  h_R_th_vs->FitSlicesY(0,0,-1,10,"QNR");

  TH1D* h_R_th_vs_mean = (TH1D*)gDirectory->Get("h_R_th_vs_1");
  h_R_th_vs_mean->GetYaxis()->SetRangeUser(scoinc_start,scoinc_end);
  h_R_th_vs_mean->Draw();

  c10->cd(7);
  h_R_ph_vs->FitSlicesY(0,0,-1,10,"QNR");

  TH1D* h_R_ph_vs_mean = (TH1D*)gDirectory->Get("h_R_ph_vs_1");
  h_R_ph_vs_mean->GetYaxis()->SetRangeUser(scoinc_start,scoinc_end);
  h_R_ph_vs_mean->Draw();

  c10->cd(8);
  h_R_x_vs->FitSlicesY(0,0,-1,10,"QNR");

  TH1D* h_R_x_vs_mean = (TH1D*)gDirectory->Get("h_R_x_vs_1");
  h_R_x_vs_mean->GetYaxis()->SetRangeUser(scoinc_start,scoinc_end);
  h_R_x_vs_mean->Draw();

  c10->cd(9);
  h_R_dp_vs->FitSlicesY(0,0,-1,10,"QNR");

  TH1D* h_R_dp_vs_mean = (TH1D*)gDirectory->Get("h_R_dp_vs_1");
  h_R_dp_vs_mean->GetYaxis()->SetRangeUser(scoinc_start,scoinc_end);
  h_R_dp_vs_mean->Draw();

  c10->cd(10);
  h_R_pad_vs->FitSlicesY(0,0,-1,10,"QNR");

  TH1D* h_R_pad_vs_mean = (TH1D*)gDirectory->Get("h_R_pad_vs_1");
  h_R_pad_vs_mean->GetYaxis()->SetRangeUser(scoinc_start,scoinc_end);
  h_R_pad_vs_mean->Draw();


  // print canvases


  c1->Print(Form("plots/coinc_time/Analyzer_coinctime_corrected_%i.png",runno));
  c2->Print(Form("plots/coinc_time/Analyzer_coinctime_uncorrected_%i.png",runno));
  c9->Print(Form("plots/coinc_time/Analyzer_coinctime_corrected_vs_L_%i.png",runno));
  c10->Print(Form("plots/coinc_time/Analyzer_coinctime_corrected_vs_R_%i.png",runno));



  
  //  plots/coinc_time




}
  

