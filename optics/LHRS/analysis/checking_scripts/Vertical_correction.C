// Script to plot several plots for the horizontal foils (LHRS)
//

#include "../Load_more_rootfiles.C"
#include "InputAPEXL.h"




// using 'triangular wave' of fast raster to conver from current to time

Double_t curr_to_time(Double_t curr, Double_t *par, Bool_t positive = true){

  Double_t period = par[0];
  Double_t amplitude = par[1];
  Double_t base = par[2];


  // positive argument gives phase
  // true if phase means that dI/dt (derivative of current wrt time) is positive

  Double_t phase_add = 0;
  Double_t phase_mult= 1;

  if(positive){
    phase_add =  (0.5*period); // was -
    phase_mult = -1;
  }

  Double_t time = phase_add + phase_mult*(period/(2*TMath::Pi())) * TMath::ASin( TMath::Sin(TMath::Pi()*((curr-base))/(2*amplitude)  )  );

  if (time < -0.75* period){
    time = 0.25*period - TMath::Abs(time + 0.75*period);
  }

  return time;
  //  return phase_add + phase_mult*(period/(2*TMath::Pi())) * TMath::ASin( TMath::Sin(TMath::Pi()*((curr-base))/(2*amplitude)  )  );
  
}

Double_t time_to_curr(Double_t time, Double_t *par){

  Double_t period = par[0];
  Double_t amplitude = par[1];
  Double_t base = par[2];
  
  
  return ( ((2*amplitude)/TMath::Pi()) * TMath::ASin( TMath::Sin(2*TMath::Pi()*(time)/period)  ) + base );
  
}


Double_t peak(Double_t *x, Double_t *par) {
  return  par[0]*TMath::Exp(-(((x[0]-par[1])*(x[0]-par[1]))/(2*par[2]*par[2])));
}

Double_t bg(Double_t *x, Double_t *par) {
  return par[0]+ (par[1]*x[0]);
}

// sum of peaks and backgorund
Double_t overall(Double_t *x, Double_t *par) {
  return peak(x,par) + peak(x,&par[3]) + bg(x,&par[6]);
}



void Vertical_correction(){


  // load in root files
  
  Int_t run_nos[] = {4766};
  // Int_t run_nos[] = {4776,4777};
  
  TChain* T;


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


  // set histogram limits

  // raster limits
  Double_t Lrby_l = 0.0014;
  Double_t Lrby_h = 0.0030;

  Double_t Lrbx_l = -0.0022;
  Double_t Lrbx_h = -0.0011;

  // BPM limits
  Double_t BPM_lim_l = 0.000;
  Double_t BPM_lim_h = 0.003;

  // raster current limits
  Double_t Lrbx_curr_l = 26500;
  Double_t Lrbx_curr_h = 33000;


    
  // reactz limits
  Double_t reactz_l = -0.4;
  Double_t reactz_h = 0.4;

  // Focal plane limits
  
  Double_t y_FP_l = -0.05;
  Double_t y_FP_h = 0.05;

  Double_t ph_FP_l = -50;
  Double_t ph_FP_h = 40;

  Double_t th_FP_l = -80;
  Double_t th_FP_h = 80;

  // target cuts

  Double_t ph_tg_l = -40;
  Double_t ph_tg_h = 40;
  
  Double_t th_tg_l = -80;
  Double_t th_tg_h = 80;
  
  // define raster cuts
  Double_t c_Lrby1_l = 0.00191;
  Double_t c_Lrby1_h = 0.00209;

  Double_t c_Lrby2_l = 0.00240;
  Double_t c_Lrby2_h = 0.00258;

  // define BPM cuts (work in addition to Raster as phase cuts)
  Double_t c_BPMB1_l = -0.0001;
  Double_t c_BPMB1_h = 0.0005;

  Double_t c_BPMB2_l = 0.0043;
  Double_t c_BPMB2_h = 0.0049;

  // Double_t c_BPMB2_l = -1;
  // Double_t c_BPMB2_h = 1;
  
  
  // define reactz cut for V2
  Double_t c_zV2_l = -0.1;
  Double_t c_zV2_h = 0.05;

  // define V2 and raster cuts

  TCut V2_zcut = Form("reactz>%f && reactz<%f",c_zV2_l,c_zV2_h);

  // lower raster cut (in y)
  TCut rast_l = Form("Lrb.y>%f && Lrb.y<%f",c_Lrby1_l,c_Lrby1_h);

  // greater raster cut (in y)
  TCut rast_h = Form("Lrb.y>%f && Lrb.y<%f",c_Lrby2_l,c_Lrby2_h);

  // lower BPMB y
  TCut BPMB_l = Form("Lrb.BPMB.y>%f && Lrb.BPMB.y<%f",c_BPMB1_l,c_BPMB1_h);

  // higher BPMB y
  TCut BPMB_h = Form("Lrb.BPMB.y>%f && Lrb.BPMB.y<%f",c_BPMB2_l,c_BPMB2_h);

  rast_l += BPMB_h;

  rast_h += BPMB_l;
  
  // cuts from InputAPEXl.h
  TCut GenrealCut = GeneralSieveCut + PID_cuts + FP_cuts;


  // plot X curent versus BPMB x (gives phase info)
  
  TCanvas *c1 = new TCanvas("c1"); 

  cout << GenrealCut << endl << endl;


  // define cut for positive and negative phase

  Double_t l_x1 = 0.0013;
  Double_t l_y1 = 32000;
  Double_t l_x2 = 0.00162;
  Double_t l_y2 = 27400;
  

  TLine* phase_l = new TLine(l_x1,l_y1,l_x2,l_y2);
  phase_l->SetLineColor(kRed);

  

  TH2F* h_rcurr = new TH2F("h_rcurr","raster current vs BPMB X",1000,BPM_lim_l,BPM_lim_h,1000,Lrbx_curr_l,Lrbx_curr_h);
  h_rcurr->GetXaxis()->SetTitle("BPMB X [m]");
  h_rcurr->GetYaxis()->SetTitle("Raster X current [ADC channels]");
  h_rcurr->GetYaxis()->SetTitleOffset(2.0);
    
  T->Draw("Lrb.Raster2.rawcur.x:Lrb.BPMB.x>>h_rcurr","","colz");

  phase_l->Draw("same");


  // form cut from line seperating postive and negative phase regions
  
  Double_t p_slope = (l_y2-l_y1)/(l_x2-l_x1);

  Double_t p_off = l_y1 - p_slope*l_x1;


  cout << "p_slope = " << p_slope << endl;
  cout << "p_off = " << p_off << endl;


  // positive phase
  TCut phase_p = Form("(Lrb.Raster2.rawcur.x - %f*Lrb.BPMB.x)>%f",p_slope,p_off);

  // negative phase
  TCut phase_n = Form("(Lrb.Raster2.rawcur.x - %d*Lrb.BPMB.x)<%d",p_slope,p_off);

  cout << "phase_p = " << phase_p << endl;


  TCanvas* c2 = new TCanvas("c2","positve phase Current");

  TH1F* h_xcurr_p = new TH1F("h_xcurr_p","raster current X",500,Lrbx_curr_l,Lrbx_curr_h);
  h_xcurr_p->GetXaxis()->SetTitle("Raster X current [ADC channels]");
   
  T->Draw("Lrb.Raster2.rawcur.x>>h_xcurr_p",phase_p);

  // fit positive raster current to get peak position


  Int_t max_bin = h_xcurr_p->GetMaximumBin();
    
  cout << "Max bin is at " << h_xcurr_p->GetBinCenter(max_bin) << endl;

  Double_t mean = 29000;
  
  Double_t min = mean - 1000;
  Double_t max = mean + 1000;

  h_xcurr_p->Fit("gaus","Q","",min,max);

  TF1* gaus1 = h_xcurr_p->GetFunction("gaus");

  min = gaus1->GetParameter(1) - 0.5*gaus1->GetParameter(2);
  max = gaus1->GetParameter(1) + 0.5*gaus1->GetParameter(2);

  h_xcurr_p->Fit("gaus","Q","",min,max);
  gaus1 = h_xcurr_p->GetFunction("gaus");

  Double_t mean_1 = gaus1->GetParameter(1);
  cout << "mean 1 = " << gaus1->GetParameter(1) << endl;
  
  
  TCanvas* c3 = new TCanvas("c3","negatve phase Current");

  TH1F* h_xcurr_n = new TH1F("h_xcurr_n","raster current X",500,Lrbx_curr_l,Lrbx_curr_h);
  h_xcurr_n->GetXaxis()->SetTitle("Raster X current [ADC channels]");
   
  T->Draw("Lrb.Raster2.rawcur.x>>h_xcurr_n",phase_n);

  // fit negative-phase raster current

  mean = 29000;
  
  min = mean - 1000;
  max = mean + 1000;


  h_xcurr_n->Fit("gaus","Q","",min,max);

  TF1* gaus2 = h_xcurr_n->GetFunction("gaus");

  min = gaus2->GetParameter(1) - 0.5*gaus2->GetParameter(2);
  max = gaus2->GetParameter(1) + 0.5*gaus2->GetParameter(2);

  h_xcurr_n->Fit("gaus","Q","",min,max);
  gaus2 = h_xcurr_n->GetFunction("gaus");

  Double_t mean_2 = gaus2->GetParameter(1);
  cout << "mean_2 = " << mean_2 << endl;
  

  TCanvas* c4 = new TCanvas("c4","test");

  T->Draw(Form("Lrb.Raster2.rawcur.x - %f*Lrb.BPMB.x>>h1(1000,4e4,6.2e4)",p_slope));


  // convert found difference in peaks to time difference

  
  Double_t curr_corr = TMath::Abs(mean_2 - mean_1)/2.0;

  Double_t period = 1/(25e3); // 25khz

  Double_t amplitude = 2250;
  Double_t base = 2.975e4; // baseline current

  
  // set wave parameters for converting currents to time (and vice versa)
  Double_t wave_pars[3] = {period,amplitude,base};

  Double_t time_1 = curr_to_time(mean_1,wave_pars,false);

  Double_t time_1_corr = curr_to_time(mean_1+curr_corr,wave_pars,false);

  Double_t time_1_diff = time_1 - time_1_corr;


  
  Double_t time_2 = curr_to_time(mean_2,wave_pars,false);

  Double_t time_2_corr = curr_to_time(mean_2+curr_corr,wave_pars,false);

  Double_t time_2_diff = time_2 - time_2_corr;

  
  cout << "time_1 = " << time_1 << endl;
  cout << "time_2 = " << time_2 << endl;

  
  cout << "time_1_diff = " << time_1_diff << endl;
  cout << "time_2_diff = " << time_2_diff << endl;



  
   
   
}
  
