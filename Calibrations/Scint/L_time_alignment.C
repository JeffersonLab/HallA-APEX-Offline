/*
*************************************************************
17/12/20 John Williamson

Performs S2 timing alignment for LHRS.
Uses hits in adjacent panels, and presumes average time difference between hits should be zero. 


Set 8th left PMT offset to zero (arbitary) and base other offsets as difference from this. 

Use difference between left and right pmts for one paddle and its relation to offsets and the hit position in the paddle, y (obtained from track projection)m to obtain difference in left and right offsets for one paddle.  
*************************************************************
 */

#include "file_def.h"
#include "Load_more_rootfiles.C"
#include <tuple>
#include "InputAPEXL.h"


void L_time_alignment(Int_t runno){

  
  TChain* t = Load_more_rootfiles(runno);
  
  //  TChain* t = new TChain("T");
  //  t->Add("/w/work3/home/johnw/Rootfiles/apex_online_20_1_21_tof_4668.root");// TOF + TW
  

  const Int_t NS2Pad = 16;
  
  // Histograms for difference betwen adjacent paddle L-times and R-times
  // establish cuts for adjacent paddles 
  // also plot difference between paddles of sum of left and right time
  // (this should eliminate 

  TH1F* h_adj_l[NS2Pad-1];
  TH1F* h_adj_r[NS2Pad-1];
  TH1F* h_adj_lr[NS2Pad-1];
  
  TString adj_l_name[NS2Pad-1];
  TString adj_r_name[NS2Pad-1];
  TString adj_lr_name[NS2Pad-1];

  TString l_adj_cut[NS2Pad-1];
  TString r_adj_cut[NS2Pad-1];
  TString lr_adj_cut[NS2Pad-1];

  TString lr_adj_alt_cut[NS2Pad-1];
  TString r_adj_alt_cut[NS2Pad-1];


  // cross-talk cuts based on s2_crosstalk.C script
  // Try and distinguish particles that caused signal in both adjacent paddles from those that formed signal in one and through cross-talk resulted in signals in both
  // cut based on projected VDC track in S2 paddle: this should also ensure good events (in which S2 signals are same as those that created tracks in VDC)

  TString l_ct_cut[NS2Pad-1];
  TString r_ct_cut[NS2Pad-1];
  TString lr_ct_cut[NS2Pad-1];



  
  Double_t L_ct_x_l[NS2Pad-1] = {-1.05,-0.84,-0.74,-0.58,-0.42,-1,-0.16,-0.03,0.1,0.24,0.38,0.53,0.66,0.78,0.93};
  Double_t L_ct_x_h[NS2Pad-1] = {-0.95,-0.78,-0.66,-0.45,-0.38,1,-0.-0.075,0.07,0.2,0.32,0.45,0.61,0.73,0.88,0.99};

  
  Double_t R_ct_x_l[NS2Pad-1] = {-1.03,-0.85,-0.73,-0.58,-0.43,-0.3,-0.13,-0.03,0.12,0.26,0.39,0.53,0.66,0.8,0.93};
  Double_t R_ct_x_h[NS2Pad-1] = {-0.87,-0.73,-0.63,-0.50,-0.35,-0.2,-0.08,0.07,0.22,0.34,0.45,0.6,0.74,0.9,1.01};

  
  for(int i=0;i<NS2Pad-1;i++){
    

    adj_l_name[i] = Form("adj_hist_l_%i",i);
    adj_r_name[i] = Form("adj_hist_r_%i",i);
    adj_lr_name[i] = Form("adj_hist_lr_%i",i);

    h_adj_l[i] = new TH1F(adj_l_name[i],Form("L-Arm S2, L-paddles time difference %d and channel %d",i,i+1),40,-20,20);
    h_adj_l[i]->GetXaxis()->SetTitle(Form("(L paddle-%i) - (L paddle-%i), TDC channel",i,i+1));

    
    h_adj_r[i] = new TH1F(adj_r_name[i],Form("L-Arm S2, R-paddles time difference %d and channel %d",i,i+1),40,-20,20);
    h_adj_r[i]->GetXaxis()->SetTitle(Form("(R paddle-%i) - (R paddle-%i), TDC channel",i,i+1));

    h_adj_lr[i] = new TH1F(adj_lr_name[i],Form("L-Arm S2, L+R-paddles time difference %d and channel %d",i,i+1),40,-20,20);
    h_adj_lr[i]->GetXaxis()->SetTitle(Form("(L+R paddle-%i) - (L+R paddle-%i), TDC channel",i,i+1));

    l_adj_cut[i] = Form("L.s2.lt[%i]>500&&L.s2.la_p[%i]>200&&L.s2.lt[%i]>500&&L.s2.la_p[%i]>200",i,i,i+1,i+1);

    r_adj_cut[i] = Form("L.s2.rt[%i]>500&&L.s2.ra_p[%i]>200&&L.s2.rt[%i]>500&&L.s2.ra_p[%i]>200",i,i,i+1,i+1);
    
    r_adj_alt_cut[i] = Form("(L.s2.rt[%i]>500&&L.s2.ra_p[%i]>200)||(L.s2.rt[%i]>500&&L.s2.ra_p[%i]>200)",i,i,i+1,i+1);

    lr_adj_cut[i] = Form("L.s2.lt[%i]>500&&L.s2.la_p[%i]>200&&L.s2.lt[%i]>500&&L.s2.la_p[%i]>200 && L.s2.rt[%i]>500&&L.s2.ra_p[%i]>200&&L.s2.rt[%i]>500&&L.s2.ra_p[%i]>200",i,i,i+1,i+1,i,i,i+1,i+1);

    lr_adj_alt_cut[i] = Form("(L.s2.lt[%i]>500&&L.s2.la_p[%i]>200)||(L.s2.lt[%i]>500&&L.s2.la_p[%i]>200) || (L.s2.rt[%i]>500&&L.s2.ra_p[%i]>200)||(L.s2.rt[%i]>500&&L.s2.ra_p[%i]>200)",i,i,i+1,i+1,i,i,i+1,i+1);


    //    l_ct_cut[i] = Form("L.s2.trx>%f&&L.s2.trx<%f",L_ct_x_l[i],L_ct_x_h[i]);
    l_ct_cut[i] = "abs(L.s2.trx)<1 && abs(L.s2.trdx)<0.07" + PID_cuts + GeneralSieveCut;
    
    //    cout << l_ct_cut[i] << endl;
    r_ct_cut[i] = Form("R.s2.trx>%f&&R.s2.trx<%f",R_ct_x_l[i],R_ct_x_h[i]);
    //    cout << r_ct_cut[i] << endl<<endl;

    lr_ct_cut[i] = Form("L.s2.trx>%f && L.s2.trx<%f && R.s2.trx>%f && R.s2.trx<%f",L_ct_x_l[i],L_ct_x_h[i],R_ct_x_l[i],R_ct_x_h[i]);

   
    }


  // cuts for intercept plots (plot of difference between pmts for one paddle)

  TString l_int_cut[NS2Pad];
  TString r_int_cut[NS2Pad];
  TString lr_int_cut[NS2Pad];

  TString l_int_alt_cut[NS2Pad];  
  

  for(int i=0;i<NS2Pad;i++){
  
    l_int_cut[i] = Form("L.s2.nthit==1 && L.s2.t_pads==%i && L.s2.la_c[%i]>100 && L.s2.ra_c[%i]>100 && abs(L.s2.trx)<1 && abs(L.s2.trdx)<0.07 && L.tr.n==1",i,i,i) + PID_cuts + GeneralSieveCut;

    if( i == 4 || i ==5){
      l_int_alt_cut[i] = Form("L.s2.nthit==1 && L.s2.t_pads==%i && (L.s2.la_c[%i]>100 || L.s2.ra_c[%i]>100) && abs(L.s2.trx)<1 && abs(L.s2.trdx)<0.07 && L.tr.n==1",i,i,i);

    }
    else{
      l_int_alt_cut[i] = Form("L.s2.nthit==1 && L.s2.t_pads==%i && L.s2.la_c[%i]>100 && L.s2.ra_c[%i]>100 && abs(L.s2.trx)<1 && abs(L.s2.trdx)<0.07 && L.tr.n==1",i,i,i);
      // TString r_int_cut[NS2Pad-1];
      // TString lr_int_cut[NS2Pad-1];
    }

    
  }
  



  // perform timing alignment for left and right paddles

  // offsets here record offset between channels (later will set 8th paddle to have zero correction and use offsets to get corrections for other paddles)
  Double_t offset_l[NS2Pad-1], offset_l_err[NS2Pad-1], offset_r[NS2Pad-1], offset_r_err[NS2Pad-1], offset_lr[NS2Pad-1], offset_lr_err[NS2Pad-1];


  TCanvas* c_l_adj =  new TCanvas("c_l_adj","Alignment between left pmts");
  c_l_adj->Divide(4,4);
  TCanvas* c_r_adj =  new TCanvas("c_r_adj","Alignment between right pmts");
  c_r_adj->Divide(4,4);
  TCanvas* c_lr_adj = new TCanvas("c_lr_adj","Alignment between left + right pmts");
  c_lr_adj->Divide(4,4);

  
  // store gaussian used to fit difference
  TF1 *f1_l;
  TF1 *f1_r;
  TF1 *f1_lr;
  Double_t lmax,lpar[3],rmax,rpar[3],lrmax,lrpar[3];

  
  
  for(int i=0;i<NS2Pad-1;i++){

    c_l_adj->cd(i+1);
    
    t->Draw(Form("L.s2.lt[%i]-L.s2.lt[%i]>>adj_hist_l_%i",i,i+1,i),lr_adj_cut[i] + "&&" + l_ct_cut[i],"");
    
    f1_l = new TF1("f1_l","gaus",-10,10);
    
    h_adj_l[i]->Fit(f1_l,"RQ");
    f1_l->GetParameters(&lpar[0]);
    f1_l->Draw("same");
    offset_l[i] = lpar[1];
    offset_l_err[i] = f1_l->GetParError(1);
    
    c_r_adj->cd(i+1);

    //    t->Draw(Form("L.s2.rt[%i]-L.s2.rt[%i]>>adj_hist_r_%i",i,i+1,i),r_adj_cut[i] + "&&" + l_ct_cut[i],"");
    
    if(i==1 || i==5){
      t->Draw(Form("L.s2.rt[%i]-L.s2.rt[%i]>>adj_hist_r_%i",i,i+1,i),r_adj_alt_cut[i] + "&&" + l_ct_cut[i],"");
    }
    else{
      t->Draw(Form("L.s2.rt[%i]-L.s2.rt[%i]>>adj_hist_r_%i",i,i+1,i),r_adj_cut[i] + "&&" + l_ct_cut[i],"");
    }

    f1_r = new TF1("f1_r","gaus",-10,10);
    h_adj_r[i]->Fit(f1_r,"RQ");
    f1_r->GetParameters(&rpar[0]);
    f1_r->Draw("same");
    offset_r[i] = rpar[1];
    offset_r_err[i] = f1_r->GetParError(1);
    


    c_lr_adj->cd(i+1);

    if(i == 0 || i == 3 || i==5){
      t->Draw(Form("((L.s2.lt[%i]+L.s2.rt[%i])-(L.s2.lt[%i]+L.s2.rt[%i]))>>adj_hist_lr_%i",i,i,i+1,i+1,i),lr_adj_alt_cut[i] + "&&" + l_ct_cut[i],"");
    }
    else{
      t->Draw(Form("((L.s2.lt[%i]+L.s2.rt[%i])-(L.s2.lt[%i]+L.s2.rt[%i]))>>adj_hist_lr_%i",i,i,i+1,i+1,i),lr_adj_cut[i] + "&&" + l_ct_cut[i],"");
    }

    f1_lr = new TF1("f1_lr","gaus",-10,10);
    h_adj_lr[i]->Fit(f1_lr,"RQ");
    f1_lr->GetParameters(&rpar[0]);
    f1_lr->Draw("same");
    offset_lr[i] = rpar[1];
    offset_lr_err[i] = f1_lr->GetParError(1);
           
  }

  // Print adjacent paddle plots

  c_l_adj->Print(Form("plots/s2_alignment/L_s2_l_alignment_%i.pdf",runno));
  
  c_r_adj->Print(Form("plots/s2_alignment/L_s2_r_alignment_%i.pdf",runno));
  
  c_lr_adj->Print(Form("plots/s2_alignment/L_s2_lr_alignment_%i.pdf",runno));

   

  /* Use difference in time between left and right pmts of same paddle with y position of hit in paddle (taken from vdc track projection) to get difference in timing offsets between left and right pmts
   */


  TProfile *pdiff[NS2Pad];
  TH2F *hdiff[NS2Pad];


  Int_t s2l_low=2600;
  Int_t s2l_high=3000;
  Int_t s2r_low=2600;
  Int_t s2r_high=3000;
  Int_t s2l_bin=(s2l_high-s2l_low);
  Int_t s2r_bin=(s2r_high-s2r_low);
  Int_t s2diff_low=s2l_low-s2r_high;//s2 diff: Left-Right
  Int_t s2diff_high=s2l_high-s2r_low;
  Int_t s2diff_bin=s2diff_high-s2diff_low;


  Double_t intercept[NS2Pad];

  TCanvas* c_lr_intercept = new TCanvas("c_lr_intercept","Difference between left and right pmts for paddle");
  c_lr_intercept->Divide(4,4);

  for(Int_t i = 0; i<NS2Pad; i++){

    c_lr_intercept->cd(i+1);
    
    hdiff[i] = new TH2F(Form("hdiff[%d]",i),Form("S2 Left - Right vs Y_s2: Paddle %d",i+1),300,-0.3,0.3,s2diff_bin,s2diff_low,s2diff_high);

    pdiff[i] = new TProfile(Form("pdiff[%i]",i),Form("S2 Left - Right vs. Track Projection Y: Paddle %d",i+1),90,-0.3,0.3,s2diff_low,s2diff_high);

    pdiff[i]->GetXaxis()->SetTitle("Y_s2 [m]"); hdiff[i]->GetXaxis()->CenterTitle();
    pdiff[i]->GetYaxis()->SetTitle("TDC Channel"); hdiff[i]->GetYaxis()->CenterTitle();


    //    cout << "l_int_cut[" << i << "] = " << l_int_cut[i] << endl;
    t->Draw(Form("L.s2.lt[%i]-L.s2.rt[%i]:L.s2.try>>pdiff[%i]",i,i,i),l_int_cut[i],"prof");

    pdiff[i]->Fit("pol1","RQ","",-0.125,0.125);
    if(pdiff[i]->GetFunction("pol1")){
      intercept[i] = pdiff[i]->GetFunction("pol1")->GetParameter(0);
    }
    else{
           
      t->Draw(Form("L.s2.lt[%i]-L.s2.rt[%i]:L.s2.try>>pdiff[%i]",i,i,i),l_int_alt_cut[i],"prof");
      pdiff[i]->Fit("pol1","RQ","",-0.125,0.125);
      intercept[i] = pdiff[i]->GetFunction("pol1")->GetParameter(0);
    }
      
    
    pdiff[i]->GetFunction("pol1")->Draw("same");

    //    cout << "For " << i << "th paddle, intercept = " << intercept[i] << endl;

  }

  c_lr_intercept->Print(Form("plots/s2_alignment/L_s2_lr_intercept_%i.pdf",runno));



  // Combine information from time difference between adjacent paddles and between left and right pmts of paddles to get time offsets for all S2 paddles

  // first set 8th S2 paddle left offset to zero
  // this choice is arbitary and is premised on difference between left and right spectrometer S2 being of interest

  Double_t l_coeff[NS2Pad] = {0.0};
  Double_t r_coeff[NS2Pad] = {0.0};


  l_coeff[7] = 0.0;
  r_coeff[7] = l_coeff[7] - intercept[7];


  

  // now loop forward and backward from the 8th paddle
  // can use knowledge of difference between adjacent paddles and between left and right paddles to determine left and right calibration coefficients

  // loop from 9th paddle to 16th
  for(Int_t i = 8; i<16; i++ ){


    Double_t sum = -offset_lr[i-1] + l_coeff[i-1] + r_coeff[i-1];
    Double_t diff = intercept[i];
    
    l_coeff[i] = .5*(sum+diff);
    r_coeff[i] = .5*(sum-diff);

  }
  

  // loop from 7th paddle to 1st
  for(Int_t i = 6; i>-1; i-- ){


    Double_t sum = offset_lr[i] + l_coeff[i+1] + r_coeff[i+1];
    Double_t diff = intercept[i];
    
    l_coeff[i] = .5*(sum+diff);
    r_coeff[i] = .5*(sum-diff);
       
  }


  Double_t lr_coeff[NS2Pad] = {0.0};
  

  cout << endl << endl;

  cout << "| Paddle # \t " << "| L coeff \t " << "| R coeff \t " << "| Intercept \t " << "| diff to next paddle | " << endl;


  ofstream oofile(Form("DB/LHRS_S2_align_t_%i.txt",runno));

  oofile<<"LHRS Scintillator timing offsets: "<<endl<<endl;
  oofile<<"S2 Offsets: "<<endl;
  oofile<<"L.s2.L.off =  ";

  for(Int_t i = 0; i<NS2Pad; i++){
    oofile<<l_coeff[i]<<" ";
  }
  oofile<<endl;
  oofile<<"L.s2.R.off =  ";
  for(Int_t i = 0; i<NS2Pad; i++){
    oofile<<r_coeff[i]<<" ";
  }

  oofile.close();
  
  for(Int_t i = 0; i<NS2Pad; i++){

    if(i<NS2Pad-1){
      cout << "| " << i+1 << "|\t" << l_coeff[i] << "|\t" << r_coeff[i] << "|\t" << intercept[i] << "|\t" << offset_lr[i] << "|" << endl;
    }
    else{
      cout << "| " << i+1 << "|\t" << l_coeff[i] << "|\t" << r_coeff[i] << "|\t" << intercept[i] << "|\t" << "-" << "|" << endl;
    }
  }

  cout << endl << endl << endl;
  for(Int_t i = 0; i<NS2Pad; i++){

    if(i<NS2Pad-1){
      cout << "| " << i+1 << "|\t" << l_coeff[i] << "|\t" << r_coeff[i] << "|\t" << intercept[i] << "|\t" << offset_lr[i] << "|" << endl;
      cout << "| " << i+1 << "|\t" << Form("%.9g",l_coeff[i]) << "|\t" << Form("%.9g",r_coeff[i]) << "|\t" << Form("%.9g",intercept[i]) << "|\t" << Form("%.9g",offset_lr[i]) << "|" << endl;
    }
    else{
      cout << "| " << i+1 << "|\t" << l_coeff[i] << "|\t" << r_coeff[i] << "|\t" << intercept[i] << "|\t" << "-" << "|" << endl;
    }
  }
  

  cout << endl << endl;



  
  ofstream oofile_csv(Form("DB/LHRS_S2_align_t_%i.csv",runno));


  oofile_csv<<"L.s2.L.off" << "," << "L.s2.R.off" << endl;

  for(Int_t i = 0; i<NS2Pad; i++){
    oofile_csv<<l_coeff[i]<<","<<r_coeff[i]<<endl;
  }
  oofile_csv.close();
  


  // draw time distributions before and after correction

  TH1F* h_time_l[NS2Pad];
  TH1F* h_time_r[NS2Pad];
  TH1F* h_time_lr[NS2Pad];

  TH1F* h_time_un_l[NS2Pad];
  TH1F* h_time_un_r[NS2Pad];
  TH1F* h_time_un_lr[NS2Pad];

  Double_t l_coeff_mean = 0.0;
  Double_t r_coeff_mean = 0.0;

  Double_t lr_coeff_mean = 0.0;

  for(int i=0;i<NS2Pad;i++){
    l_coeff_mean += l_coeff[i];
    r_coeff_mean += r_coeff[i];
    lr_coeff[i] = (l_coeff[i] + r_coeff[i])/2;
    lr_coeff_mean += lr_coeff[i];
  }

  l_coeff_mean = l_coeff_mean/NS2Pad;
  r_coeff_mean = r_coeff_mean/NS2Pad;
  //  lr_coeff_mean = (l_coeff_mean + r_coeff_mean)/2.;
  lr_coeff_mean = lr_coeff_mean/NS2Pad;


  cout << "Printing time means: " << endl;
  cout << "l_coeff_mean = " << l_coeff_mean << endl;
  cout << "r_coeff_mean = " << r_coeff_mean << endl;
  cout << "lr_coeff_mean = " << lr_coeff_mean << endl;
  
  Double_t lr_mean_time[NS2Pad] = {0.0};
  Double_t l_mean_time[NS2Pad] = {0.0};
  Double_t r_mean_time[NS2Pad] = {0.0};

  Double_t lr_mean_time_un[NS2Pad] = {0.0};
  Double_t l_mean_time_un[NS2Pad] = {0.0};
  Double_t r_mean_time_un[NS2Pad] = {0.0};

  Double_t pad[NS2Pad] = {0};
  
  for(int i=0;i<NS2Pad;i++){
      
    h_time_l[i] = new TH1F(Form("h_time_l_%i",i),Form("h_time_l_%i",i),57,2753-l_coeff[i],2810-l_coeff[i]);
    t->Draw(Form("L.s2.lt[%i]-%.9g >>h_time_l_%i",i,l_coeff[i],i),Form("L.s2.t_pads[0]==%i",i));    
    h_time_l[i]->Fit("gaus","RQ","",2753-l_coeff[i],2810-l_coeff[i]);
    l_mean_time[i] = h_time_l[i]->GetFunction("gaus")->GetParameter(1);

    h_time_r[i] = new TH1F(Form("h_time_r_%i",i),Form("h_time_r_%i",i),25,2725-r_coeff[i],2750-r_coeff[i]);
    t->Draw(Form("L.s2.rt[%i]-%.9g >>h_time_r_%i",i,r_coeff[i],i),Form("L.s2.t_pads[0]==%i",i));
    //    10 width

    if (i == 13){
      h_time_r[i]->Fit("gaus","RQ","",2725-r_coeff[i],2740-r_coeff[i]);
    }
    else{
      h_time_r[i]->Fit("gaus","RQ","",2725-r_coeff[i],2750-r_coeff[i]);
    }
    // Double_t r_mean = h_time_r[i]->GetMean();
    // h_time_r[i]->Fit("gaus","RQ","",r_mean-20,r_mean+20);
    r_mean_time[i] = h_time_r[i]->GetFunction("gaus")->GetParameter(1);

    h_time_lr[i] = new TH1F(Form("h_time_lr_%i",i),Form("h_time_lr_%i",i),27,2753-lr_coeff[i],2780-lr_coeff[i]);
    //    t->Draw(Form("((L.s2.lt[%i] - %.9g ) + (L.s2.rt[%i] - %.9g ))/2. >>h_time_lr_%i",i,l_coeff[i],i,r_coeff[i],i),Form("L.s2.t_pads[0]==%i",i));
    t->Draw(Form("(((L.s2.lt[%i] ) + (L.s2.rt[%i] )))/2. - %.9g >>h_time_lr_%i",i,i,lr_coeff[i],i),Form("L.s2.t_pads[0]==%i",i));    
    h_time_lr[i]->Fit("gaus","RQ","",2753-lr_coeff[i],2780-lr_coeff[i]);
    lr_mean_time[i] = h_time_lr[i]->GetFunction("gaus")->GetParameter(1);

    
    h_time_un_l[i] = new TH1F(Form("h_time_un_l_%i",i),Form("h_time_un_l_%i",i),57,2753,2810);
    t->Draw(Form("L.s2.lt[%i]>>h_time_un_l_%i",i,i),Form("L.s2.t_pads[0]==%i",i));    
    h_time_un_l[i]->Fit("gaus","RQ","",2753,2810);
    l_mean_time_un[i] = h_time_un_l[i]->GetFunction("gaus")->GetParameter(1);

    h_time_un_r[i] = new TH1F(Form("h_time_un_r_%i",i),Form("h_time_un_r_%i",i),25,2725,2750);
    t->Draw(Form("L.s2.rt[%i] >>h_time_un_r_%i",i,i),Form("L.s2.t_pads[0]==%i",i));    
    h_time_un_r[i]->Fit("gaus","RQ","",2725,2750);
    r_mean_time_un[i] = h_time_un_r[i]->GetFunction("gaus")->GetParameter(1);

    h_time_un_lr[i] = new TH1F(Form("h_time_un_lr_%i",i),Form("h_time_un_lr_%i",i),27,2753,2780);
    t->Draw(Form("(L.s2.lt[%i] + L.s2.rt[%i])/2. >>h_time_un_lr_%i",i,i,i),Form("L.s2.t_pads[0]==%i",i));    
    h_time_un_lr[i]->Fit("gaus","RQ","",2753,2780);
    lr_mean_time_un[i] = h_time_un_lr[i]->GetFunction("gaus")->GetParameter(1);   
    // h_time_r[i] = new TH1F(Form("h_time_r_%i",i),Form("h_time_r_%i",i),25,2725-r_coeff_mean,2750-r_coeff_mean);
    // h_time_lr[i] = new TH1F(Form("h_time_lr_%i",i),Form("h_time_lr_%i",i)27,2753-lr_coeff_mean,2780-lr_coeff_mean);

    // h_time_un_l[i] = new TH1F(Form("h_time_un_l_%i",i),Form("h_time_un_l_%i",i),57,2753,2810);
    // h_time_un_r[i] = new TH1F(Form("h_time_un_r_%i",i),Form("h_time_un_r_%i",i),25,2725,2750);
    // h_time_un_lr[i] = new TH1F(Form("h_time_un_lr_%i",i),Form("h_time_un_lr_%i",i)27,2753,2780);


    pad[i] = i;

  }

  // plot tgraphs of pre and post correction s2 paddle time averages

  
  
  TCanvas* g_canv = new TCanvas("g_canv","Corrections");

  g_canv->Divide(3,4);
  
  g_canv->cd(1);
  TGraph* g_l_t = new TGraph(NS2Pad,pad,l_mean_time);
  g_l_t->SetTitle("Corrected left times");
  g_l_t->Draw("AC*");

  g_canv->cd(2);
  TGraph* g_r_t = new TGraph(NS2Pad,pad,r_mean_time);
  g_r_t->SetTitle("Corrected right times");
  g_r_t->Draw("AC*");

  g_canv->cd(3);
  TGraph* g_lr_t = new TGraph(NS2Pad,pad,lr_mean_time);
  g_lr_t->SetTitle("Corrected left + right times");
  g_lr_t->Draw("AC*");


  g_canv->cd(4);
  TGraph* g_l_t_un = new TGraph(NS2Pad,pad,l_mean_time_un);
  g_l_t_un->SetTitle("Uncorrected left times");
  g_l_t_un->SetLineColor(kRed);
  g_l_t_un->Draw("AC*");
  g_canv->cd(1);
  g_l_t_un->Draw("C* SAME");

  g_canv->cd(5);
  TGraph* g_r_t_un = new TGraph(NS2Pad,pad,r_mean_time_un);
  g_r_t_un->SetTitle("Uncorrected right times");
  g_r_t_un->SetLineColor(kRed);
  g_r_t_un->Draw("AC*");

  g_canv->cd(6);
  TGraph* g_lr_t_un = new TGraph(NS2Pad,pad,lr_mean_time_un);
  g_lr_t_un->SetTitle("Uncorrected left + right times");
  g_lr_t_un->SetLineColor(kRed);
  g_lr_t_un->Draw("AC*");


  g_canv->cd(7);
		  
  TGraph* g_l_corr = new TGraph(NS2Pad,pad,l_coeff);
  g_l_corr->SetTitle("Left corrections");
  g_l_corr->Draw("AC*");

  g_canv->cd(8);
  TGraph* g_r_corr = new TGraph(NS2Pad,pad,r_coeff);
  g_r_corr->SetTitle("Right corrections");
  g_r_corr->Draw("AC*");

  g_canv->cd(9);
  TGraph* g_lr_corr = new TGraph(NS2Pad,pad,lr_coeff);
  g_lr_corr->SetTitle("Left + Right corrections");
  g_lr_corr->Draw("AC*");



  Double_t off[NS2Pad-1] = {0};


  for(int i=0;i<NS2Pad;i++){
    off[i] = i;
  }

  g_canv->cd(10);
		  
  TGraph* g_l_off = new TGraph(NS2Pad-1,off,offset_l);
  g_l_off->SetTitle("Left offsets");
  g_l_off->Draw("AC*");

  g_canv->cd(11);
  TGraph* g_r_off = new TGraph(NS2Pad-1,off,offset_r);
  g_r_off->SetTitle("Right offsets");
  g_r_off->Draw("AC*");

  g_canv->cd(12);
  TGraph* g_lr_off = new TGraph(NS2Pad-1,off,offset_lr);
  g_lr_off->SetTitle("Left + Right offsets");
  g_lr_off->Draw("AC*");
  


  TCanvas* g_canv_intercept = new TCanvas("g_canv_intercept","g_canv_intercept");
  g_canv_intercept->cd(1);
  TGraph* g_intercept = new TGraph(NS2Pad,off,intercept);
  g_intercept->SetTitle("Intercepts");
  g_intercept->Draw("AC*");
		  


  
  /// plot lt, rt and lt + rt for all channels

  TCanvas* time_check_l = new TCanvas("time_check_l","time_check_l");
  time_check_l->Divide(4,4);

  TCanvas* time_check_r = new TCanvas("time_check_r","time_check_r");
  time_check_r->Divide(4,4);

  TCanvas* time_check_lr = new TCanvas("time_check_lr","time_check_lr");
  time_check_lr->Divide(4,4);

  
    for(int i=0;i<NS2Pad;i++){

     
      time_check_l->cd(i+1);
      h_time_l[i]->Draw();

      time_check_r->cd(i+1);
      h_time_r[i]->Draw();

      time_check_lr->cd(i+1);
      h_time_lr[i]->Draw();

    }
  


    
  TCanvas* time_check_un_l = new TCanvas("time_check_un_l","time_check_l (uncorrected)");
  time_check_un_l->Divide(4,4);

  TCanvas* time_check_un_r = new TCanvas("time_check_un_r","time_check_r (uncorrected)");
  time_check_un_r->Divide(4,4);

  TCanvas* time_check_un_lr = new TCanvas("time_check_un_lr","time_check_lr (uncorrected)");
  time_check_un_lr->Divide(4,4);

  
    for(int i=0;i<NS2Pad;i++){

     
      time_check_un_l->cd(i+1);
      h_time_un_l[i]->Draw();

      time_check_un_r->cd(i+1);
      h_time_un_r[i]->Draw();

      time_check_un_lr->cd(i+1);
      h_time_un_lr[i]->Draw();

    }

    

 
    time_check_un_l->Print(Form("plots/s2_alignment/check/L_s2_l_align_uncorr_t_l_%i.pdf",runno));
        
    time_check_un_r->Print(Form("plots/s2_alignment/check/L_s2_r_align_uncorr_t_r_%i.pdf",runno));

    time_check_un_lr->Print(Form("plots/s2_alignment/check/L_s2_lr_align_uncorr_t_lr_%i.pdf",runno));

    time_check_l->Print(Form("plots/s2_alignment/check/L_s2_l_align_t_l_%i.pdf",runno));
        
    time_check_r->Print(Form("plots/s2_alignment/check/L_s2_r_align_t_r_%i.pdf",runno));

    time_check_lr->Print(Form("plots/s2_alignment/check/L_s2_lr_align_t_lr_%i.pdf",runno));


    
    
    //    time_check_un_lr = new TCanvas("time_check_un_lr","time_check_lr (uncorrected)");



    
}
