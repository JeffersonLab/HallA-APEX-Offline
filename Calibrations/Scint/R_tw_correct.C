/*
*************************************************************
10/1/21 John Williamson

Performs timewalk corrections for RHRS.

Uses common Hall A technique for timewalk calibration as stated in Hall A analyzer which cites:
  // Traditional correction according to
  // J.S.Brown et al., NIM A221, 503 (1984); T.Dreyer et al., NIM A309, 184 (1991)
  // Assume that for a MIP (peak ~2000 ADC channels) the correction is 0

This correction involves two coeffecients: one for a 'Minimising Ionising Particle' (MIP) whose timewalk offset is assumed to be zero, this is arbitary but is set at reasonable (and potentially different) value for S0 and S2. The second coeffecients is a mutliplicative factor and this is what must be calibrated.
*************************************************************
 */

#include "file_def.h"
#include "Load_more_rootfiles.C"


TF1* f_cb;
TF1* f_line;

Double_t finter(Double_t *x, Double_t *par) {
  return TMath::Abs(f_cb->EvalPar(x,par) - f_line->EvalPar(x,par));
}


void R_tw_correct(Int_t runno){

  TChain* T = Load_more_rootfiles(runno);


  // set levels of MIP for s0 and s2 (this must match definition in DB)

  // const Double_t s0_MIP = 500.0;
  // const Double_t s2_MIP = 200.0;



  
  const Int_t NS2Pad = 16;

  Double_t s0_L = 0.0;
  Double_t s0_R = 0.0;
  
  Double_t s2_L[NS2Pad] = {0.0};
  Double_t s2_R[NS2Pad] = {0.0};

  TCanvas *c1 = new TCanvas("c1","s0 timing");
  c1->Divide(2,1);
  c1->cd(1);



  TH1F* h_R_s0_t = new TH1F("h_R_s0_t","s0 left + right time",25,4.98e3,5.03e3);

  
  T->Draw("(R.s0.lt + R.s0.rt)>>h_R_s0_t","R.tr.n==1 && (R.ps.asum_c+0.9**R.sh.asum_c)>800 && R.ps.asum_c>350 && R.cer.asum_c>200 && R.s0.nthit==1","colz");

  
  Double_t s0_t_mean = h_R_s0_t->GetMean();

  cout << "s0_t_mean = " << s0_t_mean << endl;

  c1->cd(2);
  TH1F* h_R_s0_ADC = new TH1F("h_R_s0_ADC","s0 left + right TW effect",300,0,0.16);

  T->Draw("((1/TMath::Sqrt(R.s0.la_p)) + (1/TMath::Sqrt(R.s0.ra_p)))>>h_R_s0_ADC","R.tr.n==1 && (R.ps.asum_c+0.9**R.sh.asum_c)>800 && R.ps.asum_c>350 && R.cer.asum_c>200 && R.s0.nthit==1");

  
  TF1* f_cb_s0 = new TF1("f_cb_s0","crystalball",0,0.16);

  Double_t max_ADC_s0 = h_R_s0_ADC->GetBinCenter(h_R_s0_ADC->GetMaximumBin());
    
  Int_t max_entriess_s0 = h_R_s0_ADC->GetBinContent(h_R_s0_ADC->GetMaximumBin());
  
  f_cb_s0->SetParameters(max_entriess_s0,max_ADC_s0,0.0048,0.736,2.7e7);

  h_R_s0_ADC->Fit("f_cb_s0","QR","",0,0.16);

  Double_t s0_mean = f_cb_s0->GetParameter(1);
  Double_t s0_sigma = f_cb_s0->GetParameter(2);

  // check 2 sigma above and below mean and check which has more bin entries

  Double_t s0_lower = s0_mean-1.5*s0_sigma;
  Double_t s0_higher = s0_mean+1.5*s0_sigma;
    
  Int_t l_check = h_R_s0_ADC->GetBinContent(h_R_s0_ADC->FindBin(s0_lower));
    
  Int_t h_check = h_R_s0_ADC->GetBinContent(h_R_s0_ADC->FindBin(s0_higher));

  cout << "s0_lower = " << s0_lower << endl;
  cout << "s0_higher = " << s0_higher << endl;

    
  Double_t s0_ADC_cutoff = (l_check>h_check) ? l_check : h_check;

  TF1* f_line_ADC_s0= new TF1("f_line_ADC_s0","pol1",0,0.16);

  f_line_ADC_s0->SetParameters(s0_ADC_cutoff,0);


  f_cb = f_cb_s0;
  f_line = f_line_ADC_s0;
    

  TF1* s0_ADC_intersection = new TF1("s0_ADC_intersection","abs(f_line_ADC_s0 - f_cb_s0)",s0_lower,s0_higher);


  s0_ADC_intersection->SetLineColor(kOrange);
  s0_ADC_intersection->Draw("same");
    
  f_line_ADC_s0->SetLineColor(kGreen+2);
  f_line_ADC_s0->Draw("same");


  Double_t s0_ADC_l_cut_val = s0_ADC_intersection->GetMinimumX(s0_lower,s0_mean);
  Double_t s0_ADC_h_cut_val = s0_ADC_intersection->GetMinimumX(s0_mean,s0_higher);

  //s0_ADC_cut = new TLine(s0_lower,0,s0_lower,s0_ADC_cutoff);
  TLine*  s0_ADC_l_cut = new TLine(s0_ADC_l_cut_val,0,s0_ADC_l_cut_val,s0_ADC_cutoff);
  TLine* s0_ADC_h_cut = new TLine(s0_ADC_h_cut_val,0,s0_ADC_h_cut_val,s0_ADC_cutoff);

  TLine* s0_ADC_l_cut_t = new TLine(s0_ADC_l_cut_val,-0.05e-7,s0_ADC_l_cut_val,0.05e-7);
  TLine* s0_ADC_h_cut_t = new TLine(s0_ADC_h_cut_val,-0.05e-7,s0_ADC_h_cut_val,0.05e-7);
 

  s0_ADC_l_cut->SetLineColor(kMagenta);
  s0_ADC_h_cut->SetLineColor(kMagenta);

  s0_ADC_l_cut->Draw("same");
  s0_ADC_h_cut->Draw("same");
    
  cout << " Low cut = " << s0_lower << ", high cut = " << s0_ADC_h_cut_val << endl;
  cout << " Low cut = " << s0_ADC_l_cut_val << ", high cut = " << s0_ADC_h_cut_val << endl;
    
  
  

  
  TCanvas *c2 = new TCanvas("c2","s0 timing vs TW effects");

  c2->Divide(2,1);
  c2->cd(1);
  
  TH2F* h_R_s0_t_ADC = new TH2F("h_R_s0_t_ADC","s0 left + right time vs TW effect",50,0,0.16,60,-35,35);

  TString draw_s = Form("((R.s0.lt + R.s0.rt) - %.9g):((1/TMath::Sqrt(R.s0.la_p)) + (1/TMath::Sqrt(R.s0.ra_p)))>>h_R_s0_t_ADC",s0_t_mean);

  cout << "draw string = " << draw_s << endl;
  
  T->Draw(Form("((R.s0.lt + R.s0.rt) - %.9g):((1/TMath::Sqrt(R.s0.la_p)) + (1/TMath::Sqrt(R.s0.ra_p)))>>h_R_s0_t_ADC",s0_t_mean),"R.tr.n==1 && (R.ps.asum_c+0.9*R.sh.asum_c)>800 && R.ps.asum_c>350 && R.cer.asum_c>200 && R.s0.nthit==1","colz");


  c2->cd(2);
  TProfile* hprof_R_s0_t_ADC = new TProfile("hprof_R_s0_t_ADC","s0 left + right time vs TW effect",50,0,0.16,-35,35);

  T->Draw(Form("((R.s0.lt + R.s0.rt) - %.9g):((1/TMath::Sqrt(R.s0.la_p)) + (1/TMath::Sqrt(R.s0.ra_p)))>>hprof_R_s0_t_ADC",s0_t_mean),"R.tr.n==1 && (R.ps.asum_c+0.9*R.sh.asum_c)>800 && R.ps.asum_c>350 && R.cer.asum_c>200 && R.s0.nthit==1","colz");


    
  hprof_R_s0_t_ADC->Fit("pol1","QR","",0,0.16);
  
  Double_t s0_offset =  hprof_R_s0_t_ADC->GetFunction("pol1")->GetParameter(0);
  Double_t s0_k = hprof_R_s0_t_ADC->GetFunction("pol1")->GetParameter(1);

  // can now derive ADC_MIP for s0 from these parameters

  Double_t s0_MIP = (4*TMath::Power(s0_k,2))/TMath::Power(s0_offset,2);


  cout << "s0 k = " << s0_k << ", s0 ADC_MIP = " << s0_MIP << endl;
  
  
  
  
  //  T->Draw(Form("(R.s0.rt - %f):(1/TMath::Sqrt(R.s0.la_p))>>h2(300,0,0.14,300,-0.25e-7,0.25e-7)","R.tr.n==1 && (R.ps.asum_c+0.9**R.sh.asum_c)>800 && R.ps.asum_c>350 && R.cer.asum_c>200 && R.s0.nthit==1","colz",s0_t_mean));




  TH1F* h_R_s2_t[NS2Pad];
  TH1F* h_R_s2_ADC[NS2Pad];
  TF1* f_cb_ADC[NS2Pad];
  TH2F* h_R_s2_t_ADC[NS2Pad];
  TProfile* hprof_R_s2_t_ADC[NS2Pad];


  TH1F* h_R_s2_t_l[NS2Pad];
  TH1F* h_R_s2_ADC_l[NS2Pad];
  TH1F* h_R_s2_ADC_l_check[NS2Pad];
  TF1* f_cb_ADC_l[NS2Pad];
  TF1* f_line_ADC_l[NS2Pad];
  TF1* s2_l_ADC_intersection[NS2Pad];
  TF1* s2_l_ADC_intersection_alt[NS2Pad];
  TLine* s2_ADC_l_cut[NS2Pad];
  TLine* s2_ADC_h_cut[NS2Pad];
  TLine* s2_ADC_l_cut_t[NS2Pad];
  TLine* s2_ADC_h_cut_t[NS2Pad];
  Double_t s2_offset[NS2Pad];
  Double_t s2_k[NS2Pad];
  Double_t s2_MIP[NS2Pad];
  TH2F* h_R_s2_t_ADC_l[NS2Pad];
  TProfile* hprof_R_s2_t_ADC_l[NS2Pad];

  TH1F* h_R_s2_t_r[NS2Pad];
  TH1F* h_R_s2_ADC_r[NS2Pad];
  TF1* f_cb_ADC_r[NS2Pad];
  TH2F* h_R_s2_t_ADC_r[NS2Pad];
  TProfile* hprof_R_s2_t_ADC_r[NS2Pad];
  

  TCanvas* c3[NS2Pad];
  TCanvas* c4[NS2Pad];

  Double_t s2_mean[NS2Pad];
  Double_t s2_mean_l[NS2Pad];
  Double_t s2_mean_r[NS2Pad];

  for(Int_t i = 0; i<NS2Pad; i++){

    h_R_s2_t[i] = new TH1F(Form("h_R_s2_t_%i",i),Form("s2 left + right time, %i",i),30,5.06e3,5.12e3);

  h_R_s2_t_l[i] = new TH1F(Form("h_R_s2_t_l_%i",i),Form("s2 left time, %i",i),30,2.52e3,2.57e3);

    h_R_s2_t_r[i] = new TH1F(Form("h_R_s2_t_r_%i",i),Form("s2 right time, %i",i),30,2.52e3,2.57e3);


    h_R_s2_ADC[i] = new TH1F(Form("h_R_s2_ADC_%i",i),Form("s2 left + right TW effect, %i",i),300,0,0.16);

    h_R_s2_ADC_l[i] = new TH1F(Form("h_R_s2_ADC_l_%i",i),Form("s2 left TW effect, %i",i),300,0,0.16);

    h_R_s2_ADC_r[i] = new TH1F(Form("h_R_s2_ADC_r_%i",i),Form("s2 right TW effect, %i",i),300,0,0.16);


    //    f_cb_ADC[i] = new TF1(Form("f_cb_ADC_%i",i),"crystalball",0,0.16);    
    
    f_cb_ADC_l[i] = new TF1(Form("f_cb_ADC_l_%i",i),"crystalball",0,0.16);

    f_line_ADC_l[i] = new TF1(Form("f_line_ADC_l_%i",i),"pol1",0,0.16);
    

    c3[i] = new TCanvas(Form("c3_%i",i),Form("S2 timing, %i",i));
    c3[i]->Divide(2,3);

    c3[i]->cd(1);

    T->Draw(Form("(R.s2.lt[%i] + R.s2.rt[%i])>>h_R_s2_t_%i",i,i,i),Form("R.tr.n==1 && (R.ps.asum_c+0.9**R.sh.asum_c)>800 && R.ps.asum_c>350 && R.cer.asum_c>200 && R.s2.nthit==1 && R.s2.la_c[%i]>100 && R.s2.ra_c[%i]>100",i,i),"colz");

    s2_mean[i] = h_R_s2_t[i]->GetMean();

    
    c3[i]->cd(2);

    T->Draw(Form("((1/TMath::Sqrt(R.s2.la_p[%i])) + (1/TMath::Sqrt(R.s2.ra_p[%i])))>>h_R_s2_ADC_%i",i,i,i),Form("R.tr.n==1 && (R.ps.asum_c+0.9**R.sh.asum_c)>800 && R.ps.asum_c>350 && R.cer.asum_c>200 && R.s2.nthit==1 && R.s2.la_c[%i]>100 && R.s2.ra_c[%i]>100",i,i),"colz");

    // f_cb_ADC[i]->SetParameters(1, 0, 0.1, 0.02, 0.5);

    // h_R_s2_ADC[i]->Fit(Form("f_cb_ADC_%i",i),"QR","",0,0.16);
    
    
    c3[i]->cd(3);

    T->Draw(Form("(R.s2.lt[%i])>>h_R_s2_t_l_%i",i,i),Form("R.tr.n==1 && (R.ps.asum_c+0.9**R.sh.asum_c)>800 && R.ps.asum_c>350 && R.cer.asum_c>200 && R.s2.nthit==1 && R.s2.la_c[%i]>100",i),"colz");

    s2_mean_l[i] = h_R_s2_t_l[i]->GetMean();


    c3[i]->cd(4);

    T->Draw(Form("((1/TMath::Sqrt(R.s2.la_p[%i])))>>h_R_s2_ADC_l_%i",i,i),Form("R.tr.n==1 && (R.ps.asum_c+0.9**R.sh.asum_c)>800 && R.ps.asum_c>350 && R.cer.asum_c>200 && R.s2.nthit==1 && R.s2.la_c[%i]>100",i),"colz");

    Double_t max_ADC_l = h_R_s2_ADC_l[i]->GetBinCenter(h_R_s2_ADC_l[i]->GetMaximumBin());
    
    Int_t max_entriess_l = h_R_s2_ADC_l[i]->GetBinContent(h_R_s2_ADC_l[i]->GetMaximumBin());


    cout << "For " << i << " max ADC = " << max_ADC_l << ", with " << max_entriess_l << endl;
      
    f_cb_ADC_l[i]->SetParameters(max_entriess_l,max_ADC_l,0.0048,0.736,2.7e7);

    h_R_s2_ADC_l[i]->Fit(Form("f_cb_ADC_l_%i",i),"QR","",0,0.16);

    Double_t s2_l_mean = f_cb_ADC_l[i]->GetParameter(1);
    Double_t s2_l_sigma = f_cb_ADC_l[i]->GetParameter(2);

    // check 2 sigma above and below mean and check which has more bin entries

    Double_t s2_l_lower = s2_l_mean-1.5*s2_l_sigma;
    Double_t s2_l_higher = s2_l_mean+1.5*s2_l_sigma;
    
    Int_t l_check = h_R_s2_ADC_l[i]->GetBinContent(h_R_s2_ADC_l[i]->FindBin(s2_l_lower));
    
    Int_t h_check = h_R_s2_ADC_l[i]->GetBinContent(h_R_s2_ADC_l[i]->FindBin(s2_l_higher));

    cout << "s2_l_lower = " << s2_l_lower << endl;
    cout << "s2_l_higher = " << s2_l_higher << endl;

    
    Double_t s2_l_ADC_cutoff = (l_check>h_check) ? l_check : h_check;



    f_line_ADC_l[i]->SetParameters(s2_l_ADC_cutoff,0);


    f_cb = f_cb_ADC_l[i];
    f_line = f_line_ADC_l[i];
    
    s2_l_ADC_intersection[i] = new TF1(Form("s2_l_ADC_intersection_%i",i),finter,s2_l_mean,s2_l_higher,7);

    //    s2_l_ADC_intersection_alt[i] = new TF1(Form("s2_l_ADC_intersection_alt_%i",i),Form("abs(f_cb_ADC_l_%i - f_line_ADC_l_%i)",i,i),s2_l_mean,s2_l_higher);
    s2_l_ADC_intersection_alt[i] = new TF1(Form("s2_l_ADC_intersection_alt_%i",i),Form("abs(f_line_ADC_l_%i - f_cb_ADC_l_%i)",i,i),s2_l_lower,s2_l_higher);

    s2_l_ADC_intersection[i]->SetLineColor(kBlack);
    s2_l_ADC_intersection[i]->Draw("same");

    s2_l_ADC_intersection_alt[i]->SetLineColor(kOrange);
    s2_l_ADC_intersection_alt[i]->Draw("same");
    
    f_line_ADC_l[i]->SetLineColor(kGreen+2);
    f_line_ADC_l[i]->Draw("same");


    Double_t s2_l_ADC_l_cut = s2_l_ADC_intersection_alt[i]->GetMinimumX(s2_l_lower,s2_l_mean);
    Double_t s2_l_ADC_h_cut = s2_l_ADC_intersection_alt[i]->GetMinimumX(s2_l_mean,s2_l_higher);

    //s2_ADC_l_cut[i] = new TLine(s2_l_lower,0,s2_l_lower,s2_l_ADC_cutoff);
    s2_ADC_l_cut[i] = new TLine(s2_l_ADC_l_cut,0,s2_l_ADC_l_cut,s2_l_ADC_cutoff);
    s2_ADC_h_cut[i] = new TLine(s2_l_ADC_h_cut,0,s2_l_ADC_h_cut,s2_l_ADC_cutoff);

    s2_ADC_l_cut_t[i] = new TLine(s2_l_ADC_l_cut,-0.05e-7,s2_l_ADC_l_cut,0.05e-7);
    s2_ADC_h_cut_t[i] = new TLine(s2_l_ADC_h_cut,-0.05e-7,s2_l_ADC_h_cut,0.05e-7);
 

    s2_ADC_l_cut[i]->SetLineColor(kMagenta);
    s2_ADC_h_cut[i]->SetLineColor(kMagenta);

    s2_ADC_l_cut[i]->Draw("same");
    s2_ADC_h_cut[i]->Draw("same");
    
    cout << " Low cut = " << s2_l_lower << ", high cut = " << s2_l_ADC_h_cut << endl;
    cout << " Low cut = " << s2_l_ADC_l_cut << ", high cut = " << s2_l_ADC_h_cut << endl;
    
    //    Double_t h_or_l = (l_check>h_check) ? l_check : h_check;


    // if(l_check>h_check){
      
      
    // }
    // else if(l_check>h_check){

    // }
    // else{


    // }
    

    
    
    cout << "h_check = " << h_check << endl;
    cout << "l_check = " << l_check << endl;
    cout << "s2_l_ADC_cutoff = " << s2_l_ADC_cutoff << endl << endl;;

    
    
      
    c3[i]->cd(5);

    T->Draw(Form("(R.s2.rt[%i])>>h_R_s2_t_r_%i",i,i),Form("R.tr.n==1 && (R.ps.asum_c+0.9**R.sh.asum_c)>800 && R.ps.asum_c>350 && R.cer.asum_c>200 && R.s2.nthit==1 && R.s2.ra_c[%i]>100",i),"colz");

    s2_mean_r[i] = h_R_s2_t_r[i]->GetMean();

    c3[i]->cd(6);
	
    T->Draw(Form("((1/TMath::Sqrt(R.s2.ra_p[%i])))>>h_R_s2_ADC_r_%i",i,i),Form("R.tr.n==1 && (R.ps.asum_c+0.9**R.sh.asum_c)>800 && R.ps.asum_c>350 && R.cer.asum_c>200 && R.s2.nthit==1 && R.s2.ra_c[%i]>100",i));


    c4[i] = new TCanvas(Form("c4_%i",i),Form("s2 timing vs TW effects, %i",i));

    c4[i]->Divide(2,3);
    c4[i]->cd(1);
    
    h_R_s2_t_ADC[i] = new TH2F(Form("h_R_s2_t_ADC_%i",i),Form("s2 left + right time vs TW effect (with correction), %i",i),50,0,0.16,20,-0.015e3,0.015e3);

    T->Draw(Form("((R.s2.lt[%i] + R.s2.rt[%i]) - %.9g):((1/TMath::Sqrt(R.s2.la_p[%i])) + (1/TMath::Sqrt(R.s2.ra_p[%i])))>>h_R_s2_t_ADC_%i",i,i,s2_mean[i],i,i,i),Form("R.tr.n==1 && (R.ps.asum_c+0.9**R.sh.asum_c)>800 && R.ps.asum_c>350 && R.cer.asum_c>200 && R.s2.nthit==1 && R.s2.la_c[%i]>100 && R.s2.ra_c[%i]>100",i,i),"colz");

    
    c4[i]->cd(2);
    
    hprof_R_s2_t_ADC[i] = new TProfile(Form("hprof_R_s2_t_ADC_%i",i),Form("s2 left + right time vs TW effect  (with correction), %i",i),50,0,0.16,-0.015e3,0.015e3);

    T->Draw(Form("((R.s2.lt[%i] + R.s2.rt[%i]) - %.9g):((1/TMath::Sqrt(R.s2.la_p[%i])) + (1/TMath::Sqrt(R.s2.ra_p[%i])))>>hprof_R_s2_t_ADC_%i",i,i,s2_mean[i],i,i,i),Form("R.tr.n==1 && (R.ps.asum_c+0.9**R.sh.asum_c)>800 && R.ps.asum_c>350 && R.cer.asum_c>200 && R.s2.nthit==1 && R.s2.la_c[%i]>100 && R.s2.ra_c[%i]>100",i,i),"colz");

   
    

    c4[i]->cd(3);
    
    h_R_s2_t_ADC_l[i] = new TH2F(Form("h_R_s2_t_ADC_l_%i",i),Form("s2 left time vs TW effect, %i",i),50,0,0.16,20,-0.015e3,0.015e3);

    T->Draw(Form("((R.s2.lt[%i]) - %.9g):((1/TMath::Sqrt(R.s2.la_p[%i])))>>h_R_s2_t_ADC_l_%i",i,s2_mean_l[i],i,i),Form("R.tr.n==1 && (R.ps.asum_c+0.9**R.sh.asum_c)>800 && R.ps.asum_c>350 && R.cer.asum_c>200 && R.s2.nthit==1 && R.s2.la_c[%i]>100",i),"colz");

    
    c4[i]->cd(4);
    
    hprof_R_s2_t_ADC_l[i] = new TProfile(Form("hprof_R_s2_t_ADC_l_%i",i),Form("s2 left time vs TW effect, %i",i),50,0,0.16,-0.015e3,0.015e3);

    T->Draw(Form("((R.s2.lt[%i]) - %.9g):((1/TMath::Sqrt(R.s2.la_p[%i])))>>hprof_R_s2_t_ADC_l_%i",i,s2_mean_l[i],i,i),Form("R.tr.n==1 && (R.ps.asum_c+0.9**R.sh.asum_c)>800 && R.ps.asum_c>350 && R.cer.asum_c>200 && R.s2.nthit==1 && R.s2.la_c[%i]>100",i),"colz");

    s2_ADC_l_cut_t[i]->SetLineColor(kMagenta);
    s2_ADC_h_cut_t[i]->SetLineColor(kMagenta);
    
    s2_ADC_l_cut_t[i]->Draw("same");
    s2_ADC_h_cut_t[i]->Draw("same");

    

    hprof_R_s2_t_ADC_l[i]->Fit("pol1","QR","",s2_l_ADC_l_cut,s2_l_ADC_h_cut);

    if(hprof_R_s2_t_ADC_l[i]->GetFunction("pol1")){
    
      s2_offset[i] =  hprof_R_s2_t_ADC_l[i]->GetFunction("pol1")->GetParameter(0);
      s2_k[i] = hprof_R_s2_t_ADC_l[i]->GetFunction("pol1")->GetParameter(1);

      s2_MIP[i] = (TMath::Power(s2_k[i],2))/TMath::Power(s2_offset[i],2);
      
      cout << "For S2 " << i << ": s2_offset = " << s2_offset[i] << ", s2_k = " << s2_k[i] << endl;
    }
    else{
      s2_offset[i] =  0;
      s2_k[i] = 0;
      s2_MIP[i] = 0;
      cout << "For S2 " << i << ": fit did not work" << endl;
    }
      
    
    // perform fit 
            
    c4[i]->cd(5);
    
    h_R_s2_t_ADC_r[i] = new TH2F(Form("h_R_s2_t_ADC_r_%i",i),Form("s2 right time vs TW effect, %i",i),50,0,0.16,20,-0.015e3,0.015e3);

    T->Draw(Form("((R.s2.rt[%i]) - %.9g):((1/TMath::Sqrt(R.s2.ra_p[%i])))>>h_R_s2_t_ADC_r_%i",i,s2_mean_r[i],i,i),Form("R.tr.n==1 && (R.ps.asum_c+0.9**R.sh.asum_c)>800 && R.ps.asum_c>350 && R.cer.asum_c>200 && R.s2.nthit==1 && R.s2.ra_c[%i]>100",i),"colz");

    
    c4[i]->cd(6);
    
    hprof_R_s2_t_ADC_r[i] = new TProfile(Form("hprof_R_s2_t_ADC_r_%i",i),Form("s2 right time vs TW effect, %i",i),50,0,0.16,-0.015e3,0.015e3);

    T->Draw(Form("((R.s2.rt[%i]) - %.9g):((1/TMath::Sqrt(R.s2.ra_p[%i])))>>hprof_R_s2_t_ADC_r_%i",i,s2_mean_r[i],i,i),Form("R.tr.n==1 && (R.ps.asum_c+0.9**R.sh.asum_c)>800 && R.ps.asum_c>350 && R.cer.asum_c>200 && R.s2.nthit==1 && R.s2.ra_c[%i]>100",i),"colz");    
    
    
  }
  
  //  = new TH1F("h_R_s0_t","s0 left + right time",300,9.4e-7,11e-7);





  // save canvases to pdf

  c1->Print(Form("time_walk/plots/R_s0_tw_%i.pdf(",runno));
  c2->Print(Form("time_walk/plots/R_s0_tw_%i.pdf)",runno));


  
  for(int i=0;i<NS2Pad;i++){
    
    if(i==0){
      c3[i]->Print(Form("time_walk/plots/R_s2_tw_%i.pdf(",runno));
      c4[i]->Print(Form("time_walk/plots/R_s2_tw_%i.pdf",runno));
    }
    else if(i==NS2Pad-1){
      c3[i]->Print(Form("time_walk/plots/R_s2_tw_%i.pdf",runno));
      c4[i]->Print(Form("time_walk/plots/R_s2_tw_%i.pdf)",runno));
    }
    else{
      c3[i]->Print(Form("time_walk/plots/R_s2_tw_%i.pdf",runno));
      c4[i]->Print(Form("time_walk/plots/R_s2_tw_%i.pdf",runno));
    }
  }
 
  
  
  // Print all corrections coeffecients and print to DB

  ofstream oofile(Form("time_walk/DB/RHRS_tw_%i.txt",runno));

  oofile << "S0: "  << endl;
  oofile << "R.s0.MIP = " << s0_MIP << endl;
  oofile << "R.s0.timewalk_params = " << s0_k << endl;
  oofile << endl << endl;

  cout << "S0: "  << endl;
  cout << "R.s0.MIP = " << s0_MIP << endl;
  cout << "R.s0.timewalk_params = " << s0_k << endl;
  cout << endl << endl;


  oofile << "S2: "  << endl;
  oofile << "R.s2.MIP = ";

  cout << "S2: "  << endl;
  cout << "R.s2.MIP = ";

  // calculate mean MIP
  // MIP is arbitary (slope of each channel is not)
  // and analyzer class sets all PMTs to same value
  // mean of all channel MIPs seems reasonable choice for this single MIP value

  Double_t mean_s2_MIP = 0.0; 
  
  for(Int_t i = 0; i<NS2Pad; i++){
    oofile << " " << s2_MIP[i];
    cout << " " << s2_MIP[i];
    mean_s2_MIP += s2_MIP[i];
  }

  mean_s2_MIP /= NS2Pad;

  oofile << endl;
  oofile << "mean_s2_MIP = " << mean_s2_MIP << endl;  
  cout << endl;
  cout << "mean_s2_MIP = " << mean_s2_MIP << endl;
  
  cout << endl << endl;
  oofile << endl << endl;



  
  oofile << "S2: "  << endl;
  oofile << "R.s2.timewalk_params = ";

  cout << "S2: "  << endl;
  cout << "R.s2.timewalk_params = ";

  for(Int_t i = 0; i<NS2Pad; i++){
    oofile << " " << s2_k[i];
    cout << " " << s2_k[i];
  }

  cout << endl << endl;
  oofile << endl << endl;
  oofile.close();
 
}
