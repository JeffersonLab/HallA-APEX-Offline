/*
*************************************************************
7/1/21 John Williamson

Performs timewalk corrections for LHRS.

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


void L_tw_correct(Int_t runno, Bool_t Coinc = true){

  TChain* T = Load_more_rootfiles(runno);

  // set levels of MIP for s0 and s2 (this must match definition in DB)

  const Double_t s0_MIP = 500.0;
  //  const Double_t s2_MIP = 200.0;

  
  const Int_t NS2Pad = 16;
  const Int_t MaxS2Hit = 15;

  Double_t s0_L = 0.0;
  Double_t s0_R = 0.0;
  
  Double_t s2_L[NS2Pad] = {0.0};
  Double_t s2_R[NS2Pad] = {0.0};


  // set-up branches to be read from TTree


  //Define Variables
  Double_t L_tr_n,L_cer_asum_c,L_ps_e,L_sh_e;
  Double_t L_tr_p[100],L_s0_trx[100],L_s2_try[100];
  Double_t L_s0_lt[10],L_s0_rt[10];
  Double_t L_s0_la_p[10],L_s0_ra_p[10];
  Double_t L_s0_nthit;
  Double_t L_s2_lt[16],L_s2_rt[16];
  Double_t L_s2_la_p[16],L_s2_ra_p[16];
  Double_t L_s2_la_c[16],L_s2_ra_c[16];
  Double_t L_s2_nthit;
  Double_t L_s2_t_pads[16];
  Double_t L_s2_trx[100];
  Double_t L_s2_trdx[100];
  Double_t evtypebits;
  Double_t L_s0_trpath[100],L_s2_trpath[100];


  Double_t R_s2_la_c[16],R_s2_ra_c[16];
  Double_t R_s2_t_pads[16];
  Double_t R_tr_n,R_cer_asum_c,R_ps_e,R_sh_e;
  Double_t R_ps_asum_c, R_sh_asum_c;
  Double_t R_s2_nthit;

  // TW variables
  Double_t L_TW_l = 0.0;
  Double_t L_TW_r = 0.0;
  Double_t L_TW_b = 0.0;
  Double_t L_TW_p = 0.0;

  Double_t R_TW_l = 0.0;
  Double_t R_TW_r = 0.0;
  Double_t R_TW_b = 0.0;
  Double_t R_TW_p = 0.0;

  
  
  //trigger
  Double_t Trig_type;

  Double_t Trig_time_T2[MaxS2Hit];
  Double_t Trig_No_T2;
  Double_t Trig_time_T5[MaxS2Hit];
  Double_t Trig_No_T5;


  // record which paddle is hit for an entry
  Int_t LHRS_pad = -1;
  Int_t RHRS_pad = -1;

  
  //Define Branch Status/Addresses
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("L.tr.n",1);
  T->SetBranchStatus("L.tr.p",1);
  T->SetBranchStatus("L.cer.asum_c",1);
  T->SetBranchStatus("L.prl1.e",1);
  T->SetBranchStatus("L.prl2.e",1);
  T->SetBranchStatus("L.s0.lt",1);
  T->SetBranchStatus("L.s0.rt",1);
  T->SetBranchStatus("L.s0.trx",1);
  T->SetBranchStatus("L.s0.la_p",1);
  T->SetBranchStatus("L.s0.ra_p",1);
  T->SetBranchStatus("L.s0.nthit",1);
  T->SetBranchStatus("L.s2.lt",1);
  T->SetBranchStatus("L.s2.rt",1);
  T->SetBranchStatus("L.s2.la_p",1);
  T->SetBranchStatus("L.s2.ra_p",1);
  T->SetBranchStatus("L.s2.la_c",1);
  T->SetBranchStatus("L.s2.ra_c",1);
  T->SetBranchStatus("L.s2.try",1);
  T->SetBranchStatus("L.s2.trx",1);
  T->SetBranchStatus("L.s2.trdx",1);
  T->SetBranchStatus("L.s2.nthit",1);
  T->SetBranchStatus("L.s2.t_pads",1);
  T->SetBranchStatus("L.s0.trpath",1);
  T->SetBranchStatus("L.s2.trpath",1);
  T->SetBranchStatus("DL.evtypebits",1);
  T->SetBranchStatus("DL.evtype",1);
  T->SetBranchStatus("DR.rrawt2",1); 
  T->SetBranchStatus("Ndata.DR.rrawt2",1);
  T->SetBranchStatus("DR.rrawt5",1);
  T->SetBranchStatus("Ndata.DR.rrawt5",1);

  T->SetBranchAddress("L.tr.n",&L_tr_n);
  T->SetBranchAddress("L.tr.p",L_tr_p);
  T->SetBranchAddress("L.cer.asum_c",&L_cer_asum_c);
  T->SetBranchAddress("L.prl1.e",&L_ps_e);
  T->SetBranchAddress("L.prl2.e",&L_sh_e);
  T->SetBranchAddress("L.s0.lt",L_s0_lt);
  T->SetBranchAddress("L.s0.rt",L_s0_rt);
  T->SetBranchAddress("L.s0.trx",L_s0_trx);
  T->SetBranchAddress("L.s0.la_p",L_s0_la_p);
  T->SetBranchAddress("L.s0.ra_p",L_s0_ra_p);
  T->SetBranchAddress("L.s0.nthit",&L_s0_nthit);
  T->SetBranchAddress("L.s2.lt",L_s2_lt);
  T->SetBranchAddress("L.s2.rt",L_s2_rt);
  T->SetBranchAddress("L.s2.la_p",L_s2_la_p);
  T->SetBranchAddress("L.s2.ra_p",L_s2_ra_p);
  T->SetBranchAddress("L.s2.la_c",L_s2_la_c);
  T->SetBranchAddress("L.s2.ra_c",L_s2_ra_c);
  T->SetBranchAddress("L.s2.try",L_s2_try);
  T->SetBranchAddress("L.s2.trx",L_s2_trx);
  T->SetBranchAddress("L.s2.trdx",L_s2_trdx);
  T->SetBranchAddress("L.s2.nthit",&L_s2_nthit);
  T->SetBranchAddress("L.s2.t_pads",L_s2_t_pads);
  T->SetBranchAddress("DL.evtypebits",&evtypebits);
  T->SetBranchAddress("L.s0.trpath",L_s0_trpath);
  T->SetBranchAddress("L.s2.trpath",L_s2_trpath);
  T->SetBranchAddress("DL.evtype",&Trig_type);
  T->SetBranchAddress("DR.rrawt2", Trig_time_T2); 
  T->SetBranchAddress("Ndata.DR.rrawt2",&Trig_No_T2);
  T->SetBranchAddress("DR.rrawt5",Trig_time_T5);
  T->SetBranchAddress("Ndata.DR.rrawt5",&Trig_No_T5);

  T->SetBranchAddress("R.tr.n",&R_tr_n);
  T->SetBranchAddress("R.s2.la_c",R_s2_la_c);
  T->SetBranchAddress("R.s2.ra_c",R_s2_ra_c);
  T->SetBranchAddress("R.s2.t_pads",R_s2_t_pads);
  T->SetBranchAddress("R.ps.asum_c",&R_ps_asum_c);
  T->SetBranchAddress("R.sh.asum_c",&R_sh_asum_c);
  T->SetBranchAddress("R.cer.asum_c",&R_cer_asum_c);
  T->SetBranchAddress("R.ps.e",&R_ps_e);
  T->SetBranchAddress("R.sh.e",&R_sh_e);
  T->SetBranchAddress("R.s2.nthit",&R_s2_nthit);
  
  Int_t nentries = T->GetEntries();

  cout<<"Total Number of Events = "<<nentries<<endl;


  // set up histograms
  
  // s0
  
  TH1F* h_L_s0_t = new TH1F("h_L_s0_t","s0 left + right time",70,5.38e3,5.45e3);

  TH1F* h_L_s0_ADC = new TH1F("h_L_s0_ADC","s0 left + right TW effect",300,0,0.16);

  TH2F* h_L_s0_t_ADC = new TH2F("h_L_s0_t_ADC","s0 left + right time vs TW effect",50,0,0.16,60,-35,35);

  TProfile* hprof_L_s0_t_ADC = new TProfile("hprof_L_s0_t_ADC","s0 left + right time vs TW effect",50,0,0.16,-35,35);


  // S2

  
  TCut LCut = "L.tr.n==1 && abs(L.s2.trdx)<0.07 && (L.prl1.e/(L.gold.p*1000))>0.3 && ((L.prl1.e+L.prl2.e)/(L.gold.p*1000))>0.625 && ((L.prl1.e+L.prl2.e)/(L.gold.p*1000))<1.11 &&  L.cer.asum_c >650";
    
  TCut L_s2_Cut[NS2Pad] ;
  TCut L_s2_lCut[NS2Pad] ;
  TCut L_s2_rCut[NS2Pad] ;

  
  TH1F* h_L_s2_t[NS2Pad];
  TH1F* h_L_s2_ADC[NS2Pad];
  TF1* f_cb_ADC[NS2Pad];
  TF1* f_line_ADC[NS2Pad];
  TF1* s2_ADC_intersection[NS2Pad];
  TLine* s2_ADC_l_cut[NS2Pad];
  TLine* s2_ADC_h_cut[NS2Pad];
  TLine* s2_ADC_l_cut_t[NS2Pad];
  TLine* s2_ADC_h_cut_t[NS2Pad];
  Double_t s2_ADC_l_val[NS2Pad];
  Double_t s2_ADC_h_val[NS2Pad];
  Double_t s2_offset[NS2Pad];
  Double_t s2_k[NS2Pad];
  Double_t s2_MIP[NS2Pad];
  TH2F* h_L_s2_t_ADC[NS2Pad]; 
  TProfile* hprof_L_s2_t_ADC[NS2Pad];
  TH2F* h_L_s2_t_ADCl[NS2Pad]; 
  TProfile* hprof_L_s2_t_ADCl[NS2Pad];
  TH2F* h_L_s2_t_ADCr[NS2Pad]; 
  TProfile* hprof_L_s2_t_ADCr[NS2Pad];
  

  TH1F* h_L_s2_t_l[NS2Pad];
  TH1F* h_L_s2_ADC_l[NS2Pad];
  TH1F* h_L_s2_ADC_l_check[NS2Pad];
  TF1* f_cb_ADC_l[NS2Pad];
  TF1* f_line_ADC_l[NS2Pad];
  TF1* s2_l_ADC_intersection[NS2Pad];
  TLine* s2_l_ADC_l_cut[NS2Pad];
  TLine* s2_l_ADC_h_cut[NS2Pad];
  TLine* s2_l_ADC_l_cut_t[NS2Pad];
  TLine* s2_l_ADC_h_cut_t[NS2Pad];
  Double_t s2_l_ADC_l_val[NS2Pad];
  Double_t s2_l_ADC_h_val[NS2Pad];
  Double_t s2_offset_l[NS2Pad];
  Double_t s2_k_l[NS2Pad];
  Double_t s2_MIP_l[NS2Pad];
  TH2F* h_L_s2_tl_ADC[NS2Pad];
  TProfile* hprof_L_s2_tl_ADC[NS2Pad];
  TH2F* h_L_s2_tl_ADCl[NS2Pad];
  TProfile* hprof_L_s2_tl_ADCl[NS2Pad];
  TH2F* h_L_s2_tl_ADCr[NS2Pad];
  TProfile* hprof_L_s2_tl_ADCr[NS2Pad];
  Double_t tw_stop[NS2Pad];

  
  
  TH1F* h_L_s2_t_r[NS2Pad];
  TH1F* h_L_s2_ADC_r[NS2Pad];
  TF1* f_cb_ADC_r[NS2Pad];
  // TH2F* h_L_s2_t_ADC_r[NS2Pad];
  // TProfile* hprof_L_s2_t_ADC_r[NS2Pad];
  TH2F* h_L_s2_tr_ADC[NS2Pad];
  TProfile* hprof_L_s2_tr_ADC[NS2Pad];
  TH2F* h_L_s2_tr_ADCl[NS2Pad];
  TProfile* hprof_L_s2_tr_ADCl[NS2Pad];
  TH2F* h_L_s2_tr_ADCr[NS2Pad];
  TProfile* hprof_L_s2_tr_ADCr[NS2Pad];
  

  TCanvas* c3[NS2Pad];
  TCanvas* c4[NS2Pad];

  Double_t s2_mean[NS2Pad];
  Double_t s2_mean_l[NS2Pad];
  Double_t s2_mean_r[NS2Pad];


  // variable limits
  
  const Double_t fTdc2T = 0.5e-9;

  
  // cuts on difference of timewalk effect between L-PMT and R-PMT

  Double_t TW_lim = 0.007;



  // time limits change depdending on if run was coincidence or single-arm
  // single-arm is self-timing peak (I think)
  
  // time limits (difference from mean)
  Double_t tdiff_l = -10;
  Double_t tdiff_h = 10;
  Int_t tdiff_nbins = (tdiff_h - tdiff_l);
  
  // time limits (no correction)

  
  Double_t time_l = 0.0; // left + right pmts
  Double_t time_h = 0.0; // left + right pmts
  
  // left time
  Double_t timel_l = 0.0; // left pmts
  Double_t timel_h = 0.0; // left pmts


  // right time 
  Double_t timer_l = 0.0; // right pmts
  Double_t timer_h = 0.0; // right pmts

  if(Coinc){

    // coinc peak (between LHRS and RHRS)

    time_l = 5.56e3;
    time_h = 5.66e3; 
  
    timel_l = 2.80e3; 
    timel_h = 2.86e3; 

    timer_l = 2.74e3; 
    timer_h = 2.80e3;
    
  }
  else{

    // for single arm LHRS run
    
    time_l = 5.5e3;
    time_h = 5.56e3; 

    timel_l = 2.76e3; 
    timel_h = 2.82e3;

    timer_l = 2.72e3; 
    timer_h = 2.75e3;     

  }


  
  Int_t time_nbins = (time_h - time_l);
  Int_t timel_nbins = (timel_h - timel_l);
  Int_t timer_nbins = (timer_h - timer_l);
  

  
  
  // single pmt TW lims
  Double_t tw_s_low = 0.0;
  Double_t tw_s_high = 0.08;
  Int_t tw_s_nbins = 50;

  Double_t tw_l_low = tw_s_low;
  Double_t tw_l_high = tw_s_high;
  Int_t tw_l_nbins = tw_s_nbins;

  Double_t tw_r_low = tw_s_low;
  Double_t tw_r_high = tw_s_high;
  Int_t tw_r_nbins = tw_s_nbins;

  
  // both arm pmt TW lims
  Double_t tw_b_low = 0.05;
  Double_t tw_b_high = 0.16;
  Int_t tw_b_nbins = 69;


  
  // parameters for S2 and trigger timing cuts

  Double_t T2LowCut = 1683;
  Double_t T2HighCut = 1710;

  Double_t T5LowCut = 1625;
  Double_t T5HighCut = 1634;

  Bool_t T2Pass = false;
  Bool_t T5Pass = false;
		     		     
  Bool_t CoincPass = false;

  for(Int_t i = 0; i<NS2Pad; i++){



    // define cuts
    L_s2_Cut[i] = Form("L.s2.nthit==1 && L.s2.t_pads==%i && L.s2.la_c[%i]>100 && L.s2.ra_c[%i]>100",i,i,i);

    L_s2_lCut[i] = Form("L.s2.nthit==1 && L.s2.t_pads==%i && L.s2.la_c[%i]>100",i,i);

    L_s2_rCut[i] = Form("L.s2.nthit==1 && L.s2.t_pads==%i && L.s2.ra_c[%i]>100",i,i);
    

    h_L_s2_t[i] = new TH1F(Form("h_L_s2_t_%i",i),Form("s2 left + right time, %i",i),time_nbins,time_l,time_h);

    h_L_s2_t_l[i] = new TH1F(Form("h_L_s2_t_l_%i",i),Form("s2 left time, %i",i),timel_nbins,timel_l,timel_h);
 
    h_L_s2_t_r[i] = new TH1F(Form("h_L_s2_t_r_%i",i),Form("s2 right time, %i",i),timer_nbins,timer_l,timer_h);


    h_L_s2_ADC[i] = new TH1F(Form("h_L_s2_ADC_%i",i),Form("s2 left + right TW effect, %i",i),tw_b_nbins,tw_b_low,tw_b_high);

    h_L_s2_ADC_l[i] = new TH1F(Form("h_L_s2_ADC_l_%i",i),Form("s2 left TW effect, %i",i),tw_l_nbins,tw_l_low,tw_l_high);

    h_L_s2_ADC_r[i] = new TH1F(Form("h_L_s2_ADC_r_%i",i),Form("s2 right TW effect, %i",i),tw_r_nbins,tw_r_low,tw_r_high);


    // (left + right time)
    
    h_L_s2_t_ADC[i] = new TH2F(Form("h_L_s2_t_ADC_%i",i),Form("s2 left + right time vs (L+R) TW effect, %i",i),tw_b_nbins,tw_b_low,tw_b_high,tdiff_nbins,tdiff_l,tdiff_h);

    hprof_L_s2_t_ADC[i] = new TProfile(Form("hprof_L_s2_t_ADC_%i",i),Form("s2 left + right time vs (L+R) TW effect, %i",i),tw_b_nbins,tw_b_low,tw_b_high,tdiff_l,tdiff_h);

    h_L_s2_t_ADCl[i] = new TH2F(Form("h_L_s2_t_ADCl_%i",i),Form("s2 left + right time vs L TW effect, %i",i),tw_l_nbins,tw_l_low,tw_l_high,tdiff_nbins,tdiff_l,tdiff_h);

    hprof_L_s2_t_ADCl[i] = new TProfile(Form("hprof_L_s2_t_ADCl_%i",i),Form("s2 left + right time vs L TW effect, %i",i),tw_l_nbins,tw_l_low,tw_l_high,tdiff_l,tdiff_h);

    h_L_s2_t_ADCr[i] = new TH2F(Form("h_L_s2_t_ADCr_%i",i),Form("s2 left + right time vs R TW effect, %i",i),tw_r_nbins,tw_r_low,tw_r_high,tdiff_nbins,tdiff_l,tdiff_h);

    hprof_L_s2_t_ADCr[i] = new TProfile(Form("hprof_L_s2_t_ADCr_%i",i),Form("s2 left + right time vs R TW effect, %i",i),tw_r_nbins,tw_r_low,tw_r_high,tdiff_l,tdiff_h);


    // (left time)

    h_L_s2_tl_ADC[i] = new TH2F(Form("h_L_s2_tl_ADC_%i",i),Form("s2 left time vs (L+R) TW effect, %i",i),tw_b_nbins,tw_b_low,tw_b_high,tdiff_nbins,tdiff_l,tdiff_h);
	
    hprof_L_s2_tl_ADC[i] = new TProfile(Form("hprof_L_s2_tl_ADC_%i",i),Form("s2 left time vs (L+R) TW effect, %i",i),tw_b_nbins,tw_b_low,tw_b_high,tdiff_l,tdiff_h);
    
    h_L_s2_tl_ADCl[i] = new TH2F(Form("h_L_s2_tl_ADCl_%i",i),Form("s2 left time vs L TW effect, %i",i),tw_l_nbins,tw_l_low,tw_l_high,tdiff_nbins,tdiff_l,tdiff_h);

    hprof_L_s2_tl_ADCl[i] = new TProfile(Form("hprof_L_s2_tl_ADCl_%i",i),Form("s2 left time vs L TW effect, %i",i),tw_l_nbins,tw_l_low,tw_l_high,tdiff_l,tdiff_h);

    h_L_s2_tl_ADCr[i] = new TH2F(Form("h_L_s2_tl_ADCr_%i",i),Form("s2 left time vs R TW effect, %i",i),tw_r_nbins,tw_r_low,tw_r_high,tdiff_nbins,tdiff_l,tdiff_h);

    hprof_L_s2_tl_ADCr[i] = new TProfile(Form("hprof_L_s2_tl_ADCr_%i",i),Form("s2 left time vs R TW effect, %i",i),tw_r_nbins,tw_r_low,tw_r_high,tdiff_l,tdiff_h);


    // (right time)

    h_L_s2_tr_ADC[i] = new TH2F(Form("h_L_s2_tr_ADC_%i",i),Form("s2 right time vs (L+R) TW effect, %i",i),tw_b_nbins,tw_b_low,tw_b_high,tdiff_nbins,tdiff_l,tdiff_h);

    hprof_L_s2_tr_ADC[i] = new TProfile(Form("hprof_L_s2_tr_ADC_%i",i),Form("s2 right time vs (L+R) TW effect, %i",i),tw_b_nbins,tw_b_low,tw_b_high,tdiff_l,tdiff_h);

    h_L_s2_tr_ADCl[i] = new TH2F(Form("h_L_s2_tr_ADCl_%i",i),Form("s2 right time vs L TW effect, %i",i),tw_l_nbins,tw_l_low,tw_l_high,tdiff_nbins,tdiff_l,tdiff_h);

    hprof_L_s2_tr_ADCl[i] = new TProfile(Form("hprof_L_s2_tr_ADCl_%i",i),Form("s2 right time vs L TW effect, %i",i),tw_l_nbins,tw_l_low,tw_l_high,tdiff_l,tdiff_h);
    
    h_L_s2_tr_ADCr[i] = new TH2F(Form("h_L_s2_tr_ADCr_%i",i),Form("s2 right time vs R TW effect, %i",i),tw_r_nbins,tw_r_low,tw_r_high,tdiff_nbins,tdiff_l,tdiff_h);

    hprof_L_s2_tr_ADCr[i] = new TProfile(Form("hprof_L_s2_tr_ADCr_%i",i),Form("s2 right time vs R TW effect, %i",i),tw_r_nbins,tw_r_low,tw_r_high,tdiff_l,tdiff_h);

    
    f_cb_ADC[i] = new TF1(Form("f_cb_ADC_%i",i),"crystalball",tw_b_low,tw_b_high);    
    
    f_cb_ADC_l[i] = new TF1(Form("f_cb_ADC_l_%i",i),"crystalball",tw_l_low,tw_l_high);

    f_line_ADC[i] = new TF1(Form("f_line_ADC_%i",i),"pol1",tw_b_low,tw_b_high);

    f_line_ADC_l[i] = new TF1(Form("f_line_ADC_l_%i",i),"pol1",tw_l_low,tw_l_high);
    
  
  }
  
  // marker

  
  for(Int_t i=0;i<nentries;i++){
    
    T2Pass = false;
    T5Pass = false;
    CoincPass = false;

    LHRS_pad = -1;
    RHRS_pad = -1;
    
    
    if(i%100000==0) cout << " events processed = " << i << endl;
    T->GetEntry(i);
    
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
    
    if((T2Pass && T5Pass && Trig_type == 6) || !Coinc){
      CoincPass = true;
    }


    
    
    if(L_tr_n==1 && L_cer_asum_c>1500 && (L_ps_e+L_sh_e)/(1000.*L_tr_p[0])>0.8 && L_s0_nthit==1 && CoincPass){
      h_L_s0_t->Fill(L_s0_lt[0] + L_s0_rt[0]);
      h_L_s0_ADC->Fill((1/TMath::Sqrt(L_s0_la_p[0])) + (1/TMath::Sqrt(L_s0_ra_p[0])));
      //      h_L_s0_t_ADC->Fill();

    }


    
    if(L_s2_nthit ==1 &&  L_tr_n==1 && L_cer_asum_c>1500 && (L_ps_e+L_sh_e)/(1000.*L_tr_p[0])>0.8  && R_tr_n==1 && R_cer_asum_c>200 &&  R_s2_nthit==1 && (R_ps_asum_c + 0.9*R_sh_asum_c)>800 && R_ps_asum_c>350 && CoincPass){
    
      for(int j=0;j<NS2Pad;j++){                                                                              
	if ( L_s2_lt[j]>0 && L_s2_rt[j]>0 && abs(L_s2_trdx[0])<0.07 && L_s2_la_c[j]>160 && L_s2_ra_c[j]>160 && L_s2_t_pads[0]==j){    
	  LHRS_pad = j;
	}    
	
	
	
	if (R_s2_la_c[j]>160 && R_s2_ra_c[j]>160 && R_s2_t_pads[0]==j){    
	  RHRS_pad = j;
	}	
      }
    }


    // condition on LHRS pad being found for single-arm run
    // or both LHRS and RHRS pad found for coincidence run
    if((Coinc && (LHRS_pad == -1 || RHRS_pad == -1)) || (!Coinc && (LHRS_pad == -1)) ){
      continue;
    }
    

    L_TW_l = 1/TMath::Sqrt(L_s2_la_c[LHRS_pad]);
    L_TW_r = 1/TMath::Sqrt(L_s2_ra_c[LHRS_pad]);
    L_TW_b = L_TW_l - L_TW_r;
    L_TW_p = L_TW_l + L_TW_r;
	
    if(Coinc){
      R_TW_l = 1/TMath::Sqrt(R_s2_la_c[RHRS_pad]);
      R_TW_r = 1/TMath::Sqrt(R_s2_ra_c[RHRS_pad]);
      R_TW_b = R_TW_l - R_TW_r;
      R_TW_p = R_TW_l + R_TW_r;
    }


    if(Coinc && !(abs(L_TW_b) < TW_lim && abs(R_TW_b) < TW_lim) ){
      // cut on difference on left and right PMTs being too large (for LHRS and RHRS)
      continue;
    }
    
    h_L_s2_t[LHRS_pad]->Fill(L_s2_lt[LHRS_pad] + L_s2_rt[LHRS_pad]);
    h_L_s2_ADC[LHRS_pad]->Fill((1/TMath::Sqrt(L_s2_la_p[LHRS_pad])) + (1/TMath::Sqrt(L_s2_ra_p[LHRS_pad])));

	  
    h_L_s2_t_l[LHRS_pad]->Fill(L_s2_lt[LHRS_pad]);
    h_L_s2_ADC_l[LHRS_pad]->Fill((1/TMath::Sqrt(L_s2_la_p[LHRS_pad])));


    h_L_s2_t_r[LHRS_pad]->Fill(L_s2_rt[LHRS_pad]);
    h_L_s2_ADC_r[LHRS_pad]->Fill((1/TMath::Sqrt(L_s2_ra_p[LHRS_pad])));


  }
  


  
  

  TCanvas *c1 = new TCanvas("c1","s0 timing");
  c1->Divide(2,1);
  c1->cd(1);
  
  h_L_s0_t->Draw();
  
  
  Double_t s0_t_mean = h_L_s0_t->GetMean();

  cout << "s0_t_mean = " << s0_t_mean << endl;

  c1->cd(2); 
  h_L_s0_ADC->Draw();


  
  for(Int_t i = 0; i<NS2Pad; i++){

    // define cut

    c3[i] = new TCanvas(Form("c3_%i",i),Form("S2 timing, %i",i));
    c3[i]->Divide(2,3);

    c3[i]->cd(1);

    h_L_s2_t[i]->Draw();
    
    s2_mean[i] = h_L_s2_t[i]->GetMean();

    
    c3[i]->cd(2);

    h_L_s2_ADC[i]->Draw();

    Double_t max_ADC = h_L_s2_ADC[i]->GetBinCenter(h_L_s2_ADC[i]->GetMaximumBin());
    Int_t max_entriess = h_L_s2_ADC[i]->GetBinContent(h_L_s2_ADC[i]->GetMaximumBin());

    f_cb_ADC[i]->SetParameters(max_entriess,max_ADC,0.0048,0.736,2.7e7);

    h_L_s2_ADC[i]->Fit(Form("f_cb_ADC_%i",i),"QR","",0,0.16);

    Double_t s2_lr_mean = f_cb_ADC[i]->GetParameter(1);
    Double_t s2_lr_sigma = f_cb_ADC[i]->GetParameter(2);

    // check 2 sigma above and below mean and check which has more bin entries

    Double_t s2_lr_lower = s2_lr_mean-2.5*s2_lr_sigma;
    Double_t s2_lr_higher = s2_lr_mean+2.5*s2_lr_sigma;
    
    Int_t l_check = h_L_s2_ADC[i]->GetBinContent(h_L_s2_ADC[i]->FindBin(s2_lr_lower));
    
    Int_t h_check = h_L_s2_ADC[i]->GetBinContent(h_L_s2_ADC[i]->FindBin(s2_lr_higher));

    Double_t s2_ADC_cutoff = (l_check>h_check) ? l_check : h_check;

    cout << "S2: L+R check " << endl;
    cout << "s2_lr_mean = " << s2_lr_mean << endl;
    cout << "s2_lr_sigma = " << s2_lr_sigma << endl;
    cout << "l_check = " << l_check << endl;
    cout << "h_check = " << h_check << endl;
    cout << "s2_ADC_cutoff = " << s2_ADC_cutoff << endl;
    
    

    f_line_ADC[i]->SetParameters(s2_ADC_cutoff,0);

    f_cb = f_cb_ADC[i];
    f_line = f_line_ADC[i];
    
    s2_ADC_intersection[i] = new TF1(Form("s2_ADC_intersection_%i",i),Form("abs(f_line_ADC_%i - f_cb_ADC_%i)",i,i),s2_lr_lower,s2_lr_higher);


    s2_ADC_intersection[i]->SetLineColor(kOrange);
    s2_ADC_intersection[i]->Draw("same");
    
    f_line_ADC[i]->SetLineColor(kGreen+2);
    f_line_ADC[i]->Draw("same");
    
    s2_ADC_l_val[i] = s2_ADC_intersection[i]->GetMinimumX(s2_lr_lower,s2_lr_mean);
    s2_ADC_h_val[i] = s2_ADC_intersection[i]->GetMinimumX(s2_lr_mean,s2_lr_higher);

    s2_ADC_l_cut[i] = new TLine(s2_ADC_l_val[i],0,s2_ADC_l_val[i],s2_ADC_cutoff);
    s2_ADC_h_cut[i] = new TLine(s2_ADC_h_val[i],0,s2_ADC_h_val[i],s2_ADC_cutoff);

    s2_ADC_l_cut_t[i] = new TLine(s2_ADC_l_val[i],-0.05e-7,s2_ADC_l_val[i],0.05e-7);
    s2_ADC_h_cut_t[i] = new TLine(s2_ADC_h_val[i],-0.05e-7,s2_ADC_h_val[i],0.05e-7);
 

    s2_ADC_l_cut[i]->SetLineColor(kMagenta);
    s2_ADC_h_cut[i]->SetLineColor(kMagenta);

    s2_ADC_l_cut[i]->Draw("same");
    s2_ADC_h_cut[i]->Draw("same");

    cout << "S2: L+R " << endl;
    cout << " Low cut = " << s2_lr_lower << ", high cut = " << s2_ADC_h_val[i] << endl;
    cout << " Low cut = " << s2_ADC_l_val[i] << ", high cut = " << s2_ADC_h_val[i] << endl;
    cout << endl << endl;

               
    c3[i]->cd(3);

    h_L_s2_t_l[i]->Draw();
    s2_mean_l[i] = h_L_s2_t_l[i]->GetMean();


    c3[i]->cd(4);

    h_L_s2_ADC_l[i]->Draw();
    Double_t max_ADC_l = h_L_s2_ADC_l[i]->GetBinCenter(h_L_s2_ADC_l[i]->GetMaximumBin());
    
    Int_t max_entriess_l = h_L_s2_ADC_l[i]->GetBinContent(h_L_s2_ADC_l[i]->GetMaximumBin());


    cout << "For " << i << " max ADC = " << max_ADC_l << ", with " << max_entriess_l << endl;
      
    f_cb_ADC_l[i]->SetParameters(max_entriess_l,max_ADC_l,0.0048,0.736,2.7e7);

    h_L_s2_ADC_l[i]->Fit(Form("f_cb_ADC_l_%i",i),"QR","",0,0.16);

    Double_t s2_l_mean = f_cb_ADC_l[i]->GetParameter(1);
    Double_t s2_l_sigma = f_cb_ADC_l[i]->GetParameter(2);

    // check 2 sigma above and below mean and check which has more bin entries

    Double_t s2_l_lower = s2_l_mean-1.5*s2_l_sigma;
    Double_t s2_l_higher = s2_l_mean+1.5*s2_l_sigma;
    
    Int_t l_check_l = h_L_s2_ADC_l[i]->GetBinContent(h_L_s2_ADC_l[i]->FindBin(s2_l_lower));
    
    Int_t h_check_l = h_L_s2_ADC_l[i]->GetBinContent(h_L_s2_ADC_l[i]->FindBin(s2_l_higher));

    cout << "s2_l_lower = " << s2_l_lower << endl;
    cout << "s2_l_higher = " << s2_l_higher << endl;

    
    Double_t s2_l_ADC_cutoff = (l_check_l>h_check_l) ? l_check_l : h_check_l;

    f_line_ADC_l[i]->SetParameters(s2_l_ADC_cutoff,0);


    f_cb = f_cb_ADC_l[i];
    f_line = f_line_ADC_l[i];
    

    s2_l_ADC_intersection[i] = new TF1(Form("s2_l_ADC_intersection_%i",i),Form("abs(f_line_ADC_l_%i - f_cb_ADC_l_%i)",i,i),s2_l_lower,s2_l_higher);


    s2_l_ADC_intersection[i]->SetLineColor(kOrange);
    s2_l_ADC_intersection[i]->Draw("same");
    
    f_line_ADC_l[i]->SetLineColor(kGreen+2);
    f_line_ADC_l[i]->Draw("same");


    s2_l_ADC_l_val[i] = s2_l_ADC_intersection[i]->GetMinimumX(s2_l_lower,s2_l_mean);
    s2_l_ADC_h_val[i] = s2_l_ADC_intersection[i]->GetMinimumX(s2_l_mean,s2_l_higher);

    s2_l_ADC_l_cut[i] = new TLine(s2_l_ADC_l_val[i],0,s2_l_ADC_l_val[i],s2_l_ADC_cutoff);
    s2_l_ADC_h_cut[i] = new TLine(s2_l_ADC_h_val[i],0,s2_l_ADC_h_val[i],s2_l_ADC_cutoff);

    s2_l_ADC_l_cut_t[i] = new TLine(s2_l_ADC_l_val[i],-0.05e-7,s2_l_ADC_l_val[i],0.05e-7);
    s2_l_ADC_h_cut_t[i] = new TLine(s2_l_ADC_h_val[i],-0.05e-7,s2_l_ADC_h_val[i],0.05e-7);
 

    s2_l_ADC_l_cut[i]->SetLineColor(kMagenta);
    s2_l_ADC_h_cut[i]->SetLineColor(kMagenta);

    s2_l_ADC_l_cut[i]->Draw("same");
    s2_l_ADC_h_cut[i]->Draw("same");

    cout << "S2: L " << endl;
    cout << " Low cut = " << s2_l_lower << ", high cut = " << s2_l_ADC_h_val[i] << endl;
    cout << " Low cut = " << s2_l_ADC_l_val[i] << ", high cut = " << s2_l_ADC_h_val[i] << endl;
    cout << endl << endl;    
    
    cout << "h_check = " << h_check << endl;
    cout << "l_check = " << l_check << endl;
    cout << "s2_l_ADC_cutoff = " << s2_l_ADC_cutoff << endl << endl;;

             
    c3[i]->cd(5);

    h_L_s2_t_r[i]->Draw();
    s2_mean_r[i] = h_L_s2_t_r[i]->GetMean();

    c3[i]->cd(6);
	    
    h_L_s2_ADC_r[i]->Draw();


  }




  for(Int_t i=0;i<nentries;i++){

    T2Pass = false;
    T5Pass = false;
    CoincPass = false;

    LHRS_pad = -1;
    RHRS_pad = -1;

    
    if(i%100000==0) cout << " events processed = " << i << endl;
    T->GetEntry(i);
  


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
    
    if((T2Pass && T5Pass && Trig_type == 6) || !Coinc){
      CoincPass = true;
    }

    
    if(L_tr_n==1 && L_cer_asum_c>1500 && (L_ps_e+L_sh_e)/(1000.*L_tr_p[0])>0.8 && L_s0_nthit==1  && CoincPass){
      
      h_L_s0_t_ADC->Fill((1/TMath::Sqrt(L_s0_la_p[0])) + (1/TMath::Sqrt(L_s0_ra_p[0])), L_s0_lt[0] + L_s0_rt[0] - s0_t_mean);
      hprof_L_s0_t_ADC->Fill((1/TMath::Sqrt(L_s0_la_p[0])) + (1/TMath::Sqrt(L_s0_ra_p[0])), L_s0_lt[0] + L_s0_rt[0] - s0_t_mean);

    }



    if(L_s2_nthit ==1 &&  L_tr_n==1 && L_cer_asum_c>1500 && (L_ps_e+L_sh_e)/(1000.*L_tr_p[0])>0.8  && R_tr_n==1 && R_cer_asum_c>200 &&  R_s2_nthit==1 && (R_ps_asum_c + 0.9*R_sh_asum_c)>800 && R_ps_asum_c>350 && CoincPass){
    
      for(int j=0;j<NS2Pad;j++){                                                                              
	if ( L_s2_lt[j]>0 && L_s2_rt[j]>0 && abs(L_s2_trdx[0])<0.07 && L_s2_la_c[j]>160 && L_s2_ra_c[j]>160 && L_s2_t_pads[0]==j){    
	  LHRS_pad = j;
	}    
	
	
	
	if (R_s2_la_c[j]>160 && R_s2_ra_c[j]>160  && R_s2_t_pads[0]==j){    
	  RHRS_pad = j;
	}	
      }
    }
    
    if((Coinc && (LHRS_pad == -1 || RHRS_pad == -1)) || (!Coinc && (LHRS_pad == -1)) ){
      continue;
    }

    

    L_TW_l = 1/TMath::Sqrt(L_s2_la_c[LHRS_pad]);
    L_TW_r = 1/TMath::Sqrt(L_s2_ra_c[LHRS_pad]);
    L_TW_b = L_TW_l - L_TW_r;
    L_TW_p = L_TW_l + L_TW_r;
	
    if(Coinc){
      R_TW_l = 1/TMath::Sqrt(R_s2_la_c[RHRS_pad]);
      R_TW_r = 1/TMath::Sqrt(R_s2_ra_c[RHRS_pad]);
      R_TW_b = R_TW_l - R_TW_r;
      R_TW_p = R_TW_l + R_TW_r;
    }


    if(Coinc && !(abs(L_TW_b) < TW_lim && abs(R_TW_b) < TW_lim) ){
      // cut on difference on left and right PMTs being too large (for LHRS and RHRS)
      continue;
    }

    
    
    Double_t b_time =  L_s2_lt[LHRS_pad] + L_s2_rt[LHRS_pad] - s2_mean[LHRS_pad];
    Double_t l_time =  L_s2_lt[LHRS_pad] - s2_mean_l[LHRS_pad];
    Double_t r_time =  L_s2_rt[LHRS_pad] - s2_mean_r[LHRS_pad];
    
    
    // (left + right time)
    h_L_s2_t_ADC[LHRS_pad]->Fill(L_TW_p,b_time);
    hprof_L_s2_t_ADC[LHRS_pad]->Fill(L_TW_p,b_time);
    
    h_L_s2_t_ADCl[LHRS_pad]->Fill(L_TW_l,b_time);
    hprof_L_s2_t_ADCl[LHRS_pad]->Fill(L_TW_l,b_time);
    
    h_L_s2_t_ADCr[LHRS_pad]->Fill(L_TW_r,b_time);
    hprof_L_s2_t_ADCr[LHRS_pad]->Fill(L_TW_r,b_time);
    
	  

    // (left time)
    h_L_s2_tl_ADC[LHRS_pad]->Fill(L_TW_p,l_time);
    hprof_L_s2_tl_ADC[LHRS_pad]->Fill(L_TW_p,l_time);
	  
    h_L_s2_tl_ADCl[LHRS_pad]->Fill(L_TW_l,l_time);
    hprof_L_s2_tl_ADCl[LHRS_pad]->Fill(L_TW_l,l_time);

    h_L_s2_tl_ADCr[LHRS_pad]->Fill(L_TW_r,l_time);
    hprof_L_s2_tl_ADCr[LHRS_pad]->Fill(L_TW_r,l_time);


    // (right time)
    h_L_s2_tr_ADC[LHRS_pad]->Fill(L_TW_p,r_time);
    hprof_L_s2_tr_ADC[LHRS_pad]->Fill(L_TW_p,r_time);
		  
    h_L_s2_tr_ADCl[LHRS_pad]->Fill(L_TW_l,r_time);
    hprof_L_s2_tr_ADCl[LHRS_pad]->Fill(L_TW_l,r_time);
	  
    h_L_s2_tr_ADCr[LHRS_pad]->Fill(L_TW_r,r_time);
    hprof_L_s2_tr_ADCr[LHRS_pad]->Fill(L_TW_r,r_time);

	  
    
  }



  
  

  
  TCanvas *c2 = new TCanvas("c2","s0 timing vs TW effects");

  c2->Divide(2,1);
  c2->cd(1);
  
  h_L_s0_t_ADC->Draw("colz");

  c2->cd(2);
  hprof_L_s0_t_ADC->Draw();


    
  Double_t s0_offset = 0.0;
  Double_t s0_k = 0.0;
  

  if(!Coinc){
    hprof_L_s0_t_ADC->Fit("pol1","QR","",0,0.16);
  
    s0_offset =  hprof_L_s0_t_ADC->GetFunction("pol1")->GetParameter(0);
    s0_k = hprof_L_s0_t_ADC->GetFunction("pol1")->GetParameter(1);

  // can now derive ADC_MIP for s0 from these parameters

    //    s0_MIP = (4*TMath::Power(s0_k,2))/TMath::Power(s0_offset,2);

    cout << "s0 k = " << s0_k << ", s0 ADC_MIP = " << s0_MIP << endl;
    
  }

    
  for(Int_t i = 0; i<NS2Pad; i++){

    // left + right time
    
    c4[i] = new TCanvas(Form("c4_%i",i),Form("s2 timing vs TW effects, %i",i));

    c4[i]->Divide(6,3);
    c4[i]->cd(1);
    
    h_L_s2_t_ADC[i]->Draw("colz");

      
    c4[i]->cd(2);
    
    hprof_L_s2_t_ADC[i]->Draw();

    s2_ADC_l_cut_t[i]->SetLineColor(kMagenta);
    s2_ADC_h_cut_t[i]->SetLineColor(kMagenta);

    hprof_L_s2_t_ADC[i]->Fit("pol1","QR","",s2_ADC_l_val[i],s2_ADC_h_val[i]);
    
    if(hprof_L_s2_t_ADC[i]->GetFunction("pol1")){
    
      s2_offset[i] =  hprof_L_s2_t_ADC[i]->GetFunction("pol1")->GetParameter(0);
      s2_k[i] = hprof_L_s2_t_ADC[i]->GetFunction("pol1")->GetParameter(1);

      s2_MIP[i] = (TMath::Power(s2_k[i],2))/TMath::Power(s2_offset[i],2);
      
      cout << "For S2 " << i << ": s2_offset = " << s2_offset[i] << ", s2_k = " << s2_k[i] << endl;
    }
    else{
      s2_offset[i] =  0;
      s2_k[i] = 0;
      s2_MIP[i] = 0;
      cout << "For S2 (l+r) " << i << ": fit did not work" << endl;
    }
      

    s2_l_ADC_l_cut_t[i]->Draw("same");
    s2_l_ADC_h_cut_t[i]->Draw("same");


    c4[i]->cd(3);

    h_L_s2_t_ADCl[i]->Draw("colz");

    c4[i]->cd(4);

    hprof_L_s2_t_ADCl[i]->Draw();

    c4[i]->cd(5);

    h_L_s2_t_ADCr[i]->Draw("colz");

    c4[i]->cd(6);

    hprof_L_s2_t_ADCr[i]->Draw();


    // left time

    c4[i]->cd(11);

    h_L_s2_tl_ADCr[i]->Draw("colz");
    
    c4[i]->cd(9);
    
    h_L_s2_tl_ADCl[i]->Draw("colz");


    
    c4[i]->cd(10);
    
    hprof_L_s2_tl_ADCl[i]->Draw();

    cout << "Pre-line " << endl;
    s2_l_ADC_l_cut_t[i]->SetLineColor(kMagenta);
    s2_l_ADC_h_cut_t[i]->SetLineColor(kMagenta);
    cout << "Post-line " << endl;
  
    //    hprof_L_s2_tl_ADCl[i]->Fit("pol1","QR","",s2_l_ADC_l_val[i],s2_l_ADC_h_val[i]);
    hprof_L_s2_tl_ADCl[i]->Fit("pol1","QR","",0.056,0.063);
    //hprof_L_s2_tl_ADCl[i]->Fit("pol1","QR","",0,0.16);

    if(hprof_L_s2_tl_ADCl[i]->GetFunction("pol1")){
    
      s2_offset_l[i] =  hprof_L_s2_tl_ADCl[i]->GetFunction("pol1")->GetParameter(0);
      s2_k_l[i] = hprof_L_s2_tl_ADCl[i]->GetFunction("pol1")->GetParameter(1);

      s2_MIP_l[i] = (TMath::Power(s2_k_l[i],2))/TMath::Power(s2_offset_l[i],2);
      
      cout << "For S2 " << i << ": s2_offset = " << s2_offset_l[i] << ", s2_k = " << s2_k_l[i] << endl;
    }
    else{
      s2_offset_l[i] =  0;
      s2_k_l[i] = 0;
      s2_MIP_l[i] = 0;
      cout << "For S2 " << i << ": fit did not work" << endl;
    }
      

    s2_l_ADC_l_cut_t[i]->Draw("same");
    s2_l_ADC_h_cut_t[i]->Draw("same");
    
    

    c4[i]->cd(7);

    h_L_s2_tl_ADC[i]->Draw("colz");

    c4[i]->cd(8);

    hprof_L_s2_tl_ADC[i]->Draw();

    c4[i]->cd(11);

    h_L_s2_tl_ADCr[i]->Draw("colz");

    c4[i]->cd(12);

    hprof_L_s2_tl_ADCr[i]->Draw();

    
    
    // perform fit 

    // right time
    
    c4[i]->cd(13);
    
    h_L_s2_tr_ADC[i]->Draw("colz");

    
    c4[i]->cd(14);
    
    hprof_L_s2_tr_ADC[i]->Draw();
    
    c4[i]->cd(15);

    h_L_s2_tr_ADCl[i]->Draw("colz");

    c4[i]->cd(16);

    hprof_L_s2_tr_ADCl[i]->Draw();

    c4[i]->cd(17);

    h_L_s2_tr_ADCr[i]->Draw("colz");

    c4[i]->cd(18);

    hprof_L_s2_tr_ADCr[i]->Draw();

    
    
    
  }







  

  // save canvases to pdf


  if(!Coinc){
    c1->Print(Form("time_walk/plots/L_s0_tw_%i.pdf(",runno));
    c2->Print(Form("time_walk/plots/L_s0_tw_%i.pdf)",runno));
  }

  
  for(int i=0;i<NS2Pad;i++){
    
    if(i==0){
      c3[i]->Print(Form("time_walk/plots/L_s2_tw_%i.pdf(",runno));
      c4[i]->Print(Form("time_walk/plots/L_s2_tw_%i.pdf",runno));
    }
    else if(i==NS2Pad-1){
      c3[i]->Print(Form("time_walk/plots/L_s2_tw_%i.pdf",runno));
      c4[i]->Print(Form("time_walk/plots/L_s2_tw_%i.pdf)",runno));
    }
    else{
      c3[i]->Print(Form("time_walk/plots/L_s2_tw_%i.pdf",runno));
      c4[i]->Print(Form("time_walk/plots/L_s2_tw_%i.pdf",runno));
    }
  }
 
  
  // Print all corrections coeffecients and print to DB

  ofstream oofile(Form("time_walk/DB/LHRS_tw_%i.txt",runno));


  if(!Coinc){
    oofile << "S0: "  << endl;
    oofile << "L.s0.MIP = " << s0_MIP << endl;
    oofile << "L.s0.timewalk_params = " << s0_k << endl;
    oofile << endl << endl;

    cout << "S0: "  << endl;
    cout << "L.s0.MIP = " << s0_MIP << endl;
    cout << "L.s0.timewalk_params = " << s0_k << endl;
    cout << endl << endl;
  }


  oofile << "S2: "  << endl;
  oofile << "L.s2.MIP = ";

  cout << "S2: "  << endl;
  cout << "L.s2.MIP = ";

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
  oofile << "L.s2.timewalk_params = ";

  cout << "S2: "  << endl;
  cout << "L.s2.timewalk_params = ";

  for(Int_t i = 0; i<NS2Pad; i++){
    oofile << " " << s2_k[i];
    cout << " " << s2_k[i];
  }

  cout << endl << endl;
  oofile << endl << endl;
  oofile.close();
 


  cout << "Starting point for TW effect for all channels" << endl;
  for(Int_t i = 0; i<NS2Pad; i++){
    cout << i << "  " << tw_stop[i] << endl;
  }

  // S2 csv outfile

  ofstream oofile_csv(Form("time_walk/DB/LHRS_tw_%i_S2.csv",runno));
  
  oofile_csv<<"L_tw_param, L_MIP, L_start"<<endl;
    
  for (Int_t ii=0;ii<16;ii++)
    {
      oofile_csv<<s2_k_l[ii]<<","<<s2_MIP_l[ii]<<","<<tw_stop[ii]<< endl;       
    }


  // S0 csv outfile


  if(!Coinc){
    ofstream oofile_csv_s0(Form("time_walk/DB/LHRS_tw_%i_S0.csv",runno));
  
    oofile_csv_s0<<"L_tw_param, L_MIP"<<endl;
  
    oofile_csv_s0<<s0_k<<","<<s0_MIP<<endl;
  }
  
}

