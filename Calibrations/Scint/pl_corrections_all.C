/*
*************************************************************
15/1/21 John Williamson

Performs pathlength corrections.

Corrections are based on factors which alter pathlength and thus arrival time of particles in s2. Simlar corrections performed for various experiments including :
- (DVCS) M. Dlamini, Measurement of Hard Exclusive Electroproduction of 0 Meson Cross Section in
Hall A of JLab with CEBAF at 12 GeV, JLAB-PHY-19-3095, DOE/OR/23177-4810, 1574097 (2019)
- C. Dutta, \u201cMEASUREMENT OF SINGLE-TARGET SPIN ASYMMETRIES IN THE ELEC-
TROPRODUCTION OF NEGATIVE PIONS IN THE SEMI-INCLUSIVE DEEP INELASTIC
REACTION n(e,é)X ON A TRANSVERSELY POLARIZED 3He TARGET\u201d, University of
Kentucky (2010)

Plot time in S2 against different variables (x,theta,dp) to determine relationship


*************************************************************
 */

#include "file_def.h"
#include "Load_more_rootfiles.C"


void pl_corrections_all(Int_t runno,  TString DB_Lname /* LHRS DB name where corrections are read from*/, TString DB_Rname /* RHRS DB name where corrections are read from*/, TString Name  = "_" /*Name to be added to csv file*/){


  
  TChain* T = Load_more_rootfiles(runno);
  
  

  const Int_t NS2Pad = 16;
  const Int_t MaxS2Hit = 15;
  
  const Double_t fTdc2T = 0.5e-9;      // seconds/channel
  
  // read offset DBs for LHRS and RHRS

  CsvParser csv_L_db(DB_Lname.Data());
  
  std::vector<string> L_ls2_coeff_s = csv_L_db.GetColumn(0);
  std::vector<string> L_rs2_coeff_s = csv_L_db.GetColumn(1);

  CsvParser csv_R_db(DB_Rname.Data());
 

  std::vector<string> R_ls2_coeff_s = csv_R_db.GetColumn(0);
  std::vector<string> R_rs2_coeff_s = csv_R_db.GetColumn(1);


  double L_ls2_coeff[NS2Pad] = {0.0};
  double L_rs2_coeff[NS2Pad] = {0.0};

  double R_ls2_coeff[NS2Pad] = {0.0};
  double R_rs2_coeff[NS2Pad] = {0.0};

    
  if(L_ls2_coeff_s.size() != NS2Pad || L_rs2_coeff_s.size() != NS2Pad || R_ls2_coeff_s.size() != NS2Pad || R_rs2_coeff_s.size() != NS2Pad){
    cout << "s2 offset DB csv files do not have " << NS2Pad << " entries!!" << endl;
  }


  Double_t L_ls2_coeff_mean = 0.0;
  Double_t L_rs2_coeff_mean = 0.0;
  Double_t R_ls2_coeff_mean = 0.0;
  Double_t R_rs2_coeff_mean = 0.0;
  
  for(Int_t i = 0; i<NS2Pad ; i++){
    L_ls2_coeff[i] = stod(L_ls2_coeff_s[i]);
    L_rs2_coeff[i] = stod(L_rs2_coeff_s[i]);
    R_ls2_coeff[i] = stod(R_ls2_coeff_s[i]);
    R_rs2_coeff[i] = stod(R_rs2_coeff_s[i]);
    
    L_ls2_coeff_mean += L_ls2_coeff[i];
    L_rs2_coeff_mean += L_rs2_coeff[i];
    R_ls2_coeff_mean += R_ls2_coeff[i];
    R_rs2_coeff_mean += R_rs2_coeff[i];

  }  
  cout << "DB files read " << endl;

  L_ls2_coeff_mean = L_ls2_coeff_mean/NS2Pad;
  L_rs2_coeff_mean = L_rs2_coeff_mean/NS2Pad;
  R_ls2_coeff_mean = R_ls2_coeff_mean/NS2Pad;
  R_rs2_coeff_mean = R_rs2_coeff_mean/NS2Pad;

  // from corrections get offset from uncorrected coincidence times
  Double_t coin_off = fTdc2T*((L_ls2_coeff_mean + L_rs2_coeff_mean)/2 - (R_ls2_coeff_mean + R_rs2_coeff_mean)/2);

  cout << "coin_off = " << coin_off << endl;
  //  coin_off = 0.0;

  Double_t L_off = fTdc2T*(L_ls2_coeff_mean + L_rs2_coeff_mean)/2;
  Double_t R_off = fTdc2T*(R_ls2_coeff_mean + R_rs2_coeff_mean)/2;
  


  // read TW DBs

  
  TString Tw_1P_Lname = "time_walk/DB/LHRS_diff_tw_1P_4668_S2.csv";
  CsvParser Tw_1P_L_csv(Tw_1P_Lname.Data());

  std::vector<string> L_s2_Tw_1P_slope_s = Tw_1P_L_csv.GetColumn(0);
  std::vector<string> L_s2_Tw_1P_off_s = Tw_1P_L_csv.GetColumn(1);

  
  TString Tw_3P_Lname = "time_walk/DB/LHRS_diff_tw_3P_4668_S2.csv";
  CsvParser Tw_3P_L_csv(Tw_3P_Lname.Data());

  std::vector<string> L_s2_Tw_3P_3p_s = Tw_3P_L_csv.GetColumn(0);
  std::vector<string> L_s2_Tw_3P_2p_s = Tw_3P_L_csv.GetColumn(1);
  std::vector<string> L_s2_Tw_3P_1p_s = Tw_3P_L_csv.GetColumn(2);
  std::vector<string> L_s2_Tw_3P_0p_s = Tw_3P_L_csv.GetColumn(3);


  TString Tw_1P_Rname = "time_walk/DB/RHRS_diff_tw_1P_4668_S2.csv";
  CsvParser Tw_1P_R_csv(Tw_1P_Rname.Data());

  std::vector<string> R_s2_Tw_1P_slope_s = Tw_1P_R_csv.GetColumn(0);
  std::vector<string> R_s2_Tw_1P_off_s = Tw_1P_R_csv.GetColumn(1);

  
  TString Tw_3P_Rname = "time_walk/DB/RHRS_diff_tw_3P_4668_S2.csv";
  CsvParser Tw_3P_R_csv(Tw_3P_Lname.Data());

  std::vector<string> R_s2_Tw_3P_3p_s = Tw_3P_R_csv.GetColumn(0);
  std::vector<string> R_s2_Tw_3P_2p_s = Tw_3P_R_csv.GetColumn(1);
  std::vector<string> R_s2_Tw_3P_1p_s = Tw_3P_R_csv.GetColumn(2);
  std::vector<string> R_s2_Tw_3P_0p_s = Tw_3P_R_csv.GetColumn(3);


  if(L_s2_Tw_1P_slope_s.size() != NS2Pad || L_s2_Tw_1P_off_s.size() != NS2Pad || L_s2_Tw_3P_3p_s.size() != NS2Pad || L_s2_Tw_3P_2p_s.size() != NS2Pad || L_s2_Tw_3P_1p_s.size() != NS2Pad || L_s2_Tw_3P_0p_s.size() != NS2Pad || R_s2_Tw_1P_slope_s.size() != NS2Pad || R_s2_Tw_1P_off_s.size() != NS2Pad || R_s2_Tw_3P_3p_s.size() != NS2Pad || R_s2_Tw_3P_2p_s.size() != NS2Pad || R_s2_Tw_3P_1p_s.size() != NS2Pad || R_s2_Tw_3P_0p_s.size() != NS2Pad){
    cout << "s2 Timewalk DB csv files do not have " << NS2Pad << " entries!!" << endl;
  }



  double L_s2_Tw_1P_slope[NS2Pad] = {0.0};
  double L_s2_Tw_1P_off[NS2Pad] = {0.0};
  
  double L_s2_Tw_3P_3p[NS2Pad] = {0.0};
  double L_s2_Tw_3P_2p[NS2Pad] = {0.0};
  double L_s2_Tw_3P_1p[NS2Pad] = {0.0};
  double L_s2_Tw_3P_0p[NS2Pad] = {0.0};

  double R_s2_Tw_1P_slope[NS2Pad] = {0.0};
  double R_s2_Tw_1P_off[NS2Pad] = {0.0};
  
  double R_s2_Tw_3P_3p[NS2Pad] = {0.0};
  double R_s2_Tw_3P_2p[NS2Pad] = {0.0};
  double R_s2_Tw_3P_1p[NS2Pad] = {0.0};
  double R_s2_Tw_3P_0p[NS2Pad] = {0.0};

  
  
  for(Int_t i = 0; i<NS2Pad ; i++){
    L_s2_Tw_1P_slope[i] = stod(L_s2_Tw_1P_slope_s[i]);
    L_s2_Tw_1P_off[i] = stod(L_s2_Tw_1P_off_s[i]);

    L_s2_Tw_3P_3p[i] = stod(L_s2_Tw_3P_3p_s[i]);
    L_s2_Tw_3P_2p[i] = stod(L_s2_Tw_3P_2p_s[i]);
    L_s2_Tw_3P_1p[i] = stod(L_s2_Tw_3P_1p_s[i]);
    L_s2_Tw_3P_0p[i] = stod(L_s2_Tw_3P_0p_s[i]);

    R_s2_Tw_1P_slope[i] = stod(R_s2_Tw_1P_slope_s[i]);
    R_s2_Tw_1P_off[i] = stod(R_s2_Tw_1P_off_s[i]);

    R_s2_Tw_3P_3p[i] = stod(R_s2_Tw_3P_3p_s[i]);
    R_s2_Tw_3P_2p[i] = stod(R_s2_Tw_3P_2p_s[i]);
    R_s2_Tw_3P_1p[i] = stod(R_s2_Tw_3P_1p_s[i]);
    R_s2_Tw_3P_0p[i] = stod(R_s2_Tw_3P_0p_s[i]);
  }


  cout << "Test reading of LHRS tw" << endl;

  for(Int_t i = 0; i<NS2Pad ; i++){
    cout << "Channel " << i << endl;
    cout << "slope = " << L_s2_Tw_1P_slope[i] << ", offset = " << L_s2_Tw_1P_off[i] << endl;
    cout << "3P = " << L_s2_Tw_3P_3p[i] << ", 2P = " << L_s2_Tw_3P_2p[i] << ", 1P = " << L_s2_Tw_3P_1p[i] << ", 0P = " << L_s2_Tw_3P_0p[i] << endl;
  } 
  
  cout << endl << endl;
  cout << "Test reading of RHRS tw" << endl;

  for(Int_t i = 0; i<NS2Pad ; i++){
    cout << "Channel " << i << endl;
    cout << "slope = " << R_s2_Tw_1P_slope[i] << ", offset = " << R_s2_Tw_1P_off[i] << endl;
    cout << "3P = " << R_s2_Tw_3P_3p[i] << ", 2P = " << R_s2_Tw_3P_2p[i] << ", 1P = " << R_s2_Tw_3P_1p[i] << ", 0P = " << R_s2_Tw_3P_0p[i] << endl;
  }

  cout << endl << endl;



  

  // define cuts


  
  // variabls used for cutting
  Double_t L_tr_n,L_cer_asum_c,L_ps_e,L_sh_e;
  Double_t L_tr_p[100],L_s0_trx[100],L_s2_try[100];

  Double_t R_tr_n,R_cer_asum_c,R_ps_e,R_sh_e;
  Double_t R_ps_asum_c, R_sh_asum_c;
  Double_t R_tr_p[100],R_s0_trx[100],R_s2_try[100];
  
  
  Double_t L_s2_trdx[NS2Pad];
  Double_t L_s2_lt[NS2Pad],L_s2_rt[NS2Pad];
  Double_t L_s2_lt_c[NS2Pad],L_s2_rt_c[NS2Pad];
  Double_t L_s2_lt_o[NS2Pad],L_s2_rt_o[NS2Pad];
  Double_t L_s2_lt_t[NS2Pad],L_s2_rt_t[NS2Pad];
  Double_t L_s2_nthit;
  Double_t L_s2_t_pads[NS2Pad];
  Double_t L_s2_la_p[NS2Pad], L_s2_ra_p[NS2Pad];
  Double_t L_s2_la_c[NS2Pad], L_s2_ra_c[NS2Pad];;

  // pl variables
  Double_t L_r_x, L_r_y, L_tg_dp, L_r_th, L_r_ph;

  

  Double_t R_s2_lt[NS2Pad],R_s2_rt[NS2Pad];
  Double_t R_s2_lt_c[NS2Pad],R_s2_rt_c[NS2Pad];
  Double_t R_s2_lt_o[NS2Pad],R_s2_rt_o[NS2Pad];
  Double_t R_s2_lt_t[NS2Pad],R_s2_rt_t[NS2Pad];
  Double_t R_s2_nthit;
  Double_t R_s2_t_pads[NS2Pad];
  Double_t R_s2_la_p[NS2Pad], R_s2_ra_p[NS2Pad];
  Double_t R_s2_la_c[NS2Pad], R_s2_ra_c[NS2Pad];;

  // pl variables
  Double_t R_r_x, R_r_y, R_tg_dp, R_r_th, R_r_ph;

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
  
  T->SetBranchStatus("R.tr.n",1);
  T->SetBranchStatus("R.tr.p",1);
  T->SetBranchStatus("R.cer.asum_c",1);
  T->SetBranchStatus("R.ps.e",1);
  T->SetBranchStatus("R.sh.e",1);

  
  
  T->SetBranchStatus("L.s2.nthit",1);
  T->SetBranchStatus("L.s2.t_pads",1);
  T->SetBranchStatus("L.s2.lt",1);
  T->SetBranchStatus("L.s2.rt",1);
  T->SetBranchStatus("L.s2.lt_c",1);
  T->SetBranchStatus("L.s2.rt_c",1);
  T->SetBranchStatus("L.s2.lt_o",1);
  T->SetBranchStatus("L.s2.rt_o",1);
  T->SetBranchStatus("L.s2.lt_t",1);
  T->SetBranchStatus("L.s2.rt_t",1);
  T->SetBranchStatus("L.s2.la_p",1);
  T->SetBranchStatus("L.s2.la_c",1);
  T->SetBranchStatus("L.s2.ra_p",1);
  T->SetBranchStatus("L.s2.ra_c",1);
  T->SetBranchStatus("L.s2.nthit",1);
  T->SetBranchStatus("L.s2.trdx",1);

  T->SetBranchStatus("L.tr.r_x",1);
  T->SetBranchStatus("L.tr.r_y",1);
  T->SetBranchStatus("L.tr.r_th",1);
  T->SetBranchStatus("L.tr.r_ph",1);
  T->SetBranchStatus("L.tr.tg_dp",1);
  
  T->SetBranchStatus("R.s2.nthit",1);
  T->SetBranchStatus("R.s2.t_pads",1);
  T->SetBranchStatus("R.s2.lt",1);
  T->SetBranchStatus("R.s2.rt",1);
  T->SetBranchStatus("R.s2.lt_c",1);
  T->SetBranchStatus("R.s2.rt_c",1);
  T->SetBranchStatus("R.s2.lt_o",1);
  T->SetBranchStatus("R.s2.rt_o",1);
  T->SetBranchStatus("R.s2.lt_t",1);
  T->SetBranchStatus("R.s2.rt_t",1);
  T->SetBranchStatus("R.s2.la_p",1);
  T->SetBranchStatus("R.s2.ra_p",1);


  T->SetBranchStatus("R.tr.r_x",1);
  T->SetBranchStatus("R.tr.r_y",1);
  T->SetBranchStatus("R.tr.r_th",1);
  T->SetBranchStatus("R.tr.r_ph",1);
  T->SetBranchStatus("R.tr.tg_dp",1);

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
  T->SetBranchAddress("R.ps.asum_c",&R_ps_asum_c);
  T->SetBranchAddress("R.sh.asum_c",&R_sh_asum_c);
  T->SetBranchAddress("R.ps.e",&R_ps_e);
  T->SetBranchAddress("R.sh.e",&R_sh_e);
  
  T->SetBranchAddress("L.s2.lt",L_s2_lt);
  T->SetBranchAddress("L.s2.rt",L_s2_rt);
  T->SetBranchAddress("L.s2.lt_c",L_s2_lt_c);
  T->SetBranchAddress("L.s2.rt_c",L_s2_rt_c);
  T->SetBranchAddress("L.s2.lt_o",L_s2_lt_o);
  T->SetBranchAddress("L.s2.rt_o",L_s2_rt_o);
  T->SetBranchAddress("L.s2.lt_t",L_s2_lt_t);
  T->SetBranchAddress("L.s2.rt_t",L_s2_rt_t);
  T->SetBranchAddress("L.s2.nthit",&L_s2_nthit);
  T->SetBranchAddress("L.s2.trdx",L_s2_trdx);
  T->SetBranchAddress("L.s2.t_pads",L_s2_t_pads);
  T->SetBranchAddress("L.s2.la_p",L_s2_la_p);
  T->SetBranchAddress("L.s2.ra_p",L_s2_ra_p);
  T->SetBranchAddress("L.s2.la_c",L_s2_la_c);
  T->SetBranchAddress("L.s2.ra_c",L_s2_ra_c);

  T->SetBranchAddress("L.tr.r_x",&L_r_x);
  T->SetBranchAddress("L.tr.r_y",&L_r_y);
  T->SetBranchAddress("L.tr.r_th",&L_r_th);
  T->SetBranchAddress("L.tr.r_ph",&L_r_ph);
  T->SetBranchAddress("L.tr.tg_dp",&L_tg_dp);
  
  T->SetBranchAddress("R.s2.lt",R_s2_lt);
  T->SetBranchAddress("R.s2.rt",R_s2_rt);
  T->SetBranchAddress("R.s2.lt_c",R_s2_lt_c);
  T->SetBranchAddress("R.s2.rt_c",R_s2_rt_c);
  T->SetBranchAddress("R.s2.lt_o",R_s2_lt_o);
  T->SetBranchAddress("R.s2.rt_o",R_s2_rt_o);
  T->SetBranchAddress("R.s2.lt_t",R_s2_lt_t);
  T->SetBranchAddress("R.s2.rt_t",R_s2_rt_t);
  T->SetBranchAddress("R.s2.nthit",&R_s2_nthit);
  T->SetBranchAddress("R.s2.t_pads",R_s2_t_pads);
  T->SetBranchAddress("R.s2.la_p",R_s2_la_p);
  T->SetBranchAddress("R.s2.ra_p",R_s2_ra_p);
  T->SetBranchAddress("R.s2.la_c",R_s2_la_c);
  T->SetBranchAddress("R.s2.ra_c",R_s2_ra_c);

  T->SetBranchAddress("R.tr.r_x",&R_r_x);
  T->SetBranchAddress("R.tr.r_y",&R_r_y);
  T->SetBranchAddress("R.tr.r_th",&R_r_th);
  T->SetBranchAddress("R.tr.r_ph",&R_r_ph);
  T->SetBranchAddress("R.tr.tg_dp",&R_tg_dp);


  T->SetBranchAddress("DR.rrawt2", Trig_time_T2); 
  T->SetBranchAddress("Ndata.DR.rrawt2", &Trig_No_T2);

  T->SetBranchAddress("DR.rrawt5",Trig_time_T5);
  T->SetBranchAddress("Ndata.DR.rrawt5",&Trig_No_T5);

  
  Int_t nentries = T->GetEntries();
  //  nentries = 1e5;

  
  TCut LCut = "L.tr.n==1 && abs(L.s2.trdx)<0.07 && (L.prl1.e/(L.gold.p*1000))>0.3 && ((L.prl1.e+L.prl2.e)/(L.gold.p*1000))>0.625 && ((L.prl1.e+L.prl2.e)/(L.gold.p*1000))<1.11 &&  L.cer.asum_c >650";

  TCut L_s2_Cut_Overall;
  
  TCut RCut = "R.tr.n==1 && (R.ps.asum_c+0.9**R.sh.asum_c)>800 && R.ps.asum_c>350 && R.cer.asum_c>200";

  
  // parameters for S2 and trigger timing cuts

  Double_t T2LowCut = 1683;
  Double_t T2HighCut = 1710;

  Double_t T5LowCut = 1625;
  Double_t T5HighCut = 1634;

  Bool_t T2Pass = false;
  Bool_t T5Pass = false;
  

  TCut L_s2_Cut[NS2Pad] ;
  TCut R_s2_Cut[NS2Pad] ;

  // cuts on difference of timewalk effect between L-PMT and R-PMT

  Double_t TW_lim = 0.007;
  

  for(Int_t i = 0; i<NS2Pad; i++){

    L_s2_Cut[i] = Form("L.s2.t_pads==%i &&  L.s2.la_c[%i]>100 && L.s2.ra_c[%i]>100",i,i,i);
    if (i ==0){
      L_s2_Cut_Overall =  L_s2_Cut[i];
    }
    else{
      L_s2_Cut_Overall = L_s2_Cut_Overall || L_s2_Cut[i];
    }
    //    LCut = Form("%s || %s",LCut,L_s2_Cut[i]);

    R_s2_Cut[i] = Form("R.s2.nthit==1 && R.s2.la_c[%i]>100 && R.s2.ra_c[%i]>100",i,i);
  }
  
  LCut += L_s2_Cut_Overall;
  
  cout << "Overall Lcut = " << LCut << endl;


  Int_t coinc_bins = 20;
  // Double_t coinc_start = 0;
  // Double_t coinc_end = 2e-7;
  // Double_t coinc_start = 1.5e-7; // tof
  // Double_t coinc_end = 1.6e-7; // tof

  // Double_t coinc_start = 1.6e-7; // align
  // Double_t coinc_end = 1.7e-7; // align

  Double_t coinc_start = 1.45e-7; //uncorrected
  Double_t coinc_end = 1.55e-7; // uncorrected


  

  

  // variables lims

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

  

  Double_t coinc_start_corr = (coinc_start-coin_off)*1e9;
  Double_t coinc_end_corr = (coinc_end-coin_off)*1e9;

  
  
  TH1F *h1 = new TH1F("h1"," S2-Time difference (corrected)",coinc_bins,coinc_start_corr,coinc_end_corr);

  
  TH2F* h_coinc_s2_pad = new TH2F("h_coinc_s2_pad","s2 left + right time vs paddle",16,-0.5,15.5,coinc_bins,coinc_start_corr,coinc_end_corr);
  
  TH2F* h_coinc_s2_x = new TH2F("h_coinc_s2_x","s2 left + right time vs X",20,L_x_low,L_x_up,coinc_bins,coinc_start_corr,coinc_end_corr);

  TH2F* h_coinc_s2_th = new TH2F("h_coinc_s2_th","s2 left + right time vs L_th",20,L_th_low,L_th_up,coinc_bins,coinc_start_corr,coinc_end_corr);

  TH2F* h_coinc_s2_ph = new TH2F("h_coinc_s2_ph","s2 left + right time vs L_ph",20,L_ph_low,L_ph_up,coinc_bins,coinc_start_corr,coinc_end_corr);

  TH2F* h_coinc_s2_dp = new TH2F("h_coinc_s2_dp","s2 left + right time vs L_dp",20,L_dp_low,L_dp_up,coinc_bins,coinc_start_corr,coinc_end_corr);
  
  // both arm pmt (difference) TW lims
  Double_t tw_b_low = -0.0075;
  Double_t tw_b_high = 0.0075;
  Int_t tw_b_nbins = 50;

  // both arm pmt (difference) TW lims
  Double_t tw_bs_low = -0.05;
  Double_t tw_bs_high = 0.05;
  Int_t tw_bs_nbins = 50;

  TH2F* h_coinc_s2_TW = new TH2F("h_coinc_s2_TW","time vs L_TW (L-R)",tw_b_nbins,tw_b_low,tw_b_high,coinc_bins,coinc_start_corr,coinc_end_corr);
  h_coinc_s2_TW->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}} - 1/#sqrt{ADC_{R}})_{LHRS}");
  h_coinc_s2_TW->GetYaxis()->SetTitle("t_{L,corr} - t_{R,corr} (ns)");


  TH2F* h_coinc_s2_TW_both = new TH2F("h_coinc_s2_TW_both","time vs L_TW (L-R)",tw_bs_nbins,tw_bs_low,tw_bs_high,coinc_bins,coinc_start_corr,coinc_end_corr);
  h_coinc_s2_TW_both->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}} - 1/#sqrt{ADC_{R}})_{LHRS} - (1/#sqrt{ADC_{L}} - 1/#sqrt{ADC_{R}})_{RHRS}");
  h_coinc_s2_TW_both->GetYaxis()->SetTitle("t_{L,corr} - t_{R,corr} (ns)");


  // single arm pmt (sum) TW lims
  Double_t tw_p_low = 0.05;
  Double_t tw_p_high = 0.15;
  Int_t tw_p_nbins = 69;
    

  // both arm pmt (sum) TW lims
  Double_t tw_ps_low = -0.05;
  Double_t tw_ps_high = 0.05;
  Int_t tw_ps_nbins = 69;

  TH2F* h_coinc_s2_TW_p = new TH2F("h_coinc_s2_TW_p","time vs L_TW (L+R)",tw_p_nbins,tw_p_low,tw_p_high,coinc_bins,coinc_start_corr,coinc_end_corr);
  h_coinc_s2_TW_p->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}} + 1/#sqrt{ADC_{R}})_{LHRS}");
  h_coinc_s2_TW_p->GetYaxis()->SetTitle("t_{L,corr} - t_{R,corr} (ns)");


  TH2F* h_coinc_s2_TW_p_both = new TH2F("h_coinc_s2_TW_p_both","time vs L_TW (L+R)",tw_ps_nbins,tw_ps_low,tw_ps_high,coinc_bins,coinc_start_corr,coinc_end_corr);
  h_coinc_s2_TW_p_both->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}} + 1/#sqrt{ADC_{R}})_{LHRS} - (1/#sqrt{ADC_{L}} + 1/#sqrt{ADC_{R}})_{RHRS}");
  h_coinc_s2_TW_p_both->GetYaxis()->SetTitle("t_{L,corr} - t_{R,corr} (ns)");

  // single pmt TW lims
  Double_t tw_s_low = 0.03;
  Double_t tw_s_high = 0.07;
  Int_t tw_s_nbins = 50;
  
  TH2F* h_coinc_s2_TW_L = new TH2F("h_coinc_s2_TW_L","time vs L_TW (L)",tw_s_nbins,tw_s_low,tw_s_high,coinc_bins,coinc_start_corr,coinc_end_corr);
  h_coinc_s2_TW_L->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}})");
  h_coinc_s2_TW_L->GetYaxis()->SetTitle("t_{L,corr} - t_{R,corr} (ns)");

  TH2F* h_coinc_s2_TW_R = new TH2F("h_coinc_s2_TW_R","time vs L_TW (R)",tw_s_nbins,tw_s_low,tw_s_high,coinc_bins,coinc_start_corr,coinc_end_corr);
  h_coinc_s2_TW_R->GetXaxis()->SetTitle("(1/#sqrt{ADC_{R}})");
  h_coinc_s2_TW_R->GetYaxis()->SetTitle("t_{L,corr} - t_{R,corr} (ns)");

  TH2F* h_coinc_TW_vs = new TH2F("h_coinc_TW_vs","L_TW (L-R) vs L_TW (L+R)",tw_p_nbins,tw_p_low,tw_p_high,tw_b_nbins,tw_b_low,tw_b_high);
  h_coinc_TW_vs->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}} + 1/#sqrt{ADC_{R}})_{LHRS}");
  h_coinc_TW_vs->GetYaxis()->SetTitle("(1/#sqrt{ADC_{L}} - 1/#sqrt{ADC_{R}})_{LHRS}");
  
  TProfile* hprof_coinc_s2_pad = new TProfile("hprof_coinc_s2_pad","s2 left + right time vs paddle",16,-0.5,15.5,coinc_start_corr,coinc_end_corr);
  
  TProfile* hprof_coinc_s2_x = new TProfile("hprof_coinc_s2_x","s2 left + right time vs X",20,L_x_low,L_x_up,coinc_start_corr,coinc_end_corr);

  TProfile* hprof_coinc_s2_th = new TProfile("hprof_coinc_s2_th","time vs L_th",20,L_th_low,L_th_up,coinc_start_corr,coinc_end_corr);

  TProfile* hprof_coinc_s2_ph = new TProfile("hprof_coinc_s2_ph","time vs L_ph",20,L_ph_low,L_ph_up,coinc_start_corr,coinc_end_corr);

  TProfile* hprof_coinc_s2_dp = new TProfile("hprof_coinc_s2_dp","time vs L_dp",20,L_dp_low,L_dp_up,coinc_start_corr,coinc_end_corr);
  
  TProfile* hprof_coinc_s2_TW = new TProfile("hprof_coinc_s2_TW","time vs L_Tw (L-R)",20,tw_b_low,tw_b_high,coinc_start_corr,coinc_end_corr);
  hprof_coinc_s2_TW->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}} - 1/#sqrt{ADC_{R}})_{LHRS}");
  hprof_coinc_s2_TW->GetYaxis()->SetTitle("t_{L,corr} - t_{R,corr} (ns)");

  TProfile* hprof_coinc_s2_TW_both = new TProfile("hprof_coinc_s2_TW_both","time vs L_TW (L-R)",20,tw_bs_low,tw_bs_high,coinc_start_corr,coinc_end_corr);
  hprof_coinc_s2_TW_both->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}} - 1/#sqrt{ADC_{R}})_{LHRS} - (1/#sqrt{ADC_{L}} - 1/#sqrt{ADC_{R}})_{RHRS}");
  hprof_coinc_s2_TW_both->GetYaxis()->SetTitle("t_{L,corr} - t_{R,corr} (ns)");
  

  TProfile* hprof_coinc_s2_TW_p = new TProfile("hprof_coinc_s2_TW_p","time vs L_TW (p)",20,tw_p_low,tw_p_high,coinc_start_corr,coinc_end_corr);
  hprof_coinc_s2_TW_p->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}} + 1/#sqrt{ADC_{R}})_{LHRS}");
  hprof_coinc_s2_TW_p->GetYaxis()->SetTitle("t_{L,corr} - t_{R,corr} (ns)");

  TProfile* hprof_coinc_s2_TW_p_both = new TProfile("hprof_coinc_s2_TW_p_both","time vs L_TW (L+R)",20,tw_ps_low,tw_ps_high,coinc_start_corr,coinc_end_corr);
  hprof_coinc_s2_TW_p_both->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}} + 1/#sqrt{ADC_{R}})_{LHRS} - (1/#sqrt{ADC_{L}} + 1/#sqrt{ADC_{R}})_{RHRS}");
  hprof_coinc_s2_TW_p_both->GetYaxis()->SetTitle("t_{L,corr} - t_{R,corr} (ns)");

  
  TProfile* hprof_coinc_s2_TW_L = new TProfile("hprof_coinc_s2_TW_L","time vs L_TW (L)",20,tw_s_low,tw_s_high,coinc_start_corr,coinc_end_corr);
  hprof_coinc_s2_TW_L->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}})");
  hprof_coinc_s2_TW_L->GetYaxis()->SetTitle("t_{L,corr} - t_{R,corr} (ns)");
  

  TProfile* hprof_coinc_s2_TW_R = new TProfile("hprof_coinc_s2_TW_R","time vs L_TW (R)",20,tw_s_low,tw_s_high,coinc_start_corr,coinc_end_corr);
  hprof_coinc_s2_TW_R->GetXaxis()->SetTitle("(1/#sqrt{ADC_{R}})");
  hprof_coinc_s2_TW_R->GetYaxis()->SetTitle("t_{L,corr} - t_{R,corr} (ns)");


  TProfile* hprof_coinc_TW_vs = new TProfile("hprof_coinc_TW_vs","L_TW (L-R) vs L_TW (L+R)",20,tw_p_low,tw_p_high,tw_b_low,tw_b_high);
  hprof_coinc_TW_vs->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}} + 1/#sqrt{ADC_{R}})_{LHRS}");
  hprof_coinc_TW_vs->GetYaxis()->SetTitle("(1/#sqrt{ADC_{L}} - 1/#sqrt{ADC_{R}})_{LHRS}");
  
  
  // RHRS vars against correlation time
  
  TH2F* h_coincR_s2_pad = new TH2F("h_coincR_s2_pad","s2 left + right time vs R paddle",16,-0.5,15.5,coinc_bins,coinc_start_corr,coinc_end_corr);
  
  TH2F* h_coincR_s2_x = new TH2F("h_coincR_s2_x","s2 left + right time vs R_X",20,L_x_low,L_x_up,coinc_bins,coinc_start_corr,coinc_end_corr);

  TH2F* h_coincR_s2_th = new TH2F("h_coincR_s2_th","time vs R_th",20,R_th_low,R_th_up,coinc_bins,coinc_start_corr,coinc_end_corr);

  TH2F* h_coincR_s2_ph = new TH2F("h_coincR_s2_ph","time vs R_ph",20,R_ph_low,R_ph_up,coinc_bins,coinc_start_corr,coinc_end_corr);

  TH2F* h_coincR_s2_dp = new TH2F("h_coincR_s2_dp","xtime vs R_dp",20,R_dp_low,R_dp_up,coinc_bins,coinc_start_corr,coinc_end_corr);

  TH2F* h_coincR_s2_TW = new TH2F("h_coincR_s2_TW","time vs R_TW",tw_b_nbins,tw_b_low,tw_b_high,coinc_bins,coinc_start_corr,coinc_end_corr);
  h_coincR_s2_TW->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}} - 1/#sqrt{ADC_{R}})_{RHRS}");
  h_coincR_s2_TW->GetYaxis()->SetTitle("t_{L,corr} - t_{R,corr} (ns)");

  TH2F* h_coincR_s2_TW_p = new TH2F("h_coincR_s2_TW_p","time vs R_TW (p)",tw_p_nbins,tw_p_low,tw_p_high,coinc_bins,coinc_start_corr,coinc_end_corr);
  h_coincR_s2_TW_p->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}} + 1/#sqrt{ADC_{R}})_{RHRS}");
  h_coincR_s2_TW_p->GetYaxis()->SetTitle("t_{L,corr} - t_{R,corr} (ns)");
  
  TH2F* h_coincR_s2_TW_L = new TH2F("h_coincR_s2_TW_L","time vs R_TW (L)",tw_s_nbins,tw_s_low,tw_s_high,coinc_bins,coinc_start_corr,coinc_end_corr);
  h_coincR_s2_TW_L->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}})");
  h_coincR_s2_TW_L->GetYaxis()->SetTitle("t_{L,corr} - t_{R,corr} (ns)");

  TH2F* h_coincR_s2_TW_R = new TH2F("h_coincR_s2_TW_R","time vs R_TW (R)",tw_s_nbins,tw_s_low,tw_s_high,coinc_bins,coinc_start_corr,coinc_end_corr);
  h_coincR_s2_TW_R->GetXaxis()->SetTitle("(1/#sqrt{ADC_{R}})");
  h_coincR_s2_TW_R->GetYaxis()->SetTitle("t_{L,corr} - t_{R,corr} (ns)");

  TH2F* h_coincR_TW_vs = new TH2F("h_coincR_TW_vs","R_TW (L-R) vs R_TW (L+R)",tw_p_nbins,tw_p_low,tw_p_high,tw_b_nbins,tw_b_low,tw_b_high);
  h_coincR_TW_vs->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}} + 1/#sqrt{ADC_{R}})_{RHRS}");
  h_coincR_TW_vs->GetYaxis()->SetTitle("(1/#sqrt{ADC_{L}} - 1/#sqrt{ADC_{R}})_{RHRS}");
  
  
  TProfile* hprof_coincR_s2_pad = new TProfile("hprof_coincR_s2_pad","s2 left + right time vs R paddle",16,-0.5,15.5,coinc_start_corr,coinc_end_corr);
  
  TProfile* hprof_coincR_s2_x = new TProfile("hprof_coincR_s2_x","s2 left + right time vs R_X",20,L_x_low,L_x_up,coinc_start_corr,coinc_end_corr);

  TProfile* hprof_coincR_s2_th = new TProfile("hprof_coincR_s2_th","time vs R_th",20,R_th_low,R_th_up,coinc_start_corr,coinc_end_corr);

  TProfile* hprof_coincR_s2_ph = new TProfile("hprof_coincR_s2_ph","time vs R_ph",20,R_ph_low,R_ph_up,coinc_start_corr,coinc_end_corr);

  TProfile* hprof_coincR_s2_dp = new TProfile("hprof_coincR_s2_dp","time vs R_dp",20,R_dp_low,R_dp_up,coinc_start_corr,coinc_end_corr);

  TProfile* hprof_coincR_s2_TW = new TProfile("hprof_coincR_s2_TW","time vs R_TW",20,tw_b_low,tw_b_high,coinc_start_corr,coinc_end_corr);
  hprof_coincR_s2_TW->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}} - 1/#sqrt{ADC_{R}})_{RHRS}");
  hprof_coincR_s2_TW->GetYaxis()->SetTitle("t_{L,corr} - t_{R,corr} (ns)");

  TProfile* hprof_coincR_s2_TW_p = new TProfile("hprof_coincR_s2_TW_p","time vs R_TW (p)",20,tw_p_low,tw_p_high,coinc_start_corr,coinc_end_corr);
  hprof_coincR_s2_TW_p->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}} + 1/#sqrt{ADC_{R}})_{RHRS}");
  hprof_coincR_s2_TW_p->GetYaxis()->SetTitle("t_{L,corr} - t_{R,corr} (ns)");  
  
  TProfile* hprof_coincR_s2_TW_L = new TProfile("hprof_coincR_s2_TW_L","time vs R_TW (L)",20,tw_s_low,tw_s_high,coinc_start_corr,coinc_end_corr);
  hprof_coincR_s2_TW_L->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}})");
  hprof_coincR_s2_TW_L->GetYaxis()->SetTitle("t_{L,corr} - t_{R,corr} (ns)");

  TProfile* hprof_coincR_s2_TW_R = new TProfile("hprof_coincR_s2_TW_R","time vs R_TW (R)",20,tw_s_low,tw_s_high,coinc_start_corr,coinc_end_corr);
  hprof_coincR_s2_TW_R->GetXaxis()->SetTitle("(1/#sqrt{ADC_{R}})");
  hprof_coincR_s2_TW_R->GetYaxis()->SetTitle("t_{L,corr} - t_{R,corr} (ns)");
  
  TProfile* hprof_coincR_TW_vs = new TProfile("hprof_coincR_TW_vs","R_TW (L-R) vs R_TW (L+R)",20,tw_p_low,tw_p_high,tw_b_low,tw_b_high);
  hprof_coincR_TW_vs->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}} + 1/#sqrt{ADC_{R}})_{RHRS}");
  hprof_coincR_TW_vs->GetYaxis()->SetTitle("(1/#sqrt{ADC_{L}} - 1/#sqrt{ADC_{R}})_{RHRS}");  

  TH2F* h_coincCorr_s2_pad = new TH2F("h_coincCorr_s2_pad","s2 left + right time vs paddle",16,-0.5,15.5,coinc_bins,coinc_start_corr,coinc_end_corr);
  
  TH2F* h_coincCorr_s2_x = new TH2F("h_coincCorr_s2_x","s2 left + right time vs X",150,L_x_low,L_x_up,coinc_bins,coinc_start_corr,coinc_end_corr);

  TH2F* h_coincCorr_s2_th = new TH2F("h_coincCorr_s2_th","s2 left + right time vs L_th",50,L_th_low,L_th_up,coinc_bins,coinc_start_corr,coinc_end_corr);

  TH2F* h_coincCorr_s2_ph = new TH2F("h_coincCorr_s2_ph","s2 left + right time vs L_ph",50,L_ph_low,L_ph_up,coinc_bins,coinc_start_corr,coinc_end_corr);

  TH2F* h_coincCorr_s2_dp = new TH2F("h_coincCorr_s2_dp","s2 left + right time vs L_dp",150,L_dp_low,L_dp_up,coinc_bins,coinc_start_corr,coinc_end_corr);

  
  TProfile* hprof_coincCorr_s2_pad = new TProfile("hprof_coincCorr_s2_pad","s2 left + right time vs paddle",16,-0.5,15.5,coinc_start_corr,coinc_end_corr);
  
  TProfile* hprof_coincCorr_s2_x = new TProfile("hprof_coincCorr_s2_x","s2 left + right time vs X",20,L_x_low,L_x_up,coinc_start_corr,coinc_end_corr);

  TProfile* hprof_coincCorr_s2_th = new TProfile("hprof_coincCorr_s2_th","s2 left + right time vs L_th",20,L_th_low,L_th_up,coinc_start_corr,coinc_end_corr);

  TProfile* hprof_coincCorr_s2_ph = new TProfile("hprof_coincCorr_s2_ph","s2 left + right time vs L_ph",20,L_ph_low,L_ph_up,coinc_start_corr,coinc_end_corr);

  TProfile* hprof_coincCorr_s2_dp = new TProfile("hprof_coincCorr_s2_dp","s2 left + right time vs L_dp",20,L_dp_low,L_dp_up,coinc_start_corr,coinc_end_corr);


  // RHRS vars against correlation time

  
  TH2F* h_coincCorrR_s2_pad = new TH2F("h_coincCorrR_s2_pad","s2 left + right time vs R paddle",16,-0.5,15.5,coinc_bins,coinc_start_corr,coinc_end_corr);
  
  TH2F* h_coincCorrR_s2_x = new TH2F("h_coincCorrR_s2_x","s2 left + right time vs R_X",20,L_x_low,L_x_up,coinc_bins,coinc_start_corr,coinc_end_corr);

  TH2F* h_coincCorrR_s2_th = new TH2F("h_coincCorrR_s2_th","time vs R_th",20,R_th_low,R_th_up,coinc_bins,coinc_start_corr,coinc_end_corr);

  TH2F* h_coincCorrR_s2_ph = new TH2F("h_coincCorrR_s2_ph","time vs R_ph",20,R_ph_low,R_ph_up,coinc_bins,coinc_start_corr,coinc_end_corr);

  TH2F* h_coincCorrR_s2_dp = new TH2F("h_coincCorrR_s2_dp","xtime vs R_dp",20,R_dp_low,R_dp_up,coinc_bins,coinc_start_corr,coinc_end_corr);

  
  TProfile* hprof_coincCorrR_s2_pad = new TProfile("hprof_coincCorrR_s2_pad","s2 left + right time vs R paddle",16,-0.5,15.5,coinc_start_corr,coinc_end_corr);
  
  TProfile* hprof_coincCorrR_s2_x = new TProfile("hprof_coincCorrR_s2_x","s2 left + right time vs R_X",20,L_x_low,L_x_up,coinc_start_corr,coinc_end_corr);

  TProfile* hprof_coincCorrR_s2_th = new TProfile("hprof_coincCorrR_s2_th","time vs R_th",20,R_th_low,R_th_up,coinc_start_corr,coinc_end_corr);

  TProfile* hprof_coincCorrR_s2_ph = new TProfile("hprof_coincCorrR_s2_ph","time vs R_ph",20,R_ph_low,R_ph_up,coinc_start_corr,coinc_end_corr);

  TProfile* hprof_coincCorrR_s2_dp = new TProfile("hprof_coincCorrR_s2_dp","time vs R_dp",20,R_dp_low,R_dp_up,coinc_start_corr,coinc_end_corr);



  // 1D variable plots

  
  TH1F* h_L_th = new TH1F("h_L_th","Theta",50,L_th_low,L_th_up);

  TH1F* h_L_ph = new TH1F("h_L_ph","Phi",50,L_ph_low,L_ph_up);

  TH1F* h_L_x = new TH1F("h_L_x","x",50,L_x_low,L_x_up);

  TH1F* h_L_dp = new TH1F("h_L_dp","dp",50,L_dp_low,L_dp_up);
  
  TH1F* h_L_TW = new TH1F("h_L_TW","TW (L-R)",tw_b_nbins,tw_b_low,tw_b_high);
  h_L_TW->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}} - 1/#sqrt{ADC_{R}})_{LHRS}");

  TH1F* h_L_TW_p = new TH1F("h_L_TW_p","TW (L+R)",tw_p_nbins,tw_p_low,tw_p_high);
  h_L_TW_p->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}} + 1/#sqrt{ADC_{R}})_{LHRS}");
  
  TH1F* h_L_TW_l = new TH1F("h_L_TW_l","TW (l)",tw_s_nbins,tw_s_low,tw_s_high);
  h_L_TW_l->GetXaxis()->SetTitle("1/#sqrt{ADC_{L}}");

  TH1F* h_L_TW_r = new TH1F("h_L_TW_r","TW (r)",tw_s_nbins,tw_s_low,tw_s_high);
  h_L_TW_r->GetXaxis()->SetTitle("1/#sqrt{ADC_{R}}");

  TH1F* h_TW_both = new TH1F("h_TW_both","TW (L-R, both arms)",tw_bs_nbins,tw_bs_low,tw_bs_high);
  h_TW_both->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}} - 1/#sqrt{ADC_{R}})_{LHRS} - (1/#sqrt{ADC_{L}} - 1/#sqrt{ADC_{R}})_{RHRS}");
  
  TH1F* h_TW_p_both = new TH1F("h_TW_p_both","TW (L+R, both arms)",tw_ps_nbins,tw_ps_low,tw_ps_high);
  h_TW_p_both->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}} + 1/#sqrt{ADC_{R}})_{LHRS} - (1/#sqrt{ADC_{L}} + 1/#sqrt{ADC_{R}})_{RHRS}");

  
  

  TH1F* h_R_th = new TH1F("h_R_th","Theta",50,R_th_low,R_th_up);

  TH1F* h_R_ph = new TH1F("h_R_ph","Phi",50,R_ph_low,R_ph_up);

  TH1F* h_R_x = new TH1F("h_R_x","x",50,R_x_low,R_x_up);

  TH1F* h_R_dp = new TH1F("h_R_dp","dp",50,R_dp_low,R_dp_up);

  TH1F* h_R_TW = new TH1F("h_R_TW","TW (L-R)",tw_b_nbins,tw_b_low,tw_b_high);
  h_R_TW->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}} - 1/#sqrt{ADC_{R}})_{RHRS}");

  TH1F* h_R_TW_p = new TH1F("h_R_TW_p","TW (L+R)",tw_p_nbins,tw_p_low,tw_p_high);
  h_R_TW_p->GetXaxis()->SetTitle("(1/#sqrt{ADC_{L}} + 1/#sqrt{ADC_{R}})_{RHRS}");
  
  TH1F* h_R_TW_l = new TH1F("h_R_TW_l","TW (l)",tw_s_nbins,tw_s_low,tw_s_high);
  h_R_TW_l->GetXaxis()->SetTitle("1/#sqrt{ADC_{L}}");
  
  TH1F* h_R_TW_r = new TH1F("h_R_TW_r","TW (r)",tw_s_nbins,tw_s_low,tw_s_high);
  h_R_TW_r->GetXaxis()->SetTitle("1/#sqrt{ADC_{R}}");
  
  Double_t L_start = 1.375e3;
  Double_t L_end = 1.42e3;
  Int_t L_bins = 90;
  
  Double_t L_start_corr = (L_start - L_off);
  Double_t L_end_corr = (L_end - L_off);
  
  
  TH1F* h_L_s2_t = new TH1F("h_L_s2_t","s2 left + right time",L_bins,L_start_corr,L_end_corr);
  

  TH2F* h_L_s2_pad = new TH2F("h_L_s2_pad","s2 left vs paddle",16,-0.5,15.5,L_bins,L_start_corr,L_end_corr);
  
  TH2F* h_L_s2_x = new TH2F("h_L_s2_x","s2 left time vs X",300,L_x_low,L_x_up,L_bins,L_start_corr,L_end_corr);

  TH2F* h_L_s2_th = new TH2F("h_L_s2_th","L time vs L_th",100,L_th_low,L_th_up,L_bins,L_start_corr,L_end_corr);

  TH2F* h_L_s2_dp = new TH2F("h_L_s2_dp","L time vs L_dp",100,R_dp_low,R_dp_up,L_bins,L_start_corr,L_end_corr);




  TProfile* hprof_L_s2_x = new TProfile("hprof_L_s2_x","s2 left + right time vs X",50,L_x_low,L_x_up,L_start_corr,L_end_corr);

  TProfile* hprof_L_s2_pad = new TProfile("hprof_L_s2_pad","s2 left + right time vs paddle",16,-0.5,15.5,L_start_corr,L_end_corr);
  
  TProfile* hprof_L_s2_th = new TProfile("hprof_L_s2_th","L time vs L_th",50,L_th_low,L_th_up,L_start_corr,L_end_corr);

  TProfile* hprof_L_s2_dp = new TProfile("hprof_L_s2_dp","L time vs L_dp",50,R_dp_low,R_dp_up,L_start_corr,L_end_corr);
  

  

  
  //300,1.378e-6,1.388e-6
  TH1F* h_L_s2_un_t = new TH1F("h_L_s2_un_t","s2 left + right time",L_bins,L_start,L_end);

  TH2F* h_L_s2_un_pad = new TH2F("h_L_s2_un_pad","s2_un left + right time vs paddle",16,-0.5,15.5,L_bins,L_start,L_end);
  
  TH2F* h_L_s2_un_x = new TH2F("h_L_s2_un_x","s2 left + right time vs X",300,L_x_low,L_x_up,L_bins,L_start,L_end);

  TH2F* h_L_s2_un_th = new TH2F("h_L_s2_un_th","L time vs L_th",100,L_th_low,L_th_up,L_bins,L_start,L_end);

  TH2F* h_L_s2_un_dp = new TH2F("h_L_s2_un_dp","L time vs L_dp",100,R_dp_low,R_dp_up,L_bins,L_start,L_end);


  TProfile* hprof_L_s2_un_x = new TProfile("hprof_L_s2_un_x","s2 left + right time vs X",50,L_x_low,L_x_up,L_start,L_end);

  TProfile* hprof_L_s2_un_pad = new TProfile("hprof_L_s2_un_pad","s2 left + right time vs paddle",16,-0.5,15.5,L_start,L_end);
  
  TProfile* hprof_L_s2_un_th = new TProfile("hprof_L_s2_un_th","R time vs R_th",50,L_th_low,L_th_up,L_start,L_end);

  TProfile* hprof_L_s2_un_dp = new TProfile("hprof_L_s2_un_dp","R time vs R_dp",50,R_dp_low,R_dp_up,L_start,L_end);

  
  // marker
  // 90,1.24e-6,1.31e-6
  Double_t R_start = 1.24e3;
  Double_t R_end = 1.31e3;
  Int_t R_bins = 140;

  Double_t R_start_corr = (R_start - R_off);
  Double_t R_end_corr = (R_end - R_off);

  
  TH1F* h_R_s2_t = new TH1F("h_R_s2_t","s2 left + right time",R_bins,R_start_corr,R_end_corr);
  

  TH2F* h_R_s2_pad = new TH2F("h_R_s2_pad","s2 left vs paddle",16,-0.5,15.5,R_bins,R_start_corr,R_end_corr);
  
  TH2F* h_R_s2_x = new TH2F("h_R_s2_x","s2 right time vs X",300,L_x_low,L_x_up,R_bins,R_start_corr,R_end_corr);

  TH2F* h_R_s2_th = new TH2F("h_R_s2_th","R time vs R_th",100,R_th_low,R_th_up,R_bins,R_start_corr,R_end_corr);

  TH2F* h_R_s2_dp = new TH2F("h_R_s2_dp","R time vs R_dp",100,R_dp_low,R_dp_up,R_bins,R_start_corr,R_end_corr);




  TProfile* hprof_R_s2_x = new TProfile("hprof_R_s2_x","s2 left + right time vs X",50,L_x_low,L_x_up,R_start_corr,R_end_corr);

  TProfile* hprof_R_s2_pad = new TProfile("hprof_R_s2_pad","s2 left + right time vs paddle",16,-0.5,15.5,R_start_corr,R_end_corr);
  
  TProfile* hprof_R_s2_th = new TProfile("hprof_R_s2_th","R time vs R_th",50,R_th_low,R_th_up,R_start_corr,R_end_corr);

  TProfile* hprof_R_s2_dp = new TProfile("hprof_R_s2_dp","R time vs R_dp",50,R_dp_low,R_dp_up,R_start_corr,R_end_corr);
  

  
  //300,1.378e-6,1.388e-6
  TH1F* h_R_s2_un_t = new TH1F("h_R_s2_un_t","s2 left + right time",R_bins,R_start,R_end);

  TH2F* h_R_s2_un_pad = new TH2F("h_R_s2_un_pad","s2_un left + right time vs paddle",16,-0.5,15.5,R_bins,R_start,R_end);
  
  TH2F* h_R_s2_un_x = new TH2F("h_R_s2_un_x","s2 left + right time vs X",300,L_x_low,L_x_up,R_bins,R_start,R_end);

  TH2F* h_R_s2_un_th = new TH2F("h_R_s2_un_th","R time vs R_th",100,R_th_low,R_th_up,R_bins,R_start,R_end);

  TH2F* h_R_s2_un_dp = new TH2F("h_R_s2_un_dp","R time vs R_dp",100,R_dp_low,R_dp_up,R_bins,R_start,R_end);


  TProfile* hprof_R_s2_un_x = new TProfile("hprof_R_s2_un_x","s2 left + right time vs X",50,L_x_low,L_x_up,R_start,R_end);

  TProfile* hprof_R_s2_un_pad = new TProfile("hprof_R_s2_un_pad","s2 left + right time vs paddle",16,-0.5,15.5,R_start,R_end);
  
  TProfile* hprof_R_s2_un_th = new TProfile("hprof_R_s2_un_th","R time vs R_th",50,R_th_low,R_th_up,R_start,R_end);

  TProfile* hprof_R_s2_un_dp = new TProfile("hprof_R_s2_un_dp","R time vs R_dp",50,R_dp_low,R_dp_up,R_start,R_end);

  

  
  
  Double_t LTime;
  Double_t LTime_un;
  
  Double_t RTime;
  Double_t RTime_un;


  Double_t L_TW_b = 0.0; // difference between left and right paddle ADCs TW effects
  
  Double_t LTime_tw_1p = 0.0;
  Double_t L_tw_1p_corr = 0.0;

  Double_t LTime_tw_3p = 0.0;
  Double_t L_tw_3p_corr = 0.0;

  Double_t R_TW_b = 0.0; // difference between left and right paddle ADCs TW effects

  Double_t RTime_tw_1p = 0.0;
  Double_t R_tw_1p_corr = 0.0;

  Double_t RTime_tw_3p = 0.0;
  Double_t R_tw_3p_corr = 0.0;


  // test pl coeffecient corrections


  // s2 align

  //  Double_t L_x_coeff =  -1.77234e-09;
  Double_t L_x_coeff =  -3.25953e-09;
  Double_t L_th_coeff = 2.62758e-08;
  Double_t L_ph_coeff = -2.50384e-08;

 
  //  Double_t R_x_coeff =  -2.72124e-09;
  Double_t R_x_coeff =  -3.04253e-09;
  Double_t R_th_coeff = 2.19428e-08;
  Double_t R_ph_coeff = 7.81744e-09;


  
  
  //  Double_t R_ph_coeff = 0;
  



  
  //        f_line_ADC[i] = new TF1(Form("f_line_ADC_%i",i),"pol1",0,0.16);
  
  
  // Record on no of hits in each paddle
  Int_t L_hits[NS2Pad] = {0};
  Int_t R_hits[NS2Pad] = {0};


  // record which paddle is hit for an entry
  Int_t LHRS_pad = -1;
  Int_t RHRS_pad = -1;
    for(Int_t i=0;i<nentries;i++){
  // for(Int_t i=0;i<100;i++){
    T->GetEntry(i);
    
    LTime = 0.0;
    LTime_un = 0.0;
    LHRS_pad = -1;

    RTime = 0.0;
    RTime_un = 0.0;
    RHRS_pad = -1;

    T2Pass = false;
    T5Pass = false;

    L_TW_b = 0.0;
    
    L_tw_1p_corr = 0.0;
    LTime_tw_1p = 0.0;

    R_TW_b = 0.0;
    
    R_tw_1p_corr = 0.0;
    RTime_tw_1p = 0.0;
    
    
    
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
    
    
    
    if(L_s2_nthit==1 &&  L_tr_n==1 && L_cer_asum_c>1500 && (L_ps_e+L_sh_e)/(1000.*L_tr_p[0])>0.8  && R_tr_n==1 && R_cer_asum_c>200 &&  R_s2_nthit==1 && (R_ps_asum_c + 0.9*R_sh_asum_c)>800 && R_ps_asum_c>350 && T2Pass && T5Pass && Trig_type == 6){
    
      for(int j=0;j<NS2Pad;j++){                                                                              



	if (L_s2_t_pads[0]==j && (L_s2_lt_c[j]+L_s2_rt_c[j]) != 0){    
	  LTime = fTdc2T*(L_s2_lt[j]-L_ls2_coeff[j] + L_s2_rt[j]-L_rs2_coeff[j])/2.;	  
	  LTime_un = (L_s2_lt_t[j]+L_s2_rt_t[j])/2.;
	  LTime *= 1e9;
	  LTime_un *= 1e9;	  
	  LHRS_pad = j;
	}
	
	

	if (j == 4 || j==5 || j ==7){
	  //if(L_s2_t_pads[0]==j){
	  if(L_s2_t_pads[0]==j && (L_s2_lt[j]!= 0 || L_s2_rt[j] != 0)  && (L_s2_la_c[j]>100 || L_s2_ra_c[j]>100) ){
	    LTime = fTdc2T*(L_s2_lt[j]-L_ls2_coeff[j] + L_s2_rt[j]-L_rs2_coeff[j])/2.;	  
	    LTime_un = (L_s2_lt_t[j]+L_s2_rt_t[j])/2.;
	    LTime *= 1e9;
	    LTime_un *= 1e9;	  
	    LHRS_pad = j;
	  }
	}
      	//   else{
	//     if(L_s2_t_pads[0]==j &&  (L_s2_lt_c[j]+L_s2_rt_c[j]) != 0 && L_s2_la_c[j]>100 && L_s2_ra_c[j]>100){    
	//     LTime = (L_s2_lt_o[j]+L_s2_rt_o[j])/2.;
	//     // cout << "LTime == " << LTime << endl;
	//     // cout << "L_s2_lt_o[" << j << "] = " << L_s2_lt_o[j] << endl;
	//     // cout << "L_s2_rt_o[" << j << "] = " << L_s2_rt_o[j] << endl;
	//     LTime_un = (L_s2_lt_t[j]+L_s2_rt_t[j])/2.;
	    
	//     LHRS_pad = j;
	//     // cout << "LHRS_pad = " << LHRS_pad << endl;

	//     }
	//   }


      
      

      
	if (R_s2_t_pads[0]==j && R_s2_lt_c[j]+R_s2_rt_c[j] != 0){    
	  //	  RTime = (R_s2_lt_o[j]+R_s2_rt_o[j])/2.;
	  RTime = fTdc2T*(R_s2_lt[j]-R_ls2_coeff[j] + R_s2_rt[j]-R_rs2_coeff[j])/2.;
	  RTime_un = (R_s2_lt_t[j]+R_s2_rt_t[j])/2.;
	  RTime *= 1e9;
	  RTime_un *= 1e9;
	  RHRS_pad = j;
	}
	

      }
	



	// if (R_s2_t_pads[0]==j && R_s2_lt_c[j]+R_s2_rt_c[j] != 0){    
	//   RTime = (R_s2_lt_c[j]+R_s2_rt_c[j])/2.;
	//   RTime_un = (fTdc2T*(R_s2_lt[j]+R_s2_rt[j]))/2.;
	//   RHRS_pad = j+1;
	// }

      
    }
    else{
      continue;
    }
    
    
    
    
    
    if( (LTime) > 0 ){

	L_hits[LHRS_pad] += 1;
	
	// cout << "Filling: LHRS_pad = " << LHRS_pad << endl;
	
	Double_t LPad = LHRS_pad;
	
	Double_t RPad = RHRS_pad;
	// cout << "LPad = " << LPad << ": LTime = " << LTime << endl;


	if(RTime > 0 ){

	  Double_t L_TW_l = 1/TMath::Sqrt(L_s2_la_c[LHRS_pad]);
	  Double_t L_TW_r = 1/TMath::Sqrt(L_s2_ra_c[LHRS_pad]);
	  L_TW_b = L_TW_l - L_TW_r;
	  Double_t L_TW_p = L_TW_l + L_TW_r;
	  
	  
	  Double_t R_TW_l = 1/TMath::Sqrt(R_s2_la_c[RHRS_pad]);
	  Double_t R_TW_r = 1/TMath::Sqrt(R_s2_ra_c[RHRS_pad]);
	  R_TW_b = R_TW_l - R_TW_r;
	  Double_t R_TW_p = R_TW_l + R_TW_r;


	  //  if( abs(L_TW_b) < TW_lim && abs(R_TW_b) < TW_lim){


	  L_tw_1p_corr = L_s2_Tw_1P_slope[LHRS_pad]*(L_TW_b +  L_s2_Tw_1P_off[LHRS_pad]);
	  LTime_tw_1p = LTime - L_tw_1p_corr;
	  
	  R_tw_1p_corr = R_s2_Tw_1P_slope[RHRS_pad]*(R_TW_b +  R_s2_Tw_1P_off[RHRS_pad]);
	  RTime_tw_1p = RTime - 0.5*R_tw_1p_corr;
	  
	  Double_t Coinc_time = LTime-RTime;
	  
	  Double_t Coinc_timeCorr = (LTime-L_x_coeff*L_r_x-L_ph_coeff*L_r_ph-L_th_coeff*L_r_th)-(RTime-R_x_coeff*R_r_x-R_ph_coeff*R_r_ph-R_th_coeff*R_r_th);
	  
	  
	  h1->Fill(LTime-RTime);
	  h_coinc_s2_pad->Fill(LPad,LTime-RTime);
	  h_coinc_s2_x->Fill(L_r_x,LTime-RTime);	  
	  h_coinc_s2_th->Fill(L_r_th,LTime-RTime);
	  h_coinc_s2_ph->Fill(L_r_ph,LTime-RTime);	  
	  h_coinc_s2_dp->Fill(L_tg_dp,LTime-RTime);
	  h_coinc_s2_TW->Fill(L_TW_b,Coinc_timeCorr);
	  h_coinc_s2_TW_p->Fill(L_TW_p,Coinc_timeCorr);
	  h_coinc_s2_TW_L->Fill(L_TW_l,Coinc_timeCorr);
	  h_coinc_s2_TW_R->Fill(L_TW_r,Coinc_timeCorr);
	  h_coinc_TW_vs->Fill(L_TW_p,L_TW_b);
	  hprof_coinc_s2_pad->Fill(LPad,LTime-RTime);
	  hprof_coinc_s2_x->Fill(L_r_x,LTime-RTime);
	  hprof_coinc_s2_th->Fill(L_r_th,LTime-RTime);
	  hprof_coinc_s2_ph->Fill(L_r_ph,LTime-RTime);
	  hprof_coinc_s2_dp->Fill(L_tg_dp,LTime-RTime);
	  hprof_coinc_s2_TW->Fill(L_TW_b,Coinc_timeCorr);
	  hprof_coinc_s2_TW_p->Fill(L_TW_p,Coinc_timeCorr);
	  hprof_coinc_s2_TW_L->Fill(L_TW_l,Coinc_timeCorr);
	  hprof_coinc_s2_TW_R->Fill(L_TW_r,Coinc_timeCorr);
	  hprof_coinc_TW_vs->Fill(L_TW_p,L_TW_b);

	    
	  h_coincCorr_s2_pad->Fill(LPad,(LTime-L_x_coeff*L_r_x-L_ph_coeff*L_r_ph-L_th_coeff*L_r_th)-(RTime-R_x_coeff*R_r_x-R_ph_coeff*R_r_ph-R_th_coeff*R_r_th));
	  h_coincCorr_s2_pad->Fill(LPad,(LTime-L_x_coeff*L_r_x-L_ph_coeff*L_r_ph-L_th_coeff*L_r_th)-(RTime-R_x_coeff*R_r_x-R_ph_coeff*R_r_ph-R_th_coeff*R_r_th));
	  h_coincCorr_s2_x->Fill(L_r_x,(LTime-L_x_coeff*L_r_x-L_ph_coeff*L_r_ph-L_th_coeff*L_r_th)-(RTime-R_x_coeff*R_r_x-R_ph_coeff*R_r_ph-R_th_coeff*R_r_th));	  
	  h_coincCorr_s2_th->Fill(L_r_th,(LTime-L_x_coeff*L_r_x-L_ph_coeff*L_r_ph-L_th_coeff*L_r_th)-(RTime-R_x_coeff*R_r_x-R_ph_coeff*R_r_ph-R_th_coeff*R_r_th));
	  h_coincCorr_s2_ph->Fill(L_r_ph,(LTime-L_x_coeff*L_r_x-L_ph_coeff*L_r_ph-L_th_coeff*L_r_th)-(RTime-R_x_coeff*R_r_x-R_ph_coeff*R_r_ph-R_th_coeff*R_r_th));	  
	  h_coincCorr_s2_dp->Fill(L_tg_dp,(LTime-L_x_coeff*L_r_x-L_ph_coeff*L_r_ph-L_th_coeff*L_r_th)-(RTime-R_x_coeff*R_r_x-R_ph_coeff*R_r_ph-R_th_coeff*R_r_th));
	  hprof_coincCorr_s2_pad->Fill(LPad,(LTime-L_x_coeff*L_r_x-L_ph_coeff*L_r_ph-L_th_coeff*L_r_th)-(RTime-R_x_coeff*R_r_x-R_ph_coeff*R_r_ph-R_th_coeff*R_r_th));
	  hprof_coincCorr_s2_x->Fill(L_r_x,(LTime-L_x_coeff*L_r_x-L_ph_coeff*L_r_ph-L_th_coeff*L_r_th)-(RTime-R_x_coeff*R_r_x-R_ph_coeff*R_r_ph-R_th_coeff*R_r_th));
	  hprof_coincCorr_s2_th->Fill(L_r_th,(LTime-L_x_coeff*L_r_x-L_ph_coeff*L_r_ph-L_th_coeff*L_r_th)-(RTime-R_x_coeff*R_r_x-R_ph_coeff*R_r_ph-R_th_coeff*R_r_th));
	  hprof_coincCorr_s2_ph->Fill(L_r_ph,(LTime-L_x_coeff*L_r_x-L_ph_coeff*L_r_ph-L_th_coeff*L_r_th)-(RTime-R_x_coeff*R_r_x-R_ph_coeff*R_r_ph-R_th_coeff*R_r_th));
	  hprof_coincCorr_s2_dp->Fill(L_tg_dp,(LTime-L_x_coeff*L_r_x-L_ph_coeff*L_r_ph-L_th_coeff*L_r_th)-(RTime-R_x_coeff*R_r_x-R_ph_coeff*R_r_ph-R_th_coeff*R_r_th));



	  h_L_th->Fill(L_r_th);
	  h_L_ph->Fill(L_r_ph);
	  h_L_x->Fill(L_r_x);
	  h_L_dp->Fill(L_tg_dp);
	  h_L_TW->Fill(L_TW_b);
	  h_L_TW_p->Fill(L_TW_p);
	  h_L_TW_l->Fill(L_TW_l);
	  h_L_TW_r->Fill(L_TW_r);

	  h_TW_both->Fill(L_TW_b-R_TW_b);
	  h_TW_p_both->Fill(L_TW_p-R_TW_p);

	  h_coinc_s2_TW_both->Fill(L_TW_b-R_TW_b,Coinc_timeCorr);
	  h_coinc_s2_TW_p_both->Fill(L_TW_p-R_TW_p,Coinc_timeCorr);

	  hprof_coinc_s2_TW_both->Fill(L_TW_b-R_TW_b,Coinc_timeCorr);
	  hprof_coinc_s2_TW_p_both->Fill(L_TW_p-R_TW_p,Coinc_timeCorr);

	  h_R_th->Fill(R_r_th);
	  h_R_ph->Fill(R_r_ph);
	  h_R_x->Fill(R_r_x);
	  h_R_dp->Fill(R_tg_dp);
	  h_R_TW->Fill(R_TW_b);
	  h_R_TW_p->Fill(R_TW_p);
	  h_R_TW_l->Fill(R_TW_l);
	  h_R_TW_r->Fill(R_TW_r);

	  
	  
	  h_coincR_s2_pad->Fill(RPad,LTime-RTime);
	  h_coincR_s2_x->Fill(R_r_x,LTime-RTime);	  
	  h_coincR_s2_th->Fill(R_r_th,LTime-RTime);
	  h_coincR_s2_ph->Fill(R_r_ph,LTime-RTime);	  
	  h_coincR_s2_dp->Fill(R_tg_dp,LTime-RTime);
	  h_coincR_s2_TW->Fill(R_TW_b,Coinc_timeCorr);
	  h_coincR_s2_TW_p->Fill(R_TW_p,Coinc_timeCorr);
	  h_coincR_s2_TW_L->Fill(R_TW_l,Coinc_timeCorr);
	  h_coincR_s2_TW_R->Fill(R_TW_r,Coinc_timeCorr);
	  h_coincR_TW_vs->Fill(R_TW_p,R_TW_b);
	  hprof_coincR_s2_pad->Fill(RPad,LTime-RTime);
	  hprof_coincR_s2_x->Fill(R_r_x,LTime-RTime);
	  hprof_coincR_s2_th->Fill(R_r_th,LTime-RTime);
	  hprof_coincR_s2_ph->Fill(R_r_ph,LTime-RTime);
	  hprof_coincR_s2_dp->Fill(R_tg_dp,LTime-RTime);
	  hprof_coincR_s2_TW->Fill(R_TW_b,Coinc_timeCorr);
	  hprof_coincR_s2_TW_p->Fill(R_TW_p,Coinc_timeCorr);
	  hprof_coincR_s2_TW_L->Fill(R_TW_l,Coinc_timeCorr);
	  hprof_coincR_s2_TW_R->Fill(R_TW_r,Coinc_timeCorr);
	  hprof_coincR_TW_vs->Fill(R_TW_p,R_TW_b);
	  
	  h_coincCorrR_s2_pad->Fill(RPad,(LTime-L_x_coeff*L_r_x-L_ph_coeff*L_r_ph-L_th_coeff*L_r_th)-(RTime-R_x_coeff*R_r_x-R_ph_coeff*R_r_ph-R_th_coeff*R_r_th));
	  h_coincCorrR_s2_x->Fill(R_r_x,(LTime-L_x_coeff*L_r_x-L_ph_coeff*L_r_ph-L_th_coeff*L_r_th)-(RTime-R_x_coeff*R_r_x-R_ph_coeff*R_r_ph-R_th_coeff*R_r_th));	  
	  h_coincCorrR_s2_th->Fill(R_r_th,(LTime-L_x_coeff*L_r_x-L_ph_coeff*L_r_ph-L_th_coeff*L_r_th)-(RTime-R_x_coeff*R_r_x-R_ph_coeff*R_r_ph-R_th_coeff*R_r_th));
	  h_coincCorrR_s2_ph->Fill(R_r_ph,(LTime-L_x_coeff*L_r_x-L_ph_coeff*L_r_ph-L_th_coeff*L_r_th)-(RTime-R_x_coeff*R_r_x-R_ph_coeff*R_r_ph-R_th_coeff*R_r_th));	  
	  h_coincCorrR_s2_dp->Fill(R_tg_dp,(LTime-L_x_coeff*L_r_x-L_ph_coeff*L_r_ph-L_th_coeff*L_r_th)-(RTime-R_x_coeff*R_r_x-R_ph_coeff*R_r_ph-R_th_coeff*R_r_th));
	  hprof_coincCorrR_s2_pad->Fill(RPad,(LTime-L_x_coeff*L_r_x-L_ph_coeff*L_r_ph-L_th_coeff*L_r_th)-(RTime-R_x_coeff*R_r_x-R_ph_coeff*R_r_ph-R_th_coeff*R_r_th));
	  hprof_coincCorrR_s2_x->Fill(R_r_x,(LTime-L_x_coeff*L_r_x-L_ph_coeff*L_r_ph-L_th_coeff*L_r_th)-(RTime-R_x_coeff*R_r_x-R_ph_coeff*R_r_ph-R_th_coeff*R_r_th));
	  hprof_coincCorrR_s2_th->Fill(R_r_th,(LTime-L_x_coeff*L_r_x-L_ph_coeff*L_r_ph-L_th_coeff*L_r_th)-(RTime-R_x_coeff*R_r_x-R_ph_coeff*R_r_ph-R_th_coeff*R_r_th));
	  hprof_coincCorrR_s2_ph->Fill(R_r_ph,(LTime-L_x_coeff*L_r_x-L_ph_coeff*L_r_ph-L_th_coeff*L_r_th)-(RTime-R_x_coeff*R_r_x-R_ph_coeff*R_r_ph-R_th_coeff*R_r_th));
	  hprof_coincCorrR_s2_dp->Fill(R_tg_dp,(LTime-L_x_coeff*L_r_x-L_ph_coeff*L_r_ph-L_th_coeff*L_r_th)-(RTime-R_x_coeff*R_r_x-R_ph_coeff*R_r_ph-R_th_coeff*R_r_th));

	  R_hits[RHRS_pad] += 1;
	}
    
	
	h_L_s2_t->Fill(LTime);
	

	
	h_L_s2_pad->Fill(LPad*1.0,LTime);
	h_L_s2_x->Fill(L_r_x,LTime);
	h_L_s2_th->Fill(L_r_th,LTime);
	h_L_s2_dp->Fill(L_tg_dp,LTime);
	hprof_L_s2_pad->Fill(LPad*1.0,LTime);
	hprof_L_s2_x->Fill(L_r_x,LTime);
	hprof_L_s2_th->Fill(L_r_th,LTime);
	hprof_L_s2_dp->Fill(L_tg_dp,LTime);
	
	h_L_s2_un_t->Fill(LTime_un);
	h_L_s2_un_pad->Fill(LPad*1.0,LTime_un);
	h_L_s2_un_x->Fill(L_r_x,LTime_un);
	h_L_s2_un_th->Fill(L_r_th,LTime_un);
	h_L_s2_un_dp->Fill(L_tg_dp,LTime_un);
	hprof_L_s2_un_pad->Fill(LPad*1.0,LTime_un);
	hprof_L_s2_un_x->Fill(L_r_x,LTime_un);
	hprof_L_s2_un_th->Fill(L_r_th,LTime_un);
	hprof_L_s2_un_dp->Fill(L_tg_dp,LTime_un);

	// marker
	
	h_R_s2_pad->Fill(RPad*1.0,RTime);
	h_R_s2_x->Fill(R_r_x,RTime);
	h_R_s2_th->Fill(R_r_th,RTime);
	h_R_s2_dp->Fill(R_tg_dp,RTime);
	hprof_R_s2_pad->Fill(RPad*1.0,RTime);
	hprof_R_s2_x->Fill(R_r_x,RTime);
	hprof_R_s2_th->Fill(R_r_th,RTime);
	hprof_R_s2_dp->Fill(R_tg_dp,RTime);
	
	h_R_s2_un_t->Fill(RTime_un);
	h_R_s2_un_pad->Fill(RPad*1.0,RTime_un);
	h_R_s2_un_x->Fill(R_r_x,RTime_un);
	h_R_s2_un_th->Fill(R_r_th,RTime_un);
	h_R_s2_un_dp->Fill(R_tg_dp,RTime_un);
	hprof_R_s2_un_pad->Fill(RPad*1.0,RTime_un);
	hprof_R_s2_un_x->Fill(R_r_x,RTime_un);
	hprof_R_s2_un_th->Fill(R_r_th,RTime_un);
	hprof_R_s2_un_dp->Fill(R_tg_dp,RTime_un);

	
	
      }
      
  }
    
   
      

  

  // first draw timing distrib 

  // TCanvas* c1 = new TCanvas("c1","s2 L Timing");

  // c1->cd(1);
  // h_L_s2_t->Draw();



  TCanvas* c2 = new TCanvas("c2","s2 L Timing vs pl variables");

  c2->Divide(4,2);
  c2->cd(1);
  h_L_s2_x->Draw("colz");

  TF1* f_L_x_sin = new TF1("f_L_x_sin","pol1",L_x_low,L_x_up);
  h_L_s2_x->Fit("f_L_x_sin","QR","",L_x_low,L_x_up);
  Double_t L_x_sin_slope = f_L_x_sin->GetParameter(1);

  cout << "L x slope = " << L_x_sin_slope << endl;


  c2->cd(2);
  h_L_s2_th->Draw("colz");

  c2->cd(3);
  h_L_s2_dp->Draw("colz");

  c2->cd(4);
  h_L_s2_pad->Draw("colz");

  c2->cd(5);
  hprof_L_s2_x->Draw();

  c2->cd(6);
  hprof_L_s2_th->Draw();

  c2->cd(7);
  hprof_L_s2_dp->Draw();

  c2->cd(8);
  hprof_L_s2_pad->Draw();

  TCanvas* c3 = new TCanvas("c3","s2 L Timing vs pl variables (uncorrected)");

  c3->Divide(4,2);
  c3->cd(1);
  h_L_s2_un_x->Draw("colz");

  c3->cd(2);
  h_L_s2_un_th->Draw("colz");

  c3->cd(3);
  h_L_s2_un_dp->Draw("colz");

  c3->cd(4);
  h_L_s2_un_pad->Draw("colz");
  
  c3->cd(5);
  hprof_L_s2_un_x->Draw();

  c3->cd(6);
  hprof_L_s2_un_th->Draw();

  c3->cd(7);
  hprof_L_s2_un_dp->Draw();

  c3->cd(8);
  hprof_L_s2_un_pad->Draw();
  


  
  TCanvas* c2_R = new TCanvas("c2_R","s2 R Timing vs pl variables");

  c2_R->Divide(4,2);
  c2_R->cd(1);
  h_R_s2_x->Draw("colz");

  TF1* f_R_x_sin = new TF1("f_R_x_sin","pol1",L_x_low,L_x_up);
  h_R_s2_x->Fit("f_R_x_sin","QR","",L_x_low,L_x_up);
  Double_t R_x_sin_slope = f_R_x_sin->GetParameter(1);

  cout << "R x slope = " << R_x_sin_slope << endl;
  
  c2_R->cd(2);
  h_R_s2_th->Draw("colz");

  c2_R->cd(3);
  h_R_s2_dp->Draw("colz");

  c2_R->cd(4);
  h_R_s2_pad->Draw("colz");

  c2_R->cd(5);
  hprof_R_s2_x->Draw();

  c2_R->cd(6);
  hprof_R_s2_th->Draw();

  c2_R->cd(7);
  hprof_R_s2_dp->Draw();

  c2_R->cd(8);
  hprof_R_s2_pad->Draw();

  TCanvas* c3_R = new TCanvas("c3_R","s2 R Timing vs pl variables (uncorrected)");

  c3_R->Divide(4,2);
  c3_R->cd(1);
  h_R_s2_un_x->Draw("colz");

  c3_R->cd(2);
  h_R_s2_un_th->Draw("colz");

  c3_R->cd(3);
  h_R_s2_un_dp->Draw("colz");

  c3_R->cd(4);
  h_R_s2_un_pad->Draw("colz");
  
  c3_R->cd(5);
  hprof_R_s2_un_x->Draw();

  c3_R->cd(6);
  hprof_R_s2_un_th->Draw();

  c3_R->cd(7);
  hprof_R_s2_un_dp->Draw();

  c3_R->cd(8);
  hprof_R_s2_un_pad->Draw();
  


  TCanvas* c4_a  = new TCanvas("c4_a","1D variables");
  c4_a->Divide(5,2);
  c4_a->cd(1);
  h_L_th->Draw();

  c4_a->cd(2);
  h_L_ph->Draw();

  c4_a->cd(3);
  h_L_x->Draw();

  c4_a->cd(4);
  h_L_dp->Draw();

  c4_a->cd(5);
  h_L_TW->Draw();
  

  c4_a->cd(6);
  h_R_th->Draw();

  c4_a->cd(7);
  h_R_ph->Draw();

  c4_a->cd(8);
  h_R_x->Draw();

  c4_a->cd(9);
  h_R_dp->Draw();

  c4_a->cd(10);
  h_R_TW->Draw();

  // plots coinc time vs path-length variables and fit to extract corrections
  
  TCanvas* c4 = new TCanvas("c4","Coinc time vs L pl variables");

  c4->Divide(6,3);

  c4->cd(1);
  h_coinc_s2_th->Draw("colz");
  
  c4->cd(2);
  h_coinc_s2_ph->Draw("colz");

  c4->cd(3);
  h_coinc_s2_x->Draw("colz");
  
  c4->cd(4);
  h_coinc_s2_dp->Draw("colz");

  c4->cd(5);
  h_coinc_s2_pad->Draw("colz");

  c4->cd(6);
  h_coinc_s2_TW->Draw("colz");


  c4->cd(7);
  hprof_coinc_s2_th->Draw();
  
  TF1* f_L_th = new TF1("f_L_th","pol1",-0.01,0.014);
  hprof_coinc_s2_th->Fit("f_L_th","QR","",-0.01,0.014);
  Double_t L_th_slope = f_L_th->GetParameter(1);
  
  c4->cd(8);
  hprof_coinc_s2_ph->Draw();

  TF1* f_L_ph = new TF1("f_L_ph","pol1",-0.02,0.02);
  hprof_coinc_s2_ph->Fit("f_L_ph","QR","",-0.02,0.02);
  Double_t L_ph_slope = f_L_ph->GetParameter(1);

  c4->cd(9);
  hprof_coinc_s2_x->Draw();
  hprof_coinc_s2_x->GetYaxis()->SetRangeUser(coinc_start_corr,coinc_end_corr);

  TF1* f_L_x = new TF1("f_L_x","pol1",-0.5,0.5);
  hprof_coinc_s2_x->Fit("f_L_x","QR","",-0.5,0.5);
  Double_t L_x_slope = f_L_x->GetParameter(1);
  
  
  c4->cd(10);
  hprof_coinc_s2_dp->Draw();

  c4->cd(11);
  hprof_coinc_s2_pad->Draw();

  c4->cd(12);
  hprof_coinc_s2_TW->Draw();

  c4->cd(13);
  h_coinc_s2_th->FitSlicesY(0,0,-1,10,"QNR");

  TH1D* h_coinc_s2_th_mean = (TH1D*)gDirectory->Get("h_coinc_s2_th_1");
  h_coinc_s2_th_mean->GetYaxis()->SetRangeUser(coinc_start_corr,coinc_end_corr);
  h_coinc_s2_th_mean->Draw();

  TF1* f_L_th_sl = new TF1("f_L_th_sl","pol1",-0.01,0.014);
  h_coinc_s2_th_mean->Fit("f_L_th_sl","QR","",-0.01,0.014);
  Double_t L_th_slope_sl = f_L_th_sl->GetParameter(1);


  c4->cd(14);
  h_coinc_s2_ph->FitSlicesY(0,0,-1,10,"QN");

  TH1D* h_coinc_s2_ph_mean = (TH1D*)gDirectory->Get("h_coinc_s2_ph_1");
  h_coinc_s2_ph_mean->GetYaxis()->SetRangeUser(coinc_start_corr,coinc_end_corr);
  h_coinc_s2_ph_mean->Draw();
  
  TF1* f_L_ph_sl = new TF1("f_L_ph_sl","pol1",-0.02,0.02);
  h_coinc_s2_ph_mean->Fit("f_L_ph_sl","QR","",-0.02,0.02);
  Double_t L_ph_slope_sl = f_L_ph_sl->GetParameter(1);

  
  c4->cd(15);
  h_coinc_s2_x->FitSlicesY(0,0,-1,10,"QN");

  TH1D* h_coinc_s2_x_mean = (TH1D*)gDirectory->Get("h_coinc_s2_x_1");
  h_coinc_s2_x_mean->GetYaxis()->SetRangeUser(coinc_start_corr,coinc_end_corr);
  h_coinc_s2_x_mean->Draw();

  TF1* f_L_x_sl = new TF1("f_L_x_sl","pol1",-0.5,0.5);
  h_coinc_s2_x_mean->Fit("f_L_x_sl","QR","",-0.5,0.5);
  Double_t L_x_slope_sl = f_L_x_sl->GetParameter(1);

    
  c4->cd(16);
  h_coinc_s2_dp->FitSlicesY(0,0,-1,10,"QN");

  TH1D* h_coinc_s2_dp_mean = (TH1D*)gDirectory->Get("h_coinc_s2_dp_1");
  h_coinc_s2_dp_mean->GetYaxis()->SetRangeUser(coinc_start_corr,coinc_end_corr);
  h_coinc_s2_dp_mean->Draw();

  c4->cd(17);
  h_coinc_s2_pad->FitSlicesY(0,0,-1,10,"QN");

  TH1D* h_coinc_s2_pad_mean = (TH1D*)gDirectory->Get("h_coinc_s2_pad_1");
  h_coinc_s2_pad_mean->GetYaxis()->SetRangeUser(coinc_start_corr,coinc_end_corr);
  h_coinc_s2_pad_mean->Draw();

  c4->cd(18);
  h_coinc_s2_TW->FitSlicesY(0,0,-1,10,"QN");

  TH1D* h_coinc_s2_TW_mean = (TH1D*)gDirectory->Get("h_coinc_s2_TW_1");
  h_coinc_s2_TW_mean->GetYaxis()->SetRangeUser(coinc_start_corr,coinc_end_corr);
  h_coinc_s2_TW_mean->Draw();



  

  

  TCanvas* c5 = new TCanvas("c5","Coinc time vs R pl variables");

  c5->Divide(6,3);
  c5->cd(1);
  h_coincR_s2_th->Draw("colz");
  
  c5->cd(2);
  h_coincR_s2_ph->Draw("colz");

  c5->cd(3);
  h_coincR_s2_x->Draw("colz");
  
  c5->cd(4);
  h_coincR_s2_dp->Draw("colz");

  c5->cd(5);
  h_coincR_s2_pad->Draw("colz");

  c5->cd(6);
  h_coincR_s2_TW->Draw("colz");

  c5->cd(7);
  hprof_coincR_s2_th->Draw();

  TF1* f_R_th = new TF1("f_R_th","pol1",-0.017,0.017);
  hprof_coincR_s2_th->Fit("f_R_th","QR","",-0.017,0.017);
  Double_t R_th_slope = f_R_th->GetParameter(1);
  
  c5->cd(8);
  hprof_coincR_s2_ph->Draw();

  TF1* f_R_ph = new TF1("f_R_ph","pol1",-0.017,0.02);
  hprof_coincR_s2_ph->Fit("f_R_ph","QR","",-0.017,0.02);
  Double_t R_ph_slope = f_R_ph->GetParameter(1);

  c5->cd(9);
  hprof_coincR_s2_x->Draw();

  TF1* f_R_x = new TF1("f_R_x","pol1",-0.45,0.4);
  hprof_coincR_s2_x->Fit("f_R_x","QR","",-0.45,0.4);
  Double_t R_x_slope = f_R_x->GetParameter(1);
  
  c5->cd(10);
  hprof_coincR_s2_dp->Draw();

  c5->cd(11);
  hprof_coincR_s2_pad->Draw();

  c5->cd(12);
  hprof_coincR_s2_TW->Draw();

  c5->cd(13);
  h_coincR_s2_th->FitSlicesY(0,0,-1,10,"QN");

  TH1D* h_coincR_s2_th_mean = (TH1D*)gDirectory->Get("h_coincR_s2_th_1");
  h_coincR_s2_th_mean->GetYaxis()->SetRangeUser(coinc_start_corr,coinc_end_corr);
  h_coincR_s2_th_mean->Draw();

  TF1* f_R_th_sl = new TF1("f_R_th_sl","pol1",-0.017,0.017);
  h_coincR_s2_th_mean->Fit("f_R_th_sl","QR","",-0.017,0.017);
  Double_t R_th_slope_sl = f_R_th_sl->GetParameter(1);


  c5->cd(14);
  h_coincR_s2_ph->FitSlicesY(0,0,-1,10,"QN");

  TH1D* h_coincR_s2_ph_mean = (TH1D*)gDirectory->Get("h_coincR_s2_ph_1");
  h_coincR_s2_ph_mean->GetYaxis()->SetRangeUser(coinc_start_corr,coinc_end_corr);
  h_coincR_s2_ph_mean->Draw();

  TF1* f_R_ph_sl = new TF1("f_R_ph_sl","pol1",-0.017,0.02);
  h_coincR_s2_ph_mean->Fit("f_R_ph_sl","QR","",-0.017,0.02);
  Double_t R_ph_slope_sl = f_R_ph_sl->GetParameter(1);

  
  c5->cd(15);
  h_coincR_s2_x->FitSlicesY(0,0,-1,10,"QN");

  TH1D* h_coincR_s2_x_mean = (TH1D*)gDirectory->Get("h_coincR_s2_x_1");
  h_coincR_s2_x_mean->GetYaxis()->SetRangeUser(coinc_start_corr,coinc_end_corr);
  h_coincR_s2_x_mean->Draw();

  TF1* f_R_x_sl = new TF1("f_R_x_sl","pol1",-0.45,0.4);
  h_coincR_s2_x_mean->Fit("f_R_x_sl","QR","",-0.45,0.4);
  Double_t R_x_slope_sl = f_R_x_sl->GetParameter(1);

    
  c5->cd(16);
  h_coincR_s2_dp->FitSlicesY(0,0,-1,10,"QN");

  TH1D* h_coincR_s2_dp_mean = (TH1D*)gDirectory->Get("h_coincR_s2_dp_1");
  h_coincR_s2_dp_mean->GetYaxis()->SetRangeUser(coinc_start_corr,coinc_end_corr);
  h_coincR_s2_dp_mean->Draw();

  c5->cd(17);
  h_coincR_s2_pad->FitSlicesY(0,0,-1,10,"QN");

  TH1D* h_coincR_s2_pad_mean = (TH1D*)gDirectory->Get("h_coincR_s2_pad_1");
  h_coincR_s2_pad_mean->GetYaxis()->SetRangeUser(coinc_start_corr,coinc_end_corr);
  h_coincR_s2_pad_mean->Draw();

  c5->cd(18);
  h_coincR_s2_TW->FitSlicesY(0,0,-1,10,"QN");

  TH1D* h_coincR_s2_TW_mean = (TH1D*)gDirectory->Get("h_coincR_s2_TW_1");
  h_coincR_s2_TW_mean->GetYaxis()->SetRangeUser(coinc_start_corr,coinc_end_corr);
  h_coincR_s2_TW_mean->Draw();







  
  // plots coinc time vs different TW variables
  
  TCanvas* c6 = new TCanvas("c6","TW variables",1600,1000);

  c6->Divide(4,3);

  c6->cd(1);
  h_L_TW->Draw();

  c6->cd(2);
  h_L_TW_p->Draw();

  c6->cd(3);
  h_L_TW_l->Draw();
  
  c6->cd(4);
  h_L_TW_r->Draw();

  
  c6->cd(5);
  h_R_TW->Draw();
  
  c6->cd(6);
  h_R_TW_p->Draw();
  
  c6->cd(7);
  h_R_TW_l->Draw();
  
  c6->cd(8);
  h_R_TW_r->Draw();
 
  c6->cd(9);
  h_TW_both->Draw();
  
  c6->cd(10);
  h_TW_p_both->Draw();
  




    
  TCanvas* c7 = new TCanvas("c7","Coinc time vs L/R TW variables",1600,1000);

  c7->Divide(4,3);

  c7->cd(1);
  h_L_TW->Draw();

  c7->cd(2);
  h_L_TW_p->Draw();

  c7->cd(3);
  h_R_TW->Draw();
  
  c7->cd(4);
  h_R_TW_p->Draw();

  
  c7->cd(5);
  h_coinc_s2_TW->Draw("colz");
  
  c7->cd(6);
  h_coinc_s2_TW_p->Draw("colz");
  
  c7->cd(7);
  h_coincR_s2_TW->Draw("colz");
  
  c7->cd(8);
  h_coincR_s2_TW_p->Draw("colz");

  c7->cd(9);
  hprof_coinc_s2_TW->Draw();
  
  c7->cd(10);
  hprof_coinc_s2_TW_p->Draw();
  
  c7->cd(11);
  hprof_coincR_s2_TW->Draw();
  
  c7->cd(12);
  hprof_coincR_s2_TW_p->Draw();


  TCanvas* c8 = new TCanvas("c8","Coinc time vs L&R TW variables",1600,1000);

  c8->Divide(2,3);

  c8->cd(1);
  h_TW_both->Draw();

  c8->cd(2);
  h_TW_p_both->Draw();

  c8->cd(3);
  h_coinc_s2_TW_both->Draw("colz");

  c8->cd(4);
  h_coinc_s2_TW_p_both->Draw("colz");

  c8->cd(5);
  hprof_coinc_s2_TW_both->Draw("colz");

  c8->cd(6);
  hprof_coinc_s2_TW_p_both->Draw("colz");


  
  TCanvas* c9 = new TCanvas("c9","TW variables Vs",1600,1000);

  c9->Divide(2,2);

  c9->cd(1);
  h_coinc_TW_vs->Draw("colz");

  c9->cd(2);
  h_coincR_TW_vs->Draw("colz");

  c9->cd(3);
  hprof_coinc_TW_vs->Draw("colz");

  c9->cd(4);
  hprof_coincR_TW_vs->Draw("colz");

  
  cout << " L_hits \t Rhits" << endl;

  for(Int_t i = 0; i<NS2Pad; i++){
    
    // cout << "L_hits[" << i << "] = " << L_hits[i] << endl;
    // cout << "R_hits[" << i << "] = " << R_hits[i] << endl;
    cout << L_hits[i] << "\t" << R_hits[i] << endl;
    
  }
					

  cout << endl << endl;

  cout << "Left pl parameters: " << endl;
  cout << "th slope = " << L_th_slope << ", th slope (sl) = " << L_th_slope_sl << endl;  
  cout << "ph slope = " << L_ph_slope << ", ph slope (sl) = " << L_ph_slope_sl << endl;  
  cout << "x slope = " << L_x_slope   << ", x slope (sl) = " << L_x_slope_sl << endl;  


  cout << endl << endl;

  cout << "Right pl parameters: " << endl;
  cout << "th slope = " << R_th_slope << ", th slope (sl) = " << R_th_slope_sl << endl;  
  cout << "ph slope = " << R_ph_slope << ", ph slope (sl) = " << R_ph_slope_sl << endl;  
  cout << "x slope = " << R_x_slope   << ", x slope (sl) = " << R_x_slope_sl << endl;  

  // save to DB file

  //  CsvParser csv_L_db(DB_Lname.Data());
  
  TString output_DB = DB_Lname.Remove(DB_Lname.Last('_'),9);
  output_DB.Remove(0,8);
  
  
  ofstream oofile(Form("pl_corr/DB/%s_Coinc_vs_%i.txt",output_DB.Data(),runno));


  oofile << "Left pl parameters: " << endl;
  oofile << "th slope = " << L_th_slope << endl;
  oofile << "ph slope = " << L_ph_slope << endl;
  oofile << "x slope = " << L_x_slope << endl;


  oofile << endl << endl;

  oofile << "Right pl parameters: " << endl;
  oofile << "th slope = " << R_th_slope << endl;
  oofile << "ph slope = " << R_ph_slope << endl;
  oofile << "x slope = " << R_x_slope << endl;

  oofile.close();


  ofstream oofile_csv(Form("pl_corr/DB/%s_Coinc_vs_%i.csv",output_DB.Data(),runno));


  oofile_csv << "L th slope" << "," << " L ph slope" << "," << "L x slope" << "," << "R th slope" << "," << " R ph slope" << "," << "R x slope" << endl;
  oofile_csv << L_th_slope << "," << L_ph_slope << "," << L_x_slope << "," << R_th_slope << "," << R_ph_slope << "," << R_x_slope << endl;

  oofile_csv.close();


  ofstream oofile_csv_sl(Form("pl_corr/DB/%s_Coinc_vs_%i_sl.csv",output_DB.Data(),runno));

  oofile_csv_sl << "L th slope" << "," << " L ph slope" << "," << "L x slope" << "," << "R th slope" << "," << " R ph slope" << "," << "R x slope" << endl;
  oofile_csv_sl << L_th_slope_sl << "," << L_ph_slope_sl << "," << L_x_slope_sl << "," << R_th_slope_sl << "," << R_ph_slope_sl << "," << R_x_slope_sl << endl;

  oofile_csv_sl.close();



  // save plots

  /*
  c6->SaveAs(Form("pl_corr/plots/L_TW_corr_%i.png",runno));
  c7->SaveAs(Form("pl_corr/plots/R_TW_corr_%i.png",runno));
  */

  // pl_corr/plots/
    
  
}
