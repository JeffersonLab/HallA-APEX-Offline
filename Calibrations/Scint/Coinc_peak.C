/*
*************************************************************
13/12/20 John Williamson
Script that measures width of coincidence peak (defined by difference of left and right S2 scintillator times)

Adds width of fit to chosen csv file

*************************************************************
*/

#include "file_def.h"
#include "Load_more_rootfiles.C"
#include "CsvParser.h"

#include <iostream>

void s2TW(Int_t runno=-1);


void Coinc_peak(Int_t runno, TString DB_Lname /* LHRS DB name where corrections are read from*/, TString DB_Rname /* RHRS DB name where corrections are read from*/, TString Name  = "_" /*Name to be added to csv file*/,TString csv_name = "_", Int_t uncor_flag = 1){


  const Int_t NS2Pad = 16;
  const Int_t MaxS2Hit = 15;
  
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
  Double_t coin_off = (L_ls2_coeff_mean + L_rs2_coeff_mean) - (R_ls2_coeff_mean + R_rs2_coeff_mean); 


  // read offset DBs for LHRS and RHRS (TW versions)

  TString DB_Lname_tw = "DB/LHRS_S2_TW_align_t_4668.csv";
  
  CsvParser csv_L_db_tw(DB_Lname_tw.Data());
  
  std::vector<string> L_ls2_coeff_tw_s = csv_L_db_tw.GetColumn(0);
  std::vector<string> L_rs2_coeff_tw_s = csv_L_db_tw.GetColumn(1);

  TString DB_Rname_tw = "DB/RHRS_S2_TW_align_t_4668.csv";
  
  CsvParser csv_R_db_tw(DB_Rname_tw.Data());

  std::vector<string> R_ls2_coeff_tw_s = csv_R_db_tw.GetColumn(0);
  std::vector<string> R_rs2_coeff_tw_s = csv_R_db_tw.GetColumn(1);


  double L_ls2_coeff_tw[NS2Pad] = {0.0};
  double L_rs2_coeff_tw[NS2Pad] = {0.0};

  double R_ls2_coeff_tw[NS2Pad] = {0.0};
  double R_rs2_coeff_tw[NS2Pad] = {0.0};


  
  if(L_ls2_coeff_tw_s.size() != NS2Pad || L_rs2_coeff_tw_s.size() != NS2Pad || R_ls2_coeff_tw_s.size() != NS2Pad || R_rs2_coeff_tw_s.size() != NS2Pad){
    cout << "s2 offset DB csv files do not have " << NS2Pad << " entries!!" << endl;
  }


  Double_t L_ls2_coeff_tw_mean = 0.0;
  Double_t L_rs2_coeff_tw_mean = 0.0;
  Double_t R_ls2_coeff_tw_mean = 0.0;
  Double_t R_rs2_coeff_tw_mean = 0.0;
  
  for(Int_t i = 0; i<NS2Pad ; i++){
    L_ls2_coeff_tw[i] = stod(L_ls2_coeff_tw_s[i]);
    L_rs2_coeff_tw[i] = stod(L_rs2_coeff_tw_s[i]);
    R_ls2_coeff_tw[i] = stod(R_ls2_coeff_tw_s[i]);
    R_rs2_coeff_tw[i] = stod(R_rs2_coeff_tw_s[i]);
    
    L_ls2_coeff_tw_mean += L_ls2_coeff_tw[i];
    L_rs2_coeff_tw_mean += L_rs2_coeff_tw[i];
    R_ls2_coeff_tw_mean += R_ls2_coeff_tw[i];
    R_rs2_coeff_tw_mean += R_rs2_coeff_tw[i];

  }  
  cout << "(TW offset) DB files read " << endl;

  L_ls2_coeff_tw_mean = L_ls2_coeff_tw_mean/NS2Pad;
  L_rs2_coeff_tw_mean = L_rs2_coeff_tw_mean/NS2Pad;
  R_ls2_coeff_tw_mean = R_ls2_coeff_tw_mean/NS2Pad;
  R_rs2_coeff_tw_mean = R_rs2_coeff_tw_mean/NS2Pad;

  // from corrections get offset from uncorrected coincidence times
  //  Double_t coin_off = (L_ls2_coeff_tw_mean + L_rs2_coeff_tw_mean) - (R_ls2_coeff_tw_mean + R_rs2_coeff_tw_mean); 

  

  // read in pl-corrections if available
  
  TString output_DB = DB_Lname.Remove(DB_Lname.Last('_'),9);
  output_DB.Remove(0,8);

  //  ofstream pl_csv(Form("pl_corr/DB/%s_Coinc_vs_%i.csv",output_DB.Data(),runno));

  TString pl_name = Form("pl_corr/DB/%s_Coinc_vs_%i.csv",output_DB.Data(),runno);
  CsvParser pl_csv(pl_name.Data());

  double L_th_slope = stod((pl_csv.GetColumn(0))[0]);
  cout << "L_th_slope = " << L_th_slope << endl;
  double L_ph_slope = stod((pl_csv.GetColumn(1))[0]);
  cout << "L_ph_slope = " << L_ph_slope << endl;
  double L_x_slope = stod((pl_csv.GetColumn(2))[0]);
  cout << "L_x_slope = " << L_x_slope << endl;

  double R_th_slope = -1.*stod((pl_csv.GetColumn(3))[0]);
  cout << "R_th_slope = " << R_th_slope << endl;
  double R_ph_slope = -1.*stod((pl_csv.GetColumn(4))[0]);
  cout << "R_ph_slope = " << R_ph_slope << endl;
  double R_x_slope = -1.*stod((pl_csv.GetColumn(5))[0]);
  cout << "R_x_slope = " << R_x_slope << endl;

  cout << "PL files read" << endl;


  // read in TW correction

  TString Tw_1P_Lname = "time_walk/DB/LHRS_diff_tw_1P_4771_S2.csv";
  CsvParser Tw_1P_L_csv(Tw_1P_Lname.Data());

  std::vector<string> L_s2_Tw_1P_slope_s = Tw_1P_L_csv.GetColumn(0);
  std::vector<string> L_s2_Tw_1P_off_s = Tw_1P_L_csv.GetColumn(1);

  
  TString Tw_3P_Lname = "time_walk/DB/LHRS_diff_tw_3P_4771_S2.csv";
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


  
  //  time_walk/DB/LHRS_diff_tw_1P_4771_S2.csv
  

  
  
  std::ofstream csvFile;

  
  gStyle->SetOptFit(0011);

  gStyle->SetOptStat(0);
  
  // test if csv file exists (if not create column headers
  if ( gSystem->AccessPathName(Form("coinc_csv/%s",csv_name.Data()))) {
    
    cout << csv_name << " does not exist" << endl;
    csvFile.open(Form("coinc_csv/%s",csv_name.Data()),ios::app);
    // creat headers for csv file
    // characterise by RunNumber, Description (type of calibration performed, Width (width of coincidence peak)
    csvFile << "RunNumber" << "," << "Description" << "," << "Width" << "," << "Mean" <<  endl;
  }
  else{
    
    cout << csv_name << " does exist" << endl;
    csvFile.open(Form("coinc_csv/%s",csv_name.Data()),ios::app);

    // here check if csv file already has result from specified run and calibration description (if so then delete old and replace with new)
    //stf::ifstream ifile(Form("coinc_csv/%s",csv_name.Data()),ios::app);
    
  }


 
  // calculate width of coincidence peak
  
  TChain* T = Load_more_rootfiles(runno);
  
  
  Double_t tmin = 1690, tmax = 1705;  // the limits for 'peak'
  
  TCut TimingCut  = Form("DR.rrawt2<%.4f&&DR.rrawt2>%.4f", tmax, tmin);


  
  // Read cuts for s2 paddles based on Left ADC ped subtracted values (portion that is TW correction calibrated)

  //((1/TMath::Sqrt(R.s2.la_p[%i])))


  /*
  
  string csvname_l = ("time_walk/s2_lims/L_4771_s2.csv");
  // string csvname_l = ("time_walk/s2_lims/L_477dfsdf1_s2.csv");
  CsvParser csv_l(csvname_l);

  std::vector<string> s2_l_lower_s = csv_l.GetColumn(0);
  std::vector<string> s2_l_upper_s = csv_l.GetColumn(1);

  double s2_l_lower[NS2Pad] = {0.0};
  double s2_l_upper[NS2Pad] = {0.0};

  if(s2_l_lower_s.size() != NS2Pad || s2_l_upper_s.size() != NS2Pad){
    cout << "s2 csv files do not have " << NS2Pad << " entries!!" << endl;
  }
  
  for(Int_t i = 0; i<NS2Pad ; i++){
    s2_l_lower[i] = stod(s2_l_lower_s[i]);
    s2_l_upper[i] = stod(s2_l_upper_s[i]);
  }

  

  
  string csvname_r = ("time_walk/s2_lims/R_4651_s2.csv");
  // string csvname_r = ("time_walk/s2_lims/R_4651fasdfa_s2.csv");
  CsvParser csv_r(csvname_r);

  std::vector<string> s2_r_lower_s = csv_r.GetColumn(0);
  std::vector<string> s2_r_upper_s = csv_r.GetColumn(1);

  double s2_r_lower[NS2Pad] = {0.0};
  double s2_r_upper[NS2Pad] = {0.0};

  if(s2_r_lower_s.size() != NS2Pad || s2_r_upper_s.size() != NS2Pad){
    cout << "s2 csv files do not have " << NS2Pad << " entries!!" << endl;
  }
  
  for(Int_t i = 0; i<NS2Pad ; i++){
    s2_r_lower[i] = stod(s2_r_lower_s[i]);
    s2_r_upper[i] = stod(s2_r_upper_s[i]);
  }  

  */

  const Double_t fTdc2T = 0.5e-9;      // seconds/channel


  // variabls used for cutting
  Double_t L_tr_n,L_cer_asum_c,L_ps_e,L_sh_e;
  Double_t L_tr_p[100],L_s0_trx[100],L_s2_try[100];

  Double_t R_tr_n,R_cer_asum_c,R_ps_e,R_sh_e;
  Double_t R_tr_p[100],R_s0_trx[100],R_s2_try[100];

  // pl variables
  Double_t L_r_x[100], L_r_y[100], L_tg_dp[100], L_r_th[100], L_r_ph[100];
  Double_t R_r_x[100], R_r_y[100], R_tg_dp[100], R_r_th[100], R_r_ph[100];
  
  
  Double_t L_s2_lt[NS2Pad],L_s2_rt[NS2Pad];
  Double_t L_s2_lt_c[NS2Pad],L_s2_rt_c[NS2Pad];
  Double_t L_s2_nthit;
  Double_t L_s2_t_pads[NS2Pad];
  Double_t L_s2_la_p[NS2Pad];
  Double_t L_s2_la_c[NS2Pad];
  Double_t L_s2_ra_p[NS2Pad];
  Double_t L_s2_ra_c[NS2Pad];

  Double_t R_s2_lt[NS2Pad],R_s2_rt[NS2Pad];
  Double_t R_s2_lt_c[NS2Pad],R_s2_rt_c[NS2Pad];
  Double_t R_s2_nthit;
  Double_t R_s2_t_pads[NS2Pad];
  Double_t R_s2_la_p[NS2Pad];
  Double_t R_s2_la_c[NS2Pad];
  Double_t R_s2_ra_p[NS2Pad];
  Double_t R_s2_ra_c[NS2Pad];

  
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
  T->SetBranchStatus("L.s2.la_p",1);
  T->SetBranchStatus("L.s2.ra_p",1);
  T->SetBranchStatus("L.s2.la_c",1);
  T->SetBranchStatus("L.s2.ra_c",1);
  T->SetBranchStatus("R.s2.nthit",1);
  T->SetBranchStatus("R.s2.t_pads",1);
  T->SetBranchStatus("R.s2.lt",1);
  T->SetBranchStatus("R.s2.rt",1);
  T->SetBranchStatus("R.s2.lt_c",1);
  T->SetBranchStatus("R.s2.rt_c",1);
  T->SetBranchStatus("R.s2.la_p",1);
  T->SetBranchStatus("R.s2.ra_p",1);
  T->SetBranchStatus("R.s2.la_c",1);
  T->SetBranchStatus("R.s2.ra_c",1);

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
  T->SetBranchAddress("L.s2.nthit",&L_s2_nthit);
  T->SetBranchAddress("L.s2.t_pads",L_s2_t_pads);
  T->SetBranchAddress("L.s2.la_p",L_s2_la_p);
  T->SetBranchAddress("L.s2.ra_p",L_s2_ra_p);
  T->SetBranchAddress("L.s2.la_c",L_s2_la_c);
  T->SetBranchAddress("L.s2.ra_c",L_s2_ra_c);
  
  T->SetBranchAddress("R.s2.lt",R_s2_lt);
  T->SetBranchAddress("R.s2.rt",R_s2_rt);
  T->SetBranchAddress("R.s2.lt_c",R_s2_lt_c);
  T->SetBranchAddress("R.s2.rt_c",R_s2_rt_c);
  T->SetBranchAddress("R.s2.nthit",&R_s2_nthit);
  T->SetBranchAddress("R.s2.t_pads",R_s2_t_pads);
  T->SetBranchAddress("R.s2.la_p",R_s2_la_p);
  T->SetBranchAddress("R.s2.ra_p",R_s2_ra_p);
  T->SetBranchAddress("R.s2.la_c",R_s2_la_c);
  T->SetBranchAddress("R.s2.ra_c",R_s2_ra_c);
  
  T->SetBranchAddress("DR.rrawt2", Trig_time_T2); 
  T->SetBranchAddress("Ndata.DR.rrawt2", &Trig_No_T2);
  T->SetBranchAddress("DR.rrawt5",Trig_time_T5);
  T->SetBranchAddress("Ndata.DR.rrawt5",&Trig_No_T5);
  

  Int_t nentries = T->GetEntries();
  //  nentries = 10e4;
  //  Int_t nentries = 100;

  // coincidence histogram widths
  //  Int_t coinc_bins = 140;
  Int_t coinc_bins = 300;
  Double_t coinc_start = 0;
  Double_t coinc_end = 2e-7;

  
  // intitialise timing coincidence histogram
  //  TH1F *h1 = new TH1F("h1"," S2-Time difference (corrected)",coinc_bins,coinc_start-coin_off,coinc_end-coin_off);
  TH1F *h1 = new TH1F("h1"," S2-Time difference (corrected)",coinc_bins,coinc_start,coinc_end);    

  // intitialise timing coincidence histogram
  TH1F *h2 = new TH1F("h2"," S2-Time difference (uncorrected)",coinc_bins,coinc_start,coinc_end);

  // intitialise timing coincidence histogram with offset + pl corrections
  TH1F *h3 = new TH1F("h3"," S2-Time difference (offset + pl corrected)",coinc_bins,coinc_start,coinc_end);

  // intitialise timing coincidence histogram with offset + TW (1P) corrections
  TH1F *h4 = new TH1F("h4"," S2-Time difference (offset, pl + TW (1P) corrected)",coinc_bins,coinc_start,coinc_end);

  // intitialise timing coincidence histogram with offset + TW (3P) corrections
  TH1F *h5 = new TH1F("h5"," S2-Time difference (offset, pl + TW (3P) corrected)",coinc_bins,coinc_start,coinc_end);    

  

  
  // itnitialise paddle versus coincidence histograms
  TH2F* h_L_lvT = new TH2F("h_L_lvT","S2: Coincidence time versus l-paddle (LHRS corrected)",NS2Pad,1,NS2Pad,coinc_bins,coinc_start,coinc_end);
  h_L_lvT->GetXaxis()->SetTitle("LHRS S2 l-paddle #");
  h_L_lvT->GetYaxis()->SetTitle("Coincidence time (s)");
  

  TH2F* h_R_lvT = new TH2F("h_R_lvT","S2: Coincidence time versus l-paddle (RHRS corrected)",NS2Pad,1,NS2Pad,coinc_bins,coinc_start,coinc_end);
  h_R_lvT->GetXaxis()->SetTitle("RHRS S2 l-paddle #");
  h_R_lvT->GetYaxis()->SetTitle("Coincidence time (s)");
      

  TH2F* h_L_lvT_un = new TH2F("h_L_lvT_un","S2: Coincidence time versus l-paddle (LHRS uncorrected)",NS2Pad,1,NS2Pad,coinc_bins,coinc_start,coinc_end);
  h_L_lvT_un->GetXaxis()->SetTitle("LHRS S2 l-paddle #");
  h_L_lvT_un->GetYaxis()->SetTitle("Coincidence time (s)");

  TH2F* h_R_lvT_un = new TH2F("h_R_lvT_un","S2: Coincidence time versus l-paddle (RHRS uncorrected)",NS2Pad,1,NS2Pad,coinc_bins,coinc_start,coinc_end);
  h_R_lvT_un->GetXaxis()->SetTitle("RHRS S2 l-paddle #");
  double t = 0.0;
  h_R_lvT_un->GetYaxis()->SetTitle("Coincidence time (s)");

 
  TH2F* h_L_lvT_pl = new TH2F("h_L_lvT_pl","S2: Coincidence time versus l-paddle (LHRS pl-corrected)",NS2Pad,1,NS2Pad,coinc_bins,coinc_start,coinc_end);
  h_L_lvT_pl->GetXaxis()->SetTitle("LHRS S2 l-paddle #");
  h_L_lvT_pl->GetYaxis()->SetTitle("Coincidence time (s)");
  

  TH2F* h_R_lvT_pl = new TH2F("h_R_lvT_pl","S2: Coincidence time versus l-paddle (RHRS pl-corrected)",NS2Pad,1,NS2Pad,coinc_bins,coinc_start,coinc_end);
  h_R_lvT_pl->GetXaxis()->SetTitle("RHRS S2 l-paddle #");
  h_R_lvT_pl->GetYaxis()->SetTitle("Coincidence time (s)");

   
  TH2F* h_L_lvT_tw_1p = new TH2F("h_L_lvT_tw_1p","S2: Coincidence time versus l-paddle (LHRS pl + Tw (1P)-corrected)",NS2Pad,1,NS2Pad,coinc_bins,coinc_start,coinc_end);
  h_L_lvT_tw_1p->GetXaxis()->SetTitle("LHRS S2 l-paddle #");
  h_L_lvT_tw_1p->GetYaxis()->SetTitle("Coincidence time (s)");
  

  TH2F* h_R_lvT_tw_1p = new TH2F("h_R_lvT_tw_1p","S2: Coincidence time versus l-paddle (RHRS pl + Tw (1P)-corrected) ",NS2Pad,1,NS2Pad,coinc_bins,coinc_start,coinc_end);
  h_R_lvT_tw_1p->GetXaxis()->SetTitle("RHRS S2 l-paddle #");
  h_R_lvT_tw_1p->GetYaxis()->SetTitle("Coincidence time (s)");

  TH2F* h_L_lvT_tw_3p = new TH2F("h_L_lvT_tw_3p","S2: Coincidence time versus l-paddle (LHRS pl + Tw (3P)-corrected)",NS2Pad,1,NS2Pad,coinc_bins,coinc_start,coinc_end);
  h_L_lvT_tw_3p->GetXaxis()->SetTitle("LHRS S2 l-paddle #");
  h_L_lvT_tw_3p->GetYaxis()->SetTitle("Coincidence time (s)");
  

  TH2F* h_R_lvT_tw_3p = new TH2F("h_R_lvT_tw_3p","S2: Coincidence time versus l-paddle (RHRS pl + Tw (3P)-corrected) ",NS2Pad,1,NS2Pad,coinc_bins,coinc_start,coinc_end);
  h_R_lvT_tw_3p->GetXaxis()->SetTitle("RHRS S2 l-paddle #");
  h_R_lvT_tw_3p->GetYaxis()->SetTitle("Coincidence time (s)");
 
  
  

  
  Double_t LTime = 0.0;
  Double_t RTime = 0.0;

  Double_t LTime_un = 0.0;
  Double_t RTime_un = 0.0;

  Double_t LTime_pl = 0.0;
  Double_t RTime_pl = 0.0;
  Double_t L_pl_corr = 0.0;
  Double_t R_pl_corr = 0.0;


  Double_t L_tw_b = 0.0; // difference between left and right paddle ADCs TW effects
  
  Double_t LTime_tw_1p = 0.0;
  Double_t L_tw_1p_corr = 0.0;

  Double_t LTime_tw_3p = 0.0;
  Double_t L_tw_3p_corr = 0.0;

  Double_t R_tw_b = 0.0; // difference between left and right paddle ADCs TW effects

  Double_t RTime_tw_1p = 0.0;
  Double_t R_tw_1p_corr = 0.0;

  Double_t RTime_tw_3p = 0.0;
  Double_t R_tw_3p_corr = 0.0;

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

    LTime_pl = 0.0;
    RTime_pl = 0.0;

    L_pl_corr = 0.0;
    R_pl_corr = 0.0;

    L_tw_b = 0.0;
    
    L_tw_1p_corr = 0.0;
    LTime_tw_1p = 0.0;

    R_tw_b = 0.0;
    
    R_tw_1p_corr = 0.0;
    RTime_tw_1p = 0.0;
    
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
	
	if (L_s2_t_pads[0]==j && (L_s2_lt_c[j]+L_s2_rt_c[j]) != 0){ //&& L_s2_la_c[j]>50 && L_s2_ra_c[j]>50){

	  
	  //LTime = (L_s2_lt_c[j]+L_s2_rt_c[j])/2.;
	  //	  LTime = fTdc2T*(L_s2_lt[j]+L_ls2_coeff[j] + L_s2_rt[j]+L_rs2_coeff[j])/2.;
	  LTime = fTdc2T*(L_s2_lt[j]-L_ls2_coeff[j] + L_s2_rt[j]-L_rs2_coeff[j])/2.;
	  LTime_un = (fTdc2T*(L_s2_lt[j]+L_s2_rt[j]))/2.;

          L_pl_corr = L_th_slope * L_r_th[0] + L_ph_slope * L_r_ph[0] + L_x_slope * L_r_x[0];
	  LTime_pl = (fTdc2T*(L_s2_lt[j]-L_ls2_coeff[j]+L_s2_rt[j]-L_rs2_coeff[j]))/2. - L_pl_corr;


          if(j==4 || j==5 || L_s2_la_c[j]<100 || L_s2_ra_c[j]<100){
            L_tw_b = 0.0;
          }
          else{
	    L_tw_b = ((1/TMath::Sqrt(L_s2_la_c[j])) - (1/TMath::Sqrt(L_s2_ra_c[j])));
          }

	  L_tw_1p_corr = L_s2_Tw_1P_slope[j]*L_tw_b;
	  //	  L_tw_1p_corr = L_s2_Tw_1P_slope[j]*(L_tw_b +  L_s2_Tw_1P_off[j]);	  
	  L_tw_1p_corr = 0.0;

	  
          LTime_tw_1p = (fTdc2T*(L_s2_lt[j]-L_ls2_coeff[j]+L_s2_rt[j]-L_rs2_coeff[j] - L_tw_1p_corr ))/2. - L_pl_corr;


	  LHRS_pad = j+1;
	
	}

      }


      for(int j=0;j<NS2Pad;j++){
	
	if (R_s2_t_pads[0]==j && (R_s2_lt_c[j]+R_s2_rt_c[j]) != 0 && R_s2_la_c[j]>50 && R_s2_ra_c[j]>50){
	  //	  	if (R_s2_t_pads[0]==j && (R_s2_lt_c[j]+R_s2_rt_c[j]) != 0){// && R_s2_la_c[j]>50 && R_s2_ra_c[j]>50){

	  
	  //	  RTime = (R_s2_lt_c[j]+R_s2_rt_c[j])/2.;
	  //	  RTime = fTdc2T*(R_s2_lt[j]+R_ls2_coeff[j] + R_s2_rt[j]+R_rs2_coeff[j])/2.;
	  RTime = fTdc2T*(R_s2_lt[j]-R_ls2_coeff[j] + R_s2_rt[j]-R_rs2_coeff[j])/2.;
	  RTime_un = (fTdc2T*(R_s2_lt[j]+R_s2_rt[j]))/2.;
	  R_pl_corr = R_th_slope * R_r_th[0] + R_ph_slope * R_r_ph[0] + R_x_slope * R_r_x[0];
	  RTime_pl = (fTdc2T*(R_s2_lt[j]-R_ls2_coeff[j] + R_s2_rt[j]-R_rs2_coeff[j]))/2. - R_pl_corr;

	  if(R_s2_la_c[j]<100 || R_s2_ra_c[j]<100 || !T5Pass || !T2Pass){
	    R_tw_b = 0.0;
	    R_tw_1p_corr  = 0.0;
	  }
	  else{
	    R_tw_b = ((1/TMath::Sqrt(R_s2_la_c[j])) - (1/TMath::Sqrt(R_s2_ra_c[j])));
	    R_tw_1p_corr = R_s2_Tw_1P_slope[j]*(R_tw_b -  R_s2_Tw_1P_off[j]);
	    //	  R_tw_1p_corr = R_s2_Tw_1P_slope[j]*R_tw_b;	    
	  }


	  //	  R_tw_1p_corr = 0.0;

	  
	  
	  
          RTime_tw_1p = (fTdc2T*(R_s2_lt[j]-R_ls2_coeff[j] - R_tw_1p_corr +R_s2_rt[j]-R_rs2_coeff[j]))/2. - R_pl_corr;



	  RHRS_pad = j+1;

	  if(i%10000 == 0){
	    cout << "Event " << i << ": R_pl_corr = " << R_pl_corr << ", RTime_pl = " << RTime_pl << ", RTime = " << RTime << endl;
	    cout << "R_r_th = " << R_r_th[0] << ", R_r_ph = " << R_r_ph[0] << ", R_r_x = " << R_r_x[0] << endl;
	    cout << endl;
	  }

	}

      }


      if( (LTime-RTime) != 0 ){
      
	h1->Fill(LTime-RTime);
	h2->Fill(LTime_un-RTime_un);
	h3->Fill(LTime_pl-RTime_pl);
	h4->Fill(LTime_pl - RTime_tw_1p);
      
	h_L_lvT->Fill(LHRS_pad,LTime-RTime);
	h_R_lvT->Fill(RHRS_pad,LTime-RTime);

	h_L_lvT_un->Fill(LHRS_pad,LTime_un-RTime_un);
	h_R_lvT_un->Fill(RHRS_pad,LTime_un-RTime_un);

	h_L_lvT_pl->Fill(LHRS_pad,LTime_pl-RTime_pl);
	h_R_lvT_pl->Fill(RHRS_pad,LTime_pl-RTime_pl);

	h_L_lvT_tw_1p->Fill(LHRS_pad,LTime_tw_1p-RTime_tw_1p);
	h_R_lvT_tw_1p->Fill(RHRS_pad,LTime_tw_1p-RTime_tw_1p);
      }
      
    }
    

      
  }
  


  // fit peak of corrected coincidence timing

  TCanvas* c1 = new TCanvas("c1","c1",1000,800);
  
  h1->Draw();
  cout << "h1 entries = " << h1->GetEntries() << endl;;
  
  Double_t max= h1->GetBinCenter(h1->GetMaximumBin());


  
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


  TCanvas* c3 = new TCanvas("c3","c3",1000,800);
  
  h3->Draw();
  
  Double_t max3= h3->GetBinCenter(h3->GetMaximumBin());


  
  TF1* f1_pl = new TF1("f1_pl","pol1",max3-(6e-8),max3-(1e-8));
  h3->Fit(f1_pl,"FRQ");
  
  Double_t par3[5];   
  f1_pl->GetParameters(&par3[0]);
  
  h3->Fit(f1_pl,"FRQ");
  
  f1_pl->GetParameters(&par3[0]);
  
  TF1* f2_pl = new TF1("f2_pl","gaus",max3-(2.5e-9),max3+(2.5e-9));
 
  h3->Fit(f2_pl,"FRQ");
 
  
  f2_pl->GetParameters(&par3[2]);
 
 

  TF1* f3_pl = new TF1("f3_pl","pol1(0)+gaus(2)",max3-(6e-8),max3+(1e-8));
 
  f3_pl->SetParameters(par3);
  h3->Fit(f3_pl,"FRQ");

  f3_pl->GetParameters(&par3[0]);
  
  Double_t plcorr_width = par3[4];
  Double_t plcorr_mean = par3[3];


  gPad->Update();

  
  TPaveStats *st_3 = (TPaveStats*)h3->FindObject("stats");
  st_3->SetX1NDC(0.15); 
  st_3->SetX2NDC(0.45); 
  st_3->SetY1NDC(0.45); 
  st_3->SetY2NDC(0.90); 

  
  TCanvas* c4 = new TCanvas("c4","c4",1000,800);
  
  h4->Draw();
  
  Double_t max4= h4->GetBinCenter(h4->GetMaximumBin());


  
  TF1* f1_tw_p1 = new TF1("f1_tw_p1","pol1",max4-(6e-8),max4-(1e-8));
  h4->Fit(f1_tw_p1,"FRQ");
  
  Double_t par4[5];   
  f1_tw_p1->GetParameters(&par4[0]);
  
  h4->Fit(f1_tw_p1,"FRQ");
  
  f1_tw_p1->GetParameters(&par4[0]);
  
  TF1* f2_tw_p1 = new TF1("f2_tw_p1","gaus",max4-(2.5e-9),max4+(2.5e-9));
 
  h4->Fit(f2_tw_p1,"FRQ");
 
  
  f2_tw_p1->GetParameters(&par4[2]); 
 

  TF1* f3_tw_p1 = new TF1("f3_tw_p1","pol1(0)+gaus(2)",max4-(6e-8),max4+(1e-8));
 
  f3_tw_p1->SetParameters(par4);
  h4->Fit(f3_tw_p1,"FRQ");

  f3_tw_p1->GetParameters(&par4[0]);
  
  Double_t tw_1p_width = par4[4];
  Double_t tw_1p_mean = par4[3];


  gPad->Update();

  
  TPaveStats *st_4 = (TPaveStats*)h4->FindObject("stats");
  st_4->SetX1NDC(0.15); 
  st_4->SetX2NDC(0.45); 
  st_4->SetY1NDC(0.45); 
  st_4->SetY2NDC(0.90); 




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

  h_L_lvT_pl->Draw("colz");


  TCanvas* c10 = new TCanvas("c10","c10",1000,800);

  h_R_lvT_pl->Draw("colz");


  TCanvas* c11 = new TCanvas("c11","c11",1000,800);

  h_L_lvT_tw_1p->Draw("colz");


  TCanvas* c12 = new TCanvas("c12","c12",1000,800);

  h_R_lvT_tw_1p->Draw("colz");




  // print canvases

  c1->Print(Form("plots/coinc_time/%s.pdf(",Name.Data()));
  c2->Print(Form("plots/coinc_time/%s.pdf",Name.Data()));
  c3->Print(Form("plots/coinc_time/%s.pdf",Name.Data()));
  c4->Print(Form("plots/coinc_time/%s.pdf",Name.Data()));
  c5->Print(Form("plots/coinc_time/%s.pdf",Name.Data()));
  c6->Print(Form("plots/coinc_time/%s.pdf",Name.Data()));
  c7->Print(Form("plots/coinc_time/%s.pdf",Name.Data()));
  c8->Print(Form("plots/coinc_time/%s.pdf",Name.Data()));
  c9->Print(Form("plots/coinc_time/%s.pdf",Name.Data()));
  c10->Print(Form("plots/coinc_time/%s.pdf",Name.Data()));
  c11->Print(Form("plots/coinc_time/%s.pdf",Name.Data()));
  c12->Print(Form("plots/coinc_time/%s.pdf)",Name.Data()));
  
  //  plots/coinc_time



  // write information about fit to csv file (and print to screen)

  if(uncor_flag == 1){

    csvFile << runno << "," << Name << "," << corr_width << "," << corr_mean << endl;
    
    csvFile.close();
    
    cout << "Runnumber = " << runno << ", Description: " << Name << " , width = " << corr_width << ", corr_mean = " << corr_mean << endl;
  }
  else{
    cout << "WARNING: uncor_flag is set to true, printing and saving results to csv file for uncorrected coincidence peak width" << endl;

    csvFile << runno << "," << Name << "," << uncorr_width << "," << uncorr_mean << endl;
    
    csvFile.close();
    
    cout << "Runnumber = " << runno << ", Description: " << Name << " , width = " << uncorr_width << ", uncorr_mean = " << uncorr_mean << endl;

  }
  

}
