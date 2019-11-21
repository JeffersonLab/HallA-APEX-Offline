#include <iostream>
#include <cassert>

#include "TROOT.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TString.h"

using namespace std;

class ROpticsOpt;

ROpticsOpt * opt;

UInt_t NPara = 0;
Double_t OldMatrixArray[10000] = {-99}; //NPara
Bool_t freed[10000] = {kFALSE}; //NPara

UInt_t MaxDataPerGroup = 10000000;

void replay(TString OutputFile, TString DataBase){

  //Simplified replay code. Takes some functions for ROpticsOpt class and uses them to calculate new target variables with optimized matrix.
  //Currently relies on old data being written in an ascii text file, should probably simplify that step later

  TString rootfiles = "/home/sean/Grad/Research/APEX/Rootfiles/";
  TFile* f_old = new TFile(rootfiles + "apex_4647.root","open");
  TTree* t;
  f_old->GetObject("T",t);

  int entries = t->GetEntries();

  TFile* f_new = new TFile(rootfiles + OutputFile,"recreate");
  TTree* t_new = new TTree("T","");
  
  

  opt = new ROpticsOpt();
  //Need to create a .dat file with all of the old root file information before running this. Will fix later
  opt->LoadRawData(rootfiles + "apex_4647.dat", (UInt_t) - 1, MaxDataPerGroup);
  opt->LoadDataBase(DataBase);
  //opt->PrepareSieve();

  double R_tr_n;
  double R_tr_x_fp[100];
  double R_tr_y_fp[100];
  double R_tr_th_fp[100];
  double R_tr_ph_fp[100];
  double R_tr_p[100];
  double R_cer_asum_c;
  double R_s0_nthit;
  double beam_x[100];
  double beam_y[100];
  double R_tr_x_rot[100];
  double R_tr_y_rot[100];
  double R_tr_th_rot[100];
  double R_tr_ph_rot[100];
  double R_tr_vz[100];
  double R_tr_tg_y[100];
  double R_tr_tg_th[100];
  double R_tr_tg_ph[100];
  double R_tr_tg_dp[100];
  double sieve_x[100];
  double sieve_y[100];
  
  t->SetBranchStatus("*",0);
  t->SetBranchStatus("R.tr.n",1);
  t->SetBranchStatus("R.tr.x",1);
  t->SetBranchStatus("R.tr.y",1);
  t->SetBranchStatus("R.tr.th",1);
  t->SetBranchStatus("R.tr.ph",1);
  t->SetBranchStatus("R.tr.p",1);
  t->SetBranchStatus("R.cer.asum_c",1);
  t->SetBranchStatus("R.s0.nthit",1);
  t->SetBranchStatus("Rrb.x",1);
  t->SetBranchStatus("Rrb.y",1);
  t->SetBranchStatus("R.tr.r_x",1);
  t->SetBranchStatus("R.tr.r_y",1);
  t->SetBranchStatus("R.tr.r_th",1);
  t->SetBranchStatus("R.tr.r_ph",1);
  t->SetBranchStatus("R.tr.vz",1);
  
  t->SetBranchAddress("R.tr.n",&R_tr_n);
  t->SetBranchAddress("R.tr.x",R_tr_x_fp);
  t->SetBranchAddress("R.tr.y",R_tr_y_fp);
  t->SetBranchAddress("R.tr.th",R_tr_th_fp);
  t->SetBranchAddress("R.tr.ph",R_tr_ph_fp);
  t->SetBranchAddress("R.tr.p",R_tr_p);
  t->SetBranchAddress("R.cer.asum_c",&R_cer_asum_c);
  t->SetBranchAddress("R.s0.nthit",&R_s0_nthit);
  t->SetBranchAddress("Rrb.x",beam_x);
  t->SetBranchAddress("Rrb.y",beam_y);
  t->SetBranchAddress("R.tr.r_x",R_tr_x_rot);
  t->SetBranchAddress("R.tr.r_y",R_tr_y_rot);
  t->SetBranchAddress("R.tr.r_th",R_tr_th_rot);
  t->SetBranchAddress("R.tr.r_ph",R_tr_ph_rot);
  t->SetBranchAddress("R.tr.vz",R_tr_vz);

  t_new->Branch("R.tr.n",&R_tr_n);
  t_new->Branch("R.tr.x",R_tr_x_fp);
  t_new->Branch("R.tr.y",R_tr_y_fp);
  t_new->Branch("R.tr.th",R_tr_th_fp);
  t_new->Branch("R.tr.ph",R_tr_ph_fp);
  t_new->Branch("R.tr.p",R_tr_p);
  t_new->Branch("R.cer.asum_c",&R_cer_asum_c);
  t_new->Branch("R.s0.nthit",&R_s0_nthit);
  t_new->Branch("Rrb.x",beam_x);
  t_new->Branch("Rrb.y",beam_y);
  t_new->Branch("R.tr.r_x",R_tr_x_rot);
  t_new->Branch("R.tr.r_y",R_tr_y_rot);
  t_new->Branch("R.tr.r_th",R_tr_th_rot);
  t_new->Branch("R.tr.r_ph",R_tr_ph_rot);
  t_new->Branch("R.tr.vz",R_tr_vz);
  t_new->Branch("R.tr.tg_y",R_tr_tg_y);
  t_new->Branch("R.tr.tg_th",R_tr_tg_th);
  t_new->Branch("R.tr.tg_ph",R_tr_tg_ph);
  t_new->Branch("R.tr.tg_dp",R_tr_tg_dp);
  t_new->Branch("Sieve.x",sieve_x);
  t_new->Branch("Sieve.y",sieve_y);
    
  for(int i=0; i<entries; i++){
  
    t->GetEntry(i);

    R_tr_tg_y[0] = opt->calc_tgy(i);
    R_tr_tg_th[0] = opt->calc_tgth(i);
    R_tr_tg_ph[0] = opt->calc_tgph(i);
    R_tr_tg_dp[0] = opt->calc_tgdp(i);
    sieve_x[0] = opt->sieve_x(i);
    sieve_y[0] = opt->sieve_y(i);

    
    t_new->Fill();
    
  }

  f_new->Write();
  
  

}

