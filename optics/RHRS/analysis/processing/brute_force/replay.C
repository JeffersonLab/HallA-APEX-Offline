#include <iostream>
#include <cassert>

#include "TROOT.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TString.h"

using namespace std;

class ROpticsOpt;

ROpticsOpt * opt;

void replay(TString run,TString order){

  //Simplified replay code. Takes some functions for ROpticsOpt class and uses them to calculate new target variables with different matrix.

  //REMEMBER to change to unraster beam in ROpticsOpt.C script for certain runs

  TString rootfiles = "/home/sean/Grad/Research/APEX/Rootfiles/";

  TChain *t = new TChain("T");
  
  t->Add(rootfiles + "apex_" + run + ".root");

 
  int entries = t->GetEntries();
  
  TFile* f_new = new TFile(rootfiles + "apex_" + run + "_opt_"+order+"_xfp_full_V_wires.root","recreate");
  TTree* t_new = new TTree("T","");

  run = "V_wires";       //Set equal to run if not using multiple wires
  //Load focal plane data with new matrix
  opt = new ROpticsOpt();
  
  opt->LoadRawData(t);
  opt->LoadDataBase("../DB/" + run + "/db_R.vdc.dat_y_"+order+"_xfp","-50_-30");
  opt->LoadDataBase("../DB/" + run + "/db_R.vdc.dat_y_"+order+"_xfp","-30_-10");
  opt->LoadDataBase("../DB/" + run + "/db_R.vdc.dat_y_"+order+"_xfp","-10_10");
  opt->LoadDataBase("../DB/" + run + "/db_R.vdc.dat_y_"+order+"_xfp","10_30");
  opt->LoadDataBase("../DB/" + run + "/db_R.vdc.dat_y_"+order+"_xfp","30_50");
   
  //Define all the variables we want our output tree to have
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
  double ubeam_x[100];
  double ubeam_y[100];
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
  double BPMA_x[100];
  double BPMA_y[100];
  double BPMB_x[100];
  double BPMB_y[100];
  
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
  t->SetBranchStatus("Rurb.y",1);
  t->SetBranchStatus("Rurb.y",1);
  t->SetBranchStatus("R.tr.r_x",1);
  t->SetBranchStatus("R.tr.r_y",1);
  t->SetBranchStatus("R.tr.r_th",1);
  t->SetBranchStatus("R.tr.r_ph",1);
  t->SetBranchStatus("R.tr.vz",1);
  t->SetBranchStatus("Rrb.BPMA.x",1);
  t->SetBranchStatus("Rrb.BPMA.y",1);
  t->SetBranchStatus("Rrb.BPMB.x",1);
  t->SetBranchStatus("Rrb.BPMB.y",1);
  
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
  t->SetBranchAddress("Rrb.BPMA.x",BPMA_x);
  t->SetBranchAddress("Rrb.BPMA.y",BPMA_y);
  t->SetBranchAddress("Rrb.BPMB.x",BPMB_x);
  t->SetBranchAddress("Rrb.BPMB.y",BPMB_y);

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
  t_new->Branch("Rurb.x",ubeam_x);
  t_new->Branch("Rurb.y",ubeam_y);
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


  //Calculate target variables
  for(int i=0; i<entries; i++){
  
    t->GetEntry(i);
    
    R_tr_tg_y[0] = opt->calc_tgy(i);
    R_tr_tg_th[0] = opt->calc_tgth(i);
    R_tr_tg_ph[0] = opt->calc_tgph(i);
    R_tr_tg_dp[0] = opt->calc_tgdp(i);
    R_tr_vz[0] = opt->calc_vz(i, R_tr_tg_y[0]);
    sieve_x[0] = opt->sieve_x(i);
    sieve_y[0] = opt->sieve_y(i);
    
    t_new->Fill();

    if(i%10000 == 0) cout<<std::setprecision(2)<<i*1.0/entries*100<<"%"<<endl;
  }

  f_new->Write();

}

