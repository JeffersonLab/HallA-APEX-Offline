/*
*************************************************************
11/12/20 John Williamson
Script that performs pedestal and gain calibration for s0 for the LHRS.

Adapted from tritium experiment scripts.

*************************************************************
*/


#include "file_def.h"
#include "Load_more_rootfiles.C"


void Rs0_cali(Int_t runno)
{

  
  TChain* T = Load_more_rootfiles(runno);
  

  
  // choose value for gains to align 1 PE peak to
  Double_t PE_peak = 400.0;


  
  cout << T->GetEntries() << endl;
  
  Double_t min, max, temp;
  Double_t ped_val,ped_wid;
  Int_t max_bin;
  Double_t lfgain;
  ofstream myfile;
  myfile.open(Form("DB/R_s0_ped_%i.txt",runno));//,fstream::app);

  TCut tritype = "";//"((DR.evtypebits>>1)&1)";
 
  TCanvas *c1 = new TCanvas("c1","c1",1200,1200);
  c1->Divide(2,1);
  c1->cd(1);
  TH1F *tt1 = new TH1F("tt1","FADC R.s0.la",500,8000,14000);
  T->Draw("R.s0.la>>tt1",tritype);
  tt1->SetXTitle("R.s0.la");
  
  max_bin = tt1->GetMaximumBin();
  min = tt1->GetBinCenter(max_bin) - 30;
  max = tt1->GetBinCenter(max_bin) + 30;
  tt1->Fit("gaus","Q","",min,max);
  TF1* gaus1 = tt1->GetFunction("gaus");
  min = gaus1->GetParameter(1) - 1.5*gaus1->GetParameter(2);
  max = gaus1->GetParameter(1) + 1.5*gaus1->GetParameter(2);

  tt1->Fit("gaus","Q","",min,max);
  gaus1 = tt1->GetFunction("gaus");
  ped_val = gaus1->GetParameter(1);
  ped_wid = gaus1->GetParameter(2);    

  
  myfile<<"s0la"<<endl;
  myfile<<"ped: "<<ped_val <<endl;
  myfile<<"ped wid: "<< ped_wid <<endl;

  cout << endl;
  cout<<"s0la"<<endl;
  cout<<"ped: "<<ped_val <<endl;
  cout<<"ped wid: "<< ped_wid <<endl;
  


  c1->cd(2);
  TH1F *tt2 = new TH1F("tt2","FADC Rs0.la_p",300,100,2300);
  T->Draw("R.s0.la_p>>tt2",tritype);
  tt2->SetXTitle("R.s0.la_p");
  max_bin = tt2->GetMaximumBin();
  // min = tt2->GetBinCenter(max_bin) - 100;
  // max = tt2->GetBinCenter(max_bin) + 300;
  min = 100;
  max = 1300;
  tt2->Fit("landau","Q","",min,max);
  
  TF1* landau1 = tt2->GetFunction("landau");
  min = landau1->GetParameter(1) - 1.5*landau1->GetParameter(2);
  max = landau1->GetParameter(1) + 1.5*landau1->GetParameter(2);
  tt2->Fit("landau","Q","",min,max);

  lfgain=PE_peak*1.0/landau1->GetParameter(1);
  // this parameter of landau is MPV (Most Probable Value), peak of landau

  
  myfile<<"gain: "<<lfgain<<endl;
  cout<<"gain: "<<lfgain<<endl<<endl;
 


  TCanvas *c2 = new TCanvas("c2","c2",1200,1200);
  c2->Divide(2,1);
  c2->cd(1);
  TH1F *tt4 = new TH1F("tt4","FADC R.s0.ra",500,8000,14000);
  T->Draw("R.s0.ra>>tt4",tritype);
  tt4->SetXTitle("R.s0.ra");
  max_bin = tt4->GetMaximumBin();
  min = tt4->GetBinCenter(max_bin) - 30;
  max = tt4->GetBinCenter(max_bin) + 30;
  tt4->Fit("gaus","Q","",min,max);
  gaus1 = tt4->GetFunction("gaus");
  min = gaus1->GetParameter(1) - 1*gaus1->GetParameter(2);
  max = gaus1->GetParameter(1) + 1*gaus1->GetParameter(2);

  tt4->Fit("gaus","Q","",min,max);
  gaus1 = tt4->GetFunction("gaus");

  ped_val = gaus1->GetParameter(1);
  ped_wid = gaus1->GetParameter(2);


  myfile<<"s0ra"<<endl;
  myfile<<"ped: "<<ped_val <<endl;
  myfile<<"ped wid: "<< ped_wid <<endl;

  cout<<"s0ra"<<endl;
  cout<<"ped: "<<ped_val <<endl;
  cout<<"ped wid: "<< ped_wid <<endl;

  c2->cd(2);
  TH1F *tt5 = new TH1F("tt5","FADC R.s0.ra_p",300,100,1300);
  T->Draw("R.s0.ra_p>>tt5",tritype);
  tt5->SetXTitle("R.s0.ra_p");
  max_bin = tt5->GetMaximumBin();
  // min = tt5->GetBinCenter(max_bin) - 100;
  // max = tt5->GetBinCenter(max_bin) + 300;
  min = 100;
  max = 1300;
  tt5->Fit("landau","Q","",min,max);
  landau1 = tt5->GetFunction("landau");
  min = landau1->GetParameter(1) - 1.5*landau1->GetParameter(2);
  max = landau1->GetParameter(1) + 1.5*landau1->GetParameter(2);
  tt5->Fit("landau","Q","",min,max);
  landau1 = tt5->GetFunction("landau");

  lfgain=PE_peak*1.0/landau1->GetParameter(1);
  // this parameter of landau is MPV (Most Probable Value), peak of landau


  myfile<<"gain: "<<lfgain<<endl;
  cout<<"gain: "<<lfgain<<endl;

  

  c1->Print(Form("plots/Rs0_%i.pdf[",runno));
  c1->Print(Form("plots/Rs0_%i.pdf",runno));
  c2->Print(Form("plots/Rs0_%i.pdf",runno));
  c2->Print(Form("plots/Rs0_%i.pdf]",runno));

  myfile.close();

}

