/*
*************************************************************
11/12/20 John Williamson
Script that performs pedestal and gain calibration for s2 for the LHRS.

Adapted from tritium experiment scripts.

*************************************************************
*/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <stdexcept>
#include <cassert>

#include "TTree.h"
#include "TFile.h"
#include "TString.h"

#include "file_def.h"
#include "Load_more_rootfiles.C"

using namespace std;
void Lscin_cali(Int_t runno)
{


  TChain* T = Load_more_rootfiles(runno);

  cout << T->GetEntries() << endl;


  const Int_t NS2Pad = 16;
  
  Double_t min, max, temp;
  Double_t ped_val,ped_wid;
  Int_t max_bin;
  Int_t i, ii;
  Double_t Lped[NS2Pad],Lped_width[NS2Pad];

  Double_t Rped[NS2Pad],Rped_width[NS2Pad];
  
  ofstream myfile;
  myfile.open(Form("DB/Lscin_ped_%i.txt",runno));//,fstream::app);

  // myfile << "Values for Left paddles" << endl << endl;
  
  // myfile << setiosflags(ios::left) << setw(2) << "#" << "   "; 
  // myfile << setiosflags(ios::left) << setw(5) << "Ped" << "   ";
  // myfile << setiosflags(ios::left) << setw(8) << "ped width" << "   "<<endl;


  TCut tritype = "";//"((DR.evtypebits>>1)&1)";


  TF1* gaus_la[NS2Pad]; // gaussian function for l ADC fitting
  TF1* gaus_ra[NS2Pad]; // gaussian function for r ADC fitting

  TF1* landau_la[NS2Pad]; // landau function for l ADC (ped subtracted) fitting
  TF1* landau_ra[NS2Pad]; // landau function for r ADC (ped subtracted) fitting
  
 
  TCanvas *c1 = new TCanvas("c1","c1",1200,1200);
  c1->Divide(4,4); 
  for(i=0;i<NS2Pad;i++){
    c1->cd(i+1);
    TH1F *tt1 = new TH1F("tt1",Form("FADC_Ls2_la%d",i+1),500,4000,8000);
    T->Draw(Form("L.s2.la[%d]>>tt1",i),tritype);
    tt1->SetXTitle(Form("L.s2.la[%d]",i));
    max_bin = tt1->GetMaximumBin();
    min = tt1->GetBinCenter(max_bin) - 40;
    max = tt1->GetBinCenter(max_bin) + 40;
    tt1->Fit("gaus","Q","",min,max);
    gaus_la[i] = tt1->GetFunction("gaus");
    min = gaus_la[i]->GetParameter(1) - 2*gaus_la[i]->GetParameter(2);
    max = gaus_la[i]->GetParameter(1) + 2*gaus_la[i]->GetParameter(2);

    tt1->Fit("gaus","Q","",min,max);
    gaus_la[i] = tt1->GetFunction("gaus");

    ped_val = gaus_la[i]->GetParameter(1);
    ped_wid = gaus_la[i]->GetParameter(2);    

    gPad->SetLogy();

    Lped[i] = ped_val;
    Lped_width[i] = ped_wid; 


    // myfile << setiosflags(ios::left) << setw(2) << setiosflags(ios::fixed) << setprecision(1) << i << "   ";
    // myfile << setiosflags(ios::left) << setw(5) << setiosflags(ios::fixed) << setprecision(2) << ped_val << "   ";
    // myfile << setiosflags(ios::left) << setw(8) << setiosflags(ios::fixed) << setprecision(1) << ped_wid <<endl;
  } 

  Double_t lfgain[NS2Pad];

  // myfile << setiosflags(ios::left) << setw(2) << "#" << "   ";
  // myfile << setiosflags(ios::left) << setw(5) << "gain" << endl;

  Double_t fscale=150;
  TCanvas *c2 = new TCanvas("c2","c2",1200,1200);
  c2->Divide(4,4);
  for(i=0;i<NS2Pad;i++){
    c2->cd(i+1);
    TH1F *tt2 = new TH1F("tt2",Form("Ls2_lap_%d",i+1),300,-200,1200);
    T->Draw(Form("L.s2.la_p[%d]>>tt2",i),tritype);
    tt2->SetXTitle(Form("L.s2.la_p[%d]",i));
    tt2->GetXaxis()->SetRangeUser(0,190);
    Int_t min_bin = tt2->GetMinimumBin();
    min = tt2->GetXaxis()->GetBinCenter(min_bin);
    tt2->GetXaxis()->SetRangeUser(-200,1200);

    TF1 *f1 = new TF1("f1", "landau",min, min+3.5*fscale);
    tt2->Fit("f1", "Rq");
    
    min = f1->GetParameter(1) - 1.5*f1->GetParameter(2);
    max = f1->GetParameter(1) + 1.5*f1->GetParameter(2);

    tt2->Fit("landau","Q","",min,max);
    landau_la[i] = tt2->GetFunction("landau");
    lfgain[i]=300*1.0/landau_la[i]->GetParameter(1);


    gPad->SetLogy();

    // myfile << setiosflags(ios::left) << setw(2) << setiosflags(ios::fixed) << setprecision(1) << i << "   ";
    // myfile << setiosflags(ios::left) << setw(8) << setiosflags(ios::fixed) << setprecision(5) << lfgain[i] <<endl;
  }


  // myfile << endl << endl;
  // myfile<< "Values for Right paddles" << endl << endl;
  
  TCanvas *c3 = new TCanvas("c3","c3",1200,1200);
  c3->Divide(4,4);
  for(i=0;i<NS2Pad;i++){
    c3->cd(i+1);
    TH1F *tt3 = new TH1F("tt3",Form("FADC_Ls2_ra%d",i+1),500,4000,8000);
    T->Draw(Form("L.s2.ra[%d]>>tt3",i,i),tritype);
    tt3->SetXTitle(Form("L.s2.ra[%d]",i));
    max_bin = tt3->GetMaximumBin();
    min = tt3->GetBinCenter(max_bin) - 40;
    max = tt3->GetBinCenter(max_bin) + 40;
    tt3->Fit("gaus","Q","",min,max);
    gaus_ra[i] = tt3->GetFunction("gaus");
    min = gaus_ra[i]->GetParameter(1) - 2*gaus_ra[i]->GetParameter(2);
    max = gaus_ra[i]->GetParameter(1) + 2*gaus_ra[i]->GetParameter(2);

    tt3->Fit("gaus","Q","",min,max);
    gaus_ra[i] = tt3->GetFunction("gaus");

    ped_val = gaus_ra[i]->GetParameter(1);
    ped_wid = gaus_ra[i]->GetParameter(2);

    gPad->SetLogy();

    Rped[i] = ped_val;
    Rped_width[i] = ped_wid;


    // myfile << setiosflags(ios::left) << setw(2) << setiosflags(ios::fixed) << setprecision(1) << i << "   ";
    // myfile << setiosflags(ios::left) << setw(5) << setiosflags(ios::fixed) << setprecision(2) << ped_val << "   ";
    // myfile << setiosflags(ios::left) << setw(8) << setiosflags(ios::fixed) << setprecision(1) << ped_wid <<endl;
  }

  Double_t rfgain[NS2Pad];

  // myfile << setiosflags(ios::left) << setw(2) << "#" << "   ";
  // myfile << setiosflags(ios::left) << setw(5) << "gain" << endl;

  fscale=150;
  TCanvas *c4 = new TCanvas("c4","c4",1200,1200);
  c4->Divide(4,4);
  for(i=0;i<NS2Pad;i++){
    c4->cd(i+1);
    TH1F *tt4 = new TH1F("tt4",Form("Ls2_rap_%d",i+1),300,-200,1200);
    T->Draw(Form("L.s2.ra_p[%d]>>tt4",i,i),tritype);
    tt4->SetXTitle(Form("L.s2.ra_p[%d]",i));
    tt4->GetXaxis()->SetRangeUser(0,190);
    Int_t min_bin = tt4->GetMinimumBin();
    min = tt4->GetXaxis()->GetBinCenter(min_bin);
    tt4->GetXaxis()->SetRangeUser(-200,1200);

    TF1 *f1 = new TF1("f1", "landau",min, min+3.5*fscale);
    tt4->Fit("f1", "Rq");

    min = f1->GetParameter(1) - 1.5*f1->GetParameter(2);
    max = f1->GetParameter(1) + 1.5*f1->GetParameter(2);

    tt4->Fit("landau","Q","",min,max);
    landau_ra[i] = tt4->GetFunction("landau");
    rfgain[i]=300*1.0/landau_ra[i]->GetParameter(1);


    gPad->SetLogy();


    // myfile << setiosflags(ios::left) << setw(2) << setiosflags(ios::fixed) << setprecision(1) << i << "   ";
    // myfile << setiosflags(ios::left) << setw(8) << setiosflags(ios::fixed) << setprecision(5) << rfgain[i] <<endl;
  }


  //loop to write pedestal and gain values to file (in style of DB format for convenience)

  //LHRS
  myfile << "Values for Left paddles" << endl << endl;

  myfile << "L.s2.L.ped = " ;
  //write pedestal values 
  for(i=0;i<NS2Pad;i++){            
    myfile << Lped[i] << " ";
  }
  // gain values
  myfile << endl;
  myfile << "L.s2.L.gain = " ;
  for(i=0;i<NS2Pad;i++){            
    myfile << lfgain[i] << " ";
  }
  // ped wid values (not in DB, written in DP-like format)
  myfile << endl;
  myfile << "L.s2.L.pedwid = ";
  for(i=0;i<NS2Pad;i++){            
    myfile << Lped_width[i] << " ";
  }


  //RHRS
  myfile << endl << endl;
  myfile << "Values for Right paddles" << endl << endl;

  myfile << "L.s2.R.ped = " ;
  //write pedestal values 
  for(i=0;i<NS2Pad;i++){            
    myfile << Rped[i] << " ";
  }
  // gain values
  myfile << endl;
  myfile << "L.s2.R.gain = " ;
  for(i=0;i<NS2Pad;i++){            
    myfile << rfgain[i] << " ";
  }
  // ped wid values (not in DB, written in DP-like format)
  myfile << endl;
  myfile << "L.s2.R.pedwid = " ;
  for(i=0;i<NS2Pad;i++){            
    myfile << Rped_width[i] << " ";
  }
  


  c1->Print(Form("plots/Lscin_%i.pdf[",runno));
  c1->Print(Form("plots/Lscin_%i.pdf",runno));
  c2->Print(Form("plots/Lscin_%i.pdf",runno));
  c3->Print(Form("plots/Lscin_%i.pdf",runno));
  c4->Print(Form("plots/Lscin_%i.pdf",runno));
  c4->Print(Form("plots/Lscin_%i.pdf]",runno));
  myfile.close();

}

