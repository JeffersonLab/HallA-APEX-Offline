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
#include "Load_more_rootfiles.C"

//bool fexists(const char *filename); 

using namespace std;
void APEX_Rcer_cali()
{
  Int_t right_or_left, run_number;    
  cout << "enter run_number: ";
  cin >> run_number;
  //  TChain *T = new TChain("T"); //first parameter for Tchain is the name of the tree
//  T->Add(Form("/adaqfs/home/a-onl/tritium_work/MARANAL0/replay/Rootfiles/tritium_%d*.root",run_number));
//  T->Add(Form("/chafs1/work1/tritium/Rootfiles/tritium_%d*.root",run_number));
//  T->Add(Form("/adaqfs/home/a-onl/apex/HallA-APEX-Online/replay/rootfiles/tritium_online_%d*.root",run_number));
  //Associate T with a specific rootfile 

  TChain* T = Load_more_rootfiles(run_number);

  cout << T->GetEntries() << endl; 
  
  Double_t min, max, temp;
  Double_t ped_val,ped_wid;
  Int_t max_bin;
  Int_t i, ii;
  Double_t ped[10],ped_width[10];
  Double_t gain[10];
  ofstream myfile;
  myfile.open("R_cer_ped.txt");//,fstream::app);
  myfile << setiosflags(ios::left) << setw(2) << "#" << "   "; 
  myfile << setiosflags(ios::left) << setw(5) << "Ped" << "   ";
  myfile << setiosflags(ios::left) << setw(8) << "gain" << "   "<<endl;


  TCut tritype = "";//"((DR.evtypebits>>1)&1)";
 
  TCanvas *c1 = new TCanvas("c1","c1",1200,1200);
  c1->Divide(5,2); 
  for(i=0;i<10;i++){
    c1->cd(i+1);
    TH1F *tt1 = new TH1F("tt1",Form("FADC_Rcer_%d",i+1),100,5000,6000);
    T->Draw(Form("R.cer.a[%d]>>tt1",i),tritype);
    tt1->SetXTitle(Form("R.cer.a[%d]",i));
    max_bin = tt1->GetMaximumBin();
    min = tt1->GetBinCenter(max_bin) - 40;
    max = tt1->GetBinCenter(max_bin) + 40;
    cout << "min = " << min << " and max = " << max << endl;
    tt1->Fit("gaus","Q","",min,max);
    TF1* gaus1 = tt1->GetFunction("gaus");    
    
    min = gaus1->GetParameter(1) - 3*gaus1->GetParameter(2);
    max = gaus1->GetParameter(1) + 3*gaus1->GetParameter(2);

    tt1->Fit("gaus","Q","",min,max);
    TF1* gaus2 = tt1->GetFunction("gaus");

//    delete tt1;
    ped_val = gaus2->GetParameter(1);
    ped_wid = gaus2->GetParameter(2);    

    gPad->SetLogy();

    ped[i] = ped_val;
    ped_width[i] = ped_wid; 

 //w   cout << "Lcer" << i << endl;
    myfile << setiosflags(ios::left) << setw(2) << setiosflags(ios::fixed) << setprecision(1) << i << "   ";
    myfile << setiosflags(ios::left) << setw(5) << setiosflags(ios::fixed) << setprecision(2) << ped_val << "   ";
    myfile << setiosflags(ios::left) << setw(8) << setiosflags(ios::fixed) << setprecision(1) << ped_wid <<endl;
  } 

  Double_t lfgain[10];

  myfile << setiosflags(ios::left) << setw(2) << "#" << "   ";
  myfile << setiosflags(ios::left) << setw(5) << "gain" << endl;

  Double_t fscale=100;
  TCanvas *c2 = new TCanvas("c2","c2",1200,1200);
  c2->Divide(5,2);

  Double_t ped_sub = 0;
  for(i=0;i<10;i++){
    c2->cd(i+1);
    TH1F *tt2 = new TH1F("tt2",Form("Rcera_p_%d",i+1),300,-100,1200);
    // temporarily usin gped as calculated from previous part
    //    T->Draw(Form("R.cer.a_p[%d]>>tt2",i),tritype);
    ped_sub = ped[i];
    T->Draw(Form("(R.cer.a[%d] - %f)>>tt2",i,ped[i]),tritype);
    tt2->SetXTitle(Form("R.cer.a_p[%d]",i));
    tt2->GetXaxis()->SetRangeUser(0,110);
    Int_t min_bin_no = tt2->GetMinimumBin();
    Int_t min_bin = tt2->GetXaxis()->GetBinCenter(min_bin_no);
    tt2->GetXaxis()->SetRangeUser(-100,1200);

    TF1 *f1 = new TF1("f1", "gaus",min_bin, min_bin+4*fscale);
    tt2->Fit("f1", "Rq");

    min = f1->GetParameter(1) - 1*f1->GetParameter(2);
    max = f1->GetParameter(1) + 1*f1->GetParameter(2);

    tt2->Fit("gaus","Q","",min,max);
    TF1* gaus1 = tt2->GetFunction("gaus");

    lfgain[i]=300*1.0/gaus1->GetParameter(1);


    gPad->SetLogy();

    myfile << setiosflags(ios::left) << setw(2) << setiosflags(ios::fixed) << setprecision(1) << i << "   ";
    myfile << setiosflags(ios::left) << setw(8) << setiosflags(ios::fixed) << setprecision(5) << lfgain[i] <<endl;
  }

  TCut peakcut="";
 
  myfile << setiosflags(ios::left) << setw(2) << "#" << "   ";
  myfile << setiosflags(ios::left) << setw(5) << "Ped" << "   ";
  myfile << setiosflags(ios::left) << setw(8) << "gain" << "   "<<endl;

 
  Double_t scale =100;
//   TCanvas *c3 = new TCanvas("c3","c3",1200,1200);
//   c3->Divide(5,2);
//   for(i=0;i<10;i++){
//     c3->cd(i+1);
//     TH1F *tt3 = new TH1F("tt3",Form("Rcer_%d",i+1),600,200,2000);
//     T->Draw(Form("FbusR.cer.a[%d]>>tt3",i),tritype);
//     tt3->SetXTitle(Form("FbusR.cer.a[%d]",i));
//     max_bin = tt3->GetMaximumBin();
//     min = tt3->GetBinCenter(max_bin) - 20;
//     max = tt3->GetBinCenter(max_bin) + 20;
//     tt3->Fit("gaus","NQ","",min,max);
//     min = gaus->GetParameter(1) - 3*gaus->GetParameter(2);
//     max = gaus->GetParameter(1) + 3*gaus->GetParameter(2);

//     tt3->Fit("gaus","Q","",min,max);

// //    delete tt1;
//     ped_val = gaus->GetParameter(1);
//     ped_wid = gaus->GetParameter(2);

//     gPad->SetLogy();

//     ped[i] = ped_val;
//     ped_width[i] = ped_wid;

//     tt3->GetXaxis()->SetRangeUser(ped[i],ped[i]+scale*1);
//     min_bin = tt3->GetMinimumBin();
//     min_bin = tt3->GetXaxis()->GetBinCenter(min_bin);
//     tt3->GetXaxis()->SetRangeUser(ped[i]-1*scale,ped[i]+5*scale);

//     TF1 *f1;
//     if(i==0)
//        f1 = new TF1("f1", "gaus",min_bin+10, min_bin+1*scale);
//     else f1 = new TF1("f1", "gaus",min_bin, min_bin+1*scale);
//     tt3->Fit("f1", "Rq");

//     min = f1->GetParameter(1) - 1*f1->GetParameter(2);
//     max = f1->GetParameter(1) + 1*f1->GetParameter(2);

//     tt3->Fit("gaus","Q","",min,max);
//     gain[i]=100*1.0/(gaus->GetParameter(1)-ped[i]);

//     TLine *l1 = new TLine(ped[i],0,ped[i],tt3->GetMaximum(););
//     l1->SetLineColor(2);
//     l1->Draw();

//  //w   cout << "Lcer" << i << endl;
//     myfile << setiosflags(ios::left) << setw(2) << setiosflags(ios::fixed) << setprecision(1) << i << "   ";
//     myfile << setiosflags(ios::left) << setw(5) << setiosflags(ios::fixed) << setprecision(2) << ped_val << "   ";
//     myfile << setiosflags(ios::left) << setw(8) << setiosflags(ios::fixed) << setprecision(5) << gain[i]<<endl;

//   }

  c1->Print("Rcer.pdf(");
  //  c1->Print("Rcer.pdf");
  c2->Print("Rcer.pdf)");
  // c3->Print("Rcer.pdf");
  // c3->Print("Rcer.pdf]");
  myfile.close();


}

