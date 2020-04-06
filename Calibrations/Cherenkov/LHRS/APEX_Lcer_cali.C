#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <stdexcept>
#include <cassert>

#include "TTree.h"
//#include "TFile.h"
#include "TString.h"
#include "Load_more_rootfiles.C"


using namespace std;
void APEX_Lcer_cali()
{
  Int_t right_or_left, run_number;    
  cout << "enter run_number: ";
  cin >> run_number;

  TChain* T = Load_more_rootfiles(run_number);


  cout << T->GetEntries() << endl;
  
  Double_t min, max, temp;
  Double_t ped_val,ped_wid;
  Int_t max_bin;
  Int_t i, ii;
  Double_t ped[10],ped_width[10];
  Double_t gain[10];
  ofstream myfile;
  myfile.open("L_cer_ped.txt");//,fstream::app);
  myfile << setiosflags(ios::left) << setw(2) << "#" << "   "; 
  myfile << setiosflags(ios::left) << setw(5) << "Ped" << "   ";
  myfile << setiosflags(ios::left) << setw(8) << "gain" << "   "<<endl;


  TCut tritype = "";//"((DR.evtypebits>>1)&1)";
 
  TCanvas *c1 = new TCanvas("c1","c1",1200,1200);
  c1->Divide(5,2); 
  for(i=0;i<10;i++){
    c1->cd(i+1);
    TH1F *tt1 = new TH1F("tt1",Form("Lcer_%d",i+1),1000,3000,8000);
    T->Draw(Form("L.cer.a[%d]>>tt1",i),tritype);
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

    //gPad->SetLogy();

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

  Double_t fscale=50;
  TCanvas *c2 = new TCanvas("c2","c2",1200,1200);
  c2->Divide(5,2);
  for(i=0;i<10;i++){
    c2->cd(i+1);
    TH1F *tt2 = new TH1F("tt2",Form("Lcera_p_%d",i+1),300,-100,1200);
    T->Draw(Form("L.cer.a_p[%d]>>tt2",i),tritype);

    tt2->GetXaxis()->SetRangeUser(0,200); // 110
    Int_t min_bin_no = tt2->GetMinimumBin();
    Double_t min_bin = tt2->GetXaxis()->GetBinCenter(min_bin_no);
    tt2->GetXaxis()->SetRangeUser(-100,1200);

    cout << "Min bin = " << min_bin << " with bin_no " << min_bin_no << endl;
    
    TF1 *f1 = new TF1("f1", "gaus",min_bin, min_bin+4*fscale);
    tt2->Fit("f1", "Rq");

    min = f1->GetParameter(1) - 1*f1->GetParameter(2);
    max = f1->GetParameter(1) + 1*f1->GetParameter(2);

    cout << "min = " << min << ", max = " << max << endl;
    cout << "Mean = " << f1->GetParameter(1) << ", width = " << f1->GetParameter(2) << endl;
    
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
//     TH1F *tt3 = new TH1F("tt3",Form("Lcer_%d",i+1),600,200,2000);
//     T->Draw(Form("FbusL.cer.a[%d]>>tt3",i),tritype);
//     max_bin = tt3->GetMaximumBin();
//     min = tt3->GetBinCenter(max_bin) - 20;
//     max = tt3->GetBinCenter(max_bin) + 20;
//     tt3->Fit("gaus","Q","",min,max);
//     TF1* gaus1 = tt3->GetFunction("gaus");

//     min = gaus1->GetParameter(1) - 3*gaus1->GetParameter(2);
//     max = gaus1->GetParameter(1) + 3*gaus1->GetParameter(2);

//     tt3->Fit("gaus","RQ","",min,max);
//     TF1* gaus2 = tt3->GetFunction("gaus");

// //    delete tt1;
//     ped_val = gaus2->GetParameter(1);
//     ped_wid = gaus2->GetParameter(2);

//     gPad->SetLogy();

//     ped[i] = ped_val;
//     ped_width[i] = ped_wid;

//     tt3->GetXaxis()->SetRangeUser(ped[i],ped[i]+scale*1);
//     Int_t min_bin_no = tt3->GetMinimumBin();
//     Double_t min_bin = tt3->GetXaxis()->GetBinCenter(min_bin_no);
//     tt3->GetXaxis()->SetRangeUser(ped[i]-1*scale,ped[i]+5*scale);

//     TF1 *f1;
//     if(i==0)
//        f1 = new TF1("f1", "gaus",min_bin+10, min_bin+1*scale);
//     else f1 = new TF1("f1", "gaus",min_bin, min_bin+1*scale);
//     tt3->Fit("f1", "Rq");

//     min = f1->GetParameter(1) - 1*f1->GetParameter(2);
//     max = f1->GetParameter(1) + 1*f1->GetParameter(2);

//     tt3->Fit("gaus","Q","",min,max);
//     TF1* gaus3 = tt3->GetFunction("gaus");

//     gain[i]=100*1.0/(gaus3->GetParameter(1)-ped[i]);

//     TLine *l1 = new TLine(ped[i],0,ped[i],tt3->GetMaximum());
//     l1->SetLineColor(2);
//     l1->Draw();

//  //w   cout << "Lcer" << i << endl;
//     myfile << setiosflags(ios::left) << setw(2) << setiosflags(ios::fixed) << setprecision(1) << i << "   ";
//     myfile << setiosflags(ios::left) << setw(5) << setiosflags(ios::fixed) << setprecision(2) << ped_val << "   ";
//     myfile << setiosflags(ios::left) << setw(8) << setiosflags(ios::fixed) << setprecision(5) << gain[i]<<endl;

//   }

  c1->Print("Lcer.pdf(");
  //  c1->Print("Lcer.pdf");
  c2->Print("Lcer.pdf)");
  // c3->Print("Lcer.pdf");
  // c3->Print("Lcer.pdf]");
  myfile.close();

}

