/**************************************************************
  APEX_VDCtable.C
  John Williamson    
  22nd February, 2021

  This script creates a lookup table for all 4 VDC planes for a chosen arm and run. This is done by integrating over the timnig TDC spectrum (converted to real time with t0 subtracted). The table is normalised such that the maximum distance is the size of the VDC drift cell. 
  
*************************************************************/
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

#include "Load_more_rootfiles.C"
#include "file_def.h"
#include "TTD_namespace.h"




#define NPLANE 4
#define MAX_HIT 1000
#define MAX_ENTRIES 1000000

#define DEG_TO_RAD 0.017453278

const int kBUFLEN = 150;


const Double_t kLongestDist = 35.12e-3; // (m)
const Double_t kTimeRes = 0.5e-9;    // (s)

const int kTrue = 1;
const int kFalse = 0;

std::vector<Double_t> wtime[NPLANE];
std::vector<Double_t> tanTh[NPLANE];
std::vector<Double_t> trdist[NPLANE];

 



/*
  Creates a dynamically sized TTD lookup table based on
  a low and high time supplied

  planename should be a string like "L.u1.vdc"
  nentries is the number of entries in the supplied treefile to
    use in making the table. a value of 0 means to use all of the
    availible entries
  low and high hive the smallest and largest values of time to be considered in bulding the table (default values give full scope)
 */
void APEX_VDCtable(const char *arm, Int_t runnumber = -1,
		   Double_t low = -10.0e-9, Double_t high = 450.0e-9,  Int_t nentries=0)
{


  TChain* T = new TChain("T");

  if(!strcmp(arm,"L")){
  T = Load_more_rootfiles(runnumber);
  }
  else if (!strcmp(arm,"R")){  
    T = Load_more_rootfiles(runnumber);
  }
  else{
    cout << "arm must be L or R, " << arm << " not acceptable" << endl;
    return;
  }



  // PID cuts (different for Left and Right arms)
  Double_t Cer_cut = 0.0;
  Double_t Ps_cut = 0.0;
  Double_t Ps_Sh_cut_l = 0.0;
  Double_t Ps_Sh_cut_h = 0.0;
    
  Int_t trg = 0;
  // here is where trigger is defined
  if(!strcmp(arm,"L")){
    trg = 1;
    Cer_cut = 1500.0;
    Ps_cut = 0.3;
    Ps_Sh_cut_l = 0.625;
    Ps_Sh_cut_h = 1.1;
    cout << "left trigger" << endl;
  }
  else if (!strcmp(arm,"R")){
    cout << "right trigger" << endl;
    trg = 5;
    Cer_cut = 650.0;
    Ps_cut = 0.2;
    Ps_Sh_cut_l = 0.51;
    Ps_Sh_cut_h = 1.11;
    // trg = 6;
  }
  else{
    cout << "arm must be L or R, " << arm << " not acceptable" << endl;
    return;

  }

  
  
  const char plane[NPLANE][8] = {"u1", "u2", "v1", "v2"};
  Double_t ang[NPLANE] = {-45.0, -45.0, 45.0, 45.0};

  Double_t nhit[NPLANE], ntr;
  Double_t hittime[NPLANE][MAX_HIT], hittrknum[NPLANE][MAX_HIT],
    hittrdist[NPLANE][MAX_HIT];
  Double_t d_th[MAX_HIT], d_ph[MAX_HIT];
  Double_t cer_sum, ps_e, sh_e;
  Double_t tr_p[100];
  THaEvent* evt = 0;

  Int_t    nent[NPLANE];
  
  Int_t i, j, hit;

  Double_t evttype;


    // Set up branches

  T->SetBranchStatus("Event_Branch*", kTRUE);
  T->SetBranchAddress("Event_Branch", &evt);

  T->SetBranchStatus("DR.evtypebits", kTRUE);
  T->SetBranchAddress("DR.evtypebits", &evttype);

  T->SetBranchStatus(Form("%s.tr.n", arm), kTRUE);
  T->SetBranchAddress(Form("%s.tr.n", arm), &ntr);

  T->SetBranchStatus(Form("%s.tr.d_th", arm), kTRUE);
  T->SetBranchAddress(Form("%s.tr.d_th", arm), &d_th);
  T->SetBranchStatus(Form("%s.tr.d_ph", arm), kTRUE);
  T->SetBranchAddress(Form("%s.tr.d_ph", arm), &d_ph);

  T->SetBranchStatus(Form("%s.tr.p",arm),kTRUE);
  T->SetBranchAddress(Form("%s.tr.p",arm),tr_p);
  
  T->SetBranchStatus(Form("%s.cer.asum_c",arm),kTRUE);
  T->SetBranchAddress(Form("%s.cer.asum_c",arm),&cer_sum);
  
  if(!strcmp(arm,"L")){
    T->SetBranchStatus("L.prl1.e",kTRUE);
    T->SetBranchAddress("L.prl1.e",&ps_e);
    T->SetBranchStatus("L.prl2.e",kTRUE);
    T->SetBranchAddress("L.prl2.e",&sh_e);
  }
  else if(!strcmp(arm,"R")){
    T->SetBranchStatus("R.ps.e",kTRUE);
    T->SetBranchAddress("R.ps.e",&ps_e);
    T->SetBranchStatus("R.ps.e",kTRUE);
    T->SetBranchAddress("R.sh.e",&sh_e);
  }
  

  for( i = 0; i < NPLANE; i++ ){
    T->SetBranchStatus(Form("%s.vdc.%s.nhit", arm, plane[i]), kTRUE);
    T->SetBranchAddress(Form("%s.vdc.%s.nhit", arm, plane[i]), &nhit[i]);

    T->SetBranchStatus(Form("%s.vdc.%s.time", arm, plane[i]), kTRUE);
    T->SetBranchAddress(Form("%s.vdc.%s.time", arm, plane[i]), hittime[i]);

    T->SetBranchStatus(Form("%s.vdc.%s.trknum", arm, plane[i]), kTRUE);
    T->SetBranchAddress(Form("%s.vdc.%s.trknum", arm, plane[i]), hittrknum[i]);

//     T->SetBranchStatus(Form("%s.vdc.%s.ltrdist", arm, plane[i]), kTRUE);
//     T->SetBranchAddress(Form("%s.vdc.%s.ltrdist", arm, plane[i]), hittrdist[i]);
    T->SetBranchStatus(Form("%s.vdc.%s.trdist", arm, plane[i]), kTRUE);
    T->SetBranchAddress(Form("%s.vdc.%s.trdist", arm, plane[i]), hittrdist[i]);

    nent[i] = 0;
  }

  // Double_t high = 450e-9;
  // Double_t low = -10.0e-9;

  Int_t num_items=0;
  Int_t num_bins = (high - low)/kTimeRes;

  cout<<"num bins = "<<num_bins<<endl;

  
  
  //  vector<Float_t> table(num_bins, 0);
  // histogram of hit times
  TH1F *hist[NPLANE];


  // histogram of slope vs hit times
  TH2F *hTimeVSlope[NPLANE];

  
  for( i = 0; i < NPLANE; i++ ){
    hist[i] = new TH1F(Form("TDC_%s_%s", arm, plane[i]),Form("TDC spectrum %s %s", arm, plane[i]), num_bins, low, high);
    hTimeVSlope[i] = new TH2F(Form("hTimeVSlope_%s_%s", arm, plane[i]),Form("Time vs Slope %s %s", arm, plane[i]),100,1,2, num_bins, low, high);
  }

  

  // Fill Timing spectrum

  if(nentries == 0){
    nentries = T->GetEntries();
  }


  Double_t this_slope;

  for( i = 0; i < nentries; i++ ){
    T->GetEntry(i);
    if( (i%5000)==0 ) { cout << "Entry " << i << endl; }

   
    //    if( ntr == 1 && TTD_func::passtrg(Int_t(evttype), trg) && cer_sum > Cer_cut && ps_e/(1e3*tr_p[0]) > Ps_cut && (ps_e+sh_e)/(1e3*tr_p[0]) > Ps_Sh_cut_l &&  (ps_e+sh_e)/(1e3*tr_p[0]) < Ps_Sh_cut_h ){
    if( ntr == 1 && TTD_func::passtrg(Int_t(evttype), trg)  && cer_sum > Cer_cut){
      for( j = 0; j < NPLANE; j++ ){
	this_slope = d_th[0]*cos(ang[j]*DEG_TO_RAD) + d_ph[0]*sin(ang[j]*DEG_TO_RAD);
	for( hit = 0; hit < nhit[j] && nent[j] < MAX_ENTRIES; hit++ ){
	  if( low < hittime[j][hit] &&
	      hittime[j][hit]  < high &&
	      hittrdist[j][hit]< kLongestDist &&
	      hittrknum[j][hit] == 1 ){

	    wtime[j].push_back(hittime[j][hit]);
	    //	    tanTh[j].push_back(this_slope);
	    trdist[j].push_back(hittrdist[j][hit]);
	    
	    hist[j]->Fill( hittime[j][hit]);
	    hTimeVSlope[j]->Fill(1/this_slope,hittime[j][hit]);
	    nent[j]++;
	  }
	}
      }
    }
  }



  // cut noise based on 'real' distance vs time
  NoiseCut* Real_cut[NPLANE];
  TH2D *htime_dist_real[NPLANE];

  for( i = 0; i < NPLANE; i++ ){
  
    htime_dist_real[i] = new TH2D(Form("htime_dist_real_%s", plane[i]), Form("%s TTD", plane[i]), 290, -10,280, 200, 0.0, 0.015);
    htime_dist_real[i]->GetXaxis()->SetTitle("Drift Time(ns)");
    htime_dist_real[i]->GetXaxis()->CenterTitle();
    htime_dist_real[i]->GetYaxis()->SetTitle("Track Dist (m)");
    htime_dist_real[i]->GetYaxis()->CenterTitle();

    for( j = 0; j < nent[i]; j++ ){
      htime_dist_real[i]->Fill(1e9*(wtime[i][j]),trdist[i][j]);
    }

    Real_cut[i] =  new NoiseCut(htime_dist_real[i]);
    Real_cut[i]->PassNoiseCut(wtime[i],trdist[i]);
    nent[i] = wtime[i].size();
  }
          
  

  TCanvas *c1[NPLANE];

  for( i = 0; i < NPLANE; i++ ){
    c1[i] = new TCanvas(plane[i], plane[i], 640, 480);
    c1[i]->Divide(2,1);

    c1[i]->cd(1);
    hist[i]->Draw();

    c1[i]->cd(2);
    hTimeVSlope[i]->Draw("colz");
  }


  // find first and last bins in time with entries for each plane    
  Int_t low_bin[4] = {0};
  Int_t high_bin[4] = {0};

  
  Double_t low_val[4] = {0};
  Double_t high_val[4] = {0};


  

  Bool_t first[4] = {false};

  for(Int_t i = 0; i<NPLANE; i++ ){
    first[i] = true;
  }
  


  for(Int_t i = 0; i<NPLANE; i++ ){
    
    Int_t j = 0;


    while(hist[i]->GetBinContent(j) == 0){
      j++;
    }

    
    low_bin[i] = j;
    low_val[i] = (low_bin[i])*kTimeRes + low;

    
    Int_t k = 0; // check against gaving small gap in time spectrum at smaller times
    
    while(hist[i]->GetBinContent(j) > 0 || k < 10){
      high_bin[i] = j;      
      j++;
      k++;
    }

    cout << "passed second while loop" << endl;

    
    cout << "lowbin[" << i << "] = " << low_bin[i] << endl;
    cout << "highbin[" << i << "] = " << high_bin[i] << endl;

  }
    
    
    // for(Int_t j = 0; j<hist[i]->GetNbinsX(); j++){

    //   if (hist[i]->GetBinContent(j) > 0){
	  
    // 	high_bin[i] = j;

    // 	if (first[i] ){
    // 	  low_bin[i] = j;
    // 	  low_val[i] = (low_bin[i])*kTimeRes + low;
    // 	  first[i] = false;
    // 	}
    //   }                  
    // }

    


  
  // set provisional value of K normlisation parameter (determined by size of drift cell)
  Double_t K[4] = {1.0};
  

  
  Int_t NBins[4] = {0}; // per-plane entries in look-up tables   
  for( i = 0; i < NPLANE; i++ ){
    NBins[i] = high_bin[i] - low_bin[i];
    cout << "Plane " << plane[i] << ", high_bin = " << high_bin[i] << " (time = " << high_bin[i]*kTimeRes + low << "), low_bin = " << low_bin[i] << " (time = " << low_bin[i]*kTimeRes + low << "), NBins = " << NBins[i] << endl;
    K[i] = 1.0;
  }

  vector<Float_t> tables[NPLANE];
  

  for( i = 0; i < NPLANE; i++ ){
    tables[i].push_back(hist[i]->GetBinContent(low_bin[i]));
    cout << "table[" << i << "][0] = " << tables[i][0] << endl;
    

    for(Int_t j=1; j<NBins[i]; j++) {
      tables[i].push_back(tables[i][j-1] + K[i]*hist[i]->GetBinContent(j+low_bin[i]));
    }
    
    cout << "Completed table for " << plane[i] << endl;
  }
  
  cout << "Filled in tables for all planes" << endl;

  
  
  // autocalculate the scale factor, if the user'd like
  
  
  // for( i = 0; i < NPLANE; i++ ){
  //   K[i] = kLongestDist/tables[i][NBins[i]-1];
  //   cout<<"K["<<i<<"] estimate: "<<K[i]<<endl;
  // }

  // for( i = 0; i < NPLANE; i++ ){
  //   cout << "Table for " << plane[i] << endl << endl;
  //   for(Int_t j=0; j<NBins[i]; j++){
  //     tables[i][j] *= K[i];
  //     cout << tables[i][j] << " ";
  //     if (j%10 == 0 && j>0){
  // 	cout << endl;
  //     }
  //   }
  //   cout << endl << endl;
  // }


  // ask whether to replace the values in the database
  char input[kBUFLEN];
  cout<<"Do you want to rebuild the database using these values? [y/n] ";
  input[0] = '\0';
  fgets(input, kBUFLEN, stdin);

  if(input[0] != 'y') {
    cout<<"Exiting without rebuilding database."<<endl;
    //goto cleanup;
    exit(1);
  }


  cout<<"Rebuilding database..."<<endl;



  for( i = 0; i < NPLANE; i++ ){
    if(TTD_func::SaveNewTTDData(tables[i], NBins[i], low_val[i], arm, plane[i], runnumber, "unnormalised"))
	cout<<"Done for "<<plane[i]<<endl;
    else{
      cout<<"Failed for "<<plane[i]<<endl;      
    }
    
  }
  

}
