#include "TString.h"
#include "Load_more_rootfiles.C"


void raster_calib_check(){

  Int_t run_number;
  cout << "Enter run number: ";
  cin >> run_number;

  TString arm;
  TString arm_name;
  cout << "Which arm (L or R)?    ";
  cin >> arm;
  cout << endl << endl;
  
  if (arm == "L" || arm == "l"){
    arm = "L";
    arm_name = "LHRS";
  }
  else if(arm == "R" || arm == "r"){
    arm = "R";
    arm_name = "RHRS";
  }
  else{
    cout << "Choose L or R, you muppet" << endl;
    return 1;
  }




  // TString filename;
  //  filename = new TString(Form("../../apex_tools/rootfiles/apex_%d.0.root",run_number));
   //filename = new TString("../replay/output.root");
  //  TChain *T = new TChain("T");
  //  T->Add(filename);


  TChain* T = Load_more_rootfiles(run_number);


  Double_t rasavg,unrasavg;
  TCut cut;
  cut = "";
  TText *tex;

//   TCanvas *c1 = new TCanvas("c1","urb + rb check",1000,1000);
//   c1->Divide(2,2);
  
//   TH1F *h1 = new TH1F("h1","X of beam on target from raster current",1000,-10,10);
//   TH1F *h2 = new TH1F("h2","X of beam on target from BPMs",1000,-10,10);
//   h1->SetLineColor(2);
//   h2->SetLineColor(4);
//   h1->SetLineWidth(2.5);
//   h2->SetLineWidth(2.5);
//   h1->SetMaximum(8000);
//   h2->SetMaximum(8000);

//   c1->cd(1);
//   T->Draw("1000*rb.Raster.bpma.x>>h1");
//   h1->GetXaxis()->SetTitle("X of beam on target (mm)");
//   h1->GetXaxis()->SetTitleSize(0.04);
//   rasavg = h1->GetMean();
//   tex = new TText(0.40,0.8,Form("Mean Value: %4.4f",rasavg));
//   tex->SetNDC();
//   tex->SetTextSize(0.035);
//   tex->SetTextColor(2);
//   tex->Draw();

//   c1->cd(2);
//   T->Draw("1000*urb.BPMA.x>>h2","","same");
//   h2->GetXaxis()->SetTitle("X of beam on target (mm)");
//   h2->GetXaxis()->SetTitleSize(0.04);
//   unrasavg = h2->GetMean();
//   tex = new TText(0.40,0.8,Form("Mean value: %4.4f",unrasavg));
//   tex->SetNDC();
//   tex->SetTextSize(0.035);
//   tex->SetTextColor(4);
//   tex->Draw();

//   gStyle->SetPalette(1);
//   c1->cd(3);
//   TH2F *h3 = new TH2F("h3","X of beam on target from raster current vs. BPMs",1000,-10,10,1000,-10,10);
//   T->Draw("1000*rb.Raster.bpma.x:1000*urb.BPMA.x>>h3","","collz");
//   h3->GetXaxis()->SetTitle("X of beam on target from raster current (mm)");
//   h3->GetXaxis()->SetTitleSize(0.04);
//   h3->GetYaxis()->SetTitle("X of beam on target from BPMs (mm)");
//   h3->GetYaxis()->SetTitleSize(0.04);

//   c1->cd(4);
//   tex = new TText(0.2,0.7,Form("Run# %d",run_number));
//   tex->SetNDC();
//   tex->SetTextSize(0.05);
//   tex->Draw();
//   tex = new TText(0.2,0.65,"Raster ON");
//   tex->SetNDC();
//   tex->SetTextSize(0.05);
//   tex->Draw();

  TCanvas *c2 = new TCanvas("c2",arm_name + ":  "  + arm + "urb + rb check "+ Form("for run %d",run_number),1000,1000);
  c2->Divide(2,3);

  TH1F *rb[6];
  TH1F *urb[6];

  gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);

  Double_t rbavg[6],urbavg[6],rbrms[6],urbrms[6];

  c2->cd(1);
  rb[0] = new TH1F("rb[0]","X position of beam at BPMA",1000,-10,10);
  urb[0] = new TH1F("urb[0]","",1000,-10,10);
  rb[0]->SetLineColor(2);
  urb[0]->SetLineColor(4);
  rb[0]->SetMaximum(30000);
  T->Draw("1000*" + arm + "rb.Raster2.bpma.x>>rb[0]");
  T->Draw("1000*" + arm + "urb.BPMA.x>>urb[0]","","same");
  rbavg[0] = rb[0]->GetMean();
  urbavg[0] = urb[0]->GetMean();
  rbrms[0] = rb[0]->GetRMS();
  urbrms[0] = urb[0]->GetRMS();
  rb[0]->GetXaxis()->SetTitle("x position (mm)");
  rb[0]->GetXaxis()->CenterTitle();
  tex = new TText(0.17,0.8,"As determined by raster 2 current");
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(2);
  tex->Draw();
  tex = new TText(0.17,0.75,Form("Mean value: %4.4f",rbavg[0]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(2);
  tex->Draw();
  tex = new TText(0.17,0.7,Form("RMS: %4.4f",rbrms[0]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(2);
  tex->Draw();
  tex = new TText(0.65,0.8,"As determined by BPMs");
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(4);
  tex->Draw();
  tex = new TText(0.65,0.75,Form("Mean value: %4.4f",urbavg[0]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(4);
  tex->Draw();
  tex = new TText(0.65,0.7,Form("RMS: %4.4f",urbrms[0]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(4);
  tex->Draw();

  c2->cd(2);
  rb[1] = new TH1F("rb[1]","Y position of beam at BPMA",1000,-10,10);
  urb[1] = new TH1F("urb[1]","",1000,-10,10);
  rb[1]->SetLineColor(2);
  urb[1]->SetLineColor(4);
  rb[1]->SetMaximum(30000);
  T->Draw("1000*" + arm + "rb.Raster2.bpma.y>>rb[1]");
  T->Draw("1000*" + arm + "urb.BPMA.y>>urb[1]","","same");
  rbavg[1] = rb[1]->GetMean();
  urbavg[1] = urb[1]->GetMean();
  rbrms[1] = rb[1]->GetRMS();
  urbrms[1] = urb[1]->GetRMS();
  rb[1]->GetXaxis()->SetTitle("y position (mm)");
  rb[1]->GetXaxis()->CenterTitle();
  tex = new TText(0.17,0.8,"As determined by raster 2 current");
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(2);
  tex->Draw();
  tex = new TText(0.17,0.75,Form("Mean value: %4.4f",rbavg[1]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(2);
  tex->Draw();
  tex = new TText(0.17,0.7,Form("RMS: %4.4f",rbrms[1]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(2);
  tex->Draw();
  tex = new TText(0.65,0.8,"As determined by BPMs");
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(4);
  tex->Draw();
  tex = new TText(0.65,0.75,Form("Mean value: %4.4f",urbavg[1]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(4);
  tex->Draw();
  tex = new TText(0.65,0.7,Form("RMS: %4.4f",urbrms[1]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(4);
  tex->Draw();
  
  c2->cd(3);
  rb[2] = new TH1F("rb[2]","X position of beam at BPMB",1000,-10,10);
  urb[2] = new TH1F("urb[2]","",1000,-10,10);
  rb[2]->SetLineColor(2);
  urb[2]->SetLineColor(4);
  rb[2]->SetMaximum(30000);
  T->Draw("1000*" + arm + "rb.Raster2.bpmb.x>>rb[2]");
  T->Draw("1000*" + arm + "urb.BPMB.x>>urb[2]","","same");
  rbavg[2] = rb[2]->GetMean();
  urbavg[2] = urb[2]->GetMean();
  rbrms[2] = rb[2]->GetRMS();
  urbrms[2] = urb[2]->GetRMS();
  rb[2]->GetXaxis()->SetTitle("x position (mm)");
  rb[2]->GetXaxis()->CenterTitle();
  tex = new TText(0.17,0.8,"As determined by raster 2 current");
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(2);
  tex->Draw();
  tex = new TText(0.17,0.75,Form("Mean value: %4.4f",rbavg[2]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(2);
  tex->Draw();
  tex = new TText(0.17,0.7,Form("RMS: %4.4f",rbrms[2]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(2);
  tex->Draw();
  tex = new TText(0.65,0.8,"As determined by BPMs");
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(4);
  tex->Draw();
  tex = new TText(0.65,0.75,Form("Mean value: %4.4f",urbavg[2]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(4);
  tex->Draw();
  tex = new TText(0.65,0.7,Form("RMS: %4.4f",urbrms[2]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(4);
  tex->Draw();

  c2->cd(4);
  rb[3] = new TH1F("rb[3]","Y position of beam at BPMB",1000,-10,10);
  urb[3] = new TH1F("urb[3]","",1000,-10,10);
  rb[3]->SetLineColor(2);
  urb[3]->SetLineColor(4);
  rb[3]->SetMaximum(30000);
  T->Draw("1000*" + arm + "rb.Raster2.bpmb.y>>rb[3]");
  T->Draw("1000*" + arm + "urb.BPMB.y>>urb[3]","","same");
  rbavg[3] = rb[3]->GetMean();
  urbavg[3] = urb[3]->GetMean();
  rbrms[3] = rb[3]->GetRMS();
  urbrms[3] = urb[3]->GetRMS();
  rb[3]->GetXaxis()->SetTitle("y position (mm)");
  rb[3]->GetXaxis()->CenterTitle();
  tex = new TText(0.17,0.8,"As determined by raster 2 current");
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(2);
  tex->Draw();
  tex = new TText(0.17,0.75,Form("Mean value: %4.4f",rbavg[3]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(2);
  tex->Draw();
  tex = new TText(0.17,0.7,Form("RMS: %4.4f",rbrms[3]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(2);
  tex->Draw();
  tex = new TText(0.65,0.8,"As determined by BPMs");
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(4);
  tex->Draw();
  tex = new TText(0.65,0.75,Form("Mean value: %4.4f",urbavg[3]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(4);
  tex->Draw();
  tex = new TText(0.65,0.7,Form("RMS: %4.4f",urbrms[3]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(4);
  tex->Draw();
  
  c2->cd(5);
  rb[4] = new TH1F("rb[4]","X position of beam at target",1000,-10,10);
  urb[4] = new TH1F("urb[4]","",1000,-10,10);
  rb[4]->SetLineColor(2);
  urb[4]->SetLineColor(4);
  rb[4]->SetMaximum(30000);
  T->Draw("1000*" + arm + "rb.x>>rb[4]");
  T->Draw("1000*" + arm + "urb.x>>urb[4]","","same");
  rbavg[4] = rb[4]->GetMean();
  urbavg[4] = urb[4]->GetMean();
  rbrms[4] = rb[4]->GetRMS();
  urbrms[4] = urb[4]->GetRMS();
  rb[4]->GetXaxis()->SetTitle("x position (mm)");
  rb[4]->GetXaxis()->CenterTitle();
  tex = new TText(0.17,0.8,"As determined by raster 2 current");
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(2);
  tex->Draw();
  tex = new TText(0.17,0.75,Form("Mean value: %4.4f",rbavg[4]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(2);
  tex->Draw();
  tex = new TText(0.17,0.7,Form("RMS: %4.4f",rbrms[4]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(2);
  tex->Draw();
  tex = new TText(0.65,0.8,"As determined by BPMs");
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(4);
  tex->Draw();
  tex = new TText(0.65,0.75,Form("Mean value: %4.4f",urbavg[4]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(4);
  tex->Draw();
  tex = new TText(0.65,0.7,Form("RMS: %4.4f",urbrms[4]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(4);
  tex->Draw();

  c2->cd(6);
  rb[5] = new TH1F("rb[5]","Y position of beam at target",1000,-10,10);
  urb[5] = new TH1F("urb[5]","",1000,-10,10);
  rb[5]->SetLineColor(2);
  urb[5]->SetLineColor(4);
  rb[5]->SetMaximum(30000);
  T->Draw("1000*" + arm + "rb.y>>rb[5]");
  T->Draw("1000*" + arm + "urb.y>>urb[5]","","same");
  rbavg[5] = rb[5]->GetMean();
  urbavg[5] = urb[5]->GetMean();
  rbrms[5] = rb[5]->GetRMS();
  urbrms[5] = urb[5]->GetRMS();
  rb[5]->GetXaxis()->SetTitle("y position (mm)");
  rb[5]->GetXaxis()->CenterTitle();
  tex = new TText(0.17,0.8,"As determined by raster 2 current");
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(2);
  tex->Draw();
  tex = new TText(0.17,0.75,Form("Mean value: %4.4f",rbavg[5]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(2);
  tex->Draw();
  tex = new TText(0.17,0.7,Form("RMS: %4.4f",rbrms[5]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(2);
  tex->Draw();
  tex = new TText(0.65,0.8,"As determined by BPMs");
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(4);
  tex->Draw();
  tex = new TText(0.65,0.75,Form("Mean value: %4.4f",urbavg[5]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(4);
  tex->Draw();
  tex = new TText(0.65,0.7,Form("RMS: %4.4f",urbrms[5]));
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->SetTextColor(4);
  tex->Draw();
}
