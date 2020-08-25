/* Script designed to check position of vertical foils according to according to analyzer raster and compare to survey

Author: John Williamson
Date: 29/4/20
*/

#include "InputAPEXL.h"
#include "Load_more_rootfiles.C"
#include "file_def.h"


void beam_x_check(){


  // load in root files for Vertical foil runs


  gStyle->SetOptStat(1101);
  
  Int_t V1_r1 = 4648; // RHRS V1 run
  Int_t V1_r2 = 4766; // LHRS V1 run
  Int_t V1_r3 = 4179; // LHRS V1 run
  
  
  
  TChain* T1 = Load_more_rootfiles(V1_r1); //4648 V1
  TChain* T1_2 = Load_more_rootfiles(V1_r2); //4766 V1
  TChain* T1_3 = Load_more_rootfiles(V1_r3); //4179 V1

  
  Int_t V2_r1 = 4647;
  Int_t V2_r2 = 4768;
  Int_t V2_r3 = 4181;
  
  TChain* T2 = Load_more_rootfiles(4647); // 4647 V2
  TChain* T2_2 = Load_more_rootfiles(4768); // 4768 V2
  TChain* T2_3 = Load_more_rootfiles(4181); // 4181 V2

  
  Int_t V3_r1 = 4650;
  Int_t V3_r2 = 4769;
  Int_t V3_r3 = 4180;
  
  TChain* T3 = Load_more_rootfiles(4650); // 4650 V3
  TChain* T3_2 = Load_more_rootfiles(4769); // 4769 V3
  TChain* T3_3 = Load_more_rootfiles(4180); // 4180 V3
  

  // format histograms

  TH1F* hv1 = new TH1F(Form("V1_run_%d",V1_r1),"V1 foil",1000,-3,3);
  hv1->GetXaxis()->SetTitle("X Position [mm]");
  TH1F* hv1_2 = new TH1F(Form("V1_run_%d",V1_r2),"V1 foil",1000,-3,3);
  hv1_2->GetXaxis()->SetTitle("X Position [mm]");
  TH1F* hv1_3 = new TH1F(Form("V1_run_%d",V1_r3),"V1 foil",1000,-3,3);
  hv1_3->GetXaxis()->SetTitle("X Position [mm]");

  
  
  TH1F* hv2 = new TH1F(Form("V2_run_%d",V2_r1),"V2 foil",1000,-3,3);
  hv2->GetXaxis()->SetTitle("X Position [mm]");
  TH1F* hv2_2 = new TH1F(Form("V2_run_%d",V2_r2),"V2 foil",1000,-3,3);
  hv2_2->GetXaxis()->SetTitle("X Position [mm]");
  TH1F* hv2_3 = new TH1F(Form("V2_run_%d",V2_r3),"V2 foil",1000,-3,3);
  hv2_3->GetXaxis()->SetTitle("X Position [mm]");
  //  hv2_2->SetMinimum(1);
  
  TH1F* hv3 = new TH1F(Form("V3_run_%d",V3_r1),"V3 foil",1000,-3,3);
  hv3->GetXaxis()->SetTitle("X Position [mm]");
  TH1F* hv3_2 = new TH1F(Form("V3_run_%d",V3_r2),"V3 Foil",1000,-3,3);
  hv3_2->GetXaxis()->SetTitle("X Position [mm]");
  TH1F* hv3_3 = new TH1F(Form("V3_run_%d",V3_r3),"V3 Foil",1000,-3,3);
  hv3_3->GetXaxis()->SetTitle("X Position [mm]");
  

  
  TCanvas *c1 = new TCanvas("c1","Beam x comparison",1000,1000);
  
  c1->Divide(3,1);

  c1->cd(1);

  Double_t Z_dist = 0.0;

  Z_dist = targetfoils[8];
  
  Double_t BeamZDir_average_  = 5.131;

  hv1->SetMarkerSize(0.5);
  hv1->SetMarkerStyle(20);
  hv1->SetMarkerColor(kBlue);  
  T1->Draw(Form("(Lrb.x + ((%f/%f)*Lrb.dir.x))*1000>>V1_run_%d",targetfoils[8],BeamZDir_average,V1_r1),"","P");
  
  hv1_2->SetMarkerSize(0.5);
  hv1_2->SetMarkerStyle(20);
  hv1_2->SetMarkerColor(kRed);
  hv1_2->SetLineColor(kRed);
  hv1_2->SetLineWidth(1);
  T1_2->Draw(Form("(Lrb.x + ((%f/%f)*Lrb.dir.x))*1000>>V1_run_%d",targetfoils[8],BeamZDir_average,V1_r2),"","sames P");

  hv1_3->SetMarkerSize(0.5);
  hv1_3->SetMarkerStyle(20);
  hv1_3->SetMarkerColor(kGreen + 3);
  hv1_3->SetLineColor(kGreen + 3);
  hv1_3->SetLineWidth(1);
  T1_3->Draw(Form("(Lrb.x + ((%f/%f)*Lrb.dir.x))*1000>>V1_run_%d",targetfoils[8],BeamZDir_average,V1_r3),"","sames P");
  
  TLine *l1 = new TLine(targetfoilsX[8]*1000,0,targetfoilsX[8]*1000,hv1->GetMaximum());
  l1->SetLineWidth(2);
  l1->Draw("same");

  Double_t xdiff_V1_1 = hv1->GetMean() - targetfoilsX[8]*1000;
  cout << "xdiff_V1_1 = " << xdiff_V1_1 << endl;

  Double_t xdiff_V1_2 = hv1_2->GetMean() - targetfoilsX[8]*1000;
  cout << "xdiff_V1_2 = " << xdiff_V1_2 << endl;

  Double_t xdiff_V1_3 = hv1_3->GetMean() - targetfoilsX[8]*1000;
  cout << "xdiff_V1_3 = " << xdiff_V1_3 << endl;

  Double_t av_xdiff_V1 = (xdiff_V1_1 + xdiff_V1_2 + xdiff_V1_3)/3.0;
  
  TLegend* leg_1 = new TLegend(.1,.65,.40,.9,"Key");
  leg_1->SetFillColor(0);
  leg_1->SetTextSize(0.02); 
  leg_1->AddEntry(hv1,Form("Beam x at V1 (run %d)",V1_r1));
  leg_1->AddEntry(hv1_2,Form("Beam x at V1 (run %d)",V1_r2));
  leg_1->AddEntry(hv1_3,Form("Beam x at V1 (run %d)",V1_r3));
  leg_1->AddEntry(l1,"V1 survey","l");
  leg_1->AddEntry((TObject*)0,Form("Mean Diff = %.3f mm",av_xdiff_V1),"");
  leg_1->AddEntry((TObject*)0,"(analyzer - survey)","");

  c1->Update();

  TPaveStats *ps_V1_1 = (TPaveStats*)hv1->GetListOfFunctions()->FindObject("stats");

  ps_V1_1->SetY1NDC(0.50);
  ps_V1_1->SetY2NDC(0.63);
  ps_V1_1->SetX1NDC(0.1);
  ps_V1_1->SetX2NDC(0.3);

  TPaveStats *ps_V1_2 = (TPaveStats*)hv1_2->GetListOfFunctions()->FindObject("stats");

  ps_V1_2->SetY1NDC(0.32);
  ps_V1_2->SetY2NDC(0.45);
  ps_V1_2->SetX1NDC(0.1);
  ps_V1_2->SetX2NDC(0.3);

  TPaveStats *ps_V1_3 = (TPaveStats*)hv1_3->GetListOfFunctions()->FindObject("stats");
    
  ps_V1_3->SetY1NDC(0.14);
  ps_V1_3->SetY2NDC(0.27);
  ps_V1_3->SetX1NDC(0.1);
  ps_V1_3->SetX2NDC(0.3);
      
  c1->Modified();

  
  leg_1->Draw("same");
  
  c1->cd(2);
  
  Z_dist = targetfoils[9];

  hv2->SetMarkerSize(0.5);
  hv2->SetMarkerStyle(20);
  hv2->SetMarkerColor(kBlue);  
  T2->Draw(Form("(Lrb.x + ((%f/%f)*Lrb.dir.x))*1000>>V2_run_%d",targetfoils[9],BeamZDir_average,V2_r1),"","P");
  TLine *l2 = new TLine(targetfoilsX[9]*1000,0,targetfoilsX[9]*1000,hv2->GetMaximum());
  l2->SetLineWidth(2);
  l2->Draw("same");


  hv2_2->SetMarkerSize(0.5);
  hv2_2->SetMarkerStyle(20);
  hv2_2->SetLineColor(kRed);
  hv2_2->SetMarkerColor(kRed);
  T2_2->Draw(Form("(Lrb.x + ((%f/%f)*Lrb.dir.x))*1000>>V2_run_%d",targetfoils[9],BeamZDir_average,V2_r2),"","sames P");


  hv2_3->SetMarkerSize(0.5);
  hv2_3->SetMarkerStyle(20);
  hv2_3->SetLineColor(kGreen + 3);
  hv2_3->SetMarkerColor(kGreen + 3);
  T2_3->Draw(Form("(Lrb.x + ((%f/%f)*Lrb.dir.x))*1000>>V2_run_%d",targetfoils[9],BeamZDir_average,V2_r3),"","sames P");
  

  Double_t xdiff_V2_1 = hv2->GetMean() - targetfoilsX[9]*1000;
  cout << "xdiff_V2_1 = " << xdiff_V2_1 << endl;

  Double_t xdiff_V2_2 = hv2_2->GetMean() - targetfoilsX[9]*1000;
  cout << "xdiff_V2_2 = " << xdiff_V2_2 << endl;

  Double_t xdiff_V2_3 = hv2_3->GetMean() - targetfoilsX[9]*1000;
  cout << "xdiff_V2_3 = " << xdiff_V2_3 << endl;
  
  Double_t av_xdiff_V2 = (xdiff_V2_1 + xdiff_V2_2 + xdiff_V2_3)/3.0;
  
  TLegend* leg_2 = new TLegend(.53,.65,.9,.9,"Key");
  leg_2->SetFillColor(0);
  leg_2->SetTextSize(0.02); 
  leg_2->AddEntry(hv2,Form("Beam x at V2 (run %d)",V2_r1));
  leg_2->AddEntry(hv2_2,Form("Beam x at V2 (run %d)",V2_r2));
  leg_2->AddEntry(hv2_3,Form("Beam x at V2 (run %d)",V2_r3));
  leg_2->AddEntry(l2,"V2 survey","l");
  leg_2->AddEntry((TObject*)0,Form("Mean Diff = %.3f mm",av_xdiff_V2),"");
  leg_2->AddEntry((TObject*)0,"(analyzer - survey)","");
  
  leg_2->Draw("same");

  c1->Update();

  TPaveStats *ps_V2_1 = (TPaveStats*)hv2->GetListOfFunctions()->FindObject("stats");

  ps_V2_1->SetY1NDC(0.50);
  ps_V2_1->SetY2NDC(0.63);
  ps_V2_1->SetX1NDC(0.7);
  ps_V2_1->SetX2NDC(0.9);

  TPaveStats *ps_V2_2 = (TPaveStats*)hv2_2->GetListOfFunctions()->FindObject("stats");

  ps_V2_2->SetY1NDC(0.32);
  ps_V2_2->SetY2NDC(0.45);
  ps_V2_2->SetX1NDC(0.7);
  ps_V2_2->SetX2NDC(0.9);

  TPaveStats *ps_V2_3 = (TPaveStats*)hv2_3->GetListOfFunctions()->FindObject("stats");

  ps_V2_3->SetY1NDC(0.14);
  ps_V2_3->SetY2NDC(0.27);
  ps_V2_3->SetX1NDC(0.7);
  ps_V2_3->SetX2NDC(0.9);

  
  c1->Modified();


  // V3 foil
  
  c1->cd(3);
  
  Z_dist = targetfoils[10];

  hv3->SetMarkerSize(0.5);
  hv3->SetMarkerStyle(20);
  hv3->SetMarkerColor(kBlue);  
  T3->Draw(Form("(Lrb.x + ((%f/%f)*Lrb.dir.x))*1000>>V3_run_%d",targetfoils[10],BeamZDir_average,V3_r1),"","P");
  TLine *l3 = new TLine(targetfoilsX[10]*1000,0,targetfoilsX[10]*1000,hv3->GetMaximum());
  l3->SetLineWidth(2);
  l3->Draw("same");

  hv3_2->SetMarkerSize(0.5);
  hv3_2->SetMarkerStyle(20);
  hv3_2->SetLineColor(kRed);
  hv3_2->SetMarkerColor(kRed);
  T3_2->Draw(Form("(Lrb.x + ((%f/%f)*Lrb.dir.x))*1000>>V3_run_%d",targetfoils[10],BeamZDir_average,V3_r2),"","sames P");


  hv3_3->SetMarkerSize(0.5);
  hv3_3->SetMarkerStyle(20);
  hv3_3->SetLineColor(kGreen + 3);
  hv3_3->SetMarkerColor(kGreen + 3);
  T3_3->Draw(Form("(Lrb.x + ((%f/%f)*Lrb.dir.x))*1000>>V3_run_%d",targetfoils[10],BeamZDir_average,V3_r3),"","sames P");

  
  Double_t xdiff_V3_1 = hv3->GetMean() - targetfoilsX[10]*1000;
  cout << "xdiff_V3_1 = " << xdiff_V3_1 << endl;

  Double_t xdiff_V3_2 = hv3_2->GetMean() - targetfoilsX[10]*1000;
  cout << "xdiff_V3_2 = " << xdiff_V3_2 << endl;

  Double_t xdiff_V3_3 = hv3_3->GetMean() - targetfoilsX[10]*1000;
  cout << "xdiff_V3_3 = " << xdiff_V3_3 << endl;

  Double_t av_xdiff_V3 = (xdiff_V3_1 + xdiff_V3_2 + xdiff_V3_3)/2.0;

  
  TLegend* leg_3 = new TLegend(.53,.65,.9,.9,"Key");
  leg_3->SetFillColor(0);
  leg_3->SetTextSize(0.02); 
  leg_3->AddEntry(hv3,Form("Beam x at V3 (run %d)",V3_r1));
  leg_3->AddEntry(hv3_2,Form("Beam x at V3 (run %d)",V3_r2));
  leg_3->AddEntry(hv3_3,Form("Beam x at V3 (run %d)",V3_r3));
  leg_3->AddEntry(l3,"V3 survey","l");
  leg_3->AddEntry((TObject*)0,Form("Mean Diff = %.3f mm",av_xdiff_V3),"");
  leg_3->AddEntry((TObject*)0,"(analyzer - survey)","");
  
  leg_3->Draw("same");

  c1->Update();

  TPaveStats *ps_V3_1 = (TPaveStats*)hv3->GetListOfFunctions()->FindObject("stats");

  ps_V3_1->SetY1NDC(0.50);
  ps_V3_1->SetY2NDC(0.63);
  ps_V3_1->SetX1NDC(0.7);
  ps_V3_1->SetX2NDC(0.9);

  TPaveStats *ps_V3_2 = (TPaveStats*)hv3_2->GetListOfFunctions()->FindObject("stats");

  ps_V3_2->SetY1NDC(0.32);
  ps_V3_2->SetY2NDC(0.45);
  ps_V3_2->SetX1NDC(0.7);
  ps_V3_2->SetX2NDC(0.9);

  TPaveStats *ps_V3_3 = (TPaveStats*)hv3_3->GetListOfFunctions()->FindObject("stats");

  ps_V3_3->SetY1NDC(0.14);
  ps_V3_3->SetY2NDC(0.27);
  ps_V3_3->SetX1NDC(0.7);
  ps_V3_3->SetX2NDC(0.9);
      
  c1->Modified();




  // create same diagrams but for BPMX values

  
  // format histograms

  TH1F* hbpm1 = new TH1F(Form("BPMB_V1_run_%d",V1_r1),"BPMB for a V1 run",1000,-3,3);
  hbpm1->GetXaxis()->SetTitle("X Position [mm]");
  hbpm1->SetMinimum(3);
  TH1F* hbpm1_2 = new TH1F(Form("BPMB_V1_run_%d",V1_r2),"BPMB for a V1 run",1000,-3,3);
  hbpm1_2->GetXaxis()->SetTitle("X Position [mm]");
  hbpm1_2->SetMinimum(3);
  TH1F* hbpm1_3 = new TH1F(Form("BPMB_V1_run_%d",V1_r3),"BPMB for a V1 run",1000,-3,3);
  hbpm1_3->GetXaxis()->SetTitle("X Position [mm]");
  hbpm1_3->SetMinimum(3);

  TH1F* hbpm2 = new TH1F(Form("BPMB_V2_run_%d",V2_r1),"BPMB for a V2 run",1000,-3,3);
  hbpm2->GetXaxis()->SetTitle("X Position [mm]");
  hbpm2->SetMinimum(3);
  TH1F* hbpm2_2 = new TH1F(Form("BPMB_V2_run_%d",V2_r2),"BPMB for a V2 run",1000,-3,3);
  hbpm2_2->GetXaxis()->SetTitle("X Position [mm]");
  hbpm2_2->SetMinimum(3);
  TH1F* hbpm2_3 = new TH1F(Form("BPMB_V2_run_%d",V2_r3),"BPMB for a V2 run",1000,-3,3);
  hbpm2_3->GetXaxis()->SetTitle("X Position [mm]");
  hbpm2_3->SetMinimum(3);

  
  TH1F* hbpm3 = new TH1F(Form("BPMB_V3_run_%d",V3_r1),"BPMB for a V3 run",1000,-3,3);
  hbpm3->GetXaxis()->SetTitle("X Position [mm]");
  hbpm3->SetMinimum(3);
  TH1F* hbpm3_2 = new TH1F(Form("BPMB_V3_run_%d",V3_r2),"BPMB for a V3 run",1000,-3,3);
  hbpm3_2->GetXaxis()->SetTitle("X Position [mm]");
  hbpm3_2->SetMinimum(3);
  TH1F* hbpm3_3 = new TH1F(Form("BPMB_V3_run_%d",V3_r3),"BPMB for a V3 run",1000,-3,3);
  hbpm3_3->GetXaxis()->SetTitle("X Position [mm]");
  hbpm3_3->SetMinimum(3);


  

  TCanvas *c2 = new TCanvas("c2","BPMB Beam x comparison",1000,1000);

  c2->Divide(3,1);
  
  c2->cd(1);

  hbpm1->SetMarkerSize(0.5);
  hbpm1->SetMarkerStyle(20);
  hbpm1->SetMarkerColor(kBlue);
  hbpm1->SetLineColor(kBlue);  
  T1->Draw(Form("Lurb.BPMB.x*1000>>BPMB_V1_run_%d",V1_r1),"","P");
  //  l1->Draw("same");

  hbpm1_2->SetMarkerSize(0.5);
  hbpm1_2->SetMarkerStyle(20);
  hbpm1_2->SetMarkerColor(kRed);
  hbpm1_2->SetLineColor(kRed);
  T1_2->Draw(Form("Lurb.BPMB.x*1000>>BPMB_V1_run_%d",V1_r2),"","sames P");


  hbpm1_3->SetMarkerSize(0.5);
  hbpm1_3->SetMarkerStyle(20);
  hbpm1_3->SetMarkerColor(kGreen + 3);
  hbpm1_3->SetLineColor(kGreen + 3);
  T1_3->Draw(Form("Lurb.BPMB.x*1000>>BPMB_V1_run_%d",V1_r3),"","sames P");

  
  TLegend* leg_bpm_1 = new TLegend(.1,.65,.47,.9,"Key");
  leg_bpm_1->SetTextSize(0.02);
  leg_bpm_1->SetFillColor(0);
  leg_bpm_1->AddEntry(hbpm1,Form("BPMB x for V1 (run %d)",V1_r1));
  leg_bpm_1->AddEntry(hbpm1_2,Form("BPMB x for V1 (run %d)",V1_r2));
  leg_bpm_1->AddEntry(hbpm1_3,Form("BPMB x for V1 (run %d)",V1_r3));
  //  leg_bpm_1->AddEntry(l1,"V1 survey","l");
  //  leg_bpm_1->AddEntry((TObject*)0,Form("Diff = %.3f mm",xdiff_1),"");
  
  leg_bpm_1->Draw("same");

  c2->Update();
  
  TPaveStats *ps_BPM_V1_1 = (TPaveStats*)hbpm1->GetListOfFunctions()->FindObject("stats");

  ps_BPM_V1_1->SetY1NDC(0.50);
  ps_BPM_V1_1->SetY2NDC(0.63);
  ps_BPM_V1_1->SetX1NDC(0.1);
  ps_BPM_V1_1->SetX2NDC(0.3);

  TPaveStats *ps_BPM_V1_2 = (TPaveStats*)hbpm1_2->GetListOfFunctions()->FindObject("stats");

  ps_BPM_V1_2->SetY1NDC(0.32);
  ps_BPM_V1_2->SetY2NDC(0.45);
  ps_BPM_V1_2->SetX1NDC(0.1);
  ps_BPM_V1_2->SetX2NDC(0.3);

  TPaveStats *ps_BPM_V1_3 = (TPaveStats*)hbpm1_3->GetListOfFunctions()->FindObject("stats");

  ps_BPM_V1_3->SetY1NDC(0.14);
  ps_BPM_V1_3->SetY2NDC(0.27);
  ps_BPM_V1_3->SetX1NDC(0.1);
  ps_BPM_V1_3->SetX2NDC(0.3);
  
  c2->Modified();

  
  c2->cd(2);

  hbpm2->SetLineColor(kBlue);
  hbpm2->SetMarkerSize(0.5);
  hbpm2->SetMarkerStyle(20);
  hbpm2->SetMarkerColor(kBlue);
  T2->Draw(Form("Lurb.BPMB.x*1000>>BPMB_V2_run_%d",V2_r1),"","P");
  //  l2->Draw("same");

  hbpm2_2->SetLineColor(kRed);
  hbpm2_2->SetMarkerSize(0.5);
  hbpm2_2->SetMarkerStyle(20);
  hbpm2_2->SetMarkerColor(kRed);
  T2_2->Draw(Form("Lurb.BPMB.x*1000>>BPMB_V2_run_%d",V2_r2),"","sames P");
  //  l2->Draw("same");

  
  hbpm2_3->SetLineColor(kGreen + 3);
  hbpm2_3->SetMarkerSize(0.5);
  hbpm2_3->SetMarkerStyle(20);
  hbpm2_3->SetMarkerColor(kGreen + 3);
  T2_3->Draw(Form("Lurb.BPMB.x*1000>>BPMB_V2_run_%d",V2_r3),"","sames P");
  

  TLegend* leg_bpm_2 = new TLegend(.53,.65,.9,.9,"Key");
  leg_bpm_2->SetTextSize(0.02);
  leg_bpm_2->SetFillColor(0);
  leg_bpm_2->AddEntry(hbpm2,Form("BPMB x for V2 (run %d)",V2_r1));
  leg_bpm_2->AddEntry(hbpm2_2,Form("BPMB x for V2 (run %d)",V2_r2));
  leg_bpm_2->AddEntry(hbpm2_3,Form("BPMB x for V2 (run %d)",V2_r3)); 
  //  leg_bpm_2->AddEntry(l2,"V2 survey","l");
  //  leg_bpm_2->AddEntry((TObject*)0,Form("Diff = %.3f mm",xdiff_2),"");
  
  leg_bpm_2->Draw("same");

  c2->Update();
  
  TPaveStats *ps_BPM_V2_1 = (TPaveStats*)hbpm2->GetListOfFunctions()->FindObject("stats");

  ps_BPM_V2_1->SetY1NDC(0.50);
  ps_BPM_V2_1->SetY2NDC(0.63);
  ps_BPM_V2_1->SetX1NDC(0.7);
  ps_BPM_V2_1->SetX2NDC(0.9);

  TPaveStats *ps_BPM_V2_2 = (TPaveStats*)hbpm2_2->GetListOfFunctions()->FindObject("stats");

  ps_BPM_V2_2->SetY1NDC(0.32);
  ps_BPM_V2_2->SetY2NDC(0.45);
  ps_BPM_V2_2->SetX1NDC(0.7);
  ps_BPM_V2_2->SetX2NDC(0.9);

  TPaveStats *ps_BPM_V2_3 = (TPaveStats*)hbpm2_3->GetListOfFunctions()->FindObject("stats");

  ps_BPM_V2_3->SetY1NDC(0.14);
  ps_BPM_V2_3->SetY2NDC(0.27);
  ps_BPM_V2_3->SetX1NDC(0.7);
  ps_BPM_V2_3->SetX2NDC(0.9);
  
  c2->Modified();

  c2->cd(3);

  hbpm3->SetLineColor(kBlue);
  hbpm3->SetMarkerSize(0.5);
  hbpm3->SetMarkerStyle(20);
  hbpm3->SetMarkerColor(kBlue);
  T3->Draw(Form("Lurb.BPMB.x*1000>>BPMB_V3_run_%d",V3_r1),"","P");
  //  l3->Draw("same");

  hbpm3_2->SetLineColor(kRed);
  hbpm3_2->SetMarkerSize(0.5);
  hbpm3_2->SetMarkerStyle(20);
  hbpm3_2->SetMarkerColor(kRed);
  T3_2->Draw(Form("Lurb.BPMB.x*1000>>BPMB_V3_run_%d",V3_r2),"","sames P");

  hbpm3_3->SetLineColor(kGreen + 3);
  hbpm3_3->SetMarkerSize(0.5);
  hbpm3_3->SetMarkerStyle(20);
  hbpm3_3->SetMarkerColor(kGreen + 3);
  T3_3->Draw(Form("Lurb.BPMB.x*1000>>BPMB_V3_run_%d",V3_r3),"","sames P");

  TLegend* leg_bpm_3 = new TLegend(.53,.65,.9,.9,"Key");
  leg_bpm_3->SetTextSize(0.02);
  leg_bpm_3->SetFillColor(0);
  leg_bpm_3->AddEntry(hbpm3,Form("BPMB x for V3 (run %d)",V3_r1));
  leg_bpm_3->AddEntry(hbpm3_2,Form("BPMB x for V3 (run %d)",V3_r2));
  leg_bpm_3->AddEntry(hbpm3_3,Form("BPMB x for V3 (run %d)",V3_r3));

  //  leg_bpm_3->AddEntry(hbpm3,"BPMB x for V3");
  
  //  leg_bpm_3->AddEntry(l3,"V3 survey","l");
  //  leg_bpm_3->AddEntry((TObject*)0,Form("Diff = %.3f mm",xdiff_3),"");
  
  leg_bpm_3->Draw("same");

  c2->Update();

  TPaveStats *ps_BPM_V3_1 = (TPaveStats*)hbpm3->GetListOfFunctions()->FindObject("stats");

  ps_BPM_V3_1->SetY1NDC(0.50);
  ps_BPM_V3_1->SetY2NDC(0.63);
  ps_BPM_V3_1->SetX1NDC(0.7);
  ps_BPM_V3_1->SetX2NDC(0.9);

  TPaveStats *ps_BPM_V3_2 = (TPaveStats*)hbpm3_2->GetListOfFunctions()->FindObject("stats");

  ps_BPM_V3_2->SetY1NDC(0.32);
  ps_BPM_V3_2->SetY2NDC(0.45);
  ps_BPM_V3_2->SetX1NDC(0.7);
  ps_BPM_V3_2->SetX2NDC(0.9);

  TPaveStats *ps_BPM_V3_3 = (TPaveStats*)hbpm3_3->GetListOfFunctions()->FindObject("stats");

  ps_BPM_V3_3->SetY1NDC(0.14);
  ps_BPM_V3_3->SetY2NDC(0.27);
  ps_BPM_V3_3->SetX1NDC(0.7);
  ps_BPM_V3_3->SetX2NDC(0.9);
  c2->Modified();
  
  
}
