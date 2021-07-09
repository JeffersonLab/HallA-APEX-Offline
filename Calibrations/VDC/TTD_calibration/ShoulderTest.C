/**************************************************************
  ShoulderTest.C
  John Williamson
  April 16th 2021

  Plot timing offset against S2 coincidence times and trigger times

     
*************************************************************/



#include "Load_more_rootfiles.C"
#include "file_def.h"
#include "TTD_namespace.h"

#define NPLANE 4
#define MAX_ENTRIES 100000
#define MAX_HIT 1000

#define DEG_TO_RAD 0.017453278


void ShoulderTest(Int_t runnumber){


  TChain* T_3 = new TChain("T");
  //  T_3->Add(Form("/w/work3/home/johnw/Rootfiles/apex_online_3P_analytic_%i.root",runnumber));
  //T_3->Add(Form("/w/work3/home/johnw/Rootfiles/apex_online_3P_new_GOM_%i.root",runnumber));
  //   T_3->Add(Form("/w/work3/home/johnw/Rootfiles/apex_online_3P_Analytic_new2_%i.root",runnumber));
  //  T_3->Add(Form("/w/work3/home/johnw/Rootfiles/apex_online_2PA_Analytic_noTrig_%i.root",runnumber));
  //  T_3->Add(Form("/w/work3/home/johnw/Rootfiles/apex_online_3P_Analytic_noTrig_%i*root",runnumber));
  //  T_3->Add(Form("/w/work3/home/johnw/Rootfiles/apex_online_2PA_Analytic_Test_appRead_%i*root",runnumber));
  T_3->Add(Form("/w/work3/home/johnw/Rootfiles/apex_online_3P_NewAnalytic_Allhit_%i*root",runnumber));

  cout << "No entries in tchain = " << T_3->GetEntries() << endl;

  TCut T6Cut = "DR.evtypebits&(1<<6)";
  

  TCanvas* c_Coinc = new TCanvas("c_Coinc","c_Coinc", 1000, 1000);
  c_Coinc->Divide(2,1);

  c_Coinc->cd(1);

  Int_t coinc_bins = 300;
  Double_t coinc_start = 0;
  Double_t coinc_end = 200;
  
  TH1F *hCoinc = new TH1F("hCoinc"," S2-Time difference",coinc_bins,coinc_start,coinc_end);
  hCoinc->GetXaxis()->SetTitle("Coinc time (ns)");

  T_3->Draw("((L.s2.lt_c[L.s2.t_pads]+L.s2.rt_c[L.s2.t_pads])/2-(R.s2.lt_c[R.s2.t_pads]+R.s2.rt_c[R.s2.t_pads])/2)*1e9>>hCoinc",T6Cut,"hist");

  // S2 cut

  Double_t CutCen = 163.5;
  Double_t CutWidth = 7.5;
  
  Double_t TLowCut = 156.0;
  Double_t THighCut = 171.0;

  TLine* TLowLine = new TLine(TLowCut,0,TLowCut,1.0*hCoinc->GetMaximum());
  TLowLine->SetLineColor(kMagenta);
  TLowLine->SetLineStyle(kDashed);

  TLine* THighLine = new TLine(THighCut,0,THighCut,1.0*hCoinc->GetMaximum());
  THighLine->SetLineColor(kMagenta);
  THighLine->SetLineStyle(kDashed);

  TLowLine->Draw("same");
  THighLine->Draw("same");


  TCut S2Cut = Form("abs(((L.s2.lt_c[L.s2.t_pads]+L.s2.rt_c[L.s2.t_pads])/2-(R.s2.lt_c[R.s2.t_pads]+R.s2.rt_c[R.s2.t_pads])/2)*1e9-%f) < %f", CutCen, CutWidth);
  
  c_Coinc->cd(2);

  TH1F *hCoincC = new TH1F("hCoincC"," S2-Time difference with cut",coinc_bins,coinc_start,coinc_end);
  hCoincC->GetXaxis()->SetTitle("Coinc time (ns)");


  T_3->Draw("((L.s2.lt_c[L.s2.t_pads]+L.s2.rt_c[L.s2.t_pads])/2-(R.s2.lt_c[R.s2.t_pads]+R.s2.rt_c[R.s2.t_pads])/2)*1e9>>hCoincC",T6Cut + S2Cut);    

  TLowLine->Draw("same");
  THighLine->Draw("same");
  
  
  // LHRS time offsets
  
  TCanvas* c_Loff =  new TCanvas ("c_Loff", "Time offset (LHRS)", 1000, 1000);

  c_Loff->Divide(2,1);

  c_Loff->cd(1);
    
  TH1F* hLoff = new TH1F("hLoff","LHRS timing offset",100,-60,60);
  hLoff->GetXaxis()->SetTitle("Timing offset (ns)");

  T_3->Draw("L.vdc.u1.t0*1e9>>hLoff",T6Cut,"hist");

  c_Loff->cd(2);

  TLegend *leg_LtoffCut =  new TLegend(.1,.65,.45,.9,"Key");
  leg_LtoffCut->SetFillColor(0);
  leg_LtoffCut->SetTextSize(0.025);
  
  THStack* hLoffStack = new THStack("hLoffStack","LHRS Timing offsets ; Timing offset (ns)");

  TF1* fGaus_L; // fit 'shoulder' for LHRS

  // for gaus fits
  Int_t max_bin;
  Double_t min, max, temp;
  
  TH1F* hLoffC = new TH1F("hLoffC","LHRS timing offset with S2 cut",100,-60,60);
  hLoffC->GetXaxis()->SetTitle("time (ns)");
  hLoffC->SetLineColor(kRed);

  T_3->Draw("L.vdc.u1.t0*1e9>>hLoffC",T6Cut + S2Cut,"0");


  max_bin = hLoffC->GetMaximumBin();
  min = hLoffC->GetBinCenter(max_bin) - 30;
  max = hLoffC->GetBinCenter(max_bin) + 30;
  hLoffC->Fit("gaus","Q0","",min,max);
  
  fGaus_L = hLoffC->GetFunction("gaus");
  min = fGaus_L->GetParameter(1) - 1.5*fGaus_L->GetParameter(2);
  max = fGaus_L->GetParameter(1) + 1.5*fGaus_L->GetParameter(2);

  hLoffC->Fit("gaus","Q0","",min,max);
  fGaus_L = hLoffC->GetFunction("gaus");
  fGaus_L->SetLineColor(kBlack);  

  hLoffStack->Add(hLoffC);
  
  
  TH1F* hLoffU = new TH1F("hLoffU","LHRS timing offset without Cut",100,-60,60);
  hLoffU->SetLineColor(kGreen+2);

  T_3->Draw("L.vdc.u1.t0*1e9>>hLoffU",T6Cut + !S2Cut,"0");
  
  hLoffStack->Add(hLoffU);
  
  hLoffStack->Draw("nostack");

  leg_LtoffCut->AddEntry(hLoffC,"with S2 cut","l");
  leg_LtoffCut->AddEntry(hLoffU,"outside of S2 cut","l");
  leg_LtoffCut->Draw();


  // RHRS time offsets
  
  TCanvas* c_Roff =  new TCanvas ("c_Roff", "Time offset (RHRS)", 1000, 1000);

  c_Roff->Divide(2,1);

  c_Roff->cd(1);
    
  TH1F* hRoff = new TH1F("hRoff","RHRS timing offset",100,-60,60);
  hRoff->GetXaxis()->SetTitle("Timing offset (ns)");

  T_3->Draw("R.vdc.u1.t0*1e9>>hRoff",T6Cut,"hist");

  c_Roff->cd(2);

  TLegend *leg_RtoffCut =  new TLegend(.1,.65,.45,.9,"Key");
  leg_RtoffCut->SetFillColor(0);
  leg_RtoffCut->SetTextSize(0.025);
  
  THStack* hRoffStack = new THStack("hRoffStack","RHRS Timing offsets ; Timing offset (ns)");

  TF1* fGaus_R; // fit 'shoulder' for RHRS
  
  TH1F* hRoffC = new TH1F("hRoffC","RHRS timing offset with S2 cut",100,-60,60);
  hRoffC->GetXaxis()->SetTitle("time (ns)");
  hRoffC->SetLineColor(kRed);

  T_3->Draw("R.vdc.u1.t0*1e9>>hRoffC",T6Cut + S2Cut,"0");

  max_bin = hRoffC->GetMaximumBin();
  min = hRoffC->GetBinCenter(max_bin) - 30;
  max = hRoffC->GetBinCenter(max_bin) + 30;
  hRoffC->Fit("gaus","Q0","",min,max);
  
  fGaus_R = hRoffC->GetFunction("gaus");
  min = fGaus_R->GetParameter(1) - 1.5*fGaus_R->GetParameter(2);
  max = fGaus_R->GetParameter(1) + 1.5*fGaus_R->GetParameter(2);

  hRoffC->Fit("gaus","Q0","",min,max);
  fGaus_R = hRoffC->GetFunction("gaus");
  fGaus_R->SetLineColor(kBlack);

  
  hRoffStack->Add(hRoffC);
  
  
  TH1F* hRoffU = new TH1F("hRoffU","RHRS timing offset without Cut",100,-60,60);
  hRoffU->SetLineColor(kGreen+2);

  T_3->Draw("R.vdc.u1.t0*1e9>>hRoffU",T6Cut + !S2Cut,"0");
  
  hRoffStack->Add(hRoffU);
  
  hRoffStack->Draw("nostack");

  leg_RtoffCut->AddEntry(hRoffC,"with S2 cut","l");
  leg_RtoffCut->AddEntry(hRoffU,"outside of S2 cut","l");
  leg_RtoffCut->Draw();


  // fitted LHRS and RHRS shoulders

  TCanvas* c_LRFit =  new TCanvas ("c_LRFit", "Time offsets fitted (LHRS)", 1000, 1000);
  c_LRFit->Divide(2,1);

  c_LRFit->cd(1);
  hLoffC->Draw();
  fGaus_L->Draw("same");

  TLatex* texGausL =  new TLatex( 0.5, 0.4, Form("#splitline{#mu = %.3fs ns}{#sigma = %.3f ns}",fGaus_L->GetParameter(1),fGaus_L->GetParameter(2)));
  texGausL->SetNDC(1);
  texGausL->SetTextFont(42);
  texGausL->SetTextColor(1);
  texGausL->SetTextSize(0.065);
  texGausL->Draw();

  c_LRFit->cd(2);
  hRoffC->Draw();
  fGaus_R->Draw("same");

  TLatex* texGausR =  new TLatex( 0.1,0.4, Form("#splitline{#mu = %.3fs ns}{#sigma = %.3f ns}",fGaus_R->GetParameter(1),fGaus_R->GetParameter(2)));
  texGausR->SetNDC(1);
  texGausR->SetTextFont(42);
  texGausR->SetTextColor(1);
  texGausR->SetTextSize(0.065);
  texGausR->Draw();

  // LHRS vs RHRS toff

  TCanvas* c_LRoff = new TCanvas ("c_LRoff", "Time offset (LHRS vs RHRS)", 1000, 1000);
  
  c_LRoff->Divide(3,1);
  c_LRoff->SetLeftMargin(0.15);
  
  c_LRoff->cd(1);

  TH2F* hLRoff = new TH2F("hLRoff","LHRS vs RHRS timing offset",100,-60,60,100,-60,60);
  hLRoff->GetXaxis()->SetTitle("RHRS Timing offset (ns)");
  hLRoff->GetYaxis()->SetTitle("LHRS Timing offset (ns)");
  hLRoff->GetYaxis()->SetTitleOffset(1.2);
  
  T_3->Draw("L.vdc.u1.t0*1e9:R.vdc.u1.t0*1e9>>hLRoff","","colz");

  
  c_LRoff->cd(2);

  TH2F* hLRoffC = (TH2F*) hLRoff->Clone("hLRoffC");
  hLRoffC->SetTitle("Time offset (LHRS vs RHRS) with cut");

  T_3->Draw("L.vdc.u1.t0*1e9:R.vdc.u1.t0*1e9>>hLRoffC",T6Cut + S2Cut,"colz");
  

  c_LRoff->cd(3);

  TH2F* hLRoffU = (TH2F*) hLRoff->Clone("hLRoffU");
  hLRoffU->SetTitle("Time offset (LHRS vs RHRS) outside of cut");

  T_3->Draw("L.vdc.u1.t0*1e9:R.vdc.u1.t0*1e9>>hLRoffU",T6Cut + !S2Cut,"colz");
  
  
  // timing offset vs S2 timing coinc

  TCanvas* c_t0_S2 = new TCanvas("c_t0_s2","Time offset vs Coinc time",1000, 1000);
  c_t0_S2->Divide(2,1);

  c_t0_S2->cd(1);
  TH2F* h_t0_S2_L = new TH2F("h_t0_S2_L","LHRS timing offset vs Coinc time",120,80,200,100,-60,60);
  h_t0_S2_L->GetXaxis()->SetTitle("Coinc time (ns)");
  h_t0_S2_L->GetYaxis()->SetTitle("LHRS timing offset (ns)");
  h_t0_S2_L->GetYaxis()->SetTitleOffset(1.2);
    
  T_3->Draw("L.vdc.u1.t0*1e9:((L.s2.lt_c[L.s2.t_pads]+L.s2.rt_c[L.s2.t_pads])/2-(R.s2.lt_c[R.s2.t_pads]+R.s2.rt_c[R.s2.t_pads])/2)*1e9>>h_t0_S2_L",T6Cut,"colz");

  TLine* TLowLineVs = new TLine(TLowCut,-60,TLowCut,60);
  TLowLineVs->SetLineWidth(5);
  TLowLineVs->SetLineColor(kBlack);
  TLowLineVs->SetLineStyle(kDashed);

  TLine* THighLineVs = new TLine(THighCut,-60,THighCut,60);
  THighLineVs->SetLineWidth(5);
  THighLineVs->SetLineColor(kBlack);
  THighLineVs->SetLineStyle(kDashed);

  TLowLineVs->Draw("same");
  THighLineVs->Draw("same");


  c_t0_S2->cd(2);
  TH2F* h_t0_S2_R = new TH2F("h_t0_S2_R","RHRS timing offset vs Coinc time",120,80,200,100,-60,60);
  h_t0_S2_R->GetXaxis()->SetTitle("Coinc time (ns)");
  h_t0_S2_R->GetYaxis()->SetTitle("RHRS timing offset (ns)");
  h_t0_S2_R->GetYaxis()->SetTitleOffset(1.2);
  T_3->Draw("R.vdc.u1.t0*1e9:((L.s2.lt_c[L.s2.t_pads]+L.s2.rt_c[L.s2.t_pads])/2-(R.s2.lt_c[R.s2.t_pads]+R.s2.rt_c[R.s2.t_pads])/2)*1e9>>h_t0_S2_R",T6Cut,"colz");

  TLowLineVs->Draw("same");
  THighLineVs->Draw("same");

  
  // T2 and T5 trigger times
  
  //T->Draw("DR.rrawt2>>h1(600,1400,2000)")

    
  // Trigger time cuts

  Double_t T2LowCut = 1683;
  Double_t T2HighCut = 1710;

  Double_t T5LowCut = 1625;
  Double_t T5HighCut = 1634;




  TCanvas* c_T2 =  new TCanvas ("c_T2", "T2 trigger times", 1000, 1000);
  c_T2->Divide(2,1);

  c_T2->cd(1);
  TH1F* hT2 = new TH1F("hT2","T2 (LHRS trigger) timing",200,1600,1800);
  hT2->GetXaxis()->SetTitle("T2 time (ADC channels)");
  T_3->Draw("DR.rrawt2>>hT2",T6Cut,"hist");


  c_T2->cd(2);

  THStack* hT2Stack = new THStack("hT2Stack","T2 (LHRS trigger) timing; T2 time (ADC channels)");
  
  TH1F* hT2C = (TH1F*) hT2->Clone("hT2C");
  hT2C->SetTitle("T2 (LHRS trigger) timing Cut");
  T_3->Draw("DR.rrawt2>>hT2C",T6Cut + S2Cut,"0");
  hT2C->SetLineColor(kRed);
  hT2Stack->Add(hT2C);

  TH1F* hT2U = (TH1F*) hT2->Clone("hT2U");
  hT2U->SetTitle("T2 (LHRS trigger) timing Cut");
  T_3->Draw("DR.rrawt2>>hT2U",T6Cut + !S2Cut,"0");
  hT2U->SetLineColor(kGreen+2);
  hT2Stack->Add(hT2U);

  hT2Stack->Draw("nostack hist");


  TLine* T2LowCutL = new TLine(T2LowCut,0,T2LowCut,0.5*hT2Stack->GetMaximum());
  T2LowCutL->SetLineColor(kMagenta);
  T2LowCutL->SetLineStyle(kDashed);
  T2LowCutL->Draw();

  TLine* T2HighCutL = new TLine(T2HighCut,0,T2HighCut,0.5*hT2Stack->GetMaximum());
  T2HighCutL->SetLineColor(kMagenta);
  T2HighCutL->SetLineStyle(kDashed);
  T2HighCutL->Draw();

  TLegend *leg_T_2_5 =  new TLegend(.1,.65,.45,.9,"Key");

  leg_T_2_5->AddEntry(hRoffC,"with S2 cut","l");
  leg_T_2_5->AddEntry(hRoffU,"outside of S2 cut","l");
  leg_T_2_5->AddEntry(T2LowCutL,"Trigger timing cut","l");
  leg_T_2_5->Draw();


  
  TCanvas* c_T5 =  new TCanvas ("c_T5", "T5 trigger times", 1000, 1000);
  c_T5->Divide(2,1);

  c_T5->cd(1);
  TH1F* hT5 = new TH1F("hT5","T5 (RHRS trigger) timing",200,1600,1800);
  hT5->GetXaxis()->SetTitle("T5 time (channels)");
  T_3->Draw("DR.rrawt5>>hT5",T6Cut,"hist");


  c_T5->cd(2);

  THStack* hT5Stack = new THStack("hT5Stack","T5 (RHRS trigger) timing; T5 time (channels)");
  
  TH1F* hT5C = (TH1F*) hT5->Clone("hT5C");
  hT5C->SetTitle("T5 (LHRS trigger) timing Cut");
  T_3->Draw("DR.rrawt5>>hT5C",T6Cut + S2Cut,"0");
  hT5C->SetLineColor(kRed);
  hT5Stack->Add(hT5C);

  TH1F* hT5U = (TH1F*) hT5->Clone("hT5U");
  hT5U->SetTitle("T5 (LHRS trigger) timing Cut");
  T_3->Draw("DR.rrawt5>>hT5U",T6Cut + !S2Cut,"0");
  hT5U->SetLineColor(kGreen+2);
  hT5Stack->Add(hT5U);

  hT5Stack->Draw("nostack hist");


  TLine* T5LowCutL = new TLine(T5LowCut,0,T5LowCut,0.5*hT5Stack->GetMaximum());
  T5LowCutL->SetLineColor(kMagenta);
  T5LowCutL->SetLineStyle(kDashed);
  T5LowCutL->Draw();

  TLine* T5HighCutL = new TLine(T5HighCut,0,T5HighCut,0.5*hT5Stack->GetMaximum());
  T5HighCutL->SetLineColor(kMagenta);
  T5HighCutL->SetLineStyle(kDashed);
  T5HighCutL->Draw();


  leg_T_2_5->Draw();


  

  // T2 vs T5

  TCanvas* c_T2_T5 =  new TCanvas ("c_T2_T5", "T2 vs T5 trigger times", 1000, 1000);
  c_T2_T5->Divide(2,1);    


  c_T2_T5->cd(1);
  TH2F* h_T2_T5 = new TH2F("h_T2_T5","T2 vs T5",200,1600,1800,200,1600,1800);
  h_T2_T5->GetXaxis()->SetTitle("T5 (ADC channels)");
  h_T2_T5->GetYaxis()->SetTitle("T2 (ADC channels)");
  h_T2_T5->GetYaxis()->SetTitleOffset(1.5);
  T_3->Draw("DR.rrawt2:DR.rrawt5>>h_T2_T5",T6Cut,"colz");

  c_T2_T5->cd(2);
  TH2F* h_T2_T5C = new TH2F("h_T2_T5C","T2 vs T5 with S2 Cut",200,1600,1800,200,1600,1800);
  h_T2_T5C->GetXaxis()->SetTitle("T5 (ADC channels)");
  h_T2_T5C->GetYaxis()->SetTitle("T2 (ADC channels)");
  h_T2_T5C->GetYaxis()->SetTitleOffset(1.5);
  T_3->Draw("DR.rrawt2:DR.rrawt5>>h_T2_T5C",T6Cut + S2Cut,"colz");


  
  
  // Triggers vs coincidence time

  TCanvas* c_T_coinc = new TCanvas("c_T_coinc","Trigger vs Coinc Time", 1000, 1000);
  c_T_coinc->Divide(2,1);
  
  c_T_coinc->cd(1);
  TH2F* hT2_S2 = new TH2F("hT2_S2","T2 vs Coinc Time",120,80,200,200,1600,1800);
  hT2_S2->GetXaxis()->SetTitle("Coinc time (ns)");
  hT2_S2->GetYaxis()->SetTitle("T2 (ADC channels)");
  hT2_S2->GetYaxis()->SetTitleOffset(1.5);
  T_3->Draw("DR.rrawt2:((L.s2.lt_c[L.s2.t_pads]+L.s2.rt_c[L.s2.t_pads])/2-(R.s2.lt_c[R.s2.t_pads]+R.s2.rt_c[R.s2.t_pads])/2)*1e9>>hT2_S2",T6Cut,"colz");
  
  TLine* TLowLineTVs = new TLine(TLowCut,1600,TLowCut,1800);
  TLowLineTVs->SetLineWidth(5);
  TLowLineTVs->SetLineColor(kBlack);
  TLowLineTVs->SetLineStyle(kDashed);

  TLine* THighLineTVs = new TLine(THighCut,1600,THighCut,1800);
  THighLineTVs->SetLineWidth(5);
  THighLineTVs->SetLineColor(kBlack);
  THighLineTVs->SetLineStyle(kDashed);

  TLowLineTVs->Draw("same");
  THighLineTVs->Draw("same");

  c_T_coinc->cd(2);
  TH2F* hT5_S2 = new TH2F("hT5_S2","T5 vs Coinc Time",120,80,200,200,1600,1800);
  hT5_S2->GetXaxis()->SetTitle("Coinc time (ns)");
  hT5_S2->GetYaxis()->SetTitle("T5 (ADC channels)");
  hT5_S2->GetYaxis()->SetTitleOffset(1.4);
  T_3->Draw("DR.rrawt5:((L.s2.lt_c[L.s2.t_pads]+L.s2.rt_c[L.s2.t_pads])/2-(R.s2.lt_c[R.s2.t_pads]+R.s2.rt_c[R.s2.t_pads])/2)*1e9>>hT5_S2",T6Cut,"colz");
  
  TLowLineTVs->Draw("same");
  THighLineTVs->Draw("same");

  
  TCanvas* c_toff_t = new TCanvas("c_toff_t","Timing offset vs trigger times", 1000, 1000);
  c_toff_t->Divide(2,1);

  c_toff_t->cd(1);
  TH2F* h_t0L_T2 = new TH2F("h_t0L_T2","LHRS timing offset vs T2 time",200,1600,1800,100,-60,60);
  h_t0L_T2->GetXaxis()->SetTitle("T2 (ADC channels)");
  h_t0L_T2->GetYaxis()->SetTitle("LHRS timing offset (ns)");
  h_t0L_T2->GetYaxis()->SetTitleOffset(1.2);
  T_3->Draw("L.vdc.u1.t0*1e9:DR.rrawt2>>h_t0L_T2",T6Cut,"colz");

  c_toff_t->cd(2);
  TH2F* h_t0L_T2C = new TH2F("h_t0L_T2C","LHRS timing offset vs T2 time (Cut on Coinc time)",200,1600,1800,100,-60,60);
  h_t0L_T2C->GetXaxis()->SetTitle("T2 (ADC channels)");
  h_t0L_T2C->GetYaxis()->SetTitle("LHRS timing offset (ns)");
  h_t0L_T2C->GetYaxis()->SetTitleOffset(1.2);
  T_3->Draw("L.vdc.u1.t0*1e9:DR.rrawt2>>h_t0L_T2C",T6Cut + S2Cut,"colz");



  // save chosen plots

  c_Coinc->SaveAs(Form("plots/shoulder/%i_coinc.png",runnumber));
  c_Loff->SaveAs(Form("plots/shoulder/%i_Loff.png",runnumber));
  c_Roff->SaveAs(Form("plots/shoulder/%i_Roff.png",runnumber));
  c_LRoff->SaveAs(Form("plots/shoulder/%i_LRoff.png",runnumber));
  c_t0_S2->SaveAs(Form("plots/shoulder/%i_c_t0_S2.png",runnumber));
  c_T2->SaveAs(Form("plots/shoulder/%i_c_T2.png",runnumber));
  c_T5->SaveAs(Form("plots/shoulder/%i_c_T5.png",runnumber));
  c_T2_T5->SaveAs(Form("plots/shoulder/%i_c_T2_T5.png",runnumber));
  c_T_coinc->SaveAs(Form("plots/shoulder/%i_c_T_coinc.png",runnumber));
  c_toff_t->SaveAs(Form("plots/shoulder/%i_c_toff_t.png",runnumber));
  

}
