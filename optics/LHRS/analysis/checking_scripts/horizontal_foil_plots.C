// Script to plot several plots for the horizontal foils (LHRS)
//

#include "../Load_more_rootfiles.C"
#include "InputAPEXL.h"

void horizontal_foil_plots(){


  // load in root files
  
  Int_t run_nos[] = {4775,4776,4777};
  
  TChain* T;


  Int_t run_po = 0;
  for(auto run_no : run_nos){

    TChain* T_add = Load_more_rootfiles(run_no);

       
    // *T_ind[run_po] = *T_add;
    
    if(run_no == run_nos[0]){
      T = T_add;
    }
    else{
      T->Add(T_add);
    }

    cout << "T->GetEntries() = " << T->GetEntries() << endl;

    run_po++;
  }


  // set histogram limits

  // raster limits
  Double_t Lrby_l = 0.0014;
  Double_t Lrby_h = 0.0030;

  Double_t Lrbx_l = -0.0022;
  Double_t Lrbx_h = -0.0011;


  // reactz limits
  Double_t reactz_l = -0.4;
  Double_t reactz_h = 0.4;

  // Focal plane limits
  
  Double_t y_FP_l = -0.05;
  Double_t y_FP_h = 0.05;

  Double_t ph_FP_l = -50;
  Double_t ph_FP_h = 40;

  Double_t th_FP_l = -80;
  Double_t th_FP_h = 80;

  // target cuts

  Double_t ph_tg_l = -40;
  Double_t ph_tg_h = 40;
  
  Double_t th_tg_l = -80;
  Double_t th_tg_h = 80;
  
  // define raster cuts
  Double_t c_Lrby1_l = 0.00191;
  Double_t c_Lrby1_h = 0.00209;

  Double_t c_Lrby2_l = 0.00240;
  Double_t c_Lrby2_h = 0.00258;

  
  // define reactz cut for V2
  Double_t c_zV2_l = -0.1;
  Double_t c_zV2_h = 0.05;

  // define V2 and raster cuts

  TCut V2_zcut = Form("reactz>%f && reactz<%f",c_zV2_l,c_zV2_h);

  // lower raster cut (in y)
  TCut rast_l = Form("Lrb.y>%f && Lrb.y<%f",c_Lrby1_l,c_Lrby1_h);

  // greater raster cut (in y)
  TCut rast_h = Form("Lrb.y>%f && Lrb.y<%f",c_Lrby2_l,c_Lrby2_h);

  // cuts from InputAPEXl.h
  TCut GenrealCut = GeneralSieveCut + PID_cuts + FP_cuts;
  
  TCanvas *c1 = new TCanvas("c1"); 

  cout << GenrealCut << endl << endl;
  
  TH2F* h_raster = new TH2F("h_raster","Raster", 400, Lrbx_l, Lrbx_h, 400, Lrby_l, Lrby_h);
  h_raster->GetYaxis()->SetTitle("Raster y [m]");
  h_raster->GetXaxis()->SetTitle("Raster x [m]");
  
  T->Draw("Lrb.y:Lrb.x>>h_raster",GenrealCut,"colz");

  // draw raster cuts
  TLine* rast1_l = new TLine(Lrbx_l,c_Lrby1_l,Lrbx_h,c_Lrby1_l);
  rast1_l->SetLineWidth(1.5);
  rast1_l->SetLineColor(kRed);
  TLine* rast1_h = new TLine(Lrbx_l,c_Lrby1_h,Lrbx_h,c_Lrby1_h);
  rast1_h->SetLineWidth(1.5);
  rast1_h->SetLineColor(kRed);

  TLine* rast2_l = new TLine(Lrbx_l,c_Lrby2_l,Lrbx_h,c_Lrby2_l);
  rast2_l->SetLineWidth(1.5);
  rast2_l->SetLineColor(kRed);
  TLine* rast2_h = new TLine(Lrbx_l,c_Lrby2_h,Lrbx_h,c_Lrby2_h);
  rast2_h->SetLineWidth(1.5);
  rast2_h->SetLineColor(kRed);


  rast1_l->Draw("same");
  rast1_h->Draw("same");
  rast2_l->Draw("same");
  rast2_h->Draw("same");

  TCanvas *c2 = new TCanvas("c2");

  TH2F* h_z_rasty = new TH2F("h_z_rasty","Reactz vs. raster y",400,Lrby_l,Lrby_h,400,reactz_l,reactz_h);
  h_z_rasty->GetYaxis()->SetTitle("Reactz [m]");
  h_z_rasty->GetXaxis()->SetTitle("Raster y [m]");

  T->Draw("reactz:Lrb.y>>h_z_rasty",GenrealCut,"colz");

  TLine* V2_rz_l = new TLine(Lrby_l,c_zV2_l,Lrby_h,c_zV2_l);
  V2_rz_l->SetLineWidth(1.5);
  V2_rz_l->SetLineColor(kRed);
  TLine* V2_rz_h = new TLine(Lrby_l,c_zV2_h,Lrby_h,c_zV2_h);
  V2_rz_h->SetLineWidth(1.5);
  V2_rz_h->SetLineColor(kRed);

  V2_rz_l->Draw("same");
  V2_rz_h->Draw("same");
  

  TCanvas *c3 = new TCanvas("c3");

  TH2F* h_z_rastx = new TH2F("h_z_rastx","Reactz vs. raster x",400,Lrbx_l,Lrbx_h,400,reactz_l,reactz_h);
  h_z_rastx->GetYaxis()->SetTitle("Reactz [m]");
  h_z_rastx->GetXaxis()->SetTitle("Raster x [m]");

  T->Draw("reactz:Lrb.x>>h_z_rastx",GenrealCut,"colz");


  TLine* V2_rzx_l = new TLine(Lrbx_l,c_zV2_l,Lrbx_h,c_zV2_l);
  V2_rzx_l->SetLineWidth(1.5);
  V2_rzx_l->SetLineColor(kRed);
  TLine* V2_rzx_h = new TLine(Lrbx_l,c_zV2_h,Lrbx_h,c_zV2_h);
  V2_rzx_h->SetLineWidth(1.5);
  V2_rzx_h->SetLineColor(kRed);
  
  V2_rzx_l->Draw("same");
  V2_rzx_h->Draw("same");

  // Focal plane plots

  TCanvas *c4 = new TCanvas("c4");

  TH2F* h_FP_py_1 = new TH2F("h_FP_py_1","Focal Plane #phi vs Y, lower beam Y",400,y_FP_l,y_FP_h,400,ph_FP_l,ph_FP_h);
  h_FP_py_1->GetXaxis()->SetTitle("FP Y [m]");
  h_FP_py_1->GetYaxis()->SetTitle("FP #phi [mrad]");
  
  TH2F* h_FP_py_2 = new TH2F("h_FP_py_2","Focal Plane #phi vs Y, greater beam Y",400,y_FP_l,y_FP_h,400,ph_FP_l,ph_FP_h);
  h_FP_py_2->GetXaxis()->SetTitle("FP Y [m]");
  h_FP_py_2->GetYaxis()->SetTitle("FP #phi [mrad]");

  TH2F* h_FP_py_3 = new TH2F("h_FP_py_3","Focal Plane #phi vs Y, V2 cut",400,y_FP_l,y_FP_h,400,ph_FP_l,ph_FP_h);
  h_FP_py_3->GetXaxis()->SetTitle("FP Y [m]");
  h_FP_py_3->GetYaxis()->SetTitle("FP #phi [mrad]");

  c4->Divide(3,1);
  c4->cd(1);

  cout << "rast_l + !rast_h = " << rast_l + !rast_h << endl;


  T->Draw("1000*L.tr.r_ph:L.tr.r_y>>h_FP_py_1",GenrealCut + rast_l + !rast_h + !V2_zcut,"colz");  
  c4->cd(2);  
  T->Draw("1000*L.tr.r_ph:L.tr.r_y>>h_FP_py_2",GenrealCut + !rast_l + rast_h + !V2_zcut,"colz");
  c4->cd(3);
  T->Draw("1000*L.tr.r_ph:L.tr.r_y>>h_FP_py_3",GenrealCut + !rast_l + !rast_h + V2_zcut,"colz");

  TCanvas *c5 = new TCanvas("c5");

  TH2F* h_FP_ty_1 = new TH2F("h_FP_ty_1","Focal Plane #theta vs Y, lower beam Y",400,y_FP_l,y_FP_h,400,th_FP_l,th_FP_h);
  h_FP_ty_1->GetXaxis()->SetTitle("FP Y [m]");
  h_FP_ty_1->GetYaxis()->SetTitle("FP #theta [mrad]");
  
  TH2F* h_FP_ty_2 = new TH2F("h_FP_ty_2","Focal Plane #theta vs Y, greater beam Y",400,y_FP_l,y_FP_h,400,th_FP_l,th_FP_h);
  h_FP_ty_2->GetXaxis()->SetTitle("FP Y [m]");
  h_FP_ty_2->GetYaxis()->SetTitle("FP #theta [mrad]");

  TH2F* h_FP_ty_3 = new TH2F("h_FP_ty_3","Focal Plane #theta vs Y, V2 cut",400,y_FP_l,y_FP_h,400,th_FP_l,th_FP_h);
  h_FP_ty_3->GetXaxis()->SetTitle("FP Y [m]");
  h_FP_ty_3->GetYaxis()->SetTitle("FP #theta [mrad]");

  c5->Divide(3,1);
  c5->cd(1);
  T->Draw("1000*L.tr.r_th:L.tr.r_y>>h_FP_ty_1",GenrealCut + rast_l + !rast_h + !V2_zcut,"colz");
  c5->cd(2);
  T->Draw("1000*L.tr.r_th:L.tr.r_y>>h_FP_ty_2",GenrealCut + !rast_l + rast_h + !V2_zcut,"colz");
  c5->cd(3);
  T->Draw("1000*L.tr.r_th:L.tr.r_y>>h_FP_ty_3",GenrealCut + !rast_l + !rast_h + V2_zcut,"colz");

  // target sieve distributions

  TCanvas *c6 = new TCanvas("c6");

  TH2F* h_tg_tp_1 = new TH2F("h_tg_tp_1","Target #theta vs #phi, lower beam Y",400,ph_tg_l,ph_tg_h,400,th_tg_l,th_tg_h);
  h_tg_tp_1->GetXaxis()->SetTitle("tg #phi [mrad]");
  h_tg_tp_1->GetYaxis()->SetTitle("tg #theta [mrad]");

  TH2F* h_tg_tp_2 = new TH2F("h_tg_tp_2","Target #theta vs #phi, greater beam Y",400,ph_tg_l,ph_tg_h,400,th_tg_l,th_tg_h);
  h_tg_tp_2->GetXaxis()->SetTitle("tg #phi [mrad]");
  h_tg_tp_2->GetYaxis()->SetTitle("tg #theta [mrad]");

  TH2F* h_tg_tp_3 = new TH2F("h_tg_tp_3","Target #theta vs #phi, V2 cut",400,ph_tg_l,ph_tg_h,400,th_tg_l,th_tg_h);
  h_tg_tp_3->GetXaxis()->SetTitle("tg #phi [mrad]");
  h_tg_tp_3->GetYaxis()->SetTitle("tg #theta [mrad]");

  c6->Divide(3,1);
  c6->cd(1);
  T->Draw("1000*th_tgt:1000*ph_tgt>>h_tg_tp_1",GenrealCut + rast_l + !rast_h + !V2_zcut,"colz");
  c6->cd(2);
  T->Draw("1000*th_tgt:1000*ph_tgt>>h_tg_tp_2",GenrealCut + !rast_l + rast_h + !V2_zcut,"colz");
  c6->cd(3);
  T->Draw("1000*th_tgt:1000*ph_tgt>>h_tg_tp_3",GenrealCut + !rast_l + !rast_h + V2_zcut,"colz");


  // make plot of raster current vs position (for y)

  TCanvas *c7 = new TCanvas("c7");


  TH2F* h_rcurr = new TH2F("h_rcurr","BPMP Y vs raster current",1000,0,0.007,1000,18000,45000);
  T->Draw("Lrb.Raster2.rawcur.y:Lrb.BPMA.y>>h_rcurr","","colz");

    

  // save results

  c1->Print("horizontal_foils/hor_Raster_plot.pdf");
  c2->Print("horizontal_foils/hor_Rasty_Reactz.pdf");
  c3->Print("horizontal_foils/hor_Rastx_Reactz.pdf");
  c4->Print("horizontal_foils/hor_FP_PY.pdf");
  c5->Print("horizontal_foils/hor_FP_TY.pdf");
  c6->Print("horizontal_foils/hor_tg_TP.pdf");
  
  

  
}
