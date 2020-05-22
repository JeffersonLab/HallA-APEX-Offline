/*
 * Author: John Williamson
 * Date 27th April 2020
 *
 * Script designed to compare raster of both arms
 * - both arms have info for same raster so should be indentical
 */

#include "Load_more_rootfiles.C"
#include "file_def.h"



void arm_compare(){


  Int_t run = 0;
  cout << "What run number would you like to perfrom raster comparison with?    ";

  cin >> run;
  cout << endl << endl;  


  if(run<=0){
    cout << "Invalid run number. Exiting." << endl << endl;
    return;
  }

  
  TChain* T = Load_more_rootfiles(run);


  // used for determining axes size of left and right raster plots
  
  Double_t Lx_wid, Lx_mean, Ly_wid, Ly_mean, Rx_wid, Rx_mean, Ry_wid, Ry_mean;
  Double_t x_lower, x_upper, y_lower, y_upper;

  
  // used for determining axes size of difference of left and right raster plots

  Double_t Xdiff_wid, Xdiff_mean, Ydiff_wid, Ydiff_mean; 
  Double_t Xdiff_lower, Xdiff_upper, Ydiff_lower, Ydiff_upper;
  
  
  TH1F *Left_x = new TH1F("Left_x", "Raster X ", 1000, -1, 1);
  Left_x->GetXaxis()->SetTitle("rast X [m]");

    
  TH1F *Left_y = new TH1F("Left_y", "Raster Y ", 1000, -1, 1);
  Left_y->GetXaxis()->SetTitle("rast Y [m]");
  
  
  TH1F *Right_x = new TH1F("Right_x", "Raster X ", 1000, -1, 1);
  Right_x->GetXaxis()->SetTitle("rast X [m]");
  
  TH1F *Right_y = new TH1F("Right_y", "Raster Y ", 1000, -1, 1);
  Right_x->GetXaxis()->SetTitle("rast Y [m]");
  
  
  TH1F *X_diff = new TH1F("X_diff", "LHRS - RHRS Raster X", 1000,-1,1);
  X_diff->GetXaxis()->SetTitle("LHRS - RHRS rast X [m]");
  TH1F *Y_diff = new TH1F("Y_diff", "LHRS - RHRS Raster Y", 1000,-1,1);
  Y_diff->GetXaxis()->SetTitle("LHRS - RHRS rast Y [m]");



  TCanvas *X_canv = new  TCanvas("X_canv","Raster X");
  X_canv->Divide(2,1);

  X_canv->cd(1);
  Left_x->SetLineColor(kBlue);
  T->Draw("Lrb.x>>Left_x");
  Right_x->SetLineColor(kRed);
  T->Draw("Rrb.x>>Right_x","","same");

  Lx_wid = Left_x->GetRMS();
  Lx_mean = Left_x->GetMean();

  Rx_wid = Right_x->GetRMS();
  Rx_mean = Right_x->GetMean();

  if( (Rx_mean-2*Rx_wid) < (Lx_mean-2*Lx_wid)){
    x_lower = Rx_mean-2*Rx_wid;
  }
  else{
    x_lower = Lx_mean-2*Lx_wid;
  }
  
  if( (Rx_mean+2*Rx_wid) > (Lx_mean+2*Lx_wid)){
    x_upper = Rx_mean+2*Rx_wid;
  }
  else{
    x_upper = Lx_mean+2*Lx_wid;
  }

  cout << "Lx_mean = " << Lx_mean << ", Lx_wid = " << Lx_wid << endl;
  cout << "Rx_mean = " << Rx_mean << ", Rx_wid = " << Rx_wid << endl;
  cout << "x_lower = " << x_lower << ", x_upper = " << x_upper << endl;
  Left_x->GetXaxis()->SetLimits(x_lower,x_upper);
  Right_x->GetXaxis()->SetLimits(x_lower,x_upper);

  
  // TPaveStats *ps_Rx = (TPaveStats*)X_canv->GetPrimitive("stats");
  // ps_Rx->PaintBox(0.9,0.5,1.1,0.7);

  
  T->Draw("Lrb.x>>Left_x");  
  T->Draw("Rrb.x>>Right_x","","sames");

  X_canv->Update();

  TPaveStats *ps_Rx = (TPaveStats*)Right_x->GetListOfFunctions()->FindObject("stats");
  cout << "List of functions: " << Right_x->GetListOfFunctions() << endl;
  //  ps_Rx->PaintBox(0.9,0.5,1.1,0.7);
  ps_Rx->SetY1NDC(0.57);
  ps_Rx->SetY2NDC(0.73);
      
  X_canv->Modified();

  TLegend* leg_L = new TLegend(.1,.65,.37,.9,"Key");
  leg_L->SetFillColor(0);
  leg_L->AddEntry(Left_x,"LHRS rast X","l");
  leg_L->AddEntry(Right_x,"RHRS rast X","l");
  leg_L->Draw("same");

  X_canv->cd(2);

  T->Draw("Lrb.x-Rrb.x>>X_diff");

  Xdiff_wid = X_diff->GetRMS();
  Xdiff_mean = X_diff->GetMean();

  Xdiff_lower = Xdiff_mean-2*Xdiff_wid;
  Xdiff_upper = Xdiff_mean+2*Xdiff_wid;

  X_diff->GetXaxis()->SetLimits(Xdiff_lower,Xdiff_upper);
  T->Draw("Lrb.x-Rrb.x>>X_diff");
  
  X_canv->Update();

  



  // same calculations and plots as above but for y component of raster

  TCanvas *Y_canv = new  TCanvas("Y_canv","Raster Y");
  Y_canv->Divide(2,1);

  Y_canv->cd(1);
  Left_y->SetLineColor(kBlue);
  T->Draw("Lrb.y>>Left_y");
  Right_y->SetLineColor(kRed);
  T->Draw("Rrb.y>>Right_y","","same");

  Ly_wid = Left_y->GetRMS();
  Ly_mean = Left_y->GetMean();

  Ry_wid = Right_y->GetRMS();
  Ry_mean = Right_y->GetMean();

  if( (Ry_mean-2*Ry_wid) < (Ly_mean-2*Ly_wid)){
    y_lower = Ry_mean-2*Ry_wid;
  }
  else{
    y_lower = Ly_mean-2*Ly_wid;
  }
  
  if( (Ry_mean+2*Ry_wid) > (Ly_mean+2*Ly_wid)){
    y_upper = Ry_mean+2*Ry_wid;
  }
  else{
    y_upper = Ly_mean+2*Ly_wid;
  }

  cout << "Ly_mean = " << Ly_mean << ", Ly_wid = " << Ly_wid << endl;
  cout << "Ry_mean = " << Ry_mean << ", Ry_wid = " << Ry_wid << endl;
  cout << "y_lower = " << y_lower << ", y_upper = " << y_upper << endl;
  Left_y->GetXaxis()->SetLimits(y_lower,y_upper);
  Right_y->GetXaxis()->SetLimits(y_lower,y_upper);

  
  // TPaveStats *ps_Ry = (TPaveStats*)Y_canv->GetPrimitive("stats");
  // ps_Ry->PaintBox(0.9,0.5,1.1,0.7);

  
  T->Draw("Lrb.y>>Left_y");  
  T->Draw("Rrb.y>>Right_y","","sames");

  Y_canv->Update();

  TPaveStats *ps_Ry = (TPaveStats*)Right_y->GetListOfFunctions()->FindObject("stats");
  cout << "List of functions: " << Right_y->GetListOfFunctions() << endl;
  //  ps_Ry->PaintBox(0.9,0.5,1.1,0.7);
  ps_Ry->SetY1NDC(0.57);
  ps_Ry->SetY2NDC(0.73);
      
  Y_canv->Modified();

  TLegend* leg_R = new TLegend(.1,.65,.37,.9,"Key");
  leg_R->SetFillColor(0);
  leg_R->AddEntry(Left_y,"LHRS rast Y","l");
  leg_R->AddEntry(Right_y,"RHRS rast Y","l");
  leg_R->Draw("same");

  Y_canv->cd(2);

  T->Draw("Lrb.y-Rrb.y>>Y_diff");

  Ydiff_wid = Y_diff->GetRMS();
  Ydiff_mean = Y_diff->GetMean();

  Ydiff_lower = Ydiff_mean-2*Ydiff_wid;
  Ydiff_upper = Ydiff_mean+2*Ydiff_wid;

  Y_diff->GetXaxis()->SetLimits(Ydiff_lower,Ydiff_upper);
  T->Draw("Lrb.y-Rrb.y>>Y_diff");
  
  Y_canv->Update();
  

  
}

