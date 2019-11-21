////////////////////////////////////////////////////
//          hole_grid.C
//
//  Simple script designed to print grid of sieve holes
//  and display x,y and theta, phi distribution of 
//  chosen run and database
///////////////////////////////////////////////////





#include "../InputAPEXL.h"
#include "../Load_rootfile.C"
#include "../Load_more_rootfiles.C"
//#include "LOpticsOpt.C"
#include "../file_def.h"
#include "../APEX_Sieve.h"


#include "Load_new_replay.C"


//TCut GeneralSieveCut ="L.tr.n==1 && L.tr.chi2<0.003 && abs(L.tr.x)<0.75 && abs(L.tr.y)<0.55 && abs(L.tr.th)<0.15 && abs(L.tr.ph)<0.045";

//TCut GenrealCut = TCut("L.tr.n==1 && L.tr.tg_dp>-0.06 && L.tr.tg_dp<0.06") && GeneralSieveCut;

TCut GeneralSieveCut ="L.tr.n==1 && L.tr.chi2<0.003 && abs(L.gold.th)<0.08 && L.gold.ph>-0.07 && L.gold.ph<0.025 && abs(L.tr.r_x)<0.1 && L.vdc.u1.nclust==1 && L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.u1.nclust==1";

TCut PID_cuts = "(L.prl1.e/(L.gold.p*1000))>0.2 && ((L.prl1.e+L.prl2.e)/(L.gold.p*1000))>0.51 &&  L.cer.asum_c >400";

TCut beam_cut = "Lrb.x>0.007 && Lrb.x<0.0016";

//TCut GenrealCut = GeneralSieveCut;
TCut GenrealCut = GeneralSieveCut + PID_cuts;





int nfoil = 3;



TString CutFileName = *SoureRootFile + ".FullCut.root";
TString CutDescFileName = "./Cuts_Foil/" + RootFileName  + ".VertexCut.cut";
TString CutDescFileNameSieve = "./Cuts_Sieve/" + RootFileName + ".SieveCut.cut";
TString CutDescFileNameSieve_dir = "./Cuts_Sieve/" + RootFileName + ".SieveCut_dir.cut";
TString CutDescFileNameSieve_ell = "./Cuts_Sieve/" + RootFileName + ".SieveCut_ell.cut";

TString CutDescFileNameDp = *SoureRootFile + ".DpCut.cut";
TString CutDescFileNameCol = *SoureRootFile + ".ColCut.cut";


// sieve cut csv file name

TString Sieve_CSV_name = "sieve_csv/";





void hole_grid(Int_t runnumber, TString DB_name){



  TChain* T;  
 


  T = Load_new_replay(DB_name,runnumber);

  
  TCanvas *c1 = new TCanvas("c1","PlotSieve Angles",1000,1000);

  TH2F* h2 = new TH2F("h2", "theta_target vs. phi_target", 300, -0.04, 0.04, 300, -0.08, 0.08);
  h2->SetMinimum(2);
  
  
  c1->cd(0);

    
  T->Draw("th_tgt:ph_tgt>>h2", GenrealCut, "COLZ");
  //    T->Draw("th_tgt:ph_tgt>>h2","", "COLZ");
  c1->Update();
  



  // also draw red crosses where holes 'should' be



  
  TRotation fTCSInHCS;
  TVector3 TCSX(0,-1,0);
  TVector3 TCSZ(TMath::Sin(HRSAngle),0,TMath::Cos(HRSAngle));
  TVector3 TCSY = TCSZ.Cross(TCSX);
  fTCSInHCS.RotateAxes(TCSX,TCSY,TCSZ);
  
  fPointingOffset.SetXYZ(-MissPointZ*TMath::Sin(HRSAngle)*TMath::Cos(HRSAngle),(Double_t)MissPointY,MissPointZ*TMath::Sin(HRSAngle)*TMath::Sin(HRSAngle));
  
  
  
  Int_t FoilID = 0;
  
  
  Double_t BeamX_average = 0;
  Double_t BeamY_average = 0;
    


    
  if (FoilID == 0){
    BeamX_average = BeamX_average_V1;
    BeamY_average = BeamY_average_V1;
  }
  else if (FoilID == 1){
    BeamX_average = BeamX_average_V2;
    BeamY_average = BeamY_average_V2;
  }    
  else if (FoilID == 2){
    BeamX_average = BeamX_average_V3;
    BeamY_average = BeamY_average_V3;
  }
    


  TVector3 BeamSpotHCS_average(BeamX_average,BeamY_average,targetfoils[FoilID]);
  TVector3 BeamSpotTCS_average = fTCSInHCS.Inverse()*(BeamSpotHCS_average-fPointingOffset);
  

  
  const Double_t plotwidth = 0.0015;
    
  for(UInt_t Hole = 0; Hole < NHoles; Hole++){


    Color_t color = kBlack;
    Int_t width = 1;
      
    if (Hole == 112){
      color = kRed;
      width = 2;
    }
    if( Hole == 160){
      color = kRed;
      width = 2;
    }
    
      
    TVector3 Hole_pos = GetSieveHoleTCS(Hole);
      
    
    TVector3 MomDirectionTCS_hole = Hole_pos - BeamSpotTCS_average;
      
    Double_t theta_hole = MomDirectionTCS_hole.X()/MomDirectionTCS_hole.Z();
    Double_t phi_hole = MomDirectionTCS_hole.Y()/MomDirectionTCS_hole.Z();
      
    
      
    // + type crosses
    // TLine *lh = new TLine(posx-plotwidth,posy,posx+plotwidth,posy);
    // TLine *lv = new TLine(posx,posy-plotwidth,posx,posy+plotwidth);
    
    
    // Saltire-type crosses
    TLine *lc1 = new TLine(phi_hole-plotwidth,theta_hole-plotwidth,phi_hole+plotwidth,theta_hole+plotwidth);
    lc1->SetLineWidth(width);
    TLine *lc2 = new TLine(phi_hole-plotwidth,theta_hole+plotwidth,phi_hole+plotwidth,theta_hole-plotwidth);
    lc2->SetLineWidth(width);
      


    lc1->SetLineColor(color);
    lc2-> SetLineColor(color);
      
    lc1 -> Draw("same");
    lc2 -> Draw("same");
      
  }

  c1->Update();

  Double_t x_lim = 1.3*TMath::Max(TMath::Abs(SieveYbyCol[0]),TMath::Abs(SieveYbyCol[NSieveCol-1]));

  Double_t y_lim = 1.3*TMath::Max(TMath::Abs(SieveXbyRow[0]),TMath::Abs(SieveXbyRow[NSieveRow-1]));


  TCanvas *c2 = new TCanvas("c2","PlotSieve XY",1000,1000);

  TH2F* h3 = new TH2F("h3", "Sieve X vs Y", 300, -x_lim, x_lim, 300, -y_lim, y_lim);
  h3->SetMinimum(2);

  
  T->Draw("x_sieve:y_sieve>>h3", GenrealCut, "COLZ");

  
  // draw markers where sieve holes 'should' be 
  for(UInt_t Hole = 0; Hole < NHoles; Hole++){

  
    Color_t color = kBlack;
    Int_t width = 1;
    
    if (Hole == 112){
      color = kRed;
      width = 2;
    }
    if( Hole == 160){
      color = kRed;
      width = 2;
    }
    

    
    TVector3 Hole_pos = GetSieveHoleTCS(Hole);

    Double_t posy = Hole_pos.X();
    Double_t posx = Hole_pos.Y();

    TLine *lc1 = new TLine(posx-plotwidth,posy-plotwidth,posx+plotwidth,posy+plotwidth);
    lc1->SetLineWidth(width);
    TLine *lc2 = new TLine(posx-plotwidth,posy+plotwidth,posx+plotwidth,posy-plotwidth);
    lc2->SetLineWidth(width);

    lc1->SetLineColor(color);
    lc2->SetLineColor(color);

    lc1->Draw("same");
    lc2->Draw("same");

  }
  

  c2->Update();

 
}


