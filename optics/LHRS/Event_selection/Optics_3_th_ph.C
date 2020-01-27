#include "SaveCanvas.C"
#include "TPad.h"
#include "InputAPEXL.h"
//#include "Load_rootfile.C"
#include "Load_more_rootfiles.C"
//#include "LOpticsOpt.C"
#include "file_def.h"

#include "APEX_Sieve.h"


//JW: commented here temporarily:
// const UInt_t NHoles = 225;
// const UInt_t NoinEachRow[] = {15, 12, 15, 11, 15, 11, 15, 11, 15, 11, 15, 11, 15, 11, 15, 12, 15};


//TCut GeneralSieveCut ="L.tr.n==1 && L.tr.chi2<0.003 && abs(L.tr.x)<0.75 && abs(L.tr.y)<0.55 && abs(L.tr.th)<0.15 && abs(L.tr.ph)<0.045";

//TCut GenrealCut = TCut("L.tr.n==1 && L.tr.tg_dp>-0.06 && L.tr.tg_dp<0.06") && GeneralSieveCut;

TCut GeneralSieveCut ="L.tr.n==1 && L.tr.chi2<0.003 && abs(L.gold.th)<0.08 && L.gold.ph>-0.07 && L.gold.ph<0.025 && abs(L.tr.r_x)<0.1 && L.vdc.u1.nclust==1 && L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.u1.nclust==1";

TCut PID_cuts = "(L.prl1.e/(L.gold.p*1000))>0.2 && ((L.prl1.e+L.prl2.e)/(L.gold.p*1000))>0.51 &&  L.cer.asum_c >400";


//TCut fid_cut = "abs(th_tgt)<0.03 && abs(ph_tgt)<0.02 && abs(L.tr.r_x)<0.1 ";
TCut fid_cut = "abs(L.tr.r_x)<0.1 ";

TCut beam_cut = "Lrb.x>0.007 && Lrb.x<0.0016";

//TCut GenrealCut = GeneralSieveCut;
TCut GenrealCut = GeneralSieveCut + PID_cuts + fid_cut;





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





void ReLoadcuts(){

  

  
  //  GenrealCut = TCut("abs(L.tr.tg_dp)<0.01") && GeneralSieveCut;

 
  //  GenrealCut = TCut(Form("L.tr.n==1 && L.tr.tg_dp>%.2f && L.tr.tg_dp<%.2f",dp_lim1,dp_lim2)) && GeneralSieveCut;

  //std::cout << "GenrealCut = " << GenrealCut << std::endl;

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
 
// TString SoureRootFile_dp = "/w/work3/home/johnw/Rootfiles/apex_4179_1_7_2019_dp_%d.root";
  
  // SoureRootFile_dp = 
  
 

  CutFileName = *SoureRootFile + ".FullCut.root";
  CutDescFileName = *SoureRootFile + ".VertexCut.cut";
  CutDescFileNameSieve = *SoureRootFile + ".SieveCut.cut";
  CutDescFileNameSieve_dir = *SoureRootFile + ".SieveCut_dir.cut";


  CutDescFileNameDp = *SoureRootFile + ".DpCut.cut";

  CutDescFileNameCol = *SoureRootFile + ".ColCut.cut";
  
}


void Optics_3_th_ph() {

  TChain* T = Load_more_rootfiles(Run_number, Run_number_2);

  ReLoadcuts();

  cout << "RootFileName = " << RootFileName << endl;
  cout << "Run_number = " << Run_number << endl;
  cout << "General cut = " << GenrealCut << endl << endl;

  // creating and filling expectation value of holes (theta and phi)
  Double_t sieve_ph[NSieveCol], sieve_th[NSieveRow];

  for(Int_t r_count = 0; r_count < NSieveRow; r_count++){
    sieve_th[r_count] = 0;
  }

  for(Int_t c_count = 0; c_count < NSieveCol; c_count++){
    sieve_ph[c_count] = 0;
  }




  
  
  TCanvas *c1 = new TCanvas("c1","PlotSieve",1000,1000);

  TH2F* h2 = new TH2F("h2", "theta_target vs. phi_target", 300, -0.04, 0.04, 300, -0.08, 0.08);
   h2->SetMinimum(2);
  
   





   // JW: added new canvas to create second cut for each hole on FP theta vs y.


   
    TCanvas* c2 = new TCanvas("c2","FP sieve",1000,1000);

    //    c2->Divide(2,1);
    
    
    TH2D* thfp_v_yfp = new TH2D("thfp_v_yfp","thfp_v_yfp",300,-0.05,0.05,300,-35,30);

    thfp_v_yfp->GetYaxis()->SetTitleOffset(1.0);
    thfp_v_yfp->GetXaxis()->SetTitleSize(0.05);
    thfp_v_yfp->GetYaxis()->SetTitleSize(0.05);
    thfp_v_yfp->GetXaxis()->SetTitle("y (FP) [m]");
    thfp_v_yfp->GetYaxis()->SetTitle("th (FP) [mrad]");

   
   




   // JW: add these and improve
   // TPaveText *pt1 = new TPaveText(0.12,0.75,0.32,0.88,"nbNDC");
   // pt1->AddText("Cerenkov signal sum > 500");
   // pt1->AddText("Single track");
   // pt1->AddText("Single hit in scintillator");
   // pt1->AddText("Run 4647");
   // pt1->AddText("|#delta| < 0.01");
   // pt1->SetFillColor(0);


    c1->cd(0);

    
    // T->Draw("L.tr.tg_th:L.tr.tg_ph>>h2", GenrealCut
    // 	  + TCut(Form("fcut_L_%d", FoilID)), "COLZ");

    T->Draw("th_tgt:ph_tgt>>h2", GenrealCut
	    + TCut(Form("fcut_L_%d", FoilID)), "COLZ");


    c1->Update();




  // also draw red crosses where holes 'should' be




    TRotation fTCSInHCS;
    TVector3 TCSX(0,-1,0);
    TVector3 TCSZ(TMath::Sin(HRSAngle),0,TMath::Cos(HRSAngle));
    TVector3 TCSY = TCSZ.Cross(TCSX);
    fTCSInHCS.RotateAxes(TCSX,TCSY,TCSZ);
    
    fPointingOffset.SetXYZ(-MissPointZ*TMath::Sin(HRSAngle)*TMath::Cos(HRSAngle),(Double_t)MissPointY,MissPointZ*TMath::Sin(HRSAngle)*TMath::Sin(HRSAngle));




    Double_t BeamX_average = 0;
    Double_t BeamY_average = 0;
    
    if (FoilID == 0){
      BeamX_average = BeamX_average_V1;
      BeamY_average = BeamY_average_V1;
    }
    if (FoilID == 1){
      BeamX_average = BeamX_average_V2;
      BeamY_average = BeamY_average_V2;
    }
    if (FoilID == 2){
      BeamX_average = BeamX_average_V3;
      BeamY_average = BeamY_average_V3;
    }
    
    BeamX_average = -0.0004606;
    BeamY_average = 0.002448;    


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
  

    



    // John W: attempt to add seconday cut on FP theta vs y. This acts on a secondary check on accidentaly including events wrongly identified from a different hole.


    // c2->cd(2);

    // T->Draw("L.tr.r_th*1000:L.tr.r_y>>thfp_v_yfp",GenrealCut + TCut(Form("fcut_L_%d", FoilID)),"colz");

  

}





