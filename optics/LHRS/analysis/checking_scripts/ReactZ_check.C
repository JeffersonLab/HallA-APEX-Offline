#include "TString.h"
#include "../Load_more_rootfiles.C"

void ReactZ_check(Int_t run_no = 4180){

  // ----------------------------------------------------------------------
  // Script designed to print different forms of Vertex Z to display the discrepancies
  // L.tr.vz is used for making cuts in optimisation process and is obtained from THaReactionPoint.C
  // is almost identical to rpl.z (reactionpoint z)
  // referenced to here as 'AnaZ'

  // `reactz' used in optical reconstruction and referenced here as 'OpZ'
  
  // a more simple conversion from y_tg to z-vertex of y_tg/sin(theta) 
  // referred to here as `SimpZ'

  TChain* T = Load_more_rootfiles(run_no);


  //Variables to be set 
  const Double_t D2R = TMath::Pi() / 180.; // degrees to radians conversion
  const Double_t HRSAngle = 5.366 * D2R;

  
  const Double_t MissPointZ = 1.690e-3;
  const Double_t MissPointY = -1.790e-3;



  // Set variables to read from tree

  Double_t Lrb_x;      //  Lrb.x rastered beam x
  Double_t L_tr_tg_y[100];  //  L.tr.tg_y (target y variable)
  Double_t L_tr_tg_ph[100]; //  L.tr.tg_ph (target phi)
  Double_t Lrb_dir_x;  //  Lrb.dir.x (x-component of beam direction vector)				
  Double_t Lrb_dir_z;  //  Lrb.dir.z (z-component of beam direction vector)				
  Double_t L_tr_n;     //  L.tr.n    (number of tracks)
  Double_t L_tr_vz[100];    //  L.tr.vz   (track reconstructed z-vertex)


  // Set-up histograms

  TH1F *hAnaZ = new TH1F("hAna",Form("%d Z-vertex as from analyzer",run_no),100,-1,1);
  TH1F *hOpZ  = new TH1F("hOpz",Form("%d Z-vertex as from Optimisation",run_no),100,-1,1);
  TH1F *hSimpZ = new TH1F("hSimpZ",Form("%d Z-vertex as from Simple Calc",run_no),100,-1,1);


  TH1F *hDiff = new TH1F("hDiff",Form("%d Difference between analyzer and optimisation z-vertices",run_no),100,-1,1);
  TH2F *hAnaVOp = new TH2F("hAnaVOp",Form("%d analyzer z-vertex against optimisation z-vertex",run_no),100,-1,1,100,-1,1);

 

  // Define branch addresses

  T->SetBranchAddress("L.tr.n",&L_tr_n);
  T->SetBranchAddress("L.tr.tg_y",L_tr_tg_y);
  T->SetBranchAddress("L.tr.tg_ph",L_tr_tg_ph);
  T->SetBranchAddress("Lrb.x",&Lrb_x);
  T->SetBranchAddress("Lrb.dir.x",&Lrb_dir_x);
  T->SetBranchAddress("Lrb.dir.z",&Lrb_dir_z);
  T->SetBranchAddress("L.tr.vz",L_tr_vz);


  Int_t nentries = T->GetEntries();

  cout<<"Total Number of Events = "<<nentries<<endl;
  
  // Prepare other vertice variables
  
  Double_t AnaZ = 0;
  Double_t OpZ = 0;
  Double_t SimpZ = 0;

  TVector3 Tg_YSpotTCS;
  TVector3 Tg_YSpotHCS;

  TVector3 MomDirectionTCS;
  TVector3 MomDirectionHCS;



  TRotation fTCSInHCS;
  TVector3 TCSX(0,-1,0);
  TVector3 TCSZ(TMath::Sin(HRSAngle),0,TMath::Cos(HRSAngle));
  TVector3 TCSY = TCSZ.Cross(TCSX);
  fTCSInHCS.RotateAxes(TCSX,TCSY,TCSZ);

  
  TVector3 fPointingOffset;
  fPointingOffset.SetXYZ(-MissPointZ*TMath::Sin(HRSAngle)*TMath::Cos(HRSAngle),(Double_t)MissPointY,MissPointZ*TMath::Sin(HRSAngle)*TMath::Sin(HRSAngle));
  

  for(Int_t i=0;i<20000;i++){

    if(i%100000==0) cout << " events processed = " << i << endl;

    AnaZ = 0;
    SimpZ = 0;
    OpZ = 0;

    T->GetEntry(i);




    AnaZ = L_tr_vz[0];

    SimpZ = L_tr_tg_y[0] / TMath::Sin(HRSAngle);



    Tg_YSpotTCS.SetXYZ(0, L_tr_tg_y[0],0);
    Tg_YSpotHCS=fTCSInHCS*Tg_YSpotTCS+fPointingOffset;

    MomDirectionTCS.SetXYZ(0,L_tr_tg_ph[0],1);
    MomDirectionHCS=fTCSInHCS*MomDirectionTCS;


    OpZ = (Tg_YSpotHCS.X() -  Lrb_dir_x - (MomDirectionHCS.X()/MomDirectionHCS.Z())*Tg_YSpotHCS.Z() )/( Lrb_dir_x/ Lrb_dir_z - (MomDirectionHCS.X()/MomDirectionHCS.Z()) );

    hAnaZ->Fill(AnaZ);
    hOpZ->Fill(OpZ);
    hSimpZ->Fill(SimpZ);

    hDiff->Fill(AnaZ-OpZ);
    hAnaVOp->Fill(AnaZ,OpZ);

    if(i%100000==0) {
      cout << "Entry " << i << ":  AnaZ = " << AnaZ << ", SimpZ = " << SimpZ << " and OpZ = " << OpZ << endl;
    }

  }
  


  TCanvas *c1 = new TCanvas(Form("c1 for %d",run_no));
  c1->Divide(2,2);
  c1->cd(1);
  hAnaZ->Draw();
  c1->cd(2);
  hOpZ->Draw();
  c1->cd(3);
  hDiff->Draw();
  c1->cd(4);
  hAnaVOp->Draw("colz");
  //  hSimpZ->Draw();

  
  // TCanvas *c2 = new TCanvas("c2");
  // c2->Divide(2,1);
  // c2->cd(1);
  // hDiff->Draw();
  // c2->cd(2);
  // hAnaVOp->Draw("colz");




  //Double_t reactz = (Tg_YSpotHCS.X() - eventdata.Data[kBeamX] - (MomDirectionHCS.X()/MomDirectionHCS.Z())*Tg_YSpotHCS.Z() )/( eventdata.Data[kBeamDirX]/eventdata.Data[kBeamDirZ] - (MomDirectionHCS.X()/MomDirectionHCS.Z()) );

  // TVector3 TCSX(0,-1,0);
  // TVector3 TCSZ(TMath::Sin(HRSAngle),0,TMath::Cos(HRSAngle));
  // TVector3 TCSY = TCSZ.Cross(TCSX);
  // fTCSInHCS.RotateAxes(TCSX,TCSY,TCSZ);


}



