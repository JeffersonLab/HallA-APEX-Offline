#include "InputR.h"
#include "SaveCanvas.C"
#include "TH1.h"
#include "TH2.h"
#include "TCut.h"
#include "TSpectrum.h"

//////////////////////////////////////////////////////////////////////////////
// Work Directory
//////////////////////////////////////////////////////////////////////////////
TString WorkDir = "/home/sean/Grad/Research/APEX/optimization/Sieve/test/";



//TString WorkDir = "./Dp1/";
//TString PlotDir = "./Plot1/";

//////////////////////////////////////////////////////////////////////////////
// Define Files
//////////////////////////////////////////////////////////////////////////////
TString RootDir = "/home/sean/Grad/Research/APEX/Rootfiles/";

// Dp & Sieve
TString RootFile_Dp_m4 = "right_gmp_22771.root";
TString RootFile_Dp_m3;
TString RootFile_Dp_m2 = "right_gmp_22772.root";
TString RootFile_Dp_m1;
TString RootFile_Dp_0 = "right_gmp_22775.root";
TString RootFile_Dp_p1;
TString RootFile_Dp_p2 = "right_gmp_22776.root";
TString RootFile_Dp_p3;
string RootFile_Dp_p4 = "right_gmp_22778.root" ;
string RootFile_cen = "apex_4647_opt_sieve_plane_5th.root" ;
string RootFile_up = "apex_4648.root apex_4648_1.root apex_4648_2.root" ;
string RootFile_dn = "apex_4650.root apex_4650_1.root apex_4650_2.root" ;

 
TString RootFile_Sieve = "right_gmp_22828.root right_gmp_22829.root right_gmp_22830.root right_gmp_22831.root right_gmp_22832.root right_gmp_22833.root right_gmp_22834.root";


string RootFile_Dp = RootFile_Dp_p4;

string RootFile_test = RootFile_cen;

//TString RootFile_Vertex = RootFile_Sieve;

//TString SourceRootFile = RootFile_Sieve;
//TString SourceRootFile = RootFile_Vertex;
string SourceRootFile = RootFile_test;

//////////////////////////////////////////////////////////////////////////////
// Define Cuts
//////////////////////////////////////////////////////////////////////////////
//Bool_t UseVertexCut = kTRUE;
Bool_t UseVertexCut = kFALSE;
//Bool_t UseDpCut = kTRUE;
Bool_t UseDpCut = kFALSE;
Bool_t UseHDpCut = kFALSE;

// Switch for Plot Cuts
// 0  : general cuts
// 1  : foil cuts
// 2  : dp cuts
// 3  : dp cuts + foil cuts

UInt_t PlotCut = 0;

//TCut GeneralSieveCut = "abs(R.tr.tg_th)<0.15 && abs(R.tr.tg_ph)<0.15&& R.vdc.u1.nclust==1&& R.vdc.v1.nclust==1&& R.vdc.u2.nclust==1&& R.vdc.v2.nclust==1"; // && urb.y<0.006

//TCut GeneralCut = "R.tr.n==1&&(R.sh.e+R.ps.e)/(1000*R.tr.p)>0.6&&(R.cer.asum_c>600) ";
TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && abs(R.tr.r_x) < 0.10";
//TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > -0.45 && R.tr.r_x < -0.25)";
//TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > 0.25 && R.tr.r_x < 0.45)";

//////////////////////////////////////////////////////////////////////////////                    
// Settings                                                                                      
//////////////////////////////////////////////////////////////////////////////                   

      
Float_t dplowlimit = -0.12, dpuplimit = 0.1;
Float_t thlowlimit = -0.07, thuplimit = 0.06;
Float_t phlowlimit = -0.03, phuplimit = 0.03;
Float_t vzlowlimit = -0.2, vzuplimit = 0.15;

int nFoil = 1;

TString RootFileName;
TString CutSuf = ".FullCut.root";
TString CutDescFileSufVertex = ".VertexCut.cut";
TString CutDescFileSufDp = ".DpCut.cut";
TString CutDescFileSufSieve = ".SieveCut.%d_%d.cut";
TString CutDat = ".cuts.csv";




//////////////////////////////////////////////////////////////////////////////                   
TChain* LoadRootFiles()
{
  gStyle->SetOptStat("ne");

  double xuplimit, xlowlimit;

  if(SourceRootFile == RootFile_cen){
    dplowlimit=-0.05;dpuplimit=0.06;
    vzlowlimit=-0.2;vzuplimit=0.15;
    phlowlimit=-0.065;phuplimit=0.065;
    thlowlimit=-0.065;thuplimit=0.065;
  }
  if(SourceRootFile == RootFile_up){
    dplowlimit=-0.05;dpuplimit=0.06;
    vzlowlimit=-0.1;vzuplimit=0.25;
    phlowlimit=-0.01;phuplimit=0.035;
    thlowlimit=-0.08;thuplimit=0.05;
  }
  if(SourceRootFile == RootFile_dn){
    dplowlimit=-0.05;dpuplimit=0.06;
    vzlowlimit=-0.6;vzuplimit=0.0;
    phlowlimit=-0.04;phuplimit=0.01;
    thlowlimit=-0.08;thuplimit=0.07;
  }
  if(SourceRootFile == RootFile_Dp_m4){
    dplowlimit=-0.05;dpuplimit=-0.01;
  }
  if(SourceRootFile == RootFile_Dp_m3){
    dplowlimit=-0.04;dpuplimit=-0.025;
    xlowlimit=-0.14-0.40;xuplimit=0.06-0.40;
  }
  if(SourceRootFile == RootFile_Dp_m2){
    dplowlimit=-0.05;dpuplimit=0.01;
    //        xlowlimit=-0.14-0.24;xuplimit=0.06-0.24;                                       
  }
  if(SourceRootFile == RootFile_Dp_m1){
    dplowlimit=-0.02;dpuplimit=-0.005;
    xlowlimit=-0.14-0.14;xuplimit=0.06-0.14;
  }
  if(SourceRootFile == RootFile_Dp_0){
    dplowlimit=-0.04;dpuplimit=0.03;
  }
  if(SourceRootFile == RootFile_Dp_p1){
    dplowlimit=0.00;dpuplimit=0.015;
    xlowlimit=-0.14+0.13;xuplimit=0.06+0.13;
  }
  if(SourceRootFile == RootFile_Dp_p2){
    dplowlimit=-0.02;dpuplimit=0.04;
    //        xlowlimit=-0.14+0.26;xuplimit=0.06+0.26;                                       
  }
  if(SourceRootFile == RootFile_Dp_p3){
    dplowlimit=0.02;dpuplimit=0.035;
    xlowlimit=-0.14+0.38;xuplimit=0.06+0.38;
  }
  if(SourceRootFile == RootFile_Dp_p4){
    dplowlimit=0.0;dpuplimit=0.055;
  }

  //char s[1000] = SourceRootFile.Data();
  char s[1000];
  strcpy(s,SourceRootFile.c_str());
  const char* d = " ";
  char* sub;
  TList Files;

  sub=strtok(s,d);
  if(sub){
    RootFileName = TString(sub);
  }
  while(sub){
    Files.Add(new TObjString(sub));
    sub=strtok(NULL,d);
  }

  TIter next(&Files);
  TObjString* rootfile;
  TChain *T = new TChain("T");
  Int_t n=0;

  while((rootfile = (TObjString*)next())){
    T->Add(RootDir + rootfile->GetString());
    cout << "Open " <<rootfile->GetString() << endl;
  }

  return T;
}






//////////////////////////////////////////////////////////////////////////////
// Cut Tools
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
void CutVertex(int overwrite = 0) {


  TChain *T = LoadRootFiles();


    TString CutFileName = WorkDir + RootFileName + CutSuf;
    TString CutDescName = WorkDir + RootFileName + CutDescFileSufVertex;

    cerr << "cp -vf " + CutFileName + " " + CutFileName + ".old" << endl;
    gSystem->Exec("cp -vf " + CutFileName + " " + CutFileName + ".old");

    fstream cutdesc(CutDescName, ios_base::out);
    assert(cutdesc.is_open());

    TFile *f1 = new TFile(CutFileName, "UPDATE");
    assert(f1);

    // Define Canvas
    TCanvas* c3 = new TCanvas("c3", "ReactZ Cuts", 900, 900);
    TH2F* h3 = new TH2F("h3", "Tg_Z vs. Tg ph", 400, phlowlimit, phuplimit, 400, vzlowlimit, vzuplimit);
    //    TH2F* h3 = new TH2F("h3", "Tg_ph vs. R.tr.vz",  400, vzlowlimit, vzuplimit, 400, phlowlimit, phuplimit);
    assert(h3);

    // Choose the foil you want to make cut
    //nFoil
    for (int FoilID = 1; FoilID < 2; FoilID++) {
        TCut DrawCut = GeneralCut;

	T->Draw("R.tr.vz:R.tr.tg_ph>>h3", DrawCut, "COLZ");
	//	T->Draw("R.tr.tg_ph:R.tr.vz>>h3", DrawCut, "COLZ");
	c3->Update();
        
        cout << "Testing " << Form("fcut_R_%d", FoilID) << endl;
        TCutG* cutg = (TCutG*)gROOT->FindObject(Form("fcut_R_%d", FoilID)); //looking for old cut definition
        
        if(!cutg || overwrite){
            cout << "Making cut " << Form("fcut_R_%d", FoilID) << endl;

            cutg = (TCutG*)(TVirtualPad::Pad()->WaitPrimitive("CUTG", "CutG")); // making cut, store to CUTG
            c3->Update();
	    
            cutg->SetName(Form("fcut_R_%d", FoilID)); //
	    cutg->SetVarX("R.tr.tg_ph");
	    cutg->SetVarY("R.tr.vz");
	    //	    cutg->SetVarY("R.tr.tg_ph");
	    //           cutg->SetVarX("R.tr.vz");

            cout << "done!" << endl;
        }
        else{
            cout << Form("fcut_R_%d", FoilID) << " is found, using old one" << endl;
	    }

        cutg->SetLineColor(kMagenta);
        cutg->SetLineWidth(2);
        cutg->Draw("PL");
        c3->Update();

        // output cut to disk
        cutg->Write("", TObject::kOverwrite); // Overwrite old cut

        cout << "Log to " << CutDescName << endl;

        cutdesc << Form("fcut_R_%d", FoilID) << " && " << (const char*)GeneralCut << endl;

        SaveCanvas(c3, WorkDir + RootFileName + Form(".fcut_R_%d", FoilID), kFALSE);
    }

    f1->Write();
    f1->ls();
    cutdesc.close();
}
/////////////////////////////////////////////////////////////////////////////
void CutVertex_Dp(int overwrite = 0) {


  TChain *T = LoadRootFiles();
  

    TString CutFileName = WorkDir + RootFileName + CutSuf;
    TString CutDescName = WorkDir + RootFileName + CutDescFileSufVertex;

    cerr << "cp -vf " + CutFileName + " " + CutFileName + ".old" << endl;
    gSystem->Exec("cp -vf " + CutFileName + " " + CutFileName + ".old");

    fstream cutdesc(CutDescName, ios_base::out);
    assert(cutdesc.is_open());

    TFile *f1 = new TFile(CutFileName, "UPDATE");
    assert(f1);

    // Define Canvas
    TCanvas* c3 = new TCanvas("c3", "ReactZ Cuts", 900, 900);
    TH2F* h3 = new TH2F("h3", "Tg_Z vs. Tg ph", 400, phlowlimit, phuplimit, 400, vzlowlimit, vzuplimit);
    //    TH1F *h3 = new TH1F("h3","Tg_Z",400,vzlowlimit,vzuplimit);
    assert(h3);

    // Choose the foil you want to make cut
    //nFoil
    for (int FoilID = 1; FoilID < 2; FoilID++) {
        TCut DrawCut = GeneralCut;

	T->Draw("R.tr.vz:R.tr.tg_ph>>h3", DrawCut, "COLZ");
	//	T->Draw("R.tr.vz>>h3", DrawCut, "COLZ");
	c3->Update();
        
        cout << "Testing " << Form("fcut_R_%d", FoilID) << endl;
        TCutG* cutg = (TCutG*)gROOT->FindObject(Form("fcut_R_%d", FoilID)); //looking for old cut definition
        
        if(!cutg || overwrite){
            cout << "Making cut " << Form("fcut_R_%d", FoilID) << endl;

            cutg = (TCutG*)(TVirtualPad::Pad()->WaitPrimitive("CUTG", "CutG")); // making cut, store to CUTG
            c3->Update();
	    
            cutg->SetName(Form("fcut_R_%d", FoilID)); //
	    cutg->SetVarX("R.tr.tg_ph");
	    cutg->SetVarY("R.tr.vz");

            cout << "done!" << endl;
        }
        else{
            cout << Form("fcut_R_%d", FoilID) << " is found, using old one" << endl;
	    }

        cutg->SetLineColor(kMagenta);
        cutg->SetLineWidth(2);
        cutg->Draw("PL");
        c3->Update();

        // output cut to disk
        cutg->Write("", TObject::kOverwrite); // Overwrite old cut

        cout << "Log to " << CutDescName << endl;

        cutdesc << Form("fcut_R_%d", FoilID) << " && " << (const char*)GeneralCut << endl;

        SaveCanvas(c3, WorkDir + RootFileName + Form(".fcut_R_%d", FoilID), kFALSE);
    }

    f1->Write();
    f1->ls();
    cutdesc.close();
}
///////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
void CutDp(int overwrite = 0)
{
  TChain *T = LoadRootFiles();
  

    TString CutFileName = WorkDir + RootFileName + CutSuf;
    TString CutDescName = WorkDir + RootFileName + CutDescFileSufDp;

    cerr << "cp -vf " + CutFileName + " " + CutFileName + ".old" << endl;
    gSystem->Exec("cp -vf " + CutFileName + " " + CutFileName + ".old");

    fstream cutdesc(CutDescName, ios_base::out);
    assert(cutdesc.is_open());

    TFile *f1 = new TFile(CutFileName, "UPDATE");
    assert(f1);

    TCanvas *c4 = new TCanvas("c4","dp Cuts",900,900);
    TH2F *h4 = new TH2F("h4", "Tg th vs. Tg dp", 400, dplowlimit, dpuplimit, 400, thlowlimit, thuplimit);
    assert(h4);

    //nFoil
    for(int FoilID = 1; FoilID < 2; FoilID++){
        TCut DrawCut = GeneralCut;
        if(UseVertexCut){
            DrawCut = TCut(Form("fcut_R_%d", FoilID));
            (TCutG*)gROOT->FindObject(Form("fcut_R_%d", FoilID));
        }
        
	T->Draw("R.tr.tg_th:R.tr.tg_dp>>h4", DrawCut, "COLZ");
        c4->Update();

        cout << "Testing " << Form("dpcut_R_%d", FoilID) << endl;
        TCutG *cutg = (TCutG*)gROOT->FindObject(Form("dpcut_R_%d", FoilID));

        if(!cutg || overwrite){
            cutg = (TCutG*)(TVirtualPad::Pad()->WaitPrimitive("CUTG", "CutG")); // making cut, store to CUTG
            c4->Update();

            cutg->SetName(Form("dpcut_R_%d", FoilID));
            // set axises' name
            cutg->SetVarX("R.tr.tg_dp");
            cutg->SetVarY("R.tr.tg_th");
            cout << "done!" << endl;
        }
        else {
            cout << Form("dpcut_R_%d", FoilID) << " is found, using old one" << endl;
        }

        cutg->SetLineColor(kMagenta);
        cutg->SetLineWidth(2);
        cutg->Draw("PL");
        c4->Update();

        cutg->Write("",TObject::kOverwrite); // Overwrite old cut

        cout << "Log to " << CutDescName << endl;

        if(UseVertexCut){
            cutdesc << Form("fcut_R_%d", FoilID) << " && ";
        }
        cutdesc << Form("dpcut_R_%d", FoilID) << " && " << (const char*)GeneralCut << endl;

        SaveCanvas(c4, WorkDir + RootFileName + Form(".dpcut_R_%d", FoilID), kFALSE);
    }

    f1->Write();
    f1->ls();
    cutdesc.close();
}

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
void CutDp_Dp(int overwrite = 0)
{
  TChain *T = LoadRootFiles();
  

    TString CutFileName = WorkDir + RootFileName + CutSuf;
    TString CutDescName = WorkDir + RootFileName + CutDescFileSufDp;

    cerr << "cp -vf " + CutFileName + " " + CutFileName + ".old" << endl;
    gSystem->Exec("cp -vf " + CutFileName + " " + CutFileName + ".old");

    fstream cutdesc(CutDescName, ios_base::out);
    assert(cutdesc.is_open());

    TFile *f1 = new TFile(CutFileName, "UPDATE");
    assert(f1);

    TCanvas *c4 = new TCanvas("c4","dp Cuts",1000,1000);
    TH2F *h4 = new TH2F("h4", "Tg ph vs. Tg dp", 400, dplowlimit, dpuplimit, 400, phlowlimit, phuplimit);
    //   TH2F *h4 = new TH2F("h4", "Tg ph vs. Tg dp", 400, dplowlimit, 0.1, 400, phlowlimit, phuplimit);
    assert(h4);

    //nFoil
    for(int FoilID = 1; FoilID < 2; FoilID++){
        TCut DrawCut = GeneralCut;
        if(UseVertexCut){
            DrawCut = TCut(Form("fcut_R_%d", FoilID));
            (TCutG*)gROOT->FindObject(Form("fcut_R_%d", FoilID));
        }
        
	T->Draw("R.tr.tg_ph:R.tr.tg_dp>>h4", DrawCut, "COLZ");
        c4->Update();

        cout << "Testing " << Form("dpcut_R_%d", FoilID) << endl;
        TCutG *cutg = (TCutG*)gROOT->FindObject(Form("dpcut_R_%d", FoilID));

        if(!cutg || overwrite){
            cutg = (TCutG*)(TVirtualPad::Pad()->WaitPrimitive("CUTG", "CutG")); // making cut, store to CUTG
            c4->Update();

            cutg->SetName(Form("dpcut_R_%d", FoilID));
            // set axises' name
            cutg->SetVarX("R.tr.tg_dp");
            cutg->SetVarY("R.tr.tg_ph");
            cout << "done!" << endl;
        }
        else {
            cout << Form("dpcut_R_%d", FoilID) << " is found, using old one" << endl;
        }

        cutg->SetLineColor(kMagenta);
        cutg->SetLineWidth(2);
        cutg->Draw("PL");
        c4->Update();

        cutg->Write("",TObject::kOverwrite); // Overwrite old cut

        cout << "Log to " << CutDescName << endl;

        if(UseVertexCut){
            cutdesc << Form("fcut_R_%d", FoilID) << " && ";
        }
        cutdesc << Form("dpcut_R_%d", FoilID) << " && " << (const char*)GeneralCut << endl;

        SaveCanvas(c4, WorkDir + RootFileName + Form(".dpcut_R_%d", FoilID), kFALSE);
    }

    f1->Write();
    f1->ls();
    cutdesc.close();
}

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
void CutSieve(int FoilID = 0, int col = 6, int overwrite = 0) {

  double sieve_ph[27], sieve_th[17];
  
  //Use for angles
  //Calculate expected positions of rows and columns
  for(int i = 0; i<27; i++){
    if(i < 25) sieve_ph[i] = (14-i)*2.99 - 8.64;
    else sieve_ph[i] = sieve_ph[24] - (i-24)*2.99*2;
  }
  
  for(int j = 0; j<NSieveRow; j++){
    sieve_th[j] = (j-8)*7.254 + 9.89;
  }
  

  /*
  //Use for sieve plane
  //Calculate expected positions of rows and columns
  for(int i = 0; i<27; i++){
    if(i < 25) sieve_ph[i] = (14-i)*0.238 - 1.094;
    else sieve_ph[i] = sieve_ph[24] - (i-24)*0.238*2;
  }
  
  for(int j = 0; j<17; j++){
    sieve_th[j] = -(j-8)*0.575 - 0.0155;
  }
  */
  
  
  TChain *T = LoadRootFiles();
  

    TString CutFileName = WorkDir + RootFileName + CutSuf;
    TString TempString(Form(CutDescFileSufSieve.Data(), FoilID, col));
    TString PlotDir(RootFileName + Form(".hcut_R_%d_%d/", FoilID, col));
    TString CutDescName = WorkDir + RootFileName + TempString;
    TString CutDataName = WorkDir + RootFileName + CutDat;

    cerr << "cp -vf " + CutFileName + " " + CutFileName + ".old" << endl;
    gSystem->Exec("cp -vf " + CutFileName + " " + CutFileName + ".old");
    gSystem->Exec("mkdir -vp " + WorkDir + PlotDir);

    fstream cutdesc(CutDescName, ios_base::out);
    assert(cutdesc.is_open());

    //The following code creates a csv file that has all of the cut information
    ifstream cutcsvtest;
    cutcsvtest.open(WorkDir + "apex_4647.root.cuts_full.csv");

    TDatime* date = new TDatime();  //Get Current date

    //If file does not exists then create a new template
    if(!cutcsvtest.good()){

      ofstream createcsv;
      createcsv<<fixed<<setprecision(2);
      createcsv.open(WorkDir + "apex_4647.root.cuts_full.csv");
      createcsv<<date->GetDay()<<"/"<<date->GetMonth()<<"/"<<date->GetYear()<<" (dd/mm/yyyy)"<<endl;
      createcsv << "Hole ID (col:row),Hole Exists,Included in opt,Ellipse ph cen,Expected ph,Ellipse th cen,Expected th,Semi axis ph,Semi axis th,Ellipse Tilt (deg),Ellipse ph rms,Ellipse th rms,Statistics,All angles in mrad except ellipse tilt which is positive counterclockwise from vertical axis"<<endl;

      //Add expected positions to csv file
    for(int i = 0; i < 27; i++)
      for(int j = 0; j < 17; j++) 
	createcsv<<i<<":"<<j<<",0,0,0,"<<sieve_ph[i]<<",0,"<<sieve_th[j]<<",0,0,0"<<endl;
    
      
      createcsv.close();
    }
    
    cutcsvtest.close();

    //Open old file
    ifstream cutcsvold;
    cutcsvold.open(WorkDir + "apex_4647.root.cuts_full.csv");

    //Open file where new variables will be written
    ofstream cutcsv;
    cutcsv.open(CutDataName);
    
    //Set output to only 2 digits
    cutcsv<<fixed<<setprecision(2);
    cout<<fixed<<setprecision(2);
    
    
    cutcsv<<date->GetDay()<<"/"<<date->GetMonth()<<"/"<<date->GetYear()<<" (dd/mm/yyyy)"<<endl;
    cutcsv << "Hole ID (col:row),Hole Exists,Included in opt,Ellipse ph cen,Expected ph,Ellipse th cen,Expected th,Semi axis ph,Semi axis th,Ellipse Tilt (deg),Ellipse ph rms,Ellipse th rms,Statistics,All angles in mrad except ellipse tilt which is positive counterclockwise from vertical axis"<<endl;
    
    TFile *f1 = new TFile(CutFileName, "UPDATE");
    assert(f1);

    
    TCanvas *c5 = new TCanvas("c5","PlotSieve",1000,1000);
    TH2F *h5 = new TH2F("h5", "Tg th vs. Tg ph", 400, phlowlimit*1000, phuplimit*1000, 400, thlowlimit*1000, thuplimit*1000);  //Use for angles opt
    //TH2F *h5 = new TH2F("h5", "Sieve Plane x vs y", 400, -5, 5, 400, -6, 4);  //Use for phi opt
    //h5->SetTitle("Sieve Plane x vs y;y (cm);x (cm)");
    h5->SetTitle("Tg th vs. Tg ph;#phi (mrad);#theta (mrad)");
    h5->GetZaxis()->SetRangeUser(0,50);
    assert(h5);
    
    TPaveText *pt1 = new TPaveText(0.12,0.75,0.32,0.88,"nbNDC");
    pt1->AddText("Run 4647");
    pt1->AddText("Cerenkov signal sum > 500");
    pt1->AddText("Single track");
    pt1->AddText("|x_{fp}| < 10 cm ");
    pt1->SetFillColor(0);

    TText *text = pt1->GetLineWith("Run");
    text->SetTextColor(kRed);
    text->SetTextFont(23);
    text->SetTextSize(23);
    
    

    TCut DrawCut = GeneralCut;
    if(UseVertexCut){
        DrawCut = DrawCut + TCut(Form("fcut_R_%d", FoilID));
        (TCutG*)gROOT->FindObject(Form("fcut_R_%d", FoilID));
    }
    if(UseDpCut){
        DrawCut = DrawCut + TCut(Form("dpcut_R_%d", FoilID));
        (TCutG*)gROOT->FindObject(Form("dpcut_R_%d", FoilID));
    }

    T->Draw("R.tr.tg_th*1000:R.tr.tg_ph*1000>>h5", DrawCut, "COLZ"); //Angles cut
    //T->Draw("Sieve.x*100:Sieve.y*100>>h5", DrawCut, "COLZ");  //Sieve plane cut
    pt1->Draw("same");
    c5->Update();
    //c5->Modified();

    //    TCanvas *cth = new TCanvas("cth","cth");
    //    T->Draw("L.tr.tg_th>>hth", DrawCut);

    //    TCanvas *cph = new TCanvas("cph","cph");
    //    T->Draw("L.tr.tg_ph>>hph", DrawCut);

    int nhol = 0;
    cout << "How many holes in this No." << col << " column?" << endl;
    cin >> nhol;
    if(nhol < 0)return;
    nhol = nhol*2; //double because of row holes that don't exists
    cout << "min hole id : ";
    int rmin = 0;    
    cin >> rmin;
    if(rmin < 0)return;

    Double_t ph_peak[7]={.0, 2.e-2, 1.e-2, -0.03e-2, -1.07e-2, -2.08e-2};
    Double_t th_peak[7]={.0, -4.01e-2, -1.93e-2, 0.3e-2, 2.52e-2, 4.55e-2};
    Double_t ph_width = 1.1*0.003;
    Double_t th_width = 1.3*0.005;

    string line;
    getline(cutcsvold,line);
    getline(cutcsvold,line);
    for(int i = 0; i<col; i++){
      for(int j = 0; j<NSieveRow;j++){
	getline(cutcsvold,line);
	cutcsv << line << endl;
      }
    }
    
    bool hole_exist = false;

    for(int i = 0; i < rmin; i++){
      //check if holes exist or not
      //cutdesc << "fEvtHdr.fRun==0" << endl;
      cutdesc << "R.tr.n > 1000" << endl;
      if((col%2==0 && i%2==0) || (col%2==1 && i%2==1)) hole_exist = true;
      if((col%2==0 && i%2==1) || (col%2==1 && i%2==0)) hole_exist = false;
      if(col == 13 && (i != 1 || i !=15)) hole_exist = false;
      if(col == 25 && i%2==0) hole_exist = true;
      if(col == 25 && i%2==1) hole_exist = false;
      getline(cutcsvold,line);
      cutcsv << line << endl;  //Write old file for holes not cut
      //cutcsv<<col<<":"<<i<<","<<hole_exist<<",0,0,"<<sieve_ph[col]<<",0,"<<sieve_th[i]<<",0,0,0"<<endl;
	}
    //make cuts for each holes from begin number
    for(int row = rmin; row < rmin + nhol; row++){      
      //check if holes exist or not
      if((col%2==0 && row%2==0) || (col%2==1 && row%2==1)) hole_exist = true;
      if((col%2==0 && row%2==1) || (col%2==1 && row%2==0)) hole_exist = false;
      if(col == 13 && (row != 1 || row !=15)) hole_exist = false;
      if(col == 25 && row%2==0) hole_exist = true;
      if(col == 25 && row%2==1) hole_exist = false;
      if(!hole_exist){
	//cutdesc << "fEvtHdr.fRun==0" << endl;
	cutdesc << "R.tr.n > 1000" << endl;
	getline(cutcsvold,line);
	//cout<<line<<endl;
	//cutcsv << line << endl;

	//Write variables for current cuts happening 
	cutcsv<<col<<":"<<row<<","<<hole_exist<<",0,0,"<<sieve_ph[col]<<",0,"<<sieve_th[row]<<",0,0,0"<<endl;
	
	continue;
      }
      cout << "Testing " << Form("hcut_R_%d_%d_%d", FoilID, col, row) << endl;
	TCutG* cutg = (TCutG*)gROOT->FindObject(Form("hcut_R_%d_%d_%d", FoilID, col, row));

	//Create an ellipse that we will use to make our cuts
	TEllipse* Ellipse = new TEllipse(0,0,.02,.04,0,360,0);
	Ellipse->SetFillStyle(0);
	Ellipse->SetLineWidth(3);
	Ellipse->Draw("same");
	
        if(!cutg || overwrite){
            if(cutg)cutg->Delete();
	     
	    
	    //Move ellipse around and then draw any polygon using the toolbar to save the ellipse and move onto the next one
	    TGraph* g = (TGraph*)(TVirtualPad::Pad()->WaitPrimitive("Graph"));
	    
	    double kPI = 3.14159265358979323846;
	    Double_t x1 = Ellipse->GetX1(); Double_t y1 = Ellipse->GetY1();
	    Double_t r1 = Ellipse->GetR1(); Double_t r2 = Ellipse->GetR2();
	    Double_t theta = Ellipse->GetTheta();

	    //Write out ellipse values
	    getline(cutcsvold,line);	 
	    cutcsv<<col<<";"<<row<<","<<hole_exist<<",1,"<<x1<<","<<sieve_ph[col]<<","<<y1<<","<<sieve_th[row]<<","<<r1<<","<<r2<<","<<theta<<endl;
	    
	    //Write the ellipse into a bunch of points
	    const Int_t n = 200;
	    Double_t x[n], y[n];
	    Double_t circ = kPI*(r1+r2);
	    Double_t angle,dx,dy;
	    Double_t dphi = (360)*kPI/(180*n);
	    Double_t ct   = TMath::Cos(kPI*theta/180);
	    Double_t st   = TMath::Sin(kPI*theta/180);
	    for (Int_t i=0;i<n;i++) {
	      angle = Double_t(i)*dphi;
	      dx    = r1*TMath::Cos(angle);
	      dy    = r2*TMath::Sin(angle);
	      x[i]  = x1 + dx*ct - dy*st;
	      y[i]  = y1 + dx*st + dy*ct;
	      
	    }
	    g->Delete();
	    Ellipse->Delete();
	    //Save graphical cut as the ellipse
	    cutg = new TCutG("",n,x,y);
	    
	    //cutg = new TCutG("cutg",4);
	    //cutg->SetVarX("x");
	    //cutg->SetVarY("y");
	   
	    //cutg->SetPoint(0,ph_peak[col]-ph_width+0.001*(row-3)/2.,th_peak[row]-th_width);
	    //cutg->SetPoint(1,ph_peak[col]+ph_width+0.001*(row-3)/2,th_peak[row]-th_width);
	    //cutg->SetPoint(2,ph_peak[col]+ph_width+0.001*(row-3)/2,th_peak[row]+th_width);
	    //cutg->SetPoint(3,ph_peak[col]-ph_width+0.001*(row-3)/2,th_peak[row]+th_width);
	    //cutg->SetPoint(4,ph_peak[col]-ph_width+0.001*(row-3)/2,th_peak[row]-th_width);
	   
	      


            c5->Update();
	    
            cutg->SetName(Form("hcut_R_%d_%d_%d", FoilID, col, row));
	    
            // set axises' name
	    
            cutg->SetVarX("R.tr.tg_ph");
            cutg->SetVarY("R.tr.tg_th");
	    //cutg->SetVarX("Sieve.y");
	    //cutg->SetVarY("Sieve.x");
            cout << "done!" << endl;
        }
        else {
            cout << Form("hcut_R_%d_%d_%d", FoilID, col, row) << " is found, using old one" << endl;
        }
	
        cutg->SetLineColor(kMagenta);	
        cutg->SetLineWidth(2);	
        cutg->Draw("PL");
        c5->Update();

        cutg->Write("", TObject::kOverwrite); // Overwrite old cut

        if(UseHDpCut){
            TCanvas *c6 = new TCanvas("c6", Form("dp in col %d row %d", col+1, row+1), 800, 600);
            TH1F *h6 = new TH1F("h6", "dp", 200, dplowlimit, dpuplimit);
            assert(h6);

            c6->cd();

            TCut DpDrawCut = DrawCut + Form("hcut_R_%d_%d_%d", FoilID, col, row);
            (TCutG*)gROOT->FindObject(Form("hcut_R_%d_%d_%d", FoilID, col, row));
            T->Draw("R.tr.tg_dp>>h6", DpDrawCut);
            c6->Update();

            h6->Smooth(10);

            TSpectrum *s = new TSpectrum(5);
            Int_t nfound = s->Search(h6, 2, "", 0.10);
            c6->Update();

            Double_t *xpeaks = s->GetPositionX();
            Double_t *ypeaks = s->GetPositionY();

            cout << "Found " << nfound << " peaks in the dp plots:" << endl;
            for (int k=0; k<nfound; k++) {
                cout << k + 1 << "\t" << xpeaks[k] << "\t" << ypeaks[k] << endl;
            }
            cout << "Which peak is the C elastic peak?" << endl;

            Int_t npeak;
            cin >> npeak;

            if(npeak==0){
                cout << "No Data in this hole" << endl;
                //cutdesc << "fEvtHdr.fRun==0" << endl;
                cutdesc << "R.tr.n > 1000" << endl;
                SaveCanvas(c6, WorkDir + PlotDir + RootFileName + Form(".hdpcut_R_%d_%d_%d", FoilID, col, row), kFALSE);
                c5->cd();

                continue;
            }
            else{
                cout << "Will use peak No " << npeak << " to generate cut file" << endl;

                Int_t idx;
                h6->GetBinWithContent(ypeaks[npeak-1],idx,0,400,ypeaks[npeak-1]*0.1);

                Double_t lowthr, hithr;
                for (Int_t k=idx; k<400; k++) {
                    if (h6->GetBinContent(k)<ypeaks[npeak-1]*0.5) {
                        hithr = h6->GetBinCenter(k);
                        break;
                    }
                }
                for (Int_t k=idx; k>0; k--) {
                    if (h6->GetBinContent(k)<ypeaks[npeak-1]*0.5) {
                        lowthr = h6->GetBinCenter(k);
                        break;
                    }
                }

                cout << lowthr << "<R.tr.tg_dp<" << hithr << endl;

                TLine *l1 = new TLine(lowthr,0,lowthr,h6->GetMaximum());
                TLine *l2 = new TLine(hithr,0,hithr,h6->GetMaximum());

                l1->SetLineStyle(1);
                l1->SetLineColor(2);
                l2->SetLineStyle(1);
                l2->SetLineColor(2);
                
                l1->Draw();
                l2->Draw();

                c6->Update();
                SaveCanvas(c6, WorkDir + PlotDir + RootFileName + Form(".hdpcut_R_%d_%d_%d", FoilID, col, row), kFALSE);

                c5->cd();

                cutdesc << Form("(R.tr.tg_dp>%7.5f && R.tr.tg_dp<%7.5f) && ", lowthr, hithr);
            }
        }
        
        //cout << "Log to " << CutDescName << endl;

        cutdesc << Form("hcut_R_%d_%d_%d", FoilID, col, row) << " && ";
        if(UseVertexCut){
            cutdesc << Form("fcut_R_%d", FoilID) << " && ";
        }
        if(UseDpCut)
            cutdesc << Form("dpcut_R_%d", FoilID) << " && ";
        cutdesc << (const char*)GeneralCut << endl;
    }

    //    for(int i = rmin + nhol; i < 7; i++)
    for(int i = rmin + nhol; i < NSieveRow; i++){
      //cutdesc << "fEvtHdr.fRun==0" << endl;
      cutdesc << "R.tr.n > 1000" << endl;
      if((col%2==0 && i%2==0) || (col%2==1 && i%2==1)) hole_exist = true;
      if((col%2==0 && i%2==1) || (col%2==1 && i%2==0)) hole_exist = false;
      if(col == 13 && !(i==1 || i==15)) hole_exist = false;
      if(col == 25 && i%2==0) hole_exist = true;
      if(col == 25 && i%2==1) hole_exist = false;      
      getline(cutcsvold,line);
      cutcsv<<line<<endl;
      //cutcsv<<col<<":"<<i<<","<<hole_exist<<",0,0,"<<sieve_ph[col]<<",0,"<<sieve_th[i]<<",0,0,0"<<endl;
    }
    
    gSystem->Exec("rm -f " + WorkDir + PlotDir + RootFileName + Form(".hcut_R_%d_%d_*", FoilID, col) + ".*");
    SaveCanvas(c5, WorkDir + PlotDir + RootFileName + Form(".hcut_R_%d_%d_%d", FoilID, col, nhol), kFALSE);

    while(!cutcsvold.eof()){
      getline(cutcsvold,line);
      cutcsv<<line<<endl;
    }

       
    f1->Write();
    f1->ls();
    cutcsv.close();
    cutdesc.close();

    
    //Copy file to new template
    ofstream dest;
    dest.open(WorkDir + "apex_4647.root.cuts_full.csv");
    ifstream src;
    src.open(CutDataName);

    dest << src.rdbuf();
    
}
  //}

//#include "plot_R.C"
