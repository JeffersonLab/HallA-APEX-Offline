#include "SaveCanvas.C"
#include "TPad.h"
#include "InputAPEXL.h"
#include "Load_more_rootfiles.C"
#include "file_def.h"

#include "APEX_Sieve.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCutG.h"
#include "TLine.h"
#include "TRotation.h"
#include "TEllipse.h"
#include <iostream>
#include <fstream>
#include "InputAPEXL.h"
#include "TPaveText.h"

#include <vector>
#include <algorithm>



TCut GenrealCut = GeneralSieveCut + PID_cuts + FP_cuts;




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

  Beam_info.insert(std::pair<int,std::pair<int,int>>(1,std::pair<int,int>(1,1)));


  CutFileName = *SoureRootFile + ".FullCut.root";
  CutDescFileName = *SoureRootFile + ".VertexCut.cut";
  CutDescFileNameSieve = *SoureRootFile + ".SieveCut.cut";
  CutDescFileNameSieve_dir = *SoureRootFile + ".SieveCut_dir.cut";


  CutDescFileNameDp = *SoureRootFile + ".DpCut.cut";

  CutDescFileNameCol = *SoureRootFile + ".ColCut.cut";


  GenrealCut += get_Beamcut(Run_number);
  
}


void test_load(Int_t runnum, Int_t runnum_2 = 0){


  //~~~~~~~~ Test of Load_rootfile ~~~~~~~~
  // cout << "Start of test_load function" << endl;
  // TChain* T = Load_rootfile(runnum);

  // 		//johnw
  // cout << "After loading file !! " << endl;
  // assert(T);


  // ~~~~~~~ Test of Load_more_rootfiles ~~~~~~~~~~
  
  cout << "Start of test_load function (testing Load_more_rootfiles" << endl;

  TChain* T = Load_more_rootfiles(runnum, runnum_2);

		//johnw
  cout << "After loading file !! " << endl;
  assert(T);
  

}


void cut_Vertex(int overwrite = 0, int nfoils = 3, int FoilID = -1, int append = 0) {


  cout << "cut being used " << GenrealCut << endl << endl;

  cout << "open " << *SoureRootFile << endl;

  ReLoadcuts();
	
  //	TChain* T = Load_rootfile(Run_number);
  TChain* T = Load_more_rootfiles(Run_number, Run_number_2);


 
  gSystem->Exec("cp -vf " + CutFileName + " " + CutFileName + ".old");


  fstream cutdesc;
  
  if(append){
    cutdesc.open(CutDescFileName, ios_base::out |ios::app);
  }
  if(!append){
    cutdesc.open(CutDescFileName, ios_base::out);
  }

  //	fstream cutdesc(CutDescFileName, ios_base::out);
  assert(cutdesc.is_open());
  
  TFile *f1 = new TFile(CutFileName, "UPDATE");
  assert(f1);
  // 	f->cd(); // goback to f directory

  // Define Canvas
  TCanvas* c1 = new TCanvas("Cuts", "Cuts", 900, 900);
  gStyle->SetPalette(1, 0);



  //	TH2F* h1 = new TH2F("h1", "ReactZ vs. Target Phi", 400, -0.05, 0.035, 200,-0.4, 0.4);
  TH2F* h1 = new TH2F("h1", Form("Vertex for Foil #%d",FoilID), 400, -0.05, 0.05, 200,-0.5, 0.5);
  //	h1->CenterTitle();


  //	gStyle->SetTitleAlign(13);
  // 1st number is horizontal placement 1 = left, 

  h1->GetXaxis()->SetTitle("\\phi_{tg} [rad]");
  h1->GetYaxis()->SetTitle("Reactz [m]");
	
  // alternative to cutting on reactz vs target phi plot
  // -> use ph (FP) vs y (FP) plot instead 
	
  //	TH2F* h1 = new TH2F("h1", "ph vs y (FP)", 400, -0.05, 0.05, 400,-0.05, 0.04);
	
	
  //	TH2F* h1 = new TH2F("h1", "ReactZ vs. Target Phi", 900, -0.1, 0.06, 900,-0.5, 0.5);
  assert(h1);
  // 	TH2F* h2 = new TH2F ("h2","theta_target vs. phi_target", 900, -0.06,0.06,900,-0.07,0.07);

  // Draw ReactZ vs. Phi_rotate
  // T->Draw("L.tr.vz:L.tr.tg_ph>>h1", GenrealCut, "COLZ"); // need finer delta cut later

  // use 'reactz' as calculated in analyzer

  cout << "GenrealCut = " << GenrealCut << endl;
  T->Draw("reactz:ph_tgt>>h1", GenrealCut, "COLZ"); // need finer delta cut later



	
  // alternative to cutting on reactz vs target phi plot
  // -> use ph (FP) vs y (FP) plot instead 

  //	T->Draw("L.tr.r_ph:L.tr.r_y>>h1", GenrealCut, "COLZ"); // need finer delta cut later


  c1->Update();

  // Set y target number

  // lines added to start from foil other than zeroth
  int fmin = 0;
  if(FoilID > -1){
    fmin = FoilID;	  
  }


  // output ps files
  //c1->Print(Form("./cutfiles/gcut_L_%d.ps[",nrun));

  //	T->Draw("L.tr.vz:L.tr.tg_ph>>h1", GenrealCut, "COLZ");
  // Choose the foil you want to make cut)

  Int_t nfoil = nfoils; // Added to process seperate runs with one foil each)
  for (int i = fmin; i < nfoil + fmin; i++) {
    cout << "new foil " << i << endl;
    TCutG* cutg, *tmpcut;
    // 		if (i==fmin)  c1->Print(Form("./cutfiles/gcut_L_%d.ps(",nrun));
    // 		else c1->Print(Form("./cutfiles/gcut_L_%d.ps",nrun)); // insert page into ps file
    // 		f1->Write();
    c1->Update();
    TVirtualPad *cpad = gPad;
	  
    cout << "Get list of files: " << gROOT->GetListOfFiles() << endl;
    cout << "Testing " << Form("fcut_L_%d", i) << endl;
    // 		f1->cd();
    tmpcut = (TCutG*) gROOT->FindObject(Form("fcut_L_%d", i)); //looking for old cut definition
    // 		f->cd();
    if (tmpcut &&  !overwrite) {
      cout << Form("fcut_L_%d", i) << " is found, using old one" << endl;
      tmpcut->Draw("PL");
      // 		delete tmpcut; //delete old cut
    } else {
	    
      cout << "making cut for foil No." << i << ", waiting ..." << endl;
	    
      cutg = (TCutG*) (cpad->WaitPrimitive("CUTG", "CutG")); // making cut, store to CUTG
      c1->Update();
      cout << "done!" << endl;
      cutg->SetName(Form("fcut_L_%d", i)); //
      // set axises' nameif (i==fmin)  c1->Print(Form("./cutfiles/gcut_L_%d.ps(",nrun));
      // 			else c1->Print(Form("./cutfiles/gcut_L_%d.ps",nrun)); // insert page into ps file
      // 			f1->Write();
      cutg->SetVarX("ph_tgt");
      cutg->SetVarY("reactz");
	    


      // cut on regular analyzer replay variables
	    
      // cutg->SetVarX("L.tr.tg_ph");
      // cutg->SetVarY("L.tr.vz");

	    
      // alternative to cutting on reactz vs target phi plot
      // -> use ph (FP) vs y (FP) plot instead 
	    
      // cutg->SetVarX("L.tr.r_y");
      // cutg->SetVarY("L.tr.r_ph");

	    
      // output cut to disk
      // 			f1->cd();
      cutg->Write("", TObject::kWriteDelete); // Overwrite old cut
	    
      // output ps file
	    
	    
      // 			if (i==fmin)  c1->Print(Form("./cutfiles/gcut_L_%d.ps(",nrun));
      // 			else c1->Print(Form("./cutfiles/gcut_L_%d.ps",nrun)); // insert page into ps file
      // 			f1->Write();
    }
	  
    cout << "log to " << CutDescFileName << endl;
    cutdesc << Form("fcut_L_%d", i) << " && " << (const char *) GenrealCut
	    << endl;
  }
	
  cutdesc.close();
  f1->Write();
  cout << " --> " << CutDescFileName << endl;
  cout << " --> " << CutFileName << endl;

  SaveCanvas(c1, *SoureRootFile + ".PlotVertexCut", kFALSE);
}




// version of FP cut for single foil runs

void cut_Vertex_FP_SF(int overwrite = 0, int nfoils = 3, int FoilID = -1, int append = 0) {

  

  cout << "cut being used " << GenrealCut << endl << endl;

  cout << "open " << *SoureRootFile << endl;

  ReLoadcuts();


  TChain* T = Load_more_rootfiles(Run_number, Run_number_2);

  
  gSystem->Exec("cp -vf " + CutFileName + " " + CutFileName + ".old");


  fstream cutdesc;
  
  if(append){
    cutdesc.open(CutDescFileName, ios_base::out |ios::app);
  }
  if(!append){
    cutdesc.open(CutDescFileName, ios_base::out);
  }


  assert(cutdesc.is_open());

  TFile *f1 = new TFile(CutFileName, "UPDATE");
  assert(f1);



  TCanvas* c1 = new TCanvas("c1", "Foil cuts on ph vs y (FP)", 900, 900);
  gStyle->SetPalette(1, 0);
  


  
  // c1->Divide(2,1);		
  
  // first plot of vs Y FP with no foil cuts to show general distribution


  c1->cd(1);

  
  TH2F* h1 = new TH2F("h1", "ph vs y (FP) (no foil cut)", 400, -0.05, 0.05, 400,-0.05, 0.04);

  TH2F* h2 = new TH2F("h2", "ph vs y (FP) (with foil cut)", 400, -0.05, 0.05, 400,-0.05, 0.04);

  T->Draw("L.tr.r_ph:L.tr.r_y>>h1", GenrealCut, "COLZ"); // need finer delta cut later




  // second plot showing FP distrib with cut

  //  c1->cd(2);
  TCanvas* c2 = new TCanvas("c2", "Foil cuts on ph vs y (FP)", 900, 900);
  c2->cd(1);
  
  (TCutG*) gROOT->FindObject(Form("fcut_L_%d", FoilID));

  T->Draw("L.tr.r_ph:L.tr.r_y>>h2", GenrealCut + Form("fcut_L_%d", FoilID), "COLZ"); // need finer delta cut later

   	      

  Int_t nfoil = nfoils; // Added to process seperate runs with one foil each)

  int fmin = 0;
  if(FoilID > -1){
    fmin = FoilID;	  
  }


  for (int i = fmin; i < nfoil + fmin; i++) {

    TCutG* cutg, *tmpcut;

    //    TVirtualPad *cpad = gPad;
    //    TVirtualPad *cpad = c1->GetPad(2);
    

    tmpcut = (TCutG*) gROOT->FindObject(Form("fcut_L_FP_%d", i)); //looking for old cut definition

    
  if (tmpcut &&  !overwrite) {


    tmpcut->SetLineColor(kMagenta);	
    tmpcut->SetLineWidth(2);	

    cout << Form("fcut_L_FP_%d", i) << " is found, using old one" << endl;
    tmpcut->Draw("PL");
      // 		delete tmpcut; //delete old cut
  } else {
    
    cout << "making cut for foil No." << FoilID << ", waiting ..." << endl;
    
    //    c1->cd(2);
    //    cutg = (TCutG*) (cpad->WaitPrimitive("FPCUTG", "CutG")); // making cut, store to CUTG
    c2->cd(1);
    cutg = (TCutG*) (gPad->WaitPrimitive("CUTG", "CutG")); // making cut, store to CUTG
    //    c1->Update();

    cutg->SetLineColor(kMagenta);	
    cutg->SetLineWidth(2);	
    

    cout << "done!" << endl;
    cutg->SetName(Form("fcut_L_FP_%d", i)); //
    
    
    // alternative to cutting on reactz vs target phi plot
    // -> use ph (FP) vs y (FP) plot instead 
    
    cutg->SetVarX("L.tr.r_y");
    cutg->SetVarY("L.tr.r_ph");
    
    cutg->Write("", TObject::kWriteDelete); // Overwrite old cut
    
  }
  
  cout << "log to " << CutDescFileName << endl;
  cutdesc <<  Form("fcut_L_%d", i) << " && " <<  Form("fcut_L_FP_%d", i) << " && " << (const char *) GenrealCut << endl;
  }
  
  
  cutdesc.close();
  f1->Write();
  cout << " --> " << CutDescFileName << endl;
  cout << " --> " << CutFileName << endl;
  

}



// version of FP cut for multiple foils runs

void cut_Vertex_FP_MF(int overwrite = 0, int nfoils = 3, int FoilID = -1, int append = 0) {


  cout << "cut being used " << GenrealCut << endl << endl;

  cout << "open " << *SoureRootFile << endl;

  ReLoadcuts();


  TChain* T = Load_more_rootfiles(Run_number, Run_number_2);

  
  gSystem->Exec("cp -vf " + CutFileName + " " + CutFileName + ".old");


  fstream cutdesc;
  
  if(append){
    cutdesc.open(CutDescFileName, ios_base::out |ios::app);
  }
  if(!append){
    cutdesc.open(CutDescFileName, ios_base::out);
  }


  assert(cutdesc.is_open());

  TFile *f1 = new TFile(CutFileName, "UPDATE");
  assert(f1);
  
  
  // load foil cuts based on reactz vs phi_tgt plot

  Int_t offset = 0;
  Int_t first_foil = 0;
  Int_t last_foil = 8;



  if( Run_number == 4771){
    first_foil = 1;
    offset = 1;
    cout << "Run_number == 4771 condition" << endl;
  }
  else if( Run_number == 4774){
    first_foil = 0;
    cout << "Run_number == 4774 condition" << endl;
  }


  

  


  for(Int_t i = first_foil; i < last_foil; i++){
  
    (TCutG*) f1->GetObjectChecked(Form("fcut_L_%d", i), "TCutG"); //looking for old cut definition

  }


  //  TCanvas* c2 = new TCanvas("c2", "New Foil cut on ph vs y (FP)", 900, 900);


  TCanvas* c1 = new TCanvas("c1", "Foil cuts on ph vs y (FP)", 900, 900);
  gStyle->SetPalette(1, 0);
  


  c1->Divide(2,1);
  
  // first plot of vs Y FP with no foil cuts to show general distribution


  c1->cd(1);


  TH2F* h1 = new TH2F("h1", "ph vs y (FP)", 400, -0.05, 0.05, 400,-0.05, 0.04);

  T->Draw("L.tr.r_ph:L.tr.r_y>>h1", GenrealCut, "COLZ"); // need finer delta cut later




  // second plot showing all

  c1->cd(2);
	
  //  TH2F* h2 = new TH2F("h2", "ph vs y (FP) overlap", 400, -0.05, 0.05, 400,-0.05, 0.04);

  TH2F* h_FP_comp[7] = {NULL};


  for (Int_t i = 0; i<7; i++){

    h_FP_comp[i] = new  TH2F(Form("h_FP_comp_%d",i),"Foil cuts on ph vs y (FP)", 60, -0.05, 0.05, 60,-0.04, 0.04);
    
  }


 
  h_FP_comp[0]->SetMarkerColor(kBlue);

  cout << "cut 1 = " << TCut(Form("fcut_L_%d", 0+offset)) + TCut(Form("!fcut_L_%d", 2+offset)) << endl;

  T->Draw(Form("L.tr.r_ph:L.tr.r_y>>h_FP_comp_%d",0), GenrealCut + TCut(Form("fcut_L_%d", 0+offset)) +  TCut(Form("!fcut_L_%d", 2+offset)) + TCut(Form("!fcut_L_%d", 4+offset)) + TCut(Form("!fcut_L_%d", 6+offset)) ,  "");



  h_FP_comp[1]->SetMarkerColor(kMagenta+3);

  //  T->Draw(Form("L.tr.r_ph:L.tr.r_y>>h_FP_comp_%d",0), GenrealCut + TCut(Form("!fcut_L_%d", 0+offset)) +  TCut(Form("fcut_L_%d", 2+offset)),  "SCAT");


  //  cout << "cut 2 = " << Form("!fcut_L_%d", 0+offset) + Form("fcut_L_%d", 2+offset) << endl;

  T->Draw(Form("L.tr.r_ph:L.tr.r_y>>h_FP_comp_%d",1), GenrealCut + TCut(Form("fcut_L_%d", 0+offset)) +  TCut(Form("fcut_L_%d", 2+offset)) + TCut(Form("!fcut_L_%d", 4+offset)) + TCut(Form("!fcut_L_%d", 6+offset)),  "SCAT SAME");


  h_FP_comp[2]->SetMarkerColor(kRed);


  T->Draw(Form("L.tr.r_ph:L.tr.r_y>>h_FP_comp_%d",2), GenrealCut + TCut(Form("!fcut_L_%d", 0+offset)) +  TCut(Form("fcut_L_%d", 2+offset)) + TCut(Form("!fcut_L_%d", 4+offset)) + TCut(Form("!fcut_L_%d", 6+offset)),  "SCAT SAME");



  //  draw overlap between 2nd and 3rd foils in focal plane

  h_FP_comp[3]->SetMarkerColor(kOrange);  
  
  T->Draw(Form("L.tr.r_ph:L.tr.r_y>>h_FP_comp_%d",3), GenrealCut + TCut(Form("!fcut_L_%d", 0+offset)) +  TCut(Form("fcut_L_%d", 2+offset)) + TCut(Form("fcut_L_%d", 4+offset)) + TCut(Form("!fcut_L_%d", 6+offset)),  "SCAT SAME");


  // 3rd foil

  h_FP_comp[4]->SetMarkerColor(kYellow);  
  
  T->Draw(Form("L.tr.r_ph:L.tr.r_y>>h_FP_comp_%d",4), GenrealCut + TCut(Form("!fcut_L_%d", 0+offset)) +  TCut(Form("!fcut_L_%d", 2+offset)) + TCut(Form("fcut_L_%d", 4+offset)) + TCut(Form("!fcut_L_%d", 6+offset)),  "SCAT SAME");


  
  // 3rd and 4th foil overlap
  
  h_FP_comp[5]->SetMarkerColor(kGreen + 3);  
  
  T->Draw(Form("L.tr.r_ph:L.tr.r_y>>h_FP_comp_%d",5), GenrealCut + TCut(Form("!fcut_L_%d", 0+offset)) +  TCut(Form("!fcut_L_%d", 2+offset)) + TCut(Form("fcut_L_%d", 4+offset)) + TCut(Form("fcut_L_%d", 6+offset)),  "SCAT SAME");


  // 4th foil


  h_FP_comp[6]->SetMarkerColor(kBlue + 3);  
  
  T->Draw(Form("L.tr.r_ph:L.tr.r_y>>h_FP_comp_%d",6), GenrealCut + TCut(Form("!fcut_L_%d", 0+offset)) +  TCut(Form("!fcut_L_%d", 2+offset)) + TCut(Form("!fcut_L_%d", 4+offset)) + TCut(Form("fcut_L_%d", 6+offset)),  "SCAT SAME");




  // draw ph vs y again but with cut on particular foil (from reactz vs ph_tgt) and from this determine a second cut for that foil
  

  TCanvas* c2 = new TCanvas("c2", "New Foil cut on ph vs y (FP)", 900, 900);



  Int_t nfoil = nfoils; // Added to process seperate runs with one foil each)


  c2->cd(1);

  TH2F* h2 = new TH2F("h2", "ph vs y (FP)", 400, -0.05, 0.05, 400,-0.05, 0.04);

  
  
  
  int fmin = 0;
  if(FoilID > -1){
    fmin = FoilID; 
  }


  c2->cd(1);
  

  

  for (int i = fmin; i < nfoil + fmin; i++) {
    cout << "new foil " << i << endl;
    c2->cd(1);

    c2->Update();

    //    T->Draw("L.tr.r_ph:L.tr.r_y>>h2", GenrealCut + TCut(Form("fcut_L_%d",i)), "COLZ"); // need finer delta cut later


    T->Draw("L.tr.r_ph:L.tr.r_y:th_tgt>>h2", GenrealCut + TCut(Form("fcut_L_%d",i)), "COLZ"); 

    

    // following lines are finding TCutG for other foils not being cut
            




    
    
    
    TCutG* other_foils[3];

    
    //    for(Int_t FoilID = first_foil; FoilID<last_foil; FoilID=FoilID+2){
    
    Int_t count = 0;
    
    for(Int_t j = 0; j<4; j++){

      if((2*j + offset) != i){
	
	cout << " Begin foil " << 2*j + offset << endl;
	
	// cout << "count == " << count << endl;
	// cout << "Foil Printed: " << 2*j + offset << endl;

	// test if other foils have FP cut
	
	TCutG* tmpcut =  (TCutG*) gROOT->FindObject(Form("fcut_L_FP_%d", 2*j+offset));

	if (tmpcut){
	  other_foils[count] =  (TCutG*) gROOT->FindObject(Form("fcut_L_FP_%d", 2*j+offset));

	  other_foils[count]->SetLineColor(kRed);
	// cout << "Foil Printed: " << 2*j + offset << endl;
	  other_foils[count]->Draw("PL same");
	// cout << "Foil Printed: " << 2*j + offset << endl;

	}
	count++;

	cout << " end foil " << 2*j + offset << endl;
      }
      else{

	cout << "Foil Not Printed: " << 2*j + offset << endl;
      }
 
    }


    cout << "Reached this stage" << endl;

    (TCutG*) gROOT->FindObject(Form("fcut_L_FP_%d", i));
    



    TCutG* cutg, *tmpcut;
    // 		if (i==fmin)  c1->Print(Form("./cutfiles/gcut_L_%d.ps(",nrun));
    // 		else c1->Print(Form("./cutfiles/gcut_L_%d.ps",nrun)); // insert page into ps file
		// 		f1->Write();
    c2->Update();
    //    TVirtualPad *cpad = gPad;
    TVirtualPad *cpad = c2->cd(1);;
		
    cout << "Get list of files: " << gROOT->GetListOfFiles() << endl;
    cout << "Testing " << Form("fcut_L_FP_%d", i) << endl;
    // 		f1->cd();
    tmpcut = (TCutG*) gROOT->FindObject(Form("fcut_L_FP_%d", i)); //looking for old cut definition
    // 		f->cd();
    if (tmpcut &&  !overwrite) {
      cout << Form("fcut_L_FP_%d", i) << " is found, using old one" << endl;
      tmpcut->Draw("PL");
      // 		delete tmpcut; //delete old cut
    } else {


      c2->cd(1);
      cout << "making cut for foil No." << i << ", waiting ..." << endl;

      cutg = (TCutG*) (cpad->WaitPrimitive("CUTG", "CutG")); // making cut, store to CUTG
      c2->Update();
      cout << "done!" << endl;
      cutg->SetName(Form("fcut_L_FP_%d", i)); //


      // alternative to cutting on reactz vs target phi plot
      // -> use ph (FP) vs y (FP) plot instead 

      cutg->SetVarX("L.tr.r_y");
      cutg->SetVarY("L.tr.r_ph");
    
      cutg->Write("", TObject::kWriteDelete); // Overwrite old cut

    }
	
    
    cout << "log to " << CutDescFileName << endl;
    cutdesc <<  Form("fcut_L_%d", i) << " && " <<  Form("fcut_L_FP_%d", i) << " && " << (const char *) GenrealCut	    << endl;
  }
  
  cutdesc.close();
  f1->Write();
  cout << " --> " << CutDescFileName << endl;
  cout << " --> " << CutFileName << endl;
  
	//SaveCanvas(c1, *SoureRootFile + ".PlotVertexCut", kFALSE);
	// }


}





void cut_Vertex_FP(int overwrite = 0, int nfoils = 3, int FoilID = -1, int append = 0) {


  Bool_t multifoil = IsMultiFoil(Run_number);

  if(multifoil){

    cut_Vertex_FP_MF(overwrite, nfoils, FoilID, append);
      
  }
  else if(!multifoil){

    cut_Vertex_FP_SF(overwrite, nfoils, FoilID, append);
      }

}




void display_Foils(){

  // function displays foil tg-phi z plots

  TChain* T = Load_more_rootfiles(Run_number, Run_number_2);

 ReLoadcuts();

 TCanvas *c1 = new TCanvas("c1","PlotFoils",1000,1000);
 
 TH2F* h1 = new TH2F("h1",Form("Vertex plot for run #%d",Run_number), 400, -0.05, 0.05, 200,-0.4, 0.4);
 //	h1->CenterTitle();

 h1->GetXaxis()->SetTitle("\\phi_{tg} [rad]");
 h1->GetYaxis()->SetTitle("Reactz [m]");
	
 T->Draw("reactz:ph_tgt>>h1", GenrealCut, "COLZ"); // need finer delta cut later

 TCanvas *c2 = new TCanvas("c2","PlotFoils FP",1000,1000);

 TH2F* h2 = new TH2F("h2",Form("FP Vertex plot for run #%d",Run_number), 400, -0.05, 0.05, 200,-0.05, 0.04);


 h2->GetXaxis()->SetTitle("y_{FP} [m]");
 h2->GetYaxis()->SetTitle("#phi_{FP} [rad]");
	
 T->Draw("L.tr.r_ph:L.tr.r_y>>h2", GenrealCut, "COLZ"); 
 
}


void display_Sieve(int FoilID = 0, int col = 0){

  // function plots sieve x-y and theta-phi for a chosen foil

 TChain* T = Load_more_rootfiles(Run_number, Run_number_2);

 ReLoadcuts();
 
 cout << "RootFileName = " << RootFileName << endl;
 cout << "Run_number = " << Run_number << endl;
 cout << "General cut = " << GenrealCut << endl << endl;


 TFile *f1 = new TFile(CutFileName, "READ");
 assert(f1);
 (TCutG*) f1->GetObjectChecked(Form("fcut_L_%d", FoilID), "TCutG"); //looking for old cut definition
 
 // line added to extract previous FP vertex cut if existing
 (TCutG*) f1->GetObjectChecked(Form("fcut_L_FP_%d", FoilID), "TCutG"); //looking for old cut definition

 // default vertex cuts (designed to pass all cuts)

 // TCutG *fcut = new TCutG(Form("fcut_L_%d", FoilID),4);

 // TCutG *fcut_FP = new TCutG(Form("fcut_L_FP_%d", FoilID),4);

 // fcut->SetVarX("L.tr.r_x");
 // fcut->SetVarY("L.tr.r_y");

 // fcut->SetPoint(0,-1e10,-1e10);
 // fcut->SetPoint(1,-1e10,+1e10);
 // fcut->SetPoint(2,+1e10,+1e10);
 // fcut->SetPoint(3,+1e10,-1e10);

 
 // fcut_FP->SetVarX("L.tr.r_x");
 // fcut_FP->SetVarY("L.tr.r_y");

 // fcut_FP->SetPoint(0,-1e10,-1e10);
 // fcut_FP->SetPoint(1,-1e10,+1e10);
 // fcut_FP->SetPoint(2,+1e10,+1e10);
 // fcut_FP->SetPoint(3,+1e10,-1e10);

 
 

 TCanvas *c1 = new TCanvas("c1","PlotSieve",1000,1000);
 
 TH2F* h2 = new TH2F("h2", Form("Sieve plot for Foil #%d",FoilID), 400,-30,30,400,-65,65);
 
 h2->GetXaxis()->SetTitle("#phi_{tg} [mrad]");
 h2->GetYaxis()->SetTitle("#theta_{tg} [mrad]");
 //  TH2F* h2 = new TH2F("h2", "theta_target vs. phi_target", 300, -0.06, 0.04, 300, -0.08, 0.08);
 // h2->SetMinimum(2);
 

 c1->cd(0);
 
 cout << "Pre Tree-Draw" << endl;
 
 
 // T->Draw("L.tr.tg_th:L.tr.tg_ph>>h2", GenrealCut + TCut(Form("fcut_L_%d", FoilID)), "COLZ");
  



 T->Draw("(1000*th_tgt):(1000*ph_tgt)>>h2", GenrealCut + TCut(Form("fcut_L_%d", FoilID)) + TCut(Form("fcut_L_FP_%d", FoilID)), "COLZ");



  // test of how many events pass each cut
  cout << "How many pass overall events: " << T->GetEntries() << endl;
  cout << "How many pass general cut: " << T->GetEntries(GenrealCut) << endl;
  cout << "How many pass foil cut: " << T->GetEntries(TCut(Form("fcut_L_%d", FoilID))) << endl;
  cout << "How many pass foil FP cut: " << T->GetEntries(TCut(Form("fcut_L_FP_%d", FoilID))) << endl;
  cout << "How many pass all cuts: " << T->GetEntries(GenrealCut + TCut(Form("fcut_L_FP_%d", FoilID)) + TCut(Form("fcut_L_%d", FoilID))) << endl;
  
  TRotation fTCSInHCS;
  TVector3 TCSX(0,-1,0);
  TVector3 TCSZ(TMath::Sin(HRSAngle),0,TMath::Cos(HRSAngle));
  TVector3 TCSY = TCSZ.Cross(TCSX);
  fTCSInHCS.RotateAxes(TCSX,TCSY,TCSZ);
  
  fPointingOffset.SetXYZ(-MissPointZ*TMath::Sin(HRSAngle)*TMath::Cos(HRSAngle),(Double_t)MissPointY,MissPointZ*TMath::Sin(HRSAngle)*TMath::Sin(HRSAngle));



  const TVector3 BeamSpotHCS_average(BeamX_average[FoilID] + (targetfoils[FoilID]/BeamZDir_average)*BeamXDir_average[FoilID], BeamY_average[FoilID] + (targetfoils[FoilID]/BeamZDir_average)*BeamYDir_average[FoilID], targetfoils[FoilID]);

  TVector3 BeamSpotTCS_average = fTCSInHCS.Inverse()*(BeamSpotHCS_average-fPointingOffset);
  

  cout << "BeamX_average + (targetfoils[FoilID]/BeamZDir_average)*BeamXDir_average = " << BeamX_average[FoilID] + (targetfoils[FoilID]/BeamZDir_average)*BeamXDir_average[FoilID] << endl;
  cout << "BeamX_average = " << BeamX_average[FoilID] << ", targetfoils[FoilID] = " << targetfoils[FoilID] << ", BeamZDir_average = " << BeamZDir_average << ", BeamXDir_average = " << BeamXDir_average[FoilID] << ", BeamY_average = " << BeamY_average[FoilID] << ", BeamYDir_average = " << 
    BeamYDir_average << endl;
  cout << "BeamSpotHCS_average = [" << BeamSpotHCS_average.X() << ", " << BeamSpotHCS_average.Y() << ", " << BeamSpotHCS_average.Z() << "]" << endl;
  cout << "BeamSpotTCS_average = [" << BeamSpotTCS_average.X() << ", " << BeamSpotTCS_average.Y() << ", " << BeamSpotTCS_average.Z() << "]" << endl;
  
  //  const Double_t plotwidth = 0.0015;
  Double_t plotwidth = 0.0030*1000;
  Double_t plotheight = 0.0080*1000;
  
  
  for(UInt_t Hole = 0; Hole < NHoles; Hole++){


    Color_t color = kRed;
    Int_t width = 1;
    
    if (Hole == 112){
      color = kBlack;
      width = 3;
    }
    if( Hole == 160){
      color = kBlack;
      width = 3;
    }
    
      
    //    TVector3 Hole_pos = GetSieveHoleTCS(Hole);

    TVector3 Hole_pos = GetSieveHoleCorrectionTCS(FoilID,Hole);
      
      
    TVector3 MomDirectionTCS_hole = Hole_pos - BeamSpotTCS_average;

      
    Double_t theta_hole = MomDirectionTCS_hole.X()/MomDirectionTCS_hole.Z();
    Double_t phi_hole = MomDirectionTCS_hole.Y()/MomDirectionTCS_hole.Z();

    theta_hole *= 1000;
    phi_hole *=  1000;
    
    // + type crosses
    
    TLine *lh = new TLine(phi_hole-plotwidth,theta_hole,phi_hole+plotwidth,theta_hole);
    TLine *lv = new TLine(phi_hole,theta_hole-plotheight,phi_hole,theta_hole+plotheight);


    lh->SetLineColor(color);
    lv-> SetLineColor(color);


    lh->SetLineWidth(width);
    lv->SetLineWidth(width);
    
    lh -> Draw("same");
    lv -> Draw("same");

    // adding text  lable to columns and rows

    
    std::vector<int> x_y;
  
    x_y = Get_Col_Row(Hole);    

    Int_t col =  x_y[0];
    Int_t row =  x_y[1];
    
    
    if( col == 0 || col ==1){
      TPaveText *row_text = new TPaveText(phi_hole-10,theta_hole-2,phi_hole-5,theta_hole+2);
      row_text->AddText(Form("row_%d",row));
      row_text->Draw("same");
    }

    if( row == 0 || row ==1){
      TPaveText *col_text = new TPaveText(phi_hole-2,theta_hole+2,phi_hole+2,theta_hole-5);
      col_text->AddText(Form("col_%d",col));
      col_text->Draw("same");
    }

    
    
    
    
    // Saltire-type crosses
    // TLine *lc1 = new TLine(phi_hole-plotwidth,theta_hole-plotwidth,phi_hole+plotwidth,theta_hole+plotwidth);
    // lc1->SetLineWidth(width);
    // TLine *lc2 = new TLine(phi_hole-plotwidth,theta_hole+plotwidth,phi_hole+plotwidth,theta_hole-plotwidth);
    // lc2->SetLineWidth(width);
    


    // lc1->SetLineColor(color);
    // lc2-> SetLineColor(color);
    
    // lc1 -> Draw("same");
    // lc2 -> Draw("same");
    
  }


  
  TCanvas *c2 = new TCanvas("c2","PlotSieve",1000,1000);


  
  TH2F* h3 = new TH2F("h3", Form("Sieve plot for Foil #%d",FoilID), 300, -0.04, 0.04, 300, -0.08, 0.08);
 

  h3->GetXaxis()->SetTitle("y_sieve_alt [m]");
  h3->GetYaxis()->SetTitle("x_sieve [m]");

  // h3->GetXaxis()->SetTitle("y_tg [m]");
  // h3->GetYaxis()->SetTitle("x_tg [m]");
 //  TH2F* h2 = new TH2F("h2", "theta_target vs. phi_target", 300, -0.06, 0.04, 300, -0.08, 0.08);
  h3->SetMinimum(2);



  T->Draw("x_sieve:y_sieve_alt>>h3", GenrealCut + TCut(Form("fcut_L_FP_%d", FoilID)) + TCut(Form("fcut_L_%d", FoilID)), "COLZ");



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
      

    plotwidth = 0.005;
    plotheight = 0.005;
    
    // + type crosses
    
    TLine *lh = new TLine(posx-plotwidth,posy,posx+plotwidth,posy);
    TLine *lv = new TLine(posx,posy-plotheight,posx,posy+plotheight);


    lh->SetLineColor(color);
    lv-> SetLineColor(color);


    lh->SetLineWidth(width);
    lv->SetLineWidth(width);
    
    lh -> Draw("same");
    lv -> Draw("same");



    
    // Saltire-type crosses
    // TLine *lc1 = new TLine(posx-plotwidth,posy-plotwidth,posx+plotwidth,posy+plotwidth);
    // lc1->SetLineWidth(width);
    // TLine *lc2 = new TLine(posx-plotwidth,posy+plotwidth,posx+plotwidth,posy-plotwidth);
    // lc2->SetLineWidth(width);
	
 	
    // lc1->SetLineColor(color);
    // lc2->SetLineColor(color);
      
    // lc1->Draw("same");
    // lc2->Draw("same");
      
  }
   


  // extra comparison canvas

  TCanvas *c2_b = new TCanvas("c2_b","PlotSieve",1000,1000);

  TH2F* h3_b = new TH2F("h3_b", Form("Sieve plot for Foil #%d",FoilID), 300, -0.04, 0.04, 300, -0.08, 0.08);
 

  h3_b->GetXaxis()->SetTitle("y_sieve [m]");
  h3_b->GetYaxis()->SetTitle("x_sieve [m]");

  h3_b->SetMinimum(2);

  T->Draw("x_sieve:y_sieve>>h3_b", GenrealCut + TCut(Form("fcut_L_FP_%d", FoilID)) + TCut(Form("fcut_L_%d", FoilID)), "COLZ");

  
  // third plot of y_tg against Zpos*tan (check for correlation)

  TCanvas *c3 = new TCanvas("c3","PlotSieve - Corr",1000,1000);
  
  TH2F* h4 = new TH2F("h4", Form("Corr Sieve plot for Foil #%d",FoilID), 300, -0.04, 0.04, 300, -0.08, 0.08);

  h4->GetXaxis()->SetTitle("Z * tan(\\phi_tg) ");
  h4->GetYaxis()->SetTitle("y_tg []");

 //  TH2F* h2 = new TH2F("h2", "theta_target vs. phi_target", 300, -0.06, 0.04, 300, -0.08, 0.08);
  h4->SetMinimum(2);

  const Double_t ZPos = 31.23 * 25.4e-3;

  TH2F *h2_cp = (TH2F*)h2->Clone("h2_cp");
  


  T->Draw("1000*th_tgt:1000*ph_tgt>>h2_cp", GenrealCut + TCut(Form("fcut_L_%d", FoilID)) + TCut(Form("fcut_L_FP_%d", FoilID)), "COLZ");




  // fourth plot of y-beam vs x-beam

  TCanvas *c4 = new TCanvas("c4","x sieve vs x beam",1000,1000);
  
  TH2F* h5 = new TH2F("h5", Form("x sieve vs x beam (foil #%d)",FoilID), 300, 0.00, 0.01, 300, -0.08, 0.08);

  h5->GetXaxis()->SetTitle("x_beam [m] ");
  h5->GetYaxis()->SetTitle("x_sieve [m]");


 //  TH2F* h2 = new TH2F("h2", "theta_target vs. phi_target", 300, -0.06, 0.04, 300, -0.08, 0.08);
  h5->SetMinimum(2);

  

  T->Draw("1000*th_tgt:1000*ph_tgt>>h2", GenrealCut + TCut(Form("fcut_L_FP_%d", FoilID)) + TCut(Form("fcut_L_%d", FoilID)), "COLZ");



  // fifth plot of x-sieve vs y_tgt_alt

  TCanvas *c5 = new TCanvas("c5","c5",1000,1000);
  
  TH2F* h6 = new TH2F("h6", Form("#theta_{tg} vs y_{FP}(foil #%d)",FoilID), 300, -0.04, 0.04, 300, -65,65);

  h6->GetXaxis()->SetTitle("y_{FP} [m] ");
  h6->GetYaxis()->SetTitle("#theta_{tg} [mrad]");


 //  TH2F* h2 = new TH2F("h2", "theta_target vs. phi_target", 300, -0.06, 0.04, 300, -0.08, 0.08);
  h6->SetMinimum(2);

  

  T->Draw("(1000*th_tgt):(L.tr.r_y)>>h6", GenrealCut + TCut(Form("fcut_L_FP_%d", FoilID)) + TCut(Form("fcut_L_%d", FoilID)), "COLZ");



  // sixth plot of theta_FP vs y_FP 

  TCanvas *c6 = new TCanvas("c6","c6",1000,1000);

  c6->Divide(2);
  
  TH2F* h7 = new TH2F("h7", Form("FP theta vs y (foil #%d)",FoilID), 300, -0.04, 0.04, 300, -0.08, 0.08);

  h7->GetXaxis()->SetTitle("y FP [m] ");
  h7->GetYaxis()->SetTitle("theta FP [rad]");


 //  TH2F* h2 = new TH2F("h2", "theta_target vs. phi_target", 300, -0.06, 0.04, 300, -0.08, 0.08);
  h7->SetMinimum(2);


  TH2F* h7_b = new TH2F("h7_b", "theta vs y FP", 300, -0.06, 0.04, 300, -0.08, 0.08);
  


  c6->cd(1);
  //  T->Draw("(L.tr.r_th):(L.tr.r_y)>>h7_b", GenrealCut + TCut(Form("fcut_L_FP_%d", FoilID)) + TCut(Form("fcut_L_%d", FoilID)), "COLZ");
  T->Draw("1000*th_tgt:1000*ph_tgt>>h2", GenrealCut + TCut(Form("fcut_L_FP_%d", FoilID)) + TCut(Form("fcut_L_%d", FoilID)), "COLZ");
  c6->cd(2);
  T->Draw("(1000*th_tgt):(L.tr.r_y)>>h6", GenrealCut + TCut(Form("fcut_L_FP_%d", FoilID)) + TCut(Form("fcut_L_%d", FoilID)), "COLZ");


  
  c1->Update();


  gSystem->Exec(Form("mkdir plots/%d",Run_number));

  c1->Print(Form("plots/%d/sieve_angles.pdf",Run_number));

  c2->Print(Form("plots/%d/sieve_XY.pdf",Run_number));

}



// function designed for optics targets to display all sieve for all foils
// creates similar plots to display_Sieve:
//  theta_tg vs phi_tg and theta_FP vs y_FP
void display_all_Sieve(TString plot_name = ""){

  // check if runs are optics target runs
  if(!IsMultiFoil(Run_number) || !IsMultiFoil(Run_number_2)){
    cout << "display_all_Sieve() must be used with Optics target runs" << endl;
    cout << "Either " << Run_number << ", or " << Run_number_2 << " are not optics target runs" <<  endl;
    std::exit(0);       
  }

  Int_t foil_nos[4];

  Int_t foil_no = GetFoilID(Run_number);

  for(UInt_t i = 0; i<(sizeof(foil_nos)/sizeof(foil_nos[0])); i++ ){
    foil_nos[i] = foil_no + i*2;
  }

   
  TChain* T = Load_more_rootfiles(Run_number, Run_number_2);

  ReLoadcuts();
  cout << "General cut = " << GenrealCut << endl << endl;
  
  TFile *f1 = new TFile(CutFileName, "READ");
  assert(f1);

  TCanvas* c_tg_foils[4];
  TCanvas* c_FP_foils[4];
  

  TH2F* h2[4];
  
  TH2F* thfp_v_yfp[4];
  

  TCutG* foil_tg_gcuts[4];
  TCutG* foil_FP_gcuts[4];



  // TCutG *cutg = new TCutG("mycut",5);
  // cutg->SetVarX("L.tr.tg_ph");
  // cutg->SetVarY("L.tr.tg_y");
  // cutg->SetPoint(0, -0.028, 0.025);
  // cutg->SetPoint(1, 0.010, 0.025);
  // cutg->SetPoint(2, 0.030, -0.030);
  // cutg->SetPoint(3, 0.009, -0.030);
  // cutg->SetPoint(4, -0.028, 0.025);
  
  for(UInt_t i = 0; i<(sizeof(foil_nos)/sizeof(foil_nos[0])); i++ ){
    Int_t FoilID = foil_nos[i];
    
    foil_tg_gcuts[i] = (TCutG*) f1->GetObjectChecked(Form("fcut_L_%d", FoilID), "TCutG"); //looking for foil cut definition
    foil_FP_gcuts[i] = (TCutG*) f1->GetObjectChecked(Form("fcut_L_FP_%d", FoilID), "TCutG"); //looking for foil cut definition


    c_tg_foils[i] = new TCanvas(Form("foil Tg  %d", FoilID),Form("foil Tg  %d", FoilID),1000,1000);

    h2[i] = new TH2F(Form("h2_%d",FoilID), Form("theta_target vs. phi_target, Foil %d", FoilID),  400,-30,30,400,-65,65);

    h2[i]->GetXaxis()->SetTitle("#phi_{tg} [mrad]");
    h2[i]->GetYaxis()->SetTitle("#theta_{tg} [mrad]");
    //    h2[i]->SetMinimum(2);
    
    T->Draw(Form("(1000*th_tgt):(1000*ph_tgt)>>h2_%d",FoilID), GenrealCut + TCut(Form("fcut_L_%d", FoilID)) + TCut(Form("fcut_L_FP_%d", FoilID)) + "mycut", "COLZ");
       
    c_tg_foils[i]->Print(Form("plots/Optics/tg_sieve_Foil_%d",FoilID) + plot_name + ".pdf");
    
    // test of how many events pass each cut
    // cout << "How many pass overall events: " << T->GetEntries() << endl;
    // cout << "How many pass general cut: " << T->GetEntries(GenrealCut) << endl;
    // cout << "How many pass foil cut: " << T->GetEntries(TCut(Form("fcut_L_%d", FoilID))) << endl;
    // cout << "How many pass foil FP cut: " << T->GetEntries(TCut(Form("fcut_L_FP_%d", FoilID))) << endl;
    // cout << "How many pass all cuts: " << T->GetEntries(GenrealCut + TCut(Form("fcut_L_FP_%d", FoilID)) + TCut(Form("fcut_L_%d", FoilID))) << endl << endl;
    
    c_FP_foils[i] = new TCanvas(Form("foil FP  %d", FoilID),Form("foil FP  %d", FoilID),1000,1000);

    thfp_v_yfp[i] = new TH2F(Form("thfp_v_yfp_%d",FoilID),Form("#theta_{FP} vs y_{FP}, foil %d",FoilID),300,-0.05,0.05,300,-35,30);

    thfp_v_yfp[i]->GetYaxis()->SetTitleOffset(1.0);
    thfp_v_yfp[i]->GetXaxis()->SetTitleSize(0.05);
    thfp_v_yfp[i]->GetYaxis()->SetTitleSize(0.05);
    thfp_v_yfp[i]->GetXaxis()->SetTitle("y (FP) [m]");
    thfp_v_yfp[i]->GetYaxis()->SetTitle("th (FP) [mrad]");
    
    T->Draw(Form("L.tr.r_th*1000:L.tr.r_y>>thfp_v_yfp_%d",FoilID),GenrealCut + TCut(Form("fcut_L_%d", FoilID)) + TCut(Form("fcut_L_FP_%d", FoilID)) + "mycut","colz");

    c_FP_foils[i]->Print(Form("plots/Optics/FP_sieve_Foil_%d",FoilID) + plot_name + ".pdf");

  }

  
  // draw target and focal plane foil cuts used on seperate canvases

  
  TH2F* h_f_tg = new TH2F("h_f_tg", "Vertex Tg", 400, -0.05, 0.05, 200,-0.5, 0.5);
  h_f_tg->GetXaxis()->SetTitle("\\phi_{tg} [rad]");
  h_f_tg->GetYaxis()->SetTitle("Reactz [m]");
  
  TCanvas* c_tg_all_foils = new TCanvas("foil Tg cuts","Tg Sieve",1000,1000);
 
  T->Draw("reactz:ph_tgt>>h_f_tg", GenrealCut + "mycut", "COLZ");

  for(UInt_t i = 0; i<(sizeof(foil_nos)/sizeof(foil_nos[0])); i++ ){
    
    foil_tg_gcuts[i]->SetLineColor(kRed);
    foil_tg_gcuts[i]->Draw("PL same");
  }

  
  TH2F* h_f_FP = new TH2F("h_f_FP", "#phi vs y (FP)", 400, -0.05, 0.05, 400,-0.05, 0.04);
  h_f_FP->GetXaxis()->SetTitle("y_{FP} [m]");
  h_f_FP->GetYaxis()->SetTitle("#phi [rad]");

  c_tg_all_foils->Print("plots/Optics/tg_sieve_all" + plot_name + ".pdf");
  
  
  
  TCanvas* c_fp_all_foils = new TCanvas("foil FP cuts","FP Sieve",1000,1000);

  T->Draw("L.tr.r_ph:L.tr.r_y>>h1", GenrealCut + "mycut", "COLZ");

  for(UInt_t i = 0; i<(sizeof(foil_nos)/sizeof(foil_nos[0])); i++ ){    
    foil_FP_gcuts[i]->SetLineColor(kRed);
    foil_FP_gcuts[i]->Draw("PL same");
  }

  c_fp_all_foils->Print("plots/Optics/fp_sieve_all" + plot_name + ".pdf");
  
}



void CutSieve_ellipse(int FoilID = 0, uint col = 0, int overwrite = 0, int append = 1){

  TChain* T = Load_more_rootfiles(Run_number, Run_number_2);

  ReLoadcuts();

  cout << "RootFileName = " << RootFileName << endl;
  cout << "Run_number = " << Run_number << endl;
  cout << "General cut = " << GenrealCut << endl << endl;

  // creating and filling expectation value of holes (theta and phi)
  Double_t sieve_ph[NSieveCol], sieve_th[NSieveRow];

  for(UInt_t r_count = 0; r_count < NSieveRow; r_count++){
    sieve_th[r_count] = 0;
  }

  for(UInt_t c_count = 0; c_count < NSieveCol; c_count++){
    sieve_ph[c_count] = 0;
  }




  // file that stores TCutGs (from foil and sieve)
  gSystem->Exec("cp -vf " + CutFileName + " " + CutFileName + ".old");
  TFile *f1 = new TFile(CutFileName, "UPDATE");
  assert(f1);
   (TCutG*) f1->GetObjectChecked(Form("fcut_L_%d", FoilID), "TCutG"); //looking for old cut definition

  // line added to extract previous FP vertex cut if existing
   (TCutG*) f1->GetObjectChecked(Form("fcut_L_FP_%d", FoilID), "TCutG"); //looking for old cut definition


  // Text file containing text names of all cuts (lines have TCutG names combined with general cuts used)
  fstream cutdesc;
    
  gSystem->Exec("cp -vf " + CutDescFileNameSieve_ell + " " +  CutDescFileNameSieve_ell + ".old");
  if(append){
    cutdesc.open(CutDescFileNameSieve_ell, ios_base::out |ios::app);
  }
  if(!append){
    cutdesc.open(CutDescFileNameSieve_ell, ios_base::out);
  }
  assert(cutdesc.is_open());
  

  // could add 'PlotDir' so save picture of cuts



  // set-up old and new csv files (to store names of ellipse cuts (and properties ie r1,r2,x and y of each ellipse)

  TString csv_name = Sieve_CSV_name + RootFileName + ".csv";


  gSystem->Exec("cp -vf " + csv_name + " " + csv_name + ".old");
  ifstream cutcsvold;
  cutcsvold.open(csv_name + ".old");  
    
  ofstream cutcsv;
  cutcsv.open(csv_name);
  
  
  cutcsv<<fixed<<setprecision(1);
  cout<<fixed<<setprecision(10);

  TDatime* date = new TDatime();
  cutcsv<<date->GetDay()<<"/"<<date->GetMonth()<<"/"<<date->GetYear()<<" (dd/mm/yyyy)"<<endl;
  //  cutcsv << "Hole ID (col:row),Hole Exists,Included in opt,Ellipse ph cen,Expected ph,Ellipse th cen,Expected th,Semi axis ph,Semi axis th,Ellipse Tilt (deg),Ellipse ph rms,Ellipse th rms,Statistics,All angles in mrad except ellipse tilt which is positive counterclockwise from vertical axis"<<endl;
  cutcsv << "Hole ID: col, row, Included in opt, Ellipse ph cen, Expected ph, Ellipse th cen, Expected th, Semi axis ph, Semi axis th, Ellipse Tilt (deg), -  All angles in mrad except ellipse tilt which is positive counterclockwise from vertical axis"<<endl;
  
  TCanvas *c1 = new TCanvas("c1","PlotSieve",1000,1000);
  
  TH2F* h2 = new TH2F("h2", "theta_target vs. phi_target", 300, -0.04, 0.04, 300, -0.08, 0.08);


  h2->GetXaxis()->SetTitle("\\phi_{tg}");
  h2->GetYaxis()->SetTitle("\\theta_{tg}");
   //  TH2F* h2 = new TH2F("h2", "theta_target vs. phi_target", 300, -0.06, 0.04, 300, -0.08, 0.08);
  h2->SetMinimum(2);
  
   

  
  
  c1->cd(0);

  cout << "Pre Tree-Draw" << endl;
  
  
  T->Draw("th_tgt:ph_tgt>>h2", GenrealCut + Form("fcut_L_%d", FoilID) + TCut(Form("fcut_L_FP_%d", FoilID)) , "COLZ");


  
  cout << "Post Tree-Draw" << endl;
  
  c1->Update();
  
  
  
  
  // also draw red crosses where holes 'should' be




  TRotation fTCSInHCS;
  TVector3 TCSX(0,-1,0);
  TVector3 TCSZ(TMath::Sin(HRSAngle),0,TMath::Cos(HRSAngle));
  TVector3 TCSY = TCSZ.Cross(TCSX);
  fTCSInHCS.RotateAxes(TCSX,TCSY,TCSZ);
  
  fPointingOffset.SetXYZ(-MissPointZ*TMath::Sin(HRSAngle)*TMath::Cos(HRSAngle),(Double_t)MissPointY,MissPointZ*TMath::Sin(HRSAngle)*TMath::Sin(HRSAngle));



  const TVector3 BeamSpotHCS_average(BeamX_average[FoilID] + (targetfoils[FoilID]/BeamZDir_average)*BeamXDir_average[FoilID], BeamY_average[FoilID] + (targetfoils[FoilID]/BeamZDir_average)*BeamYDir_average[FoilID], targetfoils[FoilID]);

  TVector3 BeamSpotTCS_average = fTCSInHCS.Inverse()*(BeamSpotHCS_average-fPointingOffset);
  

  cout << "BeamX_average + (targetfoils[FoilID]/BeamZDir_average)*BeamXDir_average = " << BeamX_average[FoilID] + (targetfoils[FoilID]/BeamZDir_average)*BeamXDir_average[FoilID] << endl;
  cout << "BeamX_average = " << BeamX_average[FoilID] << ", targetfoils[FoilID] = " << targetfoils[FoilID] << ", BeamZDir_average = " << BeamZDir_average << ", BeamXDir_average = " << BeamXDir_average[FoilID] << ", BeamY_average = " << BeamY_average[FoilID] << ", BeamYDir_average = " << 
    BeamYDir_average[FoilID] << endl;
  cout << "BeamSpotHCS_average = [" << BeamSpotHCS_average.X() << ", " << BeamSpotHCS_average.Y() << ", " << BeamSpotHCS_average.Z() << "]" << endl;
  cout << "BeamSpotTCS_average = [" << BeamSpotTCS_average.X() << ", " << BeamSpotTCS_average.Y() << ", " << BeamSpotTCS_average.Z() << "]" << endl;
  
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
    
      
    //    TVector3 Hole_pos = GetSieveHoleTCS(Hole);

    TVector3 Hole_pos = GetSieveHoleCorrectionTCS(FoilID,Hole);
      
      
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
  

  // hole cutting stage

  int nhol = 0;
  cout << "How many holes in this No." << col << " column?" << endl;
  cin >> nhol;
  if (nhol < 0){
    exit(EXIT_FAILURE);
  }
    
  cout << "min hole in this column? : ";
  int rmin = -1;
  cin >> rmin;
  if (rmin < 0){
    exit(EXIT_FAILURE);
  }



  // parameter to keep track of 'current hole number': next hole entry to be added to csv file
  // this starts with hole = 0, then 1 and so on
  Int_t cur_hole = 0;

  string line; // used to gather csv file
  getline(cutcsvold,line); // these lines are for first two line of csv file (explenatory)
  getline(cutcsvold,line);

  //make cuts for each holes 
  for (int cnt = 0; cnt < nhol; cnt++) {


   c1->cd(0);

    int k = rmin + cnt * 2; // line as holes in column are in every other row
    
    if (k < 0)
      continue;

    Int_t hole_no;
    vector<int> col_row;
    
    for(UInt_t ii = 0; ii<col; ii++){
      for(UInt_t jj = 0; jj<NSieveRow; jj++){

	hole_no = Get_Hole(ii,jj);
	col_row = Get_Col_Row(hole_no);
	
	if(col_row[0] == ii && col_row[1] == jj){	
	  getline(cutcsvold,line);		
	  cutcsv << line << endl;
	  cout << "next csv line for hole " << hole_no << ", col =  " << ii << ", row = " << jj << endl;
	  cout << line << endl;
	}
	
      }
    }
	
  
  
    // ensure that old file is cycled through
    getline(cutcsvold,line);	 

  
    // try to find existing hole cut
    // use name e_hcut for ellipse_holecut
    cout << "Testing " << Form("e_hcut_L_%d_%d_%d", FoilID, col, k) << endl;
  
    TCutG* cutg = (TCutG*)gROOT->FindObject(Form("e_hcut_L_%d_%d_%d", FoilID, col, k));

  
    
    // if new ellipse is being drawn
    if(!cutg || overwrite){
      if(cutg)cutg->Delete();

  
      TVector3 Hole_pos = GetSieveHoleCorrectionTCS(FoilID,col, k);
            
      TVector3 MomDirectionTCS_hole = Hole_pos - BeamSpotTCS_average;
      
      Double_t theta_hole = MomDirectionTCS_hole.X()/MomDirectionTCS_hole.Z();
      Double_t phi_hole = MomDirectionTCS_hole.Y()/MomDirectionTCS_hole.Z();
      
      
      // TEllipse* Ellipse = new TEllipse(phi_hole,theta_hole,0.0008,0.0014,0.,360.,0); // larger sieve hole cuts
      TEllipse* Ellipse = new TEllipse(phi_hole,theta_hole,0.00038,0.0009,0.,360.,0); // smaller sieve hole cuts    
      Ellipse->SetFillStyle(0);
      Ellipse->SetLineWidth(5);
      Ellipse->Draw("same");
      

    
    
      // wait for graph object to be made in editor
      TGraph* g = (TGraph*)(TVirtualPad::Pad()->WaitPrimitive("Graph"));


      // get parameters of new ellipse
      
      double kPI = 3.14159265358979323846;
      Double_t x1 = Ellipse->GetX1(); Double_t y1 = Ellipse->GetY1();
      Double_t r1 = Ellipse->GetR1(); Double_t r2 = Ellipse->GetR2();
      Double_t theta = Ellipse->GetTheta();
      
	    
      // file in csv file with ellipse/ cut parameters
      cutcsv<<cur_hole<<","<<col<<","<<k<<","<<1<<","<<x1*1000<<","<<sieve_ph[col]<<","<<y1*1000<<","<<sieve_th[k]<<","<<r1*1000<<","<<r2*1000<<","<<theta<<endl;

      
    
      // Use parameters of ellipse to create new TCutG

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
      
      cutg = new TCutG("",n,x,y);
	    

      c1->Update();
	    
      cutg->SetName(Form("e_hcut_L_%d_%d_%d", FoilID, col, k));
	    
      // set axises' name
	    
      // cutg->SetVarX("L.tr.tg_ph");
      // cutg->SetVarY("L.tr.tg_th");

      // alternative use theta and phi-target from newer replay

      cutg->SetVarX("ph_tgt");
      cutg->SetVarY("th_tgt");

      // cutg->SetVarX("y_sieve");
      // cutg->SetVarY("x_sieve");


      cutg->Write("", TObject::kWriteDelete); // Overwrite old cut
      cout << "done!" << endl;
	    
    }    
    else {
      cout << Form("e_hcut_L_%d_%d_%d", FoilID, col, k) << " is found, using old one" << endl;

      
      
      cout << "Using line : " << endl;      
      cout << line << endl << endl;      
      cutcsv<<line<<endl;
      

      cutg->SetVarX("ph_tgt");
      cutg->SetVarY("th_tgt");

    }


    // draw TCutG (new or old if not overwritten)
    cutg->SetLineColor(kMagenta);	
    cutg->SetLineWidth(2);	
    cutg->Draw("PL");
    c1->Update();

    //    cutg->Write("", TObject::kWriteDelete); // Overwrite old cut
  

    // Log cuts names to text file   
 
    cutdesc << Form("e_hcut_L_%d_%d_%d", FoilID, col, k) << " && " << Form("fcut_L_%d", FoilID) << " && " << Form("fcut_L_FP_%d", FoilID) << " && " << (const char *) GenrealCut << endl;
  }

  cout << " added to cutdesc file " << endl;

 
  // puts remainder of old csv file into new
  if(cutcsvold){
    while(!cutcsvold.eof()){
      getline(cutcsvold,line);
      cutcsv<<line<<endl;
    }
  }
       
  f1->Write();
  f1->ls();
  cutcsv.close();
  cutdesc.close();
  

}


// used to make FP cuts where sieve hole cuts have already been made
void CutSieve_FP(int FoilID = 0, uint col = 0, int overwrite = 0, int append = 1){

  TChain* T = Load_more_rootfiles(Run_number, Run_number_2);

  ReLoadcuts();

  cout << "RootFileName = " << RootFileName << endl;
  cout << "Run_number = " << Run_number << endl;
  cout << "General cut = " << GenrealCut << endl << endl;

  // creating and filling expectation value of holes (theta and phi)
  Double_t sieve_ph[NSieveCol], sieve_th[NSieveRow];

  for(UInt_t r_count = 0; r_count < NSieveRow; r_count++){
    sieve_th[r_count] = 0;
  }

  for(UInt_t c_count = 0; c_count < NSieveCol; c_count++){
    sieve_ph[c_count] = 0;
  }




  // file that stores TCutGs (from foil and sieve)
  gSystem->Exec("cp -vf " + CutFileName + " " + CutFileName + ".old");
  TFile *f1 = new TFile(CutFileName, "UPDATE");
  assert(f1);
   (TCutG*) f1->GetObjectChecked(Form("fcut_L_%d", FoilID), "TCutG"); //looking for old cut definition

  // line added to extract previous FP vertex cut if existing
   (TCutG*) f1->GetObjectChecked(Form("fcut_L_FP_%d", FoilID), "TCutG"); //looking for old cut definition


  // Text file containing text names of all cuts (lines have TCutG names combined with general cuts used)
  fstream cutdesc;
    
  gSystem->Exec("cp -vf " + CutDescFileNameSieve_ell + " " +  CutDescFileNameSieve_ell + ".old");
  if(append){
    cutdesc.open(CutDescFileNameSieve_ell, ios_base::out |ios::app);
  }
  if(!append){
    cutdesc.open(CutDescFileNameSieve_ell, ios_base::out);
  }
  assert(cutdesc.is_open());
  



  // JW: added new canvas to create second cut for each hole on FP theta vs y.
   
   TCanvas* c2 = new TCanvas("c2","FP sieve",1000,1000);
     
    
  TH2D* thfp_v_yfp = new TH2D("thfp_v_yfp","thfp_v_yfp",300,-0.05,0.05,300,-35,30);
  // // TH2D* thfp_v_yfp = new TH2D("thfp_v_yfp","thfp_v_yfp",300,-0.2,0.2,300,-100,100);
  
  thfp_v_yfp->GetYaxis()->SetTitleOffset(1.0);
  thfp_v_yfp->GetXaxis()->SetTitleSize(0.05);
  thfp_v_yfp->GetYaxis()->SetTitleSize(0.05);
  thfp_v_yfp->GetXaxis()->SetTitle("y (FP) [m]");
  thfp_v_yfp->GetYaxis()->SetTitle("th (FP) [mrad]");
  
  // John W: attempt to add secondary cut on FP theta vs y. This acts on a secondary check on accidentaly including events wrongly identified from a different hole (or background).

    c2->cd(0);

    T->Draw("L.tr.r_th*1000:L.tr.r_y>>thfp_v_yfp",GenrealCut + TCut(Form("fcut_L_%d", FoilID)) + TCut(Form("fcut_L_FP_%d", FoilID)),"colz");


    c2->Update();

    int nhol = 0;
    cout << "How many holes in this No." << col << " column?" << endl;
    cin >> nhol;
    if (nhol < 0){
      exit(EXIT_FAILURE);
    }


    cout << "min hole in this column? : ";
    int rmin = -1;
    cin >> rmin;
    if (rmin < 0){
      exit(EXIT_FAILURE);
    }
    

    //make cuts for each holes 
    for (int cnt = 0; cnt < nhol; cnt++) {


      c2->cd(0);
      
      int k = rmin + cnt * 2; // line as holes in column are in every other row
    
      if (k < 0)
	continue;

      T->Draw("L.tr.r_th*1000:L.tr.r_y>>thfp_v_yfp",GenrealCut + TCut(Form("fcut_L_%d", FoilID)) + TCut(Form("fcut_L_FP_%d", FoilID)) + TCut(Form("e_hcut_L_%d_%d_%d", FoilID, col, k)),"colz");

           
      TCutG* tmpcut_FP = (TCutG*) gROOT->FindObject(Form("FPhcut_L_%d_%d_%d",FoilID, col, k));

      
      if (tmpcut_FP && !overwrite) {
      
	cout << Form("FPhcut_L_%d_%d_%d",FoilID, col, k)
	     << " is found, using old one" << endl;
	tmpcut_FP->Draw("PL");

	tmpcut_FP->SetVarX("L.tr.r_y[0]");
	tmpcut_FP->SetVarY("L.tr.r_th[0]*1000");
	//      tmpcut_FP->Write("",TObject::kWriteDelete);
      

      }
      else{

	cout << "Looking for FP hole cut... " << endl;

	c2->cd(0);

	TCutG* cutg_FP = (TCutG*) (gPad->WaitPrimitive("CUTG", "CutG"));

	c2->Update();

	cutg_FP->SetName(Form("FPhcut_L_%d_%d_%d",FoilID, col, k));

      
	cutg_FP->SetVarX("L.tr.r_y[0]");
	cutg_FP->SetVarY("L.tr.r_th[0]*1000");

	cout << "Pre file directory change" << endl;
	f1->cd();
	cutg_FP->Write("",TObject::kWriteDelete);
	tmpcut_FP = cutg_FP;
	cout << "reached end of loop" << endl;

      }



      //    tmpcut_FP =  (TCutG*) f1->GetObjectChecked(Form("FPhcut_L_%d_%d_%d",FoilID, col, k), "TCutG");
      assert(tmpcut_FP);
      tmpcut_FP->SetLineWidth(2);
      tmpcut_FP->SetLineColor(kMagenta);
      tmpcut_FP->Draw("PL");
      c2->Update();
 
    
      cout << "updated C2" << endl;


  

      // Log cuts names to text file   
 
      cutdesc << Form("FPhcut_L_%d_%d_%d",FoilID, col, k) << " && " << Form("e_hcut_L_%d_%d_%d", FoilID, col, k) << " && " << Form("fcut_L_%d", FoilID) << " && " << Form("fcut_L_FP_%d", FoilID) << " && " << (const char *) GenrealCut << endl;

    }

    cout << " added to cutdesc file " << endl;
    
    


    

    f1->Write();
    f1->ls();
    cutdesc.close();
  

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




  //  gSystem->Exec("cp -vf " + CutFileName + " " + CutFileName + ".old");
  
  cout << "CutFileName = " << CutFileName << endl;



  // CutFileName = "/w/work3/home/johnw/Rootfiles/apex_4771_27_11_2019.root.FullCut.root_13_12_19";

  // cout << "For this function only change to -> " << CutFileName << endl;

  TFile *f1 = new TFile(CutFileName, "READ");
  assert(f1);


  Int_t first_foil = 0;
  Int_t last_foil = 8;
  Int_t offset = 0;

  if( Run_number == 4771){
    first_foil = 1;
    offset = 1;
    cout << "Run_number == 4771 condition" << endl;
  }
  else if( Run_number == 4774){
    first_foil = 0;
    cout << "Run_number == 4774 condition" << endl;
  }
  



  TCutG* Foil_cuts[4] = {NULL};

  Int_t Foil_it = 0; // iterate through 4 chosen foils
  //  for(Int_t FoilID = 1; FoilID<8; FoilID=FoilID+2){
  for(Int_t FoilID = first_foil; FoilID<last_foil; FoilID=FoilID+2){
    cout << "FoilID = " << FoilID << endl;
    Foil_cuts[Foil_it] = (TCutG*) f1->GetObjectChecked(Form("fcut_L_%d", FoilID), "TCutG");
    (TCutG*) f1->GetObjectChecked(Form("fcut_L_FP_%d", FoilID), "TCutG");
    Foil_it++;
  }

  
  
  //  TCanvas *c1 = new TCanvas("c1","PlotSieve",1000,1000);

  // first canvases to show theta-phi plots for each foil
  TCanvas *c1[4];

  // second canvases to show X-Y plots for each foil
  TCanvas *c2[4];


  // third canvas to show FP distribution of sieve cuts
  TCanvas *c3[4];


  // fourth canvas to show FP distribution of sieve cuts with added z-scale of theta-target value
  
  TCanvas *c_FP_theta[4];


  // comparison FP canvas (attempt to show overlap between adjacent foils: 0(1) and 2(3) etc)

  TCanvas *c4 = new TCanvas("c4","FP comparison between foils");
  
  
  // show reactz vs theta_tg foil cuts

  TCanvas *c5 = new TCanvas("c5","FP comparison between foils");
  
  



  TRotation fTCSInHCS;
  TVector3 TCSX(0,-1,0);
  TVector3 TCSZ(TMath::Sin(HRSAngle),0,TMath::Cos(HRSAngle));
  TVector3 TCSY = TCSZ.Cross(TCSX);
  fTCSInHCS.RotateAxes(TCSX,TCSY,TCSZ);
  
  fPointingOffset.SetXYZ(-MissPointZ*TMath::Sin(HRSAngle)*TMath::Cos(HRSAngle),(Double_t)MissPointY,MissPointZ*TMath::Sin(HRSAngle)*TMath::Sin(HRSAngle));
  



  // Double_t BeamX_average = 0;
  // Double_t BeamY_average = 0;
  
  
  // BeamX_average = -0.0004606;
  // BeamY_average = 0.002448;  

  

  // 4774 beamx and beamy unrastered

   // BeamX_average = -0.0006391;
   // BeamY_average = 0.002405;  
  

  // 4771 beamx and beamy unrastered
  // BeamX_average = -0.0006361;
  // BeamY_average = 0.002419; 

  


  TH2F* h2[4] = {NULL};

  TH2F* h3[4] = {NULL};

  Double_t x_lim = 1.3*TMath::Max(TMath::Abs(SieveYbyCol[0]),TMath::Abs(SieveYbyCol[NSieveCol-1]));
  Double_t y_lim = 1.3*TMath::Max(TMath::Abs(SieveXbyRow[0]),TMath::Abs(SieveXbyRow[NSieveRow-1]));
 


  TH2F* h_FP[4] = {NULL};

  TH2F* h_FP_theta[4] = {NULL};


  TH2F* h_FP_comp[7] = {NULL};



  TH2F* h_Foil_comp = new TH2F("h_Foil_comp", "ReactZ vs. Target Phi", 400, -0.05, 0.105, 200,-0.4, 0.4);



 
  

  for(Int_t i = 0; i<4; i++){

    
    Int_t foil_no = (i*2)+offset;
    h2[i] = new TH2F(Form("h2_%d",i), Form("theta_target vs. phi_target Foil_%d",foil_no), 300, -0.04, 0.04, 300, -0.08, 0.08);
    h2[i]->SetMinimum(2);

    c1[i] = new TCanvas(Form("c1_%d",i),Form("PlotSieve foil_%i",foil_no),1000,1000);


    h3[i] = new TH2F(Form("h3_%d",i),Form("Sieve Plane Proj. (tg_X vs tg_Y) for Foil #%d",foil_no), 300,-x_lim,x_lim, 300,-y_lim,y_lim);

    c2[i] = new TCanvas(Form("c2_%d",i),Form("PlotXYSieve foil_%i",foil_no),1000,1000);



    h_FP[i] = new TH2F(Form("h_FP_%d",i),Form("Foil cuts on ph vs y (FP) for Foil #%d",foil_no), 400, -0.05, 0.05, 200,-0.04, 0.04);
    
    c3[i] = new TCanvas(Form("c3_%d",i),Form("FP Plot_%i",foil_no),1000,1000);


    //    h_FP_theta[i] = new TH2F(Form("h_FP_theta_%d",i),Form("Foil cuts on ph vs y (FP) for Foil #%d with theta_tg",foil_no), 400, -0.05, 0.05, 200,-0.04, 0.04, 200,-0.06,0.06);


    h_FP_theta[i] = new TH2F(Form("h_FP_theta_%d",i),Form("Foil cuts on ph vs y (FP) for Foil #%d with theta_tg",foil_no), 400, -0.05, 0.05, 200,-0.04, 0.04);
    
    c_FP_theta[i] = new TCanvas(Form("c_FP_ theta_%d",i),Form("FP Plot on Foil #%d with theta_tg",foil_no),1000,1000);
    

    
  }

  

  //  h_FP_comp[0] = new  TH2F(Form("h_FP_comp_%d",0),Form("Foil cuts on ph vs y (FP)",0+offset,2+offset), 400, -0.05, 0.05, 200,-0.04, 0.04);

  for (Int_t i = 0; i<7; i++){
    
    h_FP_comp[i] = new  TH2F(Form("h_FP_comp_%d",i),"Foil cuts on ph vs y (FP)", 60, -0.05, 0.05, 60,-0.04, 0.04);
    
  }



  for(Int_t i = 0; i<4; i++){

    c1[i]->cd(0);
    
    Int_t foil_no = (i*2)+offset;
    
    //    T->Draw(Form("th_tgt:ph_tgt>>h2_%d",i), GenrealCut + TCut(Form("fcut_L_%d", foil_no)) + TCut(Form("fcut_L_FP_%d", foil_no)), "COLZ");
    //    T->Draw(Form("th_tgt:ph_tgt>>h2_%d",i), GenrealCut + TCut(Form("fcut_L_FP_%d", foil_no)), "COLZ");


    // new optimisation
    T->Draw(Form("th_tgt:ph_tgt>>h2_%d",i), GenrealCut + TCut(Form("fcut_L_%d", foil_no)), "COLZ");

    // old optimisation
    //    T->Draw(Form("L.tr.tg_th:L.tr.tg_ph>>h2_%d,i"), GenrealCut + TCut(Form("fcut_L_%d",foil_no)), "COLZ");

    // // also draw red crosses where holes 'should' be


      

      
    //    TVector3 BeamSpotHCS_average(BeamX_average,BeamY_average,targetfoils[foil_no]);

    const TVector3 BeamSpotHCS_average(BeamX_average[foil_no] + (targetfoils[foil_no]/BeamZDir_average)*BeamXDir_average[foil_no], BeamY_average[foil_no] + (targetfoils[foil_no]/BeamZDir_average)*BeamYDir_average[foil_no], targetfoils[foil_no]);


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
      
      TVector3 Hole_pos = GetSieveHoleCorrectionTCS(i,Hole);
      
      //      TVector3 Hole_pos = GetSieveHoleTCS(Hole);
      
      
      TVector3 MomDirectionTCS_hole = Hole_pos - BeamSpotTCS_average;
      
      Double_t theta_hole = MomDirectionTCS_hole.X()/MomDirectionTCS_hole.Z();
      Double_t phi_hole = MomDirectionTCS_hole.Y()/MomDirectionTCS_hole.Z();
      

      // cout << "For Hole " << Hole << " : Hole_pos = [" << Hole_pos.X() << "," << Hole_pos.Y() << "," << Hole_pos.Z() << "]" << endl;
      
      // cout << "For BeamSpotTCS_average " << Hole << " : BeamSpotTCS_average = [" << BeamSpotTCS_average.X() << "," << BeamSpotTCS_average.Y() << "," << BeamSpotTCS_average.Z() << "]" << endl;

      //      cout << "For Theta " << Hole << " : theta_hole = [" << theta_hole.X() << "," << theta_hole.Y() << "," << theta_hole.Z() << "]" << endl;


      
      
      // + type crosses
      // TLine *lh = new TLine(posx-plotwidth,posy,posx+plotwidth,posy);
      // TLine *lv = new TLine(posx,posy-plotwidth,posx,posy+plotwidth);
	
      
      //Saltire-type crosses
      TLine *lc1 = new TLine(phi_hole-plotwidth,theta_hole-plotwidth,phi_hole+plotwidth,theta_hole+plotwidth);
      lc1->SetLineWidth(width);
      TLine *lc2 = new TLine(phi_hole-plotwidth,theta_hole+plotwidth,phi_hole+plotwidth,theta_hole-plotwidth);

      // TLine *lc1 = new TLine(phi_hole - plotwidth, theta_hole, phi_hole + plotwidth, theta_hole);
      // lc1->SetLineWidth(width);
      // TLine *lc2 = new TLine(phi_hole, theta_hole - plotwidth, phi_hole, theta_hole + plotwidth);
      // lc2->SetLineWidth(width);



      lc2->SetLineWidth(width);
      
      
	
      lc1->SetLineColor(color);
      lc2->SetLineColor(color);
      
      lc1->Draw("same");
      lc2->Draw("same");
	
      
    }
    
    c2[i]->cd(0);
      

    //    T->Draw(Form("x_sieve:y_sieve_alt>>h3_%d",i), GenrealCut + TCut(Form("fcut_L_%d", foil_no)) + TCut(Form("fcut_L_FP_%d", foil_no)), "COLZ");

    T->Draw(Form("x_sieve:y_sieve_alt>>h3_%d",i), GenrealCut + TCut(Form("fcut_L_%d", foil_no)), "COLZ");


    // T->Draw(Form("x_sieve:y_sieve_alt>>h3_%d",i), GenrealCut + TCut(Form("fcut_L_FP_%d", foil_no)), "COLZ");


      
    for(UInt_t Hole = 0; Hole < NHoles; Hole++){

	
      // Bool_t redr = Row % 2 ==0;
      // Bool_t redc = Col==0 or Col==2 or Col==5 or Col==8 or Col==11;
      
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
      
      // if (redr and redc) color=kRed;
      // else if (!redr and !redc) color=kBlue;
      // else continue;
      
      
      // std::vector<int> x_y = {};
      // x_y = Get_Col_Row(j);
      
      // UInt_t Col = x_y[0];
      // UInt_t Row = x_y[1];
      
	
	
      // const Double_t posx = SieveOffY + SieveYbyCol[Col];
      // const Double_t posy = SieveOffX + SieveXbyRow[Row];
	
      TVector3 Hole_pos = GetSieveHoleTCS(Hole);
	
      Double_t posy = Hole_pos.X();
      Double_t posx = Hole_pos.Y();
      
      //				cout << "For Hole: " << Hole << ": posy = " << posy << " and posx = " << posx << endl;
	
	
      // + type crosses
      // TLine *lh = new TLine(posx-plotwidth,posy,posx+plotwidth,posy);
      // TLine *lv = new TLine(posx,posy-plotwidth,posx,posy+plotwidth);
	
	
      // Saltire-type crosses
      TLine *lc1 = new TLine(posx-plotwidth,posy-plotwidth,posx+plotwidth,posy+plotwidth);
      lc1->SetLineWidth(width);
      TLine *lc2 = new TLine(posx-plotwidth,posy+plotwidth,posx+plotwidth,posy-plotwidth);
      lc2->SetLineWidth(width);
	
	
	
	
      // lh -> SetLineColor(color);
      // lv -> SetLineColor(color);
	
      // lh -> Draw();
      // lv -> Draw();
	
      lc1->SetLineColor(color);
      lc2->SetLineColor(color);
      
      lc1->Draw("same");
      lc2->Draw("same");
      
    }
    

    c3[i]->cd(0);
    
    T->Draw(Form("L.tr.r_ph:L.tr.r_y>>h_FP_%d",i), GenrealCut + TCut(Form("fcut_L_%d", foil_no)), "COLZ");
    // T->Draw(Form("L.tr.r_ph:L.tr.r_y>>h_FP_%d",i), GenrealCut + TCut(Form("fcut_L_FP_%d", foil_no)), "COLZ");


    c_FP_theta[i]->cd(0);
    
    
    T->Draw(Form("L.tr.r_ph:L.tr.r_y:th_tgt>>h_FP_theta_%d",i), GenrealCut + TCut(Form("fcut_L_%d", foil_no)), "COLZ");
      
      
  }
    


  // perform comparison between adjacent foil FP plots

  //  c4->Divide(2,1);

  c4->cd(1);
  
  
  //  T->SetMarkerColor(kRed);

  
  h_FP_comp[0]->SetMarkerColor(kBlue);

  cout << "cut 1 = " << TCut(Form("fcut_L_%d", 0+offset)) + TCut(Form("!fcut_L_%d", 2+offset)) << endl;

  T->Draw(Form("L.tr.r_ph:L.tr.r_y>>h_FP_comp_%d",0), GenrealCut + TCut(Form("fcut_L_%d", 0+offset)) +  TCut(Form("!fcut_L_%d", 2+offset)) + TCut(Form("!fcut_L_%d", 4+offset)) + TCut(Form("!fcut_L_%d", 6+offset)) ,  "");



  h_FP_comp[1]->SetMarkerColor(kMagenta+3);

  //  T->Draw(Form("L.tr.r_ph:L.tr.r_y>>h_FP_comp_%d",0), GenrealCut + TCut(Form("!fcut_L_%d", 0+offset)) +  TCut(Form("fcut_L_%d", 2+offset)),  "SCAT");


  //  cout << "cut 2 = " << Form("!fcut_L_%d", 0+offset) + Form("fcut_L_%d", 2+offset) << endl;

  T->Draw(Form("L.tr.r_ph:L.tr.r_y>>h_FP_comp_%d",1), GenrealCut + TCut(Form("fcut_L_%d", 0+offset)) +  TCut(Form("fcut_L_%d", 2+offset)) + TCut(Form("!fcut_L_%d", 4+offset)) + TCut(Form("!fcut_L_%d", 6+offset)),  "SCAT SAME");


  h_FP_comp[2]->SetMarkerColor(kRed);


  T->Draw(Form("L.tr.r_ph:L.tr.r_y>>h_FP_comp_%d",2), GenrealCut + TCut(Form("!fcut_L_%d", 0+offset)) +  TCut(Form("fcut_L_%d", 2+offset)) + TCut(Form("!fcut_L_%d", 4+offset)) + TCut(Form("!fcut_L_%d", 6+offset)),  "SCAT SAME");



  //  draw overlap between 2nd and 3rd foils in focal plane

  h_FP_comp[3]->SetMarkerColor(kOrange);  
  
  T->Draw(Form("L.tr.r_ph:L.tr.r_y>>h_FP_comp_%d",3), GenrealCut + TCut(Form("!fcut_L_%d", 0+offset)) +  TCut(Form("fcut_L_%d", 2+offset)) + TCut(Form("fcut_L_%d", 4+offset)) + TCut(Form("!fcut_L_%d", 6+offset)),  "SCAT SAME");


  // 3rd foil

  h_FP_comp[4]->SetMarkerColor(kYellow);  
  
  T->Draw(Form("L.tr.r_ph:L.tr.r_y>>h_FP_comp_%d",4), GenrealCut + TCut(Form("!fcut_L_%d", 0+offset)) +  TCut(Form("!fcut_L_%d", 2+offset)) + TCut(Form("fcut_L_%d", 4+offset)) + TCut(Form("!fcut_L_%d", 6+offset)),  "SCAT SAME");


  
  // 3rd and 4th foil overlap
  
  h_FP_comp[5]->SetMarkerColor(kGreen + 3);  
  
  T->Draw(Form("L.tr.r_ph:L.tr.r_y>>h_FP_comp_%d",5), GenrealCut + TCut(Form("!fcut_L_%d", 0+offset)) +  TCut(Form("!fcut_L_%d", 2+offset)) + TCut(Form("fcut_L_%d", 4+offset)) + TCut(Form("fcut_L_%d", 6+offset)),  "SCAT SAME");


  // 4th foil


  h_FP_comp[6]->SetMarkerColor(kBlue + 3);  
  
  T->Draw(Form("L.tr.r_ph:L.tr.r_y>>h_FP_comp_%d",6), GenrealCut + TCut(Form("!fcut_L_%d", 0+offset)) +  TCut(Form("!fcut_L_%d", 2+offset)) + TCut(Form("!fcut_L_%d", 4+offset)) + TCut(Form("fcut_L_%d", 6+offset)),  "SCAT SAME");



  // draw Focal-Plane ph vs y with all foil cuts as shapes

  ///  c4->cd(2);




  // perform comparison between adjacent foil reactz vs phi_tg cuts

  c5->cd(1);


  T->Draw("reactz:ph_tgt>>h_Foil_comp", GenrealCut, "COLZ"); 
  

  Int_t colours[] = {kRed,kBlue,kGreen,kOrange};

  for(Int_t i = 0; i<4; i++){
    
    Foil_cuts[i]->Draw("same");
    
    Foil_cuts[i]->SetLineColor(colours[i]);	
    Foil_cuts[i]->SetLineWidth(4);	
    Foil_cuts[i]->Draw("PL");
    

  }
  
  



    
  TString save_file = ROOTFILE_DIR;
  
  save_file.Remove(save_file.Last('/'),15);
  save_file.Remove(save_file.First('/'),43);
  
  cout << "savefile = " << save_file << endl;
    
  save_file = "../Record/" + save_file;
  
  cout << "savefile = " << save_file << endl;
    
  //    gSystem->Exec("mkdir Record/" + save_file); 
  gSystem->Exec("mkdir " +  save_file); 


    


  for(Int_t i = 0; i<4; i++){
    
    Int_t foil_no = (i*2)+offset;
    c1[i]->Print(save_file + Form("/Foil_%d_th_ph.png",foil_no));
    c2[i]->Print(save_file + Form("/Foil_%d_sieve_x_y.png",foil_no));

    c3[i]->Print(save_file + Form("/Foil_%d_FP_ph_y.png",foil_no));


  }
    

  c4->Print(save_file + "/Foil_all_foils_FP_ph_y.png");


  c5->Print(save_file + "/Foil_all_foils_zreact_v_phi_tg.png");
	

    
      
}




void CutSieve(int i = 0 /*Foil ID*/, int cmin = 0 /*Column ID*/, int overwrite = 0, bool append = true) {

  int ncol = 1;

  ReLoadcuts();

  //  gROOT->ProcessLine(".!display Example_Sieve.png;");

  cout << "Line reached" << endl;

  //	TCanvas * c1 = new TCanvas(Form("PlotSieve%d_%d", i, cmin), "PlotSieve",1800, 1150);
  TCanvas * c1 = new TCanvas(Form("PlotSieve%d_%d", i, cmin), "PlotSieve",1000, 1000);
  c1->Update();
  
  // TFile* f = new TFile(SoureRootFile);
  // TTree * T = (TTree *) f->GetObjectChecked("T", "TTree");
  TChain* T = Load_more_rootfiles(Run_number, Run_number_2);

  gSystem->Exec("cp -vf " + CutFileName + " " + CutFileName + ".old");
  TFile *f1 = new TFile(CutFileName, "UPDATE");
  assert(f1);
  (TCutG*) f1->GetObjectChecked(Form("fcut_L_%d", i), "TCutG"); //looking for old cut definition
  
  // make cuts for each holes
   // TH2F* h2 = new TH2F("h2", "theta_target vs. phi_target", 300, -0.045, 0.01, 300, 0.0, 0.07);
   // h2->SetMinimum(3);
   TH2F* h2 = new TH2F("h2", "theta_target vs. phi_target", 900, -0.04, 0.04, 900, -0.08, 0.08);
   h2->SetMinimum(2);
   

  T->Draw("L.tr.tg_th:L.tr.tg_ph>>h2", GenrealCut
	  + TCut(Form("fcut_L_%d", i)), "COLZ");
  c1->Update();

  //	int ncol = 1, cmin=1;


  fstream cutdesc;


  if(append){
    cutdesc.open(CutDescFileNameSieve_dir, ios_base::out |ios::app);
  }
  if(!append){
    cutdesc.open(CutDescFileNameSieve_dir, ios_base::out);
  }
  assert(cutdesc.is_open());



  
  // make cuts from the first available column
  for (int j = cmin; j < ncol + cmin; j++) //column number j
    {
      
      int nhol = 0;
      cout << "How many holes in this No." << j << " column?" << endl;
      cin >> nhol;
      if (nhol < 0)
	continue;
      cout << "min hole id : ";
      int rmin = -1;
      cin >> rmin;
      if (rmin < 0)
	continue;
      
      //make cuts for each holes from begin number
      for (int cnt = 0; cnt < nhol; cnt++) {
	int k = rmin + cnt * 2;
	//			cout<<"Hole ID : "<<k<<endl;
	//			cin>>k;
	if (k < 0)
	  continue;
	
	TCutG* tmpcut = (TCutG*) gROOT->FindObject(Form("hhcut_L_%d_%d_%d",i, j, k));
	//WARNING: Overwriting
	//			tmpcut = 0;
	//			cout<<"WARNING: Overwiting!"<<endl;
	//WARNING: Overwriting
	if (tmpcut && !overwrite) {
	  
	  cout << Form("hhcut_L_%d_%d_%d", i, j, k)
	       << " is found, using old one" << endl;
	} else {
	  //				TPad * gPad = TVirtualPad::Pad();
	  
	  // define cuts
	  cout << "making cut for foil=" << i << ", column=" << j
	       << ", No." << k << " hole: waiting ..." << endl;
	  TCutG* cutg = (TCutG*) gPad->WaitPrimitive("CUTG", "CutG");
	  c1->Update();
	  cout << "done!" << endl;
	  cutg->SetName(Form("hhcut_L_%d_%d_%d", i, j, k));
	  
	  // set axises
	  cutg->SetVarX("L.tr.tg_ph[0]");
	  cutg->SetVarY("L.tr.tg_th[0]");
	  
	  // output
	  f1->cd();
	  cutg->Write("", TObject::kWriteDelete);
	  tmpcut = cutg;
	  //			f->cd();
	}
	
	tmpcut = (TCutG*) f1->GetObjectChecked(Form("hhcut_L_%d_%d_%d", i, j, k), "TCutG");
	assert(tmpcut);
	tmpcut->SetLineWidth(2);
	tmpcut->SetLineColor(kMagenta);
	tmpcut->Draw("PL");
	c1->Update();


	cout << "log to " << CutDescFileNameSieve_dir << endl;
	cutdesc << Form("fcut_L_%d && hhcut_L_%d_%d_%d", i, i, j, k) << " && " << (const char *) GenrealCut
		<< endl;
      }
    }


  cutdesc.close();
  f1->Write();
  SaveCanvas(c1, *SoureRootFile + Form("hhcut_L_%d_%d", i, cmin), kFALSE);
  //	SaveCanvas(c1, SoureRootFile + Form("hcut_L_%d_%d_%d", i, cmin), kFALSE);
  //	f1->Close();
  f1->ls();
  
}

void cut_L_dp(int overwrite = 0)
{

//	fstream cutdesc(CutDescFileName, ios_base::out);
//	assert(cutdesc.is_open());
//
//	fstream cutdescSieve(CutDescFileName, ios_base::in);
//	assert(cutdescSieve.is_open());


	TFile* f = new TFile(*SoureRootFile);
	TTree * t = (TTree *) f->GetObjectChecked("T", "TTree");

//	gSystem->Exec("cp -vf ./cutfiles/gcut_L_dpfull.root ./cutfiles/gcut_L_dpfull.root.old2");
	TFile* fcut = new TFile(CutFileName,"UPDATE");fcut->cd();

	TCanvas * c1 = new TCanvas("dpCuts","asdf",1700,900);

	TCutG * cutg = NULL;

	for(UInt_t FoilID = 0;FoilID<10;FoilID++)
	{

		TString VertexCut = Form("fcut_L_%d",FoilID);
		cout<<"Working on Vertex Cut "<<VertexCut<<endl;

		TCutG *cVertexCut = (TCutG*)fcut->GetObjectChecked(VertexCut,"TCutG");
		if( !cVertexCut) break;



		t->Draw(Form("L.tr.tg_th:L.tr.tg_dp>>hDpCutFoil%d(1000,-.01,.05,1000,-.1,.1)",FoilID)
				,TCut(VertexCut) && GenrealCut,"COLZ");
		c1->Update();

		cutg = (TCutG *)fcut->GetObjectChecked(Form("L_dp_%d",FoilID),"TCutG");
		assert(cutg);

		if (!cutg || overwrite)
		{

			cutg = (TCutG*) (TVirtualPad::Pad()->WaitPrimitive("CUTG","CutG")); // making cut, store to CUTG
			c1->Update();

			cutg->SetName(Form("L_dp_%d",FoilID)); //
		  // set axises' name
			cutg->SetVarX("L.tr.tg_dp[0]");
			cutg->SetVarY("L.tr.tg_th[0]");
		}
		else {cout<<"Found "<<Form("L_dp_%d",FoilID)<<endl;}

		cutg->SetLineColor(kMagenta);
		cutg->SetLineWidth(2);
		cutg->Draw("PL");

		cutg->Write("",TObject::kWriteDelete); // Overwrite old cut

		SaveCanvas(c1, *SoureRootFile + Form("L_dp_%d",FoilID), kFALSE);
	}

	fcut->Write();

//	cutdesc.close();
//	cutdescSieve.close(); sd
}

