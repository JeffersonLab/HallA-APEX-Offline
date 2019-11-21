#include "SaveCanvas.C"
#include "TPad.h"
#include "InputAPEXL.h"
#include "Load_rootfile.C"
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





void ReLoadcuts(){

  

  
  //  GenrealCut = TCut("abs(L.tr.tg_dp)<0.01") && GeneralSieveCut;


  Double_t Dp_values[11];

  for (Int_t i = 0; i<sizeof(Dp_values)/sizeof(Dp_values[0]); i++){
    
    Dp_values[i] = -0.05 + i*0.01;
    //    std::cout << "Test: Dp_values[" << i <<"] = " << Dp_values[i] << std::endl;
  }

  Double_t dp_lim1 = -0.01;
  Double_t dp_lim2 = 0.01;
  
  for(Int_t i = 0; i < sizeof(Dp)/sizeof(Dp[0]); i++){
    
    if(Dp[i]){
      dp_lim1 =  Dp_values[i];
      dp_lim2 =  Dp_values[i+1];
      
    }
    
  }
 
  //  GenrealCut = TCut(Form("L.tr.n==1 && L.tr.tg_dp>%.2f && L.tr.tg_dp<%.2f",dp_lim1,dp_lim2)) && GeneralSieveCut;

  //std::cout << "GenrealCut = " << GenrealCut << std::endl;

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
 
// TString SoureRootFile_dp = "/w/work3/home/johnw/Rootfiles/apex_4179_1_7_2019_dp_%d.root";
  
  // SoureRootFile_dp = 
  
  for(Int_t i = 0; i < sizeof(Dp)/sizeof(Dp[0]); i++){
    
    if(Dp[i]){
      *SoureRootFile = Form(SoureRootFile_dp,i);
    }
    
  }
 

  CutFileName = *SoureRootFile + ".FullCut.root";
  CutDescFileName = *SoureRootFile + ".VertexCut.cut";
  CutDescFileNameSieve = *SoureRootFile + ".SieveCut.cut";
  CutDescFileNameSieve_dir = *SoureRootFile + ".SieveCut_dir.cut";


  CutDescFileNameDp = *SoureRootFile + ".DpCut.cut";

  CutDescFileNameCol = *SoureRootFile + ".ColCut.cut";
  
}

void PlotTest() {

//	TCut HoleCut = "abs(L.tr.tg_ph+.009)<.0015 && abs(L.tr.tg_th-.0265)<.005";
	TCut HoleCut = "1";

//	T->Draw("L.tr.tg_th:L.tr.tg_ph>>tpasdf2(500,-.04,.04,500,-.08,.08)",
//			GenrealCut, "COLZ");

	TFile* f = new TFile(*SoureRootFile);
	assert(f);
	TTree * T = (TTree *) f->GetObjectChecked("T", "TTree");
	assert(T);


	TCanvas * c1 = new TCanvas("PlotTest", "PlotTest", 1800, 900);
	c1->Divide(4, 3);
	int idx = 1;

	c1->cd(idx++);
	c1->Update();
	T->Draw("L.tr.x:L.tr.y>>xyasdf(3000,-.1,.1,3000,-1,1)", GenrealCut && HoleCut, "COLZ");
	// 	xyasdf->SetTitle("VDC Hits x VS y",GenrealCut && HoleCut,"COLZ");

	c1->cd(idx++);
	c1->Update();
	T->Draw("L.tr.th:L.tr.ph>>xyasdfasdf(3000,-.06,.06,3000,-.2,.2)", GenrealCut && HoleCut,
			"COLZ");
	// 	xyasdfasdf->SetTitle("VDC Hits th VS ph",GenrealCut && HoleCut,"COLZ");

	c1->cd(idx++);
	c1->Update();
	T->Draw("L.tr.x:L.tr.th>>xyasdf2(3000,-.2,.2,3000,-1,1)", GenrealCut && HoleCut, "COLZ");
	// 	xyasdf->SetTitle("VDC Hits x VS y",GenrealCut && HoleCut,"COLZ");

	c1->cd(idx++);
	c1->Update();
	T->Draw("L.tr.y:L.tr.ph>>xyasdfasdf2(3000,-.06,.06,3000,-.1,.1)", GenrealCut && HoleCut,
			"COLZ");

	c1->cd(idx++);
	c1->Update();
	T->Draw("L.tr.x:L.tr.ph>>xyasdf2asdfasdf(3000,-.06,.06,3000,-1,1)", GenrealCut && HoleCut, "COLZ");

	c1->cd(idx++);
	c1->Update();
	T->Draw("L.tr.y:L.tr.th>>xyasdfasdf2asdfw2re2(3000,-.2,.2,3000,-.1,.1)", GenrealCut && HoleCut,
			"COLZ");

//	//////////////////////////////////////////////////////////////////////////
////	c1->cd(idx++);
////	c1->Update();
////	T->Draw("L.tr.tg_th:L.tr.tg_ph>>tpasdf1(500,-.04,.04,500,-.08,.08)",
////			GenrealCut && HoleCut, "COLZ");
//
//	c1->cd(idx++);
//	c1->Update();
//	T->Draw("L.tr.tg_th:L.tr.x>>tpasdf2(3000,-1,1,500,-.08,.08)",
//			GenrealCut && HoleCut, "COLZ");
//
//	c1->cd(idx++);
//	c1->Update();
//	T->Draw("L.tr.tg_th:L.tr.y>>tpasdf3(3000,-.1,.1,500,-.08,.08)",
//			GenrealCut && HoleCut, "COLZ");
//
//	c1->cd(idx++);
//	c1->Update();
//	T->Draw("L.tr.tg_th:L.tr.th>>tpasdf4(3000,-.06,.06,500,-.08,.08)",
//			GenrealCut && HoleCut, "COLZ");
//
//	c1->cd(idx++);
//	c1->Update();
//	T->Draw("L.tr.tg_th:L.tr.ph>>tpasdf5(3000,-.06,.06,500,-.08,.08)",
//			GenrealCut && HoleCut, "COLZ");
//
//	///////////////////////////////////////////////////////////////////////////
//
////	c1->cd(idx++);
////	c1->Update();
////	T->Draw("L.tr.tg_th:L.tr.tg_ph>>tpasdf1(500,-.04,.04,500,-.08,.08)",
////			GenrealCut && HoleCut, "COLZ");
//
//
//	c1->cd(idx++);
//	c1->Update();
//	T->Draw("L.tr.tg_ph:L.tr.x>>tpasdfasdf3(3000,-1,1,500,-.04,.04)",
//			GenrealCut && HoleCut, "COLZ");
//
//	c1->cd(idx++);
//	c1->Update();
//	T->Draw("L.tr.tg_ph:L.tr.y>>tpasdfaasdfsdf3(3000,-.1,.1,500,-.04,.04)",
//			GenrealCut && HoleCut, "COLZ");
//
//	c1->cd(idx++);
//	c1->Update();
//	T->Draw("L.tr.tg_ph:L.tr.th>>tpasdfasdf4(3000,-.06,.06,500,-.04,.04)",
//			GenrealCut && HoleCut, "COLZ");
//
//	c1->cd(idx++);
//	c1->Update();
//	T->Draw("L.tr.tg_ph:L.tr.ph>>tpasdfasdf5(3000,-.06,.06,500,-.04,.04)",
//			GenrealCut && HoleCut, "COLZ");



}

void PlotRun()
{

	TFile* f = new TFile(*SoureRootFile);
	assert(f);
	TTree * T = (TTree *) f->GetObjectChecked("T", "TTree");
	assert(T);


	cout << *SoureRootFile << endl << (const char *) GenrealCut << endl;
	TCanvas * c1 = new TCanvas("PlotRun", "PlotRun", 1800, 900);
	c1->Divide(3, 2);
	int idx = 1;

	c1->cd(idx++);
	c1->Update();
	T->Draw("L.tr.x:L.tr.y>>xyasdf(500,-.1,.1,500,-1,1)", GenrealCut , "COLZ");
	// 	xyasdf->SetTitle("VDC Hits x VS y",GenrealCut,"COLZ");

	c1->cd(idx++);
	c1->Update();
	T->Draw("L.tr.th:L.tr.ph>>xyasdfasdf(3000,-.06,.06,3000,-.2,.2)", GenrealCut ,
			"COLZ");
	// 	xyasdfasdf->SetTitle("VDC Hits th VS ph",GenrealCut,"COLZ");

	c1->cd(idx++);
	c1->Update();
	T->Draw("L.tr.chi2>>dpasdfasdf(6000,0,0.01)", GenrealCut );
	// 	dpasdf->SetTitle("dp");
	gPad->SetLogy();

	c1->cd(idx++);
	c1->Update();
	T->Draw("L.tr.tg_th:L.tr.tg_ph>>tpasdf(500,-.04,.04,500,-.08,.08)",
			GenrealCut , "COLZ");
	// 	tpasdf->SetTitle("Sieve");

	c1->cd(idx++);
	c1->Update();
	T->Draw("L.tr.tg_dp>>dpasdf(6000,-.06,.06)", GenrealCut );
	// 	dpasdf->SetTitle("dp");

	c1->cd(idx++);
	c1->Update();
	T->Draw("ReactPt_L.z:L.tr.tg_ph>>dasdfasdf(500,-.04,.04,500,-.6,.6)",
			GenrealCut , "COLZ");
	// 	dasdfasdf->SetTitle("VZ vs tg_ph (each strip is a foil)");

	
	TString Canvas_name = *SoureRootFile;
	Canvas_name =+  "." + TString(c1->GetName());


	Bool_t canvas_bool = kFALSE;
	//	SaveCanvas(c1, *SoureRootFile + "." + TString(c1->GetName()),kFALSE);
	SaveCanvas(c1,Canvas_name,canvas_bool);
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


void cut_Vertex(int overwrite = 0, int nfoils = 3, int FoilID = -1) {


  cout << "cut being used " << GenrealCut << endl << endl;

  cout << "open " << *SoureRootFile << endl;

	ReLoadcuts();

	// TFile* f = new TFile(SoureRootFile);
	// assert(f);
	// TTree * T = (TTree *) f->GetObjectChecked("T", "TTree");
	// assert(T);
	
	//	TChain* T = Load_rootfile(Run_number);
	TChain* T = Load_more_rootfiles(Run_number, Run_number_2);



	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// loop for alternative z-react calculation


	TH1F *hOpZ  = new TH1F("hOpz","Z-vertex as from Optimisation",100,-1,1);
       
	Double_t Lrb_x;      //  Lrb.x rastered beam x
	Double_t L_tr_tg_y[100];  //  L.tr.tg_y (target y variable)
	Double_t L_tr_tg_ph[100]; //  L.tr.tg_ph (target phi)
	Double_t Lrb_dir_x;  //  Lrb.dir.x (x-component of beam direction vector)				
	Double_t Lrb_dir_z;  //  Lrb.dir.z (z-component of beam direction vector)				
	Double_t L_tr_n;     //  L.tr.n    (number of tracks)
	Double_t L_tr_vz[100];    //  L.tr.vz   (track reconstructed z-vertex)


	// Define branch addresses

	T->SetBranchAddress("L.tr.n",&L_tr_n);
	T->SetBranchAddress("L.tr.tg_y",L_tr_tg_y);
	T->SetBranchAddress("L.tr.tg_ph",L_tr_tg_ph);
	T->SetBranchAddress("Lrb.x",&Lrb_x);
	T->SetBranchAddress("Lrb.dir.x",&Lrb_dir_x);
	T->SetBranchAddress("Lrb.dir.z",&Lrb_dir_z);
	T->SetBranchAddress("L.tr.vz",L_tr_vz);

	TVector3 fPointingOffset;
	fPointingOffset.SetXYZ(-MissPointZ*TMath::Sin(HRSAngle)*TMath::Cos(HRSAngle),(Double_t)MissPointY,MissPointZ*TMath::Sin(HRSAngle)*TMath::Sin(HRSAngle));

	Double_t OpZ = 0;

	TVector3 Tg_YSpotTCS;
	TVector3 Tg_YSpotHCS;

	TVector3 MomDirectionTCS;
	TVector3 MomDirectionHCS;


	TRotation fTCSInHCS;
	TVector3 TCSX(0,-1,0);
	TVector3 TCSZ(TMath::Sin(HRSAngle),0,TMath::Cos(HRSAngle));
	TVector3 TCSY = TCSZ.Cross(TCSX);
	fTCSInHCS.RotateAxes(TCSX,TCSY,TCSZ);


		
	Int_t No_entries = T->GetEntries();


	for(Int_t i=0;i<No_entries;i++){

	  if(i%100000==0) cout << " events processed = " << i << endl;


	  
	  
	  T->GetEntry(i);


	  
	  Tg_YSpotTCS.SetXYZ(0, L_tr_tg_y[0],0);
	  Tg_YSpotHCS=fTCSInHCS*Tg_YSpotTCS+fPointingOffset;
	  
	  MomDirectionTCS.SetXYZ(0,L_tr_tg_ph[0],1);
	  MomDirectionHCS=fTCSInHCS*MomDirectionTCS;
	  
	  
	  OpZ = (Tg_YSpotHCS.X() -  Lrb_dir_x - (MomDirectionHCS.X()/MomDirectionHCS.Z())*Tg_YSpotHCS.Z() )/( Lrb_dir_x/ Lrb_dir_z - (MomDirectionHCS.X()/MomDirectionHCS.Z()) );
	  

	  hOpZ->Fill(OpZ);

   
	}

	TCanvas* canv_op = new TCanvas("canv_op","canv_op");

	hOpZ->Draw();

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	TString OTPT = 'root';

	gSystem->Exec("cp -vf " + CutFileName + " " + CutFileName + ".old");

	fstream cutdesc(CutDescFileName, ios_base::out);
	assert(cutdesc.is_open());

	TFile *f1 = new TFile(CutFileName, "UPDATE");
	assert(f1);
	// 	f->cd(); // goback to f directory

	// Define Canvas
	TCanvas* c1 = new TCanvas("Cuts", "Cuts", 900, 900);
	gStyle->SetPalette(1, 0);

	// Define Histogramsif (i==fmin)  c1->Print(Form("./cutfiles/gcut_L_%d.ps(",nrun));
	// 	else c1->Print(Form("./cutfiles/gcut_L_%d.ps",nrun)); // insert page into ps file
	// 	f1->Write();

	TH2F* h1 = new TH2F("h1", "ReactZ vs. Target Phi", 400, -0.05, 0.035, 400,-1, 1);
	//	TH2F* h1 = new TH2F("h1", "ReactZ vs. Target Phi", 900, -0.1, 0.06, 900,-0.5, 0.5);
	assert(h1);
	// 	TH2F* h2 = new TH2F ("h2","theta_target vs. phi_target", 900, -0.06,0.06,900,-0.07,0.07);

	// Draw ReactZ vs. Phi_rotate
	T->Draw("L.tr.vz:L.tr.tg_ph>>h1", GenrealCut, "COLZ"); // need finer delta cut later
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

	nfoil = nfoils; // Added to process seperate runs with one foil each)
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
			cutg->SetVarX("L.tr.tg_ph");
			cutg->SetVarY("L.tr.vz");

			// output cut to disk
			// 			f1->cd();
			cutg->Write("", TObject::kOverwrite); // Overwrite old cut

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


void cut_Vertex_new(int overwrite = 0, int nfoils = 3, int FoilID = -1) {
	cout << "open " << *SoureRootFile << endl;

	ReLoadcuts();

	// TFile* f = new TFile(SoureRootFile);
	// assert(f);
	// TTree * T = (TTree *) f->GetObjectChecked("T", "TTree");
	// assert(T);
	
	//	TChain* T = Load_rootfile(Run_number);
	TChain* T = Load_more_rootfiles(Run_number, Run_number_2);

	TString OTPT = 'root';

	gSystem->Exec("cp -vf " + CutFileName + " " + CutFileName + ".old");

	fstream cutdesc(CutDescFileName, ios_base::out);
	assert(cutdesc.is_open());

	TFile *f1 = new TFile(CutFileName, "UPDATE");
	assert(f1);
	// 	f->cd(); // goback to f directory

	// Define Canvas
	TCanvas* c1 = new TCanvas("Cuts", "Cuts", 900, 900);
	gStyle->SetPalette(1, 0);

	// Define Histogramsif (i==fmin)  c1->Print(Form("./cutfiles/gcut_L_%d.ps(",nrun));
	// 	else c1->Print(Form("./cutfiles/gcut_L_%d.ps",nrun)); // insert page into ps file
	// 	f1->Write();

	TH2F* h1 = new TH2F("h1", "ReactZ vs. Target Phi", 600, -0.05, 0.035, 600,-0.2, 0.2);
	//	TH2F* h1 = new TH2F("h1", "ReactZ vs. Target Phi", 900, -0.1, 0.06, 900,-0.5, 0.5);
	assert(h1);
	// 	TH2F* h2 = new TH2F ("h2","theta_target vs. phi_target", 900, -0.06,0.06,900,-0.07,0.07);

	// Draw ReactZ vs. Phi_rotate
	T->Draw("L.tr.vz:L.tr.tg_ph>>h1", GenrealCut, "COLZ"); // need finer delta cut later
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

	nfoil = nfoils; // Added to process seperate runs with one foil each)
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
			cutg->SetVarX("L.tr.tg_ph");
			cutg->SetVarY("L.tr.vz");

			// output cut to disk
			// 			f1->cd();
			cutg->Write("", TObject::kOverwrite); // Overwrite old cut

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




void CutColumn(int FoilID = 0 /* FoilID*/, int colNo = 0 /*Number of columns*/, int firstCol = 0, int overwrite = 0, bool append = true) {

  /* Function designed to cut on Columns of sieve based on plot of Focal Plane y against phi
     Idea from Vardan
  
     yyp_fp->Fill(R_tr_ph_fp[0], R_tr_y_fp[0]);
     t->SetBranchAddress("R.tr.ph",R_tr_ph_fp);
     t->SetBranchAddress("R.tr.y",R_tr_y_fp);

  */

  ReLoadcuts();

  std::cout << "CutFileName = " << CutFileName << std::endl;


  TCanvas * c1 = new TCanvas(Form("PlotColSieve%d", FoilID), "PlotSieve (Columns)",1000, 1000);
  c1->Update();
  
  // TFile* f = new TFile(SoureRootFile);
  // TTree * T = (TTree *) f->GetObjectChecked("T", "TTree");
  TChain* T = Load_more_rootfiles(Run_number, Run_number_2);

  
  //	gSystem->Exec("cp -vf " + SoureRootFile + " " + CutFileName + ".old");
  TFile *f1 = new TFile(CutFileName, "UPDATE");
  assert(f1);
  (TCutG*) f1->GetObjectChecked(Form("fcut_L_%d", FoilID), "TCutG"); //looking for foil cut definition
  

  fstream cutdesc;

  if(append){
    cutdesc.open(CutDescFileNameCol, ios_base::out |ios::app);
  }
  if(!append){
    cutdesc.open(CutDescFileNameCol, ios_base::out);
  }
  assert(cutdesc.is_open());

  // plot of y against phi (focal plane)
  //  TH2F* h2 = new TH2F("h2", "y vs. phi (focal plane ", 1000, -0.03, 0.04, 900, -0.04, 0.035);
  TH2F* h2 = new TH2F("h2", "y vs. phi (focal plane) ", 1000, -0.04, 0.04, 900, -0.04, 0.04);
  //  T->Draw("L.tr.y:L.tr.ph>>h2","", "COLZ");
  T->Draw("L.tr.y:L.tr.ph>>h2", GenrealCut + TCut(Form("fcut_L_%d", FoilID)), "COLZ");
  c1->Update();
  
  //	int ncol = 1, colID=1;
  
  // make cuts from the first available column
  for (int ColID = firstCol; ColID < colNo; ColID++) //column number ColID
    {
      
      
      TCutG* tmpcut = (TCutG*) gROOT->FindObject(Form("ccut_L_%d_%d",FoilID,ColID));
      //WARNING: Overwriting
      //			tmpcut = 0;
      //			cout<<"WARNING: Overwiting!"<<endl;
      //WARNING: Overwriting
      if (tmpcut && !overwrite) {
	
	cout << Form("ccut_L_%d_cols", FoilID)
	     << " is found, using old one" << endl;
      } else {
	//				TPad * gPad = TVirtualPad::Pad();
	
	// define cuts
	cout << "making cut for foil=" << FoilID << ", column=" << ColID
	     << " ... waiting ..." << endl;
	TCutG* cutg = (TCutG*) gPad->WaitPrimitive("CUTG", "CutG");
	c1->Update();
	cout << "done!" << endl;
	cutg->SetName(Form("ccut_L_%d_%d", FoilID,ColID));
	  
	// set axises
	cutg->SetVarX("L.tr.ph[0]");
	cutg->SetVarY("L.tr.y[0]");
	
	// output
	f1->cd();
	cutg->Write("", TObject::kOverwrite);
	tmpcut = cutg;
	//			f->cd();
      }
      
      tmpcut = (TCutG*) f1->GetObjectChecked(Form("ccut_L_%d_%d", FoilID,ColID), "TCutG");
      assert(tmpcut);
      tmpcut->SetLineWidth(2);
      tmpcut->SetLineColor(kMagenta);
      tmpcut->Draw("PL");
      c1->Update();


      cout << "log to " << CutDescFileNameCol << endl;
      cutdesc << Form("ccut_L_%d_%d", FoilID, ColID) << " && " << (const char *) GenrealCut
	      << endl;


    }

  
  cutdesc.close();
  f1->Write();
  SaveCanvas(c1, *SoureRootFile + Form("ccut_L_%d_cols", FoilID), kFALSE);
  //	SaveCanvas(c1, SoureRootFile + Form("hcut_L_%d_%d_%d", i, colID), kFALSE);
  //	f1->Close();
  f1->ls();
  
}




void test_APEX(){

  
  TVector3 test_vect = GetSieveHoleTCS(4,4);
  cout << "Vector X,Y,Z = " << test_vect.X() << "," << test_vect.Y() << "," << test_vect.Z() << endl;


  std::vector<int> x_y;

  x_y = Get_Col_Row(88);
  
  cout << "Col = " << x_y[0] << endl;


  Int_t hole_no = Get_Hole(x_y[0],x_y[1]);

  cout << "hole_no = " << hole_no << endl;

  TVector3 test_vect2 = GetSieveHoleTCS(hole_no);

  cout << "Vector X,Y,Z = " << test_vect2.X() << "," << test_vect2.Y() << "," << test_vect2.Z() << endl;  


}


void CutSieve_ellipse(int FoilID = 0, int col = 0, int overwrite = 0, int append = 1) {

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




  // file that stores TCutGs (from foil and sieve)
  gSystem->Exec("cp -vf " + CutFileName + " " + CutFileName + ".old");
  TFile *f1 = new TFile(CutFileName, "UPDATE");
  assert(f1);
  (TCutG*) f1->GetObjectChecked(Form("fcut_L_%d", FoilID), "TCutG"); //looking for old cut definition



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
  cout<<fixed<<setprecision(1);

  TDatime* date = new TDatime();
  cutcsv<<date->GetDay()<<"/"<<date->GetMonth()<<"/"<<date->GetYear()<<" (dd/mm/yyyy)"<<endl;
  //  cutcsv << "Hole ID (col:row),Hole Exists,Included in opt,Ellipse ph cen,Expected ph,Ellipse th cen,Expected th,Semi axis ph,Semi axis th,Ellipse Tilt (deg),Ellipse ph rms,Ellipse th rms,Statistics,All angles in mrad except ellipse tilt which is positive counterclockwise from vertical axis"<<endl;
  cutcsv << "Hole ID: col, row, Included in opt, Ellipse ph cen, Expected ph, Ellipse th cen, Expected th, Semi axis ph, Semi axis th, Ellipse Tilt (deg), -  All angles in mrad except ellipse tilt which is positive counterclockwise from vertical axis"<<endl;
  
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

    
    T->Draw("L.tr.tg_th:L.tr.tg_ph>>h2", GenrealCut
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



  // Designed to copy over old file contents until hole of interest is reached

  // Int_t minhole = Get_Hole(col,rmin);

  // cout << "min hole id is " << minhole << endl;

  // string line;
  // getline(cutcsvold,line); // these lines are for first two line of csv file (explenatory)
  // getline(cutcsvold,line);
  // for(int i = 0; i<minhole; i++){
  //   getline(cutcsvold,line);
  //   cutcsv << line << endl;
  //  }
   

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

    while(cur_hole < Get_Hole(col,k)){
      //      cout << " Current hole = " << cur_hole << ", Get_hole = " << Get_Hole(col,k) << endl;
      getline(cutcsvold,line);	 
      //      cout << " current line (cvs) is " << endl;
      cout << line << endl << endl;
      
      cur_hole++;
      cutcsv<<line<<endl;
      
    }

	
  
  
    // ensure that old file is cycled through
    getline(cutcsvold,line);	 

  
    // try to find existing hole cut
    // use name e_hcut for ellipse_holecut
    cout << "Testing " << Form("e_hcut_L_%d_%d_%d", FoilID, col, k) << endl;
  
    TCutG* cutg = (TCutG*)gROOT->FindObject(Form("e_hcut_L_%d_%d_%d", FoilID, col, k));

  
    // TEllipse* Ellipse = new TEllipse(0,0,0.0021,0.0043,0,360,0);
    // Ellipse->SetFillStyle(0);
    // Ellipse->SetLineWidth(3);
    // Ellipse->Draw("same");
    
    // if new ellipse is being drawn
    if(!cutg || overwrite){
      if(cutg)cutg->Delete();

      //    TEllipse* Ellipse = new TEllipse(0,0,0.0006,0.0014,0,360,0);
      TEllipse* Ellipse = new TEllipse(-0.02,0,0.0006,0.0014,0,360,0);
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
	    
      cutg->SetVarX("L.tr.tg_ph");
      cutg->SetVarY("L.tr.tg_th");
      cout << "done!" << endl;
	    
    }
    else {
      cout << Form("e_hcut_L_%d_%d_%d", FoilID, col, k) << " is found, using old one" << endl;
      
      cout << "Using line : " << endl;
      cout << line << endl << endl;
      cutcsv<<line<<endl;

    }


    // draw TCutG (new or old if not overwritten)
    cutg->SetLineColor(kMagenta);	
    cutg->SetLineWidth(2);	
    cutg->Draw("PL");
    c1->Update();

    cutg->Write("", TObject::kOverwrite); // Overwrite old cut



    



    // John W: attempt to add seconday cut on FP theta vs y. This acts on a secondary check on accidentaly including events wrongly identified from a different hole.


    // c2->cd(2);

    // T->Draw("L.tr.r_th*1000:L.tr.r_y>>thfp_v_yfp",GenrealCut + TCut(Form("fcut_L_%d", FoilID)),"colz");


    c2->cd();

    T->Draw("L.tr.r_th*1000:L.tr.r_y>>thfp_v_yfp",GenrealCut + TCut(Form("fcut_L_%d", FoilID)) + TCut(Form("e_hcut_L_%d_%d_%d", FoilID, col, k)),"colz");

    c2->Update();


    //    TVirtualPad *cpad = gPad;
    
    TCutG* tmpcut_FP = (TCutG*) gROOT->FindObject(Form("FPhcut_L_%d_%d_%d",FoilID, col, k));

    
    if (tmpcut_FP && !overwrite) {
      
      cout << Form("FPhcut_L_%d_%d_%d",FoilID, col, k)
	   << " is found, using old one" << endl;
      tmpcut_FP->Draw("PL");

      tmpcut_FP->SetVarX("L.tr.r_y[0]");
      tmpcut_FP->SetVarY("L.tr.r_th[0]*1000");
      tmpcut_FP->Write("",TObject::kOverwrite);
      

    }
    else{

      cout << "Looking for FP hole cut... " << endl;

      c2->cd();

      TCutG* cutg_FP = (TCutG*) (gPad->WaitPrimitive("CUTG", "CutG"));

      c2->Update();

      cutg_FP->SetName(Form("FPhcut_L_%d_%d_%d",FoilID, col, k));

      
      cutg_FP->SetVarX("L.tr.r_y[0]");
      cutg_FP->SetVarY("L.tr.r_th[0]*1000");

      f1->cd();
      cutg_FP->Write("",TObject::kOverwrite);
      tmpcut_FP = cutg_FP;
      

    }



    tmpcut_FP =  (TCutG*) f1->GetObjectChecked(Form("FPhcut_L_%d_%d_%d",FoilID, col, k), "TCutG");
    assert(tmpcut_FP);
    tmpcut_FP->SetLineWidth(2);
    tmpcut_FP->SetLineColor(kMagenta);
    tmpcut_FP->Draw("PL");
    c2->Update();
 
    
    




    // Log cuts names to text file   
    cutdesc << Form("FPhcut_L_%d_%d_%d",FoilID, col, k) << " && " << Form("e_hcut_L_%d_%d_%d", FoilID, col, k) << " && " << Form("fcut_L_%d", FoilID) << " && " << (const char *) GenrealCut << endl;
    //    cutdesc << Form("e_hcut_L_%d_%d_%d", FoilID, col, k) << " && " << Form("fcut_L_%d", FoilID) << " && " << (const char *) GenrealCut << endl;
    
  }


 
  // puts remainder of old csv file into new
  while(!cutcsvold.eof()){
    getline(cutcsvold,line);
    cutcsv<<line<<endl;
  }


       
  f1->Write();
  f1->ls();
  cutcsv.close();
  cutdesc.close();
  

}





void CutSieve_col(int i = 0 /*Foil ID*/, int cmin = 0 /*Column ID*/, int overwrite = 0, bool append = true) {

  int ncol = 1;

  ReLoadcuts();

  //	TCanvas * c1 = new TCanvas(Form("PlotSieve%d_%d", i, cmin), "PlotSieve",1800, 1150);
  TCanvas * c1 = new TCanvas(Form("PlotSieve_%d_%d", i, cmin), "PlotSieve",1000, 1000);
  c1->Update();
  
  // TFile* f = new TFile(SoureRootFile);
  // TTree * T = (TTree *) f->GetObjectChecked("T", "TTree");
  TChain* T = Load_more_rootfiles(Run_number, Run_number_2);

  //	gSystem->Exec("cp -vf " + SoureRootFile + " " + CutFileName + ".old");
  TFile *f1 = new TFile(CutFileName, "UPDATE");
  assert(f1);

  //  TString test_string = Form("ccut_L_%d_cols",i,cmin);
  //  std::cout << "test_string is " << test_string << std::endl;
  (TCutG*) f1->GetObjectChecked(Form("ccut_L_%d_%d", i, cmin), "TCutG"); //looking for column cut definition
  (TCutG*) f1->GetObjectChecked(Form("fcut_L_%d", i), "TCutG"); //looking for foil cut definition



  fstream cutdesc;


  if(append){
    cutdesc.open(CutDescFileNameSieve, ios_base::out |ios::app);
  }
  if(!append){
    cutdesc.open(CutDescFileNameSieve, ios_base::out);
  }
  assert(cutdesc.is_open());


  // make cuts for each holes
  TH2F* h2 = new TH2F("h2", "theta_target vs. phi_target", 600, -0.04, 0.01, 600, -0.04, 0.04);
  //  TH2F* h2 = new TH2F("h2", "theta_target vs. phi_target", 900, -0.04, 0.04, 900, -0.08, 0.08);
  T->Draw("L.tr.tg_th:L.tr.tg_ph>>h2", GenrealCut
	  + TCut(Form("fcut_L_%d", i)) + TCut(Form("ccut_L_%d_%d", i, cmin)) , "COLZ");
  c1->Update();

  //	int ncol = 1, cmin=1;
  
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
	int k = rmin + cnt * 2 ;
	// JW changed cnt * 2 -> cnt, as I wasn't sure of purpose of *2
	//			cout<<"Hole ID : "<<k<<endl;
	//			cin>>k;
	if (k < 0)
	  continue;
	
	TCutG* tmpcut = (TCutG*) gROOT->FindObject(Form("hcut_L_%d_%d_%d",
							i, j, k));
	//WARNING: Overwriting
	//			tmpcut = 0;
	//			cout<<"WARNING: Overwiting!"<<endl;
	//WARNING: Overwriting
	if (tmpcut && !overwrite) {
	  
	  cout << Form("hcut_L_%d_%d_%d", i, j, k)
	       << " is found, using old one" << endl;
	} else {
	  //				TPad * gPad = TVirtualPad::Pad();
	  
	  // define cuts
	  cout << "making cut for foil=" << i << ", column=" << j
	       << ", No." << k << " hole: waiting ..." << endl;
	  TCutG* cutg = (TCutG*) gPad->WaitPrimitive("CUTG", "CutG");
	  c1->Update();
	  cout << "done!" << endl;
	  cutg->SetName(Form("hcut_L_%d_%d_%d", i, j, k));
	  
	  // set axises
	  cutg->SetVarX("L.tr.tg_ph[0]");
	  cutg->SetVarY("L.tr.tg_th[0]");
	  
	  // output
	  f1->cd();
	  cutg->Write("", TObject::kOverwrite);
	  tmpcut = cutg;
	  //			f->cd();
	}
	
	tmpcut = (TCutG*) f1->GetObjectChecked(Form("hcut_L_%d_%d_%d", i, j, k), "TCutG");
	assert(tmpcut);
	tmpcut->SetLineWidth(2);
	tmpcut->SetLineColor(kMagenta);
	tmpcut->Draw("PL");
	c1->Update();


	cout << "log to " << CutDescFileNameSieve << endl;
	cutdesc << Form("hcut_L_%d_%d_%d", i, j, k) << " && " << (const char *) GenrealCut
		<< endl;
      }
    }


  cutdesc.close();
  f1->Write();
  SaveCanvas(c1, *SoureRootFile + Form("hcut_L_%d_%d", i, cmin), kFALSE);
  //	SaveCanvas(c1, SoureRootFile + Form("hcut_L_%d_%d_%d", i, cmin), kFALSE);
  //	f1->Close();
  f1->ls();
  
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
	  cutg->Write("", TObject::kOverwrite);
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

		cutg->Write("",TObject::kOverwrite); // Overwrite old cut

		SaveCanvas(c1, *SoureRootFile + Form("L_dp_%d",FoilID), kFALSE);
	}

	fcut->Write();

//	cutdesc.close();
//	cutdescSieve.close();
}

void cut_L(/*int nrun=0*/) {

	// 	PlotRun();
	// 	cut_Vertex();
	CutSieve(0, 9);

}


void cut_VertexAdv(int overwrite = 0, TCut extracuts = "1") {
	cout << "open " << *SoureRootFile << endl;

	TFile* f = new TFile(*SoureRootFile);
	assert(f);
	TTree * T = (TTree *) f->GetObjectChecked("T", "TTree");
	assert(T);

	TString OTPT = 'root';

	gSystem->Exec("cp -vf " + CutFileName + " " + CutFileName + ".old");

	fstream cutdesc(CutDescFileName, ios_base::out);
	assert(cutdesc.is_open());

	TFile *f1 = new TFile(CutFileName, "UPDATE");
	assert(f1);
	// 	f->cd(); // goback to f directory

	// Define Canvas
	TCanvas* c1 = new TCanvas("Cuts", "Cuts", 900, 900);
	gStyle->SetPalette(1, 0);

	// Define Histogramsif (i==fmin)  c1->Print(Form("./cutfiles/gcut_L_%d.ps(",nrun));
	// 	else c1->Print(Form("./cutfiles/gcut_L_%d.ps",nrun)); // insert page into ps file
	// 	f1->Write();
	TH2F* h1 = new TH2F("h1", "ReactZ vs. Target Phi", 900, -0.06, 0.06, 900,
			-0.5, 0.5);
	assert(h1);
	// 	TH2F* h2 = new TH2F ("h2","theta_target vs. phi_target", 900, -0.06,0.06,900,-0.07,0.07);

	// Draw ReactZ vs. Phi_rotate
	c1->Update();

	// Set y target number

	int fmin = 0;

	// output ps files
	//c1->Print(Form("./cutfiles/gcut_L_%d.ps[",nrun));

		T->Draw("L.tr.vz:L.tr.tg_ph>>h1", GenrealCut, "COLZ"); // need finer delta cut later
//		T->Draw("ReactPt_L_urb.z:L.tr.tg_ph>>h1", GenrealCut, "COLZ"); // need finer delta cut later
	// Choose the foil you want to make cut)
	for (int i = fmin; i < nfoil + fmin; i++) {
		cout << "new foil " << i << endl;
		TCutG* cutg, *tmpcut;
		// 		if (i==fmin)  c1->Print(Form("./cutfiles/gcut_L_%d.ps(",nrun));
		// 		else c1->Print(Form("./cutfiles/gcut_L_%d.ps",nrun)); // insert page into ps file
		// 		f1->Write();
		c1->Update();
		TVirtualPad *cpad = gPad;

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
			cutg->SetVarX("L.tr.tg_ph");
			cutg->SetVarY("L.tr.vz");

			// output cut to disk
			// 			f1->cd();
			cutg->Write("", TObject::kOverwrite); // Overwrite old cut

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



//_____________________________________________________________________________


// std::vector<int> Get_Col_Row(Int_t Hole){


//   Int_t row_comp = 0;
//   Int_t no_col = 0;
//   Int_t col = 0;
//   Int_t row = 0;


//   for(Int_t i = 0; i<NSieveRow; i++){

//     row_comp += NoinEachRow[i];

//     if( (row_comp-1) >= Hole){
//       row = i;
//       no_col = Hole - ( row_comp - NoinEachRow[i]);
      
//       break;
//     }

//   }



//   if(row%2 == 0){
//     if(no_col==13){
//       col = 25;
//     }
//     else if (no_col==14){
//       col = 26;
//     }
//     else{

//     col = no_col *2;
//     }
//   }

//   if(row%2 == 1){

//     if(row > 1 && row < 15){
//       if(no_col >5){
// 	col = (no_col*2)+3;
//       }
//       else{
// 	col = (no_col*2) +1;
//       }
//     }
//     else{
//       col = (no_col*2) +1;
//     }

    
//   }


//   //  cout << "For Hole " << Hole << ": row = " << row << " and col = " << col << endl;
//     //  hole_no += Col;
//   std::vector<int> rowcol{col, row};

//   return rowcol;
// }




// //____________________________________________________________________________



// Int_t Get_Hole(Int_t Col, Int_t Row){


//   Int_t hole_no = 0;

//   for(Int_t i = 0; i<Row; i++){

//     hole_no += NoinEachRow[i];

//   }

//   // else if conditions here deal with columns at right edge of sieve slit (area odd rows do not have holes)

//   if(Row%2 == 0){
//     if(Col < 25){
//       hole_no += (Col/2);
//     }
//     else if(Col == 25){
//       hole_no += 13;      
//     }
//     else if(Col == 26){
//       hole_no += 14;
//     }      
//   }


//   if(Row%2 == 1){
//     if(Row > 1 && Row < 15){
//       if(Col < 13){
// 	Col = ((Col+1)/2) - 1;
//       }
//       else if (Col >= 13){
// 	Col = ((Col-3)/2);
//       }      
//     }
//     else{
//       Col = ((Col+1)/2) -1;
//     }    
//     hole_no += Col;
//   }
  

//   //  hole_no += Col;



//   return hole_no;
// }


//____________________________________________________________________________




