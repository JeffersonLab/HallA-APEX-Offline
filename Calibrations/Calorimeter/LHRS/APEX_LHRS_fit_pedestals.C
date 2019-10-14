#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include  <stdio.h>
#include  <stdlib.h>
#include <vector>
#include "TMinuit.h"
#include "TFitter.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"

// Need these up here so they can be accessed by both
// the main function and the function to be minimized

Int_t numEvents;
Int_t nBlkPRL1;
Int_t nBlkPRL2;
char arm[2]; 
//TChain *T = new TChain("T");
vector<Double_t> ped_prl1;		// prl1er pedestals
vector<Double_t> ped_prl2;		// shower pedestals

#include "Load_more_rootfiles.C"

void APEX_LHRS_fit_pedestals(Int_t run = 4179) { 

	//LHRS
	sprintf(arm, "L");
	nBlkPRL1 = 34;
	nBlkPRL2 = 34;	

	TChain* T = Load_more_rootfiles(run);

  //	T->Add("./rootfiles/apex_online_4240.root");

	numEvents = T->GetEntries();

	cout << " numevents = " << numEvents << endl;

	// fit the pedestals

	TCanvas *c1=new TCanvas("c1","",800,600);
	c1->SetLogy();
	c1->SetTickx();

	TString ped_name1=Form("%sHRS_PRL_Ped.pdf[",arm);
	TString ped_name2=Form("%sHRS_PRL_Ped.pdf",arm);
	TString ped_name3=Form("%sHRS_PRL_Ped.pdf]",arm);

	Double_t prl1Ped;
	Double_t prl2Ped;
	
	for (Int_t i=0;i<nBlkPRL1;i++){

	  cout << " ~~~~~~~~~~~~~~ \n" << " iteration " << i << "\n ~~~~~~~~~~~" << endl;

		TH1F *prl1_ped= new TH1F("prl1_ped",Form("Pedstal for preshower cal PMT %d",i),3000,-500,15000);
		
		//Preshower
	
		//		T->Draw(Form("%s.prl1.a[%d]>>prl1_ped",arm,i),Form("D%s.evtypebits>>2&1",arm)); 

		//tried removing trigger requirement

		T->Draw(Form("%s.prl1.a[%d]>>prl1_ped",arm,i)); 

		Int_t binmax=prl1_ped->GetMaximumBin();

		cout << " binmax = " << binmax << endl;

		Double_t x1=prl1_ped->GetXaxis()->GetBinCenter(binmax);

		TH1D *htmp = (TH1D*)prl1_ped->Clone("htmp");
		htmp->GetXaxis()->SetRange(binmax-50, binmax+50);// changed 20 to 50
		Double_t ped_fit_range = 2.5*htmp->GetRMS(); 

		cout << " ped_fit_range " << ped_fit_range << endl; 

		prl1_ped->Fit("gaus", "R", "", x1-ped_fit_range, x1+ped_fit_range);
		TF1* gaus1 = prl1_ped->GetFunction("gaus");
		prl1Ped = gaus1->GetParameter(1);
		prl1_ped->Fit("gaus", "R", "", gaus1->GetParameter(1)-2.*gaus1->GetParameter(2), gaus1->GetParameter(1)+2.*gaus1->GetParameter(2));
	       	TF1* gaus1_2 = prl1_ped->GetFunction("gaus");
		prl1Ped = gaus1_2->GetParameter(1);

		delete htmp;

		ped_prl1.push_back(prl1Ped);
//		ped_prl1.push_back(300.);

		gPad->SetTickx();
		gPad->SetLogy();
		gPad->SetGridx();
		if(i == 0) {
			c1->SaveAs(ped_name1);
		} else {
			c1->SaveAs(ped_name2);
		}
		
		delete prl1_ped;
	}

	for(int i = 0; i<nBlkPRL2; i++) { 

	  cout << "+++++++++++++++" << "\n iteration " << i << "\n ++++++++++++++ " << endl;
	  
		TH1F *prl2_ped= new TH1F("prl2_ped",Form("Pedstal for shower cal PMT %d",i),600,-500,15000);
		
		//Shower
	
		//		T->Draw(Form("%s.prl2.a[%d]>>prl2_ped",arm,i),Form("D%s.evtypebits>>2&1",arm));
		T->Draw(Form("%s.prl2.a[%d]>>prl2_ped",arm,i));
	
		Int_t binmax2=prl2_ped->GetMaximumBin();
		cout << " binmax2 = " << binmax2 << endl;

		Double_t x2=prl2_ped->GetXaxis()->GetBinCenter(binmax2);

		TH1D *htmp2 = (TH1D*)prl2_ped->Clone("htmp2");
		htmp2->GetXaxis()->SetRange(binmax2-20, binmax2+20);
		Double_t ped_fit_range2 = 2.5*htmp2->GetRMS();

		cout << " ped_fit_range2 " << ped_fit_range2 << endl; 


		prl2_ped->Fit("gaus", "R", "", x2-ped_fit_range2, x2+ped_fit_range2);
		TF1* gaus2 = prl2_ped->GetFunction("gaus");
		prl2Ped = gaus2->GetParameter(1);	
		prl2_ped->Fit("gaus", "R", "", gaus2->GetParameter(1)-2.*gaus2->GetParameter(2), gaus2->GetParameter(1)+2.*gaus2->GetParameter(2));
		TF1* gaus2_2 = prl2_ped->GetFunction("gaus");
		prl2Ped = gaus2_2->GetParameter(1);

		delete htmp2;
		
		ped_prl2.push_back(prl2Ped);

		gPad->SetTickx();
		gPad->SetLogy();
		gPad->SetGridx();
		if(i == nBlkPRL2-1) {
			c1->SaveAs(ped_name3);
		} else {
			c1->SaveAs(ped_name2);
		}
		delete prl2_ped;
	}

	ofstream ped1file(Form("DB_output/%d_ped_prl1.dat",run));
	ofstream ped2file(Form("DB_output/%d_ped_prl2.dat",run));
	
	for(int i = 0; i < nBlkPRL1; i++) {
		ped1file << ped_prl1[i] << "\n";
	}
	for(int i = 0; i < nBlkPRL2; i++) {
		ped2file << ped_prl2[i] << "\n";
	}
	ped1file.close();
	ped2file.close();
}



