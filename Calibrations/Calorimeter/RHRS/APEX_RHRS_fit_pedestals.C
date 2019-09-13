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

#include "Load_more_rootfiles.C"

// Need these up here so they can be accessed by both
// the main function and the function to be minimized

Int_t nBlkPreShower;
Int_t nBlkShower;
char arm[2]; 
//TChain *T = new TChain("T");
vector<Double_t> ped_preShow;		// preshower pedestals
vector<Double_t> ped_show;		// shower pedestals

void APEX_RHRS_fit_pedestals(Int_t run = 4648) { 

	//RHRS
	sprintf(arm, "R");
	nBlkPreShower = 48;
	nBlkShower = 75;	

	// add root files to chain
	//	T->Add("./rootfiles/tritium_online_*");

	// fit the pedestals

	TChain* T = Load_more_rootfiles(run);


	Int_t numEvents = T->GetEntries();

	cout << " numevents = " << numEvents << endl;	

	TCanvas *c1=new TCanvas("c1","For fitting the pedestals",800,600);
	c1->SetLogy();
	c1->SetTickx();
	c1->SetGridx();

	TString ped_name1=Form("%sHRS_PS-SH_Ped.pdf[",arm);
	TString ped_name2=Form("%sHRS_PS-SH_Ped.pdf",arm);
	TString ped_name3=Form("%sHRS_PS-SH_Ped.pdf]",arm);

	Double_t preShowPed;
	Double_t showPed;
	
	for (Int_t i=0;i<nBlkPreShower;i++){

		TH1F *preShow_ped= new TH1F("preShow_ped",Form("Pedstal for preshower cal PMT %d",i),1800,200,2000);
		
		//Preshower
	
//		c1->cd(1);
		T->Draw(Form("%s.ps.a[%d]>>preShow_ped",arm,i));

		Int_t binmax=preShow_ped->GetMaximumBin();

		cout << " binmax = " << binmax << endl;

		Double_t x1=preShow_ped->GetXaxis()->GetBinCenter(binmax);

		TH1D *htmp = (TH1D*)preShow_ped->Clone("htmp");
		htmp->GetXaxis()->SetRange(binmax-5, binmax+5);
		//		Double_t ped_fit_range = 2.5*htmp->GetRMS();
		Double_t ped_fit_range = htmp->GetRMS();
		//		Double_t ped_fit_range = 2.5*htmp->GetRMS();

		preShow_ped->Fit("gaus", "R", "", x1-ped_fit_range, x1+ped_fit_range);
		TF1* gaus1 = preShow_ped->GetFunction("gaus");
		preShowPed = gaus1->GetParameter(1);
		cout << "preShowped = " << preShowPed << ", preshowwidth = " << gaus1->GetParameter(2) << endl;
		preShow_ped->Fit("gaus", "R", "", gaus1->GetParameter(1)-2.*gaus1->GetParameter(2), gaus1->GetParameter(1)+2.*gaus1->GetParameter(2));
		delete htmp;

		ped_preShow.push_back(preShowPed);

		gPad->SetTickx();
		gPad->SetLogy();
		gPad->SetGridx();

		if(i == 0) {
			c1->SaveAs(ped_name1);
		} else {
			c1->SaveAs(ped_name2);
		}
		
		
		delete preShow_ped;
	}

	for(int i = 0; i<nBlkShower; i++) { 

		//Shower
		TH1F *show_ped= new TH1F("show_ped",Form("Pedstal for shower cal PMT %d",i),600,200,800);
		//		T->Draw(Form("%s.sh.a[%d]>>show_ped",arm,i),Form("D%s.evtypebits>>4&1",arm));
		T->Draw(Form("%s.sh.a[%d]>>show_ped",arm,i));
	
		Int_t binmax2=show_ped->GetMaximumBin();

		cout << " binmax2 = " << binmax2 << endl;

		Double_t x2=show_ped->GetXaxis()->GetBinCenter(binmax2);

		TH1D *htmp2 = (TH1D*)show_ped->Clone("htmp2");
		htmp2->GetXaxis()->SetRange(binmax2-5, binmax2+5);
		Double_t ped_fit_range2 = htmp2->GetRMS();

		show_ped->Fit("gaus", "R", "", x2-ped_fit_range2, x2+ped_fit_range2);
		TF1* gaus2 = show_ped->GetFunction("gaus");
		showPed = gaus2->GetParameter(1);	
		show_ped->Fit("gaus", "R", "", gaus2->GetParameter(1)-2.*gaus2->GetParameter(2), gaus2->GetParameter(1)+2.*gaus2->GetParameter(2));
		delete htmp2;
		
		ped_show.push_back(showPed);

		gPad->SetTickx();
		gPad->SetLogy();
		gPad->SetGridx();

		if(i == nBlkShower-1) {
			c1->SaveAs(ped_name3);
		} else {
			c1->SaveAs(ped_name2);
		}
	
		delete show_ped;
	}

	ofstream pedPSfile("output_data/ped_preshower.dat");
	ofstream pedSHfile("output_data/ped_shower.dat");
	
	for(int i = 0; i < nBlkPreShower; i++) {
		pedPSfile << ped_preShow[i] << "\n";
	}
	for(int i = 0; i < nBlkShower; i++) {
		pedSHfile << ped_show[i] << "\n";
	}
	pedPSfile.close();
	pedSHfile.close();
}


