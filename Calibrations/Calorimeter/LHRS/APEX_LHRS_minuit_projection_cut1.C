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
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "Minuit2/Minuit2Minimizer.h"
#include <TMinuit.h>


#include "Load_more_rootfiles.C"


// Need these up here so they can be accessed by both
// the main function and the function to be minimized

Int_t numEvents;
Int_t numBlockPRL1;
Int_t numBlockPRL2;
Int_t numRowPRL1;
Int_t numColPRL1;
Int_t numRowPRL2;
Int_t numColPRL2;
char arm[2]; 
TChain *T = new TChain("T");
vector<Double_t> pedestals_PRL1;		// preshower pedestals
vector<Double_t> pedestals_PRL2;		// shower pedestals
Int_t numGoodEvents = 0;
vector<pair<pair<int, double>, vector<pair<int, double> > > > goodEvents;

// Block positions for projection
Double_t prl1_xmin[17] = {1.107,0.959,0.812,0.664,0.517,0.369,0.222,0.074,-0.073,-0.221,-0.369,-0.516,-0.664,-0.811,-0.959,-1.106,-1.254};
Double_t prl1_xmax[17] = {1.254,1.107,0.959,0.812,0.664,0.517,0.369,0.222,0.074,-0.073,-0.221,-0.369,-0.516,-0.664,-0.811,-0.959,-1.106};
Double_t prl2_xmin[17] = {1.125,0.975,0.825,0.675,0.525,0.375,0.225,0.075,-0.075,-0.225,-0.375,-0.525,-0.675,-0.825,-0.975,-1.125,-1.275};
Double_t prl2_xmax[17] = {1.275,1.125,0.975,0.825,0.675,0.525,0.375,0.225,0.075,-0.075,-0.225,-0.375,-0.525,-0.675,-0.825,-0.975,-1.125};

// Find PRL1 row from projected track
Int_t prl1_x_to_row(Double_t x) {
	Int_t row = -99;
	for(Int_t i = 0; i < numRowPRL1; i++) {
		if(x > prl1_xmin[i] && x < prl1_xmax[i]) {
			row = i;
			break;
		}
	}
	return row;
}

// Find PRL2 row from projected track
Int_t prl2_x_to_row(Double_t x) {
	Int_t row = -99;
	for(Int_t i = 0; i < numRowPRL2; i++) {
		if(x > prl2_xmin[i] && x < prl2_xmax[i]) {
			row = i;
			break;
		}
	}
	return row;
}

void getGoodEvents() {
	
	Int_t edge[] = {0, 16, 17, 33};
	Int_t Nedge = sizeof(edge)/sizeof(Int_t);

	Double_t scer[100];
	Double_t Eps, Esh;
	Double_t n_tr;
	Double_t n_u1_clust;
	Double_t n_u2_clust;
	Double_t n_v1_clust;
	Double_t n_v2_clust;
	Double_t p[100];
	Double_t prl1_adc[100];
	Double_t prl2_adc[100];	
	Double_t prl1_trx[100];
	Double_t prl1_try[100];
	Double_t prl2_trx[100];
	Double_t prl2_try[100];
	Double_t evtype;
	
	T->SetBranchStatus("*",0);  

	T->SetBranchStatus(Form("%s.tr.n",arm),1);
	T->SetBranchStatus(Form("%s.cer.asum_c",arm),1);
	T->SetBranchStatus(Form("%s.prl1.e",arm),1);
	T->SetBranchStatus(Form("%s.prl2.e",arm),1);
	T->SetBranchStatus(Form("%s.prl1.a",arm),1);
	T->SetBranchStatus(Form("%s.prl2.a",arm),1);
	T->SetBranchStatus(Form("%s.tr.p",arm),1);
	T->SetBranchStatus(Form("D%s.evtypebits",arm),1);
	T->SetBranchStatus(Form("%s.vdc.u1.nclust",arm),1);
	T->SetBranchStatus(Form("%s.vdc.u2.nclust",arm),1);
	T->SetBranchStatus(Form("%s.vdc.v1.nclust",arm),1);
	T->SetBranchStatus(Form("%s.vdc.v2.nclust",arm),1);
	T->SetBranchStatus(Form("%s.prl1.trx",arm),1);
	T->SetBranchStatus(Form("%s.prl1.try",arm),1);
	T->SetBranchStatus(Form("%s.prl2.trx",arm),1);
	T->SetBranchStatus(Form("%s.prl2.try",arm),1);


	T->SetBranchAddress(Form("%s.tr.n",arm),&n_tr);
	T->SetBranchAddress(Form("%s.cer.asum_c",arm),&scer);
	T->SetBranchAddress(Form("%s.prl1.e",arm),&Eps);
	T->SetBranchAddress(Form("%s.prl2.e",arm),&Esh);
	T->SetBranchAddress(Form("%s.prl1.a",arm),&prl1_adc);
	T->SetBranchAddress(Form("%s.prl2.a",arm),&prl2_adc);
	T->SetBranchAddress(Form("%s.tr.p",arm),&p);
	T->SetBranchAddress(Form("D%s.evtypebits",arm),&evtype);
	T->SetBranchAddress(Form("%s.vdc.u1.nclust",arm),&n_u1_clust);
	T->SetBranchAddress(Form("%s.vdc.u2.nclust",arm),&n_u2_clust);
	T->SetBranchAddress(Form("%s.vdc.v1.nclust",arm),&n_v1_clust);
	T->SetBranchAddress(Form("%s.vdc.v2.nclust",arm),&n_v2_clust);
	T->SetBranchAddress(Form("%s.prl1.trx",arm),&prl1_trx);
	T->SetBranchAddress(Form("%s.prl1.try",arm),&prl1_try);
	T->SetBranchAddress(Form("%s.prl2.trx",arm),&prl2_trx);
	T->SetBranchAddress(Form("%s.prl2.try",arm),&prl2_try);

	cout << "Completed tree branch setting" << endl;

	// loop the events

	for (Int_t i_evt=0; i_evt<numEvents;i_evt++) {  
//	for (Int_t i_evt=0; i_evt<10000;i_evt++) {  
	  
		T->GetEntry(i_evt);
		// cout << "get entry done for event " << i_evt <<endl;
		
		Int_t evtypebits = evtype;
		
		if( n_tr == 1 && n_u1_clust == 1 && n_u2_clust == 1 && n_v1_clust == 1 && n_v2_clust == 1 &&
		    (Eps + 0.85*Esh) > 200. && (Eps + 0.85*Esh) < 4000. && Eps > 150. && Esh > 75. &&
		    scer[0] > 3000. && p[0] > 0. ) {

		  // 	Tritium cut levels    (Eps + 0.85*Esh) > 1500. && (Eps + 0.85*Esh) < 4000. && Eps > 250. && Esh > 750. &&
		  //	    	    (evtypebits == 14) &&
		  //		    scer[0] > 5000. && p[0] > 0. ) {

			Double_t prl1_max_adc = 0.;
			Int_t prl1_max_block = -1;
			Double_t prl2_max_adc = 0.;
			Int_t prl2_max_block = -1;

			vector<pair<int, double> > block_adc;
	
			int goodBlocks = 0;

			Int_t prl1_hit_row = prl1_x_to_row(prl1_trx[0]);
			Int_t prl2_hit_row = prl2_x_to_row(prl2_trx[0]);
			Int_t prl1_hit_col = 0;
			Int_t prl2_hit_col = 0;
			Int_t irow;
			Int_t icol; 
			Int_t drow;
			Int_t dcol;

			for(Int_t block = 0; block<numBlockPRL1; block++) {
				irow = block%numRowPRL1;
				icol = block/numRowPRL1;
				drow = abs(prl1_hit_row - irow);
				dcol = abs(prl1_hit_col - icol);
				if(drow < 2 && dcol < 2){
				  // cout << "PRL1 good loop: block = " << block << endl;
				  Double_t adc_cor = (prl1_adc[block] - pedestals_PRL1[block])/1000.;
		
					if(adc_cor > 0.) {
		
						block_adc.push_back(make_pair(block, adc_cor));
						goodBlocks++;
		
					}
					if(adc_cor > prl1_max_adc) {
						prl1_max_adc = adc_cor;	
						prl1_max_block = block;
					}
				}
			} 

			//			cout << "processed loop for prl1" << endl;
			for(Int_t block = 0; block<numBlockPRL2; block++) {
				irow = block%numRowPRL2;
				icol = block/numRowPRL2;
				drow = abs(prl2_hit_row - irow);
				dcol = abs(prl2_hit_col - icol);
				if(drow < 2 && dcol < 2){
				  //				   cout << "PRL2 good loop: block = " << block << endl;
					Double_t adc_cor = (prl2_adc[block] - pedestals_PRL2[block])/1000.;
					if(adc_cor > 0) {
					  //					  cout << "entered first loop" << endl;
					  block_adc.push_back(make_pair(block+numBlockPRL1, adc_cor));
					  goodBlocks++;										       
					}
					if(adc_cor > prl2_max_adc) {
					  //					  cout << "entered 2nd loop" << endl;
					  prl2_max_adc = adc_cor;	
					  prl2_max_block = block;
					}
					//					cout << "PRL2 end of good loop: block = " << block << endl;
				}

				//				cout << "PRL2: processed for block " << block << endl;
			}
			//			cout << "processed loop for prl2" << endl;

			Bool_t edge_max_block = false;

			for(Int_t i = 0; i < Nedge; i++) {
				if(edge[i] == prl1_max_block || edge[i] == prl2_max_block) {
					edge_max_block = true;
				}
			}

			if(!edge_max_block) {
				numGoodEvents++;
				goodEvents.push_back(make_pair(make_pair(goodBlocks, p[0]), block_adc));
			}

			block_adc.clear();
			
		}// end good event if statement

		if (i_evt%10000==0){cout << "Processed " << i_evt << "/" << numEvents << " events \t " << numGoodEvents << " good"  << endl;}

	} // end event for loop	

	cout << "DONE!  There are " << numGoodEvents << " good events!\n\n";
}


void chi2(int &nDim, double *gout, double &result, double par[], int flg) {


	double sumE = 0.;
	double deltaE = 0.;
	double sumDeltaE2 = 0.;
	
	for(int evt = 0; evt < numGoodEvents; evt++) {
//		cout << "chi2 looking at event " << evt << endl;
		sumE = 0.;
		deltaE = 0.;

		int nBlocks = goodEvents[evt].first.first;
		double p = goodEvents[evt].first.second; 		

//		cout << "nBlocks = " << nBlocks << "\tP = " << p << endl;

		for(int blk = 0; blk < nBlocks; blk++) {

			int blockNum = goodEvents[evt].second[blk].first;
			double blockADC = goodEvents[evt].second[blk].second;

//			cout << "Block num = " << blockNum << "\tADC = " << blockADC << endl;
			sumE += par[blockNum]*blockADC;

		} //block loop
		
		deltaE = sumE - p;
		sumDeltaE2 += deltaE*deltaE;

	} //event loop

	result = sumDeltaE2;

	cout << fixed << sumDeltaE2 << endl;

}

void APEX_LHRS_minuit_projection_cut1(Int_t run = 4179) {

	//LHRS
	sprintf(arm, "L");
	numBlockPRL1 = 34;
	numBlockPRL2 = 34;	
	numRowPRL1 = 17;
	numColPRL1 = 2;
	numRowPRL2 = 17;
	numColPRL2 = 2;
	
/*
	Int_t runlist[9];
	runlist[0] = 90840;	// optics
	runlist[1] = 90849;	// hydrogen
	runlist[2] = 90850;	// hydrogen
	runlist[3] = 90855;	// helium
	runlist[4] = 90856;	// helium
	runlist[5] = 90874;	// tritium
	runlist[6] = 90875;	// tritium
	runlist[7] = 90876;	// deuterium
	runlist[8] = 90877;	// deuterium
	for(int i = 0; i<9; i++) {
		T->Add(Form("/home/tyler/tritium/rootfiles/tritium_%i*.root",runlist[i]));
	}
*/	


	T = Load_more_rootfiles(run);
	//TChain* T = Load_more_rootfiles(run);
		

//	T->Add("./rootfiles/tritium_*.root");
	
	numEvents = T->GetEntries();

	// read the pedestals
	ifstream preshowerPedFile("ped_prl1.dat");
	ifstream showerPedFile("ped_prl2.dat");
	Double_t prl1Ped;
	Double_t prl2Ped;

	while(preshowerPedFile >> prl1Ped){
		pedestals_PRL1.push_back(prl1Ped);
	}

	cout << "Read PRL1 pedestals!\n\n";

	while(showerPedFile >> prl2Ped) { 
		pedestals_PRL2.push_back(prl2Ped);
	}

	cout << "Read PRL2 pedestals!\n\n";

	getGoodEvents();
	cout << "completed getGoodEvents() function!" <<endl;


	TVirtualFitter::SetDefaultFitter("Minuit"); //default is Minuit (JW)


	//	TFitter* minimizer = new TFitter(numBlockPRL1 + numBlockPRL2);
	TVirtualFitter* minimizer = TVirtualFitter::Fitter(NULL, numBlockPRL1 + numBlockPRL2);
	minimizer->SetFCN(chi2);

	for(Int_t i = 0; i < numBlockPRL1 + numBlockPRL2; i++) {
		minimizer->SetParameter(i, Form("alpha%i",i+1), 1.0, 0.05, 0., 10.);
	}
	cout << "Minimizing..." << endl; 
	//	minimizer->ExecuteCommand("MINIMIZE",0,0);
	Double_t arglist[1] = {0};

 	minimizer->ExecuteCommand("MIGRAD",arglist,0);
	
	ofstream fit_preshower("minuit_projection_cut1.prl1.dat");
	ofstream fit_shower("minuit_projection_cut1.prl2.dat");
	cout << setprecision(5);
	fit_preshower << "L.prl1.pedestals =\n";
	for(int i = 0; i < numBlockPRL1; i++) {
		fit_preshower << setw(10) << pedestals_PRL1[i] << "\t";
		if((i+1)%numColPRL1==0) {
			fit_preshower << "\n";
		}
	}	
	fit_shower << "L.prl2.pedestals =\n";
	for(int i = 0; i < numBlockPRL2; i++) {
		fit_shower << setw(10) << pedestals_PRL2[i] << "\t";
		if((i+1)%numColPRL2==0) {
			fit_shower << "\n";
		}
	}
	cout << setprecision(7);
	fit_preshower << "L.prl1.gains =\n";
	for(int i = 0; i < numBlockPRL1; i++) {
		fit_preshower << setw(15) << minimizer->GetParameter(i) << "\t";
		if((i+1)%numColPRL1==0) {
			fit_preshower << "\n";
		}
	}
	fit_shower << "L.prl2.gains =\n";
	for(int i = 0; i < numBlockPRL2; i++) {
		fit_shower << setw(15) << minimizer->GetParameter(i+numBlockPRL1) << "\t";
		if((i+1)%numColPRL2==0) {
			fit_shower << "\n";
		}
	}
	fit_preshower.close();
	fit_preshower.close();
}



