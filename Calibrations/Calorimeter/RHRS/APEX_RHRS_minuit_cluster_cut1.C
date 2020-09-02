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

// number of events, blocks, and columns (needed for loops)
Int_t numEvents;
Int_t nBlkPreShower;
Int_t nBlkShower;
Int_t nColPreShower;
Int_t nColShower;
// which arm: L or R
char arm[2]; 
// TChain for data ROOT files
TChain *T = new TChain("T");

// vector of pedestal values for shower and preshower
// should be read in from text file
vector<Double_t> ped_preShow;		
vector<Double_t> ped_show;	

// vector of gain values for shower and preshower
// should be read in from text file
vector<Double_t> gain_preShow;		
vector<Double_t> gain_show;	


// number of good events (needed for loops)	
Int_t numGoodEvents = 0;

// good event vector
vector<pair<pair<int, double>, vector<pair<int, double> > > > goodEvents;
/* JW: vector of pairs of (  (  pairs of ints and doubles = no of goodblocks and the momentum value) +  ( vectors of (pairs of ints and doubles) = block no + ADC value) )
large vector is for all events
 */


// ############################################
// #### A NOTE ABOUT THE GOOD EVENT VECTOR ####
// ############################################
// for a single event with signal in 7 blocks (2 in preshower, 3 in shower), and track momentum P_tr, the vector entry is structured as:
// 	(7, P_tr)	(psBlock1, ADC_psBlock1)
// 			(psBlock2, ADC_psBlock2)
// 			(shBlock1, ADC_shBlock1)
// 			(shBlock2, ADC_shBlock2)
// 			(shBlock3, ADC_shBlock3) 
// 	(next event)	


void getGoodEvents() {

	// histograms to keep track of how many events were accepted for each block
	TH1I* preshowerEventHist = new TH1I("preshowerEventHist","",48,0,48); 
	TH1I* showerEventHist = new TH1I("showerEventHist","",75,0,75);
	
	// block numbers for blocks on the edge
	// don't accept events where the cluster center is on the edge
	Int_t ps_edge[] = {0, 24, 23, 47};
	Int_t sh_edge[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 30, 45, 29, 44, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74};
	Int_t ps_Nedge = sizeof(ps_edge)/sizeof(Int_t);
	Int_t sh_Nedge = sizeof(sh_edge)/sizeof(Int_t);

	// all the branches required for cuts and cluster information
	Double_t scer[100];
	Double_t Eps, Esh;
	Double_t n_tr;
	Double_t n_u1_clust;
	Double_t n_u2_clust;
	Double_t n_v1_clust;
	Double_t n_v2_clust;
	Double_t nMultShower;
	Double_t nMultPreShower;
	Double_t ps_nBlk[100];
	Double_t sh_nBlk[100];
	Double_t p[100];
	Double_t ps_adc[100];
	Double_t sh_adc[100];	
	Double_t evtype;
	
	T->SetBranchStatus("*",0);  

	T->SetBranchStatus(Form("%s.tr.n",arm),1);
	T->SetBranchStatus(Form("%s.cer.asum_c",arm),1);
	T->SetBranchStatus(Form("%s.ps.e",arm),1);
	T->SetBranchStatus(Form("%s.sh.e",arm),1);
	T->SetBranchStatus(Form("%s.ps.mult",arm),1);
	T->SetBranchStatus(Form("%s.ps.nblk",arm),1);
	T->SetBranchStatus(Form("%s.ps.a",arm),1);
	T->SetBranchStatus(Form("%s.sh.mult",arm),1);
	T->SetBranchStatus(Form("%s.sh.nblk",arm),1);
	T->SetBranchStatus(Form("%s.sh.a",arm),1);
	T->SetBranchStatus(Form("%s.tr.p",arm),1);
	T->SetBranchStatus(Form("D%s.evtypebits",arm),1);
	T->SetBranchStatus(Form("%s.vdc.u1.nclust",arm),1);
	T->SetBranchStatus(Form("%s.vdc.u2.nclust",arm),1);
	T->SetBranchStatus(Form("%s.vdc.v1.nclust",arm),1);
	T->SetBranchStatus(Form("%s.vdc.v2.nclust",arm),1);

	T->SetBranchAddress(Form("%s.tr.n",arm),&n_tr);
	T->SetBranchAddress(Form("%s.cer.asum_c",arm),&scer);
	T->SetBranchAddress(Form("%s.ps.e",arm),&Eps);
	T->SetBranchAddress(Form("%s.sh.e",arm),&Esh);
	T->SetBranchAddress(Form("%s.ps.mult",arm),&nMultPreShower);
	T->SetBranchAddress(Form("%s.ps.nblk",arm),&ps_nBlk);
	T->SetBranchAddress(Form("%s.ps.a",arm),&ps_adc);
	T->SetBranchAddress(Form("%s.sh.mult",arm),&nMultShower);
	T->SetBranchAddress(Form("%s.sh.nblk",arm),&sh_nBlk);
	T->SetBranchAddress(Form("%s.sh.a",arm),&sh_adc);
	T->SetBranchAddress(Form("%s.tr.p",arm),&p);
	T->SetBranchAddress(Form("D%s.evtypebits",arm),&evtype);
	T->SetBranchAddress(Form("%s.vdc.u1.nclust",arm),&n_u1_clust);
	T->SetBranchAddress(Form("%s.vdc.u2.nclust",arm),&n_u2_clust);
	T->SetBranchAddress(Form("%s.vdc.v1.nclust",arm),&n_v1_clust);
	T->SetBranchAddress(Form("%s.vdc.v2.nclust",arm),&n_v2_clust);


	// loop the all events


	for (Int_t i_evt=0; i_evt<numEvents;i_evt++) {  
//	for (Int_t i_evt=0; i_evt<1000000;i_evt++) {  
	
		T->GetEntry(i_evt);

		Int_t evtypebits = evtype;

		// good event cuts
		// JW: first cuts based on track info: no of tracks and tracks in both planes of each VDC chamber
		// JW: 2nd line based on energy deposited in shower and post-shower detectors (characteristic of electron) 
		// JW : evtypebits are trigger conditions 
		// JW: scer is cerenkov sum (again for electron PID)
		
		if( n_tr == 1 && n_u1_clust == 1 && n_u2_clust == 1 && n_v1_clust == 1 && n_v2_clust == 1 &&
		    (Eps + 0.85*Esh) > 500. && (Eps + 0.85*Esh) < 2000. && Eps > 100. && Esh > 100. && scer[0] > 3000. && p[0] > 0. ) {


		  // tritium cut levels
		  // (Eps + 0.85*Esh) > 450. && (Eps + 0.85*Esh) < 4000. && Eps > 250. && Esh > 200. && scer[0] > 3000. && p[0] > 0. ) {

		    //	    	    (evtypebits == 48 || evtypebits == 112) &&
		    // scer[0] > 3000. && p[0] > 0. ) {


		  // below are the levels of the tritium cuts
		  //		    (Eps + 0.85*Esh) > 1500. && (Eps + 0.85*Esh) < 4000. && Eps > 250. && Esh > 750. &&
		  //		    scer[0] > 5000. && p[0] > 0. ) {

			Double_t ps_max_adc = 0.;
			Int_t ps_max_block = -1;
			Double_t sh_max_adc = 0.;
			Int_t sh_max_block = -1;

			vector<pair<int, double> > block_adc;
	
			int goodBlocks = 0;

			// loop over blocks in preshower cluster
			Int_t nMultPS = nMultPreShower; // JW: blocks in main cluster
			for(Int_t n_ps = 0; n_ps < nMultPS; n_ps++) {
			  Int_t blockNum = ps_nBlk[n_ps]; // JW: No of block in the cluster
				double adc_cor = (ps_adc[blockNum] - ped_preShow[blockNum])/1000.;
				// keep block information if pedestal-subtracted ADC signal is greater than 0
				if(adc_cor > 0.) {
					block_adc.push_back(make_pair(blockNum, adc_cor));	
					// JW: adds block number and it adc_corrected value to a list (for each event)
					goodBlocks++;
				}
				// keep track of which block saw the maximum signal
				if(adc_cor > ps_max_adc) {
					ps_max_adc = adc_cor;	
					ps_max_block = blockNum;
				}
			}	
			
			// loop over block in shower cluster
			Int_t nMultSH = nMultShower;
			for(Int_t n_sh = 0; n_sh < nMultSH; n_sh++) {
				Int_t blockNum = sh_nBlk[n_sh];
				double adc_cor = (sh_adc[blockNum] - ped_show[blockNum])/1000.;
				if(adc_cor > 0.) {
					block_adc.push_back(make_pair(blockNum+nBlkPreShower, adc_cor));	
					goodBlocks++;
				}
				if(adc_cor > sh_max_adc) {
					sh_max_adc = adc_cor;	
					sh_max_block = blockNum;
				}
			}

			// check if block with maximum signal is on the edge
			// if so, reject event
			Bool_t ps_edge_max_block = false;
			Bool_t sh_edge_max_block = false;

			for(Int_t i = 0; i < ps_Nedge; i++) {
				if(ps_edge[i] == ps_max_block) {
					ps_edge_max_block = true;
				}
			}

			for(Int_t i = 0; i < sh_Nedge; i++) {
				if(sh_edge[i] == sh_max_block) {
					sh_edge_max_block = true;
				}
			}

			// if good event, add event information to good event vector
			if(!ps_edge_max_block && !sh_edge_max_block) {
				numGoodEvents++;
				goodEvents.push_back(make_pair(make_pair(goodBlocks, p[0]), block_adc));

				// JW: this loop simply fills no of good events for pre and shower detectores
				for(Int_t i = 0; i < goodBlocks; i++) {
					Int_t thisBlock = block_adc[i].first;
					if(thisBlock < nBlkPreShower) {
						preshowerEventHist->Fill(thisBlock);						
					} else {
						showerEventHist->Fill(thisBlock - nBlkPreShower);
					}	
				}
			}

			block_adc.clear();
			
		}// end good event if statement

		if (i_evt%10000==0){cout << "Processed " << i_evt << "/" << numEvents << " events \t " << numGoodEvents << " good"  << endl;}

	} // end event for loop	

	cout << "DONE!  There are " << numGoodEvents << " good events!\n\n";
	TCanvas* can = new TCanvas("can","",800,600);
	gPad->SetLogy();
	preshowerEventHist->Draw();
	can->SaveAs("preshowerBlockEvents.pdf");
	showerEventHist->Draw();
	can->SaveAs("showerBlockEvents.pdf");

}


void chi2(int &nDim, double *gout, double &result, double par[], int flg) {

	// cumulative sum of block energies (for one event)
	double sumE = 0.;
	// difference between sum of block energies and track momentum (for one event)
	double deltaE = 0.;
	// cumulative sum of deltaE^2, for all events 
	double sumDeltaE2 = 0.;
	
	// loop over good events
	for(int evt = 0; evt < numGoodEvents; evt++) {
		sumE = 0.;
		deltaE = 0.;

		int nBlocks = goodEvents[evt].first.first;
		double p = goodEvents[evt].first.second; 		

		// loop over each block and add energy to sumE
		for(int blk = 0; blk < nBlocks; blk++) {

			int blockNum = goodEvents[evt].second[blk].first;
			double blockADC = goodEvents[evt].second[blk].second;

			sumE += par[blockNum]*blockADC;

		} //block loop
		
		// update sumDeltaE2 for this event
		deltaE = sumE - p;
		// JW: translates to par[i]*A(ij) - p(j) (i over blocks, j over events)
		// where p is momentum 
		sumDeltaE2 += deltaE*deltaE;

	} //event loop

	// return sumDeltaE2 for this set of block gains
	result = sumDeltaE2;

	// print out result (fun to watch it converge)
	cout << fixed << sumDeltaE2 << endl;

}

void APEX_RHRS_minuit_cluster_cut1(Int_t run = 4179) {

	//RHRS
	sprintf(arm, "R");
	nBlkPreShower = 48;
	nBlkShower = 75;	
	nColPreShower = 2;
	nColShower = 5;
	
	// add all rootfiles in some directory to TChain
	//	T->Add("./rootfiles/tritium_online_*");	
	T = Load_more_rootfiles(run);
	


	numEvents = T->GetEntries();

	// read the pedestals
	ifstream preshowerPedFile(Form("output_data/%d_ped_preshower.dat",run));
	ifstream showerPedFile(Form("output_data/%d_ped_shower.dat",run));
	Double_t preShowPed;
	Double_t showPed;

	// read the preshower pedestals from file
	while(preshowerPedFile >> preShowPed){
		ped_preShow.push_back(preShowPed);
	}

	// read the shower pedestals from file
	while(showerPedFile >> showPed) { 
		ped_show.push_back(showPed);
	}



	// read input gain files to use as first guess of parameters

	ifstream preshowerGainsFile("ps_gains.dat");
	ifstream showerGainsFile("sh_gains.dat");
	Double_t psGain;
	Double_t shGain;

	Int_t counter = 0;

	while(preshowerGainsFile >> psGain){
	  gain_preShow.push_back(psGain);
	}

	while(showerGainsFile >> shGain) { 
	  gain_show.push_back(shGain);
	}




	// fill good event vector with good events
	getGoodEvents();

	// create minimizer and set chi2 as the function to be minimized	
	//	TFitter* minimizer = new TFitter(nBlkPreShower + nBlkShower);
	TVirtualFitter* minimizer = TVirtualFitter::Fitter(NULL, nBlkPreShower + nBlkShower);
	minimizer->SetFCN(chi2);
	//	minimizer->SetFitOption(q);



	// use input gains (read from old DB entry) as initial guess of gain parameters
	for(Int_t i = 0; i < nBlkPreShower + nBlkShower; i++) {

	  if (i < nBlkPreShower){
	    minimizer->SetParameter(i, Form("alpha%i",i+1), gain_preShow[i], 0.05, 0., 10.);
	  }
	  else{
	    minimizer->SetParameter(i, Form("alpha%i",i+1), gain_show[i - nBlkPreShower], 0.05, 0., 10.);
	  }
	  
	  
	  
	}
	
	// complete the minimization
	cout << "Minimizing..." << endl; 
	//	minimizer->ExecuteCommand("MINIMIZE",0,0);
 	minimizer->ExecuteCommand("MIGRAD",0,0);
	
	// JW: think these should return int but allowed not to (root 5 ?)

	// write out to file in current directory
	// the output has been structured with the correct ordering of the block gains
	// i.e., you should be able to just copy and paste this into your database file
	ofstream fit_preshower(Form("output_data/%d_minuit_cluster_cut1.ps.dat",run));
	ofstream fit_shower(Form("output_data/%d_minuit_cluster_cut1.sh.dat",run));
	cout << setprecision(5);

	fit_preshower << "R.ps.pedestals =\n";
	for(int i = 0; i < nBlkPreShower; i++) {
		fit_preshower << setw(10) << ped_preShow[i] << "\t";
		if((i+1)%nColPreShower==0) {
			fit_preshower << "\n";
		}
	}	
	fit_shower << "R.sh.pedestals =\n";
	for(int i = 0; i < nBlkShower; i++) {
		fit_shower << setw(10) << ped_show[i] << "\t";
		if((i+1)%nColShower==0) {
			fit_shower << "\n";
		}
	}
	cout << setprecision(7);
	fit_preshower << "R.ps.gains =\n";
	for(int i = 0; i < nBlkPreShower; i++) {
		fit_preshower << setw(15) << minimizer->GetParameter(i) << "\t";
		if((i+1)%nColPreShower==0) {
			fit_preshower << "\n";
		}
	}

	fit_shower << "R.sh.gains =\n";
	for(int i = 0; i < nBlkShower; i++) {
		fit_shower << setw(15) << minimizer->GetParameter(i+nBlkPreShower) << "\t";
		if((i+1)%nColShower==0) {
			fit_shower << "\n";
		}
	}
	fit_preshower.close();
	fit_preshower.close();

	// done!
}



