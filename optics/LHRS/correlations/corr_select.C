////////////////////////////////////////////////////
//          corr_select
//
//   Script designed to select a hole, column or row
//   for a particular run and DB and use the corr_plot.C 
//   script to plot various correlation variables
//
//
//  John Williamson
//  28/10/2019
///////////////////////////////////////////////////




#include "TChain.h"
#include "corr_plot.C"
#include "Load_new_replay.C"

#include <iostream>
#include <cstdlib>
#include "file_def.h"
#include "TCutG.h"
#include "TCut.h"

#include "InputAPEXL.h"
#include "APEX_Sieve.h"




//void corr_select(TString DB_name, Int_t runnumber = 0, Int_t colin = -1, Int_t rowin = -1){

//void corr_select(DB_info DB_run, event_info event = {4179,1,-1,-1}){
void corr_select(Run_spec run_info){


  TString DB_name = run_info.DB_info[DB_NAME];
  TString Cut_name = run_info.DB_info[CUT_NAME];
  
  Int_t runnumber = run_info.event_info[RUN_NO];
  Int_t foil_no = run_info.event_info[FOIL_NO];
  Int_t colin = run_info.event_info[COL_NO];
  Int_t rowin = run_info.event_info[ROW_NO];


  
  if( runnumber == 0){

    cout << "Select run number: " << endl;
   
    cin >> runnumber;
  }


  Char_t choice = {0};

  

  
  if( rowin == -1 && colin == -1){
  
  cout << "Plot all holes (A), Plot all data (E) (E for Everything), select a column (C), select a row (R) or select a particular hole (H): " << endl;

  cin >> choice;
  }




  // set-up general track and PID cuts (taken from cut_L.C inthe opics scripts)


  // TCut GeneralSieveCut ="L.tr.n==1 && L.tr.chi2<0.003 && abs(L.gold.th)<0.05 && L.gold.ph>-0.07 && L.gold.ph<0.025 && abs(L.gold.dp)<0.05 && L.vdc.u1.nclust==1 && L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.u1.nclust==1";

  TCut GeneralSieveCut ="L.tr.n==1 && L.tr.chi2<0.003 && abs(L.gold.th)<0.05 && L.gold.ph>-0.07 && L.gold.ph<0.025 && abs(L.tr.r_x)<0.1 && L.vdc.u1.nclust==1 && L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.u1.nclust==1";

  TCut PID_cuts = "(L.prl1.e/(L.gold.p*1000))>0.2 && ((L.prl1.e+L.prl2.e)/(L.gold.p*1000))>0.51 && L.cer.asum_c>400";






  Int_t FoilID = 0;

  if(runnumber == 4179){
    FoilID = 0;
  }
  else if(runnumber == 4181){
    FoilID = 1;
  }
  else if(runnumber == 4180){
    FoilID = 2;
  }

  // cuts based on foil number
 

  
  //  *SoureRootFile = TString("newfile");

 

  //  TString CutFileName = *SoureRootFile + ".FullCut.root";


  TString CutFileName = Form("/w/work3/home/johnw/Rootfiles/apex_%d",runnumber) +  Cut_name + ".FullCut.root";

  cout << "~~~~~~~##CutFileName = " << CutFileName << endl;


  TFile *f1 = new TFile(CutFileName, "READ");
   //  TFile *f1 = new TFile(CutFileName, "READ"); // changed to read


  TString foil_cut = Form("fcut_L_%d", foil_no);

  TCutG* foil_cutG;

  

  //  (TCutG*) f1->GetObject(foil_cut, "TCutG");

  f1->GetObject(foil_cut, foil_cutG);


  //  TFile *new_file = new TFile("temp.root","RECREATE");


  
  TChain* T;




  // get choices from function input 

  if( rowin == -1 && colin >= 0){
    choice = 'C';
  }
  else if( rowin >= 0 && colin == -1){
    choice = 'R';
  }
  else if( rowin == 0 && colin == 0){
    choice = 'A';  
  }
  else if( rowin == 1 && colin == 1){
    choice = 'E';  
  }  
  else if( rowin >= 0 && colin >= 0){
    choice = 'H';
  }
  

  cout << " choice = " << choice << endl;


  //  TString foil_cut = Form("fcut_L_%d", foil_no);

  //  if( choice == (Char_t)'A' ){
  switch(choice) {
  case 'A':
    {
      //      Char_t all_choice = {0};
      cout << "Plotting all data ... (This plots all the holes, in future option with no hole cuts may be added)" << endl << endl;
      // cout << "Would you like to plot all hole that have been cut? (C) Or would you like to plot all data without holes cuts? (N)" << endl;
      // cin >> all_choice;   


      
      // cycle through all holes and plot ALL of them



      TString final_cut;

      TString hole_string;
      

      //      TCutG* cutg[NHoles];

      T = Load_new_replay(DB_name,runnumber);


      TCutG* cutg[NSieveRow*NSieveCol];




      Bool_t first_hole = true;
      // TCutG* cutg;

      for(Int_t row_i = 0; row_i < NSieveRow; row_i++){
  	for(Int_t col_i = 0; col_i < NSieveCol; col_i++){
  	  //	  cout << "(NSieveCol*row_i) + (col_i) = " << (NSieveCol*row_i) + (col_i) << endl;

  	  hole_string = Form("e_hcut_L_%d_%d_%d", FoilID, col_i, row_i);
	  
  	  cutg[(NSieveCol*row_i) + (col_i)] = (TCutG*)gROOT->FindObject(hole_string);
  		  //	  if( cutg[Get_Hole(col_i,row_i)] ){
  	  if( cutg[(NSieveCol*row_i) + (col_i)] ){
	  
  	    if(first_hole){
  	      final_cut += " ( "  + hole_string;
  	      first_hole = false;
  	  }
  	    else{
  	      final_cut += " || " + hole_string;
  	    }
	  
  	  }	  
  	}
	


      }


      final_cut += ") && (" + foil_cut + " && " + (TString) PID_cuts + " && " + (TString) GeneralSieveCut + ")";

      cout << "final cut for all = " << final_cut << endl;


      if(first_hole == false){
  	T = (TChain*) T->CopyTree(final_cut);
	
      
	
	
  	corr_plot(DB_name, T, runnumber, -1, -1);
      }
      else{
	
  	cout << "no holes found, exiting script" << endl;
      }

      
      break;
    }    
  case 'E':
    {
      //      Char_t all_choice = {0};
      cout << "Plotting all data with (no hole or foil cuts)" << endl << endl;


      TString final_cut;

      TString hole_string;
      

      //      TCutG* cutg[NHoles];

      T = Load_new_replay(DB_name,runnumber);



      final_cut =  (TString) PID_cuts + " && " + (TString) GeneralSieveCut + " && " + foil_cut; 

      cout << "final cut for all = " << final_cut << endl;



      T = (TChain*) T->CopyTree(final_cut);
	
      cout << "After tree copy" << endl;

      
      //      corr_plot(DB_name, T, runnumber, 1, 1);
      corr_plot(run_info, T);
      
      
      break;
    }    
  case 'C':
    {

      Int_t col_no;

      if( colin == -1){
  	cout << "Choose column number to plot: " << endl;
  	cin >> col_no;
      }
      else{
  	col_no = colin;
      }


      T = Load_new_replay(DB_name,runnumber);


      // cycle through all rows for one column and test if hole cut exists (and add if it does)

      TString final_cut;

      TString hole_string;
      

      TCutG* cutg[NSieveRow];


      Bool_t first_hole = true;
      // TCutG* cutg;

      for(Int_t row_i = 0; row_i < NSieveRow; row_i++){

  	hole_string = Form("e_hcut_L_%d_%d_%d", FoilID, col_no, row_i);

  	cutg[row_i] = (TCutG*)gROOT->FindObject(hole_string);

  	if( cutg[row_i] ){

  	  if(first_hole){

  	    final_cut += " ( "  + hole_string;
  	    first_hole = false;
  	  }
  	  else{
  	    final_cut += " || " + hole_string;
  	  }

  	}	  
      }


      final_cut += ") && (" + foil_cut + " && " + (TString) PID_cuts + " && " + (TString) GeneralSieveCut + ")";

      cout << "final cut = " << final_cut << endl;
      


      if(first_hole == false){
  	cout << "Pre-tree copy entries = " << T->GetEntries() << endl;
  	T = (TChain*) T->CopyTree(final_cut);
  	cout << "Post-tree copy entries = "  << T->GetEntries() << endl;
      	
  	corr_plot(DB_name, T, runnumber,col_no,-1);	
      }
      else{
	
  	cout << "no holes found for this column: " << col_no << ", exiting script" << endl;
      }



      break;
    }
  case 'R':
    {
      Int_t row_no;

      if( rowin == -1){
  	cout << "Choose a row number to plot: " << endl;
  	cin >> row_no;
      }
      else{
  	row_no = rowin;
      }

      T = Load_new_replay(DB_name, runnumber);



      // cycle through all columns for one row and test if hole cut exists (and add if it does)

      TString final_cut;

      TString hole_string;
      

      TCutG* cutg[NSieveCol];


      Bool_t first_hole = true;
      // TCutG* cutg;

      for(Int_t col_i = 0; col_i < NSieveCol; col_i++){

  	hole_string = Form("e_hcut_L_%d_%d_%d", FoilID, col_i, row_no);

  	cutg[col_i] = (TCutG*)gROOT->FindObject(hole_string);

  	if( cutg[col_i] ){

  	  if(first_hole){

  	    final_cut += "(" + hole_string;
  	    first_hole = false;
  	  }
  	  else{
  	    final_cut += " || " + hole_string;
  	  }

  	}	  
      }


      final_cut += ") && (" + foil_cut + " && " + (TString) PID_cuts + " && " + (TString) GeneralSieveCut + ")";

      cout << "final cut = " << final_cut << endl;



      if(first_hole == false){
  	T = (TChain*) T->CopyTree(final_cut);
      	
  	corr_plot(DB_name, T, runnumber,-1,row_no);
      }
      else{
	
  	cout << "no holes found for this row, exiting script" << endl;
      }
            


      break;
    }
  case 'H':
    {
      Int_t col_no;
      Int_t row_no;


      
      if( rowin == -1 && colin == -1){
  	cout << "Choose column number of hole to plot: " << endl;
  	cin >> col_no;
  	cout << "Choose row number of hole to plot: " << endl;
  	cin >> row_no;
  	cout << "hole row, column = (" << row_no << "," << col_no << ")" << endl;
      }
      else{
  	col_no = colin;
  	row_no = rowin;
      }

      T = Load_new_replay(DB_name, runnumber);






      // basic logic: load in cut and then use to copy version of TTree


      TString hole_string = Form("e_hcut_L_%d_%d_%d", FoilID, col_no, row_no);

      TCutG* cutg = (TCutG*)gROOT->FindObject(hole_string);

      if(cutg){
  	cout << "Found cut :)" << endl;
	

  	TString combined_cut = hole_string + " && " + foil_cut + " && " + (TString) PID_cuts + " && " + (TString) GeneralSieveCut;
  	cout << "combined cut = " << combined_cut << endl;

	
  	T = (TChain*) T->CopyTree(combined_cut);

  	corr_plot(DB_name, T, runnumber,col_no,row_no);
	

      }
      else{
  	cout << "didn't find cut :(" << endl;
      }

      
      
      break;
    }
  default:
    {
    cout << "Invalid choice" << endl;
    }
  }

  cout << "pre-file close" << endl;
  f1->Close();
  cout << "post-file close" << endl;
  //  new_file->Close();



}
