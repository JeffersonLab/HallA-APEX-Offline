
#ifndef ROOT_file_def
#define ROOT_file_def


#include "TObject.h" //trick to avoid error: `Int_t' does not name a type


// older cut files
//TString file_date = "27_11_2019";

// newer optimisation cut files
TString file_date = "16_1_2020";



//TString ROOTFILE_DIR = "/w/work3/home/johnw/Rootfiles/Correlations/16_1_20_TP/%d_replay.root";



//TString RootFileName = "apex_4771_27_11_2019.root";
//TString RootFileName = "apex_4774_27_11_2019.root";


TString RootFileName = "apex_4771_16_1_2020.root";
//TString RootFileName = "apex_4774_16_1_2020.root";


// Vertical target runs

//TString RootFileName = "apex_4770_27_11_2019.root";
 

//TString ROOTFILE_DIR = "/w/work3/home/johnw/Rootfiles/Correlations/17_1_20_TPY_new_offset/%d_replay.root";

TString ROOTFILE_DIR = "/w/work3/home/johnw/Rootfiles/apex_%d_27_1_2020.root";


TString* SoureRootFile = new TString("/w/work3/home/johnw/Rootfiles/" + RootFileName);

//TString* SoureRootFile = new TString("/w/work3/home/johnw/Rootfiles/" + RootFileName);


//TString* SoureRootFile = new TString("/w/work3/home/johnw/Rootfiles/apex_4771_27_1_2020.root");

Int_t Run_number = 4771;
Int_t Run_number_2 = 4771;






#endif
