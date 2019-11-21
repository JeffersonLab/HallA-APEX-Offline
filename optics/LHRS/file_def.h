#ifndef ROOT_file_def
#define ROOT_file_def


#include "TObject.h" //trick to avoid error: `Int_t' does not name a type



TString file_date = "2_11_2019";

//TString ROOTFILE_DIR = "/w/work3a/home/johnw/Rootfiles/apex_%d_2_11_2019.root";

/* TString DB_name = " */
TString ROOTFILE_DIR = "/w/work3/home/johnw/Rootfiles/Correlations/";



//TString ROOTFILE_DIR = "/w/work3/home/johnw/Rootfiles/apex_%d_15_10_2019.root";
//TString ROOTFILE_DIR = "/w/work3/home/johnw/Rootfiles/apex_%d_10_9_2019.root";
//TString ROOTFILE_DIR = "/w/work3/home/johnw/Rootfiles/apex_%d_3_9_2019.root";
//TString ROOTFILE_DIR = "/w/work3/home/johnw/Rootfiles/apex_%d_22_8_2019_rast.root";

//TString ROOTFILE_DIR = "/w/work3/home/johnw/Rootfiles/apex_%d_10_7_2019.root";
/* TString ROOTFILE_DIR = "/w/work3/home/johnw/Rootfiles/apex_%d_2010.root"; */

TString RootFileName = "apex_4181_2_11_2019.root";
 

TString* SoureRootFile = new TString("/w/work3/home/johnw/Rootfiles/" + RootFileName);


//TString* SoureRootFile = new TString("/w/work3/home/johnw/Rootfiles/apex_4181_3_9_2019.root");
/* TString* SoureRootFile = new TString("/w/work3/home/johnw/Rootfiles/apex_4771_2010.root"); */

Int_t Run_number = 4181;
Int_t Run_number_2 = 4181;



/* Following section designed to deal with splitting data into Dp fractions where
   Dp = (p - p_0_/ p

 */
// where the zeroth element of Dp is -5% to -4 % through to the ninth element 4% to 5%
/*
Dp[0] = -5% -> -4%
Dp[1] = -4% -> -3%
Dp[2] = -3% -> -2%
Dp[3] = -2% -> -1%
Dp[4] = -1% ->  0%
Dp[5] =  0% ->  1%
Dp[6] =  1% ->  2%
Dp[7] =  2% ->  3%
Dp[8] =  3% ->  4%
Dp[9] =  4% ->  5%




 */

//Bool_t* Dp[10] = {NULL};

TString SoureRootFile_dp = "/w/work3/home/johnw/Rootfiles/apex_4181_22_8_2019_rast_dp_%d.root";
/* TString SoureRootFile_dp = "/w/work3/home/johnw/Rootfiles/apex_4181_2010_dp_%d.root"; */

Bool_t Dp[10] = {false};

void Dp_choose(Int_t j){

  //  Bool_t Dp[10] = {0};

  for (Int_t i = 0; i<sizeof(Dp)/sizeof(Dp[0]); i++){
    
    Dp[i] = false;
  }
  
  Dp[j] = true;
  
  
  
  
  // SoureRootFile_dp = 
  
  for(Int_t i = 0; i < sizeof(Dp)/sizeof(Dp[0]); i++){
    
    if(Dp[i]){
      *SoureRootFile = Form(SoureRootFile_dp,i);
    }
    
  }
 


}

//SoureRootFile = Dp_choose();

#endif
