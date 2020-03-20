
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


//TString RootFileName = "apex_4771_30_1_2020.root";

//TString RootFileName = "apex_4771_16_1_2020.root";

//TString RootFileName = "apex_4774_27_1_2020.root";




// Optics targets cuts

/* TString RootFileName = "apex_4771_26_2_2020.root"; // Optics 3 */
// TString RootFileName = "apex_4774_26_2_2020.root"; // Optics 1


//TString RootFileName = "apex_4771_28_2_2020.root"; // Optics 3
//TString RootFileName = "apex_4774_28_2_2020.root"; // Optics 1




//Vertical Foil cuts

/* TString RootFileName = "apex_4766_27_1_2020.root"; //V1 */
/* TString RootFileName = "apex_4768_27_1_2020.root"; //V2 */
/* TString RootFileName = "apex_4769_27_1_2020.root"; //V3 */

/* TString RootFileName = "apex_4766_6_2_2020.root"; //V1 */
/* TString RootFileName = "apex_4768_6_2_2020.root"; //V2 */
/* TString RootFileName = "apex_4769_6_2_2020.root"; //V3 */

//TString RootFileName = "apex_4766_17_2_2020.root"; //V1
/* TString RootFileName = "apex_4768_17_2_2020.root"; //V2 */
/* TString RootFileName = "apex_4769_17_2_2020.root"; //V3 */

TString RootFileName = "apex_4766_29_2_2020.root"; //V1
/* TString RootFileName = "apex_4768_29_2_2020.root"; //V2  */
/* TString RootFileName = "apex_4769_29_2_2020.root"; //V3 */





// Vertical target runs

//TString RootFileName = "apex_4770_27_11_2019.root";

 
//TString ROOTFILE_DIR = "/w/work3/home/johnw/Rootfiles/Correlations/17_1_20_TPY_new_offset/%d_replay.root";


//~~~~~~~~

//TString ROOTFILE_DIR = "/w/work3/home/johnw/Rootfiles/Correlations/30_1_20_TPY/%d_replay.root";

//TString ROOTFILE_DIR = "/w/work3/home/johnw/Rootfiles/Correlations/V1_V2_TP/%d_replay.root";

//TString ROOTFILE_DIR = "/w/work3/home/johnw/Rootfiles/Correlations/V1_V2_V3_TPY/%d_replay.root";

//TString ROOTFILE_DIR = "/w/work3/home/johnw/Rootfiles/Correlations/V1_V2_V3_TPY_25_2_20/%d_replay.root";


//28/2/20
TString ROOTFILE_DIR = "/w/work3/home/johnw/Rootfiles/Correlations/all_foil_min7_TPY_3rd/%d_replay.root";


//~~~~~~~~

//TString ROOTFILE_DIR = "/w/work3/home/johnw/Rootfiles/Correlations/TP_13_1_20_Old_Y_no_col_0/%d_replay.root";


//TString ROOTFILE_DIR = "/w/work3/home/johnw/Rootfiles/apex_online_%d_22_11_2019.root";


//TString* SoureRootFile = new TString("/w/work3/home/johnw/Rootfiles/" + RootFileName);


//TString ROOTFILE_DIR = "/w/work3/home/johnw/Rootfiles/apex_%d_30_1_2020.root";

//TString ROOTFILE_DIR = "/w/work3/home/johnw/Rootfiles/apex_%d_27_1_2020.root";


TString* SoureRootFile = new TString("/w/work3/home/johnw/Rootfiles/" + RootFileName);

//TString* SoureRootFile = new TString("/w/work3/home/johnw/Rootfiles/apex_4771_2_11_2019.root");

Int_t Run_number = 4766;
Int_t Run_number_2 = 4766;






#endif
