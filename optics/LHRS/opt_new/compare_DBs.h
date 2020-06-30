///////////////////////////////////////////////////////////////////////////////
//
// Header for script designed to compare and plots two optics matrix BDs
//
//

// Author: John Williamson <jwilli@jlab.org>
// 
//  29th June 2020
//
///////////////////////////////////////////////////////////////////////////////


enum{
     T = 0,
     P,
     Y,
     D
};

struct ME{
  Int_t Order; // order of element
  Int_t type; // type of element (T (theta),P (phi), Y (y) or D (momentum deviation))
  Int_t x_pow; // indice of x_FP
  Int_t y_pow; // indice of y_FP
  Int_t th_pow; // indice of theta_FP
  Int_t ph_pow; // indice of phi_FP
  Double_t Coeff; // coeffecient of element
};


// struct designed to hold vector of MEs for particular DB and the largest order element
struct DB_MEs{
  vector <ME> MEs; // vector of MEs
  Int_t max_order; // largest order element in vector of MEs
};


TString element_string[4];


// type is vector of coeffecients and names for Matrix Elements (MEs). Usually a vector is used to store MEs of a particular order.
typedef vector<pair<Double_t, TString>> ME_Names;


// type is used to store vector of vectors of Matrix ELement (ME) coeffecient and names for diferrent orders 
typedef vector<ME_Names> ME_Ordered;


// used to read Database into DB_MEs 
DB_MEs ReadDB(LOpticsOpt* DB, Int_t variable);



// function that extracts vectors of co-effecients and names of elements for various orders
// produces 0th, 1st, 2nd etc order vectors
void extract_elems(DB_MEs &ME_vect, vector<vector<pair<Double_t, TString>>> &Coeff_names);



// function designed to compare MEs from two DBs and remove elements which are not common to both and return ME_Ordered with differences of two MEs
// flag == 0 : returns proportional difference
// flag == 1 : returns absolute difference
ME_Ordered compare_elems(ME_Ordered &Coeff_names_1, ME_Ordered &Coeff_names_2, Int_t flag = 0);



// compare_elems function used for individual orders
ME_Names compare_elems_order(ME_Names &Order_1, ME_Names &Order_2, Int_t flag = 0);



// function used to fill graphs for each order from ME_Ordered
void fill_graphs(vector<TGraph*> &graphs, ME_Ordered coeff_names);

