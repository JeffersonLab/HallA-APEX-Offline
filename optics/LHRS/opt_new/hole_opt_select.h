///////////////////////////////////////////////////////////////////////////////
//
// Allows chosen hole not to be used in optimisation
//
//

// Author: John Williamson <jwilli@jlab.org>
// Modification:
//			Jan 13th 2020
//
///////////////////////////////////////////////////////////////////////////////



#ifndef HOLE_OPT
#define HOLE_OPT


/* switches which change way hole selection works
 kcopy_old: this uses Old_DataSource in LOpticsOptScript.C to read in the holes used by a previous optimisation
The file read in has to first be produced by LOpticsOpt::Print_holes(TString DataSource)

khole_opt_select: this uses series of foils, columns, rows, holes and combinations of the above to define which holes shoudl be ignored in an optimisation

*/ 


enum hole_method{
		 kcopy_old = 0,
		 khole_opt_select
};

Int_t method = kcopy_old;


//Bool_t hole_select[NHoles];



Bool_t foil_select[11] = {true};


// create struct for holes (based on foil, column and row number)

struct holes_fcr {
  Int_t Foil; // Foil number for hole
  Int_t Col; // Col number for hole
  Int_t Row; // Row number for hole
};

struct cols_fc {
  Int_t Foil; // Foil number for col
  Int_t Col; // Col number for col
};

struct rows_fr {
  Int_t Foil; // Foil number for col
  Int_t Row; // Row number for col
};



// create struct for columns-foil combinations (based on foil and column numbers)

Bool_t hole_select[11][NSieveCol][NSieveRow] = {false};

Bool_t col_select[11][NSieveCol] = {false};

Bool_t row_select[11][NSieveRow] = {false};





void hole_opt_select();


void hole_select_copy(TString File_name);


/* for(Int_t i = 0; i<NHoles; i++){ */


/*   hole_select[i] = false; */

/*  } */







#endif
