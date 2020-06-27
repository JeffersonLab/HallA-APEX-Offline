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




//Bool_t hole_select[NHoles];



Bool_t foil_select[11];


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

Bool_t hole_select[11][NSieveCol][NSieveRow];

Bool_t col_select[11][NSieveCol];

Bool_t row_select[11][NSieveRow];



void hole_opt_select();
/* for(Int_t i = 0; i<NHoles; i++){ */


/*   hole_select[i] = false; */

/*  } */







#endif
