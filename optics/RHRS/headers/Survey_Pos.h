#include "TRotation.h"




TRotation fTCSInHCS;
TVector3 fPointingOffset;

void Initialize(void);

void Hole_exp(TString run,int n_foil, int &row_min,int &row_max,int &col_min,int &col_max);

TVector3 BeamSpotHCS_Correction(UInt_t FoilID, double beam_y, double beam_z);
const TVector3 GetSieveHoleTCS(UInt_t Col, UInt_t Row);
const TVector3 GetSieveHoleCorrectionTCS(UInt_t nfoil, UInt_t Col, UInt_t Row);


void Sieve_hole_pos(int FoilID, int Col, int Row, double sieve_ph_th[], double sieve_yx[]);
//void Sieve_hole_pos(int FoilID, int Col, int Row);
