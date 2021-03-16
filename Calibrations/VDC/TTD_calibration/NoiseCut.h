#ifndef NOISECUT_H
#define NOISECUT_H

/*
  Class desinged to read 2D historgramm, call FitSlicesY() which fits Gaussian slices in Y, and store the resulting historgrams of constants, means and sigmas as vectors. 

Can then be used to test if vector of data passes conditions of noise defined by Gaussian fit: if the data lies between the mean plus adjustable value times sigmae. Designed to remove noise that comes from unphysical noise (times that correspond to vertical distances which are too large or small to be real tracks). Can be used for generic 2D histogram, to then cut vectors (need access to same two variables as in 2D histogram).


*/


#include "TH2.h"
#include "TH1.h"
#include "TDirectory.h"

#include <iostream>

using namespace std;

class NoiseCut {

 public:
  NoiseCut(TH2D *histo, int firstxbin = 0, int lastxbin = -1, int cut = 10 /*Min num of entries in slices*/);


  void ReadHist(TH1D* hist,  std::vector<Double_t>& vec);
  

  void PassNoiseCut(std::vector<double>& InVectX, std::vector<double>& InVectY);

  bool PassNoiseCut(double X, double Y);
  
  
 private:
  std::vector<double> consts; // constants from gaussian fits of slices (peak); will scale with entries in histogram
  std::vector<double> means;  // means from gaussian fits of slices
  std::vector<double> sigmas; // sigmas from gaussian fits of slices
  double Low; // lowest value of 'x' var
  double High; // lowest value of 'x' var
  double BinWid; // width of bins
  double XSigmaUp = 3.0; // value of sigma to cut above
  double XSigmaDown = 2.0; // value of sigma to cut below
};




#endif
