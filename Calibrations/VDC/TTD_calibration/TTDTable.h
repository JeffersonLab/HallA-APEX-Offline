#ifndef TTDTable_H
#define TTDTable_H

/*
  Class designed to store TTD lookup table and perform lookup up conversion to distance given time (and slope)


 */

#include <iostream>
#include "TMath.h"

using namespace std;

class TTDTable {

 public:
  TTDTable(std::vector<Double_t> Table, double LowVal, int nbins, double *apars, double *aparsext): LTable(Table), Low(LowVal), NBins(nbins), par(apars), parExt(aparsext), BExtPars(true) {};
  TTDTable(std::vector<Double_t> Table, double LowVal, int nbins): LTable(Table),  Low(LowVal), NBins(nbins), par(NULL), parExt(NULL) {};


  double Convert(double dtime);
  double ConvertAngleCorr(double dtime, double tantheta);
  void PrintParams();


 private:
  std::vector<double> LTable; // lookup table
  const double * par; // parameters for angle correction
  const double invTanTheta0 = 1.4; // central value of inverse tan theta (inverse of slope)
  const double Low; // value of lowest time
  const int    NBins; // number of entries in Lookup table
  const double bin_res = 0.5e-9;
  const double * parExt; // parameters for times outside normal range
  Bool_t BExtPars = false; // bool to keep track if ext parametershave been intialised

};

#endif
