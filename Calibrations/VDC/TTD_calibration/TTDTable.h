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
  TTDTable(std::vector<Double_t> Table, double LowVal, double *apars): LTable(Table), Low(LowVal), par(apars) {};
  TTDTable(std::vector<Double_t> Table, double LowVal): LTable(Table),  Low(LowVal), par(NULL) {};

  double Convert(double dtime);
  double ConvertAngleCorr(double dtime, double tantheta);
  void PrintParams();


 private:
  std::vector<double> LTable; // lookup table
  const double * par; // parameters for angle correction
  const double invTanTheta0 = 1.4; // central value of inverse tan theta (inverse of slope)
  const double Low; // value of lowest time 

};

#endif
