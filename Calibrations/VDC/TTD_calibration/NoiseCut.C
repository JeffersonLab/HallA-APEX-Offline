#include "NoiseCut.h"


NoiseCut::NoiseCut(TH2D *hist, int firstxbin, int lastxbin, int cut /*Min num of entries in slices*/)
{

  // Attempt to Call FitSlicesY on histogram, this fit Gaussains on slices in Y and returns parameters of Gaussian fits
  
  hist->FitSlicesY(0,0,-1,10,"QN");

  TH1D* hconst = (TH1D*)gDirectory->Get(Form("%s_0",hist->GetName()));
  TH1D* hmean = (TH1D*)gDirectory->Get(Form("%s_1",hist->GetName()));
  TH1D* hsigma = (TH1D*)gDirectory->Get(Form("%s_2",hist->GetName()));

  ReadHist(hsigma,consts);
  ReadHist(hmean,means);
  ReadHist(hsigma,sigmas);


  // retrive min,max and bin width from histogram
  
  TAxis* XAxis = hist->GetXaxis();
  Low = XAxis->GetXmin();
  High = XAxis->GetXmax();
  BinWid = XAxis->GetBinWidth(0); // assumption here that bin width is constant
     
}



void NoiseCut::ReadHist(TH1D* hist, std::vector<double>& vec){
  
  for(Int_t j = 0; j < hist->GetNbinsX(); j++){
    vec.push_back(hist->GetBinContent(j));
  }

  
  
}



// check if members of vectors pass cuts in y, based on gaussian fitting of slices

void NoiseCut::PassNoiseCut(std::vector<double>& InVectX, std::vector<double>& InVectY){

  if(InVectX.size() != InVectY.size()){
    cout << "Means size = " << means.size() << ", InVectX size = " << InVectX.size() << ", InVectY size = " << InVectY.size() << endl;
    std::cout << "Two vectors given are either not the same size!" << std::endl;
    return;
  }

  
  for(Int_t i = 0; i < InVectX.size(); i++){
    
    double x = InVectX[i];
    double y = InVectY[i];

    Int_t Bin_no = ((x*1e9)-Low)/BinWid;

    Double_t mean = means[Bin_no];
    
    Double_t sigma = sigmas[Bin_no];


    
    if(y > mean - XSigmaDown*sigma && y < mean + XSigmaUp*sigma){

      // event passes cut
    }
    else{
      // event fails cut
    
      InVectX.erase(InVectX.begin() + i);
      InVectY.erase(InVectY.begin() + i);
      i--; // vector length is reduced by one so to get to 'next' element of original vector we need to reduce iterator by one

    }      
  }
  
     
}


bool NoiseCut::PassNoiseCut(double x, double y){


  Int_t Bin_no = ((x*1e9)-Low)/BinWid;

  Double_t mean = means[Bin_no];
    
  Double_t sigma = sigmas[Bin_no];


  bool Pass = false;
    
  if(y > mean - XSigmaDown*sigma && y < mean + XSigmaUp*sigma){

    Pass = true;
    // event passes cut
  }
  else{
    // event fails cut
    
    Pass = false;
  }      
  
  return Pass;
  
}



