#include <iostream>
#include <cassert>

#include "TROOT.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TString.h"
#include "TVirtualFitter.h"

#define th_ph_optimize true
#define draw_plots true
#define y_optimize true
#define dp_optimize false

// this parameter can turn on cutting of events not in main gaus for reconstructed phi_tg
// uses initial matrix to produce phi_tg, fit gaus to y_sieve diff plot and cut based on this
#define tail_cutting true

//#include "ROpticsOpt.h"
//#include "SaveCanvas.C"

using namespace std;

class ROpticsOpt;

ROpticsOpt * opt;
UInt_t NPara = 0;
Double_t OldMatrixArray[10000] = {-99}; //NPara
Bool_t freed[10000] = {kFALSE}; //NPara

UInt_t MaxDataPerGroup = 100;
//UInt_t MaxDataPerGroup = 100;

// Inputs for minimiser and algorithm used
char* minimiser = NULL;
char* algorithm = NULL;

// default minimiser and algorithm used
const char *min_def = "Minuit";
const char *al_def = "Migrad";


// set max calls for Minuit minimisers and GSL minimisers

const Int_t Minuit_maxCalls = 1000000;
const Int_t GSL_maxCalls = 100;



TString run;     //Don't forget to change CheckSieve(#) as well
TString range;
TString DataSource;


typedef void (*PTRFCN)(Int_t &, Double_t *, Double_t &, Double_t*, Int_t);
PTRFCN myfcn = NULL;

typedef double (*PTRFCNX)(const double *);
PTRFCNX myfcnX = NULL;

double myfcn1(const double *par)
{
  //compute the sum of squares of dth
 

  static UInt_t NCall = 0;
  NCall++;


  assert(opt);
  assert(opt->fCurrentMatrixElems);

   

  opt->Array2Matrix(par);
  double f = opt->SumSquareDTh();
  
  //  cout << "For NCall :" << NCall << ", f = " << f << endl;

  return f;
}

double myfcn2(const double *par)
{
  //compute the sum of squares of dph
  static UInt_t NCall = 0;
  NCall++;

  assert(opt);
  assert(opt->fCurrentMatrixElems);
  
  opt->Array2Matrix(par);
  
  std::pair<Double_t,Double_t> rms_ph_dy = opt->SumSquareDPhi();
  //  Double_t f = rms_ph_dy.first;
  Double_t f = rms_ph_dy.second;
  
  //double f = opt->SumSquareDPhi();
  

  //  cout << "For NCall :" << NCall << ", f = " << f << endl;

  return f;
}

double myfcn3(const double *par)
{
  //compute the sum of squares of dtgy
  static UInt_t NCall = 0;
  NCall++;

  assert(opt);
  assert(opt->fCurrentMatrixElems);

  opt->Array2Matrix(par);
  double f = opt->SumSquareDTgY();
  
  //  cout << "For NCall :" << NCall << ", f = " << f << endl;

  return f;
}

void myfcn4(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t)
{
    //compute the sum of squares of dp

    assert(opt);
    assert(opt->fCurrentMatrixElems);

    opt->Array2Matrix(par);
    f = opt->SumSquareDp();

    return;
}

void DoMinTP(TString SourceDataBase, TString DestDataBase, UInt_t MaxDataPerGroup = 200)
{
    // minimize with root

    assert(opt);
    assert(opt->fCurrentMatrixElems);

    opt->LoadDataBase(SourceDataBase);
    NPara = opt->Matrix2Array(OldMatrixArray, freed);
    
    opt->LoadRawData(DataSource, (UInt_t) - 1, MaxDataPerGroup);
    opt->PrepareSieve();


#if tail_cutting
    
    opt->Sieve_hole_diff_tail(NFoils);
    opt->Ignore_tail();
    
    
#else
#if !tail_cutting
    
    opt->Sieve_hole_diff(NFoils);
    
#endif
#endif
    

    opt->Print("");
    
#if th_ph_optimize  

    /*
    TVirtualFitter::SetDefaultFitter(); //default is Minuit
    TVirtualFitter *fitter = TVirtualFitter::Fitter(NULL, NPara);
    fitter->SetFCN(myfcn);
    */

    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minimiser, algorithm);

    min->SetTolerance(0.001);
    
    ROOT::Math::Functor f(myfcnX,NPara); 
    
    min->SetFunction(f);

    min->SetMaxIterations(GSL_maxCalls);  // for GSL
    min->SetMaxFunctionCalls(Minuit_maxCalls);   // for Minuit/Minuit2
    
    for (UInt_t i = 0; i < NPara; i++) {
      //      cout<<"i:"<<i<<endl;
        Double_t absold = TMath::Abs(OldMatrixArray[i]);
        Double_t abslimit = absold > 0 ? absold * 10000 : 10000;

        //fitter->SetParameter(i, Form("TMatrix%03d", i), OldMatrixArray[i], absold > 0 ? absold / 10 : 0.1, -abslimit, abslimit);
        // //fitter->SetParameter(1,"asdf",0,0,0,0);

	min->SetLimitedVariable(i,  Form("TMatrix%03d", i), OldMatrixArray[i], absold > 0 ? absold / 10 : 0.1, -abslimit, abslimit);

        //if (!freed[i]) fitter->FixParameter(i);
	if (!freed[i])  min->FixVariable(i);
    }

    //fitter->Print();
    //cout << fitter->GetNumberFreeParameters() << " Free  / " << fitter->GetNumberTotalParameters() << " Parameters\n";

    min->PrintLevel();
    
    assert(opt->fNRawData > 0);
    assert(NPara > 0);
    //assert(fitter->GetNumberFreeParameters() > 0);
    //assert(fitter->GetNumberTotalParameters() == NPara);

    //Double_t arglist[1] = {0};
    //fitter->ExecuteCommand("MIGRAD", arglist, 0);
       
    min->Minimize();
#endif                   
                 
    opt->Print("");
    //    opt->SaveDataBase(SourceDataBase);
    opt->SaveDataBase(DestDataBase); 
    
    opt->SumSquareDTh();
    opt->SumSquareDPhi();

    //#if draw_plots
    //opt->check_fit_qual_Th();
    //#endif
    
    TCanvas * c1 = opt->CheckSieve(-1);
    c1->Print(DestDataBase+".Sieve.Opt.png", "png");
    //c1->Print(DestDataBase+".Sieve.Opt.eps", "eps");
    
    //    TCanvas * c2 = opt->CheckSieveAccu(-1);
    //    c2->Print(DestDataBase + ".TpAccu.Opt.png", "png");
    //    c2->Print(DestDataBase + ".TpAccu.Opt.eps", "eps");


    TCanvas * c2_diff;
    
#if tail_cutting

    c2_diff = opt->Sieve_hole_diff_tail(NFoils);

#else
#if !tail_cutting
    c2_diff = opt->Sieve_hole_diff(NFoils);
    
#endif
#endif
    
    c2_diff->SaveAs(DestDataBase + ".phi.png");

#if th_ph_optimize
    //delete fitter;
    delete min;
#endif    
}

void DoMinY(TString SourceDataBase, TString DestDataBase, UInt_t MaxDataPerGroup = 200)
{
    // minimize with root

    assert(opt);
    assert(opt->fCurrentMatrixElems);

    opt->LoadDataBase(SourceDataBase);
    NPara = opt->Matrix2Array(OldMatrixArray, freed);
    cout<<"NPara:"<<NPara<<endl;
    opt->LoadRawData(DataSource, (UInt_t) - 1, MaxDataPerGroup);
    opt->PrepareVertex();

     // added for phi_tg tail exclusion
    opt->PrepareSieve();

    opt->Sieve_hole_diff_tail(NFoils);

#if tail_cutting
    opt->Ignore_tail();
#endif
    
    opt->Print("");


#if y_optimize
    /*
    TVirtualFitter::SetDefaultFitter(); //default is Minuit
    TVirtualFitter *fitter = TVirtualFitter::Fitter(NULL, NPara);
    fitter->SetFCN(myfcn);
    */

    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minimiser, algorithm);

    min->SetTolerance(0.001);
    
    ROOT::Math::Functor f(myfcnX,NPara); 
    
    min->SetFunction(f);

    min->SetMaxIterations(GSL_maxCalls);  // for GSL
    min->SetMaxFunctionCalls(Minuit_maxCalls);   // for Minuit/Minuit2

    cout << "NPara = " << NPara << endl;
    
    for (UInt_t i = 0; i < NPara; i++) {
        Double_t absold = TMath::Abs(OldMatrixArray[i]);
        Double_t abslimit = absold > 0 ? absold * 10000 : 10000;

        //fitter->SetParameter(i, Form("TMatrix%03d", i), OldMatrixArray[i], absold > 0 ? absold / 10 : 0.1, -abslimit, abslimit);
        min->SetLimitedVariable(i,  Form("TMatrix%03d", i), OldMatrixArray[i], absold > 0 ? absold / 10 : 0.1, -abslimit, abslimit);

        //if (!freed[i]) fitter->FixParameter(i);
	if (!freed[i])  min->FixVariable(i);
    }

    min->PrintLevel();
    
    assert(opt->fNRawData > 0);
    assert(NPara > 0);
        
    min->Minimize();

    /*
    fitter->Print();
    cout << fitter->GetNumberFreeParameters() << " Free  / " << fitter->GetNumberTotalParameters() << " Parameters\n";
    //    cout<<"NPara:"<<NPara<<endl;

    assert(opt->fNRawData > 0);
    assert(NPara > 0);
    assert(fitter->GetNumberFreeParameters() > 0);
    assert(fitter->GetNumberTotalParameters() == NPara);
    
    Double_t arglist[1] = {0};
    fitter->ExecuteCommand("MIGRAD", arglist, 0);
    */
#endif
                               
    opt->Print("");
    opt->SaveDataBase(DestDataBase);

    opt->SumSquareDTgY();

    TCanvas * c1 = opt->CheckVertex();
    c1->Print(DestDataBase + ".Vertex.Opt.png", "png");

#if y_optimize
    //delete fitter;
    delete min;
#endif
}

void DoMinDp(TString SourceDataBase, TString DestDataBase, UInt_t MaxDataPerGroup = 200)
{
    // minimize with root

    assert(opt);
    assert(opt->fCurrentMatrixElems);

    cout << "Optimizing for dp\n";
    opt->fCurrentMatrixElems = &(opt->fDMatrixElems);

    opt->LoadDataBase(SourceDataBase);
    NPara = opt->Matrix2Array(OldMatrixArray, freed);
    opt->LoadRawData(DataSource, (UInt_t) - 1, MaxDataPerGroup);
    opt->PrepareDp();

    //compensate bias due to dp event selections
    // 	opt->fArbitaryDpKinShift[0] = 2.786177e-05;
    // 	opt->fArbitaryDpKinShift[1] = 8.168538e-05;
    // 	opt->fArbitaryDpKinShift[2] = 5.299596e-05;
    // 	opt->fArbitaryDpKinShift[3] = 3.175602e-05;
    // 	opt->fArbitaryDpKinShift[4] = 9.519830e-05;
    
    opt->fArbitaryDpKinShift[0] = 0.;
    opt->fArbitaryDpKinShift[1] = 0.;
    opt->fArbitaryDpKinShift[2] = 0.;
    opt->fArbitaryDpKinShift[3] = 0.;
    opt->fArbitaryDpKinShift[4] = 0.;
    
    opt->Print("");

#if dp_optimize
         
    TVirtualFitter::SetDefaultFitter(); //default is Minuit
    TVirtualFitter *fitter = TVirtualFitter::Fitter(NULL, NPara);
    fitter->SetFCN(myfcn);

    for (UInt_t i = 0; i < NPara; i++) {
        Double_t absold = TMath::Abs(OldMatrixArray[i]);
        Double_t abslimit = absold > 0 ? absold * 10000 : 10000;

        fitter->SetParameter(i, Form("TMatrix%03d", i), OldMatrixArray[i], absold > 0 ? absold / 10 : 0.1, -abslimit, abslimit);
        // fitter->SetParameter(1,"asdf",0,0,0,0);

        if (!freed[i]) fitter->FixParameter(i);
    }

    fitter->Print();
    cout << fitter->GetNumberFreeParameters() << " Free  / " << fitter->GetNumberTotalParameters() << " Parameters\n";

    assert(opt->fNRawData > 0);
    assert(NPara > 0);
    assert(fitter->GetNumberFreeParameters() > 0);
    assert(fitter->GetNumberTotalParameters() == NPara);

    Double_t arglist[1] = {0};
    fitter->ExecuteCommand("MIGRAD", arglist, 0);
#endif         
    opt->Print("");
    opt->SaveDataBase(DestDataBase);

    opt->SumSquareDp();

    TCanvas * c1 = opt->CheckDp();
    c1->Print(DestDataBase + ".Dp.Opt.png", "png");
#if dp_optimze
    delete fitter;
#endif
}

void PlotDataBase(TString DatabaseFileName, UInt_t MaxDataPerGroup = 1000)
{
    opt = new ROpticsOpt();

    assert(opt);

    gStyle->SetOptStat(0);

    opt->LoadDataBase(DatabaseFileName);
    opt->Print("");

    opt->LoadRawData(DataSource, (UInt_t) - 1, MaxDataPerGroup);

    opt->PrepareSieve();
    //opt->PrepareDp();

    //opt->SumSquareDTh(kTRUE);
    //opt->SumSquareDPhi(kTRUE);
    //opt->SumSquareDp(kTRUE);

    TCanvas * c1 = opt->CheckSieve();
    c1->Print(DatabaseFileName + ".Sieve.png", "png");
 

   //TCanvas * c1 = opt->CheckDpGlobal();
    //c1->Print(DatabaseFileName + ".Dp.png","png");

    delete opt;
}

void ROpticsOptScript(TString run, TString range,TString select, TString SourceDataBase, TString DestDataBase)
{
  
  DataSource = "../Sieve/"+run+"/Sieve.full.f"+run;
  //DataSource = "../Sieve/"+run+"/xfp_"+range+"/Sieve.full.f4647";
  
  opt = new ROpticsOpt();
  
  //TString extra_dir = run + "/";
  TString extra_dir = "V_wires_test/";
  if(select != "phi") extra_dir = run + "/";
  //if(select != "phi" && select != "y") extra_dir = run + "/";
    
    SourceDataBase = "DB/" + extra_dir + SourceDataBase;
    DestDataBase = "DB/" + run + "/" + DestDataBase;


    const char* min = min_def;
    const char* algo = al_def;

    minimiser = (char*)min;
    algorithm = (char*)algo;
    
    Int_t s = 0;
    if (select == "theta") s = 1;
    if (select == "phi") s = 2;
    if (select == "y") s = 3;
    if (select == "delta") s = 4;
    if (select == "pta") s = 5;

    gStyle->SetOptStat(0);

    switch (s) {
    case 1:
        cout << "Optimizing for Theta\n";
        myfcnX = myfcn1;
        opt->fCurrentMatrixElems = &(opt->fTMatrixElems);
        DoMinTP(SourceDataBase, DestDataBase, 500);
        break;
    case 2:
        cout << "Optimizing for Phi\n";
        myfcnX = myfcn2;
        opt->fCurrentMatrixElems = &(opt->fPMatrixElems);
        DoMinTP(SourceDataBase, DestDataBase, 100);
        break;
    case 3:
        cout << "Optimizing for Y\n";
        myfcnX = myfcn3;
        opt->fCurrentMatrixElems = &(opt->fYMatrixElems);
        //DoMinY(SourceDataBase, DestDataBase, 200000);
	DoMinY(SourceDataBase, DestDataBase, 100);
        break;
    case 4:
        cout << "Optimizing for Delta\n";
        myfcn = myfcn4;
        opt->fCurrentMatrixElems = &(opt->fDMatrixElems);
        DoMinDp(SourceDataBase, DestDataBase, 200000);
        break;
    case 5:
        cout << "Optimizing for Phi\n";
        myfcnX = myfcn2;
        opt->fCurrentMatrixElems = &(opt->fPTAMatrixElems);
        DoMinTP(SourceDataBase, DestDataBase, 500);
        break;
    default:
        break;
    }

    gSystem->Exec(Form("cp -vf %s %s.source", SourceDataBase.Data(), DestDataBase.Data()));
    //    gSystem->Exec(Form("cp -vf log %s.log", DestDataBase.Data()));

    delete opt;

    return;
}
