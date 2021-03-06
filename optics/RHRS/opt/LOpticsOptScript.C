#include <iostream>
#include <cassert>

#include "TROOT.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TString.h"
#include "TVirtualFitter.h"

//#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"

#include <typeinfo>

#define th_ph_optimize true
#define draw_plots true
#define y_optimize true
#define dp_optimize false

// this parameter can turn on cutting of events not in main gaus for reconstructed phi_tg
// uses initial matrix to produce phi_tg, fit gaus to y_sieve diff plot and cut based on this
#define tail_cutting false

//#include "LOpticsOpt.h"
//#include "SaveCanvas.C"

using namespace std;

class LOpticsOpt;

LOpticsOpt * opt;
UInt_t NPara = 0;
Double_t OldMatrixArray[10000] = {-99}; //NPara
Bool_t freed[10000] = {kFALSE}; //NPara

UInt_t MaxDataPerGroup = 100;
//UInt_t MaxDataPerGroup = 100;

TString DataSource = "../Sieve/V_wires/xfp_-10_10/Sieve.full.fV_wires";


// store results of optimisations

Double_t good_theta = 1e10;
Double_t good_phi= 1e10;
Double_t good_y= 1e10;

// Inputs for minimiser and algorithm used
char* minimiser = NULL;
char* algorithm = NULL;

// default minimiser and algorithm used
const char *min_def = "Minuit";
const char *al_def = "Migrad";


// set max calls for Minuit minimisers and GSL minimisers

const Int_t Minuit_maxCalls = 1000000;
const Int_t GSL_maxCalls = 100;



// directory where rms values are saved
TString RMS_dir = "rms_csv/";


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

void DoMinTP(TString select,TString SourceDataBase, TString DestDataBase, UInt_t MaxDataPerGroup = 200, Int_t save = 0)
{
    // minimize with root

    assert(opt);
    assert(opt->fCurrentMatrixElems);
    
    opt->LoadDataBase(SourceDataBase);
    NPara = opt->Matrix2Array(OldMatrixArray, freed);
    
    // opt->LoadRawData(DataSource, (UInt_t) - 1, MaxDataPerGroup);

    cout << "DataSource = " << DataSource << endl;
    //    UInt_t NRead = opt->LoadRawData(DataSource);
    UInt_t NRead = opt->LoadRawData(DataSource, (UInt_t) - 1, MaxDataPerGroup);
    cout << "Events read = " << NRead << endl;
    

    opt->PrepareSieve();


#if tail_cutting

    opt->Sieve_hole_diff_tail(NFoils,select);
    opt->Ignore_tail();


#else
#if !tail_cutting
    
    opt->Sieve_hole_diff(NFoils,select);
    
#endif
#endif

       

    opt->Print("");
    
#if th_ph_optimize  
                                          
    // TVirtualFitter::SetDefaultFitter(); //default is Minuit
    // TVirtualFitter *fitter = TVirtualFitter::Fitter(NULL, NPara);
    // fitter->SetFCN(myfcn);


    //    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

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


      
      min->SetLimitedVariable(i,  Form("TMatrix%03d", i), OldMatrixArray[i], absold > 0 ? absold / 10 : 0.1, -abslimit, abslimit);

      if (!freed[i])  min->FixVariable(i);
    }

    
    min->PrintLevel();

    assert(opt->fNRawData > 0);
    assert(NPara > 0);
    // assert(fitter->GetNumberFreeParameters() > 0);
    // assert(fitter->GetNumberTotalParameters() == NPara);

    // Double_t arglist[1] = {0};
    // fitter->ExecuteCommand("MIGRAD", arglist, 0);
    cout << "Reached this stage" << endl;
    min->Minimize();
#endif                   
                 
    opt->Print("");
    //    opt->SaveDataBase(SourceDataBase);
    opt->SaveDataBase(DestDataBase); 
    
    Double_t rms_th = opt->SumSquareDTh();


    //Double_t rms_ph = opt->SumSquareDPhi();
    std::pair<Double_t,Double_t> rms_ph_dy = opt->SumSquareDPhi();
    Double_t rms_ph = rms_ph_dy.first;
    Double_t rms_dy = rms_ph_dy.second;


    //#if draw_plots
    //opt->check_fit_qual_Th();
    //#endif
    
    TCanvas * c1_check = opt->CheckSieve(NFoils);
    //c1_check->Print(DestDataBase+".Sieve.Opt.png", "png");
    //c1_check->Print(DestDataBase+".Sieve.Opt.eps", "eps");


    TCanvas * c1;
    
    
#if tail_cutting

    c1 = opt->Sieve_hole_diff_tail(NFoils,select);

#else
#if !tail_cutting
    c1 = opt->Sieve_hole_diff(NFoils,select);
    
#endif
#endif

    //    TCanvas * c2_diff = opt->Sieve_hole_diff(NFoils);

    if(select == "phi") c1->SaveAs(DestDataBase + ".phi.png");
    if(select == "theta") c1->SaveAs(DestDataBase + ".theta.png");

    //    TCanvas * c2 = opt->CheckSieveAccu(-1);
    //    c2->Print(DestDataBase + ".TpAccu.Opt.png", "png");
    //    c2->Print(DestDataBase + ".TpAccu.Opt.eps", "eps");

#if th_ph_optimize
    delete min;  
#endif    


    if(save){

      TString csv_name = RMS_dir + "algorithm_results.csv";

      gSystem->Exec("cp -vf " + csv_name + " " + csv_name + ".old");

      // open old csv file
      ifstream rmscsvold;
      rmscsvold.open(csv_name + ".old");  
    
      // open new csv file
      ofstream rmscsv;
      rmscsv.open(csv_name);

      rmscsv<<fixed<<setprecision(10);
      cout<<fixed<<setprecision(10);

      rmscsv << "Minimiser,Algorithm,Optimisation Variable,Rms (real xyz),Data_name" << endl;;
      
      string line; // used to gather csv file

      getline(rmscsvold,line); // these lines are for first two line of csv file (explenatory)

      
      // copy previous results to new file
      
      while(!rmscsvold.eof()){
	getline(rmscsvold,line);
	rmscsv<<line<<endl;
      }

      
      TString Data_name = DataSource;

      Data_name.Remove(0,24);

      cout << "Data_name = " << Data_name << endl;
       
      if( myfcnX == myfcn1){
	cout << "Theta!!!" << endl;
	rmscsv<<minimiser<<","<<algorithm<<",Theta,"<<rms_th<<","<<Data_name;
      }
      else if( myfcnX == myfcn2){
	cout << "Phi!!!" << endl;
	rmscsv<<minimiser<<","<<algorithm<<",Phi,"<<rms_ph<<","<<Data_name;
      }
      else{
	cout << "Not Theta or Phi!!!" << endl;
      }
      
      

    }
    
    
    good_theta = rms_th;
    good_phi = rms_ph;
    good_y = rms_dy;
    
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << ">>>>>>>>>>>> DoMinTP Finished! <<<<<<<<<<<<<<<<<<" << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;


}

void DoMinY(TString SourceDataBase, TString DestDataBase, UInt_t MaxDataPerGroup = 200, Int_t save = 0)
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


	min->SetLimitedVariable(i,  Form("TMatrix%03d", i), OldMatrixArray[i], absold > 0 ? absold / 10 : 0.1, -abslimit, abslimit);

	if (!freed[i])  min->FixVariable(i);
    }


    //    min->PrintLevel();
    

    assert(opt->fNRawData > 0);
    assert(NPara > 0);
        
    min->Minimize();
#endif
                               
    opt->Print("");
    opt->SaveDataBase(DestDataBase);

    Double_t rms_y = opt->SumSquareDTgY();

    TCanvas * c1 = opt->CheckVertex();
    c1->Print(DestDataBase + ".Vertex.Opt.png", "png");

#if y_optimize
    delete min;
#endif


    if(save){

      TString csv_name = RMS_dir + "algorithm_results.csv";

      gSystem->Exec("cp -vf " + csv_name + " " + csv_name + ".old");

      // open old csv file
      ifstream rmscsvold;
      rmscsvold.open(csv_name + ".old");  
    
      // open new csv file
      ofstream rmscsv;
      rmscsv.open(csv_name);

      rmscsv<<fixed<<setprecision(10);
      cout<<fixed<<setprecision(10);

      rmscsv << "Minimiser, Algorithm, Optimisation Variable, Rms (real x,y,z), Data_name" << endl;;
      
      string line; // used to gather csv file

      getline(rmscsvold,line); // these lines are for first two line of csv file (explenatory)

      
      // copy previous results to new file
      
      while(!rmscsvold.eof()){
	getline(rmscsvold,line);
	rmscsv<<line<<endl;
      }
      
      TString Data_name = DataSource;

      Data_name.Remove(0,24);
       
      if( myfcnX == myfcn3){
	cout << "Y!!!" << endl;
	rmscsv<<minimiser<<","<<algorithm<<",Y,"<<rms_y<<","<<Data_name;
      }
      else{
	cout << "Not Y!!!" << endl;
      }
            
    }        

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
    // opt->PrepareDp();

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

    // TCanvas * c1 = opt->CheckDp();
    // c1->Print(DestDataBase + ".Dp.Opt.png", "png");
#if dp_optimze
    delete fitter;
#endif
}

void PlotDataBase(TString DatabaseFileName, UInt_t MaxDataPerGroup = 1000)
{
    opt = new LOpticsOpt();

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

void LOpticsOptScript(TString run, TString range,TString select, TString SourceDataBase, TString DestDataBase)
{
    opt = new LOpticsOpt();

   DataSource = "../Sieve/"+run+"/xfp_"+range+"/Sieve.full.f"+run;
  //DataSource = "../Sieve/"+run+"/xfp_"+range+"/Sieve.full.f4650";
  
  //TString extra_dir = run + "/";
  TString extra_dir = "";
  if(select != "phi") extra_dir = run + "/";
  //if(select != "phi" && select != "y") extra_dir = run + "/";
    
    
    SourceDataBase = "DB/" + extra_dir + SourceDataBase;
    DestDataBase = "DB/" + run + "/" + DestDataBase;


    const char* min = min_def;
    const char* algo = al_def;
    Int_t save = 0;
    
    minimiser = (char*)min;
    algorithm = (char*)algo;

    // minimiser = strdup(min);
    // algorithm = strdup(algo);

    
    cout << endl;
    cout << "Using minimiser " << minimiser << endl;
    cout << "Using algorithm " << algorithm << endl;
    cout << endl;


    Int_t s = 0;
    if (select == "theta") s = 1;
    if (select == "phi") s = 2;
    if (select == "y") s = 3;
    if (select == "delta") s = 4;
    if (select == "pta") s = 5;
    if (select == "yta") s = 6;

    gStyle->SetOptStat(0);

    switch (s) {
    case 1:
        cout << "Optimizing for Theta\n";
	//       myfcn = myfcn1;
	myfcnX = myfcn1;
        opt->fCurrentMatrixElems = &(opt->fTMatrixElems);
	
        DoMinTP(select,SourceDataBase, DestDataBase, 500, save);
        break;
    case 2:
        cout << "Optimizing for Phi\n";
        myfcnX = myfcn2;
        opt->fCurrentMatrixElems = &(opt->fPMatrixElems);
        DoMinTP(select,SourceDataBase, DestDataBase, 100, save);
        break;
    case 3:
        cout << "Optimizing for Y\n";
        myfcnX = myfcn3;
        opt->fCurrentMatrixElems = &(opt->fYMatrixElems);
        //DoMinY(SourceDataBase, DestDataBase, 200000, save);
	DoMinY(SourceDataBase, DestDataBase, 100, save);
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
        DoMinTP(select,SourceDataBase, DestDataBase, 500);
        break;
   case 6:
        cout << "Optimizing for y\n";
        myfcnX = myfcn3;
        opt->fCurrentMatrixElems = &(opt->fYTAMatrixElems);
        DoMinY(SourceDataBase, DestDataBase, 200000);
        break;
    default:
        break;
    }

    gSystem->Exec(Form("cp -vf %s %s.source", SourceDataBase.Data(), DestDataBase.Data()));
    //    gSystem->Exec(Form("cp -vf log %s.log", DestDataBase.Data()));

    delete opt;

    return;
}


// script that processes inout ASCII scripts with different column numbers for the same run and saves optimisation outputs (rms) to text file
Int_t LOpticsOptScript_alter_ASCII(TString select, TString SourceDataBase, TString DestDataBase, const char* DataSourceName){

  /*
  DataSource = "../Tree2Ascii/text_cuts/" +  TString(DataSourceName);

  TString DataSource_orig = DataSource;

  //DataSource += "_V1_V2";
  DataSource += "_V1";
  
  // test if (base) text_cut file exists

  // if it does create text file to output results to
  if(!(gSystem->AccessPathName(DataSource))){
    cout << DataSource << " file was found" << endl;
  }
  else if(gSystem->AccessPathName(DataSource)){
    cout << DataSource << " file was not found!!!" << endl;
    return 0;
  }


  // create output file to store
  ofstream ofile;
  ofile.open("col_vary/" + TString(DataSourceName)); //col_vary

  ofile << "Phi_rms" << "," << "Y_rms" << "," << "Col_change" << endl;
  
  //  carry out optimisation with unaltered ascii file (original column and row numbers)
  
  LOpticsOptScript( select,SourceDataBase,DestDataBase);
  LOpticsOptScript( select,DestDataBase,DestDataBase);


  ofile << good_phi << "," << good_y << "," << 0 << endl;
  
  
  // loop through different text files

  //  TString DataSource_orig = DataSource;
  TString DestDataBase_orig = DestDataBase;
  
  for(Int_t i = -10; i<11; i++){

    
    //  DataSource = DataSource_orig + Form("_%i_col_change_V1_V2",i);
    DataSource = DataSource_orig + Form("_%i_col_change_V1",i);
    // DataSource = DataSource_orig + Form("_%i_col_change",i);

    if(!(gSystem->AccessPathName(DataSource))){
      cout << DataSource << " file was found" << endl;

      DestDataBase = DestDataBase_orig +  Form("_%i_col_change",i);
      LOpticsOptScript( select,SourceDataBase,DestDataBase);
      LOpticsOptScript( select,DestDataBase,DestDataBase);

      ofile << good_phi << "," << good_y << "," << i << endl;
    }



  }
*/
  return 1;
  
}
