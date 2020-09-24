////////////////////////////////////////////////////
//          apex_ME_calc
//
//   Script designed to take in ROOT Tree and DB and 
//   produce output tree with focal plane, target 
//   optics variables, PID variables and calculate
//   all ME values for all events.
//
//   1) Script should take in intial TTree and DB
//      and save new TTree with ME elements in 
//      directory related to DB name
//
//   2) Script should take result of 1) and allow
//      modification of matrix element co-effecient
//      and recalculation of ME for all events 
//      and re-save TTree with new parameters
///////////////////////////////////////////////////


#include "TString.h"
#include "Load_more_rootfiles.C"
#include "InputAPEXL.h"
#include "APEX_Sieve.h"
#include "file_def.h"

//#include "../file_def.h"
#include "TBenchmark.h"
#include "TCut.h"

#include <iostream>
#include <fstream>
#include <vector>



#include "LoadDataBase.C"




void apex_ME_calc(Int_t runnumber, TString db_name, Int_t rast = 1){
  



  TRotation fTCSInHCS;

  TVector3 TCSX(0,-1,0);
  TVector3 TCSZ(TMath::Sin(HRSAngle),0,TMath::Cos(HRSAngle));
  TVector3 TCSY = TCSZ.Cross(TCSX);
  fTCSInHCS.RotateAxes(TCSX,TCSY,TCSZ);


  // create a directory in which to save the results

  gSystem->Exec("mkdir rootfiles/" + db_name); 



  // read intitial TTree

  TChain* T = Load_more_rootfiles(runnumber);
  

  Int_t NEntries = T->GetEntries();

  

  T->SetBranchStatus("*",0);

  T->SetBranchStatus("L.tr.n",1);
  T->SetBranchStatus("L.tr.x",1);
  T->SetBranchStatus("L.tr.y",1);
  T->SetBranchStatus("L.tr.th",1);
  T->SetBranchStatus("L.tr.ph",1);
  T->SetBranchStatus("Lrb.x",1);
  T->SetBranchStatus("Lrb.y",1);
  T->SetBranchStatus("Lurb.x",1);
  T->SetBranchStatus("Lurb.y",1);
  T->SetBranchStatus("Lurb.dir.x",1);
  T->SetBranchStatus("Lurb.dir.y",1);
  T->SetBranchStatus("Lurb.dir.z",1);
  T->SetBranchStatus("Lrb.dir.x",1);
  T->SetBranchStatus("Lrb.dir.y",1);
  T->SetBranchStatus("Lrb.dir.z",1);
  T->SetBranchStatus("L.tr.tg_th",1);
  T->SetBranchStatus("L.tr.tg_ph",1);
  T->SetBranchStatus("L.tr.tg_y",1);
  T->SetBranchStatus("L.tr.tg_dp",1);


  T->SetBranchStatus("L.tr.r_x",1);
  T->SetBranchStatus("L.tr.r_y",1);
  T->SetBranchStatus("L.tr.r_th",1);
  T->SetBranchStatus("L.tr.r_ph",1);

  T->SetBranchStatus("L.tr.vz",1);

  //  T->SetBranchStatus("L.gold.*",1);
  


  // branches necessary for cuts on track quality  

  T->SetBranchStatus("L.tr.chi2",1);
  T->SetBranchStatus("L.gold.th",1);
  T->SetBranchStatus("L.gold.ph",1);
  T->SetBranchStatus("L.gold.dp",1);
  T->SetBranchStatus("L.vdc.u1.nclust",1);
  T->SetBranchStatus("L.vdc.v1.nclust",1);
  T->SetBranchStatus("L.vdc.u2.nclust",1);
  T->SetBranchStatus("L.vdc.v2.nclust",1);
  
  // branches necessary for PID cuts

  T->SetBranchStatus("L.prl1.e",1);
  T->SetBranchStatus("L.prl2.e",1);
  T->SetBranchStatus("L.gold.p",1);

  T->SetBranchStatus("L.cer.asum_c",1);
  
  
  // branches for trigger information
  T->SetBranchStatus("DL.evtype",1);
  T->SetBranchStatus("DR.evtype",1);


  // further branches for raster information
  T->SetBranchStatus("Lrb.*",1);
  

  // set-up output TTree


  //  auto newfile = TFile::Open("rootfiles/" + Dir_name + Form("/%d_replay.root",runnumber),"recreate");
  TFile* newfile = TFile::Open("rootfiles/" + db_name + Form("/%d_replay.root",runnumber),"recreate");
  

  // cut added due to problems relating to large focal plane or vertex entries
  TCut FP_cut  = "abs(L.tr.r_x)<1 && abs(L.tr.r_y)<1 && abs(L.tr.r_th)<2 && abs(L.tr.r_ph)<2 && abs(L.tr.vz)<1e+6 && abs(L.tr.tg_th)<10 && abs(L.tr.tg_ph)<10 && abs(L.tr.tg_y)<10 && abs(L.tr.tg_dp)<10 && abs(Lurb.x)<10 && abs(Lurb.y)<10 && abs(L.tr.x)<10 && abs(L.tr.y)<10 && abs(L.tr.th)<10 && abs(L.tr.ph)<10 ";

  



  //  TString final_cut = (TString) PID_cuts + " && " + (TString) GeneralSieveCut + " && " + (TString) FP_cut;

  // for the moment: do not implement cut but include tree variables such that cuts can be performed later
  TString final_cut = "";




  TTree*  out_T = T->CloneTree(0,final_cut);

  
  LoadDataBase("DB/" + db_name);



  MatrixElemExist(fTMatrixElems, TMatrix_vals);
  MatrixElemExist(fYMatrixElems, YMatrix_vals);
  MatrixElemExist(fPMatrixElems, PMatrix_vals);
  MatrixElemExist(fDMatrixElems, DMatrix_vals);



  cout << "TMatrix size = " << TMatrix_vals.values.size() << endl;
  cout << "YMatrix size = " << YMatrix_vals.values.size() << endl;
  cout << "PMatrix size = " << PMatrix_vals.values.size() << endl;
  cout << "DMatrix size = " << DMatrix_vals.values.size() << endl;





  //  out_T->Branch(Form("Test_2_%i",no_exist),fTMatrixElems,);
  out_T->Branch("T_elems",&TMatrix_vals.values);
  out_T->Branch("Y_elems",&YMatrix_vals.values);
  out_T->Branch("P_elems",&PMatrix_vals.values);
  out_T->Branch("D_elems",&DMatrix_vals.values);
  

  

  // add z-vertex and sieve positions to output tree

  // assigned foilz (makes sense for vertical foil runs)
  Double_t foilz = 0;
  // calculated foilz (based on Y elements/ tg_y reconstruction)
  Double_t reactz = 0;

  

  out_T->Branch("foilz",&foilz);
  out_T->Branch("reactz",&reactz);

  Double_t x_sieve = 0;

  Double_t y_sieve = 0;
  
  Double_t y_tgt_alt = 0;
  Double_t y_sieve_alt = 0;
  // y_sieve alt comes from using beam info to get y_tgt rather than matrix (could be useful when y_tg matrix elements are thought to be poorly calibrated)

  
  out_T->Branch("x_sieve",&x_sieve);
  out_T->Branch("y_sieve",&y_sieve);

  // create branch for 'y_tgt_alt' this is for the version of y_tgt used for optimisation of th_tgt
  out_T->Branch("y_tgt_alt",&y_tgt_alt);
  out_T->Branch("y_sieve_alt",&y_sieve_alt);



  Double_t x_tgt = 0;

  out_T->Branch("x_tgt",&x_tgt);





  // set-up focal plan variables needed to calculate matrix elements
  Double_t x_fp[100];
  Double_t y_fp[100];
  Double_t th_fp[100];
  Double_t ph_fp[100];
  Double_t th_tg[100];

  Double_t Lrb_x = 0;
  Double_t Lrb_y = 0;
  Double_t Lrb_dir_x = 0;
  Double_t Lrb_dir_y = 0;

  
  Double_t Lurb_x = 0;
  Double_t Lurb_y = 0;
  Double_t Lurb_dir_x = 0;
  Double_t Lurb_dir_y = 0;

 
  Double_t th_tgt = 0;
  Double_t y_tgt = 0;
  Double_t ph_tgt = 0;
  Double_t dp_tgt = 0;




  T->SetBranchAddress("L.tr.r_x",x_fp);
  T->SetBranchAddress("L.tr.r_y",y_fp);
  T->SetBranchAddress("L.tr.r_th",th_fp);
  T->SetBranchAddress("L.tr.r_ph",ph_fp);

  T->SetBranchAddress("L.tr.tg_th",th_tg);

  T->SetBranchAddress("Lrb.x",&Lrb_x);
  T->SetBranchAddress("Lrb.y",&Lrb_y);
  T->SetBranchAddress("Lrb.dir.x",&Lrb_dir_x);
  T->SetBranchAddress("Lrb.dir.y",&Lrb_dir_y);
  T->SetBranchAddress("Lurb.x",&Lurb_x);
  T->SetBranchAddress("Lurb.y",&Lurb_y);
  T->SetBranchAddress("Lurb.dir.x",&Lurb_dir_x);
  T->SetBranchAddress("Lurb.dir.y",&Lurb_dir_y);
  
  out_T->Branch("th_tgt",&th_tgt);
  out_T->Branch("y_tgt",&y_tgt);
  out_T->Branch("ph_tgt",&ph_tgt);
  out_T->Branch("dp_tgt",&dp_tgt);


  // set-up variables to Read-out ME 

  Double_t r_values[5][5][5][5] = {{{{0}}}};



  // Added to keep track of portion of times that the number of bytes read from the original tree exceed s a certain level
  Double_t tree_read_tracker = 0;
  Int_t read_lim = 600;
  Int_t bytes = 0;

  Double_t theta_track = 0;



  // set number of powers to calculate 
  // this sets for each event what powers of x_fp,th_fp,y_fp,ph_fp to calculate  
  
  Int_t max_pow = 5;

  
  Int_t loop_no = 0;
  
  for(Int_t i = 0; i<NEntries; i++){
    loop_no++;


    reactz = 0;
    
    bytes  = T->GetEntry(i);

    
    // if(bytes>read_lim){
    //   tree_read_tracker++;
    //   continue;
    // }


    
    //  calculate the powers we need
    for(int j=0; j<=max_pow; j++) {
      powers[j][0] = pow(x_fp[0], j);
      powers[j][1] = pow(th_fp[0], j);
      powers[j][2] = pow(y_fp[0], j);
      powers[j][3] = pow(ph_fp[0], j);
      // powers[j][4] = pow(TMath::Abs(th_fp),j);
    }


    
    CalcMatrixElem(TMatrix_vals);
    CalcMatrixElem(YMatrix_vals);
    CalcMatrixElem(PMatrix_vals);
    CalcMatrixElem(DMatrix_vals);

    
    th_tgt = TMatrix_vals.tg_v;
    y_tgt = YMatrix_vals.tg_v;
    ph_tgt = PMatrix_vals.tg_v;
    dp_tgt = DMatrix_vals.tg_v;

    
    const Int_t a = (HRSAngle > 0) ? 1 : -1;


    TVector3 BeamSpotHCS(0,0,0);

    Int_t FoilID = GetFoilID(runnumber);
    foilz = targetfoils[FoilID];

    
    if( rast!= 0){


      //      reactz = - ( y_tgt -a*MissPointZ)*TMath::Cos(ph_tgt)/TMath::Sin(HRSAngle + TMath::ATan(ph_tgt)) + Lrb_x*TMath::Cos(HRSAngle + ph_tgt)/TMath::Sin(HRSAngle + TMath::ATan(ph_tgt));

      reactz = - ( y_tgt -a*MissPointZ)*TMath::Cos(ph_tgt)/TMath::Sin(HRSAngle + ph_tgt) + Lrb_x*TMath::Cos(HRSAngle + ph_tgt)/TMath::Sin(HRSAngle + ph_tgt);
      
      
      if(IsMultiFoil(runnumber)){
	// for multiple foil (horizontal runs) runs use reactz (calculated z) for beamspot
	BeamSpotHCS.SetXYZ(Lrb_x + (Lrb_dir_x)*(reactz/BeamZDir_average),Lrb_y + (Lrb_dir_y)*(reactz/BeamZDir_average),reactz);
      }
      else if(!IsMultiFoil(runnumber)){
	// for single foil runs use foilz (kno	g_th_rc[row_count][col_v]->SetMarkerColor(col_colour[(int)x_row_column[row_count][col_v]-first_col]);wn z for said foil) for beamspot
	BeamSpotHCS.SetXYZ(Lrb_x + (Lrb_dir_x)*(foilz/BeamZDir_average),Lrb_y + (Lrb_dir_y)*(foilz/BeamZDir_average),foilz);;
      }
      
    }
    else{
      
      reactz = - ( y_tgt -a*MissPointZ)*TMath::Cos(TMath::ATan(ph_tgt))/TMath::Sin(HRSAngle + TMath::ATan(ph_tgt)) + Lurb_x*TMath::Cos(HRSAngle+TMath::ATan(ph_tgt))/TMath::Sin(HRSAngle+TMath::ATan(ph_tgt));

      BeamSpotHCS.SetXYZ(Lurb_x + (Lurb_dir_x)*(reactz/BeamZDir_average),Lurb_y + (Lurb_dir_y)*(reactz/BeamZDir_average),reactz);


      
    }



    TVector3 BeamSpotTCS=fTCSInHCS.Inverse()*(BeamSpotHCS-fPointingOffset);

    x_tgt = BeamSpotTCS.X() - BeamSpotTCS.Z() * TMath::Tan(th_tgt);



    x_sieve = x_tgt + (TMath::Tan(th_tgt) + x_tgt*ExtTarCor_ThetaCorr) * (ZPos); 
   

    y_tgt_alt = BeamSpotTCS.Y() - BeamSpotTCS.Z() * TMath::Tan(ph_tgt);    

    
    y_sieve = y_tgt + TMath::Tan(ph_tgt) * (ZPos);

    y_sieve_alt = y_tgt_alt + TMath::Tan(ph_tgt) * (ZPos);

   


    if(TMath::Abs(th_tgt) < 1e+8 &&  TMath::Abs(y_tgt) < 1e+8 && TMath::Abs(ph_tgt) < 1e+8 && TMath::Abs(dp_tgt) < 1e+8 && theta_track != th_tgt){

      
      
      out_T->Fill();



    }
    else{
      
    }


    
    
    theta_track = th_tgt;
    
  }


  cout << "final number of loop: " << loop_no << endl;
  cout << "Number of events: " << NEntries << endl << endl;


  out_T->Write();



  // extract DB elements for T,Y,P and D and save these to the rootfile

  //firstly create vectors to save DB coeffecients and corresponding powers of T,Y,P and D

  vector< Int_t > TDB_pow;
  vector< Double_t> TDB_co;

  vector< Int_t > YDB_pow;
  vector< Double_t> YDB_co;

  vector< Int_t > PDB_pow;
  vector< Double_t> PDB_co;

  vector< Int_t > DDB_pow;
  vector< Double_t> DDB_co;


  Read_Mat_DB(TMatrix_vals,TDB_pow,TDB_co);
  Read_Mat_DB(YMatrix_vals,YDB_pow,YDB_co);
  Read_Mat_DB(PMatrix_vals,PDB_pow,PDB_co);
  Read_Mat_DB(DMatrix_vals,DDB_pow,DDB_co);
  


  
  cout << "reached before DB set-up" << endl;

  // write to new DB file

  std::ofstream DB_file; 

  DB_file.open("rootfiles/" + db_name + Form("/DB_%i.dat",runnumber),std::ofstream::out);



  cout << "reached after DB set-up" << endl;
  
  
  DB_file << "Run_number = " << runnumber << endl;
  DB_file << "X Y Theta Phi Coeff0 Coeff1 ..." << endl;

  DB_file << endl << "TElements:" << endl;
  write_Mat_DB(TMatrix_vals, DB_file);

  DB_file << endl << "YElements:" << endl;
  write_Mat_DB(YMatrix_vals, DB_file);

  DB_file << endl << "PElements:" << endl;
  write_Mat_DB(PMatrix_vals, DB_file);

  DB_file << endl << "DElements:" << endl;
  write_Mat_DB(DMatrix_vals, DB_file);


  cout << "finished writing to DB" << endl;

  DB_file.close();

  

  cout << "closed DB" << endl;


  delete out_T;
  delete T;
  newfile->Close();
  cout << "closed newfile" << endl;





  
  

}

