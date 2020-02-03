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



//void CalcMatrixElem_APEX(vector<THaMatrixElement>& matrix );


void apex_ME_calc(Int_t runnumber, TString db_name, Int_t rast = 1){
  



  TRotation fTCSInHCS;

  TVector3 TCSX(0,-1,0);
  TVector3 TCSZ(TMath::Sin(HRSAngle),0,TMath::Cos(HRSAngle));
  TVector3 TCSY = TCSZ.Cross(TCSX);
  fTCSInHCS.RotateAxes(TCSX,TCSY,TCSZ);


  // create a directory in which to save the results

  // TString Dir_name;
  
  // cout << "Enter a directory name in which to save results: " << endl;
    
  // cin >> Dir_name;

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
  T->SetBranchStatus("L.tr.tg_th",1);
  T->SetBranchStatus("L.tr.tg_ph",1);
  T->SetBranchStatus("L.tr.tg_y",1);
  T->SetBranchStatus("L.tr.tg_dp",1);


  T->SetBranchStatus("L.tr.r_x",1);
  T->SetBranchStatus("L.tr.r_y",1);
  T->SetBranchStatus("L.tr.r_th",1);
  T->SetBranchStatus("L.tr.r_ph",1);

  T->SetBranchStatus("L.tr.vz",1);
  

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
  
  




  // set-up output TTree


  //  auto newfile = TFile::Open("rootfiles/" + Dir_name + Form("/%d_replay.root",runnumber),"recreate");
  TFile* newfile = TFile::Open("rootfiles/" + db_name + Form("/%d_replay.root",runnumber),"recreate");
  

  // cut added due to problems relating to large focal plane or vertex entries
  TCut FP_cut  = "abs(L.tr.r_x)<1 && abs(L.tr.r_y)<1 && abs(L.tr.r_th)<2 && abs(L.tr.r_ph)<2 && abs(L.tr.vz)<1e+6 && abs(L.tr.tg_th)<10 && abs(L.tr.tg_ph)<10 && abs(L.tr.tg_y)<10 && abs(L.tr.tg_dp)<10 && abs(Lurb.x)<10 && abs(Lurb.y)<10 && abs(L.tr.x)<10 && abs(L.tr.y)<10 && abs(L.tr.th)<10 && abs(L.tr.ph)<10 ";


  // TCut GeneralSieveCut ="L.tr.n==1 && L.tr.chi2<0.003 && abs(L.gold.th)<0.05 && L.gold.ph>-0.07 && L.gold.ph<0.025 && abs(L.gold.dp)<0.05 && L.vdc.u1.nclust==1 && L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.v2.nclust==1";

 // TCut GeneralSieveCut ="L.tr.n==1 && L.tr.chi2<0.003 && abs(L.gold.th)<0.05 && L.gold.ph>-0.07 && L.gold.ph<0.025 && abs(L.tr.r_x)<1 && L.vdc.u1.nclust==1 && L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.v2.nclust==1";


 //  TCut PID_cuts = "(L.prl1.e/(L.gold.p*1000))>0.2 && ((L.prl1.e+L.prl2.e)/(L.gold.p*1000))>0.51 &&  L.cer.asum_c >400";

  



  TString final_cut = (TString) PID_cuts + " && " + (TString) GeneralSieveCut + " && " + (TString) FP_cut;




  TTree*  out_T = T->CloneTree(0,final_cut);


  // add additional branches 

  
  // create structure for ME (matrix elements)
  // 





  //  auto newBranch = new TBranch();

  
  LoadDataBase("DB/" + db_name);


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // test-loop to print DB out



  //  std::vector<THaMatrixElement>& matrix = fTMatrixElems;

  
  
  //  MatrixElems_Vals TMatrix_vals;

  MatrixElemExist(fTMatrixElems, TMatrix_vals);
  MatrixElemExist(fYMatrixElems, YMatrix_vals);
  MatrixElemExist(fPMatrixElems, PMatrix_vals);
  MatrixElemExist(fDMatrixElems, DMatrix_vals);



  cout << "TMatrix size = " << TMatrix_vals.values.size() << endl;
  cout << "YMatrix size = " << YMatrix_vals.values.size() << endl;
  cout << "PMatrix size = " << PMatrix_vals.values.size() << endl;
  cout << "DMatrix size = " << DMatrix_vals.values.size() << endl;



  // for( int a = 0; a<5; a++){
  //   for( int b = 0; b<5; b++){
  //     for( int c = 0; c<5; c++){
  // 	for( int d = 0; d<5; d++){
	  
  // 	  cout << "Exist[" << a << "][" << b << "][" << c << "][" << d << "] = " << TMatrix_vals.Exist[a][b][c][d] << endl;
	  

  // 	  if(TMatrix_vals.Exist[a][b][c][d] != 0){
  // 	    //	    out_T->Branch(Form("TElems_%i_%i_%i_%i",a,b,c,d),values[no_exist]);
  // 	    //	    out_T->Branch(Form("Test_2_%i",no_exist),values);

  // 	    if(no_exist > 0){

  // 	      TMatrix_vals.values.push_back(0);

  // 	    }

  // 	    no_exist++;
	    
  // 									       }
  // 	  }
  // 	}
  //     }
  //   }





  //  out_T->Branch(Form("Test_2_%i",no_exist),fTMatrixElems,);
  out_T->Branch("T_elems",&TMatrix_vals.values);
  out_T->Branch("Y_elems",&YMatrix_vals.values);
  out_T->Branch("P_elems",&PMatrix_vals.values);
  out_T->Branch("D_elems",&DMatrix_vals.values);
  


  // add z-vertex and sieve positions to output tree

  Double_t reactz = 0;
    
  out_T->Branch("reactz",&reactz);

  Double_t x_sieve = 0;

  Double_t y_sieve = 0;
  

  Double_t y_sieve_alt = 0;
  // y_sieve alt comes from using beam info to get y_tgt rather than matrix (could be useful when y_tg matrix elements are thought to be poorly calibrated)

  
  out_T->Branch("x_sieve",&x_sieve);
  out_T->Branch("y_sieve",&y_sieve);

  out_T->Branch("y_sieve_alt",&y_sieve_alt);



  Double_t x_tgt = 0;

  out_T->Branch("x_tgt",&x_tgt);

  // std::cout.precision(7);

  // Int_t it_count = 0;
  // for( vector<THaMatrixElement>::iterator it=matrix.begin(); it!=matrix.end(); it++ ) {






  //   std::cout << "for ME[" << it_count << "] with x-order = " << it->order << " : " << std::endl;
  //   for(int i=0; i<=it->order-1; i++){
  //     std::cout << "it->poly[" << i << "] = " << it->poly[i] << std::endl;
  //     std::cout << "it->pw[" << i << "] = " << it->pw[i] << std::endl;
  //   }
    
  //   std::cout << std::endl << std::endl;
  //   it_count++;

  // }
  //    ~~~~~~~~~~~~~~~~~~~~~~~~


  //  std::cout << "finished " << std::endl;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Calculate new target variables (with given DB)

//  out_T->Branch("",fTMatrixelems);



  // set-up focal plan variables needed to calculate matrix elements
  Double_t x_fp[100];
  Double_t y_fp[100];
  Double_t th_fp[100];
  Double_t ph_fp[100];
  Double_t th_tg[100];

  Double_t Lrb_x = 0;
  Double_t Lrb_y = 0;
  Double_t Lurb_x = 0;
  Double_t Lurb_y = 0;

 
  Double_t th_tgt = 0;
  Double_t y_tgt = 0;
  Double_t ph_tgt = 0;
  Double_t dp_tgt = 0;



  // T->SetBranchAddress("L.tr.x",&x_fp);
  // T->SetBranchAddress("L.tr.y",&y_fp);
  // T->SetBranchAddress("L.tr.th",&th_fp);
  // T->SetBranchAddress("L.tr.ph",&ph_fp);

  T->SetBranchAddress("L.tr.r_x",x_fp);
  T->SetBranchAddress("L.tr.r_y",y_fp);
  T->SetBranchAddress("L.tr.r_th",th_fp);
  T->SetBranchAddress("L.tr.r_ph",ph_fp);

  T->SetBranchAddress("L.tr.tg_th",th_tg);

  T->SetBranchAddress("Lrb.x",&Lrb_x);
  T->SetBranchAddress("Lrb.y",&Lrb_y);
  T->SetBranchAddress("Lurb.x",&Lurb_x);
  T->SetBranchAddress("Lurb.y",&Lurb_y);



  out_T->Branch("th_tgt",&th_tgt);
  out_T->Branch("y_tgt",&y_tgt);
  out_T->Branch("ph_tgt",&ph_tgt);
  out_T->Branch("dp_tgt",&dp_tgt);


  // set-up variables to Read-out ME 

  Double_t r_values[5][5][5][5] = {{{{0}}}};





  //  out_T->Branch("X_ME",r_values,"r_values[5][5][5][5]/D");
  //  out_T->Branch("X_ME",&fTMatrixElems); 






  // Added to keep track of portion of times that the number of bytes read from the original tree exceed s a certain level
  Double_t tree_read_tracker = 0;
  Int_t read_lim = 600;
  Int_t bytes = 0;

  Double_t theta_track = 0;

  //  NEntries = 10;



  for(Int_t i = 0; i<NEntries; i++){



    bytes  = T->GetEntry(i);

    //    cout << "Out_T first check : out_T->GetEntries() = " << out_T->GetEntries() << endl;

    
    if(bytes>read_lim){
      tree_read_tracker++;
      continue;
    }

    //    cout << "Portion of events > 600 bytes = " << tree_read_tracker/i << endl;


    //    std::cout << " loop " << i << " where x_fp = " << x_fp << std::endl;

    //  calculate the powers we need
    for(int j=0; j<5; j++) {
      powers[j][0] = pow(x_fp[0], j);
      powers[j][1] = pow(th_fp[0], j);
      powers[j][2] = pow(y_fp[0], j);
      powers[j][3] = pow(ph_fp[0], j);
      // powers[j][4] = pow(TMath::Abs(th_fp),j);
    }


    //    CalcMatrix(x_fp,fTMatrixElems);
    // cout << " theta = " << CalcTargetVar(fTMatrixElems, powers) << endl;

    // cout << " Print focal plane variables: " << endl;
    // cout << "x_fp = " << x_fp << ", th_fp = " << th_fp << ", y_fp = " << y_fp << ", ph_fp = " << ph_fp << endl;
    
    CalcMatrixElem(TMatrix_vals);
    CalcMatrixElem(YMatrix_vals);
    CalcMatrixElem(PMatrix_vals);
    CalcMatrixElem(DMatrix_vals);

    
    th_tgt = TMatrix_vals.tg_v;
    y_tgt = YMatrix_vals.tg_v;
    ph_tgt = PMatrix_vals.tg_v;
    dp_tgt = DMatrix_vals.tg_v;



    // cout << "new target theta = " << thtgt << " and old target theta = " << th_tg << endl << endl;

   //  cout << "New target values: th_tgt = " << th_tgt << ", y_tgt = " << y_tgt << ", ph_tgt = " << ph_tgt << ", dp_tgt = " << dp_tgt << endl;


   //  cout << "New matrix values: T_elems (size " << TMatrix_vals.values.size() << ") = ";
   //  for(Int_t j = 0; j < TMatrix_vals.values.size(); j++){
   //    if(TMath::Abs(TMatrix_vals.values[j]) > 1e+8){
   // 	cout << TMatrix_vals.values[j] << ", ";
   //    }
   
      
   //  }

   //  cout << " " << endl;
   //  cout << "New matrix values: Y_elems (size " << YMatrix_vals.values.size() << ") = ";
   //  for(Int_t j = 0; j < YMatrix_vals.values.size(); j++){

   //    if(TMath::Abs(YMatrix_vals.values[j]) > 1e+8){
   // 	cout << YMatrix_vals.values[j] << ", ";
   //    }
   
      
   //  }

   //  cout << " " << endl;
   //  cout << "New matrix values: P_elems (size " << PMatrix_vals.values.size() << ") = ";
   //  for(Int_t j = 0; j < PMatrix_vals.values.size(); j++){

   //    if(TMath::Abs(PMatrix_vals.values[j]) > 1e+8){
   // 	cout << PMatrix_vals.values[j] << ", ";
   //    }
   
      
   //  }

   //  cout << " " << endl;
   //  cout << "New matrix values: D_elems (size " << DMatrix_vals.values.size() << ") = ";
   //  for(Int_t j = 0; j < DMatrix_vals.values.size(); j++){
   //    if(TMath::Abs(DMatrix_vals.values[j]) > 1e+8){
   // 	cout << DMatrix_vals.values[j] << ", ";
   //    }
   
      
   //  }

   

   //  cout << " " << endl;


    
   //  cout << "loop = " << i << endl;


    
   // cout << "bytes = " << bytes << endl;



    // calculate reactz


    // temp 

    
    const Int_t a = (HRSAngle > 0) ? 1 : -1;


    TVector3 BeamSpotHCS(0,0,0);


    if( rast!= 0){

      reactz = - ( y_tgt -a*MissPointZ)*TMath::Cos(TMath::ATan(ph_tgt))/TMath::Sin(HRSAngle + TMath::ATan(ph_tgt)) + Lrb_x*TMath::Cos(HRSAngle+TMath::ATan(ph_tgt))/TMath::Sin(HRSAngle+TMath::ATan(ph_tgt));

      BeamSpotHCS.SetXYZ(Lrb_x,Lrb_y,reactz);


    }
    else{
      
      reactz = - ( y_tgt -a*MissPointZ)*TMath::Cos(TMath::ATan(ph_tgt))/TMath::Sin(HRSAngle + TMath::ATan(ph_tgt)) + Lurb_x*TMath::Cos(HRSAngle+TMath::ATan(ph_tgt))/TMath::Sin(HRSAngle+TMath::ATan(ph_tgt));

      BeamSpotHCS.SetXYZ(Lurb_x,Lurb_y,reactz);

      
    }



    // calculating x_tg and sieve x and y


    
    //    reactz = -0.205;
    

    //    TVector3 BeamSpotHCS(0,0,-0.205);



    TVector3 BeamSpotTCS=fTCSInHCS.Inverse()*(BeamSpotHCS-fPointingOffset);

    x_tgt = BeamSpotTCS.X() - BeamSpotTCS.Z() * TMath::Tan(th_tgt);



    x_sieve = x_tgt + (TMath::Tan(th_tgt) + x_tgt*ExtTarCor_ThetaCorr) * (ZPos); 


    

    Double_t y_tgt_alt = BeamSpotTCS.Y() - BeamSpotTCS.Z() * TMath::Tan(ph_tgt);

    

    //    x_sieve = x_tgt + (TMath::Tan(th_tgt) + x_tgt*ExtTarCor_ThetaCorr) * (reactz); 

    //    y_sieve = y_tgt + TMath::Tan(ph_tgt) * (ZPos);
    y_sieve = y_tgt + TMath::Tan(ph_tgt) * (ZPos);

    y_sieve_alt = y_tgt_alt + TMath::Tan(ph_tgt) * (ZPos);



    //  y_sieve = y_tgt + TMath::Tan(ph_tgt) * (SieveHoleTCS.Z() - BeamSp);


    if(TMath::Abs(th_tgt) < 1e+8 &&  TMath::Abs(y_tgt) < 1e+8 && TMath::Abs(ph_tgt) < 1e+8 && TMath::Abs(dp_tgt) < 1e+8 && theta_track != th_tgt){


    // if(TMath::Abs(th_tgt) < 1e+8 &&  TMath::Abs(y_tgt) < 1e+8 && TMath::Abs(ph_tgt) < 1e+8 && TMath::Abs(dp_tgt) < 1e+8){

      //    if(bytes<read_lim){
     
      
      
      // cout << "before test: out_T->GetEntries() = " << out_T->GetEntries() << endl << endl;
      
      out_T->Fill();



      // cout << "after test: out_T->GetEntries() = " << out_T->GetEntries() << endl << endl;
      //    }
    }
    else{
      // cout << " " << endl;
    }

    theta_track = th_tgt;
    
  }


  out_T->Write();


  //#pragma link C++ class vector<DB_entry> +;
  //#pragma link C++ class vector<THaMatrixElement>+;


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

  //  DB_file.open("rootfiles/" + Dir_name + Form("/DB_%i.dat",runnumber),std::ofstream::out);
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

  // newfile->WriteObject(&TDB_co,"T_coeffs");
  // newfile->WriteObject(&TDB_pow,"T_powers");

  // newfile->WriteObject(&YDB_co,"Y_coeffs");
  // newfile->WriteObject(&YDB_pow,"Y_powers");

  // newfile->WriteObject(&PDB_co,"P_coeffs");
  // newfile->WriteObject(&PDB_pow,"P_powers");

  // newfile->WriteObject(&DDB_co,"D_coeffs");
  // newfile->WriteObject(&DDB_pow,"D_powers");


  // if(TFile* file_test = TFile::Open("rootfiles/" + db_name + Form("/%d_replay_1.root",runnumber),"read")){
  //   cout << "opened new file! " << endl;
  //   file_test->Close();    
  // }
  // else{
  //   newfile->Close();
  //   cout << "closed newfile" << endl;


  // }


  delete out_T;
  delete T;
  newfile->Close();
  cout << "closed newfile" << endl;




  //  delete newfile;
  
  

}



// void CalcMatrixElem_APEX(vector<THaMatrixElement>& matrix )
// {
//   // calculates value for each event of each combination of polynomial ie for T A B C this calculates for theta^A, y^B and phi^C and calculates for all powers of X related to this as well


//   // For:
//   // T 2 1 0  -1.575120e+01 -5.741883e+02  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
//   // values will be calculated for th^2*y^1*ph^0*x^0 and for th^2*y^1*ph^0*x^1



//   // loop through all of one kind of element
//   for( vector<THaMatrixElement>::iterator it=matrix.begin(); it!=matrix.end(); it++ ) {


//     // if(!it){
//     //   std::cout << "iterator not valid " << std::endl;
//     // }
    
//     for(int i=0; i<=it->order-1; i++){
//       // order here is power of X
//       // it->order = (max power of x) +1
//       // it->poly[i] is coeffecient of x^i
//       // it->pw[i] is For i = 0 power of theta,
//       //                  i = 1 power of y,
//       //                  i = 2 power of phi,
      
      
//       // combine powers of x, theta, y and phi to
//       // 
//       it->values[i][it->pw[0]][it->pw[1]][it->pw[2]] = it->poly[i] * powers[i][0] * powers[it->pw[0]][1] * powers[it->pw[1]][2]  * powers[it->pw[2]][3];
      
      
      
//     }

    
//   }



    



//   // }




// }


