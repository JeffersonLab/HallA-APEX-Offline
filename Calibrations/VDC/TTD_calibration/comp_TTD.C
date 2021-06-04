/**************************************************************
  comp_TTD.C
  John Williamson
  Feb 22, 2021

  This script compares different methods of TTD conversion: lookup table and analytic.
     
*************************************************************/

#include <fstream>

#include "Load_more_rootfiles.C"
#include "file_def.h"
#include "TTD_namespace.h"
#include "TTDTable.h"


#define NPLANE 4
#define MAX_ENTRIES 1000000
#define MAX_HIT 1000

#define DEG_TO_RAD 0.017453278


Double_t ReadDriftVel(std::string line);

Double_t* ReadAParams(std::string line);

void comp_TTD(const char *arm = "L", Int_t runnumber = -1  ){
  
  const char plane[NPLANE][8] = {"u1", "u2", "v1", "v2"};

  // set-up variables to read analytical DB into

  Double_t driftvel[NPLANE] = {0};

  Double_t an_pars[NPLANE][8] = {0}; // 4 a1 and 4 a3 parameters for each VDC plane

  Double_t Lookup_pars[NPLANE][2] = {0}; // 2 parameters for lookup table correction

  Double_t Ext_pars[NPLANE][4] = {0}; // 2 parameters for lookup table correction
  
  TString ana_DB_name = Form("DB/analytic_TTD/db_%s_TTD.vdc.%d.dat",arm,runnumber);

  std::ifstream ana_DB(Form("DB/analytic_TTD/db_%s_TTD.vdc.%d.dat",arm,runnumber));

  char* phrase = NULL;

  Int_t line_no = 0;
  
  std::string line;
  
  while (std::getline(ana_DB, line))
    {
      std::istringstream iss(line);

      Int_t plane_no = 0;
      
      for(auto pl : plane){
	
	phrase = Form("%s.vdc.%s.driftvel",arm,pl);
	if(!line.find(phrase)){
	  
	  std::getline(ana_DB, line);
	  
	  driftvel[plane_no] = ReadDriftVel(line);

	  line_no++;
	  
	}

	phrase = Form("%s.vdc.%s.ttd.param =",arm,pl);

	if(!line.find(phrase)){

	  std::getline(ana_DB, line);	  

	  Double_t* a1_pars = ReadAParams(line);
	  
	  for(Int_t j = 0; j<4; j++){
	    an_pars[plane_no][j] = a1_pars[j];
	   }
	  
	  std::getline(ana_DB, line);

	  Double_t* a2_pars = ReadAParams(line);
	  
	  for(Int_t j = 0; j<4; j++){
	    an_pars[plane_no][j+4] = a2_pars[j];
	  }

	  line_no++;	  

	}	

	plane_no++;
      }
      
      
      cout << endl;

      line_no++;
    }

  
  cout << "Testing analytic DB reading" << endl;
  
  for(Int_t i = 0; i < NPLANE; i++){

    cout << plane[i] << " driftvel = " << driftvel[i] << endl;

    for( Int_t j = 0; j<4; j++ ){      
      cout << "a1," << j << " = " << an_pars[i][j] << ", ";
    }
    cout << endl;
    for( Int_t j = 0; j<4; j++ ){
      cout << "a2," << j << " = " << an_pars[i][j+4] << ", ";
    }
    
    cout << endl << endl;
  }



  // read in drift velocity table

  // set-up variables to read analytical DB into

  Int_t NBins[NPLANE] = {0}; // number of bins lookup table for each plane
  Double_t Low[NPLANE] = {0};
  
  std::vector<Double_t> LTable[NPLANE]; // velocity lookup table
  TTDTable* PlaneTable[NPLANE]; 



  
  // read in extension parameters for Lookup table
  
  std::ifstream tableExt_DB(Form("DB/ALT_analytic_TTD/db_%s_TTD.vdc.%d.dat", arm, runnumber));
  
  cout << "Opened " << Form("DB/ALT_analytic_TTD/db_%s_TTD.vdc.%d.dat", arm, runnumber) << endl;




    
  line_no = 0;
  
  while (std::getline(tableExt_DB, line))
    {
      
      
      std::istringstream iss(line);
      
      for(Int_t i = 0; i < NPLANE; i++){	

	
	//L.vdc.u1.ttd.param =
	phrase = Form("%s.vdc.%s.ttd.param =",arm,plane[i]);

	if(!line.find(phrase)){
	  cout << "found phrase " << phrase << endl;

	  std::getline(tableExt_DB, line);
	  
	  //	  Double_t* expars = ReadAParams(line);       

	  

	  for(Int_t j = 0; j<4; j++){
	    Ext_pars[i][j] = TTD_func::ReadSingleVal<Double_t>(line);
	    //	    Ext_pars[i][j] = expars[j];
	    cout << "Ext_pars[" << i<< "][" << j << "] = " << Ext_pars[i][j] << endl;
	    std::getline(tableExt_DB, line);
	  }

	  line_no++;
	}
	
	
	
      }
      line_no++;
    }
  
  cout << "Succesfull read Ext pars " << endl << endl;

  
  // loop over planes

  for(Int_t i = 0; i<NPLANE; i++){

    // set values of Lookup angular correction parameters

    // Lookup_pars[i][0] = 0.0019;
    // Lookup_pars[i][1] = 1.35;
    

    std::ifstream tableAng_DB(Form("DB/lookup_tables/db_%s_%s_lookup_TTD_angleCorr.vdc.%d.dat", arm, plane[i], runnumber));

    cout << "Opened " << Form("DB/lookup_tables/db_%s_%s_lookup_TTD_angleCorr.vdc.%d.dat", arm, plane[i], runnumber) << endl;
    
    line_no = 0;

    while (std::getline(tableAng_DB, line))
      {


	std::istringstream iss(line);
	phrase = Form("%s.vdc.%s.ttd_table.R = ",arm,plane[i]);

	if(!line.find(phrase)){

	  std::getline(tableAng_DB, line);
	  
	  Lookup_pars[i][0] = TTD_func::ReadSingleVal<Double_t>(line);

	  line_no++;
	}

	phrase = Form("%s.vdc.%s.ttd_table.theta0 = ",arm,plane[i]);

	if(!line.find(phrase)){

	  cout << phrase << " found in line " << line_no << endl;
	  std::getline(tableAng_DB, line);
	  
	  Lookup_pars[i][1] = TTD_func::ReadSingleVal<Double_t>(line);

	  line_no++;
	}

	line_no++;
      }
    
    
    std::ifstream table_DB(Form("DB/lookup_tables/db_%s_%s_lookup_TTD_norm.vdc.%d.dat",arm,plane[i],runnumber));

    cout << "Reading " << Form("DB/lookup_tables/db_%s_%s_lookup_TTD_norm.vdc.%d.dat",arm,plane[i],runnumber) << endl;
    
    Int_t line_no = 0;
  
    std::string line;

    while (std::getline(table_DB, line))
      {
	std::istringstream iss(line);

	phrase = Form("%s.vdc.%s.ttd_table.nbins = ",arm,plane[i]);

	if(!line.find(phrase)){
	  
	  std::getline(table_DB, line);

	  NBins[i] = TTD_func::ReadNBins(line);

	  line_no++;

	}

	phrase = Form("%s.vdc.%s.ttd_table.low = ",arm,plane[i]);
	
	if(!line.find(phrase)){
	  
	  
	  std::getline(table_DB, line);

	  Low[i] = TTD_func::ReadSingleVal<Double_t>(line);

	  line_no++;

	}
	

	
	phrase = Form("%s.vdc.%s.ttd_table.table",arm,plane[i]);

	if(!line.find(phrase)){

	  LTable[i] = TTD_func::ReadLookupTable(table_DB, line_no, NBins[i]);
	  PlaneTable[i] = new TTDTable(LTable[i],Low[i],NBins[i],Lookup_pars[i],Ext_pars[i]);
	  
	  line_no++;
	}	

	line_no++;
      }  

    
        
  }



  cout << "Test reading of Table parameters" << endl;


  for(Int_t i = 0; i < NPLANE; i++){

    cout << "For " << plane[i] << ", NBins = " << NBins[i] << endl;

  }


  for(Int_t i = 0; i < NPLANE; i++){

    cout << "For " << plane[i] << ", Table = " << endl;
    
    Int_t j = 1;
    for( auto val : LTable[i]){
      cout << val << " ";
      if(j%10==0){
	cout << endl;
      }
      
      j++;
    }

    cout << endl << endl;
  }

  
  cout << "Print angular parameters from TTDTable objects: " << endl << endl;
  for(Int_t i = 0; i < NPLANE; i++){

    cout << "For " << plane[i] << ": " << endl;
    PlaneTable[i]->PrintParams();
    cout << endl << endl;
  }
    
  cout << endl << endl;




  // Load Root file
  
  
  TChain* T = new TChain("T");
  
  if(!strcmp(arm,"L")){
  T = Load_more_rootfiles(runnumber);
  }
  else if (!strcmp(arm,"R")){  
    //  TChain *T = new TChain("T");
  T = Load_more_rootfiles(runnumber);
  }
  else{
    cout << "arm must be L or R, " << arm << " not acceptable" << endl;
    return;
  }




  // declare cut variables
  
  // PID cuts (different for Left and Right arms)
  Double_t Cer_cut = 0.0;
  Double_t Ps_cut = 0.0;
  Double_t Ps_Sh_cut_l = 0.0;
  Double_t Ps_Sh_cut_h = 0.0;
  

  // define cuts based on arm 
  
  Int_t trg = 0;
  // here is where trigger is defined
  if(!strcmp(arm,"L")){
    trg = 1;
    Cer_cut = 1500.0;
    Ps_cut = 0.3;
    Ps_Sh_cut_l = 0.625;
    Ps_Sh_cut_h = 1.1;
    cout << "left trigger" << endl;
  }
  else if (!strcmp(arm,"R")){
    cout << "right trigger" << endl;
    trg = 5;
    Cer_cut = 650.0;
    Ps_cut = 0.2;
    Ps_Sh_cut_l = 0.51;
    Ps_Sh_cut_h = 1.11;
    // trg = 6;
  }
  else{
    cout << "arm must be L or R, " << arm << " not acceptable" << endl;
    return;

  }

  
  Double_t ang[NPLANE] = {-45.0, -45.0, 45.0, 45.0};

  Double_t nhit[NPLANE], ntr;
  Double_t hittime[NPLANE][MAX_HIT], hittrknum[NPLANE][MAX_HIT],
    hittrdist[NPLANE][MAX_HIT];
  Double_t d_th[MAX_HIT], d_ph[MAX_HIT];
  Double_t cer_sum, ps_e, sh_e;
  Double_t tr_p[100];
  THaEvent* evt = 0;

  Int_t    nent[NPLANE];
  
  Int_t i, j, hit;

  Double_t evttype;


      // cluster slope and intercept from analyzer
  Double_t slope[NPLANE][MAX_HIT];
  Double_t intercept[NPLANE][MAX_HIT]; 

  Double_t nclust[NPLANE];
  Double_t clsiz[NPLANE][MAX_HIT];
  Double_t wireno[NPLANE][MAX_HIT];
  Double_t clpivot[NPLANE][MAX_HIT];
  Double_t clbeg[NPLANE][MAX_HIT];
  Double_t clend[NPLANE][MAX_HIT];
  
  

  // Set up branches

  T->SetBranchStatus("Event_Branch*", kTRUE);
  T->SetBranchAddress("Event_Branch", &evt);

  T->SetBranchStatus("DR.evtypebits", kTRUE);
  T->SetBranchAddress("DR.evtypebits", &evttype);

  T->SetBranchStatus(Form("%s.tr.n", arm), kTRUE);
  T->SetBranchAddress(Form("%s.tr.n", arm), &ntr);

  T->SetBranchStatus(Form("%s.tr.d_th", arm), kTRUE);
  T->SetBranchAddress(Form("%s.tr.d_th", arm), &d_th);
  T->SetBranchStatus(Form("%s.tr.d_ph", arm), kTRUE);
  T->SetBranchAddress(Form("%s.tr.d_ph", arm), &d_ph);

  T->SetBranchStatus(Form("%s.tr.p",arm),kTRUE);
  T->SetBranchAddress(Form("%s.tr.p",arm),tr_p);
  
  T->SetBranchStatus(Form("%s.cer.asum_c",arm),kTRUE);
  T->SetBranchAddress(Form("%s.cer.asum_c",arm),&cer_sum);
  
  if(!strcmp(arm,"L")){
    T->SetBranchStatus("L.prl1.e",kTRUE);
    T->SetBranchAddress("L.prl1.e",&ps_e);
    T->SetBranchStatus("L.prl2.e",kTRUE);
    T->SetBranchAddress("L.prl2.e",&sh_e);
  }
  else if(!strcmp(arm,"R")){
    T->SetBranchStatus("R.ps.e",kTRUE);
    T->SetBranchAddress("R.ps.e",&ps_e);
    T->SetBranchStatus("R.ps.e",kTRUE);
    T->SetBranchAddress("R.sh.e",&sh_e);
  }
  

  for( i = 0; i < NPLANE; i++ ){
    T->SetBranchStatus(Form("%s.vdc.%s.nhit", arm, plane[i]), kTRUE);
    T->SetBranchAddress(Form("%s.vdc.%s.nhit", arm, plane[i]), &nhit[i]);

    T->SetBranchStatus(Form("%s.vdc.%s.time", arm, plane[i]), kTRUE);
    T->SetBranchAddress(Form("%s.vdc.%s.time", arm, plane[i]), hittime[i]);

    T->SetBranchStatus(Form("%s.vdc.%s.trknum", arm, plane[i]), kTRUE);
    T->SetBranchAddress(Form("%s.vdc.%s.trknum", arm, plane[i]), hittrknum[i]);

//     T->SetBranchStatus(Form("%s.vdc.%s.ltrdist", arm, plane[i]), kTRUE);
//     T->SetBranchAddress(Form("%s.vdc.%s.ltrdist", arm, plane[i]), hittrdist[i]);
    T->SetBranchStatus(Form("%s.vdc.%s.trdist", arm, plane[i]), kTRUE);
    T->SetBranchAddress(Form("%s.vdc.%s.trdist", arm, plane[i]), hittrdist[i]);

    T->SetBranchStatus(Form("%s.vdc.%s.slope", arm, plane[i]), kTRUE);
    T->SetBranchAddress(Form("%s.vdc.%s.slope", arm, plane[i]), slope[i]);
    
    nent[i] = 0;
  }




  Double_t this_slope;



  // Double_t wtime[NPLANE][MAX_ENTRIES], tanTh[NPLANE][MAX_ENTRIES], trdist[NPLANE][MAX_ENTRIES];

  std::vector<Double_t> wtime[NPLANE], tanTh[NPLANE], trdist[NPLANE], vslope[NPLANE];
  
  // Load up drift times and tangents





  // set-up histograms

  // real distributions

  // plot drift time spectra for all planes
  TH1F *htime[NPLANE];

  // plot 'real' drift distance spectra for all planes
  TH1F *hdist[NPLANE];
  TH1F *hdistTable[NPLANE];
  TH1F *hdistTableCorr[NPLANE];
  TH1F* hdistAna[NPLANE];

  // plot 'real' drift distance versus time spectra for all planes
  TH2D *htime_dist_real[NPLANE];
  TH1D *hprof_time_dist_real[NPLANE];


  // plot table drift distance versus time spectra for all planes
  TH2D *htime_dist_table[NPLANE];

  // plot table drift distance versus 'real' distances for all planes
  TH2D *hreal_dist_table[NPLANE];

  
  // plot table drift distance with angular correction versus time spectra for all planes
  TH2D *htime_dist_tableCorr[NPLANE];

  // plot table drift distance with angular correction versus 'real' distances for all planes
  TH2D *hreal_dist_tableCorr[NPLANE];
  

  // plot analytical drift distance versus time spectra for all planes
  TH2D *htime_dist_ana[NPLANE];

  // plot analytical drift distance versus 'real' distances for all planes
  TH2D *hreal_dist_ana[NPLANE];


  

  // slice distributions
  TH1D *h_means[NPLANE];
  TH1D *h_sigma[NPLANE];
  TH1D *h_const[NPLANE];

  // plot sigma limits above and below means of distribution
  TH1D *h_fit_up[NPLANE];
  TH1D *h_fit_low[NPLANE];


  // chi^2 for both TTD methods
  
  TH1D *h_chi2_table[NPLANE];
  TH1D *h_chi2_tableCorr[NPLANE]; // with angular correction
  TH1D *h_chi2_ana[NPLANE];

  TH1D *h_chi_table[NPLANE];
  TH1D *h_chi_tableCorr[NPLANE]; // with angular correction
  TH1D *h_chi_ana[NPLANE]; 

  
  // plot slopes for all planes

  TH1D *hslope[NPLANE];
  TH1D *hslopeAlt[NPLANE];
  

  for( i = 0; i < NPLANE; i++ ){
    
    htime_dist_real[i] = new TH2D(Form("htime_dist_real_%s", plane[i]), Form("%s TTD", plane[i]), 760, -30, 350, 200, 0.0, 0.020);
    htime_dist_real[i]->GetXaxis()->SetTitle("Drift Time(ns)");
    htime_dist_real[i]->GetXaxis()->CenterTitle();
    htime_dist_real[i]->GetYaxis()->SetTitle("Dist from Time (m)");
    htime_dist_real[i]->GetYaxis()->CenterTitle();

    
    htime_dist_table[i] = new TH2D(Form("htime_dist_table_%s", plane[i]), Form("%s TTD Lookup Table", plane[i]), 760, -30, 350, 200, 0.0, 0.020);
    htime_dist_table[i]->GetXaxis()->SetTitle("Drift Time(ns)");
    htime_dist_table[i]->GetXaxis()->CenterTitle();
    htime_dist_table[i]->GetYaxis()->SetTitle("Dist from Time (m)");
    htime_dist_table[i]->GetYaxis()->CenterTitle();

    htime_dist_tableCorr[i] = (TH2D*) htime_dist_table[i]->Clone(Form("htime_dist_tableCorr_%s", plane[i]));
    htime_dist_tableCorr[i]->SetTitle(Form("%s TTD Lookup table (with angular correction)", plane[i]));

    htime_dist_ana[i] = new TH2D(Form("htime_dist_ana_%s", plane[i]), Form("%s TTD Analytic function", plane[i]), 760, -30, 350, 200, 0.0, 0.020);
    htime_dist_ana[i]->GetXaxis()->SetTitle("Drift Time(ns)");
    htime_dist_ana[i]->GetXaxis()->CenterTitle();
    htime_dist_ana[i]->GetYaxis()->SetTitle("Dist from Time (m)");
    htime_dist_ana[i]->GetYaxis()->CenterTitle();
    

    htime[i] = new TH1F(Form("htime_%s", plane[i]), Form("%s TTD", plane[i]), 760, -30, 350);
    htime[i]->GetXaxis()->SetTitle("Drift Time(ns)");
    htime[i]->GetXaxis()->CenterTitle();

    hdist[i] = new TH1F(Form("hdist_%s", plane[i]), Form("%s TTD", plane[i]), 200, -0.001, 0.020);
    hdist[i]->GetXaxis()->SetTitle("Track Dist (m)");
    hdist[i]->GetXaxis()->CenterTitle();

    hdistTable[i] = (TH1F*) hdist[i]->Clone(Form("hdistTable_%s", plane[i]));
    hdistTable[i]->SetTitle(Form("%s TTD", plane[i]));

    hdistTableCorr[i] = (TH1F*) hdist[i]->Clone(Form("hdistTableCorr_%s", plane[i]));
    hdistTableCorr[i]->SetTitle(Form("%s TTD", plane[i]));

    hdistAna[i] = (TH1F*) hdist[i]->Clone(Form("hdistAna_%s", plane[i]));
    hdistAna[i]->SetTitle(Form("%s TTD", plane[i]));

    hreal_dist_table[i] = new TH2D(Form("hreal_dist_table_%s", plane[i]), Form("%s TTD Lookup table", plane[i]), 200, 0.0, 0.020, 200, 0.0, 0.020);
    hreal_dist_table[i]->GetXaxis()->SetTitle("'Real' Track Dist (m)");
    hreal_dist_table[i]->GetXaxis()->CenterTitle();
    hreal_dist_table[i]->GetYaxis()->SetTitle("Table Track Dist (m)");
    hreal_dist_table[i]->GetYaxis()->CenterTitle();

    hreal_dist_tableCorr[i] = (TH2D*) hreal_dist_table[i]->Clone(Form("hreal_dist_tableCorr_%s", plane[i]));
    hreal_dist_tableCorr[i]->SetTitle(Form("%s TTD Lookup table (with angular correction)", plane[i]));

    hreal_dist_ana[i] = new TH2D(Form("hreal_dist_ana_%s", plane[i]), Form("%s TTD analytical", plane[i]), 200, 0.0, 0.020, 200, 0.0, 0.020);
    hreal_dist_ana[i]->GetXaxis()->SetTitle("'Real' Track Dist (m)");
    hreal_dist_ana[i]->GetXaxis()->CenterTitle();
    hreal_dist_ana[i]->GetYaxis()->SetTitle("Analytic Track Dist (m)");
    hreal_dist_ana[i]->GetYaxis()->CenterTitle();



    h_chi2_table[i] = new TH1D(Form("h_chi2_table_%s", plane[i]), Form("%s #chi^{2}", plane[i]), 200, 0.0, 1e-6);
    h_chi2_table[i]->GetXaxis()->SetTitle("#chi^{2}");
    h_chi2_table[i]->GetXaxis()->CenterTitle();

    h_chi2_tableCorr[i] = (TH1D*) h_chi2_table[i]->Clone(Form("h_chi2_tableCorr_%s", plane[i]));
    h_chi2_tableCorr[i]->SetTitle(Form("%s #chi^{2}", plane[i]));

    h_chi2_ana[i] = (TH1D*) h_chi2_table[i]->Clone(Form("h_chi2_tableCorr_%s", plane[i]));
    h_chi2_ana[i]->SetTitle(Form("%s #chi^{2}", plane[i]));



    h_chi_table[i] = new TH1D(Form("h_chi_table_%s", plane[i]), Form("%s Real distance - TTD distance", plane[i]),  100, -0.0015, 0.0015);
    h_chi_table[i]->GetXaxis()->SetTitle("Real distance - TTD distance");
    h_chi_table[i]->GetXaxis()->CenterTitle();

    h_chi_tableCorr[i] = (TH1D*) h_chi_table[i]->Clone(Form("h_chi_tableCorr_%s", plane[i]));
    h_chi_tableCorr[i]->SetTitle(Form("%s Real distance - TTD distance",plane[i]));
    
    h_chi_ana[i] = (TH1D*) h_chi_table[i]->Clone(Form("h_chi_ana_%s", plane[i]));
    h_chi_ana[i]->SetTitle(Form("%s Real distance - TTD distance",plane[i]));

    hslope[i] = new TH1D(Form("hslope_%s",plane[i]), Form("slope distribution %s",plane[i]), 100,1.0,2.0);
    hslope[i]->GetXaxis()->SetTitle("Slope");
    hslope[i]->GetXaxis()->CenterTitle();


    hslopeAlt[i] = (TH1D*) hslope[i]->Clone(Form("hslopeAlt_%s",plane[i]));
    hslopeAlt[i]->SetTitle(Form("slope distribution %s",plane[i]));
					    
  }

  
  for( i = 0; i < T->GetEntries(); i++ ){
    // for( i = 0; i < 1e4; i++ ){
    T->GetEntry(i);
    if( (i%5000)==0 ) { cout << "Entry " <<  i << endl; }
    
    
    if( ntr == 1 && TTD_func::passtrg(Int_t(evttype), trg) && cer_sum > Cer_cut && ps_e/(1e3*tr_p[0]) > Ps_cut && (ps_e+sh_e)/(1e3*tr_p[0]) > Ps_Sh_cut_l &&  (ps_e+sh_e)/(1e3*tr_p[0]) < Ps_Sh_cut_h ){      
      for( j = 0; j < NPLANE; j++ ){
	this_slope = d_th[0]*cos(ang[j]*DEG_TO_RAD)
	  + d_ph[0]*sin(ang[j]*DEG_TO_RAD);;
	for( hit = 0; hit < nhit[j] && nent[j] < MAX_ENTRIES; hit++)
	  {
	    
	    //	  if( 12e-9 < hittime[j][hit]
	    if( //0 < hittime[j][hit] &&
	       // hittime[j][hit]  < 350e-9
	       hittime[j][hit]  < 350e-9
	       && hittrdist[j][hit]< 15.12e-3
	       && hittrknum[j][hit] == 1
	       // && wireno[j][hit] >= clbeg[j][0]
	       // && wireno[j][hit] <= clend[j][0]
	       
		)
	      {
				  
		  
		// fill distribution arrays	    
		wtime[j].push_back(hittime[j][hit]);
		tanTh[j].push_back(this_slope);
		trdist[j].push_back(hittrdist[j][hit]);
		vslope[j].push_back(slope[j][0]);
		  
		// wtime[j][nent[j]]  = hittime[j][hit];
		// tanTh[j][nent[j]] = this_slope;
		// trdist[j][nent[j]] = hittrdist[j][hit];
		  
		Double_t time = wtime[j][nent[j]];
		Double_t time_ns = 1e9*(time); // time converted to nanoseconds
		Double_t dist = trdist[j][nent[j]];
		Double_t slope_alt = vslope[j][nent[j]]; // slope from analyzer

		htime[j]->Fill(time_ns);
		hdist[j]->Fill(dist);
		htime_dist_real[j]->Fill(time_ns,dist);

		Double_t ana_dist = TTD_func::TTDAna(time, this_slope, driftvel[j],an_pars[j]);
		htime_dist_ana[j]->Fill(time_ns,ana_dist);
		hreal_dist_ana[j]->Fill(dist,ana_dist);
	    
		//	    Double_t table_dist = TTD_table(time, LTable[j], NBins[j],this_slope,an_pars[j]);
		Double_t table_dist = PlaneTable[j]->Convert(time);

		//	    Double_t table_dist = PlaneTable[j]->ConvertAngleCorr(time,this_slope);
		htime_dist_table[j]->Fill(time_ns,table_dist);
		hreal_dist_table[j]->Fill(dist,table_dist);

		//	    Double_t table_distCorr = PlaneTable[j]->ConvertAngleCorr(time,1/this_slope);
		Double_t table_distCorr = PlaneTable[j]->ConvertAngleCorr(time,1/this_slope);
		htime_dist_tableCorr[j]->Fill(time_ns,table_distCorr);
		hreal_dist_tableCorr[j]->Fill(dist,table_distCorr);
	    

		hdistTable[j]->Fill(table_dist);
		hdistTableCorr[j]->Fill(table_distCorr);
		hdistAna[j]->Fill(ana_dist);
		
		h_chi2_table[j]->Fill(TMath::Power((dist-table_dist),2));
		h_chi2_tableCorr[j]->Fill(TMath::Power((dist-table_distCorr),2));
		h_chi2_ana[j]->Fill(TMath::Power((dist-ana_dist),2));


		h_chi_table[j]->Fill(dist-table_dist);
		h_chi_tableCorr[j]->Fill(dist-table_distCorr);
		h_chi_ana[j]->Fill(dist-ana_dist);
		

		hslope[j]->Fill(1/this_slope);
		hslopeAlt[j]->Fill(1/slope_alt);
		
		nent[j]++;	    
		
	      }
	  }
      }
    }
  }
  


  // draw real distributions

  TCanvas *c_real[NPLANE];

  TCanvas *c_dist = new TCanvas("c_dist", "Distance distributions",640,480);
  c_dist->Divide(2,2);

  TLegend *leg_dist =  new TLegend(.05,.65,.45,.9,"Key");
  leg_dist->SetFillColor(0);
  leg_dist->SetTextSize(0.025);

  THStack* hDistStack[NPLANE];

  TCanvas *c_slices[NPLANE];

  // TTD results for table and analytic method
  TCanvas *c_TTD[NPLANE];
  TLegend *leg_TTD[NPLANE];

  // function that display perfect linear correlation between 'real' distance and distance caluclated from time (y = x)
  TF1 *Lin_cor = new TF1("Lin_cor","pol1(0)", 0.0, 0.015);
  Lin_cor->SetParameter(0,0.0);
  Lin_cor->SetParameter(1,1.0);
  Lin_cor->SetLineStyle(4); // dashed line
  

  Double_t XSigma_up = 2.0; // how many sigma to cut away from mean of 'real' distance distrib

  Double_t XSigma_down = 2.0; // how many sigma to cut away from mean of 'real' distance distrib
    
  
  TCanvas *c_Chi2[NPLANE];
  TLegend *leg_Chi2[NPLANE];


  TCanvas *c_slope = new TCanvas("c_slope", "Slopes",1400,1400);
  c_slope->Divide(2,2);
  TLegend *leg_slope =  new TLegend(.45,.65,.9,.9,"Key");
  leg_slope->SetFillColor(0);
  leg_slope->SetTextSize(0.025);
  
  for( i = 0; i < NPLANE; i++ ){
    
    c_real[i] = new TCanvas(Form("c_real_%s",plane[i]), Form("Real distributions %s",plane[i]),640,480);

    c_real[i]->Divide(2,2);

    c_real[i]->cd(1);
    htime[i]->Draw();
    c_real[i]->cd(2);
    hdist[i]->Draw();
    
    c_real[i]->cd(3);
    htime_dist_real[i]->Draw("colz");


    c_dist->cd(i+1);

    hDistStack[i]  = new THStack(Form("hDistStack_%i",i),Form("Dist %s; Dist",plane[i]));

    hdist[i]->SetLineColor(kBlack);
    hDistStack[i]->Add(hdist[i]);
    
    hdistAna[i]->SetLineColor(kRed);
    hDistStack[i]->Add(hdistAna[i]);

    hdistTable[i]->SetLineColor(kBlue);
    hDistStack[i]->Add(hdistTable[i]);

    hdistTableCorr[i]->SetLineColor(kGreen+2);
    hDistStack[i]->Add(hdistTableCorr[i]);

        
    hDistStack[i]->Draw("nostack");
    
    if(i == 0){
      // create dist legend
      leg_dist->AddEntry(hdist[i],"'Real' track dist ","l");
      leg_dist->AddEntry(hdistAna[i],"Analytic method","l");
      leg_dist->AddEntry(hdistTable[i],"Dist from Lookup table","l");
      leg_dist->AddEntry(hdistTableCorr[i],"Dist from Lookup table w corr","l");            
    }

    leg_dist->Draw("same");


    
    

    // fit slices in Y ('real' distance), over all bins in x (time), with cut off of 10 bins being filled in y)
    htime_dist_real[i]->FitSlicesY(0,0,-1,10,"QN");

    
    c_slices[i] = new TCanvas(Form("c_slices_%s",plane[i]), Form("Slices distributions %s",plane[i]),640,480);

    c_slices[i]->Divide(2,2);

    c_slices[i]->cd(1);
    htime_dist_real[i]->Draw("colz");

    c_slices[i]->cd(2);
    h_const[i] = (TH1D*)gDirectory->Get(Form("htime_dist_real_%s_0",plane[i]));
    h_const[i]->Draw();

    c_slices[i]->cd(3);
    h_means[i] = (TH1D*)gDirectory->Get(Form("htime_dist_real_%s_1",plane[i]));
    h_means[i]->Draw();

    c_slices[i]->cd(4);
    h_sigma[i] = (TH1D*)gDirectory->Get(Form("htime_dist_real_%s_2",plane[i]));
    h_sigma[i]->Draw();

    h_fit_up[i] = (TH1D*) h_means[i]->Clone();
    h_fit_up[i]->Add(h_sigma[i],XSigma_up);
    h_fit_up[i]->SetLineColor(kRed);

    h_fit_low[i] = (TH1D*) h_means[i]->Clone();
    h_fit_low[i]->Add(h_sigma[i],-XSigma_down);
    h_fit_low[i]->SetLineColor(kRed);


    c_slices[i]->cd(1);
    h_means[i]->SetLineColor(kBlack);
    h_means[i]->Draw("same");
    h_fit_up[i]->Draw("hist  L same");
    h_fit_low[i]->Draw("hist L same");


    c_TTD[i] = new TCanvas(Form("c_TTD_%s",plane[i]), Form("TTD distributions %s",plane[i]),1400,1400);
    c_TTD[i]->Divide(3,2);
    c_TTD[i]->cd(1);
    gPad->SetLeftMargin(0.15);
    htime_dist_ana[i]->GetYaxis()->SetTitleOffset(1.4);
    htime_dist_ana[i]->Draw("colz");

    c_TTD[i]->cd(2);
    gPad->SetLeftMargin(0.15);
    htime_dist_table[i]->GetYaxis()->SetTitleOffset(1.4);
    htime_dist_table[i]->Draw("colz");

    c_TTD[i]->cd(3);
    gPad->SetLeftMargin(0.15);
    htime_dist_tableCorr[i]->GetYaxis()->SetTitleOffset(1.4);
    htime_dist_tableCorr[i]->Draw("colz");

    
    c_TTD[i]->cd(4);
    gPad->SetLeftMargin(0.15);
    hreal_dist_ana[i]->GetYaxis()->SetTitleOffset(1.4);
    hreal_dist_ana[i]->Draw("colz");
    Lin_cor->Draw("same");
    leg_TTD[i] = new TLegend(.2,.65,.57,.9,"Key");
    leg_TTD[i]->SetFillColor(0);
    leg_TTD[i]->SetTextSize(0.025);
    leg_TTD[i]->AddEntry(Lin_cor,"Exact linear correlation","l");
    leg_TTD[i]->Draw("same");
    
    
    c_TTD[i]->cd(5);
    gPad->SetLeftMargin(0.15);
    hreal_dist_table[i]->GetYaxis()->SetTitleOffset(1.4);
    hreal_dist_table[i]->Draw("colz");
    Lin_cor->Draw("same");
    leg_TTD[i]->Draw("same");

    c_TTD[i]->cd(6);
    gPad->SetLeftMargin(0.15);
    hreal_dist_tableCorr[i]->GetYaxis()->SetTitleOffset(1.4);
    hreal_dist_tableCorr[i]->Draw("colz");
    Lin_cor->Draw("same");
    leg_TTD[i]->Draw("same");

    c_Chi2[i] = new TCanvas(Form("c_Chi2_%s",plane[i]), Form("chi2 %s",plane[i]),1400,1400);

    c_Chi2[i]->Divide(2,1);
    
    c_Chi2[i]->cd(1);

    gPad->SetLogy();
    
    h_chi2_ana[i]->SetLineColor(kRed);
    h_chi2_ana[i]->Draw();

    h_chi2_table[i]->SetLineColor(kBlue);
    h_chi2_table[i]->Draw("same");

    h_chi2_tableCorr[i]->SetLineColor(kGreen+2);
    h_chi2_tableCorr[i]->Draw("same");

    leg_Chi2[i] = new TLegend(.70,.65,.9,.9,"Key");
    leg_Chi2[i]->SetFillColor(0);
    leg_Chi2[i]->SetTextSize(0.02);
    leg_Chi2[i]->AddEntry(h_chi2_ana[i],"Analytic","l");
    leg_Chi2[i]->AddEntry(h_chi2_table[i],"Lookup table","l");
    leg_Chi2[i]->AddEntry(h_chi2_tableCorr[i],"#splitline{Lookup table}{(corrected)}","l");
    leg_Chi2[i]->Draw("same");


    c_Chi2[i]->cd(2);

    gPad->SetLogy();

    h_chi_ana[i]->SetLineColor(kRed);
    h_chi_ana[i]->Draw();

    h_chi_table[i]->SetLineColor(kBlue);
    h_chi_table[i]->Draw("same");

    h_chi_tableCorr[i]->SetLineColor(kGreen+2);
    h_chi_tableCorr[i]->Draw("same");

    leg_Chi2[i]->Draw("same");

    c_slope->cd(i+1);
    
    hslope[i]->SetLineColor(kBlue);
    hslope[i]->Draw();
    
    hslopeAlt[i]->SetLineColor(kRed);
    hslopeAlt[i]->Draw("same");

    if(i == 0){
      leg_slope->AddEntry(hslope[i],"Slope from calc","l");
      leg_slope->AddEntry(hslopeAlt[i],"Slope from analyzer","l");
    }
    
    leg_slope->Draw("same");
    
  }


  // save plots

  for( i = 0; i < NPLANE; i++ ){

    if(i == 0){
      c_TTD[i]->Print(Form("plots/comparison/comp_TTD_%s_%i.pdf(",arm,runnumber));
      c_Chi2[i]->Print(Form("plots/comparison/comp_Chi2_%s_%i.pdf(",arm,runnumber));
    }    
    else if (i == (NPLANE-1)){
      c_TTD[i]->Print(Form("plots/comparison/comp_TTD_%s_%i.pdf)",arm,runnumber));
      c_Chi2[i]->Print(Form("plots/comparison/comp_Chi2_%s_%i.pdf)",arm,runnumber));
    }
    else{
      c_TTD[i]->Print(Form("plots/comparison/comp_TTD_%s_%i.pdf",arm,runnumber));
      c_Chi2[i]->Print(Form("plots/comparison/comp_Chi2_%s_%i.pdf",arm,runnumber));
    }
     
      
    
    c_TTD[i]->Print(Form("plots/comparison/TTD_%s_%s_%i.png",arm,plane[i],runnumber));

    c_Chi2[i]->Print(Form("plots/comparison/Chi2_%s_%s_%i.png",arm,plane[i],runnumber));
    


    //plots/comparison/
  }
  

}





Double_t* ReadAParams(std::string line){

  std::istringstream iss(line);

  Double_t a_pars[4] = {0};
  
  for(Int_t j = 0; j<4; j++){
    iss >> a_pars[j];    
  }

  return a_pars;
  
}



Double_t ReadDriftVel(std::string line){


  Double_t DriftVel = 0.0;
  
  DriftVel = TTD_func::ReadSingleVal<Double_t>(line);

  return DriftVel;

}







// perform analytic TTD conversion

// Double_t TTDform( Double_t dtime, Double_t tanTheta, Double_t fDriftVel, const Double_t *par){
  
//   Double_t a1 = 0.0, a2 = 0.0;

//   const Double_t* fA1tdcCor = &(par[0]);
//   const Double_t* fA2tdcCor = &(par[4]);

//   if( fabs(tanTheta) > 1e-7 ){
//     tanTheta = 1.0 / tanTheta;
//   } else {
//     //    cerr << "TTDform: Invalid tanTheta = " << tanTheta << endl;
//     return 0.0;
//   }

//   for (Int_t i = 3; i >= 1; i--) {
//     a1 = tanTheta * (a1 + fA1tdcCor[i]);
//     a2 = tanTheta * (a2 + fA2tdcCor[i]);
//   }
//   a1 += fA1tdcCor[0];
//   a2 += fA2tdcCor[0];

//   Double_t dist = fDriftVel * dtime;

//   if (dist < 0) {
//     //    cerr << "ttdForm: invalid dist = " << dist << endl;
//     return 0;
//   // } else if (a2<0 || a1<0 || dist < 0) {

//   //   return 1e32;
//   } else if (dist < a1 ) { 
//     dist *= ( 1.0 + a2 / a1);
//   }  else {
//     dist +=  a2;
//   }

//   return dist;
// }


// performs lookup table ttd conversion
// Double_t TTD_table( Double_t dtime, std::vector<Double_t> &LTable, Int_t NBins, Double_t tanTheta, const Double_t *par){

//   // get relevant entry in lookup table from value of time
//   Double_t bin_res = 0.5e-9;
//   Int_t bin_no = dtime/(0.5e-9);
//   Double_t dist = LTable[bin_no];
  
//   Double_t dist_nocorr = dist;
//   //  cout << "bin_no = " << bin_no << ", time = " << dtime*1e9 << " ns, distance = " << LTable[bin_no]*1e-3 << endl;

  
  
//   Double_t a1 = 0.0, a2 = 0.0;

//   const Double_t* fA1tdcCor = &(par[0]);
//   const Double_t* fA2tdcCor = &(par[4]);


  
//   if( fabs(tanTheta) > 1e-7 ){
//     tanTheta = 1.0 / tanTheta;
//   } else {
//     //    cerr << "TTDform: Invalid tanTheta = " << tanTheta << endl;
//     return 0.0;
//   }

//   for (Int_t i = 3; i >= 1; i--) {
//     a1 = tanTheta * (a1 + fA1tdcCor[i]);
//     a2 = tanTheta * (a2 + fA2tdcCor[i]);
//   }
//   a1 += fA1tdcCor[0];
//   a2 += fA2tdcCor[0];


//     if (dist < 0) {
//     //    cerr << "ttdForm: invalid dist = " << dist << endl;
//     return 0;
//   // } else if (a2<0 || a1<0 || dist < 0) {

//   //   return 1e32;
//   } else if (dist < a1 ) { 
//     dist *= ( 1.0 + a2 / a1);
//   }  else {
//     dist +=  a2;
//   }

//     return dist_nocorr;
//     //    return dist;
  

  
// }


