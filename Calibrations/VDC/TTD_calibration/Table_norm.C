/**************************************************************
  Table_norm.C
  John Williamson    
  2nd March, 2021

  This script aims to calculte difference in vertical distance between two vertical wires in a cluster. Firstly clusters with 5 wires are selected, differences in track vertical distance between the 1st and 2nd (and 4th adn 5th) are calculated (determined from cluster/ track slope and known wire seperation). This is used to normal TTD lookup table, by dividing lookup table entries for relevant wires by this calculated distance. 
  
*************************************************************/


#include "Load_more_rootfiles.C"
#include "file_def.h"
#include "TTD_namespace.h"
#include "TTDTable.h"

#define MAX_ENTRIES 1000000
#define MAX_HIT 1000
#define NPLANE 4

#define DEG_TO_RAD 0.017453278

const Double_t wire_sep = 0.0042426; // seperation of wires in VDC plane

Double_t an_pars[NPLANE][8] = {0}; // 4 a1 and 4 a3 parameters for each VDC plane


std::vector<Double_t> wtime[NPLANE];
std::vector<Double_t> tanTh[NPLANE];
std::vector<Double_t> tanTh_alt[NPLANE];
std::vector<Double_t> trdist[NPLANE];

// cluster info (selecting for clusters with 5 wires)
std::vector<Int_t> WireNumbers[NPLANE];
std::vector<Int_t> CentralWireNo[NPLANE]; // Central wire no for cluster
vector <int> BegWireNo[NPLANE]; // Central wire no for cluster
std::vector<Int_t> EndWireNo[NPLANE]; // Central wire no for cluster


std::vector<Int_t> GoodWireNumbers[NPLANE]; // wire numbers for cluster passing conditions



void Table_norm(const char *arm, Int_t runnumber = -1)
{

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


  // PID cuts (different for Left and Right arms)
  Double_t Cer_cut = 0.0;
  Double_t Ps_cut = 0.0;
  Double_t Ps_Sh_cut_l = 0.0;
  Double_t Ps_Sh_cut_h = 0.0;
    
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


    
  const char plane[NPLANE][8] = {"u1", "u2", "v1", "v2"};
  Double_t ang[NPLANE] = {-45.0, -45.0, 45.0, 45.0};

  Double_t nhit[NPLANE], ntr;
  Double_t hittime[NPLANE][MAX_HIT], hittrknum[NPLANE][MAX_HIT],
    hittrdist[NPLANE][MAX_HIT];
  Double_t d_th[MAX_HIT], d_ph[MAX_HIT];
  Double_t cer_sum, ps_e, sh_e;
  Double_t tr_p[100];
  THaEvent* evt = 0;

  Int_t    nent[NPLANE];
  Int_t    ClustCount[NPLANE]; // count number of clusters sorted through
  
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


    T->SetBranchStatus(Form("%s.vdc.%s.trdist", arm, plane[i]), kTRUE);
    T->SetBranchAddress(Form("%s.vdc.%s.trdist", arm, plane[i]), hittrdist[i]);

    T->SetBranchStatus(Form("%s.vdc.%s.slope", arm, plane[i]), kTRUE);
    T->SetBranchAddress(Form("%s.vdc.%s.slope", arm, plane[i]), slope[i]);

    T->SetBranchStatus(Form("%s.vdc.%s.clpos", arm, plane[i]), kTRUE);
    T->SetBranchAddress(Form("%s.vdc.%s.clpos", arm, plane[i]), intercept[i]);

    T->SetBranchStatus(Form("%s.vdc.%s.nclust", arm, plane[i]), kTRUE);
    T->SetBranchAddress(Form("%s.vdc.%s.nclust", arm, plane[i]), &nclust[i]);

    T->SetBranchStatus(Form("%s.vdc.%s.clsiz", arm, plane[i]), kTRUE);
    T->SetBranchAddress(Form("%s.vdc.%s.clsiz", arm, plane[i]), clsiz[i]);
    
    T->SetBranchStatus(Form("%s.vdc.%s.wire", arm, plane[i]), kTRUE);
    T->SetBranchAddress(Form("%s.vdc.%s.wire", arm, plane[i]), wireno[i]);

    T->SetBranchStatus(Form("%s.vdc.%s.clpivot", arm, plane[i]), kTRUE);
    T->SetBranchAddress(Form("%s.vdc.%s.clpivot", arm, plane[i]), clpivot[i]);

    T->SetBranchStatus(Form("%s.vdc.%s.clbeg", arm, plane[i]), kTRUE);
    T->SetBranchAddress(Form("%s.vdc.%s.clbeg", arm, plane[i]), clbeg[i]);

    T->SetBranchStatus(Form("%s.vdc.%s.clend", arm, plane[i]), kTRUE);
    T->SetBranchAddress(Form("%s.vdc.%s.clend", arm, plane[i]), clend[i]);

    nent[i] = 0;
    ClustCount[i] = 0;
  }


  // Important to note that in Analyzer the 'slope' is 1/m where y = mx + c (here y is vertical position of track above plane and x is coord in u or v plane), and intercept is -c/b (this is crossover point in u/v plane



  Double_t this_slope;
  
  Int_t NoEntries = 1e4;
  Int_t NoTEntries = T->GetEntries();
  
  if(NoEntries > NoTEntries){
    NoEntries = NoTEntries;
  }

  
  for( i = 0; i < NoEntries; i++ ){
    T->GetEntry(i);
    if( (i%5000)==0 ) { cout << "Entry " << i << endl; }
    
   
    if( ntr == 1 && TTD_func::passtrg(Int_t(evttype), trg) && cer_sum > Cer_cut && ps_e/(1e3*tr_p[0]) > Ps_cut && (ps_e+sh_e)/(1e3*tr_p[0]) > Ps_Sh_cut_l &&  (ps_e+sh_e)/(1e3*tr_p[0]) < Ps_Sh_cut_h ){
      for( j = 0; j < NPLANE; j++ ){
	this_slope = d_th[0]*cos(ang[j]*DEG_TO_RAD) 
	  + d_ph[0]*sin(ang[j]*DEG_TO_RAD);

	Bool_t NewClust = kTRUE;
	std::vector<Double_t> Clust_times; // vector to hold custer hits times
	std::vector<Double_t> Clust_tantheta; // vector to hold custer tanetheta
	std::vector<Double_t> Clust_tantheta_alt; // vector to hold custer tantheta alt
	std::vector<Double_t> Clust_dist; // vector to hold custer dist
	std::vector<Double_t> Clust_wires; // vector to hold wire numbers
	

	Int_t Pivot = 0; // var to hold pivot number
	Int_t BegWire = 0; // var to hold first wire
	Int_t EndWire = 0; // var to hold end wire	    



	
	if(nclust[j] == 1 // one cluster
	   && clsiz[j][0] == 5 // cluster size is 5
	   )
	  {
	
	    for( hit = 0; hit < nhit[j]  && nent[j] < MAX_ENTRIES; hit++ ){


	      
	      
	      if( 0 < hittime[j][hit] 
		  //		  &&  hittime[j][hit]  < 260.0e-9
		  //		  && hittrdist[j][hit]< 0.015
		  //		  && hittrknum[j][hit] == 1
		  && wireno[j][hit] >= clbeg[j][0]
		  && wireno[j][hit] <= clend[j][0]
		  ){


		Clust_times.push_back(hittime[j][hit]);
		Clust_tantheta.push_back(this_slope);
		Clust_tantheta_alt.push_back(slope[j][0]);
		Clust_dist.push_back(hittrdist[j][hit]);		
		
		Clust_wires.push_back(wireno[j][hit]);
		
		if(NewClust){

		  
		  Pivot = clpivot[j][0];
		  BegWire = clbeg[j][0];
		  EndWire = clend[j][0];

		  NewClust = kFALSE;
		}
		
		
	      }
	    }
	    
	  }

	// check that clust vector passed conditions

      
	if(Clust_wires.size() == 5){
	  if(Clust_wires[0] >= BegWire && Clust_wires[5] <= EndWire){

	    ClustCount[j]++;
	    CentralWireNo[j].push_back(Pivot);
	    BegWireNo[j].push_back(BegWire);
	    EndWireNo[j].push_back(EndWire);
	    
	    for(Int_t k = 0; k<5; k++){
	      wtime[j].push_back(Clust_times[k]);
	      tanTh[j].push_back(Clust_tantheta[k]);
	      tanTh_alt[j].push_back(Clust_tantheta_alt[k]);
	      trdist[j].push_back(Clust_dist[k]);
	      WireNumbers[j].push_back(Clust_wires[k]);
	      nent[j]++;
	    }
	    
	    
	  }
	}
	
      }
    }
  }
  
  
  // set-up variables to read unnormalised TTD DB into

  Int_t NBins[NPLANE] = {0}; // number of bins lookup table for each plane

  Double_t Low[NPLANE] = {0}; // lowest value for time
  
  std::vector<Double_t> LTable[NPLANE]; // velocity lookup table
  TTDTable* PlaneTable[NPLANE];

  for(Int_t i = 0; i<NPLANE; i++){
    std::ifstream table_DB(Form("DB/lookup_tables/db_%s_%s_lookup_TTD_unnormalised.vdc.%d.dat",arm,plane[i],runnumber));

    cout << "Reading " << Form("DB/lookup_tables/db_%s_%s_lookup_TTD_unnormalised.vdc.%d.dat",arm,plane[i],runnumber) << endl;
    
    Int_t line_no = 0;
  
    std::string line;


    char* phrase = NULL;
    
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
	  PlaneTable[i] = new TTDTable(LTable[i],Low[i]);	  

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




  

  // set-up histograms
  
  
  TH1F *hslope[NPLANE];
  TH1F *hslope_alt[NPLANE];

  // plot wire numbers versus events

  TH2I *hWireClust[NPLANE]; //
  TH2I *hWireGoodClust[NPLANE]; //


  // plot calculated distance difference veruss analyzer distance difference
  // select for clusters with 5 wires where pivot is central wire
  // should be able to use slope of track, and knwon seperation of wires
  // to calculate difference of vertical distance between 1st and 2nd wires (same with 4th and 5th)
  // Plot against slope

  TH2D *hW1W2Diff[NPLANE];
  TH2D *hW4W5Diff[NPLANE];


  // plot calculated TTD table normalisation constant agaisnt angle/ slope

  TH2D *hW1W2TTDSlope[NPLANE];
  TH2D *hW4W5TTDSlope[NPLANE];

  TH1D *hW1W2TTD[NPLANE];
  TH1D *hW4W5TTD[NPLANE];

  TH1D *hW1W2Dist[NPLANE];
  TH1D *hW4W5Dist[NPLANE];
  
  TH2D *hW1W2DistSlope[NPLANE];
  TH2D *hW4W5DistSlope[NPLANE];

  // plot calculated difference of distance against slope calculation

  TH2D* hW1W2DistMSlope[NPLANE];
  TH2D* hW4W5DistMSlope[NPLANE];
  
  for( i = 0; i < NPLANE; i++ ){

    Int_t PlaneCount = 0;
    
    

    hslope[i] = new TH1F(Form("hslope_%s",plane[i]), Form("slope distribution %s",plane[i]), 100,0.5,1.0);
    hslope[i]->GetXaxis()->SetTitle("Slope");
    hslope[i]->GetXaxis()->CenterTitle();

    hslope_alt[i] = (TH1F*) hslope[i]->Clone(Form("hslope_%s",plane[i]));


    hWireClust[i] = new TH2I(Form("hWireClust_%s", plane[i]), Form("Wire Clusters %s", plane[i]), ClustCount[i], 0 , ClustCount[i], 450,0, 450);
    hWireClust[i]->GetXaxis()->SetTitle("Cluster Number");
    hWireClust[i]->GetYaxis()->SetTitle("Wire Number");
    

    hWireGoodClust[i] = (TH2I*) hWireClust[i]->Clone(Form("hWireGoodClust_%s", plane[i]));


    hW1W2Diff[i] = new TH2D(Form("hW1W2Diff_%s", plane[i]), Form("Vertical Dstance Difference %s", plane[i]),100,1,2,100,0.5e-2,0.82e-2);
    hW1W2Diff[i]->GetXaxis()->SetTitle("Analyzer distance difference");
    hW1W2Diff[i]->GetYaxis()->SetTitle("Calculated distance difference");
    
    hW4W5Diff[i] = (TH2D*) hW1W2Diff[i]->Clone(Form("hW4W5Diff_%s", plane[i]));



    // calc limits of correction based on distancece and TTD table entries
    Double_t norm_max = LTable[i].back()/0.013;    
    

    hW1W2TTDSlope[i] = new TH2D(Form("hW1W2TTDSlope_%s", plane[i]),Form("(Wires 1 & 2) TTD Normsalistion vs angle %s", plane[i]),100,1,2,100,norm_max/2,norm_max);
    hW1W2TTDSlope[i]->GetXaxis()->SetTitle("Slope");
    hW1W2TTDSlope[i]->GetYaxis()->SetTitle("TTD normalisation");


    hW4W5TTDSlope[i] = (TH2D*) hW1W2TTDSlope[i]->Clone(Form("hW4W5TTDSlope_%s", plane[i]));
    hW4W5TTDSlope[i]->SetTitle(Form("(Wires 4 & 5) TTD Normsalistion vs angle %s", plane[i]));
    
    
    hW1W2TTD[i] = new TH1D(Form("hW1W2TTD_%s", plane[i]), Form("(Wires 1 & 2) TTD Normsalistion %s", plane[i]),100,norm_max/2,norm_max);
    hW1W2TTD[i]->GetXaxis()->SetTitle("TTD normalisation");
    
    hW4W5TTD[i] = (TH1D*) hW1W2TTD[i]->Clone(Form("hW1W2TTD_%s", plane[i]));
    hW4W5TTD[i]->SetTitle( Form("(Wires 4 & 5) TTD Normsalistion %s", plane[i]));
    

    hW1W2Dist[i] = new TH1D(Form("hW1W2Dist_%s", plane[i]), Form("(Wires 1 & 2) Dist %s", plane[i]),100,0,3e-7);
    hW1W2Dist[i]->GetXaxis()->SetTitle("Time Diff");
    
    hW4W5Dist[i] = (TH1D*) hW1W2Dist[i]->Clone(Form("hW1W2Dist_%s", plane[i]));
    hW4W5Dist[i]->SetTitle( Form("(Wires 4 & 5) Dist %s", plane[i]));
    
    hW1W2DistSlope[i] = new TH2D(Form("hW1W2DistSlope_%s", plane[i]), Form("(Wires 1 & 2) Dist v slope %s", plane[i]),100,1,2,100,0,3e-7);
    hW1W2DistSlope[i]->GetXaxis()->SetTitle("Slope");
    hW1W2DistSlope[i]->GetYaxis()->SetTitle("Time Diff");

    hW4W5DistSlope[i] = (TH2D*) hW1W2DistSlope[i]->Clone(Form("hW4W5DistSlope_%s", plane[i]));
    hW4W5DistSlope[i]->SetTitle(Form("(Wires 4 & 5) Dist v slope %s", plane[i]));


    hW1W2DistMSlope[i] = new TH2D(Form("hW1W2DistMSlope_%s", plane[i]), Form("(Wires 1 & 2) Dist v Mslope %s", plane[i]),100,0,0.015,100,3e4,5e4);

    hW1W2DistMSlope[i]->GetXaxis()->SetTitle("m * wire_sep");
    hW1W2DistMSlope[i]->GetYaxis()->SetTitle("Dist Diff");

    hW4W5DistMSlope[i] =  (TH2D*) hW1W2DistMSlope[i]->Clone(Form("hW4W5DistMSlope_%s", plane[i]));
    hW4W5DistMSlope[i]->SetTitle(Form("(Wires 1 & 2) Dist v Mslope %s", plane[i]));
    
    for( j = 0; j <ClustCount[i]; j++ ){

      
      Bool_t GoodClust = kTRUE;


      
      for(Int_t k = 0; k<5; k++){

	hWireClust[i]->Fill(j,WireNumbers[i][j*5 + k]);
	
	// check that wires are consecutive (no gaps in cluster)
	if(k>0){
	  if(! (WireNumbers[i][j*5 + k] == (WireNumbers[i][j*5 + k - 1])+1) ){
	    GoodClust = kFALSE;
	    break;
	  }	  
	}

	// check that pivot wire in cluster ins 3rd wire (5 wire clusters, checking that pivot is central wire)
	if(k==2){	  
	  if(!(WireNumbers[i][j*5 + k] == CentralWireNo[i][j])){
	    GoodClust = kFALSE;
	    break;
	  }
	  
	}
      }
      
      if(GoodClust){
	
	  for(Int_t k = 0; k<5; k++){      
	    hWireGoodClust[i]->Fill(j,WireNumbers[i][j*5 + k]);	    
	  }

	  // difference in vertical ditance between 1st and 2nd wires
	  // AnaDiff from analyzer
	  // Calc from calulation with angle and wire seperation (should be equal)
	  Double_t AnaDiff_1_2 = -(trdist[i][j*5 + 1] - trdist[i][j*5 + 0]);
	  Double_t CalcDiff_1_2 = (1/tanTh_alt[i][j*5+1]) * (wire_sep);

	  Double_t AnaDiff_4_5 = trdist[i][j*5 + 4] - trdist[i][j*5 + 3];
	  Double_t CalcDiff_4_5 = (1/tanTh_alt[i][j*5+3]) * (wire_sep);
	  
	  
	  hW1W2Diff[i]->Fill(1/tanTh_alt[i][j*5+3],CalcDiff_1_2);
	  hW4W5Diff[i]->Fill(1/tanTh_alt[i][j*5+3],CalcDiff_4_5);


	  // calculate Difference of TTD lookup table divided by calculated seperation between wires	  	  
	  Double_t TTDTableDiff_1_2 = PlaneTable[i]->Convert(wtime[i][j*5+1])-PlaneTable[i]->Convert(wtime[i][j*5+0]);
	  
	  Double_t TTDW1W2Norm = -TTDTableDiff_1_2/CalcDiff_1_2;
	  

	  Double_t TTDTableDiff_4_5 = PlaneTable[i]->Convert(wtime[i][j*5+4])-PlaneTable[i]->Convert(wtime[i][j*5+3]);
		  
	  Double_t TTDW4W5Norm = TTDTableDiff_4_5/CalcDiff_4_5;
	  
	  hW1W2TTDSlope[i]->Fill(1/tanTh_alt[i][j*5+3],TTDW1W2Norm);	  	  
	  hW1W2TTD[i]->Fill(TTDW1W2Norm);


	  hW4W5TTDSlope[i]->Fill(1/tanTh_alt[i][j*5+3],TTDW4W5Norm);	  	  
	  hW4W5TTD[i]->Fill(TTDW4W5Norm);



	  hW1W2Dist[i]->Fill(-wtime[i][j*5+1] + wtime[i][j*5+0]);

	  hW4W5Dist[i]->Fill(wtime[i][j*5+4] - wtime[i][j*5+3]);

	  //	  CalcDiff_4_5
	  


	  hW1W2DistSlope[i]->Fill(1/tanTh_alt[i][j*5+3],-wtime[i][j*5+1] + wtime[i][j*5+0]);
	  hW4W5DistSlope[i]->Fill(1/tanTh_alt[i][j*5+3],wtime[i][j*5+4] - wtime[i][j*5+3]);


	  hW1W2DistMSlope[i]->Fill(CalcDiff_1_2,-TTDTableDiff_1_2);
	  hW4W5DistMSlope[i]->Fill(CalcDiff_4_5,TTDTableDiff_4_5);
	  
      }

    }

    
  }
  




  
  


  TCanvas* cwire[NPLANE];

  TCanvas* cdist[NPLANE];

  TCanvas* cnorm[NPLANE];


  Double_t TTDNormMeans[NPLANE];

  for( i = 0; i < NPLANE; i++ ){

    
    cwire[i] = new TCanvas(Form("cwire_%s",plane[i]),Form("Wire distib %s",plane[i]));
    cwire[i]->Divide(2,1);
    
    cwire[i]->cd(1);
    hWireClust[i]->Draw("colz");


    cwire[i]->cd(2);
    hWireGoodClust[i]->Draw("colz");

    cdist[i] = new TCanvas(Form("cdist_%s",plane[i]),Form("Dist distib %s",plane[i]));
    cdist[i]->Divide(2,1);
    
    cdist[i]->cd(1);
    hW1W2Diff[i]->Draw("colz");

    cdist[i]->cd(2);
    hW4W5Diff[i]->Draw("colz");


    cnorm[i] = new TCanvas(Form("cnorm_%s",plane[i]),Form("Normalisation %s",plane[i]));
    cnorm[i]->Divide(4,3);

    cnorm[i]->cd(1);
    hW1W2TTDSlope[i]->Draw("colz");

    cnorm[i]->cd(2);
    hW1W2TTD[i]->Draw();
    TTDNormMeans[i] = hW1W2TTD[i]->GetMean();
    

    cnorm[i]->cd(3);
    hW4W5TTDSlope[i]->Draw("colz");

    cnorm[i]->cd(4);
    hW4W5TTD[i]->Draw();
    TTDNormMeans[i] += hW4W5TTD[i]->GetMean();

    TTDNormMeans[i] /= 2.0;

    cnorm[i]->cd(5);
    hW1W2Dist[i]->Draw();

    cnorm[i]->cd(6);
    hW4W5Dist[i]->Draw();

    cnorm[i]->cd(7);
    hW1W2DistSlope[i]->Draw("colz");

    cnorm[i]->cd(8);
    hW4W5DistSlope[i]->Draw("colz");

    cnorm[i]->cd(9);
    hW1W2DistMSlope[i]->Draw("colz");

    cnorm[i]->cd(10);
    hW4W5DistMSlope[i]->Draw("colz");

    
  }



  // gather means of distributions and use to normalise TTD tables

  const int kBUFLEN = 150;
  // cout<<"Do you want to rebuild the database using these values? [y/n] ";
  // input[0] = '\0';
  // fgets(input, kBUFLEN, stdin);
  
  for( i = 0; i < NPLANE; i++ ){
    cout << "Table for " << plane[i] << endl << endl;
    for(Int_t j=0; j<NBins[i]; j++){
      LTable[i][j] /= TTDNormMeans[i];
      if (j%10 == 0 && j>0){
  	cout << endl;
      }
      cout << LTable[i][j] << " ";
      
    }
    cout << endl << endl;
  }

  
  // ask whether to replace the values in the database
  char input[kBUFLEN];
  cout<<"Do you want to rebuild the database using these values? [y/n] ";
  input[0] = '\0';
  fgets(input, kBUFLEN, stdin);

  if(input[0] != 'y') {
    cout<<"Exiting without rebuilding database."<<endl;
    // goto cleanup;
    exit(1);
  }

  cout<<"Rebuilding database..."<<endl;


  for( i = 0; i < NPLANE; i++ ){
    if(TTD_func::SaveNewTTDData(LTable[i], NBins[i], Low[i], arm, plane[i], runnumber, "norm"))
      cout<<"Done for "<<plane[i]<<endl;
    else
      cout<<"Failed for "<<plane[i]<<endl;
    
  }
  
}


