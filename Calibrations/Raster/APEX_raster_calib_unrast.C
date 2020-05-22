/*
 * Author: Tyler Hague
 * Date 2 Nov 17
 *
 * This code will calibrate the raster both in the Fastbus and fADCs.
 * Work in progress
 */



#include "TString.h"
#include "Load_more_rootfiles.C"

void APEX_raster_calib_unrast(){

  Int_t run = 0;
  cout << "What run number would you like to calibrate with?    ";
  cin >> run;
  cout << endl << endl;

  TString arm;
  cout << "Which arm (L or R)?    ";
  cin >> arm;
  cout << endl << endl;

  TString unrast_arm;
  TString arm_string;

  arm_string = arm;

  if (arm == "L" || arm == "l"){
    arm = "Lrb";
    unrast_arm = "Lurb";
  }
  else if(arm == "R" || arm == "r"){
    arm = "Rrb";
    unrast_arm = "Rurb";
  }
  else{
    cout << "Choose L or R, you muppet" << endl;
    return 1;
  }


  if(run<=0){
    cout << "Invalid run number. Exiting." << endl << endl;
    return;
  }


  cout << "Run = " << run << " and arm = " << arm << endl;




  //Set Raster correlation factors
  Double_t kx = 1.; //Horrizontal beam direction
  Double_t ky = 1.; //Vertical beam direction

  //Open Root File

  // TChain *rootfile = new TChain("T");

  // int i = 1;
  // //TFile *test_file = new TFile(Form("/volatile/halla/triton/tjhague/rootfiles/coinc_test_%d_%d.root",run,i));

  // if(!gSystem->AccessPathName(TString::Format("/adaqfs/home/a-onl/tritium_work/MARANAL0/replay/Rootfiles/tritium_%d.root",run),kFileExists)){
  //   rootfile->Add(TString::Format("/adaqfs/home/a-onl/tritium_work/MARANAL0/replay/Rootfiles/tritium_%d.root",run));
  //   cout << "Added file: tritium_" << run << ".root" << endl;
  // }else{
  //   cout << "Requested run has not been replayed. Exiting." << endl << endl;
  //   return;
  // }

  // while(!gSystem->AccessPathName(TString::Format("/adaqfs/home/a-onl/tritium_work/MARANAL0/replay/Rootfiles/tritium_%d_%d.root",run,i),kFileExists)){
  //   rootfile->Add(TString::Format("/adaqfs/home/a-onl/tritium_work/MARANAL0/replay/Rootfiles/tritium_%d_%d.root",run,i));
  //   cout << "Added file: tritium_" << run << "_" << i << ".root" << endl;
  //   i=i+1;
  // }                       


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  TChain* T = Load_more_rootfiles(run);





  //Create strings that will get the position of the beam at the target from BPM info
  //Hardcoded numbers are the position of BPMB w.r.t. the target and the distance between the two BPMs

  //  TString  BtoTarg = "2.214"; // BPMB to target distance value for tritium

  TString   BtoTarg = "1.166"; // BPMB to target distance value for APEX
  TString AtoB    = "5.131" ; // BPMA to BPMB distance

  //  TString targ_x_pos = "(" + arm + ".BPMB.x)" + " + ((" + BtoTarg + ") * ((" + arm + ".BPMB.x - " + arm + ".BPMA.x) / (" + AtoB + " )))";

  TString targ_x_pos = "((" + arm + ".BPMB.x)" + " + ((" + BtoTarg + ")/(" + AtoB + ") * (" + arm + ".BPMB.x - " + arm + ".BPMA.x))) ";

  cout << targ_x_pos << endl;

  //  TString targ_x_pos = "(" + arm + ".BPMB.x) + ((2.214) * ((" + arm + ".BPMB.x - " + arm + ".BPMA.x) / (5.131)))";


  TString targ_y_pos = "((" + arm + ".BPMB.y)" + " + ((" + BtoTarg + ")/(" + AtoB + ") * (" + arm + ".BPMB.y - " + arm + ".BPMA.y))) ";

  //  TString targ_y_pos = "(" + arm + ".BPMB.y) + ((2.214) * ((" + arm + ".BPMB.y - " + arm + ".BPMA.y) / (5.131)))";

  cout << " targ_x_pos = " << targ_x_pos << " \n targ_y_pos = " << targ_y_pos << endl;


  // Define Plots of BPM and Raster Current (and Beam Position at Target?)
  // Also get the Mean and RMS of each plot


  Double_t lim_1 = 0;
  Double_t lim_2 = 0;

  
  if(arm == "Lrb"){
    // lim_1 = 25000;
    // lim_2 = 35000;
    lim_1 = 20000;
    lim_2 = 40000;
  }
  else if (arm == "Rrb"){
    lim_1 = 18000;
    lim_2 = 28000;
    if (run > 4500){
      // lim_1 = 35000;
      // lim_2 = 50000;
      lim_1 = 25000;
      lim_2 = 55000;
      
    }
  }
  else{
    cout << "Error, 'arm' variable found to be neither 'Lrb' or 'Rrb'" << endl;

  }

  //Plot X and Y Currents for both Rasters
  TH1F *r1xcurr = new TH1F("r1xcurr", "Raster 1-X Current vs ADC Channel", 1000, lim_1, lim_2);
  TH1F *r1ycurr = new TH1F("r1ycurr", "Raster 1-Y Current vs ADC Channel", 1000, lim_1, lim_2);
  TH1F *r2xcurr = new TH1F("r2xcurr", "Raster 2-X Current vs ADC Channel", 1000, lim_1, lim_2);
  TH1F *r2ycurr = new TH1F("r2ycurr", "Raster 2-Y Current vs ADC Channel", 1000, lim_1, lim_2);

  //Plot X and Y Position for both BPMs
  TH1F *bpmaxpos = new TH1F("bpmaxpos", "BPM A-X Position (m)", 400, -0.005, 0.005);
  TH1F *bpmaypos = new TH1F("bpmaypos", "BPM A-Y Position (m)", 400, -0.005, 0.005);
  TH1F *bpmbxpos = new TH1F("bpmbxpos", "BPM B-X Position (m)", 400, -0.005, 0.005);
  TH1F *bpmbypos = new TH1F("bpmbypos", "BPM B-Y Position (m)", 400, -0.005, 0.005);

  //Plot X and Y Position at the Target
  TH1F *targxpos = new TH1F("targxpos", "Target X Position (based on BPMA and BPMB) (m)", 400, -0.005, 0.005);
  TH1F *targypos = new TH1F("targypos", "Target Y Position (based on BPMA and BPMB) (m)", 400, -0.005, 0.005);

  
  // Unrastered and rastered beam variables from the analyzer
  TH1F *An_targxpos = new TH1F("An_targxpos", "Analyzer Unrastered Target X Position (m)", 400, -0.005, 0.005);
  TH1F *An_targypos = new TH1F("An_targypos", " Analyzer Unrastered Target Y Position (m)", 400, -0.005, 0.005);

  TH1F *An_Rast_targxpos = new TH1F("An_Rast_targxpos", "Analyzer rastered Target X Position (m) (check)", 400, -0.005, 0.005);
  TH1F *An_Rast_targypos = new TH1F("An_Rast_targypos", " Analyzer rastered Target Y Position (m) (check)", 400, -0.005, 0.005);




  //Populate the plots
  //Need to start with a cut for when the beam is off
  //If the BPM goes to a very large (negative?) value, the beam is off

  //Probably add a current cut later? Maybe unnecessary

  // TString cut = "(TMath::Abs(" + arm + ".BPMA.x)<100)&&((ev";
  // if(arm == "Rrb"){
  //   cut += "Right";
  // }else if("Lrb"){
  //   cut += "Left";
  // }
  // cut += "dnew_r*0.0003299)>20)";
  // TCut beamcut = cut.Data();
  TString cut = "(TMath::Abs(" + arm + ".BPMA.x)<100)";
  TCut beamcut = cut.Data();



  cout << "cut for no beam is << " << cut << endl;



  //The plots are added to a canvas as they are populated
  TCanvas *raster_canvas = new TCanvas("raster_canvas",Form("Raster canvas (run-%d ",run) + arm_string + "-arm)");
  raster_canvas->Divide(2,2);

  raster_canvas->cd(1);
  T->Draw(arm + ".Raster.rawcur.x>>r1xcurr",beamcut);
  Double_t r1x_mean = r1xcurr->GetMean();
  Double_t r1x_rms = r1xcurr->GetRMS();

  raster_canvas->cd(2);
  T->Draw(arm + ".Raster.rawcur.y>>r1ycurr",beamcut);
  Double_t r1y_mean = r1ycurr->GetMean();
  Double_t r1y_rms = r1ycurr->GetRMS();


  raster_canvas->cd(3);
  T->Draw(arm + ".Raster2.rawcur.x>>r2xcurr",beamcut);
  Double_t r2x_mean = r2xcurr->GetMean();
  Double_t r2x_rms = r2xcurr->GetRMS();

  raster_canvas->cd(4);
  T->Draw(arm + ".Raster2.rawcur.y>>r2ycurr",beamcut);
  Double_t r2y_mean = r2ycurr->GetMean();
  Double_t r2y_rms = r2ycurr->GetRMS();




  TCanvas *position_canvas = new TCanvas("position_canvas",Form("Position canvas (run-%d ",run) + arm_string + "-arm)");
  position_canvas->Divide(2,5);
  position_canvas->cd(1);
  T->Draw(arm + ".BPMA.x>>bpmaxpos",beamcut);
  Double_t bpmax_mean = bpmaxpos->GetMean();
  Double_t bpmax_rms = bpmaxpos->GetRMS();

  position_canvas->cd(2);
  T->Draw(arm + ".BPMA.y>>bpmaypos",beamcut);
  Double_t bpmay_mean = bpmaypos->GetMean();
  Double_t bpmay_rms = bpmaypos->GetRMS();

  position_canvas->cd(3);
  T->Draw(arm + ".BPMB.x>>bpmbxpos",beamcut);
  Double_t bpmbx_mean = bpmbxpos->GetMean();
  Double_t bpmbx_rms = bpmbxpos->GetRMS();

  position_canvas->cd(4);
  T->Draw(arm + ".BPMB.y>>bpmbypos",beamcut);
  Double_t bpmby_mean = bpmbypos->GetMean();
  Double_t bpmby_rms = bpmbypos->GetRMS();

  position_canvas->cd(5);
  T->Draw(targ_x_pos + ">>targxpos",beamcut);
  Double_t targx_mean = targxpos->GetMean();
  Double_t targx_rms = targxpos->GetRMS();

  position_canvas->cd(6);
  T->Draw(targ_y_pos + ">>targypos",beamcut);
  Double_t targy_mean = targypos->GetMean();
  Double_t targy_rms = targypos->GetRMS();

  
  position_canvas->cd(7);
  T->Draw(unrast_arm + ".x>>An_targxpos",beamcut);
  Double_t An_targx_mean = An_targxpos->GetMean();
  Double_t An_targx_rms = An_targxpos->GetRMS();

  position_canvas->cd(8);
  T->Draw(unrast_arm + ".y>>An_targypos",beamcut);
  Double_t An_targy_mean = An_targypos->GetMean();
  Double_t An_targy_rms = An_targypos->GetRMS();


  position_canvas->cd(9);
  T->Draw(arm + ".x>>An_Rast_targxpos",beamcut);
  // Double_t An_targx_mean = An_targxpos->GetMean();
  // Double_t An_targx_rms = An_targxpos->GetRMS();

  position_canvas->cd(10);
  T->Draw(arm + ".y>>An_Rast_targypos",beamcut);
  // Double_t An_targy_mean = An_targypos->GetMean();
  // Double_t An_targy_rms = An_targypos->GetRMS();



//   position_canvas->cd(6);

  //Calculate and display coeffiecients
  //Give instructions on how to update database?
  Double_t UbpmAx_offset = bpmax_mean - ((r1x_mean*bpmax_rms)/(r1x_rms*kx));
  Double_t UbpmAx_slope = bpmax_rms/(r1x_rms*kx);
  Double_t UbpmAy_offset = bpmay_mean - ((r1y_mean*bpmay_rms)/(r1y_rms*ky));
  Double_t UbpmAy_slope = bpmay_rms/(r1y_rms*ky);

  Double_t UbpmBx_offset = bpmbx_mean - ((r1x_mean*bpmbx_rms)/(r1x_rms*kx));
  Double_t UbpmBx_slope = bpmbx_rms/(r1x_rms*kx);
  Double_t UbpmBy_offset = bpmby_mean - ((r1y_mean*bpmby_rms)/(r1y_rms*ky));
  Double_t UbpmBy_slope = bpmby_rms/(r1y_rms*ky);

  Double_t Utargx_offset = targx_mean - ((r1x_mean*targx_rms)/(r1x_rms*kx));
  Double_t Utargx_slope = targx_rms/(r1x_rms*kx);
  Double_t Utargy_offset = targy_mean - ((r1y_mean*targy_rms)/(r1y_rms*ky));
  Double_t Utargy_slope = targy_rms/(r1y_rms*ky);

  Double_t DbpmAx_offset = bpmax_mean - ((r2x_mean*bpmax_rms)/(r2x_rms*kx));
  Double_t DbpmAx_slope = bpmax_rms/(r2x_rms*kx);
  Double_t DbpmAy_offset = bpmay_mean - ((r2y_mean*bpmay_rms)/(r2y_rms*ky));
  Double_t DbpmAy_slope = bpmay_rms/(r2y_rms*ky);

  Double_t DbpmBx_offset = bpmbx_mean - ((r2x_mean*bpmbx_rms)/(r2x_rms*kx));
  Double_t DbpmBx_slope = bpmbx_rms/(r2x_rms*kx);
  Double_t DbpmBy_offset = bpmby_mean - ((r2y_mean*bpmby_rms)/(r2y_rms*ky));
  Double_t DbpmBy_slope = bpmby_rms/(r2y_rms*ky);

  Double_t Dtargx_offset = targx_mean - ((r2x_mean*targx_rms)/(r2x_rms*kx));
  Double_t Dtargx_slope = targx_rms/(r2x_rms*kx);
  Double_t Dtargy_offset = targy_mean - ((r2y_mean*targy_rms)/(r2y_rms*ky));
  Double_t Dtargy_slope = targy_rms/(r2y_rms*ky);

  cout << arm << ".Raster.raw2posA = " << UbpmAx_offset << " " << UbpmAy_offset << " " << UbpmAx_slope << " " << UbpmAy_slope << " 0.0 0.0" << endl; 
  cout << arm << ".Raster.raw2posB = " << UbpmBx_offset << " " << UbpmBy_offset << " " << UbpmBx_slope << " " << UbpmBy_slope << " 0.0 0.0" << endl; 
  cout << arm << ".Raster.raw2posT = " << Utargx_offset << " " << Utargy_offset << " " << Utargx_slope << " " << Utargy_slope << " 0.0 0.0" << endl << endl; 
  
  cout << arm << ".Raster2.raw2posA = " << DbpmAx_offset << " " << DbpmAy_offset << " " << DbpmAx_slope << " " << DbpmAy_slope << " 0.0 0.0" << endl; 
  cout << arm << ".Raster2.raw2posB = " << DbpmBx_offset << " " << DbpmBy_offset << " " << DbpmBx_slope << " " << DbpmBy_slope << " 0.0 0.0" << endl; 
  cout << arm << ".Raster2.raw2posT = " << Dtargx_offset << " " << Dtargy_offset << " " << Dtargx_slope << " " << Dtargy_slope << " 0.0 0.0" << endl; 
}
