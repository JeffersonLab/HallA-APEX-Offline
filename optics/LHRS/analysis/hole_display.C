#include "APEX_Sieve.h"
#include "InputAPEXL.h"
#include "../correlations/Load_new_replay.C"
#include "Load_more_rootfiles.C"

void hole_display(Int_t foil_no = 0){

  //  Int_t FoilID = 1;
  Int_t FoilID = foil_no;

  gStyle->SetOptStat(0);
  

  
  // TCut  GeneralSieveCut ="L.tr.n==1 && L.tr.chi2<0.003 && abs(L.tr.x)<0.75 && abs(L.tr.y)<0.55 && abs(L.tr.th)<0.15 && abs(L.tr.ph)<0.045 && abs(L.tr.tg_dp)<0.01";

  // TCut fid_cut = "abs(L.tr.r_x)<0.1 && abs(th_tgt)<0.03 && abs(ph_tgt)<0.02 ";

  //  TCut  GenrealCut = TCut("L.tr.n==1 && L.tr.tg_dp>-0.01 && L.tr.tg_dp<0.01") && GeneralSieveCut;

   // TCut PID_cuts = "(L.prl1.e/(L.gold.p*1000))>0.2 && ((L.prl1.e+L.prl2.e)/(L.gold.p*1000))>0.51 &&  L.cer.asum_c >400";


  
  TCut GenrealCut = GeneralSieveCut + PID_cuts + FP_cuts;
  // TCut GenrealCut = GeneralSieveCut + PID_cuts;

  TString CutFileName = *SoureRootFile + ".FullCut.root";

  //  TChain* T = Load_new_replay("V1_V2_V3_TPY",Run_number);
  TChain* T = Load_more_rootfiles(Run_number);

  gStyle->SetPalette(1);

  TFile *tcuts = new TFile(CutFileName, "READ");


  // tcuts->GetObject(Form("e_hcut_L_%d_%d_%d",FoilID,n_col,n_row), g[Get_Hole(n_col,n_row)]);


  cout << "CutFileName = " << CutFileName << endl;


  //  TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && R.s0.nthit==1  && abs(R.tr.tg_dp) < 0.01";
  double phi_rms, th_rms;
  int  stat;

  //  TH2F* xpyp = new TH2F("xpyp","",400,-0.065,0.065,400,-0.065,0.065);
  TH2F* xpyp = new TH2F("xpyp","",400,-65,65,400,-65,65);
  //  T->Draw("L.tr.tg_th:L.tr.tg_ph>>xpyp", GenrealCut,"");
  T->Draw("th_tgt*1000:ph_tgt*1000>>xpyp", GenrealCut,"");

      
  TH1F *htemp = (TH1F*)gPad->GetPrimitive("xpyp");   
  xpyp->SetTitle("Tg x' vs y'");
  // xpyp->GetXaxis()->SetTitle("Tg #phi (rad)");
  // xpyp->GetYaxis()->SetTitle("Tg #theta (rad)");
  xpyp->GetXaxis()->SetTitle("#phi_{tg} [mrad]");
  xpyp->GetYaxis()->SetTitle("#theta_{tg} [mrad]");
  
  xpyp->GetYaxis()->SetTitleOffset(1.4);
  xpyp->GetZaxis()->SetRangeUser(0,30);

  //  TImage* img = TImage::Open("Sieve.png");
  TImage* img = TImage::Open("../diagrams/LHRS_Sieve_annotated.png");

  double idx_00 = 0.9035, idy_00 = 0.8735;
  double dx = 0.0505/2, dy = 0.103/2;

  //TCanvas *c = new TCanvas("c","",1000,1000);
  TCanvas *c = new TCanvas("c","Side by Side",1300,800);
  c->Divide(2,1);
  c->cd(1);
  xpyp->Draw("colz");
  c->cd(2);
  img->Draw();

  c->Update();


  TCanvas *c_b = new TCanvas("c_b","target cuts",1300,800);

  c_b->cd(1);
  xpyp->Draw("colz");

  // create second canvas to show FP cuts also

  
  TH2D* thfp_v_yfp = new TH2D("thfp_v_yfp","thfp_v_yfp",300,-0.05,0.05,300,-35,30);

  thfp_v_yfp->GetYaxis()->SetTitleOffset(1.0);
  thfp_v_yfp->GetXaxis()->SetTitleSize(0.05);
  thfp_v_yfp->GetYaxis()->SetTitleSize(0.05);
  thfp_v_yfp->GetXaxis()->SetTitle("y_{FP} [m]");
  thfp_v_yfp->GetYaxis()->SetTitle("th_{FP} [mrad]");



  // Canvas showing target angle picture with cuts

  
  TCanvas *c2 = new TCanvas("c2","Sieve Angles",1000,1000);
  
  TH2F* thtg_v_phtg = new TH2F("thtg_v_phtg","Sieve Angles", 300, -0.04, 0.04, 300, -0.08, 0.08);

  T->Draw("L.tr.tg_th:L.tr.tg_ph>>thtg_v_phtg", GenrealCut,"colz");

  thtg_v_phtg->GetXaxis()->SetTitle("#phi_{tg} [rad]");
  thtg_v_phtg->GetYaxis()->SetTitle("#theta_{tg} [rad]");
  

  // Canvas showing target angle picture without cuts
  
  
  TCanvas *c3 = new TCanvas("c3","Sieve Angles (no ellipses)",1000,1000);

  TH2D* thtg_v_phtg_h = new TH2D("thtg_v_phtg_h","Sieve Angles", 300, -0.04, 0.04, 300, -0.08, 0.08);

  thtg_v_phtg_h->GetXaxis()->SetTitle("#phi_{tg} [rad]");
  thtg_v_phtg_h->GetYaxis()->SetTitle("#theta_{tg} [rad]");
  
  T->Draw("L.tr.tg_th:L.tr.tg_ph>>thtg_v_phtg_h", GenrealCut,"colz");
  

  
  // Canvas showing target sieve xy picture

  TCanvas *c4 = new TCanvas("c4","Sieve x-y",1000,1000);


  //    TH2F* h3 = new TH2F("h3", Form("Sieve plot for Foil #%d",FoilID), 300, -0.04, 0.04, 300, -0.08, 0.08);
  
  TH2D* x_sieve_v_y_sieve = new TH2D("x_sieve_v_y_sieve","Sieve X-Y",300, -0.04, 0.04, 300, -0.08, 0.08);

  x_sieve_v_y_sieve->GetXaxis()->SetTitle("y_sieve [m]");
  x_sieve_v_y_sieve->GetYaxis()->SetTitle("x_sieve [m]");

    
  T->Draw("x_sieve:y_sieve>>x_sieve_v_y_sieve", GenrealCut, "COLZ");
  

  // canvas showing FP cuts


  TCanvas *c5 = new TCanvas("c5","FP cuts",1000,1000);

  
  
  T->Draw("1000*L.tr.r_th:L.tr.r_y>>thfp_v_yfp",GenrealCut,"colz");

  c5->Update();




  //  Int_t FoilID = 1;


  TCutG* g[NSieveRow*NSieveCol];
  TEllipse* Ellipse[NSieveRow*NSieveCol];
  TCutG* g_FP[NSieveRow*NSieveCol];


  for(Int_t n_row = 0; n_row < NSieveRow; n_row++){
    for(Int_t n_col = 0; n_col < NSieveCol; n_col++){


      


  // while (loop){
  
  // int n_col;
  // int n_row;
  // string answer;



  
  //for(int n_col = 0; n_col < 27; n_col++){
  //for(int n_row = 0; n_row < 17; n_row++){
      g[Get_Hole(n_col,n_row)] = NULL;
      //tcuts->GetObject(Form("hcut_R_1_%d_%d",n_col,n_row), g);

      cout << " hole selection = " << Form("e_hcut_L_%d_%d_%d",FoilID,n_col,n_row) << endl;
      tcuts->GetObject(Form("e_hcut_L_%d_%d_%d",FoilID,n_col,n_row), g[Get_Hole(n_col,n_row)]);
      
      if (!g[Get_Hole(n_col,n_row)]){
	cout<<"Hole "<<n_col<<":"<<n_row<<" does not exist"<<endl;
	continue;
      }
      
      cout<<"Showing col:row = "<<n_col<<":"<<n_row<<endl;

      c->cd(1);
      if(g[Get_Hole(n_col,n_row)]){

	for(int i = 0; i<g[Get_Hole(n_col,n_row)]->GetN();i++){
	  g[Get_Hole(n_col,n_row)]->GetX()[i] *= 1000;
	  g[Get_Hole(n_col,n_row)]->GetY()[i] *= 1000;
	}             

	g[Get_Hole(n_col,n_row)]->SetLineWidth(2);
	g[Get_Hole(n_col,n_row)]->Draw("same");
	c_b->cd(0);
	g[Get_Hole(n_col,n_row)]->Draw("same");
	c2->cd(0);
	g[Get_Hole(n_col,n_row)]->Draw("same");
      }
      
      
      // next if conditions deal with end columns (25 and 26) which are not spread out evenly as columns 0-24 are

      Int_t x_col = 0;
      if(n_col == 25){
	x_col = 26;       
      }
      else if(n_col == 26){
	x_col = 28;
      }
      else{
	x_col = n_col;
      }
      


      Ellipse[Get_Hole(n_col,n_row)] = new TEllipse(idx_00-x_col*dx,idy_00-n_row*dy,0.02,0.02);
      Ellipse[Get_Hole(n_col,n_row)]->SetFillStyle(0);
      Ellipse[Get_Hole(n_col,n_row)]->SetLineColor(kRed);
      Ellipse[Get_Hole(n_col,n_row)]->SetLineWidth(2);
      c->cd(2);
      Ellipse[Get_Hole(n_col,n_row)]->Draw("same");


      c->Update();
      c_b->Update();
      c2->Update();
      
      
      // cin.get();
            
      //      Ellipse->Delete();
      // if(g){
      // g->Delete();
      // }
      


      //      TCutG* g_FP = NULL;

      g_FP[Get_Hole(n_col,n_row)] = NULL;

      tcuts->GetObject(Form("FPhcut_L_%d_%d_%d",FoilID,n_col,n_row), g_FP[Get_Hole(n_col,n_row)]);

      if(!g_FP[Get_Hole(n_col,n_row)]){
	cout<<"Hole FP cut "<<n_col<<":"<<n_row<<" does not exist"<<endl;
	cout <<  Form("FPhcut_L_%d_%d_%d",FoilID,n_col,n_row) << endl;
      }
      
      
      c5->cd();
      if(g_FP[Get_Hole(n_col,n_row)]){

	g_FP[Get_Hole(n_col,n_row)]->SetLineColor(kMagenta);	
	g_FP[Get_Hole(n_col,n_row)]->SetLineWidth(2);
	g_FP[Get_Hole(n_col,n_row)]->Draw("same");
	cout << "Showing FP hole: " << n_col << ":" << n_row << endl;
      }


      c5->Update();


      //}
//}


//      cin.get();
            
      // if(g_FP){
      // 	g_FP->Delete();
      // }
      


    }
  }


  c2->Update();


  gSystem->Exec(Form("mkdir output/%d",Run_number));
    
   // gSystem->Exec("mkdir /home/johnw/public_html/correlation_plots/" + DB_name + "/hole_cuts_used");
  


  
  //  c->Print(Form("sieve_angles_w_diagram.pdf",Run_number));

  c2->Print(Form("output/%d/sieve_angles_ellipse.pdf",Run_number));

  c3->Print(Form("output/%d/sieve_angles.pdf",Run_number));

  c4->Print(Form("output/%d/sieve_XY.pdf",Run_number));

  c5->Print(Form("output/%d/sieve_FP.pdf",Run_number));
    
    
    

  
}
