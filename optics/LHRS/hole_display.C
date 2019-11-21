#include "APEX_Sieve.h"
#include "InputAPEXL.h"

void hole_display(){


  TCut  GeneralSieveCut ="L.tr.n==1 && L.tr.chi2<0.003 && abs(L.tr.x)<0.75 && abs(L.tr.y)<0.55 && abs(L.tr.th)<0.15 && abs(L.tr.ph)<0.045";

  TCut  GenrealCut = TCut("L.tr.n==1 && L.tr.tg_dp>-0.01 && L.tr.tg_dp<0.01") && GeneralSieveCut;



  TString CutFileName = *SoureRootFile + ".FullCut.root";

  TChain* T = Load_more_rootfiles(Run_number, Run_number_2);

  gStyle->SetPalette(1);

  TFile *tcuts = new TFile(CutFileName, "READ");


  cout << "CutFileName = " << CutFileName << endl;


  //  TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && R.s0.nthit==1  && abs(R.tr.tg_dp) < 0.01";
  double phi_rms, th_rms;
  int  stat;

  TH2F* xpyp = new TH2F("xpyp","",400,-0.065,0.065,400,-0.065,0.065);
  T->Draw("L.tr.tg_th:L.tr.tg_ph>>xpyp", GenrealCut,"");

      
  TH1F *htemp = (TH1F*)gPad->GetPrimitive("xpyp");   
  xpyp->SetTitle("Tg x' vs y'");
  xpyp->GetXaxis()->SetTitle("Tg #phi (rad)");
  xpyp->GetYaxis()->SetTitle("Tg #theta (rad)");
  xpyp->GetYaxis()->SetTitleOffset(1.4);
  xpyp->GetZaxis()->SetRangeUser(0,30);

  //  TImage* img = TImage::Open("Sieve.png");
  TImage* img = TImage::Open("LHRS_Sieve_annotated.png");

  double idx_00 = 0.9035, idy_00 = 0.8735;
  double dx = 0.0505/2, dy = 0.103/2;

  //TCanvas *c = new TCanvas("c","",1000,1000);
  TCanvas *c = new TCanvas("c","",1300,800);
  c->Divide(2,1);
  c->cd(1);
  xpyp->Draw("colz");
  c->cd(2);
  img->Draw();

  c->Update();





  // create second canvas to show FP cuts also

  
  TH2D* thfp_v_yfp = new TH2D("thfp_v_yfp","thfp_v_yfp",300,-0.05,0.05,300,-35,30);

  thfp_v_yfp->GetYaxis()->SetTitleOffset(1.0);
  thfp_v_yfp->GetXaxis()->SetTitleSize(0.05);
  thfp_v_yfp->GetYaxis()->SetTitleSize(0.05);
  thfp_v_yfp->GetXaxis()->SetTitle("y (FP) [m]");
  thfp_v_yfp->GetYaxis()->SetTitle("th (FP) [mrad]");


  


  TCanvas *c2 = new TCanvas("c2","",1300,800);

  T->Draw("L.tr.r_th*1000:L.tr.r_y>>thfp_v_yfp",GenrealCut,"colz");

  c2->Update();




  Int_t FoilID = 0;


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
      
      
      c2->cd();
      if(g_FP[Get_Hole(n_col,n_row)]){

	g_FP[Get_Hole(n_col,n_row)]->SetLineColor(kMagenta);	
	g_FP[Get_Hole(n_col,n_row)]->SetLineWidth(2);
	g_FP[Get_Hole(n_col,n_row)]->Draw("same");
	cout << "Showing FP hole: " << n_col << ":" << n_row << endl;
      }


      c2->Update();


      //}
//}


//      cin.get();
            
      // if(g_FP){
      // 	g_FP->Delete();
      // }
      


    }
  }
}
