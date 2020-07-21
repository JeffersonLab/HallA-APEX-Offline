void hole_all(){

  // Macro that shows all of the holes and the graphical cuts side by side for a nice visualization

  TString run = "4647";
  TString range = "full";
  
  gStyle->SetPalette(1);
  //TFile* tcuts = new TFile("../Sieve/"+run+"/xfp_"+range+"/apex_"+run+".root.FullCut.root","read");
  TFile* tcuts = new TFile("../Sieve/"+run+"/xfp_full_brute/apex_"+run+"_opt_5th_xfp_full_brute.root.FullCut.root","read");
  TChain * t = new TChain("T");
  //t->Add("/home/sean/Grad/Research/APEX/Rootfiles/apex_"+run+".root");
  t->Add("/home/sean/Grad/Research/APEX/Rootfiles/apex_"+run+"_opt_5th_xfp_full_brute.root");
  
  TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500)";
  if(range == "-50_-30") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > -0.50 && R.tr.r_x < -0.30)";
  if(range == "-10_10") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && abs(R.tr.r_x) < 0.10";
  //TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > -0.45 && R.tr.r_x < -0.25)";
  //TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > 0.25 && R.tr.r_x < 0.45)";
  

  TH2F* xpyp = new TH2F("xpyp","",400,-65,65,400,-65,65);
  t->Draw("R.tr.tg_th*1000:R.tr.tg_ph*1000>>xpyp", GeneralCut,"");

      
  TH1F *htemp = (TH1F*)gPad->GetPrimitive("xpyp");   
  xpyp->SetTitle("Tg th vs Tg ph;#phi (mrad);#theta (mrad)");  
  xpyp->GetYaxis()->SetTitleOffset(1.4);
  xpyp->GetZaxis()->SetRangeUser(0,150);

  TCanvas *c = new TCanvas("c","",1000,1000);
  xpyp->Draw("colz");

  TPaveText *pt1 = new TPaveText(0.12,0.78,0.32,0.89,"nbNDC");
  pt1->AddText("Run "+run);
  pt1->AddText("Cerenkov signal sum > 500");
  pt1->AddText("Single track");
  //pt1->AddText("|x_{fp}| < 0.10 m");
  //pt1->AddText("-0.45 < x_{fp} < -0.25 m");
  //pt1->AddText("0.25 < x_{fp} < 0.45 m");
  pt1->SetFillColor(0);

  TText *text = pt1->GetLineWith("Run");
  text->SetTextColor(kRed);
  text->SetTextFont(23);
  text->SetTextSize(23);
  pt1->Draw("same");
  
  gStyle->SetOptStat(10);
  TPaveStats *s = (TPaveStats*) gPad->GetPrimitive("stats");
  s->SetX2NDC(0.9);
  s->SetY2NDC(0.9);
  s->SetX1NDC(0.72);
  s->SetY1NDC(0.85);
  
  for(int n_col = 0; n_col < 27; n_col++){
    for(int n_row = 0; n_row < 17; n_row++){
      TCutG* g = NULL;
  
      tcuts->GetObject(Form("hcut_R_1_%d_%d",n_col,n_row), g);
      
      if (!g){
	continue;
      }
      
      g->Draw("same");
      
    }
  }

  c->SaveAs("xfp_full_brute/apex_"+run+"_xfp_full_brute_cuts.gif");
  
}
