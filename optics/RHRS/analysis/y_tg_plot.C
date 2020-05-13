void y_tg_plot(){

  //Macro makes plots to analyze the new theta and phi after optimization

  TString order = "5th";    //Optimization order
  //example for range is -10_10 for -10 cm < x_fp <10 cm
  //use "full" for full focal plane range
  TString range = "full";   //Range in focal plane
  bool before = true;      //Are we doing before optimization plots
  bool make_plots = true;
  bool V_wires = true;

  TString rootfiles = "/home/sean/Grad/Research/APEX/Rootfiles/";
  gStyle->SetPalette(1);
  TFile* f;


  TChain* t = new TChain("T");
  
  if(before) {
    t->Add(rootfiles + "apex_4647.root");
    t->Add(rootfiles + "apex_4648.root");
    t->Add(rootfiles + "apex_4650.root");
  }
  else {
    t->Add(rootfiles + "apex_4647_opt_5th_xfp_full_V_wires.root");
    t->Add(rootfiles + "apex_4648_opt_5th_xfp_full_V_wires.root");
    t->Add(rootfiles + "apex_4650_opt_5th_xfp_full_V_wires.root");
  }
  
  

  //Cuts made for all the plots
  TCut GeneralCut;
  //TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && R.s0.nthit==1  && abs(R.tr.tg_dp) < 0.01";
  if(range == "full") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500)";
  if(range == "-10_10") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && abs(R.tr.r_x) < 0.10";
  if(range == "-45_-25") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > -0.45 && R.tr.r_x < -0.25)";
  if(range == "25_45") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > 0.25 && R.tr.r_x < 0.45)";
  if(range == "-50_-30") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > -0.50 && R.tr.r_x < -0.30)";
  if(range == "-30_-10") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > -0.30 && R.tr.r_x < -0.10)";
  if(range == "30_50") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > 0.30 && R.tr.r_x < 0.50)";
  if(range == "10_30") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > 0.10 && R.tr.r_x < 0.30)";

  TH1F * y_tg = new TH1F("Tg y", "", 400, -40, 40);
  
  //Draw tg y plot
  TCanvas* c = new TCanvas("c","c",1000,1000);
  t->Draw("R.tr.tg_y*1000>>Tg y",GeneralCut,"");
  if(before) y_tg->SetTitle("Tg y Before Opt;Tg y (mm);");
  else y_tg->SetTitle("Tg y "+order+" Order Opt;Tg y (mm);");

   
  /////Add labels for the run number and cuts ////
  TPaveText *pt1 = new TPaveText(0.12,0.78,0.32,0.89,"nbNDC");
  pt1->AddText("Vertical Wires");
  pt1->AddText("Cerenkov signal sum > 500");
  pt1->AddText("Single track");
  if(range == "-10_10") pt1->AddText("|x_{fp}| < 0.10 m");
  if(range == "-45_-25") pt1->AddText("-0.45 m < x_{fp} < -0.25 m");
  if(range == "25_45") pt1->AddText("0.25 m < x_{fp} < 0.45 m");
  pt1->SetFillColor(0);

  TText *text = pt1->GetLineWith("Vertical");
  text->SetTextColor(kRed);
  text->SetTextFont(23);
  text->SetTextSize(23);
  pt1->Draw("same");

  ///// Set Stats to just show number of events/////
  gStyle->SetOptStat(10);
  TPaveStats *s = (TPaveStats*) gPad->GetPrimitive("stats");
  s->SetX2NDC(0.9);
  s->SetY2NDC(0.9);
  s->SetX1NDC(0.72);
  s->SetY1NDC(0.85);

  
  if(make_plots){
    if(before) c->SaveAs("plots/vertex/y_tg_before_opt_xfp_"+range+".gif");
    else c->SaveAs("plots/vertex/y_tg_opt_"+order+"_xfp_"+range+".gif");
  }

  
  TH1F * z_r = new TH1F("Z React", "", 400, -350, 300);

  //Draw z react plot
  TCanvas* c2 = new TCanvas("c2","c2",1000,1000);
  t->Draw("R.tr.vz*1000>>Z React",GeneralCut,"");
  if(before) z_r->SetTitle("Z React Before Opt;Z (mm);");
  else z_r->SetTitle("Z React "+order+" Order Opt;Z(mm);");

  pt1->Draw("same");


  if(make_plots){
    if(before) c2->SaveAs("plots/vertex/z_r_before_opt_xfp_"+range+".gif");
    else c2->SaveAs("plots/vertex/z_r_opt_"+order+"_xfp_"+range+".gif");
  }
}
