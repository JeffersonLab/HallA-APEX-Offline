void correlation(){

  //Analyzes how the optimized matrix th and phi are correlated with focal plane variables

  TString run = "4647";     //Run number
  TString order = "5th";
  TString range = "-10_10";
  bool before = false; //Are we doing before optimization plots
  bool make_plots = false;
  bool brute_force = false;
  
  TString rootfiles = "/home/sean/Grad/Research/APEX/Rootfiles/";
  gStyle->SetPalette(1);
  TFile* f;

  if(before) f = new TFile(rootfiles + "apex_"+run+".root","read");
  else if(range == "full" && !brute_force) f = new TFile(rootfiles + "apex_"+run+"_opt_"+order+"_xfp_full.root","read");       //Full x_fp with single matrix
  else if(range == "full" && brute_force) f = new TFile(rootfiles + "apex_"+run+"_opt_"+order+"_xfp_full_brute.root","read");  //Full x_fp with 5 matrices
  else f = new TFile(rootfiles + "apex_"+run+"_opt_"+order+"_xfp_"+range+".root","read");
  
  TTree* t;
  f->GetObject("T",t);

  
  TCut GeneralCut;
  //TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && R.s0.nthit==1  && abs(R.tr.tg_dp) < 0.01";
  if(range == "full") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500)";
  if(range == "-10_10") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && abs(R.tr.r_x) < 0.10";
  if(range == "-45_-25") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > -0.45 && R.tr.r_x < -0.25)";
  if(range == "25_45") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > 0.25 && R.tr.r_x < 0.45)";
  


  TPaveText *pt1 = new TPaveText(0.12,0.25,0.40,0.49,"nbNDC");
  pt1->AddText("Run 4647");
  pt1->AddText("Cerenkov signal sum > 500");
  pt1->AddText("Single track");
  if(range == "-10_10") pt1->AddText("|x_{fp}| < 0.10 m");
  if(range == "-45_-25") pt1->AddText("-0.45 m < x_{fp} < -0.25 m");
  if(range == "25_45") pt1->AddText("0.25 m < x_{fp} < 0.45 m");
  //pt1->AddText("-15 cm < sieve x < 15 cm");
  pt1->SetFillColor(0);

  TText *text = pt1->GetLineWith("Run");
  text->SetTextColor(kRed);
  text->SetTextFont(23);
  text->SetTextSize(23);
 
  TPaveText *pt2 = new TPaveText(0.12,0.65,0.40,0.89,"nbNDC");
  pt2->AddText("Run 4647");
  pt2->AddText("Cerenkov signal sum > 500");
  pt2->AddText("Single track");
  if(range == "-10_10") pt2->AddText("|x_{fp}| < 0.10 m");
  if(range == "-45_-25") pt2->AddText("-0.45 m < x_{fp} < -0.25 m");
  if(range == "25_45") pt2->AddText("0.25 m < x_{fp} < 0.45 m");
  //pt2->AddText("-75 cm < sieve y < 85 cm");
  pt2->SetFillColor(0);

  TText *text2 = pt2->GetLineWith("Run");
  text2->SetTextColor(kRed);
  text2->SetTextFont(23);
  text2->SetTextSize(23);

  TPaveText *pt3 = new TPaveText(0.65,0.74,0.90,0.89,"nbNDC");
  pt3->AddText("Vertical Band");
  pt3->SetFillColor(0);

  TPaveText *pt4 = new TPaveText(0.65,0.74,0.90,0.89,"nbNDC");
  pt4->AddText("Horizontal Band");
  pt4->SetFillColor(0);
  
  
  TCut cut_phi = "(-0.15 < Sieve.x*100) && (0.15 >  Sieve.x*100)";
  TCut cut_theta = "(0.75 < Sieve.y*100) && (0.85 >  Sieve.y*100)";

  double th_min = -40, th_max = 60; // -40, 60
  double ph_min = -30, ph_max = 30; // -30, 30
   
  TH2D * th_th = new TH2D("th_th", "", 400, -30, 30, 400, th_min, th_max);
  TH2D * th_ph = new TH2D("th_ph", "", 400, -50, 50, 400, th_min, th_max);
  TH2D * th_x = new TH2D("th_x", "", 400, -60, 60, 400, th_min, th_max);
  TH2D * th_y = new TH2D("th_y", "", 400, -6, 6, 400, th_min, th_max);

  TH2D * ph_th = new TH2D("ph_th", "", 400, -30, 30, 400, ph_min, ph_max);
  TH2D * ph_ph = new TH2D("ph_ph", "", 400, -50, 50, 400, ph_min, ph_max);
  TH2D * ph_x = new TH2D("ph_x", "", 400, -60, 60, 400, ph_min, ph_max);
  TH2D * ph_y = new TH2D("ph_y", "", 400, -6, 6, 400, ph_min, ph_max);

  TH2D * th_xraster = new TH2D("th_xraster", "", 200, -5, -3, 200, th_min, th_max);
  TH2D * th_yraster = new TH2D("th_yraster", "", 200, 6, 10, 200, th_min, th_max);
  TH2D * ph_xraster = new TH2D("ph_xraster", "", 200, -5, -3, 200, ph_min, ph_max);
  TH2D * ph_yraster = new TH2D("ph_yraster", "", 200, 6, 10, 200, ph_min, ph_max); 
  
  TString title;

  title = order+" Order Opt Projection";

  th_th->SetTitle(title+";#theta_{fp} (mrad);#theta_{tg} (mrad)");
  th_ph->SetTitle(title+";#phi_{fp} (mrad);#theta_{tg} (mrad)");
  th_x->SetTitle(title+";x_{fp} (cm);#theta_{tg} (mrad)");
  th_y->SetTitle(title+";y_{fp} (cm);#theta_{tg} (mrad)");

  ph_th->SetTitle(title+";#theta_{fp} (mrad);#phi_{tg} (mrad)");
  ph_ph->SetTitle(title+";#phi_{fp} (mrad);#phi_{tg} (mrad)");
  ph_x->SetTitle(title+";x_{fp} (cm);#phi_{tg} (mrad)");
  ph_y->SetTitle(title+";y_{fp} (cm);#phi_{tg} (mrad)");

  th_xraster->SetTitle(title+";Beam x (mm);#theta_{tg} (mrad)");
  th_yraster->SetTitle(title+";Beam y (mm);#theta_{tg} (mrad)");
  ph_xraster->SetTitle(title+";Beam x (mm);#phi_{tg} (mrad)");
  ph_yraster->SetTitle(title+";Beam y (mm);#phi_{tg} (mrad)");
  
  gStyle->SetOptStat(10);
  //Make tg theta plots
  TCanvas* c = new TCanvas("c","c",1200,1200);
  c->Divide(2,2);
  c->cd(1);
  t->Draw("R.tr.tg_th*1000:R.tr.r_th*1000>>th_th",GeneralCut && cut_theta,"colz");
  pt3->Draw("same");
  
  /*
  TPaveStats *s = (TPaveStats*) gPad->GetPrimitive("stats");
  s->SetX2NDC(0.9);
  s->SetY2NDC(0.9);
  s->SetX1NDC(0.72);
  s->SetY1NDC(0.85);
  */
  
  c->cd(2);
  t->Draw("R.tr.tg_th*1000:R.tr.r_ph*1000>>th_ph",GeneralCut && cut_theta,"colz");
  pt2->Draw("same");

  c->cd(3);
  t->Draw("R.tr.tg_th*1000:R.tr.r_x*100>>th_x",GeneralCut && cut_theta,"colz");

  c->cd(4);
  t->Draw("R.tr.tg_th*1000:R.tr.r_y*100>>th_y",GeneralCut && cut_theta,"colz");

  c->Update();
  
  //Make tg phi plots
  TCanvas* c2 = new TCanvas("c2","c2",1200,1200);
  c2->Divide(2,2);
  c2->cd(1);
  t->Draw("R.tr.tg_ph*1000:R.tr.r_th*1000>>ph_th",GeneralCut && cut_phi,"colz");
  
  
  c2->cd(2);
  t->Draw("R.tr.tg_ph*1000:R.tr.r_ph*1000>>ph_ph",GeneralCut && cut_phi,"colz");
  pt1->Draw("same");

  c2->cd(3);
  t->Draw("R.tr.tg_ph*1000:R.tr.r_x*100>>ph_x",GeneralCut && cut_phi,"colz");

  c2->cd(4);
  t->Draw("R.tr.tg_ph*1000:R.tr.r_y*100>>ph_y",GeneralCut && cut_phi,"colz");
  pt4->Draw("same");

  TCanvas* c3 = new TCanvas("c3","c3",1200,1200);
  c3->Divide(2,2);

  c3->cd(1);
  t->Draw("R.tr.tg_th*1000:Rrb.x*1000>>th_xraster",GeneralCut && cut_theta,"colz");
  pt3->Draw("same");

  c3->cd(2);
  t->Draw("R.tr.tg_th*1000:Rrb.y*1000>>th_yraster",GeneralCut && cut_theta,"colz");
  pt2->Draw("same");
  pt3->Draw("same");

  c3->cd(3);
  t->Draw("R.tr.tg_ph*1000:Rrb.x*1000>>ph_xraster",GeneralCut && cut_phi,"colz");
  pt4->Draw("same");
  
  c3->cd(4);
  t->Draw("R.tr.tg_ph*1000:Rrb.y*1000>>ph_yraster",GeneralCut && cut_phi,"colz");
  pt4->Draw("same");
  
  c->Update();
  c2->Update();
  c3->Update();
  
  
  if(make_plots){
    if(before) c->SaveAs("plots/correlations/"+run+"/th_tg_before_opt_xfp_"+range+".gif");
    else c->SaveAs("plots/correlations/"+run+"/th_tg_opt_"+order+"_xfp_"+range+".gif");

    if(before) c2->SaveAs("plots/correlations/"+run+"/ph_tg_before_opt_xfp_"+range+".gif");
    else c2->SaveAs("plots/correlations/"+run+"/ph_tg_opt_"+order+"_xfp_"+range+".gif");

    if(before) c3->SaveAs("plots/correlations/"+run+"/tg_before_opt_raster_"+range+".gif");
    else c3->SaveAs("plots/correlations/"+run+"/tg_opt_"+order+"_raster_"+range+".gif");
  }


  
  
  
 
}
