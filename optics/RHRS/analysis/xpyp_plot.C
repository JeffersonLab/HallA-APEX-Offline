void xpyp_plot(){

  //Macro makes plots to analyze the new theta and phi after optimization
  
  TString order = "5th";
  TString range = "all";
  bool before = true; //Are we doing before optimization plots
  bool make_plots = false;

  TString rootfiles = "/home/sean/Grad/Research/APEX/Rootfiles/";
  gStyle->SetPalette(1);
  TFile* f;
  
  if(before) f = new TFile(rootfiles + "apex_4647.root","read");
  else if(range == "all") f = new TFile(rootfiles + "apex_4647_opt_"+order+"_xfp_full.root","read");
  else f = new TFile(rootfiles + "apex_4647_opt_"+order+"_xfp_"+range+".root","read");
  //else f = new TFile(rootfiles + "apex_4647_opt_sieve_plane_"+order+".root","read");
  TTree* t;
  f->GetObject("T",t);
  
  TH2D * xpyp = new TH2D("Tg angles", "", 400, -65, 65, 400, -65, 65);

  TCut GeneralCut;
  //TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && R.s0.nthit==1  && abs(R.tr.tg_dp) < 0.01";
  if(range == "all") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500)";
  if(range == "-10_10") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && abs(R.tr.r_x) < 0.10";
  if(range == "-45_-25") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > -0.45 && R.tr.r_x < -0.25)";
  if(range == "25_45") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > 0.25 && R.tr.r_x < 0.45)";
  if(range == "full") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && ((abs(R.tr.r_x) < 0.10) || (R.tr.r_x > 0.25 && R.tr.r_x < 0.45) || (R.tr.r_x > -0.45 && R.tr.r_x < -0.25))";
  //GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && abs(R.tr.r_x) < 0.03";

  TCanvas* c = new TCanvas("c","c",1000,1000);
  t->Draw("R.tr.tg_th*1000:R.tr.tg_ph*1000>>Tg angles",GeneralCut,"colz");

  if(range == "full" || range == "all") xpyp->GetZaxis()->SetRangeUser(0,100);
  else xpyp->GetZaxis()->SetRangeUser(0,50);
  if(before) xpyp->SetTitle("Tg x' vs y' Before Opt;Tg #phi (mrad);Tg #theta (mrad)");
  else xpyp->SetTitle("Tg x' vs y' "+order+" Order Opt;Tg #phi (mrad);Tg #theta (mrad)");


  TPaveText *pt1 = new TPaveText(0.12,0.78,0.32,0.89,"nbNDC");
  pt1->AddText("Run 4647");
  pt1->AddText("Cerenkov signal sum > 500");
  pt1->AddText("Single track");
  if(range == "-10_10" || range == "full") pt1->AddText("|x_{fp}| < 0.10 m");
  if(range == "-45_-25" || range == "full") pt1->AddText("-0.45 m < x_{fp} < -0.25 m");
  if(range == "25_45" || range == "full") pt1->AddText("0.25 m < x_{fp} < 0.45 m");
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
  
  double sieve_ph[27], sieve_th[17];

  for(int i = 0; i<27; i++){
    if(i < 25) sieve_ph[i] = (14-i)*2.99 - 8.64;
    else sieve_ph[i] = sieve_ph[24] - (i-24)*2.99*2;
  }
  
  for(int j = 0; j<17; j++){
    sieve_th[j] = (j-8)*7.254 + 9.89;
  }
  
  TLine *l[27];
  TLine *l2[17];
  
  
  for(int i=4;i<20;i++){
      l[i] = new TLine(sieve_ph[i], -65, sieve_ph[i], 65);
      l[i]->SetLineColor(2);
      l[i]->Draw("same");
    }
 
  for(int i=3;i<14;i++){
    l2[i] = new TLine(-65, sieve_th[i], 65, sieve_th[i]);
    l2[i]->SetLineColor(2);
    l2[i]->Draw("same");
  }  
  
  TLegend *leg = new TLegend (0.72, 0.8, 0.9, 0.85);
  leg->AddEntry(l[10], "Exp Hole Positions", "l");
  leg->Draw("same");

  if(make_plots){
    if(before) c->SaveAs("plots/tg_th_ph/xpyp_before_opt_xfp_"+range+".gif");
    else c->SaveAs("plots/tg_th_ph/xpyp_opt_"+order+"_xfp_"+range+".gif");
  }
  
  TH1D * y = new TH1D("y", "", 300, -3, 2);
  
  TCanvas *c3 = new TCanvas("c3","",800,600);
  t->Draw("Sieve.y*100>>y",GeneralCut && "(-0.15 < Sieve.x*100) && (0.15 >  Sieve.x*100)","hist");
  y->Scale(1.0/y->GetBinContent(y->GetMaximumBin()));
  y->SetTitle(order+" Order Opt y Projection");
  y->GetXaxis()->SetTitle("Sieve y (cm)");
  y->GetYaxis()->SetTitle("Entries/Max");

  TPaveText *pt2 = new TPaveText(0.12,0.78,0.32,0.89,"nbNDC");
  pt2->AddText("Run 4647");
  pt2->AddText("Cerenkov signal sum > 500");
  pt2->AddText("Single track");
  if(range == "-10_10" || range == "full") pt2->AddText("|x_{fp}| < 0.10 m");
  if(range == "-45_-25" || range == "full") pt2->AddText("-0.45 m < x_{fp} < -0.25 m");
  if(range == "25_45" || range == "full") pt2->AddText("0.25 m < x_{fp} < 0.45 m");
  pt2->AddText("-15 cm < x < 15 cm");
  pt2->SetY1(0.70);
  pt2->SetX2(0.35);
  pt2->SetFillColor(0);

  TText *text2 = pt2->GetLineWith("Run");
  text2->SetTextColor(kRed);
  text2->SetTextFont(23);
  text2->SetTextSize(23);
  
  pt2->Draw("same");
  
  if(make_plots) c3->SaveAs("plots/y_projection_opt_"+order+"_xfp_"+range+".gif");
  
  /*
  TH1D * x_fp = new TH1D("x_fp", "", 500, -0.7, 0.7);
  
  TCanvas *c2 = new TCanvas("c2","",800,600);
  t->Draw("R.tr.r_x>>x_fp",GeneralCut);

  x_fp->SetTitle("Focal Plane x");
  x_fp->GetXaxis()->SetTitle("x_{fp} (m)");
  x_fp->GetYaxis()->SetTitle("Entries");
  x_fp->GetYaxis()->SetRangeUser(0,3500);

  TPaveText *pt2 = new TPaveText(0.12,0.78,0.32,0.89,"nbNDC");
  pt2->AddText("Run 4647");
  pt2->AddText("Cerenkov signal sum > 500");
  pt2->AddText("Single track");
  //pt2->AddText("Single hit in scintillator");  
  pt2->SetFillColor(0);

  TText *text2 = pt2->GetLineWith("Run");
  text2->SetTextColor(kRed);
  text2->SetTextFont(23);
  text2->SetTextSize(23);
  pt2->Draw("same");

  //c2->SaveAs("plots/x_fp.gif");
  */
  /*
  TH2D * xy = new TH2D("xy", "", 500, 6, 13, 500, 10, 17);

  TCanvas *c4 = new TCanvas("c4","",800,600);
  t->Draw("Rrb.y*1000:Rrb.x*1000>>xy",GeneralCut,"colz");

  xy->SetTitle("Raster Scan New Calib");
  xy->GetXaxis()->SetTitle("x (mm)");
  xy->GetYaxis()->SetTitle("y (mm)");
  */


  double sieve_y[27], sieve_x[17];

  for(int i = 0; i<27; i++){
    if(i < 25) sieve_y[i] = (14-i)*0.238 - 1.094;
    else sieve_y[i] = sieve_y[24] - (i-24)*0.238*2;
  }
  
  for(int j = 0; j<17; j++){
    sieve_x[j] = -(j-8)*0.575 - 0.0155;
  }

  TCanvas *c5 = new TCanvas("c5","",1000,1000);
  TH2D * xy_sieve = new TH2D("xy_sieve", "", 400, -5, 5, 400, -6, 4);
  t->Draw("Sieve.x*100:Sieve.y*100>>xy_sieve",GeneralCut,"colz");
  if(range == "full" || range == "all") xy_sieve->GetZaxis()->SetRangeUser(0,100);
  else xy_sieve->GetZaxis()->SetRangeUser(0,50);
  if(before) xy_sieve->SetTitle("Sieve Plane Before Opt;y (cm);x (cm)");
  else xy_sieve->SetTitle("Sieve Plane "+order+" Order Opt;y (cm);x (cm)");

  gStyle->SetOptStat(10);
  TPaveStats *s1 = (TPaveStats*) gPad->GetPrimitive("stats");
  s1->SetX2NDC(0.9);
  s1->SetY2NDC(0.9);
  s1->SetX1NDC(0.72);
  s1->SetY1NDC(0.85);
  

  TLine* lx1 = new TLine(0.75,-6,0.75,4);
  TLine* lx2 = new TLine(0.85,-6,0.85,4);
  TLine* ly1 = new TLine(-5,-0.15,5,-0.15);
  TLine* ly2 = new TLine(-5,.15,5,0.15);

  lx1->SetLineColor(2);
  lx2->SetLineColor(2);
  ly1->SetLineColor(2);
  ly2->SetLineColor(2);
  

  lx1->Draw("same");
  lx2->Draw("same");
  ly1->Draw("same");
  ly2->Draw("same");

  TLegend *leg2 = new TLegend (0.72, 0.8, 0.9, 0.85);
  leg2->AddEntry(lx1, "Column and Row Cuts", "l");
  leg2->Draw("same");


  TLine *l3[27];
  TLine *l4[17];
  /*
  for(int i=4;i<20;i++){
      l3[i] = new TLine(sieve_y[i], -6, sieve_y[i], 4);
      l3[i]->SetLineColor(2);
      l3[i]->Draw("same");
    }
 
  for(int i=3;i<14;i++){
    l4[i] = new TLine(-5, sieve_x[i], 5, sieve_x[i]);
    l4[i]->SetLineColor(2);
    l4[i]->Draw("same");
  }  
  */
  
  
  //leg->Draw("same");
  leg2->Draw("same");
  pt1->Draw("same");

  if(make_plots){
    if(before) c5->SaveAs("plots/sieve_xy/xy_before_opt_xfp_"+range+".gif");
    else c5->SaveAs("plots/sieve_xy/xy_opt_"+order+"_xfp_"+range+".gif");
  }
}
