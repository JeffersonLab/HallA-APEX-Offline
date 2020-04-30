void xpyp_plot(){

  //Macro makes plots to analyze the new theta and phi after optimization

  TString run = "4652";     //Run number
  TString order = "5th";    //Optimization order
  //example for range is -10_10 for -10 cm < x_fp <10 cm
  //use "full" for full focal plane range
  TString range = "full";   //Range in focal plane
  bool before = true;      //Are we doing before optimization plots
  bool make_plots = false;
  bool brute_force = true;  //Are we using brute force method

  TString rootfiles = "/home/sean/Grad/Research/APEX/Rootfiles/";
  gStyle->SetPalette(1);
  TFile* f;
  
  if(before) f = new TFile(rootfiles + "apex_"+run+".root","read");
  else if(range == "full" && !brute_force) f = new TFile(rootfiles + "apex_"+run+"_opt_"+order+"_xfp_full.root","read");       //Full x_fp with single matrix
  else if(range == "full" && brute_force) f = new TFile(rootfiles + "apex_"+run+"_opt_"+order+"_xfp_full_brute.root","read");  //Full x_fp with 5 matrices
  else f = new TFile(rootfiles + "apex_"+run+"_opt_"+order+"_xfp_"+range+".root","read");

  
  TTree* t;
  f->GetObject("T",t);
  
  TH2D * xpyp = new TH2D("Tg angles", "", 400, -65, 65, 400, -65, 65);

  //Cuts made for all the plots
  TCut GeneralCut;
  //TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && R.s0.nthit==1  && abs(R.tr.tg_dp) < 0.01";
  if(range == "full") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500)";
  if(range == "-10_10") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && abs(R.tr.r_x) < 0.10";
  if(range == "-45_-25") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > -0.45 && R.tr.r_x < -0.25)";
  if(range == "25_45") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > 0.25 && R.tr.r_x < 0.45)";
  if(range == "-50_-30") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > -0.30 && R.tr.r_x < -0.10)";


  //Draw tg theta and phi plots
  TCanvas* c = new TCanvas("c","c",1000,1000);
  t->Draw("R.tr.tg_th*1000:R.tr.tg_ph*1000>>Tg angles",GeneralCut,"colz");
  if(range == "full") xpyp->GetZaxis()->SetRangeUser(0,150);
  else xpyp->GetZaxis()->SetRangeUser(0,50);
  if(before) xpyp->SetTitle("Tg x' vs y' Before Opt;Tg #phi (mrad);Tg #theta (mrad)");
  else xpyp->SetTitle("Tg x' vs y' "+order+" Order Opt;Tg #phi (mrad);Tg #theta (mrad)");

  /////Add labels for the run number and cuts ////
  TPaveText *pt1 = new TPaveText(0.12,0.78,0.32,0.89,"nbNDC");
  pt1->AddText("Run "+run);
  pt1->AddText("Cerenkov signal sum > 500");
  pt1->AddText("Single track");
  if(range == "-10_10") pt1->AddText("|x_{fp}| < 0.10 m");
  if(range == "-45_-25") pt1->AddText("-0.45 m < x_{fp} < -0.25 m");
  if(range == "25_45") pt1->AddText("0.25 m < x_{fp} < 0.45 m");
  pt1->SetFillColor(0);

  TText *text = pt1->GetLineWith("Run");
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

  ///// Get expected hole positions from csv file and draw /////
  double sieve_ph[27], sieve_th[17];

  TString csv = range;
  if(brute_force) csv = "full_brute";
  ifstream csv_file("../Sieve/"+run+"/xfp_"+csv+"/apex_"+run+".root.cuts_full.csv");
  
  string line;

  getline(csv_file,line);
  getline(csv_file,line);
  int col_n = 0, row_n = 0;
  
  while (getline(csv_file,line)){
    if(row_n == 17){
      col_n++;
      row_n = 0;
    }
    stringstream linestream(line);
    string cell;

    int n_elem = 1;

    while(getline(linestream,cell,',')){
      if(n_elem == 5) sieve_ph[col_n] = stod(cell);
      else if(n_elem == 7) sieve_th[row_n] = stod(cell);
      n_elem++;
    }
    
    row_n++;
  }


  TLine *l[27];
  TLine *l2[17];
  
  
  for(int i=4;i<20;i++){
      l[i] = new TLine(sieve_ph[i], -65, sieve_ph[i], 65);
      l[i]->SetLineColor(2);
      if(!before) l[i]->Draw("same");
    }
 
  for(int i=3;i<14;i++){
    l2[i] = new TLine(-65, sieve_th[i], 65, sieve_th[i]);
    l2[i]->SetLineColor(2);
    if(!before) l2[i]->Draw("same");
  }  
  
  TLegend *leg = new TLegend (0.72, 0.8, 0.9, 0.85);
  leg->AddEntry(l[10], "Exp Hole Positions", "l");
  if(!before) leg->Draw("same");

  ///// Draw labels for column and row number ///////
  for(int n_col = 4; n_col <=19; n_col++){
    if(n_col%2 != 0) continue;
    TText *text = new TText;
    text -> SetTextFont(1);
    text -> SetTextSize(0.02);
    text->SetTextAlign(22);
    if(!before) text -> DrawText(sieve_ph[n_col], 67, Form("%d",n_col));
  }

  for(int n_row = 3; n_row <=13; n_row++){
    TText *text = new TText;
    text -> SetTextFont(1);
    text -> SetTextSize(0.02);
    text->SetTextAlign(22);
    if(!before) text -> DrawText(-60,sieve_th[n_row], Form("%d",n_row));
  }

  if(make_plots){
    if(before) c->SaveAs("plots/tg_th_ph/"+run+"/xpyp_before_opt_xfp_"+range+".gif");
    else if(brute_force) c->SaveAs("plots/tg_th_ph/"+run+"/xpyp_opt_"+order+"_xfp_full_brute.gif");
    else c->SaveAs("plots/tg_th_ph/"+run+"/xpyp_opt_"+order+"_xfp_"+range+".gif");
  }
  
  //// Plots for x_fp ///
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


  gStyle->SetOptStat(1);
  TH2D * xy = new TH2D("xy", "", 500, -3.5, 3.5, 500, 1, 4);

  TCanvas *c4 = new TCanvas("c4","",800,600);
  if(run == "4652" || run = "4653") t->Draw("Rurb.y*1000:Rurb.x*1000>>xy",GeneralCut,"colz");
  else t->Draw("Rrb.y*1000:Rrb.x*1000>>xy",GeneralCut,"colz");
  
  xy->SetTitle("Raster Scan Calib");
  xy->GetXaxis()->SetTitle("x (mm)");
  xy->GetYaxis()->SetTitle("y (mm)");
  pt1->Draw("same");

  if(make_plots) c4->SaveAs("plots/raster/"+run+"/"+run+"raster.gif");

}
