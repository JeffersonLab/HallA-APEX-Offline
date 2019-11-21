void proj(TString DataBase){

  gStyle->SetPalette(1);
  
  //TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && R.s0.nthit==1  && abs(R.tr.tg_dp) < 0.01";
  //TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500)";
  TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && abs(R.tr.r_x) < 0.10";
  //TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > -0.45 && R.tr.r_x < -0.25)";
  //TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > 0.25 && R.tr.r_x < 0.45)";
  //TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && ((abs(R.tr.r_x) < 0.10) || (R.tr.r_x > 0.25 && R.tr.r_x < 0.45) || (R.tr.r_x > -0.45 && R.tr.r_x < -0.25))";

  TString rootfiles = "/home/sean/Grad/Research/APEX/Rootfiles/";
  
  //TFile* f = new TFile(rootfiles + "apex_4647_opt_sieve_plane_3rd.root","open");
  TFile* f = new TFile(rootfiles + "apex_4647.root","open");
  TTree* t;
  f->GetObject("T",t);

  
  TH2D * xpyp = new TH2D("Sieve Plane", "", 400, -5, 5, 400, -6, 4);

  
  opt = new ROpticsOpt();
  opt->LoadRawData(rootfiles + "apex_4647.dat", (UInt_t) - 1, MaxDataPerGroup);
  opt->LoadDataBase(DataBase);
  opt->PrepareSieve();
  //opt->CheckSieve(1);

  double R_tr_n;
  double R_tr_x_fp[100];
  double R_tr_y_fp[100];
  double R_tr_th_fp[100];
  double R_tr_ph_fp[100];
  double R_tr_p[100];
  double R_cer_asum_c;
  double R_s0_nthit;
  double beam_x[100];
  double beam_y[100];
  double R_tr_x_rot[100];
  double R_tr_y_rot[100];
  double R_tr_th_rot[100];
  double R_tr_ph_rot[100];
  double R_tr_vz[100];
  double R_tr_tg_y[100];
  double R_tr_tg_th[100];
  double R_tr_tg_ph[100];
  double R_tr_tg_dp[100];

  
  t->SetBranchStatus("*",0);
  t->SetBranchStatus("R.tr.n",1);
  t->SetBranchStatus("R.tr.x",1);
  t->SetBranchStatus("R.tr.y",1);
  t->SetBranchStatus("R.tr.th",1);
  t->SetBranchStatus("R.tr.ph",1);
  t->SetBranchStatus("R.tr.p",1);
  t->SetBranchStatus("R.cer.asum_c",1);
  t->SetBranchStatus("R.s0.nthit",1);
  t->SetBranchStatus("Rrb.x",1);
  t->SetBranchStatus("Rrb.y",1);
  t->SetBranchStatus("R.tr.r_x",1);
  t->SetBranchStatus("R.tr.r_y",1);
  t->SetBranchStatus("R.tr.r_th",1);
  t->SetBranchStatus("R.tr.r_ph",1);
  t->SetBranchStatus("R.tr.vz",1);
  
  t->SetBranchAddress("R.tr.n",&R_tr_n);
  t->SetBranchAddress("R.tr.x",R_tr_x_fp);
  t->SetBranchAddress("R.tr.y",R_tr_y_fp);
  t->SetBranchAddress("R.tr.th",R_tr_th_fp);
  t->SetBranchAddress("R.tr.ph",R_tr_ph_fp);
  t->SetBranchAddress("R.tr.p",R_tr_p);
  t->SetBranchAddress("R.cer.asum_c",&R_cer_asum_c);
  t->SetBranchAddress("R.s0.nthit",&R_s0_nthit);
  t->SetBranchAddress("Rrb.x",beam_x);
  t->SetBranchAddress("Rrb.y",beam_y);
  t->SetBranchAddress("R.tr.r_x",R_tr_x_rot);
  t->SetBranchAddress("R.tr.r_y",R_tr_y_rot);
  t->SetBranchAddress("R.tr.r_th",R_tr_th_rot);
  t->SetBranchAddress("R.tr.r_ph",R_tr_ph_rot);
  t->SetBranchAddress("R.tr.vz",R_tr_vz);


  int j = 0;
  for(int i=0; i<t->GetEntries(); i++){
  
    t->GetEntry(i);

    double theta = opt->calc_tgth(i);
    double phi = opt->calc_tgph(i);
    
    if(R_tr_n == 1 && R_cer_asum_c > 500 && abs(R_tr_x_fp[0]) < 0.10 ) xpyp->Fill(opt->sieve_y(i)*100,opt->sieve_x(i)*100);
    
    if(R_tr_n == 1 && R_cer_asum_c > 500 && abs(R_tr_x_fp[0]) < 0.10 ){
      if(j<10)cout<<j<<" "<<opt->sieve_x(i)<<" "<<opt->sieve_y(i)<<endl;
      j++;
    }
  }

  TCanvas *c = new TCanvas("c","c",1000,1000);
  xpyp->Draw("colz");
  //xpyp->GetZaxis()->SetRangeUser(0,80);
  xpyp->GetZaxis()->SetRangeUser(0,50);
  xpyp->SetTitle("Sieve Plane 3rd Order Opt");
  xpyp->GetXaxis()->SetTitle("Sieve y (cm)");
  xpyp->GetYaxis()->SetTitle("Sieve x (cm)");


  TPaveText *pt1 = new TPaveText(0.12,0.78,0.32,0.89,"nbNDC");
  pt1->AddText("Run 4647");
  pt1->AddText("Cerenkov signal sum > 500");
  pt1->AddText("Single track");
  pt1->AddText("|x_{fp}| < 0.10 m");
  pt1->AddText("-0.45 m < x_{fp} < -0.25 m");
  pt1->AddText("0.25 m < x_{fp} < 0.45 m");
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
    if(i < 25) sieve_ph[i] = (14-i)*0.238 - 1.094;
    else sieve_ph[i] = sieve_ph[24] - (i-24)*0.238*2;
  }
  
  for(int j = 0; j<17; j++){
    sieve_th[j] = -(j-8)*0.575 - 0.0155;
  }
  
  TLine *l[27];
  TLine *l2[17];
  
  
  for(int i=4;i<20;i++){
      l[i] = new TLine(sieve_ph[i], -5, sieve_ph[i], 5);
      l[i]->SetLineColor(2);
      l[i]->Draw("same");
    }
 
  for(int i=3;i<14;i++){
    l2[i] = new TLine(-5, sieve_th[i], 5, sieve_th[i]);
    l2[i]->SetLineColor(2);
    l2[i]->Draw("same");
  }  

  //c->SaveAs("../plots/sieve_plane_opt_3rd_xfp_full.gif");
}
