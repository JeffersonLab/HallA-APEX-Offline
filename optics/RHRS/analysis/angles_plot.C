void angles_plot(){

  //Macro makes plots to analyze the new theta and phi after optimization with optional cuts////

  Int_t run;
  TString after_opt;
  
  cout <<"\n: Please enter a Run Number (-1 to exit):";
  cin >> run;

  cout <<"\n: After optimization? (yes or no):";
  cin >> after_opt;

  TString rootfiles = "/home/sean/Grad/Research/APEX/Rootfiles/";
  gStyle->SetPalette(1);

  TChain *t = new TChain("T");
  TString name;
  
  if(after_opt == "no") name = rootfiles + Form("apex_%d",run);
  else name = rootfiles + Form("apex_%d_opt_5th_xfp_full_V_wires",run);  //Full x_fp with 5 matrices
  
  t->Add(name + ".root");

  //Cuts made for all the plots
  TCut GeneralCut;
  //TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && R.s0.nthit==1  && abs(R.tr.tg_dp) < 0.01";
  GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500)";


  gStyle->SetOptStat(10);

  
  TH1D * x_fp = new TH1D("x_fp", "", 500, -0.7, 0.7);
  
  TCanvas *c1 = new TCanvas("c1","",800,600);
  t->Draw("R.tr.r_x>>x_fp",GeneralCut);

  x_fp->SetTitle("Focal Plane x");
  x_fp->GetXaxis()->SetTitle("x_{fp} (m)");
  x_fp->GetYaxis()->SetTitle("Entries");
  x_fp->GetYaxis()->SetRangeUser(0,3500);


  TPaveText *pt2 = new TPaveText(0.12,0.78,0.32,0.89,"nbNDC");
  pt2->AddText(Form("Run %d",run));
  pt2->AddText("Cerenkov signal sum > 500");
  pt2->AddText("Single track");
  pt2->SetFillColor(0);

  TText *text2 = pt2->GetLineWith("Run");
  text2->SetTextColor(kRed);
  text2->SetTextFont(23);
  text2->SetTextSize(23);
  pt2->Draw("same");


 
  TH2D * xy = new TH2D("xy", "", 500, -3.5, 3.5, 500, 1, 4);

  TCanvas *c2 = new TCanvas("c2","",800,600);
  if(run == 4652 || run == 4653) t->Draw("Rurb.y*1000:Rurb.x*1000>>xy",GeneralCut,"colz");
  else t->Draw("Rrb.y*1000:Rrb.x*1000>>xy",GeneralCut,"colz");
  
  xy->SetTitle("Raster");
  xy->GetXaxis()->SetTitle("x (mm)");
  xy->GetYaxis()->SetTitle("y (mm)");



  
  //Draw tg theta and phi plots
  TH2D * xpyp = new TH2D("Tg angles", "", 400, -65, 65, 400, -65, 65);
  
  TCanvas* c3 = new TCanvas("c3","c3",1000,1000);
  t->Draw("R.tr.tg_th*1000:R.tr.tg_ph*1000>>Tg angles",GeneralCut,"colz");
  //xpyp->GetZaxis()->SetRangeUser(0,150);
  if(after_opt == "no") xpyp->SetTitle("Tg x' vs y' Before Opt;Tg #phi (mrad);Tg #theta (mrad)");
  else xpyp->SetTitle("Tg x' vs y' 5th Order Opt;Tg #phi (mrad);Tg #theta (mrad)");

  pt2->Draw("same");


  c1->Update();
  c2->Update();
  c3->Update();
  
  //// Now we make cuts and check the sieve again ////

  TString use_xfp;
  TString use_raster_x;
  TString use_raster_y;

  double xfp_min = 0;
  double xfp_max = 0;
  double raster_x_min = 0;
  double raster_y_min = 0;
  double raster_x_max = 0;
  double raster_y_max = 0;

  TCut cut_xfp = "";
  TCut cut_raster_x = "";
  TCut cut_raster_y = "";
  TCut cut_theta = "";
  
  cout <<"\n: Make cuts on x_fp? (yes or no):";
  cin >> use_xfp;


  if(use_xfp == "yes"){

    cout <<"\n: x_fp minimum value (in meters):";
    cin >> xfp_min;

    cout <<"\n: x_fp maximum value (in meters):";
    cin >> xfp_max;

    cut_xfp = Form("R.tr.r_x > %f && R.tr.r_x < %f",xfp_min,xfp_max);

    t->Draw("R.tr.r_x>>x_fp",GeneralCut && cut_xfp);

    if(run == 4652 || run == 4653) t->Draw("Rurb.y*1000:Rurb.x*1000>>xy",GeneralCut && cut_xfp,"colz");
    else t->Draw("Rrb.y*1000:Rrb.x*1000>>xy",GeneralCut && cut_xfp,"colz");

    t->Draw("R.tr.tg_th*1000:R.tr.tg_ph*1000>>Tg angles",GeneralCut && cut_xfp,"colz");

    pt2->AddText(Form("%.2f m < x_{fp} < %.2f m",xfp_min,xfp_max));
    pt2->Draw("same");

    c1->Modified();
    c2->Modified();
    c3->Modified();
    
    c1->Update();
    c2->Update();
    c3->Update();

  }


  
  cout <<"\n: Make cuts on raster x? (yes or no):";
  cin >> use_raster_x;


  if(use_raster_x == "yes"){

    cout <<"\n: raster x minimum value (in mm):";
    cin >> raster_x_min;

    cout <<"\n: raster x maximum value (in mm):";
    cin >> raster_x_max;

    if(run == 4652 || run == 4653) cut_raster_x = Form("Rurb.x > %f/1000 && Rurb.x < %f/1000",raster_x_min,raster_x_max);
    else cut_raster_x = Form("Rrb.x > %f/1000 && Rrb.x < %f/1000",raster_x_min,raster_x_max);

    t->Draw("R.tr.r_x>>x_fp",GeneralCut && cut_xfp && cut_raster_x);

    if(run == 4652 || run == 4653) t->Draw("Rurb.y*1000:Rurb.x*1000>>xy",GeneralCut && cut_xfp && cut_raster_x,"colz");
    else t->Draw("Rrb.y*1000:Rrb.x*1000>>xy",GeneralCut && cut_xfp && cut_raster_x,"colz");

    t->Draw("R.tr.tg_th*1000:R.tr.tg_ph*1000>>Tg angles",GeneralCut && cut_xfp && cut_raster_x,"colz");

    pt2->AddText(Form("%.2f mm < raster x < %.2f mm",raster_x_min,raster_x_max));
    pt2->Draw("same");

    c1->Modified();
    c2->Modified();
    c3->Modified();
    
    c1->Update();
    c2->Update();
    c3->Update();
    
  }


  cout <<"\n: Make cuts on raster y? (yes or no):";
  cin >> use_raster_y;

  if(use_raster_y == "yes"){

    cout <<"\n: raster y minimum value (in mm):";
    cin >> raster_y_min;

    cout <<"\n: raster y maximum value (in mm):";
    cin >> raster_y_max;

    if(run == 4652 || run == 4653) cut_raster_y = Form("Rurb.y > %f/1000 && Rurb.y < %f/1000",raster_y_min,raster_y_max);
    else cut_raster_y = Form("Rrb.y > %f/1000 && Rrb.y < %f/1000",raster_y_min,raster_y_max);

    t->Draw("R.tr.r_x>>x_fp",GeneralCut && cut_xfp && cut_raster_x && cut_raster_y);

    if(run == 4652 || run == 4653) t->Draw("Rurb.y*1000:Rurb.x*1000>>xy",GeneralCut && cut_xfp && cut_raster_x && cut_raster_y,"colz");
    else t->Draw("Rrb.y*1000:Rrb.x*1000>>xy",GeneralCut && cut_xfp && cut_raster_x && cut_raster_y,"colz");

    t->Draw("R.tr.tg_th*1000:R.tr.tg_ph*1000>>Tg angles",GeneralCut && cut_xfp && cut_raster_x && cut_raster_y,"colz");

    pt2->AddText(Form("%.2f mm < raster x < %.2f mm",raster_y_min,raster_y_max));
    pt2->Draw("same");

    c1->Modified();
    c2->Modified();
    c3->Modified();
    
    c1->Update();
    c2->Update();
    c3->Update();
    
  }
  


  //// Make Projection Cuts/////

  TString use_proj;
  double theta_min;
  double theta_max;
  

  cout <<"\n: Make phi projection? (yes or no):";
  cin >> use_proj;

  if(use_proj == "yes"){

    cout <<"\n: Theta minimum value (in mrad):";
    cin >> theta_min;

    cout <<"\n: Theta maximum value (in mrad):";
    cin >> theta_max;

    cut_theta = Form("R.tr.tg_th > %f/1000 && R.tr.tg_th < %f/1000",theta_min,theta_max);

    TLine* l1 = new TLine(-65,theta_min,65,theta_min);
    TLine* l2 = new TLine(-65,theta_max,65,theta_max);

    l1->SetLineColor(2);
    l2->SetLineColor(2);
    
    c3->cd();

    l1->Draw("same");
    l2->Draw("same");
    
    c3->Modified();
    c3->Update();

    TH1D* proj = new TH1D("projection","Phi Projection;#phi_{tg} (mrad)",300,-65,65);

    TCanvas *c4 = new TCanvas("c4","",800,600);
    t->Draw("R.tr.tg_ph*1000>>projection",GeneralCut && cut_xfp && cut_raster_x && cut_raster_y && cut_theta,"colz");

    
  }
  

  
}
