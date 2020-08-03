void raster_correlation(){

  TString target = "Horizontal Wire";
  bool save = false;
  TString run;
  
  TString rootfiles = "/home/sean/Grad/Research/APEX/Rootfiles/";
  gStyle->SetPalette(1);
  
  TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500)";

  if(target == "Horizontal Wire") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && abs(R.tr.r_x) < 0.10";
  
  TChain* t = new TChain("T");

  
  if(target == "Horizontal Wire") run = "4657";
  else run = "4647";
  
  if(target == "Horizontal Wire")  t->Add(rootfiles + "apex_"+run+".root");
  else t->Add(rootfiles + "apex_"+run+"_opt_5th_xfp_full_brute.root");

  
  TH2D * xpyp = new TH2D("Tg angles", ";#phi_{tg} (mrad);#theta_{tg} (mrad)", 400, -45, 45, 400, -45, 45);

  
  xpyp->GetYaxis()->SetTitleOffset(1.2);

  
  if(target == "Horizontal Wire"){
    TH2D *yx = new TH2D("yx",";y (mm); x (mm)", 200, -4, 0, 200, 1, 2.5);

    TCanvas *c2 = new TCanvas("c2","",800,800);
    t->Draw("Rrb.y*1000:Rrb.x*1000>>yx",GeneralCut,"colz");
    
    TCutG* cut_ras = (TCutG*)(TVirtualPad::Pad()->WaitPrimitive("CUTG", "CutG"));
    c2->Update();
    cut_ras->SetName("rast_cut");
    cut_ras->SetVarX("Rrb.x");
    cut_ras->SetVarY("Rrb.y");
    
    cut_ras->SetLineColor(kMagenta);
    cut_ras->SetLineWidth(2);
    cut_ras->Draw("PL");
    c2->Update();

    for(int i = 0; i<cut_ras->GetN();i++){
      cut_ras->GetX()[i] /= 1000;
      cut_ras->GetY()[i] /= 1000;
  }

    GeneralCut += "rast_cut";

    c2->Close();
  }
  
  
  TCanvas *c3 = new TCanvas("c3","",800,800);
  t->Draw("R.tr.tg_th*1000:R.tr.tg_ph*1000>>Tg angles",GeneralCut,"colz");
  
  xpyp->SetStats(0);
   
  TCutG* cutg = (TCutG*)(TVirtualPad::Pad()->WaitPrimitive("CUTG", "CutG"));

  c3->Update();
  
  cutg->SetName("hole_cut"); //
  cutg->SetVarX("R.tr.tg_ph");
  cutg->SetVarY("R.tr.tg_th");

  cutg->SetLineColor(kMagenta);
  cutg->SetLineWidth(2);
  cutg->Draw("PL");
  c3->Update();

  
  TCutG* cutg2 = new TCutG("test",cutg->GetN());

  double theta_min = 100;
  double theta_max = -100;
  double phi_min = 100;
  double phi_max = -100;

  double raster_x_min = -2.5;
  double raster_y_min = 1.5;
  double raster_x_max = 0.5;
  double raster_y_max = 4;

  if(target == "Horizontal Wire"){
    raster_x_min = -4;
    raster_y_min = 1.2;
    raster_x_max = 0;
    raster_y_max = 1.6;
  }

  for(int i = 0; i<cutg->GetN();i++){
    cutg2->SetPoint(i,cutg->GetX()[i],cutg->GetY()[i]);

    if(cutg->GetX()[i] < phi_min) phi_min = cutg->GetX()[i];
    if(cutg->GetX()[i] > phi_max) phi_max = cutg->GetX()[i];
    if(cutg->GetY()[i] < theta_min) theta_min = cutg->GetY()[i];
    if(cutg->GetY()[i] > theta_max) theta_max = cutg->GetY()[i];
    
    cutg->GetX()[i] /= 1000;
    cutg->GetY()[i] /= 1000;
  }

  cutg2->SetLineColor(kMagenta);
  cutg2->SetLineWidth(2);

  TH2D * th_x = new TH2D("Th vs x", ";raster x (mm);#theta_{tg} (mrad)", 100, raster_x_min, raster_x_max, 100, theta_min, theta_max);
  TH2D * th_y = new TH2D("Th vs y", ";raster y (mm);#theta_{tg} (mrad)", 100, raster_y_min, raster_y_max, 100, theta_min, theta_max);
  TH2D * ph_x = new TH2D("Ph vs x", ";raster x (mm);#phi_{tg} (mrad)", 100, raster_x_min, raster_x_max, 100, phi_min, phi_max);
  TH2D * ph_y = new TH2D("Ph vs y", ";raster y (mm);#phi_{tg} (mrad)", 100, raster_y_min, raster_y_max, 100, phi_min, phi_max);


  TPaveText *pt = new TPaveText(0.70,0.25,0.90,0.45,"nbNDC");
  pt->AddText(target);
  pt->SetFillColor(0);

  
  TCanvas *c = new TCanvas("c","",1400,1000);
  c->Divide(2,2);

  c->cd(1);
  t->Draw("R.tr.tg_th*1000:Rrb.x*1000>>Th vs x",GeneralCut && "hole_cut","colz");
  pt->Draw("same");
  th_x->Fit("pol1");
  
  c->cd(2);
  t->Draw("R.tr.tg_th*1000:Rrb.y*1000>>Th vs y",GeneralCut && "hole_cut","colz");
  th_y->Fit("pol1");

  c->cd(3);
  t->Draw("R.tr.tg_ph*1000:Rrb.x*1000>>Ph vs x",GeneralCut && "hole_cut","colz");
  ph_x->Fit("pol1");
  
  c->cd(4);
  t->Draw("R.tr.tg_ph*1000:Rrb.y*1000>>Ph vs y",GeneralCut && "hole_cut","colz");
  ph_y->Fit("pol1");

  if(save) c->SaveAs("plots/correlations/"+run+"/raster_correlations.gif");
  
}
