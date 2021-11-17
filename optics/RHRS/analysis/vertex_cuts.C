void vertex_cuts(){

  /////Draw cut on z vs ph plot and check projections/////

  
  TString rootfiles = "/home/sean/Grad/Research/APEX/Rootfiles/";
  gStyle->SetPalette(1);
  
  TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500)";

  TChain* t = new TChain("T");

   
  t->Add(rootfiles + "apex_4653_opt_3rd_xfp_full_V_Opt_All.root");


  double phi_min = -65;
  double phi_max = 65;
  double z_min = -350;
  double z_max = 350;
  
  
  TH2F* zph = new TH2F("zph", "Z vs. #phi_{tg}; #phi_{tg} (mrad); z (mm)", 400, -40, 40, 400, z_min, z_max);
  TH2F* FP_cut = new TH2F("FP_cut", "y_{FP} vs. #phi_{FP} with cut; #phi_{FP} (mrad);y_{FP} (mm)", 400, -40, 50, 400, -50, 50);
  
  TH2F* xpyp = new TH2F("xpyp", "Angles Without Cuts; #phi_{tg} (mrad); #theta_{tg} (mrad)", 400, phi_min, phi_max, 400, -65, 65);
  TH2F* xpyp_cut = new TH2F("xpyp_cut", "Angles With Cuts; #phi_{tg} (mrad); #theta_{tg} (mrad)", 400, phi_min, phi_max, 400, -65, 65);
  TH2F* xpyp_cut_all = new TH2F("xpyp_cut_all", "Angles With All Cuts; #phi_{tg} (mrad); #theta_{tg} (mrad)", 400, phi_min, phi_max, 400, -65, 65);


  
  zph->GetYaxis()->SetTitleOffset(1.2);
  

  TCanvas *c_zph = new TCanvas("c_zph","",800,800);
  t->Draw("R.tr.vz*1000:R.tr.tg_ph*1000>>zph",GeneralCut,"colz");
  
  zph->SetStats(0);
   
  TCutG* cutg = (TCutG*)(TVirtualPad::Pad()->WaitPrimitive("CUTG", "CutG"));

  c_zph->Update();
  
  cutg->SetName("Vertex_cut"); 
  cutg->SetVarX("R.tr.tg_ph");
  cutg->SetVarY("R.tr.vz");

  cutg->SetLineColor(kMagenta);
  cutg->SetLineWidth(2);
  cutg->Draw("PL");
  c_zph->Update();


  TCutG* cutg2 = new TCutG("test",cutg->GetN());

  
  for(int i = 0; i<cutg->GetN();i++){
    cutg2->SetPoint(i,cutg->GetX()[i],cutg->GetY()[i]);
    cutg->GetX()[i] /= 1000;
    cutg->GetY()[i] /= 1000;
  }
  cutg2->SetLineColor(kMagenta);
  cutg2->SetLineWidth(2);

 
  t->Draw("R.tr.r_ph*1000:R.tr.r_y*1000>>FP_cut",GeneralCut && "Vertex_cut","colz");
  
   
  TCutG* cutg3 = (TCutG*)(TVirtualPad::Pad()->WaitPrimitive("CUTG", "CutG"));

  cutg3->SetName("FP_Vertex_cut"); 
  cutg3->SetVarX("R.tr.r_y");
  cutg3->SetVarY("R.tr.r_ph");

  cutg3->SetLineColor(kMagenta);
  cutg3->SetLineWidth(2);
  cutg3->Draw("PL");
  
  TCutG* cutg4 = new TCutG("test2",cutg3->GetN());

  
  for(int i = 0; i<cutg3->GetN();i++){
    cutg4->SetPoint(i,cutg3->GetX()[i],cutg3->GetY()[i]);
    cutg3->GetX()[i] /= 1000;
    cutg3->GetY()[i] /= 1000;
  }
  cutg4->SetLineColor(kMagenta);
  cutg4->SetLineWidth(2);

  
  TCanvas *c = new TCanvas("c","",1200,1000);
  c->Divide(2,3);


  c->cd(1);
  zph->Draw("colz");

  c->cd(2);
  t->Draw("R.tr.tg_th*1000:R.tr.tg_ph*1000>>xpyp",GeneralCut,"colz");

  c->cd(3);
  zph->Draw("colz");
  //cutg2->Draw("PL same");
  cutg2->Draw("PL same");
  
  
  c->cd(4);
  t->Draw("R.tr.tg_th*1000:R.tr.tg_ph*1000>>xpyp_cut",GeneralCut && "Vertex_cut","colz");

  c->cd(5);
  FP_cut->Draw("colz");
  cutg4->Draw("PL same");

  c->cd(6);
  t->Draw("R.tr.tg_th*1000:R.tr.tg_ph*1000>>xpyp_cut_all",GeneralCut && "Vertex_cut" && "FP_Vertex_cut","colz");
  
}
