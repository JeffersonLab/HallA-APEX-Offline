void vertex_cuts(){

  /////Draw cut on z vs ph plot and check projections/////

  
  TString rootfiles = "/home/sean/Grad/Research/APEX/Rootfiles/";
  gStyle->SetPalette(1);
  
  TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500)";

  TChain* t = new TChain("T");

   
  t->Add(rootfiles + "apex_4653_opt_5th_xfp_full_V_wires.root");

  
  TH2F* zph = new TH2F("zph", "Z vs. Tg Phi; #phi_{tg} (mm); z (mm)", 400, -40, 40, 400, -350, 200);
  TH1F* z1d = new TH1F("z1d", "Z React;z (mm);", 400, -350, 200);
  TH1F* z1d_cut = new TH1F("z1d_cut", "Z React with cuts;z (mm)", 400, -350, 200);
  zph->GetYaxis()->SetTitleOffset(1.2);
  

  TCanvas *c3 = new TCanvas("c3","",800,800);
  t->Draw("R.tr.vz*1000:R.tr.tg_ph*1000>>zph",GeneralCut,"colz");
  
  zph->SetStats(0);
   
  TCutG* cutg = (TCutG*)(TVirtualPad::Pad()->WaitPrimitive("CUTG", "CutG"));

  c3->Update();
  
  cutg->SetName("Vertex_cut"); //
  cutg->SetVarX("R.tr.tg_ph");
  cutg->SetVarY("R.tr.vz");

  cutg->SetLineColor(kMagenta);
  cutg->SetLineWidth(2);
  cutg->Draw("PL");
  c3->Update();


  TCutG* cutg2 = new TCutG("test",cutg->GetN());

  
  for(int i = 0; i<cutg->GetN();i++){
    cutg2->SetPoint(i,cutg->GetX()[i],cutg->GetY()[i]);
    cutg->GetX()[i] /= 1000;
    cutg->GetY()[i] /= 1000;
  }
  cutg2->SetLineColor(kMagenta);
  cutg2->SetLineWidth(2);

  TCanvas *c = new TCanvas("c","",1400,1000);
  c->Divide(2,2);


  c->cd(1);
  zph->Draw("colz");

  c->cd(2);
  t->Draw("R.tr.vz*1000>>z1d",GeneralCut,"colz");

  c->cd(3);
  zph->Draw("colz");
  cutg2->Draw("PL same");
  
  
  
  c->cd(4);
  t->Draw("R.tr.vz*1000>>z1d_cut",GeneralCut && "Vertex_cut","colz");


}
