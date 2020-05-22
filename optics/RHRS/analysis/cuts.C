void cuts(){

  TString order = "5th";    //Optimization order
  //example for range is -10_10 for -10 cm < x_fp <10 cm
  //use "full" for full focal plane range
  TString range = "30_50";   //Range in focal plane
  TString target = "Optics 1";  //Vertical Wires Optics 3 or Optics 1

  bool before = false;      //Are we doing before optimization plots
  bool make_plots = true;  

  TString output = "V_wires";
  if(target == "Optics 3") output = "Opt3";
  if(target == "Optics 1") output = "Opt1";

  TString run = "4647";
  if(target == "Optics 3") run = "4652";
  if(target == "Optics 1") run = "4653";
  
  TFile* tcuts = new TFile("../Sieve/"+run+"/xfp_"+range+"/apex_"+run+"_opt_5th_xfp_full_V_wires.root.FullCut.root","read");
  
  TString rootfiles = "/home/sean/Grad/Research/APEX/Rootfiles/";
  gStyle->SetPalette(1);
  TFile* f;


  TChain* t = new TChain("T");
  
  if(before && target == "Vertical Wires") {
    t->Add(rootfiles + "apex_4647.root");
    t->Add(rootfiles + "apex_4648.root");
    t->Add(rootfiles + "apex_4650.root");
  }
  if(!before && target == "Vertical Wires") {
    t->Add(rootfiles + "apex_4647_opt_5th_xfp_full_V_wires.root");
    t->Add(rootfiles + "apex_4648_opt_5th_xfp_full_V_wires.root");
    t->Add(rootfiles + "apex_4650_opt_5th_xfp_full_V_wires.root");
  }
  if(before && target == "Optics 3") t->Add(rootfiles + "apex_4652.root");
 
  if(!before && target == "Optics 3") t->Add(rootfiles + "apex_4652_opt_5th_xfp_full_V_wires.root");
  
  if(before && target == "Optics 1") t->Add(rootfiles + "apex_4653.root");
 
  if(!before && target == "Optics 1") t->Add(rootfiles + "apex_4653_opt_5th_xfp_full_V_wires.root");

  

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


  /////Add labels for the run number and cuts ////
  TPaveText *pt1 = new TPaveText(0.12,0.78,0.32,0.89,"nbNDC");
  pt1->AddText(target);
  pt1->AddText("Cerenkov signal sum > 500");
  pt1->AddText("Single track");
  if(range == "-10_10") pt1->AddText("|x_{fp}| < 0.10 m");
  if(range == "-50_-30") pt1->AddText("-0.50 m < x_{fp} < -0.30 m");
  if(range == "30_50") pt1->AddText("0.30 < x_{fp} < 0.50 m");
  pt1->SetFillColor(0);

  TText *text = pt1->GetLineWith(target);
  text->SetTextColor(kRed);
  text->SetTextFont(23);
  text->SetTextSize(23);


  TCanvas *c_tg = new TCanvas("c_tg","",1000,1000);
  c_tg->Divide(2,2);
  TH2F *h_tg[4];
  
  for(int n_foil = 0; n_foil < 4; n_foil++){
    TCutG* g = NULL;
    
    tcuts->GetObject(Form("fcut_R_%d",n_foil), g);
    
    if (!g){
      continue;
    }
    

    
    for(int i = 0; i<g->GetN();i++){
      g->GetX()[i] /= 1000;
      g->GetY()[i] /= 1000;
    }
    
    TCut foil_cut = Form("fcut_R_%d",n_foil);
    
    TString name = Form("Foil_%d",n_foil);
    c_tg->cd(n_foil + 1);
    h_tg[n_foil] = new TH2F(name,"", 400, -65, 65,400,-65,65);
    
    t->Draw("R.tr.tg_th*1000:R.tr.tg_ph*1000>>"+name,GeneralCut && foil_cut,"colz");

    if(before) h_tg[n_foil]->SetTitle("Tg x' vs y' Before Opt;Tg #phi (mrad);Tg #theta (mrad)");
    else h_tg[n_foil]->SetTitle("Tg x' vs y' "+order+" Order Opt;Tg #phi (mrad);Tg #theta (mrad)");

    pt1->Draw("same");
    
      
  }

  if(make_plots){
      if(before) c_tg->SaveAs("plots/vertex_cuts/"+output+"/tg_th_ph_foils_before_opt_xfp_"+range+".gif");
      else c_tg->SaveAs("plots/vertex_cuts/"+output+"/tg_th_ph_foils_opt_xfp_"+range+".gif");
  }
 

  TH2F * zph = new TH2F("zph","", 400, -65, 65,400,-600,350);
  
  //Draw tg y plot
  TCanvas* c = new TCanvas("c","c",1000,800);
  t->Draw("R.tr.vz*1000:R.tr.tg_ph*1000>>zph",GeneralCut,"colz");
  if(before) zph->SetTitle("Before Y Opt;Tg #phi (mrad);z (mm)");
  else zph->SetTitle("After Tg y Opt;Tg #phi (mrad);z (mm)");
  
  pt1->Draw("same");
  

  for(int n_foil = 0; n_foil < 4; n_foil++){
    TCutG* g = NULL;
    
    tcuts->GetObject(Form("fcut_R_%d",n_foil), g);
    
    if (!g){
      continue;
    }

    g->Draw("same");
    
  }

  if(make_plots){
      if(before) c->SaveAs("plots/vertex_cuts/"+output+"/z_ph_before_opt_xfp_"+range+".gif");
      else c->SaveAs("plots/vertex_cuts/"+output+"/z_ph_opt_xfp_"+range+".gif");;
  }

 
}
