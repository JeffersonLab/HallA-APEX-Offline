void cuts(){

  ////Makes 2D angular plot for each of the foils /////
  
  TString order = "3rd";    //Optimization order
  //example for range is -10_10 for -10 cm < x_fp <10 cm
  //use "full" for full focal plane range
  TString range = "full";   //Range in focal plane
  TString target = "Vertical Wires";  //Vertical Wires, Optics 3, or Optics 1
  TString opt = "V_Opt_All";
  
  bool before = false;      //Are we doing before optimization plots
  bool make_plots = true;  

  TString output = "Vertical Wires";
  if(target == "Optics 3") output = "Opt3";
  if(target == "Optics 1") output = "Opt1";


  TString run = "4647";
  if(target == "Optics 3") run = "4652";
  if(target == "Optics 1") run = "4653";

  TString cut_name = "4647";
  if(target == "Optics 3") cut_name = "V_Opt3_test";
  if(target == "Optics 1") cut_name = "V_Opt1_test";
  //cut_name = "V_Opt_All";

  //Initialize();


  
  
  TFile* tcuts = new TFile("../Sieve/"+opt+"/xfp_"+range+"/apex_"+run+"_opt_"+order+"_xfp_" + range + "_" + opt + ".root.FullCut.root","read");
  //TFile* tcuts = new TFile("../Sieve/"+opt+"/xfp_"+range+"/apex_Vertical_Wires_opt_"+order+"_xfp_" + range + "_" + opt + ".root.FullCut.root","read");

  
  
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
    t->Add(rootfiles + "apex_4647_opt_"+order+"_xfp_full_" + opt + ".root");
    t->Add(rootfiles + "apex_4648_opt_"+order+"_xfp_full_" + opt + ".root");
    t->Add(rootfiles + "apex_4650_opt_"+order+"_xfp_full_" + opt + ".root");
  }
  if(before && target == "Optics 3") t->Add(rootfiles + "apex_4652.root");
 
  if(!before && target == "Optics 3") t->Add(rootfiles + "apex_4652_opt_" + order + "_xfp_" + range + "_" + opt +".root");
  
  if(before && target == "Optics 1") t->Add(rootfiles + "apex_4653.root");
 
  if(!before && target == "Optics 1") t->Add(rootfiles + "apex_4653_opt_" + order + "_xfp_" + range + "_" + opt +".root");

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
  if(range == "-50_-30") pt1->AddText("-0.50 m < x_{fp} < -0.30 m");
  if(range == "-30_-10") pt1->AddText("-0.30 m < x_{fp} < -0.10 m");
  if(range == "-10_10") pt1->AddText("|x_{fp}| < 0.10 m");
  if(range == "10_30") pt1->AddText("0.10 < x_{fp} < 0.30 m");
  if(range == "30_50") pt1->AddText("0.30 < x_{fp} < 0.50 m");
  pt1->SetFillColor(0);

  TText *text = pt1->GetLineWith(target);
  text->SetTextColor(kRed);
  text->SetTextFont(23);
  text->SetTextSize(23);


  TCanvas *c_tg = new TCanvas("c_tg","",1000,1000);
  c_tg->Divide(2,2);
  TH2F *h_tg[4];

  TCanvas *c_cut = new TCanvas("c_cut","",1400,1000);
  c_cut->Divide(3,2);
  TH2F *FP_cut[4];
  
  int startFoil = 3;
  if(target == "Optics 3") startFoil = 7;
  int endFoil = startFoil + 4;

  if(target == "Vertical Wires"){
    startFoil = 0;
    endFoil = startFoil + 3;
  }

  
  for(int n_foil = startFoil; n_foil < endFoil; n_foil++){
    TCutG* g = NULL;
    TCutG* g_FP = NULL;


    tcuts->GetObject(Form("fcut_R_%d",n_foil), g);
    tcuts->GetObject(Form("fcut_R_FP_%d",n_foil), g_FP);
    
    if (!g || !g_FP){
      continue;
    }
    

    
    for(int i = 0; i<g->GetN();i++){
      g->GetX()[i] /= 1000;
      g->GetY()[i] /= 1000;
    }

    for(int i = 0; i<g_FP->GetN();i++){
      g_FP->GetX()[i] /= 1000;
      g_FP->GetY()[i] /= 1000;
    }
    
    TCut foil_cut = Form("fcut_R_%d",n_foil);
    TCut foil_FP_cut = Form("fcut_R_FP_%d",n_foil);

    int foil_num = n_foil;
    if(target == "Optics 1") foil_num = n_foil - 3;
    if(target == "Optics 3") foil_num = foil_num - 7;

    TString name = Form("Foil %d",n_foil);
    c_tg->cd(foil_num + 1);
    h_tg[foil_num] = new TH2F(name,"", 400, -65, 65,400,-65,65);
    
    t->Draw("R.tr.tg_th*1000:R.tr.tg_ph*1000>>"+name,GeneralCut && foil_cut && foil_FP_cut,"colz");

    if(before) h_tg[foil_num]->SetTitle("Tg x' vs y' Before Opt;#phi_{tg} (mrad);#theta_{tg} (mrad)");
    else h_tg[foil_num]->SetTitle("Tg x' vs y' "+order+" Order Opt;#phi_{tg} (mrad);#theta_{tg} (mrad)");

    pt1->Draw("same");

    double sieve_ph[27], sieve_th[17];


    for(int i_col = 0; i_col < 27; i_col++){
      double ph_th[2];
      double yx[2];
      
      Sieve_hole_pos(n_foil,i_col,8,ph_th,yx);
      
      sieve_ph[i_col] = ph_th[0]*1000;
    }
    
    
    for(int i_row = 0; i_row < 17; i_row++){
      double ph_th[2];
      double yx[2];
      
      Sieve_hole_pos(n_foil,4,i_row,ph_th,yx);
    
      sieve_th[i_row] = ph_th[1]*1000;
    }
    
    
    TLine *l[27];
    TLine *l2[17];
    
    int row_min = 0;
    int row_max = 17;
    int col_min = 0;
    int col_max =27;
    
    Hole_exp(n_foil,row_min,row_max,col_min,col_max);
    
    for(int i=col_min;i<=col_max;i++){
      l[i] = new TLine(sieve_ph[i], -65, sieve_ph[i], 65);
      l[i]->SetLineColor(2);
      if(!before) l[i]->Draw("same");
    }
    
    for(int i=row_min;i<=row_max;i++){
      l2[i] = new TLine(-65, sieve_th[i], 65, sieve_th[i]);
      l2[i]->SetLineColor(2);
      if(!before) l2[i]->Draw("same");
    }  

    ///// Draw labels for column and row number ///////
  for(int n_col = col_min; n_col <=col_max; n_col++){
    if(n_col%2 != 0) continue;
    TText *text = new TText;
    text -> SetTextFont(1);
    text -> SetTextSize(0.02);
    text->SetTextAlign(22);
    if(!before) text -> DrawText(sieve_ph[n_col], 67, Form("%d",n_col));
  }

  for(int n_row = row_min; n_row <=row_max; n_row++){
    TText *text = new TText;
    text -> SetTextFont(1);
    text -> SetTextSize(0.02);
    text->SetTextAlign(22);
    if(!before) text -> DrawText(-60,sieve_th[n_row], Form("%d",n_row));
  }


     ///// Set Stats to just show number of events/////
  gStyle->SetOptStat(11);
  TPaveStats *s = (TPaveStats*) gPad->GetPrimitive("stats");
  s->SetX2NDC(0.9);
  s->SetY2NDC(0.9);
  s->SetX1NDC(0.72);
  s->SetY1NDC(0.82);


  }


  if(make_plots){
      if(before) c_tg->SaveAs("plots/vertex_cuts/"+output+"/tg_th_ph_foils_before_opt_xfp_"+range+".gif");
      else c_tg->SaveAs("plots/vertex_cuts/V_Opt_All/"+ output + "_foils_opt_" + order + "_xfp_"+range+"_" + opt + ".gif");
  }
 

  TH2F * zph = new TH2F("zph","", 400, -65, 65,400,-350,350);
  
  //Draw tg y plot
  c_cut->cd(1);
  
  t->Draw("R.tr.vz*1000:R.tr.tg_ph*1000>>zph",GeneralCut,"colz");
  if(before) zph->SetTitle("Before Y Opt;#phi_{tg} (mrad);z (mm)");
  else zph->SetTitle("Tg y "+order+" Order Opt;#phi_{tg} (mrad);z (mm)");
  
  zph->GetYaxis()->SetTitleOffset(1.2);
  
  pt1->Draw("same");
 

  for(int n_foil = startFoil; n_foil < endFoil; n_foil++){
    TCutG* g = NULL;
    TCutG* g_FP = NULL;
    
    tcuts->GetObject(Form("fcut_R_%d",n_foil), g);
    tcuts->GetObject(Form("fcut_R_FP_%d",n_foil), g_FP);
    
    if (!g){
      continue;
    }

    c_cut->cd(1);
    g->Draw("same");

    TCut foil_cut = Form("fcut_R_%d",n_foil);

    int foil_num = n_foil;
    if(target == "Optics 1") foil_num = n_foil - 3;
    if(target == "Optics 3") foil_num = foil_num - 7;


    TString name2 = Form("Foil  %d",n_foil);
    c_cut->cd(foil_num + 2);
    FP_cut[foil_num] = new TH2F(name2, "y_{FP} vs. #phi_{FP} with cut; #phi_{FP} (mrad);y_{FP} (mm)", 400, -40, 50, 400, -50, 50);
    
    t->Draw("R.tr.r_ph*1000:R.tr.r_y*1000>>"+name2,GeneralCut && foil_cut,"colz");
    g_FP->Draw("same");
    
  }

  if(make_plots){
      if(before) c_cut->SaveAs("plots/vertex_cuts/"+output+"/z_ph_before_opt_xfp_"+range+".gif");
      else c_cut->SaveAs("plots/vertex_cuts/V_Opt_All/" + output + "_z_ph_opt_" + order + "_xfp_"+range+"_" + opt + ".gif");
  }

 
}
