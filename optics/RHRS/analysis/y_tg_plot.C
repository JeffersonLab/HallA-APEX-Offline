void y_tg_plot(){

  //Macro makes plots to analyze the new theta and phi after optimization

  TString order = "3rd";    //Optimization order
  //example for range is -10_10 for -10 cm < x_fp <10 cm
  //use "full" for full focal plane range
  TString range = "full";   //Range in focal plane
  TString target = "Optics 3";  //Vertical Wires, Optics 3, or Optics 1
  TString opt = "V_Opt_All";
  
  bool before = false;      //Are we doing before optimization plots
  bool make_plots = true;  

  TString output = "V_wires";
  if(target == "Optics 3") output = "Opt3";
  if(target == "Optics 1") output = "Opt1";
  if(target == "Vert + Opt 1") output = "V_Opt1";
  //output = "Matrix_Choosing";
  
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
    //t->Add(rootfiles + "apex_4647_opt_"+order+"_xfp_full_V_wires.root");
    //t->Add(rootfiles + "apex_4648_opt_"+order+"_xfp_full_V_wires.root");
    //t->Add(rootfiles + "apex_4650_opt_"+order+"_xfp_full_V_wires.root");

    t->Add(rootfiles + "apex_4647_opt_"+order+"_xfp_" + range + "_" + opt + ".root");
    t->Add(rootfiles + "apex_4648_opt_"+order+"_xfp_" + range + "_" + opt + ".root");
    t->Add(rootfiles + "apex_4650_opt_"+order+"_xfp_" + range + "_" + opt + ".root");
  }

  if(before && target == "Optics 3") t->Add(rootfiles + "apex_4652.root");
 
  if(!before && target == "Optics 3") t->Add(rootfiles + "apex_4652_opt_"+order+"_xfp_" + range + "_" + opt + ".root");
  
  if(before && target == "Optics 1") t->Add(rootfiles + "apex_4653.root");
 
  if(!before && target == "Optics 1") t->Add(rootfiles + "apex_4653_opt_"+order+"_xfp_" + range + "_" + opt + ".root");

  if(before && target == "Vert + Opt 1") {
    t->Add(rootfiles + "apex_4653.root");
    //t->Add(rootfiles + "apex_4647.root");
    //t->Add(rootfiles + "apex_4648.root");
    //t->Add(rootfiles + "apex_4650.root");
  }
  if(!before && target == "Vert + Opt 1") {
    t->Add(rootfiles + "apex_4653_opt_"+order+"_xfp_full_V_Opt1.root");
    //t->Add(rootfiles + "apex_4647_opt_"+order+"_xfp_"+range+"_V_Opt1.root");
    //t->Add(rootfiles + "apex_4648_opt_"+order+"_xfp_"+range+"_V_Opt1.root");
    //t->Add(rootfiles + "apex_4650_opt_"+order+"_xfp_"+range+"_V_Opt1.root");
  }

  //Cuts made for all the plots
  TCut GeneralCut;
  //TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && R.s0.nthit==1  && abs(R.tr.tg_dp) < 0.01";
  if(range == "full") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500)";
  if(range == "-10_10") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && abs(R.tr.r_x) < 0.10";
  if(range == "-50_-30") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > -0.50 && R.tr.r_x < -0.30)";
  if(range == "-30_-10") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > -0.30 && R.tr.r_x < -0.10)";
  if(range == "30_50") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > 0.30 && R.tr.r_x < 0.50)";
  if(range == "10_30") GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > 0.10 && R.tr.r_x < 0.30)";

  TH1F * y_tg = new TH1F("Tg y", "", 400, -40, 40);
  
  //Draw tg y plot
  TCanvas* c = new TCanvas("c","c",1000,800);
  t->Draw("R.tr.tg_y*1000>>Tg y",GeneralCut,"");
  if(before) y_tg->SetTitle("Tg y Before Opt;Tg y (mm);");
  else y_tg->SetTitle("Tg y "+order+" Order Opt;Tg y (mm);");
  
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
  pt1->Draw("same");

  double y_exp[3] = {-0.0124,-0.00035,0.01396};
  double z_exp[11] = {-0.205,-0.005,0.195,-0.305,-0.155,0.07,0.214,-0.224,-0.080,0.145,0.295};

  TLine *ly[3];
  TLine *lz[11];
  TLine *lz2[11];

  c->Update();

  int startfoil = 0;
  int endfoil = 3;

  if(target == "Optics 3"){
    startfoil = 7;
    endfoil = 11;
  }
  if(target == "Optics 1"){
    startfoil = 3;
    endfoil = 7;
  }

  /*
  for(int i=startfoil; i<endfoil; i++) {
    ly[i] = new TLine(y_exp[i]*1000,0,y_exp[i]*1000,c->GetUymax());
    ly[i]->SetLineColor(2);
    if(!before && target == "Vertical Wires") ly[i]->Draw("same");
  }
  
  TLegend *leg = new TLegend (0.72, 0.8, 0.9, 0.85);
  leg->AddEntry(ly[0], "Exp Wire Positions", "l");
  if(!before && target == "Vertical Wires") leg->Draw("same");
  */
  
  ///// Set Stats to just show number of events/////
  gStyle->SetOptStat(10);
  TPaveStats *s = (TPaveStats*) gPad->GetPrimitive("stats");
  s->SetX2NDC(0.9);
  s->SetY2NDC(0.9);
  s->SetX1NDC(0.72);
  s->SetY1NDC(0.85);

  
  if(make_plots){
    if(before) c->SaveAs("plots/vertex/"+output+"/y_tg_before_opt_xfp_"+range+".gif");
    else c->SaveAs("plots/vertex/V_Opt_All/" + output + "_y_tg_opt_"+order+"_xfp_"+range+"_" + opt + ".gif");
  }

  double zmax = 400;
  double zmin = -500;

  double ymax = 40;
  double ymin = -40;

  if(before){
    zmin = -150;
    zmax = 150;
  }
  
  TH1F * z_r = new TH1F("Z React", "", 400, zmin, zmax);
  
  //Draw z react plot
  TCanvas* c2 = new TCanvas("c2","c2",1000,800);
  t->Draw("R.tr.vz*1000>>Z React",GeneralCut,"");
  if(before) z_r->SetTitle("Z React Before Opt;Z (mm);");
  else z_r->SetTitle("Z React "+order+" Order Opt;Z (mm);");
    

  pt1->Draw("same");

  

  c2->Update();

  
  TF1 *fgaus[11];
  TF1 *f_refine[11];
  TPaveText *fit_txt[11];


  for(int i=startfoil; i<endfoil; i++) {
    lz[i] = new TLine(z_exp[i]*1000,0,z_exp[i]*1000,c2->GetUymax());
    lz[i]->SetLineColor(2);
    if(!before) lz[i]->Draw("same");
  
    fgaus[i] = new TF1(Form("f%i",i),"gaus",z_exp[i]*1000 - 50,z_exp[i]*1000 + 50);
    z_r->Fit(Form("f%i",i),"qR");


    f_refine[i] = new TF1(Form("fref_%i",i),"gaus",fgaus[i]->GetParameter(1) - 1.2*fgaus[i]->GetParameter(2),fgaus[i]->GetParameter(1) + 1.2*fgaus[i]->GetParameter(2));
    z_r->Fit(Form("fref_%i",i),"qR");
    
    f_refine[i]->Draw("same");

    //fit_txt[i] = new TPaveText(0.12,0.78,0.32,0.89,"nbNDC");
    double max_val = gPad->GetUymax();
    
    fit_txt[i] = new TPaveText(z_exp[i]*1000 + 20,max_val*0.7,z_exp[i]*1000 + 180,max_val*0.7 + 0.15*max_val);
    fit_txt[i]->AddText(Form("#Delta = %g",f_refine[i]->GetParameter(1) - z_exp[i]*1000));
    fit_txt[i]->AddText(Form("#sigma = %g",f_refine[i]->GetParameter(2)));
    fit_txt[i]->SetFillColor(0);
    fit_txt[i]->Draw("same");
  
  }
  

  
  TPaveStats *s2 = (TPaveStats*) gPad->GetPrimitive("stats");
  s2->SetX2NDC(0.9);
  s2->SetY2NDC(0.9);
  s2->SetX1NDC(0.72);
  s2->SetY1NDC(0.85);
  
  TLegend *leg2 = new TLegend (0.72, 0.8, 0.9, 0.85);
  leg2->AddEntry(lz[0], "Exp Wire Positions", "l");
  if(!before && target == "Vertical Wires") leg2->Draw("same");


  if(make_plots){
    if(before) c2->SaveAs("plots/vertex/"+output+"/z_before_opt_xfp_"+range+".gif");
    else c2->SaveAs("plots/vertex/V_Opt_All/" + output + "_z_opt_"+order+"_xfp_"+range+"_" + opt + ".gif");
  }
  

  

  TH2F * zph = new TH2F("zph","", 400, -65, 65,400,zmin,zmax);
  
  //Draw z vs ph plot
  TCanvas* c3 = new TCanvas("c3","",1000,800);
  t->Draw("R.tr.vz*1000:R.tr.tg_ph*1000>>zph",GeneralCut,"colz");
  if(before) zph->SetTitle("Before Y Opt;Tg #phi (mrad);z (mm)");
  else zph->SetTitle("Tg y " + order +" Order Opt;Tg #phi (mrad);z (mm)");


  for(int i=startfoil; i<endfoil; i++) {
    lz2[i] = new TLine(-65,z_exp[i]*1000,65,z_exp[i]*1000);
    lz2[i]->SetLineColor(2);
    if(!before) lz2[i]->Draw("same");
  }
  
  pt1->Draw("same");

  if(make_plots){
    if(before) c3->SaveAs("plots/vertex/"+output+"/z_ph_before_opt_xfp_"+range+".gif");
    else c3->SaveAs("plots/vertex/V_Opt_All/" + output + "_z_ph_opt_"+order+"_xfp_"+range+"_" + opt + ".gif");
  }


  TH2F * yph = new TH2F("yph","", 400, -65, 65,400,ymin,ymax);
  
  //Draw y tg vs ph plot
  TCanvas* c4 = new TCanvas("c4","",1000,800);
  t->Draw("R.tr.tg_y*1000:R.tr.tg_ph*1000>>yph",GeneralCut,"colz");
  if(before) yph->SetTitle("Before Y Opt;Tg #phi (mrad);Tg y (mm)");
  else yph->SetTitle("Tg y " + order +" Order Opt;Tg #phi (mrad);Tg y (mm)");
  
  
  pt1->Draw("same");



  TH2F * yz = new TH2F("yz","", 400, ymin, ymax,400,zmin,zmax);
  
  //Draw y tg vs z plot
  TCanvas* c5 = new TCanvas("c5","",1000,800);
  t->Draw("R.tr.vz*1000:R.tr.tg_y*1000>>yz",GeneralCut,"colz");
  if(before) yz->SetTitle("Before Y Opt;Tg y (mm);z (mm)");
  else yz->SetTitle("Tg y " + order +" Order Opt;Tg y (mm);z (mm)");
  
  
  pt1->Draw("same");

  
  



  
}
