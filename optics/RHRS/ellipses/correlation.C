void correlation(){


  ////Makes correlation plots between theta and phi target and all the focal plane variables for each hole/////
  
  gStyle->SetPalette(1);

  TString range = "full_brute";
  bool opt = true;
  
  TFile* tcuts;
  if(opt) tcuts = new TFile("../Sieve/xfp_" + range + "/apex_4647_opt_5th_xfp_"+range+".root.FullCut.root","READ");
  else tcuts = new TFile("../Sieve/xfp_" + range + "/apex_4647.root.FullCut.root","READ");
  TChain * t = new TChain("T");
  if(opt)  t->Add("/home/sean/Grad/Research/APEX/Rootfiles/apex_4647_opt_5th_xfp_"+range+".root");
  else t->Add("/home/sean/Grad/Research/APEX/Rootfiles/apex_4647.root");
  
  TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500)";


  gStyle->SetOptStat(11);
  gStyle->SetOptFit(1);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.15);

  
  int n_col = 11;
  int n_row = 9;

  TCanvas *c[4][2];
  TString var[4] = {"x","y","th","ph"};
  double x_range[4] = {0.5,0.03,0.03,0.03};
  
  TCutG* g = NULL;
  tcuts->GetObject(Form("hcut_R_1_%d_%d",n_col,n_row), g);

  double sieve_ph[27], sieve_th[17];

  for(int i = 0; i<27; i++){
    if(i < 25) sieve_ph[i] = (14-i)*2.99 - 12.5;
    else sieve_ph[i] = sieve_ph[24] - (i-24)*2.99*2;
  }
  
  for(int j = 0; j<17; j++){
    sieve_th[j] = (j-8)*7.254 + 3.5;
  }
  
  TString name[4][2] = {{Form("ID = %d:%d",n_col,n_row),Form("ID  = %d:%d",n_col,n_row)},{Form("ID =  %d:%d",n_col,n_row),Form("ID   = %d:%d",n_col,n_row)},{Form("ID =   %d:%d",n_col,n_row),Form("ID  =  %d:%d",n_col,n_row)},{Form("ID  =   %d:%d",n_col,n_row),Form("ID   =  %d:%d",n_col,n_row)}};
  
  
  for(int i = 0; i<g->GetN();i++){
    g->GetX()[i] /= 1000;
    g->GetY()[i] /= 1000;
  }             
  
  TCut id_cut = TCut(Form("hcut_R_1_%d_%d",n_col,n_row));

  for(int i = 0; i < 4; i++){

    TString fp = var[i] + "_{fp}";
    if(i < 2) fp += " (m)";
    else fp += " (rad)";
    
    TH2F *htemp = new TH2F(name[i][0],";"+fp+";Target #phi (mrad)",100,-x_range[i],x_range[i],100,sieve_ph[n_col] - 3,sieve_ph[n_col] + 3);


    c[i][0] = new TCanvas(Form("c%d%d",i,0),"",800,600);
    t->Draw("R.tr.tg_ph*1000:R.tr.r_"+var[i]+">>"+name[i][0], GeneralCut && id_cut,"colz");
    
    htemp->Draw("colz");

    TLine *l = new TLine(-x_range[i],sieve_ph[n_col],x_range[i],sieve_ph[n_col]);
    l->SetLineColor(2);
    l->Draw("same");

    TH2F *htemp2 = new TH2F(name[i][1],";"+fp+";Target #theta (mrad)",100,-x_range[i],x_range[i],100,sieve_th[n_row] - 10,sieve_th[n_row] + 10);

    c[i][1] = new TCanvas(Form("c%d%d",i,1),"",800,600);
    t->Draw("R.tr.tg_th*1000:R.tr.r_"+var[i]+">>" + name[i][1], GeneralCut && id_cut,"colz");
    
    htemp2->Draw("colz");

    TLine *l2 = new TLine(-x_range[i],sieve_th[n_row],x_range[i],sieve_th[n_row]);
    l2->SetLineColor(2);
    l2->Draw("same");
  }

  for(int i = 0; i<4; i++){
    if(i==0) c[i][0]->Print("xfp_" + range + "/correlations/"+Form("correlation_%d_%d",n_col,n_row)+".pdf(","pdf");
    else c[i][0]->Print("xfp_" + range + "/correlations/"+Form("correlation_%d_%d",n_col,n_row)+".pdf","pdf");
  }

  for(int i = 0; i<4; i++){
    if(i==3) c[i][1]->Print("xfp_" + range + "/correlations/"+Form("correlation_%d_%d",n_col,n_row)+".pdf)","pdf");
    else c[i][1]->Print("xfp_" + range + "/correlations/"+Form("correlation_%d_%d",n_col,n_row)+".pdf","pdf");
  }


}
