void apex_optics_right(){

  ///Makes a bunch of plots of all the correlated variables///
  
  gStyle->SetPalette(1);

  TString opt = "no";
  
  Int_t run_no;
  cout <<"\n: Please enter a Run Number (-1 to exit):";
  cin >> run_no;

  cout <<"\n: Use optimization? (yes or no):";
  cin >> opt;

  // ----------------------------------------------------------------------

  Int_t    prl2_min  = 800;                   // the limit for (L.prl1.e+prl2_e*L.prl2.e)
  Double_t prl2_e    = 0.8;                   // for (L.prl1.e+prl2_e*L.prl2.e)
  Int_t    prl1_min  = 200;                   // for L.prl1.e
  Int_t    cer_min   = 2000;                  // the limit for Cherenkov
  Double_t phmin     = -0.06, phmax = 0.06;  // the limits for tg_ph
  Double_t thmin     = -0.08, thmax = 0.08;  // the limits for tg_th
  Double_t dpmin     = -0.06, dpmax = 0.06;  // the limits for tg_dp
  Double_t ytmin     = -0.06, ytmax = 0.06;  // the limits for tg_y

  // ----------------------------------------------------------------------

  TCut cut_nclust = "R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1";
  TCut cut_ntrack = "R.tr.n==1&&R.s0.nthit==1";
  TCut GeneralCut = cut_nclust + cut_ntrack; 
  TCut ct1        = Form("(R.prl1.asum_c+%.2f*R.prl2.asum_c)>%d&&R.prl1.asum_c>%d", prl2_e, prl2_min, prl1_min);
  TCut ct2        = Form("R.cer.asum_c>%d", cer_min);
  TCut ct3        = Form("R.tr.tg_ph<%.2f&&R.tr.tg_ph>%.2f", phmax, phmin);
  TCut ct4        = Form("R.tr.tg_th<%.2f&&R.tr.tg_th>%.2f", thmax, thmin);
  TCut ct5        = Form("R.tr.tg_dp<%.2f&&R.tr.tg_dp>%.2f", dpmax, dpmin);
  TCut VertexCut  = Form("R.tr.tg_y<%.2f&&R.tr.tg_y>%.2f", ytmax, ytmin);

  TCut fp1        = Form("R.tr.y<%.2f&&R.tr.y>%.2f", -0.01, -0.1);

  // ----------------------------------------------------------------------

  TString Rootfiles = "/home/sean/Grad/Research/APEX/Rootfiles/";
  
  TChain *T = new TChain("T");
  if( run_no != 9999 ){
    if(opt == "no")
      T->Add(Rootfiles + Form("apex_%d.root", run_no));
    else
      T->Add(Rootfiles + Form("apex_%d_opt_5th_xfp_full_V_wires.root", run_no));
  }
  else {
  T->Add(Rootfiles + Form("apex_4766.root"));
  T->Add(Rootfiles + Form("apex_4768.root"));
  T->Add(Rootfiles + Form("apex_4769.root"));
  }

  // ----------------------------------------------------------------------
  // 1D RHRS Focal plane variables 
  // ----------------------------------------------------------------------

  TCanvas *c1 = new TCanvas("c1","RHRS Focal Plane Variables",1000,800); 
  c1->Divide(3,2); 

  TH1F* h71 = new TH1F("h71","h71",300,-1.5,1.5); 
  TH1F* h72 = new TH1F("h72","h72",300,-0.2,0.2); 
  TH1F* h73 = new TH1F("h73","h73",300,-0.1,0.1); 
  TH1F* h74 = new TH1F("h74","h74",300,-0.2,0.2); 
  Int_t hmax, hmin;

  c1->cd(1);
  /*
  TH2F *h1 = new TH2F("h1","h1",100,0,2000,100,0,2000);
  T->Draw("R.prl1.asum_c:R.prl2.asum_c>>h1",GeneralCut,"colz");
  assert(h1);
  h1->SetTitle("RHRS Shower sum vs Preshower sum");
  h1->GetYaxis()->SetTitleOffset(1.0);
  h1->GetXaxis()->SetTitleSize(0.05);
  h1->GetYaxis()->SetTitleSize(0.05);
  h1->GetXaxis()->SetTitle("R.prl2.asum_c");
  h1->GetYaxis()->SetTitle("R.prl1.asum_c");
  TLine *line1 = new TLine(0,prl2_min,prl2_min/prl2_e,0);
  line1->SetLineColor(2);
  line1->SetLineWidth(2);
  line1->Draw();
  TLine *line2 = new TLine(0,prl1_min,2000,prl1_min);
  line2->SetLineColor(2);
  line2->SetLineWidth(2);
  line2->Draw();
  */

  c1->cd(2)->SetLogy(1);
    TH1F *h2 = new TH1F("h2","h2",500,0,5000);
    //T->Draw("R.cer.asum_c>>h2",GeneralCut&&ct1);
    T->Draw("R.cer.asum_c>>h2",GeneralCut);
  h2->SetTitle("RHRS Cherenkov sum");
  hmin = h2->GetMinimum();
  hmax = h2->GetMaximum();
  TLine *line3 = new TLine(cer_min,hmin,cer_min,hmax);
  line3->SetLineColor(2);
  line3->SetLineWidth(2);
  line3->Draw();

  //GeneralCut += ct1+ct2;
  GeneralCut += ct2;
  GeneralCut += ct3+ct4+ct5;
  GeneralCut += VertexCut;
  //  GeneralCut += fp1;

  c1->cd(3); 
  T->Draw("R.tr.y>>h72",GeneralCut,"colz"); 
  h72->SetTitle("RHRS FP y (in-plane)"); 
  hmin = h72->GetMinimum();
  hmax = h72->GetMaximum();
  TLine *line18 = new TLine(-0.15,hmin,-0.15,hmax);
  line18->SetLineColor(2);
  line18->SetLineWidth(2);
  line18->Draw();
  TLine *line19 = new TLine(0.15,hmin,0.15,hmax);
  line19->SetLineColor(2);
  line19->SetLineWidth(2);
  line19->Draw();

  c1->cd(4); 
  T->Draw("R.tr.x>>h71",GeneralCut,"colz"); 
  h71->SetTitle("RHRS FP x (out-of-plane)"); 
  hmin = h71->GetMinimum();
  hmax = h71->GetMaximum();
  TLine *line118 = new TLine(-0.7,hmin,-0.7,hmax);
  line118->SetLineColor(2);
  line118->SetLineWidth(2);
  line118->Draw();
  TLine *line119 = new TLine(0.7,hmin,0.7,hmax);
  line119->SetLineColor(2);
  line119->SetLineWidth(2);
  line119->Draw();

  c1->cd(5); 
  T->Draw("R.tr.ph>>h73",GeneralCut,"colz"); 
  h73->SetTitle("RHRS FP ph (in-plane)"); 
  hmin = h73->GetMinimum();
  hmax = h73->GetMaximum();
  TLine *line218 = new TLine(-0.06,hmin,-0.06,hmax);
  line218->SetLineColor(2);
  line218->SetLineWidth(2);
  line218->Draw();
  TLine *line219 = new TLine(0.06,hmin,0.06,hmax);
  line219->SetLineColor(2);
  line219->SetLineWidth(2);
  line219->Draw();

  c1->cd(6); 
  T->Draw("(R.tr.th)>>h74",GeneralCut,"colz"); 
  h74->SetTitle("RHRS FP th (out-of-plane)"); 
  hmin = h74->GetMinimum();
  hmax = h74->GetMaximum();
  TLine *line318 = new TLine(-0.12,hmin,-0.12,hmax);
  line318->SetLineColor(2);
  line318->SetLineWidth(2);
  line318->Draw();
  TLine *line319 = new TLine(0.12,hmin,0.12,hmax);
  line319->SetLineColor(2);
  line319->SetLineWidth(2);
  line319->Draw();

  // ----------------------------------------------------------------------
  // 2D RHRS Focal plane variables 
  // ----------------------------------------------------------------------

  TCanvas *c2 = new TCanvas("c2","RHRS Focal Plane Variables 2D",1000,800); 
  c2->Divide(2,2); 

  TH2F* h701 = new TH2F("h701","h701",100,-0.1,0.1,100,-1.2,1.2); 
  TH2F* h702 = new TH2F("h702","h702",100,-0.1,0.1,100,-0.15,0.15); 
  TH2F* h703 = new TH2F("h703","h703",100,-0.15,0.15,100,-1.2,1.2); 
  TH2F* h704 = new TH2F("h704","h704",100,-0.1,0.1,100,-0.1,0.1);

  c2->cd(1); 
  T->Draw("R.tr.x:R.tr.y>>h701",GeneralCut,"colz"); 
  h701->SetTitle("RHRS FP x vs y"); 

  c2->cd(2); 
  T->Draw("(R.tr.th):(R.tr.ph)>>h702",GeneralCut,"colz"); 
  h702->SetTitle("RHRS FP th vs ph"); 

  c2->cd(3); 
  T->Draw("(R.tr.x):(R.tr.th)>>h703",GeneralCut,"colz"); 
  h703->SetTitle("RHRS FP x vs th"); 

  c2->cd(4); 
  T->Draw("R.tr.y:R.tr.ph>>h704",GeneralCut,"colz"); 
  h704->SetTitle("RHRS FP y vs ph"); 

  // ----------------------------------------------------------------------
  // RHRS Target variables 
  // ----------------------------------------------------------------------

  TCanvas *c3 = new TCanvas("c3","RHRS Target Variables",1000,800);
  c3->Divide(2,2);

  c3->cd(1);
  TH1F *h3 = new TH1F("h3","h3",300, -0.1, 0.1);
  T->Draw("R.tr.tg_ph>>h3",GeneralCut);
  h3->SetTitle("RHRS target ph (in-plane)");
  hmax = h3->GetMaximum()*1.05;
  hmin = h3->GetMinimum();
  TLine *line4 = new TLine(phmin,hmin,phmin,hmax);
  line4->SetLineColor(2);
  line4->SetLineWidth(2);
  line4->Draw();
  TLine *line5 = new TLine(phmax,hmin,phmax,hmax);
  line5->SetLineColor(2);
  line5->SetLineWidth(2);
  line5->Draw();

  c3->cd(2);
  TH1F *h4 = new TH1F("h4","h4",300, -0.2, 0.2);
  T->Draw("R.tr.tg_th>>h4",GeneralCut);
  h4->SetTitle("RHRS target th (out-of-plane)");
  hmax = h4->GetMaximum()*1.05;
  hmin = h4->GetMinimum();
  TLine *line6 = new TLine(thmin,hmin,thmin,hmax);
  line6->SetLineColor(2);
  line6->SetLineWidth(2);
  line6->Draw();
  TLine *line7 = new TLine(thmax,hmin,thmax,hmax);
  line7->SetLineColor(2);
  line7->SetLineWidth(2);
  line7->Draw();

  c3->cd(3);
  TH1F *h5 = new TH1F("h5","h5",300, -0.1, 0.1);
  T->Draw("R.tr.tg_dp>>h5",GeneralCut);
  h5->SetTitle("RHRS target delta");
  hmax = h5->GetMaximum()*1.05;
  hmin = h5->GetMinimum();
  TLine *line8 = new TLine(dpmin,hmin,dpmin,hmax);
  line8->SetLineColor(2);
  line8->SetLineWidth(2);
  line8->Draw();
  TLine *line9 = new TLine(dpmax,hmin,dpmax,hmax);
  line9->SetLineColor(2);
  line9->SetLineWidth(2);
  line9->Draw();

  c3->cd(4);
  TH1F *h63 = new TH1F("h63","h63",300, -0.1, 0.1);
  T->Draw("R.tr.tg_y>>h63",GeneralCut);
  h63->SetTitle("RHRS target y");
  h63->GetXaxis()->SetTitleSize(0.05);
  h63->SetLineColor(4);
  hmax = h63->GetMaximum()*1.05;
  hmin = h63->GetMinimum();
  TLine *line808 = new TLine(ytmin,hmin,ytmin,hmax);
  line808->SetLineColor(2);
  line808->SetLineWidth(2);
  line808->Draw();
  TLine *line909 = new TLine(ytmax,hmin,ytmax,hmax);
  line909->SetLineColor(2);
  line909->SetLineWidth(2);
  line909->Draw();

  // ----------------------------------------------------------------------
  // RHRS Target variables 
  // ----------------------------------------------------------------------

  //GeneralCut += ct3+ct4+ct5;
  TCanvas *c4 = new TCanvas("c4","RHRS Target Variables 2D",1000,800);
  c4->Divide(3,2);

  TH2F *h33 = new TH2F("h33","h33", 100,-0.05,0.05,100,-0.05,0.05);
  TH2F *h21 = new TH2F("h21","h21", 100,-0.05,0.05,100,-0.05,0.05);
  TH2F *h22 = new TH2F("h22","h22", 100,-0.05,0.05,100,-0.05,0.05);
  TH2F *h23 = new TH2F("h23","h23", 100,-0.05,0.05,100,-0.05,0.05);
  TH2F *h34 = new TH2F("h34","h34", 100,-0.05,0.05,100,-0.05,0.05);
  TH2F *h35 = new TH2F("h35","h35", 100,-0.05,0.05,100,-0.05,0.05);

  c4->cd(1);
  T->Draw("R.tr.tg_th:R.tr.tg_ph>>h33",GeneralCut,"colz");
  h33->SetTitle("RHRS th_tgt vs ph_tgt");
   
  c4->cd(2);
  T->Draw("R.tr.tg_y:R.tr.tg_ph>>h21",GeneralCut,"colz");
  //T->Draw("R.tr.vz:R.tr.tg_ph>>h21",GeneralCut,"colz");
  h21->SetTitle("RHRS y_tgt vs ph_tgt");

  c4->cd(3);
  T->Draw("R.tr.tg_y:R.tr.tg_th>>h22",GeneralCut,"colz");
  //T->Draw("R.tr.vz:R.tr.tg_ph>>h21",GeneralCut,"colz");
  h22->SetTitle("RHRS y_tgt vs th_tgt");

  c4->cd(4);
  T->Draw("R.tr.tg_y:R.tr.tg_dp>>h23",GeneralCut,"colz");
  //T->Draw("R.tr.vz:R.tr.tg_ph>>h21",GeneralCut,"colz");
  h23->SetTitle("RHRS y_tgt vs delta");

  c4->cd(5);
  T->Draw("R.tr.tg_dp:R.tr.tg_th>>h34",GeneralCut,"colz");
  //T->Draw("R.tr.vz:R.tr.tg_ph>>h21",GeneralCut,"colz");
  h34->SetTitle("RHRS delta vs th_tgt");

  c4->cd(6);
  T->Draw("R.tr.tg_dp:R.tr.tg_ph>>h35",GeneralCut,"colz");
  //T->Draw("R.tr.vz:R.tr.tg_ph>>h21",GeneralCut,"colz");
  h35->SetTitle("RHRS delta vs ph_tgt");

  c1->Print("apex_optics_left1.pdf");
  c2->Print("apex_optics_left2.pdf");
  c3->Print("apex_optics_left3.pdf");
  c4->Print("apex_optics_left4.pdf");
  
  if(opt == "no") gSystem->Exec(Form("pdfunite apex_optics_left*.pdf RHRS%d_apex_optics.pdf", run_no));
  else gSystem->Exec(Form("pdfunite apex_optics_left*.pdf RHRS%d_apex_optics_optimized.pdf", run_no)); 
  gSystem->Exec("rm apex_optics_left*.pdf"); 


}



