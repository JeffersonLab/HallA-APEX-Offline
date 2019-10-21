void hole_id(){

  // Macro that shows all of the holes and the graphical cuts side by side for a nice visualization

  gStyle->SetPalette(1);
  TFile* tcuts = new TFile("../Sieve/xfp_-10_10/apex_4647.root.FullCut.root","read");
  TChain * t = new TChain("T");
  t->Add("/home/sean/Grad/Research/APEX/Rootfiles/apex_4647.root");
  
  
  TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && abs(R.tr.r_x) < 0.10";
  //TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > -0.45 && R.tr.r_x < -0.25)";
  //TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && (R.tr.r_x > 0.25 && R.tr.r_x < 0.45)";
  

  TH2F* xpyp = new TH2F("xpyp","",400,-65,65,400,-65,65);
  //TH2F* xpyp = new TH2F("xpyp","",400,-5,5,400,-6,4);
  t->Draw("R.tr.tg_th*1000:R.tr.tg_ph*1000>>xpyp", GeneralCut,"");
  //t->Draw("Sieve.x*100:Sieve.y*100>>xpyp", GeneralCut,"");

      
  TH1F *htemp = (TH1F*)gPad->GetPrimitive("xpyp");   
  //xpyp->SetTitle("Sieve Plane;y (cm);x (cm)");
  xpyp->SetTitle("Tg th vs Tg ph;#phi (mrad);#theta (mrad)");  
  xpyp->GetYaxis()->SetTitleOffset(1.4);
  xpyp->GetZaxis()->SetRangeUser(0,30);

  TImage* img = TImage::Open("Sieve.png");

  double idx_00 = 0.14, idy_00 = 0.888;
  double dx = 0.042/2, dy = 0.105/2;

  //TCanvas *c = new TCanvas("c","",1000,1000);
  TCanvas *c = new TCanvas("c","",1600,800);
  c->Divide(2,1);
  c->cd(1);
  xpyp->Draw("colz");

  TPaveText *pt1 = new TPaveText(0.12,0.78,0.32,0.89,"nbNDC");
  pt1->AddText("Run 4647");
  pt1->AddText("Cerenkov signal sum > 500");
  pt1->AddText("Single track");
  pt1->AddText("|x_{fp}| < 0.10 m");
  //pt1->AddText("-0.45 < x_{fp} < -0.25 m");
  //pt1->AddText("0.25 < x_{fp} < 0.45 m");
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
  
  c->cd(2);
  img->Draw();

  

  c->Update();

  bool loop = true;
  while (loop){
  
  int n_col;
  int n_row;
  string answer;
  
  cout<<"Select new hole? (y or n): ";
  cin>>answer;
  if(answer == "no" || answer == "n") break; 
  cout<<"Hole Column: ";
  cin>>n_col;
  cout<<"Hole Row: ";
  cin>>n_row;
  
  
  //for(int n_col = 0; n_col < 27; n_col++){
  //for(int n_row = 0; n_row < 17; n_row++){
      TCutG* g = NULL;
  
      tcuts->GetObject(Form("hcut_R_1_%d_%d",n_col,n_row), g);
      
      if (!g){
	cout<<"Hole "<<n_col<<":"<<n_row<<" does not exits"<<endl;
	continue;
      }
      
      cout<<"Showing col:row = "<<n_col<<":"<<n_row<<endl;
      
      c->cd(1);
      g->Draw("same");
      
      TEllipse* Ellipse = new TEllipse(idx_00+n_col*dx,idy_00-n_row*dy,0.02,0.02);
      Ellipse->SetFillStyle(0);
      Ellipse->SetLineColor(kRed);
      Ellipse->SetLineWidth(2);
      c->cd(2);
      Ellipse->Draw("same");
      

      c->Update();
      cin.get();
            
      Ellipse->Delete();
      g->Delete();
      
      
      //}
      //}
    
  }
}
