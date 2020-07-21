void correction(TH2F *hist){

  TF1 *fit = new TF1("line","pol1");
  hist->Fit("line","q0");

  double b = fit->GetParameter(0);
  double m = fit->GetParameter(1);
  double x_m = hist->GetMean();
  
  for(int i = 0; i <  hist->GetNbinsX();i++){
    int start_bin = 0;
    int end_bin = hist->GetNbinsY();

    if(hist->GetXaxis()->GetBinCenter(i) > x_m){
      start_bin = end_bin;
      end_bin = 0;

      for(int j = start_bin; j > end_bin;j--){
	double y_diff = m*(x_m - hist->GetXaxis()->GetBinCenter(i));
	double y_new = hist->GetYaxis()->GetBinCenter(j) + y_diff;

	if(hist->GetYaxis()->FindBin(y_new) == j) continue;
	
	hist->SetBinContent(i,hist->GetYaxis()->FindBin(y_new),hist->GetBinContent(i,j));
	hist->SetBinContent(i,j,0);
      }
    }
    else if(hist->GetXaxis()->GetBinCenter(i) < x_m){
      for(int j = start_bin; j < end_bin;j++){
	double y_diff = m*(hist->GetXaxis()->GetBinCenter(i) - x_m);
	double y_new = hist->GetYaxis()->GetBinCenter(j) - y_diff;
	
	if(hist->GetYaxis()->FindBin(y_new) == j) continue;
	
	hist->SetBinContent(i,hist->GetYaxis()->FindBin(y_new),hist->GetBinContent(i,j));
	hist->SetBinContent(i,j,0);
      }
      
      
    }
  }


}



double Avg(double R[3]){
  double average = 0;
  
  for(int i = 0; i<3; i++) average += R[i]/3;
    
  return average;
}

void projection(){

  //Macro takes graphical cuts and makes phi and theta projections from them and writes some results to the csv file containing all the information

  TString range = "full_brute";
  bool opt = true;
  
  TFile* tcuts;
  if(opt) tcuts = new TFile("../Sieve/4647/xfp_" + range + "/apex_4647_opt_5th_xfp_"+range+".root.FullCut.root","READ");
  else tcuts = new TFile("../Sieve/4647/xfp_" + range + "/apex_4647.root.FullCut.root","READ");
  TChain * t = new TChain("T");
  if(opt)  t->Add("/home/sean/Grad/Research/APEX/Rootfiles/apex_4647_opt_5th_xfp_"+range+".root");
  else t->Add("/home/sean/Grad/Research/APEX/Rootfiles/apex_4647.root");
  
  TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500)";
  double phi_rms = 0;
  double th_rms = 0;
  double phi_cen = 0;
  double th_cen = 0;
  int  stat;

  ifstream cutcsv("../Sieve/4647/xfp_" + range + "/apex_4647.root.cuts_full.csv");
  ofstream cutcsvnew("xfp_" + range + "/apex_4647.root.EllipseCuts.csv");
  string line;

  TDatime* date = new TDatime();  //Get Current date
  
  cutcsvnew<<fixed<<setprecision(2);
  getline(cutcsv,line);
  cutcsvnew<<date->GetDay()<<"/"<<date->GetMonth()<<"/"<<date->GetYear()<<" (dd/mm/yyyy)"<<endl;
  getline(cutcsv,line);
  cutcsvnew<<line<<endl;

  gStyle->SetOptStat(11);
  gStyle->SetOptFit(1);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.15);

  TCanvas *c[100] = {NULL};
  TCanvas *c2[100] = {NULL};
  
  int n_canvas = 0;
  for(int n_col = 0; n_col < 26; n_col++){
    for(int n_row = 0; n_row < 17; n_row++){
      TCutG* g = NULL;
      tcuts->GetObject(Form("hcut_R_1_%d_%d",n_col,n_row), g);
      getline(cutcsv,line);

      if (!g){
	cutcsvnew<<line<<",0,0,0"<<endl;
	continue;
      }


      stringstream linestream(line);
      string cell;

      string ph_exp,th_exp,ph_axis,th_axis,tilt;
      int n_elem = 1;
      while(getline(linestream,cell,',')){
	if(n_elem<=3) cutcsvnew<<cell<<",";
	if(n_elem == 5) ph_exp = cell;
	if(n_elem == 7) th_exp = cell;
	if(n_elem == 8) ph_axis = cell;
	if(n_elem == 9) th_axis = cell;
	if(n_elem == 10) tilt = cell; 
	n_elem++;
      }
      
      
      TString name = Form("ID = %d:%d",n_col,n_row);
      TString name2 = Form("ID  = %d:%d",n_col,n_row);


      double theta_min = 100;
      double theta_max = -100;
      double phi_min = 100;
      double phi_max = -100;
      
      for(int i = 0; i<g->GetN();i++){
	if(g->GetX()[i] < phi_min) phi_min = g->GetX()[i];
	if(g->GetX()[i] > phi_max) phi_max = g->GetX()[i];
	if(g->GetY()[i] < theta_min) theta_min = g->GetY()[i];
	if(g->GetY()[i] > theta_max) theta_max = g->GetY()[i];
	
	g->GetX()[i] /= 1000;
	g->GetY()[i] /= 1000;
      }             

      TCut id_cut = TCut(Form("hcut_R_1_%d_%d",n_col,n_row));

      
      c[n_canvas] = new TCanvas(Form("c_%d",n_canvas),"",800,600);
      t->Draw("R.tr.tg_ph*1000>>" + name, GeneralCut && id_cut,"");
      //t->Draw("Sieve.y*100>>" + name, GeneralCut && id_cut,"");
      
      
      TH1F *htemp = (TH1F*)gPad->GetPrimitive(name);   
      htemp->SetTitle("");
      htemp->GetXaxis()->SetTitle("Target #phi (mrad)");
      htemp->GetYaxis()->SetTitle("Entries");
      stat = htemp->GetEntries();

      TF1 *f1 = new TF1("f1","gaus");
      htemp->Fit("f1","q0");

      TF1 *fph = new TF1("fph","gaus",f1->GetParameter(1) - 2*f1->GetParameter(2),f1->GetParameter(1) + 2*f1->GetParameter(2));
      
      htemp->Fit("fph","qR");

      double y[3];
      double x1 = fph->GetParameter(1) - 1.6*fph->GetParameter(2);
      double x2 = fph->GetParameter(1) + 1.6*fph->GetParameter(2);

      y[0] = htemp->GetBinContent(htemp->GetXaxis()->FindBin(x1));
      y[1] = htemp->GetBinContent(htemp->GetXaxis()->FindBin(x1) - 1);
      y[2] = htemp->GetBinContent(htemp->GetXaxis()->FindBin(x1) - 2);
      double y1 = Avg(y);
      y[0] = htemp->GetBinContent(htemp->GetXaxis()->FindBin(x2));
      y[1] = htemp->GetBinContent(htemp->GetXaxis()->FindBin(x2) + 1);
      y[2] = htemp->GetBinContent(htemp->GetXaxis()->FindBin(x2) + 2);
      double y2 = Avg(y);
      
      double m = (y2 - y1)/(x2 - x1);
      double b = y1 - m*x1;

      fph = new TF1("fph",Form("[0]*exp(-0.5*((x-[1])/[2])^2) + %g*x + %g",m,b),fph->GetParameter(1) - 2*fph->GetParameter(2),fph->GetParameter(1) + 2*fph->GetParameter(2));
      fph->SetParameters(htemp->GetBinContent(htemp->GetMaximumBin()),htemp->GetMean(),htemp->GetRMS()/2);
      fph->SetParNames("Constant","Position","Sigma");
      
      htemp->Fit("fph","qR");
      fph->SetParameter(2,abs(fph->GetParameter(2)));
      htemp->Fit("fph","qR");
      
      htemp->Draw();
      
      phi_rms = fph->GetParameter(2);
      phi_cen = fph->GetParameter(1);

      
      TH2F * th_y = new TH2F("Th vs y", ";raster y (mm);#theta_{tg} (mrad)", 100, 1.2, 4, 100, theta_min, theta_max);
      t->Draw("R.tr.tg_th*1000:Rrb.y*1000>>Th vs y", GeneralCut && id_cut,"goff");

      correction(th_y);
      

      c2[n_canvas] = new TCanvas(Form("c2_%d",n_canvas),"",800,600);
      TH1D *htemp2 = th_y->ProjectionY(name2);

      
      htemp2->SetTitle("");
      htemp2->GetXaxis()->SetTitle("Target #theta (mrad)");
      htemp2->GetYaxis()->SetTitle("Entries");
      htemp2->GetXaxis()->SetTitleOffset(1.2);

      TF1 *f2 = new TF1("f2","gaus");
      htemp2->Fit("f2","q0");

      TF1 *fth = new TF1("fth","gaus",f2->GetParameter(1) - 2*f2->GetParameter(2),f2->GetParameter(1) + 2*f2->GetParameter(2));
      
      htemp2->Fit("fth","qR");

      x1 = fth->GetParameter(1) - 1.6*fth->GetParameter(2);
      x2 = fth->GetParameter(1) + 1.6*fth->GetParameter(2);

      y[0] = htemp2->GetBinContent(htemp2->GetXaxis()->FindBin(x1));
      y[1] = htemp2->GetBinContent(htemp2->GetXaxis()->FindBin(x1) - 1);
      y[2] = htemp2->GetBinContent(htemp2->GetXaxis()->FindBin(x1) - 2);
      y1 = Avg(y);
      y[0] = htemp2->GetBinContent(htemp2->GetXaxis()->FindBin(x2));
      y[1] = htemp2->GetBinContent(htemp2->GetXaxis()->FindBin(x2) + 1);
      y[2] = htemp2->GetBinContent(htemp2->GetXaxis()->FindBin(x2) + 2);
      y2 = Avg(y);

      m = (y2 - y1)/(x2 - x1);
      b = y1 - m*x1;

      fth = new TF1("fth",Form("[0]*exp(-0.5*((x-[1])/[2])^2) + %g*x + %g",m,b),fth->GetParameter(1) - 2*fth->GetParameter(2),fth->GetParameter(1) + 2*fth->GetParameter(2));
      fth->SetParameters(htemp2->GetBinContent(htemp2->GetMaximumBin()),htemp2->GetMean(),htemp2->GetRMS()/2);
      fth->SetParNames("Constant","Position","Sigma");
      
      htemp2->Fit("fth","qR");
      fth->SetParameter(2,abs(fth->GetParameter(2)));
      htemp2->Fit("fth","qR");
      
      htemp2->Draw();
      th_rms = fth->GetParameter(2);
      th_cen = fth->GetParameter(1);

      
      cutcsvnew<<phi_cen<<","<<ph_exp<<","<<th_cen<<","<<th_exp<<","<<ph_axis<<","<<th_axis<<","<<tilt<<","<<phi_rms<<","<<th_rms<<","<<stat<<endl;

      n_canvas++;
    }
  }

  int n = 0;
  
  while(c[n] != NULL) n++;
  
  for(int i = 0; i<n; i++){
    if(i == 0) {
      c[i]->Print("xfp_" + range + "/projections/phi.pdf(","pdf");
      c2[i]->Print("xfp_" + range + "/projections/theta.pdf(","pdf");
    }
    else if(i == n-1) {
      c[i]->Print("xfp_" + range + "/projections/phi.pdf)","pdf");
      c2[i]->Print("xfp_" + range + "/projections/theta.pdf)","pdf");
    }
    else {
      c[i]->Print("xfp_" + range + "/projections/phi.pdf","pdf");
      c2[i]->Print("xfp_" + range + "/projections/theta.pdf","pdf");
    }
  }
}
