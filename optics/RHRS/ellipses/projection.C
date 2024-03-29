TString which_file(int n_foil){

  TString run_number = "4648";

  if(n_foil%2 == 0) run_number = "4653";
  if(n_foil%2 == 1) run_number = "4652";
  
  if(n_foil == 8) run_number = "4648";
  if(n_foil == 9) run_number = "4647";
  if(n_foil == 10) run_number = "4650";



  return run_number;
}


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

  
  TString range = "-10_10";
  bool opt = true;
  
  TFile* tcuts;
  //tcuts = new TFile("../Sieve/Opt_All_Ellipses/xfp_" + range + "/Opt_All_Ellipses.root","READ");
  tcuts = new TFile("../Sieve/V_Opt_All/xfp_full/apex_optics_xfp_full_V_Opt_All.root","READ");

  TString RootDir = "/home/sean/Grad/Research/APEX/Rootfiles/";

  
  //TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && abs(R.tr.r_x) < 0.10";
  TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500)";
  double phi_rms = 0;
  double th_rms = 0;
  double phi_cen = 0;
  double th_cen = 0;
  int  stat;

  
  //ofstream cutcsvnew("Opt_All_Ellipses/xfp_" + range + "/Opt_All_Ellipses.root.EllipseCuts.csv");
  ofstream cutcsvnew("V_Opt_All/xfp_full/apex_optics_test.csv");
  string line;

  TDatime* date = new TDatime();  //Get Current date
  
  cutcsvnew<<fixed<<setprecision(2);
  cutcsvnew<<date->GetDay()<<"/"<<date->GetMonth()<<"/"<<date->GetYear()<<" (dd/mm/yyyy)"<<endl;
  //cutcsvnew<< "Hole ID (foil:col:row),Hole Exists,Included in opt,Ellipse ph cen,Expected ph,Ellipse th cen,Expected th,Ellipse ph rms,Ellipse th rms,Statistics"<<endl;
  cutcsvnew << "Used,phi_cen,ph_exp,th_cen,th_exp,phi_rms_cr,th_rms_cor,stat,n_col,n_row"<<endl;

  gStyle->SetOptStat(11);
  gStyle->SetOptFit(1);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.15);

  TCanvas *c[10][200] = {NULL};
  TCanvas *c2[10][200] = {NULL}; 

  TRotation fTCSInHCS;
  TVector3 TCSX(0,-1,0);
  TVector3 TCSZ(TMath::Sin(HRSAngle),0,TMath::Cos(HRSAngle));
  TVector3 TCSY = TCSZ.Cross(TCSX);
  fTCSInHCS.RotateAxes(TCSX,TCSY,TCSZ);

  for(int n_foil = 8; n_foil < 9; n_foil++){

     const TVector3 BeamSpotHCS_average = BeamSpotHCS_Correction(n_foil,BeamY_average[n_foil], targetfoils[n_foil]);

  TVector3 BeamSpotTCS_average = fTCSInHCS.Inverse()*(BeamSpotHCS_average-fPointingOffset);

  TVector3 Hole_posCen = GetSieveHoleCorrectionTCS(n_foil,14,8);

  TVector3 MomDirectionTCS_holeCen = Hole_posCen - BeamSpotTCS_average;

  Double_t Sieve_ZCen = MomDirectionTCS_holeCen.Z();

  Double_t hole_varEx = TMath::ATan((SieveRadius)/(2*Sieve_ZCen));


  TCutG* g_F = NULL;
  TCutG* g_FP = NULL;
  
  tcuts->GetObject(Form("fcut_R_%d",n_foil), g_F);
  tcuts->GetObject(Form("fcut_R_FP_%d",n_foil), g_FP);
  
 for(int i = 0; i<g_F->GetN();i++){
    g_F->GetX()[i] /= 1000;
    g_F->GetY()[i] /= 1000;
  }

  for(int i = 0; i<g_FP->GetN();i++){
    g_FP->GetX()[i] /= 1000;
    g_FP->GetY()[i] /= 1000;
  }     

  
  TString run_number = which_file(n_foil);
  
  //TFile *Tfile = new TFile(RootDir + "apex_" + run_number + "_opt_3rd_xfp_" + range + "_V_Opt_All.root","read");
  TFile *Tfile = new TFile(RootDir + "apex_" + run_number + "_opt_3rd_xfp_full_V_Opt_All.root","read");
  
  TTree *t = (TTree*)Tfile->Get("T");
  
  int n_canvas = 0;
  
  for(int n_col = 0; n_col < 27; n_col++){
      for(int n_row = 0; n_row < 17; n_row++){
      
	TCutG* g = NULL;
	
	
	tcuts->GetObject(Form("hcut_R_%d_%d_%d",n_foil,n_col,n_row), g);
	


      if (!g){
	//cutcsvnew<<n_foil<<":"<<n_col<<":"<<n_row<<",0,0,0,0,0,0,0,0,0"<<endl;
	cutcsvnew<<line<<"0,0,0,0,0,0,0,0,0,0"<<endl;
	continue;
      }
      
      TString name = Form("ID = %d:%d:%d",n_foil,n_col,n_row);
      TString name2 = Form("ID  = %d:%d:%d",n_foil,n_col,n_row);


      double theta_min = 100;
      double theta_max = -100;
      double phi_min = 100;
      double phi_max = -100;
      
      for(int i = 0; i<g->GetN();i++){
	g->GetX()[i] /= 1000;
	g->GetY()[i] /= 1000;
      }             
      
      TCut foil_cut = TCut(Form("fcut_R_%d",n_foil));
      TCut foil_FP_cut = TCut(Form("fcut_R_FP_%d",n_foil));
      TCut hole_cut = TCut(Form("hcut_R_%d_%d_%d",n_foil,n_col,n_row));
      TCut beam_cut = "0.00120<Rrb.x && Rrb.x<0.0024";
      
      TCut id_cut = foil_cut + foil_FP_cut + hole_cut + beam_cut;

      
      c[n_foil][n_canvas] = new TCanvas(Form("c_%d_%d",n_foil,n_canvas),"",800,600);
      
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
      

      c2[n_foil][n_canvas] = new TCanvas(Form("c2_%d_%d",n_foil,n_canvas),"",800,600);
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

      
      TVector3 Hole_pos = GetSieveHoleCorrectionTCS(n_foil,n_col,n_row);
      TVector3 MomDirectionTCS_hole = Hole_pos - BeamSpotTCS_average;
      
      Double_t th_exp = MomDirectionTCS_hole.X()/MomDirectionTCS_hole.Z();
      Double_t ph_exp = MomDirectionTCS_hole.Y()/MomDirectionTCS_hole.Z();

      //Sieve_hole_pos(n_foil,n_col,n_row,ph_th,yx);
      
      //cutcsvnew<<n_foil<<":"<<n_col<<":"<<n_row<<",1,1,"<<phi_cen<<","<<ph_exp*1000<<","<<th_cen<<","<<th_exp*1000<<","<<phi_rms<<","<<th_rms<<","<<stat<<endl;
      cutcsvnew<<1<<","<<phi_cen<<","<<ph_exp*1000<<","<<th_cen<<","<<th_exp*1000<<","<<phi_rms<<","<<th_rms<<","<<stat<<","<<n_col<<","<<n_row<<endl;    

      n_canvas++;
    }
    }
  }

  TString output = "V_Opt_All/xfp_full/projections_old";
  
  
  for(int n_foil = 0; n_foil < 10; n_foil++){

    TString run_number = which_file(n_foil);
    
    
    int n = 0;
    
    while(c[n_foil][n] != NULL) n++;
    
    for(int i = 0; i<n; i++){
      if(i == 0) {
	c[n_foil][i]->Print(output + "/phi.pdf(","pdf");
	c2[n_foil][i]->Print(output + "/theta.pdf(","pdf");
      }
      else if(i == n-1) {
	c[n_foil][i]->Print(output + "/phi.pdf)","pdf");
	c2[n_foil][i]->Print(output + "/theta.pdf)","pdf");
      }
      else {
	c[n_foil][i]->Print(output + "/phi.pdf","pdf");
	c2[n_foil][i]->Print(output + "/theta.pdf","pdf");
      }
    }

  }
  
}
