
double Avg(double R[3]){
  double average = 0;
  
  for(int i = 0; i<3; i++) average += R[i]/3;
    
  return average;
}

void projection(TString option = ""){

  
  
  //Macro takes graphical cuts and makes phi and theta projections from them and writes some results to the csv file containing all the information

  TString range = "full_brute";
  // bool opt = true;
  
  // TFile* tcuts;
  // if(opt) tcuts = new TFile("../Sieve/xfp_" + range + "/apex_4647_opt_5th_xfp_"+range+".root.FullCut.root","READ");
  // else tcuts = new TFile("../Sieve/xfp_" + range + "/apex_4647.root.FullCut.root","READ");


  TFile *tcuts = new TFile(CutFileName, "UPDATE");

  TChain* t = Load_more_rootfiles(Run_number, Run_number_2);
  // TChain * t = new TChain("T");

  
  // if(opt)  t->Add("/home/sean/Grad/Research/APEX/Rootfiles/apex_4647_opt_5th_xfp_"+range+".root");
  // else t->Add("/home/sean/Grad/Research/APEX/Rootfiles/apex_4647.root");
  
  //  TCut GenrealCut = "R.tr.n==1 && (R.cer.asum_c>500)";
  double phi_rms = 0;
  double th_rms = 0;
  double phi_cen = 0;
  double th_cen = 0;
  double sy_rms = 0;
  double sx_rms = 0;
  double sy_cen = 0;
  double sx_cen = 0;
  
  int  stat;

  // ifstream cutcsv("../Sieve/xfp_" + range + "/apex_4647.root.cuts_full.csv");

  //  string new_string = "proj" + CutFileName;

  TString new_tstring =  CutFileName;
  new_tstring.Remove(0,30);
  cout << "new_tstring = " << new_tstring << endl;

  TString file_name;
  
  if( option == ""){
    file_name = "proj_csv/"  + new_tstring + ".csv";
  }
  else{
    file_name = "proj_csv/"  + new_tstring + "_" + option + ".csv";
  }
  
  cout << "file_name " << file_name << endl;
  
  ofstream cutcsvnew(file_name);
  string line;

  TDatime* date = new TDatime();  //Get Current date
  
  cutcsvnew<<fixed<<setprecision(2);
  //  getline(cutcsv,line);
  cutcsvnew<<date->GetDay()<<"/"<<date->GetMonth()<<"/"<<date->GetYear()<<" (dd/mm/yyyy)"<<endl;
  // getline(cutcsv,line);
  // cutcsvnew<<line<<endl
  cutcsvnew << "Hole ID: col, row, Included in opt, Ellipse ph cen, Expected ph, Ellipse th cen, Expected th, Semi axis ph, Semi axis th, Ellipse Tilt (deg), -  All angles in mrad except ellipse tilt which is positive counterclockwise from vertical axis"<<endl;
    ;

  gStyle->SetOptStat(11);
  gStyle->SetOptFit(1);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.15);

  TCanvas *c[100] = {NULL};
  TCanvas *c2[100] = {NULL};
  
  int n_canvas = 0;
  // for(int n_col = 0; n_col < 26; n_col++){
  //   for(int n_row = 0; n_row < 17; n_row++){


  int foil_no = 8;



  TRotation fTCSInHCS;
  TVector3 TCSX(0,-1,0);
  TVector3 TCSZ(TMath::Sin(HRSAngle),0,TMath::Cos(HRSAngle));
  TVector3 TCSY = TCSZ.Cross(TCSX);
  fTCSInHCS.RotateAxes(TCSX,TCSY,TCSZ);
  
  
  const TVector3 BeamSpotHCS_average(BeamX_average[foil_no] + (targetfoils[foil_no]/BeamZDir_average)*BeamXDir_average[foil_no], BeamY_average[foil_no] + (targetfoils[foil_no]/BeamZDir_average)*BeamYDir_average[foil_no], targetfoils[foil_no]);


  TVector3 BeamSpotTCS_average = fTCSInHCS.Inverse()*(BeamSpotHCS_average-fPointingOffset);
     

  
  for(int n_hole = 0; n_hole < NHoles; n_hole++){
  // for(int n_hole = 56; n_hole < 58; n_hole++){
  // for(int n_hole = 50; n_hole < 60; n_hole++){

    vector<int> col_row = Get_Col_Row(n_hole);
    int n_col = col_row[0];
    int n_row = col_row[1];

        
  
    TCutG* g = NULL;
    tcuts->GetObject(Form("e_hcut_L_%d_%d_%d",foil_no,n_col,n_row), g);


    if (!g){
      cutcsvnew<<line<<"0,0,0,0"<<endl;
	continue;
    }

    // get expected values of theta/phi

    TVector3 Hole_pos = GetSieveHoleCorrectionTCS(foil_no,n_hole);
    TVector3 MomDirectionTCS_hole = Hole_pos - BeamSpotTCS_average;
      
    Double_t th_exp = MomDirectionTCS_hole.X()/MomDirectionTCS_hole.Z();
    // th_exp = 0;
    // Double_t th_exp_2 = MomDirectionTCS_hole.X()/MomDirectionTCS_hole.Z();

    //    th_exp = Hole_pos.X()// JW adjusted
    
    
    Double_t ph_exp = MomDirectionTCS_hole.Y()/MomDirectionTCS_hole.Z();
    cout << "ph_exp = " << ph_exp*1000 << endl;
    //    ph_exp = Hole_pos.Y();
    
      
    TString name = Form("ID = %d:%d",n_col,n_row);
    TString name2 = Form("ID  = %d:%d",n_col,n_row);

       // for(int i = 0; i<g->GetN();i++){
       // 	  g->GetX()[i] /= 1000;
       // 	  g->GetY()[i] /= 1000;
       // 	}             

    TCut id_cut = TCut(Form("e_hcut_L_%d_%d_%d",foil_no,n_col,n_row));


      c[n_canvas] = new TCanvas(Form("c_%d",n_canvas),"",800,600);
      t->Draw("ph_tgt*1000>>" + name, GenrealCut && id_cut,"");
      //t->Draw("Sieve.y*100>>" + name, GenrealCut && id_cut,"");
      
      gPad->ls();
      
      TH1F *htemp = (TH1F*)gPad->GetPrimitive(name);   
      htemp->SetTitle("");
      htemp->GetXaxis()->SetTitle("Target #phi (mrad)");
      htemp->GetYaxis()->SetTitle("Entries");
      stat = htemp->GetEntries();

      cout << " SSTAT = " << stat << endl;
      
      TF1 *f1 = new TF1("f1","gaus");
      htemp->Fit("f1","q0");

      cout << "First phi mean and width " << f1->GetParameter(1) << " & " << f1->GetParameter(2) << endl;

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

      c2[n_canvas] = new TCanvas(Form("c2_%d",n_canvas),"",800,600);
      //      t->Draw("x_sieve*1000>>" + name2, GenrealCut && id_cut,""); //


      
      // t->Draw(Form("(th_tgt-%f)*1000>>",th_exp) + name2, GenrealCut && id_cut,""); //  altered JW
      t->Draw("th_tgt*1000>>" + name2, GenrealCut && id_cut,""); //  altered JW
      
      

      TH1F *htemp2 = (TH1F*)gPad->GetPrimitive(name2);

      htemp2->SetTitle("");
      htemp2->GetXaxis()->SetTitle("Target #theta (mrad)");
      htemp2->GetYaxis()->SetTitle("Entries");

      TF1 *f2 = new TF1("f2","gaus");
      htemp2->Fit("f2","q0");


      cout << "First theta mean and width " << f2->GetParameter(1) << " & " << f2->GetParameter(2) << endl;
      
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


    
      
      cutcsvnew<<1<<","<<phi_cen<<","<<ph_exp*1000<<","<<th_cen<<","<<th_exp*1000<<","<<phi_rms<<","<<th_rms<<","<<stat<<","<<n_col<<","<<n_row<<endl;
      
      n_canvas++;
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
