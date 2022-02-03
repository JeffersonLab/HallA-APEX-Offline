TString which_file(int n_foil){

  TString run_number = "4648";

  if(n_foil%2 == 0) run_number = "4653";
  if(n_foil%2 == 1) run_number = "4652";
  
  if(n_foil == 8) run_number = "4648";
  if(n_foil == 9) run_number = "4647";
  if(n_foil == 10) run_number = "4650";



  return run_number;
}


double Avg(double R[3]){
  double average = 0;
  
  for(int i = 0; i<3; i++) average += R[i]/3;
    
  return average;
}


Double_t CircFunc(Double_t *x, Double_t *par)
{
   Double_t xx =x[0];
   Double_t f = 0.0;
   if(par[0]<abs(xx-par[1])){
     f = 0.0;
   }
   else{
     f = par[2]*2*(sqrt(par[0]*par[0]-(xx-par[1])*(xx-par[1])));
   }
   
   return f;
   
}


Double_t GausFunc(Double_t *x, Double_t *par)
{
   Double_t xx =x[0];
   //   Double_t f = exp(-0.5*((xx)/par[0])^2);
   Double_t f = (1/(sqrt(2*TMath::Pi())*par[0]))*exp(-0.5*pow((xx)/par[0],2));
      
   return f;   
}


void projection_LHRS2(Int_t foil_no, TString option = ""){

  
  
  //Macro takes graphical cuts and makes phi and theta projections from them and writes some results to the csv file containing all the information

  TString range = "full_brute";
  TString output = "V_Opt_All/xfp_full";
  TString Rootdir = "/home/sean/Grad/Research/APEX/Rootfiles/";
  TString run = which_file(foil_no);
  // bool opt = true;
  
  // TFile* tcuts;
  // if(opt) tcuts = new TFile("../Sieve/xfp_" + range + "/apex_4647_opt_5th_xfp_"+range+".root.FullCut.root","READ");
  // else tcuts = new TFile("../Sieve/xfp_" + range + "/apex_4647.root.FullCut.root","READ");
  TString CutFileName = "../Sieve/V_Opt_All/xfp_full/apex_optics_xfp_full_V_Opt_All.root";

  
  TFile *tcuts = new TFile(CutFileName, "read");

  
  //TChain* t = Load_more_rootfiles(Run_number, Run_number_2);
  TChain * t = new TChain("T");
  t->Add(Rootdir + "apex_" + run + "_opt_3rd_xfp_full_V_Opt_All.root");
  
  // if(opt)  t->Add("/home/sean/Grad/Research/APEX/Rootfiles/apex_4647_opt_5th_xfp_"+range+".root");
  // else t->Add("/home/sean/Grad/Research/APEX/Rootfiles/apex_4647.root");
  
  TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500)";
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
  new_tstring.Remove(0,28);
  cout << "new_tstring = " << new_tstring << endl;

  TString file_name;
  
  if( option == ""){
    file_name = output + "/"  + new_tstring + Form("_%i",foil_no) + ".csv";
  }
  else{
    file_name = output + "/"  + new_tstring + "_" + Form("_%i",foil_no) + option + ".csv";
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
  //cutcsvnew << "Hole ID: col, row, Included in opt, Ellipse ph cen, Expected ph, Ellipse th cen, Expected th, Semi axis ph, Semi axis th, Ellipse Tilt (deg), -  All angles in mrad except ellipse tilt which is positive counterclockwise from vertical axis"<<endl;
  cutcsvnew << "Used,phi_cen,ph_exp,th_cen,th_exp,phi_rms_cr,th_rms_cor,stat,n_col,n_row"<<endl;


  // gStyle->SetOptStat(11);
    //    gStyle->SetOptFit(1);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.15);

  TCanvas *c[500] = {NULL};
  TCanvas *c2[500] = {NULL};
  
  int n_canvas = 0;
  // for(int n_col = 0; n_col < 26; n_col++){
  //   for(int n_row = 0; n_row < 17; n_row++){


  //  int foil_no = 4;



  TRotation fTCSInHCS;
  TVector3 TCSX(0,-1,0);
  TVector3 TCSZ(TMath::Sin(HRSAngle),0,TMath::Cos(HRSAngle));
  TVector3 TCSY = TCSZ.Cross(TCSX);
  fTCSInHCS.RotateAxes(TCSX,TCSY,TCSZ);
  
  
  //const TVector3 BeamSpotHCS_average(BeamX_average[foil_no] + (targetfoils[foil_no]/BeamZDir_average)*BeamXDir_average[foil_no], BeamY_average[foil_no] + (targetfoils[foil_no]/BeamZDir_average)*BeamYDir_average[foil_no], targetfoils[foil_no]);
  const TVector3 BeamSpotHCS_average = BeamSpotHCS_Correction(foil_no,BeamY_average[foil_no], targetfoils[foil_no]);

  TVector3 BeamSpotTCS_average = fTCSInHCS.Inverse()*(BeamSpotHCS_average-fPointingOffset);

  // calculate distance from sieve to foil for central hole
  
  TVector3 Hole_posCen = GetSieveHoleCorrectionTCS(foil_no,14,8);

  TVector3 MomDirectionTCS_holeCen = Hole_posCen - BeamSpotTCS_average;

  Double_t Sieve_ZCen = MomDirectionTCS_holeCen.Z();

  cout << "Distance from central hole to pivot = " << Sieve_ZCen << endl;

  Double_t hole_varEx = TMath::ATan((SieveRadius)/(2*Sieve_ZCen));

  cout << "Resulting sigma^2 for survery = " << hole_varEx*hole_varEx << endl;

  
  // add beam cut
  //  GeneralCut += get_Beamcut(Run_number);

  TCutG* g_F = NULL;
  TCutG* g_FP = NULL;

  tcuts->GetObject(Form("fcut_R_%d",foil_no), g_F);
  tcuts->GetObject(Form("fcut_R_FP_%d",foil_no), g_FP);

  for(int i = 0; i<g_F->GetN();i++){
    g_F->GetX()[i] /= 1000;
    g_F->GetY()[i] /= 1000;
  }

  for(int i = 0; i<g_FP->GetN();i++){
    g_FP->GetX()[i] /= 1000;
    g_FP->GetY()[i] /= 1000;
  }          
  
  //  for(int n_hole = 0; n_hole < NHoles; n_hole++){
  for(int n_col = 0; n_col < 27; n_col++){
    for(int n_row = 0; n_row < 17; n_row++){
    
	
        
  
    TCutG* g = NULL;
    tcuts->GetObject(Form("hcut_R_%d_%d_%d",foil_no,n_col,n_row), g);


    if (!g){
      //cutcsvnew<<line<<"0,0,0,0"<<endl;
      cutcsvnew<<line<<"0,0,0,0,0,0,0,0,0,0"<<endl;
	continue;
    }

    // get expected values of theta/phi

    TVector3 Hole_pos = GetSieveHoleCorrectionTCS(foil_no,n_col,n_row);
    TVector3 MomDirectionTCS_hole = Hole_pos - BeamSpotTCS_average;
      
    Double_t th_exp = MomDirectionTCS_hole.X()/MomDirectionTCS_hole.Z();
    // th_exp = 0;
    // Double_t th_exp_2 = MomDirectionTCS_hole.X()/MomDirectionTCS_hole.Z();

    //    th_exp = Hole_pos.X()// JW adjusted
    
    
    Double_t ph_exp = MomDirectionTCS_hole.Y()/MomDirectionTCS_hole.Z();
    cout << "ph_exp(" << n_col << "," << n_row << ") = " << ph_exp*1000 << endl;
    //    ph_exp = Hole_pos.Y();
    
      
    TString name = Form("ID = %d:%d",n_col,n_row);
    TString name2 = Form("ID  = %d:%d",n_col,n_row);

    for(int i = 0; i<g->GetN();i++){
      g->GetX()[i] /= 1000;
      g->GetY()[i] /= 1000;
    }             

    TCut foil_cut = TCut(Form("fcut_R_%d",foil_no));
    TCut foil_FP_cut = TCut(Form("fcut_R_FP_%d",foil_no));
    TCut hole_cut = TCut(Form("hcut_R_%d_%d_%d",foil_no,n_col,n_row));
    TCut beam_cut = "0.00120<Rrb.x && Rrb.x<0.0024";

    TCut id_cut = foil_cut + foil_FP_cut + hole_cut + beam_cut;
    

      c[n_canvas] = new TCanvas(Form("c_%d",n_canvas),"",800,600);
      t->Draw("R.tr.tg_ph*1000>>" + name, GeneralCut && id_cut,"");
      //t->Draw("Sieve.y*100>>" + name, GeneralCut && id_cut,"");
      
      //      gPad->ls();
      
      TH1F *htemp = (TH1F*)gPad->GetPrimitive(name);   
      htemp->SetTitle("");
      htemp->GetXaxis()->SetTitle("Target #phi (mrad)");
      htemp->GetYaxis()->SetTitle("Entries");
      stat = htemp->GetEntries();

      cout << " SSTAT = " << stat << endl;

      
      Double_t hrad = 0.0;
      
      
      if ((n_col == 14 && n_row == 8) || (n_col == 6 && n_row == 12)){
	// two larger holes
	hrad = SieveRadius_c;
      }
      else{
	hrad = SieveRadius;
      }
    


      Double_t Sieve_Z = MomDirectionTCS_hole.Z(); // distance from target to sieve
      
      Double_t hole_sig = TMath::ATan(hole_sig/(Sieve_Z)); // convert to angular sigma
      cout << "hole_sig = " << hole_sig << endl;

      Double_t hrad_ang = TMath::ATan(hrad/Sieve_Z); // angular 'radius'
      cout << "hrad_ang = " << hrad_ang << endl;
      
      Double_t ph_lowerlim = htemp->GetBinCenter(htemp->GetXaxis()->GetFirst());
      Double_t ph_upperlim = htemp->GetBinCenter(htemp->GetXaxis()->GetLast());
            
      // Double_t ph_lowlim = (ph_exp-hrad_ang)*1000;
      // Double_t ph_uplim = (ph_exp+hrad_ang)*1000;
      
      //cout << "ph_low = " << ph_lowlim << ", ph_uplim = " << ph_uplim << endl;


      //      TF1* fCirc = new TF1("fCirc",Form("2*sqrt(%f^2-(x-[0])^2)",hrad_ang),ph_lowerlim,ph_upperlim);
      TF1* fCirc = new TF1("fCirc",CircFunc,ph_lowerlim,ph_upperlim,3);
      cout << "Made fCirc" << endl;
      //      TF1* fGaus = new TF1("fGaus","gaus(0)");
      TF1* fGaus = new TF1("fGaus",GausFunc,-10,10,1);
      
      
      TF1Convolution *f_conv_ph= new TF1Convolution(fCirc,fGaus,ph_lowerlim,ph_upperlim,true);
      cout << "Made convolution" << endl;
      f_conv_ph->SetRange(ph_lowerlim,ph_upperlim);
      f_conv_ph->SetNofPointsFFT(10000);

      cout << "about to make fph" << endl;
      TF1   *fph = new TF1("fph",f_conv_ph,ph_lowerlim,ph_upperlim, f_conv_ph->GetNpar());

      cout << "Number of parameters in f_conv_ph = " << f_conv_ph->GetNpar() << endl;


      fph->SetParName(0,"Hole Radius");
      fph->SetParName(1,"Hole Centre");
      fph->SetParName(2,"Height");
      fph->SetParName(3,"Gaus #sigma");
           
      cout << "Number of parameters in fph = " << fph->GetNpar() << endl;

      // guess for hole centre: take as mean of plots
      Double_t ph_guess = htemp->GetMean();
      cout << "ph_lower = " << ph_lowerlim << ", ph_upperlim = " << ph_upperlim << endl;
      cout << "ph_guess = " << ph_guess << endl;
      
      fph->FixParameter(0,hrad_ang*1000);      
      //      fph->SetParameter(1,ph_exp*1000);
      fph->SetParameter(1,ph_guess);
      fph->SetParameter(2,0.5*htemp->GetBinContent(htemp->GetMaximumBin()));
      fph->SetParameter(3,0.35);

      cout << "about to fit" << endl;
      htemp->Fit("fph","qR");

      
      fCirc->SetLineColor(kRed);
      fCirc->SetParameter(0,fph->GetParameter(0));
      fCirc->SetParameter(1,fph->GetParameter(1));
      fCirc->SetParameter(2,fph->GetParameter(2));
      fCirc->Draw("same");

      
      
      phi_rms = fph->GetParameter(3);
      phi_cen = fph->GetParameter(1);

      cout << "phi_cen = " << phi_cen << endl;
      cout << "phi_rms = " << phi_rms << endl;
      

      c2[n_canvas] = new TCanvas(Form("c2_%d",n_canvas),"",800,600);
      //      t->Draw("x_sieve*1000>>" + name2, GeneralCut && id_cut,""); //


      
      // t->Draw(Form("(th_tgt-%f)*1000>>",th_exp) + name2, GeneralCut && id_cut,""); //  altered JW
      t->Draw("R.tr.tg_th*1000>>" + name2, GeneralCut && id_cut,""); //  altered JW
      
      

      TH1F *htemp2 = (TH1F*)gPad->GetPrimitive(name2);

      htemp2->SetTitle("");
      htemp2->GetXaxis()->SetTitle("Target #theta (mrad)");
      htemp2->GetYaxis()->SetTitle("Entries");


      Double_t th_lowerlim = htemp2->GetBinCenter(htemp2->GetXaxis()->GetFirst());
      Double_t th_upperlim = htemp2->GetBinCenter(htemp2->GetXaxis()->GetLast());

      TF1* fCirc2 = new TF1("fCirc2",CircFunc,th_lowerlim,th_upperlim,3);

      TF1* fGaus2 = new TF1("fGaus2",GausFunc,-10,10,1);
      TF1Convolution *f_conv_th= new TF1Convolution(fCirc2,fGaus2,th_lowerlim,th_upperlim,true);
      f_conv_th->SetRange(th_lowerlim,th_upperlim);
      f_conv_th->SetNofPointsFFT(10000);

      TF1   *fth = new TF1("fth",f_conv_th,th_lowerlim,th_upperlim, f_conv_th->GetNpar());
      
      fth->SetParName(0,"Hole Radius");
      fth->SetParName(1,"Hole Centre");
      fth->SetParName(2,"Height");
      fth->SetParName(3,"Gaus #sigma");

      // guess for hole centre: take as mean of plots
      Double_t th_guess = htemp2->GetMean();      
      
      fth->FixParameter(0,hrad_ang*1000);
      //      fth->SetParameter(1,th_exp*1000);
      fth->SetParameter(1,th_guess);
      fth->SetParameter(2,0.5*htemp2->GetBinContent(htemp2->GetMaximumBin()));
      fth->SetParameter(3,1.5);
      htemp2->Fit("fth","qR");

      
      fCirc2->SetLineColor(kRed);
      fCirc2->SetParameter(0,fth->GetParameter(0));
      fCirc2->SetParameter(1,fth->GetParameter(1));
      fCirc2->SetParameter(2,fth->GetParameter(2));
      fCirc2->Draw("same");

      
      
      th_rms = fth->GetParameter(3);
      th_cen = fth->GetParameter(1);


      // adjust value for rms by taking into account rms of sieve hole

      // calculate rms of sieve hole





      Double_t th_rms_cor = TMath::Sqrt( (th_rms*1e-3)*(th_rms*1e-3) - hole_sig*hole_sig);
      th_rms_cor *= 1e3; // convert from rad to mrad
      Double_t phi_rms_cor = TMath::Sqrt( (phi_rms*1e-3)*(phi_rms*1e-3) - hole_sig*hole_sig);
      phi_rms_cor *= 1e3; // convert from rad to mrad
      
      
      //cutcsvnew<<1<<","<<phi_cen<<","<<ph_exp*1000<<","<<th_cen<<","<<th_exp*1000<<","<<phi_rms_cor<<","<<th_rms_cor<<","<<stat<<","<<n_col<<","<<n_row<<endl;
      cutcsvnew<<1<<","<<phi_cen<<","<<ph_exp*1000<<","<<th_cen<<","<<th_exp*1000<<","<<phi_rms<<","<<th_rms<<","<<stat<<","<<n_col<<","<<n_row<<endl;    

      n_canvas++;
      }
  }
  

  int n = 0;
  
  while(c[n] != NULL) n++;

  TString output_phi = Form("/projections/phi_foil_%i.pdf",foil_no);
  TString output_theta = Form("/projections/theta_foil_%i.pdf",foil_no);
  
  for(int i = 0; i<n; i++){
    if(i == 0) {
      c[i]->Print(output + output_phi + "(","pdf");
      c2[i]->Print(output + output_theta + "(","pdf");
    }
    else if(i == n-1) {
      c[i]->Print(output + output_phi + ")","pdf");
      c2[i]->Print(output + output_theta + ")","pdf");
    }
    else {
      c[i]->Print(output + output_phi,"pdf");
      c2[i]->Print(output + output_theta,"pdf");
    }
  }
}
