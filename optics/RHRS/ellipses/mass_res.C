
vector<TString> run = {"4647","4648","4650"};
//vector<TString> run = {"4647"};
TString range = "-10_10";

double Avg(vector<double> R){
  double average = 0;
  
  int number = 0;

  for(int i = 0; i<R.size(); i++)
    if(abs(R[i]) < 1000) number++;
  
  for(int i = 0; i<R.size(); i++) {
    if(abs(R[i]) > 1000) continue;
    average += R[i]/number;

  }
  
  return average;
}

void compare(double &ph_wid,double &th_wid,double &ph_off,double &th_off){

  vector<double> y_th_tot, y_ph_tot, ph_std_tot, th_std_tot;
  
  //  for(int i_run = 0;i_run < run.size();i_run++){
  
  
    //ifstream file(run[i_run] + "_test/xfp_" + range + "/apex_" + run[i_run] + ".root.EllipseCuts.csv");
    ifstream file("Opt_All_Ellipses/xfp_" + range + "/Opt_All_Ellipses.root.EllipseCuts.csv");
  double ph, th, ph_exp, th_exp, ph_rms, th_rms;
  int hole_n = 1, opt;
  vector<double> x,y_th, y_ph, ph_std, th_std;
  
  string line;

  getline(file,line); 
  getline(file,line);
  while (getline(file,line)){
    stringstream linestream(line);
    string cell;
    
    int n_elem = 1;
    while(getline(linestream,cell,',')){
      if(n_elem == 3) opt = stoi(cell);
      else if(n_elem == 4) ph = stod(cell);
      else if(n_elem == 5) ph_exp = stod(cell);
      else if(n_elem == 6) th = stod(cell);
      else if(n_elem == 7) th_exp = stod(cell);
      else if(n_elem == 8) ph_rms = stod(cell);
      else if(n_elem == 9) th_rms = stod(cell);
      n_elem++;
    }

    
    if(opt) {
      y_th.push_back(th - th_exp);
      y_ph.push_back(ph - ph_exp);
      ph_std.push_back(ph_rms);
      th_std.push_back(th_rms);
      x.push_back(hole_n);

      y_th_tot.push_back(th - th_exp);
      y_ph_tot.push_back(ph - ph_exp);
      ph_std_tot.push_back(ph_rms);
      th_std_tot.push_back(th_rms);
    }
    hole_n++;
  }
  
  TCanvas* c = new TCanvas("c","",800,640);
  TGraphErrors* g = new TGraphErrors(x.size(),&x[0],&y_th[0],0,&th_std[0]);
  g->Draw("AP");
  g->SetMarkerStyle(20);
  g->GetYaxis()->SetRangeUser(-2.5,2.5);
  g->GetXaxis()->SetTitle("Hole #");
  g->GetYaxis()->SetTitle("#theta_{meas} - #theta_{exp} (mrad)");
  g->SetTitle("Measured Sieve Hole Centers");

  
  TLine* l = new TLine(x[0] - 30,0,x[x.size()-1] + 30,0);
  l->Draw("same");
  l->SetLineColor(kRed);

  
  //c->SaveAs(run[i_run] + "_test/xfp_" + range + "/theta_diff.gif");
  c->SaveAs("./Opt_All_Ellipses/xfp_" + range + "/theta_diff.gif");

  TCanvas* c2 = new TCanvas("c2","",800,640);
  TGraphErrors* g2 = new TGraphErrors(x.size(),&x[0],&y_ph[0],0,&ph_std[0]);
  g2->Draw("AP");
  g2->SetMarkerStyle(20);
  g2->GetYaxis()->SetRangeUser(-2.5,2.5);
  g2->GetXaxis()->SetTitle("Hole #");
  g2->GetYaxis()->SetTitle("#phi_{meas} - #phi_{exp} (mrad)");
  g2->SetTitle("Measured Sieve Hole Centers");
  
  l->Draw("same");
  //c2->SaveAs(run[i_run] + "_test/xfp_" + range + "/phi_diff.gif");
  c2->SaveAs("./Opt_All_Ellipses/xfp_" + range + "/phi_diff.gif");


  TH1F* hphi = new TH1F("hphi","#phi Offsets;#phi_{meas} - #phi_{exp} (mrad);Entries",100,-3.5,3.5);

  for(int i=0;i<y_ph.size();i++) hphi->Fill(y_ph[i]);

  TCanvas* c5 = new TCanvas("c5","",800,640);
  hphi->Draw();
  hphi->GetYaxis()->SetTitleOffset(1.0);

  //c5->SaveAs(run[i_run] + "_test/xfp_" + range + "/phi_offset.gif");
  c5->SaveAs("./Opt_All_Ellipses/xfp_" + range + "/phi_offset.gif");

  
  TH1F* htheta = new TH1F("htheta","#theta Offsets;#theta_{meas} - #theta_{exp} (mrad);Entries",100,-3.5,3.5);

  for(int i=0;i<y_ph.size();i++) htheta->Fill(y_th[i]);

  TCanvas* c6 = new TCanvas("c6","",800,640);
  htheta->Draw();
  htheta->GetYaxis()->SetTitleOffset(1.0);
 
  //c6->SaveAs(run[i_run] + "_test/xfp_" + range + "/theta_offset.gif");
  c6->SaveAs("./Opt_All_Ellipses/xfp_" + range + "/theta_offset.gif");
  

  TH1F* hphi_std = new TH1F("hphi_std","#phi Width;#phi Width (mrad);Entries",100,0,2.0);

  for(int i=0;i<ph_std.size();i++) hphi_std->Fill(ph_std[i]);

  TCanvas* c7 = new TCanvas("c7","",800,640);
  hphi_std->Draw();
  hphi_std->GetYaxis()->SetTitleOffset(1.0);

  //c7->SaveAs(run[i_run] + "_test/xfp_" + range + "/phi_std.gif");
  c7->SaveAs("./Opt_All_Ellipses/xfp_" + range + "/phi_std.gif");
  
  TH1F* htheta_std = new TH1F("htheta_std","#theta Width;#theta Width (mrad);Entries",100,0.5,3.5);

  for(int i=0;i<th_std.size();i++) htheta_std->Fill(th_std[i]);

  TCanvas* c8 = new TCanvas("c8","",800,640);
  htheta_std->Draw();
  htheta_std->GetYaxis()->SetTitleOffset(1.0);
  
  //c8->SaveAs(run[i_run] + "_test/xfp_" + range + "/theta_std.gif");
  c8->SaveAs("./Opt_All_Ellipses/xfp_" + range + "/theta_std.gif");

  //}
  
  ph_off = Avg(y_ph_tot)/1000;
  th_off = Avg(y_th_tot)/1000;
  ph_wid = Avg(ph_std_tot)/1000;
  th_wid = Avg(th_std_tot)/1000;
  
  
}


double mass_calc(double theta_p,double theta_m,double phi_p,double phi_m,double delta_p,double delta_m){

  double p0 = 1104; ///MeV
  double theta0 = 5.0*TMath::Pi()/180; ///HRS angle = 5 deg
  

  
  
  return sqrt(p0*p0*(4*theta0*theta0 + 4*theta0*theta0*delta_p + 4*theta0*theta0*delta_m + 4*theta0*(phi_p - phi_m) + 2*theta_p*theta_m));
  //return sqrt(p0*p0*(4*theta0*theta0 + 8*theta0*(phi_p - phi_m) + 2*theta_p*theta_m));

}


void mass_res(){

  gStyle->SetOptStat("emr");
  
  double ph_wid, th_wid, ph_off, th_off;

  compare(ph_wid, th_wid, ph_off,th_off);
    
  int n_events = 1000000;

  double ph_low = 0;
  double ph_high = 20./1000;
  double th_low = -40./1000;
  double th_high = 40./1000;
  double dp_lim = 0.045;
  
  
  TH1F* hmass = new TH1F("hmass","Mass Resolution;#Delta m (MeV);Entries",100,-10,10);
  hmass->GetYaxis()->SetTitleOffset(1.6);
  
  TRandom2 *tr = new TRandom2();

  
  for(int i=0; i<n_events; i++){

    double ph_p = tr->Uniform(ph_low,ph_high);
    double ph_m = ph_p;
    double th_p = tr->Uniform(th_low,th_high);
    double th_m = -th_p;
    double dp_p = tr->Uniform(-dp_lim,dp_lim);
    //dp_p = 0;
    double dp_m = -dp_p;

  
   
    double m1 = mass_calc(th_p,th_m,ph_p,ph_m,dp_p,dp_m);

    //    cout<<m1<<endl;
    
    ph_p = tr->Gaus(ph_p,ph_wid);
    ph_m = tr->Gaus(ph_m,ph_wid);
    th_p = tr->Gaus(th_p,th_wid);
    th_m = tr->Gaus(th_m,th_wid);
    dp_p = tr->Gaus(dp_p,0.001);
    dp_m = tr->Gaus(dp_m,0.001);


   
    double m2 = mass_calc(th_p,th_m,ph_p,ph_m,dp_p,dp_m);


    //cout<<m1<<" "<<m2<<" "<<m1-m2<<endl;
    hmass->Fill(m1 - m2);
  }

  TCanvas *cmass = new TCanvas("cmass","",800,640);
  hmass->Draw();

  cmass->SaveAs("./mass_res.gif");

  cout<<"Phi Offset = "<<ph_off*1000<<" mrad"<<endl;
  cout<<"Phi Width = "<<ph_wid*1000<<" mrad"<<endl;
  cout<<"Theta Offset = "<<th_off*1000<<" mrad"<<endl;
  cout<<"Theta Width = "<<th_wid*1000<<" mrad"<<endl;
  }
