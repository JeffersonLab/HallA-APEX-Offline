
////Macro calculates the angular resolutions and the mass resolution /////

TString range;

double Avg(vector<double> R){
  double average = 0;
  
  for(int i = 0; i<R.size(); i++) average += R[i]/R.size();
    
  return average;
}

void compare(double &ph_wid,double &th_wid,double &ph_off,double &th_off){

  range = "full_brute";

  ifstream file("xfp_" + range + "/apex_4647.root.EllipseCuts.csv");
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
      else if(n_elem == 11) ph_rms = stod(cell);
      else if(n_elem == 12) th_rms = stod(cell);
      n_elem++;
    }
    
    if(opt) {
      y_th.push_back(th - th_exp);
      y_ph.push_back(ph - ph_exp);
      ph_std.push_back(ph_rms);
      th_std.push_back(th_rms);
      x.push_back(hole_n);
      hole_n++;
    }
  }
  
  TCanvas* c = new TCanvas("c","",800,640);
  TGraphErrors* g = new TGraphErrors(x.size(),&x[0],&y_th[0],0,&th_std[0]);
  g->Draw("AP");
  g->SetMarkerStyle(20);
  g->GetYaxis()->SetRangeUser(-2.5,2.5);
  g->GetXaxis()->SetTitle("Hole #");
  g->GetYaxis()->SetTitle("#theta_{meas} - #theta_{exp} (mrad)");
  g->SetTitle("Measured Sieve Hole Centers");

  TLine* l = new TLine(0,0,78,0);
  l->Draw("same");
  l->SetLineColor(kRed);

  c->SaveAs("xfp_" + range + "/theta_diff.gif");

  TCanvas* c2 = new TCanvas("c2","",800,640);
  TGraphErrors* g2 = new TGraphErrors(x.size(),&x[0],&y_ph[0],0,&ph_std[0]);
  g2->Draw("AP");
  g2->SetMarkerStyle(20);
  g2->GetYaxis()->SetRangeUser(-2.5,2.5);
  g2->GetXaxis()->SetTitle("Hole #");
  g2->GetYaxis()->SetTitle("#phi_{meas} - #phi_{exp} (mrad)");
  g2->SetTitle("Measured Sieve Hole Centers");
  
  l->Draw("same");
  c2->SaveAs("xfp_" + range + "/phi_diff.gif");


  TCanvas* c3 = new TCanvas("c3","",800,640);
  TGraph* g3 = new TGraph(x.size(),&x[0],&ph_std[0]);
  g3->Draw("AP");
  g3->SetMarkerStyle(20);
  g3->GetYaxis()->SetRangeUser(0,1.5);
  g3->GetXaxis()->SetTitle("Hole #");
  g3->GetYaxis()->SetTitle("#phi Width (mrad)");
  g3->SetTitle("Width of #phi");


  TCanvas* c4 = new TCanvas("c4","",800,640);
  TGraph* g4 = new TGraph(x.size(),&x[0],&th_std[0]);
  g4->Draw("AP");
  g4->SetMarkerStyle(20);
  g4->GetYaxis()->SetRangeUser(0.5,3.5);
  g4->GetXaxis()->SetTitle("Hole #");
  g4->GetYaxis()->SetTitle("#theta Width (mrad)");
  g4->SetTitle("Width of #theta");

  TH1F* hphi = new TH1F("hphi","#phi Offsets;#phi_{meas} - #phi_{exp} (mrad);Entries",100,-1.5,1.5);

  for(int i=0;i<y_ph.size();i++) hphi->Fill(y_ph[i]);

  TCanvas* c5 = new TCanvas("c5","",800,640);
  hphi->Draw();
  hphi->GetYaxis()->SetTitleOffset(1.0);

  c5->SaveAs("xfp_" + range + "/phi_offset.gif");

  
  TH1F* htheta = new TH1F("htheta","#theta Offsets;#theta_{meas} - #theta_{exp} (mrad);Entries",100,-1.5,1.5);

  for(int i=0;i<y_ph.size();i++) htheta->Fill(y_th[i]);

  TCanvas* c6 = new TCanvas("c6","",800,640);
  htheta->Draw();
  htheta->GetYaxis()->SetTitleOffset(1.0);
 
  c6->SaveAs("xfp_" + range + "/theta_offset.gif");
  

  TH1F* hphi_std = new TH1F("hphi_std","#phi Width;#phi Width (mrad);Entries",100,0,1.4);

  for(int i=0;i<ph_std.size();i++) hphi_std->Fill(ph_std[i]);

  TCanvas* c7 = new TCanvas("c7","",800,640);
  hphi_std->Draw();
  hphi_std->GetYaxis()->SetTitleOffset(1.0);

  c7->SaveAs("xfp_" + range + "/phi_std.gif");
  
  TH1F* htheta_std = new TH1F("htheta_std","#theta Width;#theta Width (mrad);Entries",100,0.5,3.5);

  for(int i=0;i<th_std.size();i++) htheta_std->Fill(th_std[i]);

  TCanvas* c8 = new TCanvas("c8","",800,640);
  htheta_std->Draw();
  htheta_std->GetYaxis()->SetTitleOffset(1.0);
  
  c8->SaveAs("xfp_" + range + "/theta_std.gif");

  
  ph_off = Avg(y_ph)/1000;
  th_off = Avg(y_th)/1000;
  ph_wid = Avg(ph_std)/1000;
  th_wid = Avg(th_std)/1000;
  
}


double mass_calc(double theta_p,double theta_m,double phi_p,double phi_m,double delta_p,double delta_m){

  double p0 = 1104; ///MeV
  double theta0 = 5.0*TMath::Pi()/180; ///HRS angle = 5 deg

  //cout<<delta_p<<" "<<delta_m<<endl;
  return sqrt(p0*p0*(4*theta0*theta0 + 4*theta0*theta0*delta_p + 4*theta0*theta0*delta_m + 8*theta0*(phi_p - phi_m) + 2*theta_p*theta_m));
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

  cmass->SaveAs("xfp_" + range + "/mass_res.gif");
  
}
