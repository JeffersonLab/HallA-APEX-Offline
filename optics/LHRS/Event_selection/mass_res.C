
#include "file_def.h"

TString range;

TString Sieve_CSV_pos = "proj_csv/";

double Avg(vector<double> R){
  double average = 0;
  
  for(int i = 0; i<R.size(); i++) average += R[i]/R.size();
    
  return average;
}

void compare(double &ph_wid,double &th_wid,double &ph_off,double &th_off, TString option){

  range = "full_brute";


  //  ifstream file(Sieve_CSV_pos + "/EllipseCuts.csv");

  TString new_tstring =  CutFileName;
  new_tstring.Remove(0,30);

  TString file_name;
  
  if( option == ""){
    file_name = "proj_csv/"  + new_tstring + ".csv";
  }
  else{
    file_name = "proj_csv/"  + new_tstring + "_" + option + ".csv";
  }
  
 
  ifstream file(file_name);
  
  //  ifstream file("xfp_" + range + "/apex_4647.root.EllipseCuts.csv");

  
  double ph, th, ph_exp, th_exp, ph_rms, th_rms;

  int hole_n = 1, opt, col, row;
  vector<double> x,y_th, y_ph, ph_std, th_std;
  vector<double> x_col[NSieveCol];
  vector<pair<double,double>> y_th_col[NSieveCol], y_ph_col[NSieveCol], ph_std_col[NSieveCol], th_std_col[NSieveCol]; // vector filled with row number and value
  vector<double> x_row_hole[NSieveRow], x_row_column[NSieveCol]; // keeps track of hole and columns for each row
  vector<double> y_th_row[NSieveRow], y_ph_row[NSieveRow], ph_std_row[NSieveRow], th_std_row[NSieveRow]; // vector filled with column number and value


  // histograms to store results for columns and rows

  TH1F* h_row_y[NSieveRow];
  TH1F* h_col_x[NSieveCol];

  // for (Int_t i_row; i_row<NSieveRow; i_row++){   
  //   h_row_y[i_row] = new TH1F(Form("h theta %i", i_row),"#theta;Hole #;#theta_{meas} - #theta_{exp} (mrad)",);
  // }
  


  
  // store column and row of hole
  vector<int> cols,rows;
  
  string line;

  //  getline(file,line);
  getline(file,line);
  while (getline(file,line)){
    stringstream linestream(line);
    string cell;

    int n_elem = 1;
    while(getline(linestream,cell,',')){
      if(n_elem == 1) opt = stoi(cell);
      else if(n_elem == 2) ph = stod(cell);
      else if(n_elem == 3) ph_exp = stod(cell);
      else if(n_elem == 4) th = stod(cell);
      else if(n_elem == 5) th_exp = stod(cell);
      else if(n_elem == 6) ph_rms = stod(cell);
      else if(n_elem == 7) th_rms = stod(cell);
      else if(n_elem == 9) col = stod(cell);
      else if(n_elem == 10) row = stod(cell);
      n_elem++;
    }
    
    if(opt && !(abs(ph)>1e4 || abs(ph_rms)>1e4 || abs(th)>1e4 || abs(th_rms)>1e4) ) {
      y_th.push_back(th - th_exp);
      y_ph.push_back(ph - ph_exp);
      ph_std.push_back(ph_rms);
      th_std.push_back(th_rms);
      x.push_back(hole_n);
      cols.push_back(col);
      rows.push_back(row);

      
      y_th_col[col].push_back({row, th - th_exp});
      y_ph_col[col].push_back({row, ph - ph_exp});
      ph_std_col[col].push_back({row, ph_rms});
      th_std_col[col].push_back({row, th_rms});
      x_col[col].push_back(hole_n);

      
      y_th_row[row].push_back(th - th_exp);
      y_ph_row[row].push_back(ph - ph_exp);
      ph_std_row[row].push_back(ph_rms);
      th_std_row[row].push_back(th_rms);
      x_row_column[row].push_back(col);
      x_row_hole[row].push_back(hole_n);
      
      hole_n++;

      
    }
    else if(abs(ph)>1e4 || abs(ph_rms)>1e4 || abs(th)>1e4 || abs(th_rms)>1e4){

      cout << "col = " << col << ", row = " << row << endl;
      cout << "ph = " << ph << ", ph_rms = " << ph_rms << ", th = " << th << ", th_rms = " << th_rms << endl << endl;
    }


  }


  // retrieve largest and smallest rows and columns for plotting purposes

  int last_row = *std::max_element(rows.begin(),rows.end());
  int first_row = *std::min_element(rows.begin(),rows.end());

  int last_col = *std::max_element(cols.begin(),cols.end());
  int first_col = *std::min_element(cols.begin(),cols.end());

  
  
  // for plotting purposes choose series of colours and markers to be used to distinguish rows and columns

  // marker type used to distinguish rows
  Int_t marker_type[] =  {kFullCircle, kFullTriangleUp, kFullSquare, kFullDiamond, kFullTriangleDown, kFullFourTrianglesPlus, kFullThreeTriangles, kFullCross, kFullCrossX, kFourSquaresPlus, kCircle,  kMultiply, kOpenSquare, kOpenCircle,  kPlus};

  // colour used to distinguish columns
  Int_t col_colour[11] = {kBlue,kOrange+7,kGreen+4,kYellow-6,kMagenta,kSpring+10,kGray, kRed, kRed-5, kBlue+4, kBlack};

  // plot of reconsructed theta compared to survey
  
  TCanvas* c = new TCanvas("c","",800,640);
  
  // TGraphErrors* g = new TGraphErrors(x.size(),&x[0],&y_th[0],0,&th_std[0]);
  TMultiGraph* g_th_all =  new TMultiGraph;
  TGraphErrors* g_th_rc[NSieveRow][NSieveCol];

  Int_t row_count = 0;


  // create legend
  TLegend* leg = new TLegend(.1,.65,.37,.9,"Key");
  leg->SetFillColor(0);

  //  TH1I* h_marke[];

  

  
  for(auto row_v : y_th_row){
    
    if(!row_v.empty()){
      for( Int_t col_v = 0; col_v < x_row_column[row_count].size(); col_v++){      
	g_th_rc[row_count][col_v] = new TGraphErrors(1,&x_row_hole[row_count][col_v],&y_th_row[row_count][col_v],0,&th_std_row[row_count][col_v]);
	g_th_rc[row_count][col_v]->SetMarkerStyle(marker_type[row_count-first_row]);
	g_th_rc[row_count][col_v]->SetMarkerColor(col_colour[((int)x_row_column[row_count][col_v]-first_col)]);
	g_th_rc[row_count][col_v]->SetMarkerSize(2);
	g_th_all->Add(g_th_rc[row_count][col_v],"P");

	// if(!leg_added_col[x_row_column[row_count][col_v]]){
	//   leg_added_col = true;
	//   leg->Add();
	// }
	  
	//	}
      }
    }
    row_count ++;
  }
  
  g_th_all->Draw("AP");
  

  g_th_all->GetYaxis()->SetRangeUser(-3.0,3.0);
  g_th_all->GetXaxis()->SetTitle("Hole #");
  g_th_all->GetYaxis()->SetTitle("#theta_{meas} - #theta_{exp} (mrad)");
  g_th_all->SetTitle("Measured Sieve Hole Centers");

  TLine* l = new TLine(0,0,78,0);
  l->Draw("same");
  l->SetLineColor(kRed);

  new_tstring.Remove(22,38);
  
  c->SaveAs("proj_plots/" + new_tstring + "_theta_diff.pdf");


  // TCanvas* c_b = new TCanvas("c_b","",800,640);

  
  // TGraphErrors* g = new TGraphErrors(x.size(),&x[0],&y_th[0],0,&th_std[0]);
  // g->Draw("AP");
  // g->SetMarkerStyle(20);
  // g->GetYaxis()->SetRangeUser(-2.5,2.5);
  // g->GetXaxis()->SetTitle("Hole #");
  // g->GetYaxis()->SetTitle("#theta_{meas} - #theta_{exp} (mrad)");
  // g->SetTitle("Measured Sieve Hole Centers");

  // TLine* l = new TLine(0,0,78,0);
  // l->Draw("same");
  // l->SetLineColor(kRed);

  // c->SaveAs("proj_plots/theta_diff.pdf");


  

  // Phi plots

  
  // TCanvas* c2 = new TCanvas("c2","",800,640);
  // TGraphErrors* g2 = new TGraphErrors(x.size(),&x[0],&y_ph[0],0,&ph_std[0]);
  // g2->Draw("AP");
  // g2->SetMarkerStyle(20);
  // g2->GetYaxis()->SetRangeUser(-2.5,2.5);
  // g2->GetXaxis()->SetTitle("Hole #");
  // g2->GetYaxis()->SetTitle("#phi_{meas} - #phi_{exp} (mrad)");
  // g2->SetTitle("Measured Sieve Hole Centers");
  
  // l->Draw("same");
  // c2->SaveAs("proj_plots/phi_diff.pdf");


  TCanvas* c2 = new TCanvas("c2","",800,640);
  //  TGraphErrors* g2 = new TGraphErrors(x.size(),&x[0],&y_ph[0],0,&ph_std[0]);

  // loop through row vectors to see which have non-zero entries

  TMultiGraph* g2_ph_all =  new TMultiGraph;

  TGraphErrors* g2_ph_rc[NSieveRow][NSieveCol];

  row_count = 0;
  for(auto row_v : y_ph_row){
    
    if(!row_v.empty()){
      for( Int_t col_v = 0; col_v < x_row_column[row_count].size(); col_v++){
     	g2_ph_rc[row_count][col_v] = new TGraphErrors(1,&x_row_hole[row_count][col_v],&y_ph_row[row_count][col_v],0,&ph_std_row[row_count][col_v]);
	g2_ph_rc[row_count][col_v]->SetMarkerStyle(marker_type[row_count-first_row]);
	g2_ph_rc[row_count][col_v]->SetMarkerColor(col_colour[((int)x_row_column[row_count][col_v]-first_col)%8] );
	g2_ph_rc[row_count][col_v]->SetMarkerSize(2);
	g2_ph_all->Add(g2_ph_rc[row_count][col_v],"P");
      }
    }
    row_count ++;
  }

  g2_ph_all->Draw("AP");
  
  // g2->Draw("AP");
  // g2->SetMarkerStyle(20);
  g2_ph_all->GetYaxis()->SetRangeUser(-2.5,2.5);
  g2_ph_all->GetXaxis()->SetTitle("Hole #");
  g2_ph_all->GetYaxis()->SetTitle("#phi_{meas} - #phi_{exp} (mrad)");
  g2_ph_all->SetTitle("Measured Sieve Hole Centers");
  
  l->Draw("same");
  c2->SaveAs("proj_plots/" + new_tstring + "_phi_diff.pdf");



  


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

  c5->SaveAs("proj_plots/" + new_tstring + "_phi_offset.pdf");

  
  TH1F* htheta = new TH1F("htheta","#theta Offsets;#theta_{meas} - #theta_{exp} (mrad);Entries",100,-1.5,1.5);

  for(int i=0;i<y_ph.size();i++) htheta->Fill(y_th[i]);

  TCanvas* c6 = new TCanvas("c6","",800,640);
  htheta->Draw();
  htheta->GetYaxis()->SetTitleOffset(1.0);
 
  c6->SaveAs("proj_plots/" + new_tstring + "_theta_offset.pdf");
  

  TH1F* hphi_std = new TH1F("hphi_std","#phi Width;#phi Width (mrad);Entries",100,0,1.4);

  for(int i=0;i<ph_std.size();i++) hphi_std->Fill(ph_std[i]);

  TCanvas* c7 = new TCanvas("c7","",800,640);
  hphi_std->Draw();
  hphi_std->GetYaxis()->SetTitleOffset(1.0);

  c7->SaveAs("proj_plots/" + new_tstring + "_phi_std.pdf");
  
  TH1F* htheta_std = new TH1F("htheta_std","#theta Width;#theta Width (mrad);Entries",100,0.5,3.5);

  for(int i=0;i<th_std.size();i++) htheta_std->Fill(th_std[i]);

  TCanvas* c8 = new TCanvas("c8","",800,640);
  htheta_std->Draw();
  htheta_std->GetYaxis()->SetTitleOffset(1.0);
  
  c8->SaveAs("proj_plots/" + new_tstring + "_theta_std.pdf");


  
  double max_ph = *std::max_element(y_ph.begin(),y_ph.end());
  double max_ph_off = *std::max_element(ph_std.begin(),ph_std.end());

  cout << "Largest ph = " << max_ph << endl;
  cout << "Largest ph off = " << max_ph_off << endl;
  
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



void mass_res(TString option = ""){

  gStyle->SetOptStat("emr");
  
  double ph_wid, th_wid, ph_off, th_off;
  
  compare(ph_wid, th_wid, ph_off,th_off,option);

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

  TString new_tstring =  CutFileName;
  new_tstring.Remove(0,30);
  new_tstring.Remove(22,38);

  cmass->SaveAs("proj_plots/" + new_tstring + "_mass_res.pdf");

  cout<<"Phi Offset = "<<ph_off*1000<<" mrad"<<endl;
  cout<<"Theta Offset = "<<th_off*1000<<" mrad"<<endl;
  cout<<"Phi Width = "<<ph_wid*1000<<" mrad"<<endl;
  cout<<"Theta Width = "<<th_wid*1000<<" mrad"<<endl;
  cout <<endl;
  cout<<"new Theta Offset = "<<1000*(TMath::ATan(th_off/(ZPos)))<<endl;
  cout<<"new Theta Width = "<<1000*(TMath::ATan(th_wid/(ZPos)))<<endl;
  
}
