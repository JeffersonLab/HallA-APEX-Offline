// Script to plot several plots for the horizontal foils (LHRS)
//

#include "../Load_more_rootfiles.C"
#include "InputAPEXL.h"




// using 'triangular wave' of fast raster to conver from current to time

Double_t curr_to_time(Double_t curr, Double_t *par, Bool_t positive = true){

  Double_t period = par[0];
  Double_t amplitude = par[1];
  Double_t base = par[2];


  // positive argument gives phase
  // true if phase means that dI/dt (derivative of current wrt time) is positive

  Double_t phase_add = 0;
  Double_t phase_mult= 1;

  if(positive){
    phase_add =  (0.5*period); // was -
    phase_mult = -1;
  }

  Double_t time = phase_add + phase_mult*(period/(2*TMath::Pi())) * TMath::ASin( TMath::Sin(TMath::Pi()*((curr-base))/(2*amplitude)  )  );

  if (time < -0.75* period){
    time = 0.25*period - TMath::Abs(time + 0.75*period);
  }

  return time;
  //  return phase_add + phase_mult*(period/(2*TMath::Pi())) * TMath::ASin( TMath::Sin(TMath::Pi()*((curr-base))/(2*amplitude)  )  );
  
}

Double_t time_to_curr(Double_t time, Double_t *par){

  Double_t period = par[0];
  Double_t amplitude = par[1];
  Double_t base = par[2];
  
  
  return ( ((2*amplitude)/TMath::Pi()) * TMath::ASin( TMath::Sin(2*TMath::Pi()*(time)/period)  ) + base );
  
}


Double_t peak(Double_t *x, Double_t *par) {
  return  par[0]*TMath::Exp(-(((x[0]-par[1])*(x[0]-par[1]))/(2*par[2]*par[2])));
}

Double_t bg(Double_t *x, Double_t *par) {
  return par[0]+ (par[1]*x[0]);
}

// sum of peaks and backgorund
Double_t overall(Double_t *x, Double_t *par) {
  return peak(x,par) + peak(x,&par[3]) + bg(x,&par[6]);
}



void horizontal_correction(){


  // load in root files
  
   Int_t run_nos[] = {4775,4776,4777};
  // Int_t run_nos[] = {4776,4777};
  
  TChain* T;


  Int_t run_po = 0;
  for(auto run_no : run_nos){

    TChain* T_add = Load_more_rootfiles(run_no);

       
    // *T_ind[run_po] = *T_add;
    
    if(run_no == run_nos[0]){
      T = T_add;
    }
    else{
      T->Add(T_add);
    }

    cout << "T->GetEntries() = " << T->GetEntries() << endl;

    run_po++;
  }


  // set histogram limits

  // raster limits
  Double_t Lrby_l = 0.0014;
  Double_t Lrby_h = 0.0030;

  Double_t Lrbx_l = -0.0022;
  Double_t Lrbx_h = -0.0011;

  // BPM limits
  Double_t BPM_lim_l = -0.001;
  Double_t BPM_lim_h = 0.006;

  // raster current limits
  Double_t Lrby_curr_l = 18000;
  Double_t Lrby_curr_h = 45000;


    
  // reactz limits
  Double_t reactz_l = -0.4;
  Double_t reactz_h = 0.4;

  // Focal plane limits
  
  Double_t y_FP_l = -0.05;
  Double_t y_FP_h = 0.05;

  Double_t ph_FP_l = -50;
  Double_t ph_FP_h = 40;

  Double_t th_FP_l = -80;
  Double_t th_FP_h = 80;

  // target cuts

  Double_t ph_tg_l = -40;
  Double_t ph_tg_h = 40;
  
  Double_t th_tg_l = -80;
  Double_t th_tg_h = 80;
  
  // define raster cuts
  Double_t c_Lrby1_l = 0.00191;
  Double_t c_Lrby1_h = 0.00209;

  Double_t c_Lrby2_l = 0.00240;
  Double_t c_Lrby2_h = 0.00258;

  // define BPM cuts (work in addition to Raster as phase cuts)
  Double_t c_BPMB1_l = -0.0001;
  Double_t c_BPMB1_h = 0.0005;

  Double_t c_BPMB2_l = 0.0043;
  Double_t c_BPMB2_h = 0.0049;

  // Double_t c_BPMB2_l = -1;
  // Double_t c_BPMB2_h = 1;
  
  
  // define reactz cut for V2
  Double_t c_zV2_l = -0.1;
  Double_t c_zV2_h = 0.05;

  // define V2 and raster cuts

  TCut V2_zcut = Form("reactz>%f && reactz<%f",c_zV2_l,c_zV2_h);

  // lower raster cut (in y)
  TCut rast_l = Form("Lrb.y>%f && Lrb.y<%f",c_Lrby1_l,c_Lrby1_h);

  // greater raster cut (in y)
  TCut rast_h = Form("Lrb.y>%f && Lrb.y<%f",c_Lrby2_l,c_Lrby2_h);

  // lower BPMB y
  TCut BPMB_l = Form("Lrb.BPMB.y>%f && Lrb.BPMB.y<%f",c_BPMB1_l,c_BPMB1_h);

  // higher BPMB y
  TCut BPMB_h = Form("Lrb.BPMB.y>%f && Lrb.BPMB.y<%f",c_BPMB2_l,c_BPMB2_h);

  rast_l += BPMB_h;

  rast_h += BPMB_l;
  
  // cuts from InputAPEXl.h
  TCut GenrealCut = GeneralSieveCut + PID_cuts + FP_cuts;
  
  TCanvas *c1 = new TCanvas("c1"); 

  cout << GenrealCut << endl << endl;


  TH1F* h_ycurr = new TH1F("h_ycurr","raster current Y",1000,Lrby_curr_l,Lrby_curr_h);
  h_ycurr->GetXaxis()->SetTitle("Raster Y current [ADC channels]");
  
  
  T->Draw("Lrb.Raster2.rawcur.y>>h_ycurr");


  // fit peaks in raster y current
  // rough values (min max range) are set here: could be automated with effort

  Double_t pmin_1 = 27500;
  Double_t pmax_1 = 28500;
  
  Double_t pmin_2 = 29000;
  Double_t pmax_2 = 30000;

  h_ycurr->Fit("gaus","Q","",pmin_1,pmax_1);

  TF1* gaus1 = h_ycurr->GetFunction("gaus");

  Double_t min_1 = gaus1->GetParameter(1) - 3*gaus1->GetParameter(2);
  Double_t max_1 = gaus1->GetParameter(1) + 3*gaus1->GetParameter(2);
    
  h_ycurr->Fit("gaus","Q","",min_1,max_1);

  gaus1 = h_ycurr->GetFunction("gaus");

  Double_t mean1 = gaus1->GetParameter(1);

  
  h_ycurr->Fit("gaus","Q","",pmin_2,pmax_2);

  TF1* gaus2 = h_ycurr->GetFunction("gaus");

  Double_t min_2 = gaus2->GetParameter(1) - 3*gaus2->GetParameter(2);
  Double_t max_2 = gaus2->GetParameter(1) + 3*gaus2->GetParameter(2);
    
  h_ycurr->Fit("gaus","Q","",min_2,max_2);

  gaus2 = h_ycurr->GetFunction("gaus");

  Double_t mean2 = gaus2->GetParameter(1);

  
  cout << "Mean for first peak = " << mean1 << endl;
  cout << "Mean for second peak = " << mean2 << endl;

  

  TF1* overall_fit = new TF1("Overall_fit",overall,pmin_1-2000,pmax_2+2000,8);

  overall_fit->SetParameter(1,mean1);
  overall_fit->SetParLimits(1,0.98*mean1,1.02*mean1);
  overall_fit->SetParameter(2,gaus1->GetParameter(2));
  //  overall_fit->SetParameter(3,10000);
  overall_fit->SetParameter(4,mean2);
  overall_fit->SetParLimits(4,0.98*mean2,1.02*mean2);
  overall_fit->SetParameter(5,gaus2->GetParameter(2));

  // estimate linear parameters
  overall_fit->SetParameter(6,2000);
  overall_fit->SetParameter(7,0.1);
  overall_fit->SetParLimits(7,-0.2,0.2);

  h_ycurr->Fit(overall_fit,"BR");
  
  if(overall_fit){
    mean1 = overall_fit->GetParameter(1);
    mean2 = overall_fit->GetParameter(4);
  }
  else{
    cout << "over_fit failed" << endl;
  }
  
  cout << "Mean for first peak = " << mean1 << endl;
  cout << "Mean for second peak = " << mean2 << endl;
  cout << "line params: y = " << overall_fit->GetParameter(7) << "x + " << overall_fit->GetParameter(6)  << endl;


  // plot Y curent versus BPMB y (gives phase info)

  TCanvas *c2 = new TCanvas("c2"); 

  TH2F* h_rcurr = new TH2F("h_rcurr","raster current vs BPMB Y",1000,BPM_lim_l,BPM_lim_h,1000,Lrby_curr_l,Lrby_curr_h);
  h_rcurr->GetXaxis()->SetTitle("BPMB Y [m]");
  h_rcurr->GetYaxis()->SetTitle("Raster Y current [ADC channels]");
  h_rcurr->GetYaxis()->SetTitleOffset(2.0);
    
  T->Draw("Lrb.Raster2.rawcur.y:Lrb.BPMB.y>>h_rcurr","","colz");
  

  // define line to use as cut for positive and negative phase
  Double_t l_x1 = 0.0017;
  Double_t l_y1 = 20100;
  Double_t l_x2 = 0.0035;
  Double_t l_y2 = 39400;
  
  //  TLine* phase_l = new TLine(0.0017,20100,0.0035,39400);
  TLine* phase_l = new TLine(l_x1,l_y1,l_x2,l_y2);
  phase_l->SetLineColor(kRed);
  phase_l->Draw("same");
  
  // convert this difference in current to a time delay

  Double_t curr_corr = TMath::Abs(mean2 - mean1)/2.0;

  Double_t period = 1/(25e3); // 25khz

  Double_t amplitude = 9525;
  Double_t base = 2.978e4; // baseline current
  
  cout << "amplitude = " << amplitude << endl;
  cout << "period = " << period << endl;


  // set wave parameters for converting currents to time (and vice versa)
  Double_t wave_pars[3] = {period,amplitude,base};
  
  //  Double_t time_1 = - (period/2) - (period/(2*TMath::Pi())) * TMath::ASin( TMath::Sin(TMath::Pi()*(mean1-base)/(2*amplitude)  )  );
  Double_t time_1 = curr_to_time(mean1,wave_pars,false);
  
  //  Double_t time_1_corr =  - (period/2) - (period/(2*TMath::Pi())) * TMath::ASin( TMath::Sin(TMath::Pi()*((mean1-base)+curr_corr)/(2*amplitude)  )  );
  Double_t time_1_corr = curr_to_time(mean1+curr_corr,wave_pars,false);

  Double_t time_1_diff = time_1 - time_1_corr;


  // test if adjusted current can be obtained from orginal peak
  Double_t curr_1_test =time_to_curr( curr_to_time(mean1,wave_pars,false) -time_1_diff, wave_pars);

  //  Double_t time_2 = (period/(2*TMath::Pi())) * TMath::ASin( TMath::Sin(TMath::Pi()*(mean2-base)/(2*amplitude)  )  );
  Double_t time_2 = curr_to_time(mean2,wave_pars);

  //  Double_t time_2_corr = (period/(2*TMath::Pi())) * TMath::ASin( TMath::Sin(TMath::Pi()*((mean2-base)-curr_corr)/(2*amplitude)  )  );
  Double_t time_2_corr = curr_to_time(mean2-curr_corr,wave_pars);

  Double_t time_2_diff = time_2 - time_2_corr;

  // test if adjusted current can be obtained from orginal peak
  Double_t curr_2_test =time_to_curr( curr_to_time(mean2,wave_pars,false) - time_1_diff, wave_pars);


  cout << "time_1 = " << time_1 << endl;
  cout << "time_2 = " << time_2 << endl;

  
  cout << "time_1_diff = " << time_1_diff << endl;
  cout << "time_2_diff = " << time_2_diff << endl;


  cout << "mean1+curr_corr = " << mean1+curr_corr << endl;

  cout << "mean2-curr_corr = " << mean2-curr_corr << endl;


  cout << "curr_1_test = " << curr_1_test << endl;
  cout << "curr_2_test = " << curr_2_test << endl << endl;
    



  // for loop to test saw distribution

  std::vector<Double_t> cur_vals;
  std::vector<Double_t> time_vals;

  TH2D* cur_wave = new TH2D("cur_wave","Raster signal",100,1e6*-0.6*period,1e6*0.6*period,100,base-1.1*amplitude,base+1.1*amplitude);

  
  //  for(Double_t cur_t = base; cur< (base+2*amplitude); cur_t = cur_t + (2*amplitude)/100){
  for(Double_t wave_t = -period; wave_t<period; wave_t = wave_t + period/100){


    Double_t wave_i = time_to_curr(wave_t,wave_pars);
    cur_wave->Fill(1e6*wave_t,wave_i);
    

  }

  TCanvas *c3 = new TCanvas("c3","c3",1200,1200);

  c3->SetLeftMargin(0.2);
  
  cur_wave->GetXaxis()->SetTitle("time (#mus)");
  cur_wave->GetYaxis()->SetTitle("Current (ADC channels)");
  cur_wave->GetYaxis()->SetTitleOffset(2.0);

  gStyle->SetOptStat(0);
  
  cur_wave->Draw("colz");

  // draw lines representing currents for two horizontal bands

  // lower current
  TLine* l_curr_l = new TLine(-period*1e6, mean1,period*1e6, mean1);
  l_curr_l->SetLineColor(kRed);
  l_curr_l->Draw("same");

  TLine* l_curr_h = new TLine(-period*1e6, mean2,period*1e6, mean2);
  l_curr_h->SetLineColor(kBlue);
  l_curr_h->Draw("same");


  TLine* l_curr_m = new TLine(-period*1e6, mean2-curr_corr,period*1e6, mean2-curr_corr);
  l_curr_m->SetLineColor(kGreen+2);
  l_curr_m->Draw("same");


  // time of lower current peak

  TLine* l_t_l  = new TLine(1e6*time_1, base-amplitude,1e6*time_1, base+amplitude);

  l_t_l->SetLineColor(kRed);
  l_t_l->SetLineStyle(9);
  l_t_l->Draw("same");

  // corrected time of lower current peak

  TLine* l_t_l_corr  = new TLine(1e6*time_1_corr, base-amplitude,1e6*time_1_corr, base+amplitude);

  l_t_l_corr->SetLineColor(kGreen+2);
  l_t_l_corr->SetLineStyle(9);
  l_t_l_corr->Draw("same");


  // time of higher current peak

  TLine* l_t_h  = new TLine(1e6*time_2, base-amplitude,1e6*time_2, base+amplitude);
  
  l_t_h->SetLineColor(kBlue);
  l_t_h->SetLineStyle(9);
  l_t_h->Draw("same");

  // corrected time of higher current peak

  TLine* l_t_h_corr  = new TLine(1e6*time_2_corr, base-amplitude,1e6*time_2_corr, base+amplitude);

  l_t_h_corr->SetLineColor(kGreen+2);
  l_t_h_corr->SetLineStyle(9);
  l_t_h_corr->Draw("same");


  TLegend* leg_wave = new TLegend(.20,.55,.55,.9,"Key");

  leg_wave->SetFillColor(0);
  //  leg_wave->SetFillStyle(0);
   leg_wave->AddEntry(l_curr_l,"Lower current peak","l");
   leg_wave->AddEntry(l_curr_h,"Higher current peak","l");
   leg_wave->AddEntry(l_curr_m,"Average current of peaks","l");
   leg_wave->AddEntry(l_t_l,"Time of lower peak","l");
   leg_wave->AddEntry(l_t_l_corr,"Adjusted Time of lower peak","l");
   leg_wave->AddEntry(l_t_h,"Time of higher peak","l");
   leg_wave->AddEntry(l_t_h_corr,"Adjusted Time of higher peak","l");
   leg_wave->Draw("same");


   TPaveText *tt_diff = new TPaveText(0.55,0.15,0.88,0.25,"NDC");
   tt_diff->SetFillColor(0);
   tt_diff->SetShadowColor(0);
   tt_diff->AddText(Form("time correction = %.3f #mus",time_1_diff*1e6));
   tt_diff->Draw("same");
   



   // loop over all entries and plot adjusted current (and adjusted current vs BPM Y)

   Double_t Lrb_rawcur_y = 0;

   Double_t Lrb_time = 0;
   Double_t L_BPMB_y = 0;

   // adj for adjusted in following variables
   
   Double_t Lrb_time_adj = 0;
   Double_t Lrb_rawcur_y_adj = 0;
  

   TH2F* h_rcurr_adj = new TH2F("h_rcurr_adj","raster current vs BPMB Y",1000,BPM_lim_l,BPM_lim_h,1000,Lrby_curr_l,Lrby_curr_h);

   h_rcurr_adj->GetXaxis()->SetTitle("BPMB Y [m]");
   h_rcurr_adj->GetYaxis()->SetTitle("Raster Y current [ADC channels]");
  h_rcurr_adj->GetYaxis()->SetTitleOffset(2.0);
  

   TH1F* h_ycurr_adj = new TH1F("h_ycurr_adj","Adjusted raster current Y",1000,Lrby_curr_l,Lrby_curr_h);


   TH2F* h_comp_curr = new TH2F("h_comp_curr","raster current vs adjusted raster current",1000,Lrby_curr_l,Lrby_curr_h,1000,Lrby_curr_l,Lrby_curr_h);
   h_comp_curr->GetXaxis()->SetTitle("Adjusted Y Current (ADC channels)");
   h_comp_curr->GetYaxis()->SetTitle("Recorded Y Current (ADC channels)");
   h_comp_curr->GetYaxis()->SetTitleOffset(2.0);


   TH2F* h_comp_time = new TH2F("h_comp_time","time vs adjusted time",1000,-40E-6,15E-6,1000,-40E-6,15E-6);
   h_comp_time->GetXaxis()->SetTitle("Adjusted time");
   h_comp_time->GetYaxis()->SetTitle("time");

   Int_t NEntries = T->GetEntries();

   // parameters for phase cut


   Double_t p_slope = (l_y2-l_y1)/(l_x2-l_x1);

   Double_t p_off = l_y1 - p_slope*l_x1;


   cout << "p_slope = " << p_slope << endl;
   cout << "p_off = " << p_off << endl;
   
   // Double_t p_slope = 1.103e7;
   // Double_t p_off = 1.25e3; // (cut intercept with y-axis)
   
   // Double_t p_slope = 8.698e6;
   // Double_t p_off = 7.387e3; // (cut intercept with y-axis)

   
   
   Bool_t phase; 

   T->SetBranchAddress("Lrb.Raster2.rawcur.y",&Lrb_rawcur_y);
   T->SetBranchAddress("Lrb.BPMB.y",&L_BPMB_y);

   std::vector<Double_t> v_time;
   std::vector<Double_t> v_curr;

   std::vector<Double_t> v_time_adj;
   std::vector<Double_t> v_curr_adj;
   
   
   for(Int_t i = 0; i<NEntries; i++){


     T->GetEntry(i);
     
     
     // need to check phase of raster
     
     if( (Lrb_rawcur_y -  p_slope*L_BPMB_y) > p_off ){
       // positive slope
       phase = true;

     }
     else{
       // negative slope/ phase
       phase = false;
     }
     
     Lrb_time =  curr_to_time(Lrb_rawcur_y,wave_pars,phase);
     //     cout << "Lrb_time = " << Lrb_time << endl;
     Lrb_time_adj =  Lrb_time - time_2_diff;
     



     
     //     cout << "Lrb_time_adj = " << Lrb_time_adj << endl;
     Lrb_rawcur_y_adj = time_to_curr(Lrb_time_adj,wave_pars);
     //     cout << "Lrb_rawcur_y_adj = " << Lrb_rawcur_y_adj << endl;

     h_rcurr_adj->Fill(L_BPMB_y,Lrb_rawcur_y_adj);
     h_ycurr_adj->Fill(Lrb_rawcur_y_adj);


     h_comp_curr->Fill(Lrb_rawcur_y_adj,Lrb_rawcur_y);
     h_comp_time->Fill(Lrb_time,Lrb_time_adj);
     
     v_time.push_back(Lrb_time/period);
     v_curr.push_back(Lrb_rawcur_y);

     v_time_adj.push_back(Lrb_time_adj/period);
     v_curr_adj.push_back(Lrb_rawcur_y_adj);
     
   }


   TCanvas *c4 = new TCanvas("c4","Current adjusted raster vs BPM",1200,1200);

   c4->SetLeftMargin(0.2);
   c4->Divide(2,1);
   c4->cd(1);
   h_rcurr->Draw("col");
   phase_l->Draw("same");
   c4->cd(2);
   h_rcurr_adj->Draw("col");
   phase_l->Draw("same");
     
   
   TCanvas *c5 = new TCanvas("c5","Current adjusted raster",1200,1200);
   c5->SetLeftMargin(0.2);
   c5->Divide(2,1);
   c5->cd(1);
   h_ycurr->Draw("colz");
   c5->cd(2);
   h_ycurr_adj->GetXaxis()->SetTitle("Raster Y current [ADC channels]");
   h_ycurr_adj->Draw("colz");

   TGraph* g_curr = new TGraph(NEntries,&v_time[0],&v_curr[0]);
   TGraph* g_curr_adj = new TGraph(NEntries,&v_time_adj[0],&v_curr_adj[0]);
   
   TCanvas *c6 = new TCanvas("c6","Current vs Time",1200,1200);
   c6->Divide(2,1);
   c6->cd(1);
   g_curr->Draw("A*");
   c6->cd(2);
   g_curr_adj->Draw("A*");


   TCanvas *c7 = new TCanvas("c7","Current vs Time",1200,1200);
   c7->SetLeftMargin(0.2);
   // c7->Divide(2,1);
   // c7->cd(1);
   h_comp_curr->Draw("colz");
   // c7->cd(2);
   // h_comp_time->Draw("colz");



   // save canvases as pdfs
   c3->Print("horizontal_foils/corrections/raster_shape.pdf");
   c4->Print("horizontal_foils/corrections/rast_BPM.pdf");
   c5->Print("horizontal_foils/corrections/rast_y.pdf");
   //   c3->Print("horizontal_foils/corrections/raster_shape.pdf");
   c7->Print("horizontal_foils/corrections/curr_vs_curr_adj.pdf");
   
}
