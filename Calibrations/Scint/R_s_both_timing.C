/*
*************************************************************
10/12/20 John Williamson
Script that performs timing offset calibration for s0 and s2 for the RHRS.
This uses technique described in Longwu Ou thesis (GMP experiment): https://dspace.mit.edu/handle/1721.1/
and Tong Su thesis (tritium): https://etd.ohiolink.edu/apexprod/rws_olink/r/1501/10?clear=10&p10_accession_num=kent1587680491082341.

Uses relation between difference of left and right paddles for S2 (top and bottom equivalent for S0) to get information about difference in offsets between left and right paddles. 

Can use TOF (Time-Of-Flight) and difference between sums of left and right in S2 top and bottom in S0 to establish relationship.

As only difference between S2 and S0 (TOF) matters, system has a degree of freedom (could increase offsets for top, bottom in S0 and left, right in S2 and not affect TOF). This means one of the S0 offsets can be fixed, and from the knwoledge of the difference the other S0 difference can be calculated. As described above we now have knowledge of the sum and difference of the S2 coefficients (this must be done for each paddle) and so can calculate both offsets for S2.

*************************************************************
 */

#include "file_def.h"
#include "Load_more_rootfiles.C"

void R_s_both_timing(Int_t runno){


  TChain* t = Load_more_rootfiles(runno);

  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);




  // constants for S0 and S2

  const Float_t ctos = 0.5E-9; //0.5ns per TDC channel
  const Float_t ry_s2 = 0.432; //Full length of s2 along y [m]
  const Float_t rx_s0 = 1.70; //Full length of s0 along x [m]
  const Float_t c_s2 = 1.26193e+08; //Speed of light in s2 [m/s]
  const Float_t c_s0 = 1.23858e+08; //Speed of light in s0 [m/s]
  const Float_t elec_v = 3.0E8; //Speed of electron between s0 and s2 [m/s]

  
  // alternative names for constants
  Double_t L=1.7;
  Double_t l=0.432;
  Double_t c0=1.21577e+08;

 

  
  Float_t tof; //in seconds
  Float_t tof_sum; //in channels
  Float_t  TOF;//in Channel

  
  const Int_t NS2Pad = 16;

  
  // define track cuts (from projected tracks from VDC)
  
  Double_t R_x_l[16] = {-1.1,-0.99,-0.88,-0.70,-0.55,-0.40,-0.3,-0.15,-0.01,0.12,0.27,0.40,0.53,0.68,0.82};
  Double_t R_x_h[16] = {-0.95,-0.83,-0.68,-0.55,-0.40,-0.24,-0.12,0.0,0.17,0.29,0.40,0.55,0.69,0.82,0.99,1.04};
  
  
  //Initialize S0 Histograms
  TH2F *h1a = new TH2F("h1a","RHRS Left_s0_TDC vs X: Good Events",500,-1.25,1.25,400,2200,2600);
  h1a->GetXaxis()->SetTitle("X_s0 [m]"); h1a->GetXaxis()->CenterTitle();
  h1a->GetYaxis()->SetTitle("TDC Channel"); h1a->GetYaxis()->CenterTitle();
  
  TH2F *h1b = new TH2F("h1b","RHRS Right_s0_TDC vs Track  X: Good Events",500,-1.25,1.25,400,2200,2600);
  h1b->GetXaxis()->SetTitle("X_s0 [m]"); h1b->GetXaxis()->CenterTitle();
  h1b->GetYaxis()->SetTitle("TDC Channel"); h1b->GetYaxis()->CenterTitle();

  TH2F *h2a = new TH2F("h2a","RHRS S0 : Left_TDC - Right_TDC  vs  X: Good Events",500,-1.25,1.25,180,-60,120);
  h2a->GetXaxis()->SetTitle("X_s0 [m]"); h2a->GetXaxis()->CenterTitle();
  h2a->GetYaxis()->SetTitle("Left_TDC-Right_TDC"); h2a->GetYaxis()->CenterTitle();
  
  TProfile *h2b = new TProfile("h2b","RHRS S0 :Left - Right  vs X: Good Events",100,-1.5,1.5,-60,120);
  h2b->SetLineColor(kRed);h2b->SetMarkerColor(kRed);



  //Initialise S2 histograms


  TH2 *hleft[16];
  TH2 *hright[16];
  TH2 *hdiff[16];
  TProfile *pdiff[16];
  TH1 *hsum[16];
    
  Int_t s2l_low=2300;
  Int_t s2l_high=2700;
  Int_t s2r_low=2300;
  Int_t s2r_high=2700;
  Int_t s2l_bin=(s2l_high-s2l_low);
  Int_t s2r_bin=(s2r_high-s2r_low);
  Int_t s2diff_low=s2l_low-s2r_high;//s2 diff: Left-Right
  Int_t s2diff_high=s2l_high-s2r_low;
  Int_t s2diff_bin=s2diff_high-s2diff_low;
  Int_t s2sum_high=s2l_high+s2r_high;
  Int_t s2sum_low=s2l_low+s2r_high;
  Int_t s2sum_bin=s2sum_high-s2sum_low;

  for(Int_t i=0;i<16;i++){
    
    hleft[i] = new TH2F(Form("hleft[%d]",i),Form("S2 Left PMT TDC vs Y_s2: Paddle %d",i+1),300,-0.3,0.3,s2l_bin,s2l_low,s2l_high);
    hleft[i]->GetXaxis()->SetTitle("Y_s2 [m]"); hleft[i]->GetXaxis()->CenterTitle();
    hleft[i]->GetYaxis()->SetTitle("TDC Channel"); hleft[i]->GetYaxis()->CenterTitle();
    
    hright[i] = new TH2F(Form("hright[%d]",i),Form("S2 Right PMT TDC vs. Y_s2: Paddle %d",i+1),300,-0.3,0.3,s2r_bin,s2r_low,s2r_high);
    hright[i]->GetXaxis()->SetTitle("Y_s2 [m]"); hright[i]->GetXaxis()->CenterTitle();
    hright[i]->GetYaxis()->SetTitle("TDC Channel"); hright[i]->GetYaxis()->CenterTitle();

    hdiff[i] = new TH2F(Form("hdiff[%d]",i),Form("S2 Left - Right vs Y_s2: Paddle %d",i+1),300,-0.3,0.3,s2diff_bin,s2diff_low,s2diff_high);
    hdiff[i]->GetXaxis()->SetTitle("Y_s2 [m]"); hdiff[i]->GetXaxis()->CenterTitle();
    hdiff[i]->GetYaxis()->SetTitle("TDC Channel"); hdiff[i]->GetYaxis()->CenterTitle();

    if(i==0) pdiff[i] = new TProfile(Form("pdiff[%d]",i),Form("S2 Left - Right vs Y_s2: Paddle %d",i+1),6,-0.3,0.3,s2diff_low,s2diff_high);
    else pdiff[i] = new TProfile(Form("pdiff[%d]",i),Form("S2 Left - Right vs. Track Projection Y: Paddle %d",i+1),90,-0.3,0.3,s2diff_low,s2diff_high);
    pdiff[i]->GetXaxis()->SetTitle("Y_s2 [m]"); pdiff[i]->GetXaxis()->CenterTitle();
    pdiff[i]->GetYaxis()->SetTitle("TDC Channel"); pdiff[i]->GetYaxis()->CenterTitle();
    pdiff[i]->SetLineColor(kBlack);pdiff[i]->SetMarkerColor(kBlack);

    //hsum[i] = new TH2F(Form("hsum[%d]",i),Form("S2 Left + Right vs. Track Projection Y: Paddle %d",i+1),300,-0.3,0.3,10000,0,10000);
    hsum[i] = new TH1F(Form("hsum[%d]",i),Form("S2 Left + Right: Paddle %d",i+1),200,2800,3000);
    //    hsum[i]->GetXaxis()->SetTitle("Y_s2 [m]"); hsum[i]->GetXaxis()->CenterTitle();
    hsum[i]->GetXaxis()->SetTitle("TDC Channel (sum)"); hsum[i]->GetYaxis()->CenterTitle();

  }




  
    
  //Define Variables
  Double_t R_tr_n,R_cer_asum_c,R_ps_e,R_sh_e;
  Double_t R_tr_p[100],R_s0_trx[100],R_s2_try[100],R_s2_trx[100];
  Double_t R_s0_lt[10],R_s0_rt[10];
  Double_t R_s0_nthit;
  Double_t R_s2_lt[16],R_s2_rt[16];
  Double_t R_s2_nthit;
  Double_t R_s2_t_pads[16];
  Double_t evtypebits;
  Double_t R_s0_trpath[100],R_s2_trpath[100];

  
  //Define Branch Status/Addresses
  t->SetBranchStatus("*",0);
  t->SetBranchStatus("R.tr.n",1);
  t->SetBranchStatus("R.tr.p",1);
  t->SetBranchStatus("R.cer.asum_c",1);
  t->SetBranchStatus("R.ps.e",1);
  t->SetBranchStatus("R.sh.e",1);
  t->SetBranchStatus("R.s0.lt",1);
  t->SetBranchStatus("R.s0.rt",1);
  t->SetBranchStatus("R.s0.trx",1);
  t->SetBranchStatus("R.s0.nthit",1);
  t->SetBranchStatus("R.s2.try",1);
  t->SetBranchStatus("R.s2.trx",1);
  t->SetBranchStatus("R.s2.nthit",1);
  t->SetBranchStatus("R.s2.t_pads",1);
  t->SetBranchStatus("DR.evtypebits",1);
  t->SetBranchStatus("R.s0.trpath",1);
  t->SetBranchStatus("R.s2.trpath",1);

  t->SetBranchAddress("R.tr.n",&R_tr_n);
  t->SetBranchAddress("R.tr.p",R_tr_p);
  t->SetBranchAddress("R.cer.asum_c",&R_cer_asum_c);
  t->SetBranchAddress("R.ps.e",&R_ps_e);
  t->SetBranchAddress("R.sh.e",&R_sh_e);
  t->SetBranchAddress("R.s0.lt",R_s0_lt);
  t->SetBranchAddress("R.s0.rt",R_s0_rt);
  t->SetBranchAddress("R.s0.trx",R_s0_trx);
  t->SetBranchAddress("R.s0.nthit",&R_s0_nthit);
  t->SetBranchAddress("R.s2.lt",R_s2_lt);
  t->SetBranchAddress("R.s2.rt",R_s2_rt);
  t->SetBranchAddress("R.s2.try",R_s2_try);
  t->SetBranchAddress("R.s2.trx",R_s2_trx);
  t->SetBranchAddress("R.s2.nthit",&R_s2_nthit);
  t->SetBranchAddress("R.s2.t_pads",R_s2_t_pads);
  t->SetBranchAddress("DR.evtypebits",&evtypebits);
  t->SetBranchAddress("R.s0.trpath",R_s0_trpath);
  t->SetBranchAddress("R.s2.trpath",R_s2_trpath);

  Int_t nentries = t->GetEntries();

  cout<<"Total Number of Events = "<<nentries<<endl;


  // Use first Loop over events to get offsets for S0 and
  // to extract speed of light in both scitntillators
  // needed to extracting S2 offsets later on


  
  for(Int_t i=0;i<nentries;i++){
    
    if(i%100000==0) cout << " events processed = " << i << endl;
    t->GetEntry(i);

   
    //    if(R_tr_n==1 && R_cer_asum_c>1500 && (R_ps_e+R_sh_e)/(1000.*R_tr_p[0])>0.8 && R_s0_nthit==1 && ((Int_t)evtypebits>>1&1)){

    // S0
    if(R_tr_n==1 && R_cer_asum_c>1500 && (R_ps_e+R_sh_e)/(1000.*R_tr_p[0])>0.8 && R_s0_nthit==1){
     
    
      h1a->Fill(R_s0_trx[0],R_s0_lt[0]);
      h1b->Fill(R_s0_trx[0],R_s0_rt[0]);
      h2a->Fill(R_s0_trx[0],R_s0_lt[0]-R_s0_rt[0]);
      h2b->Fill(R_s0_trx[0],R_s0_lt[0]-R_s0_rt[0]);
      
    }


    //S2
    Int_t x=evtypebits;
    for(Int_t j=0;j<16;j++){
      if(R_tr_n==1 && R_cer_asum_c>2500&&R_s2_t_pads[0]==j&&R_s0_nthit==1&&(R_ps_e+R_sh_e)/R_tr_p[0]/1000>0.7 && R_s2_trx[0]>R_x_l[j] && R_s2_trx[0]<R_x_h[j]){

	tof=(R_s2_trpath[0]-R_s0_trpath[0])/elec_v;
        TOF=tof/ctos;
	//        Double_t SUM=(R_s0_lt[0]+R_s0_rt[0])-(R_s2_rt[j]+R_s2_lt[j])+(L/c0/ctos)-(l/c_s2/ctos)+3018.26-2*TOF;
        
        hleft[j]->Fill(R_s2_try[0],R_s2_lt[j]);
        hright[j]->Fill(R_s2_try[0],R_s2_rt[j]);
        hdiff[j]->Fill(R_s2_try[0],R_s2_lt[j]-R_s2_rt[j]);
        pdiff[j]->Fill(R_s2_try[0],R_s2_lt[j]-R_s2_rt[j]);
        //hsum[j]->Fill(SUM);

	
      }
    }
    
  }




  
  //Make Plots for S0
  TCanvas *c1 = new TCanvas("c1");
  c1->Divide(1,2);
  c1->cd(1);
  h1a->Draw("col");
  c1->cd(2);
  h1b->Draw("col");
  
  TCanvas *c2 = new TCanvas("c2");
  h2b->Draw();
  h2b->Fit("pol1","R","",-0.75,0.75);
  c2->Update();
  h2a->Draw("col");
  h2b->Draw("same");
  c2->Modified();

  Float_t difference =  h2b->GetFunction("pol1")->GetParameter(0);
  Float_t temp_speed = h2b->GetFunction("pol1")->GetParameter(1);
  Float_t c_speed = abs((2./temp_speed)/ctos);
  
 

  c1->Print(Form("plots/R_s0_timing_%i.pdf[",runno));
  c1->Print(Form("plots/R_s0_timing_%i.pdf",runno));
  c2->Print(Form("plots/R_s0_timing_%i.pdf",runno));
  c2->Print(Form("plots/R_s0_timing_%i.pdf]",runno));


  // t_stop value is used to convert from 'real time' to TDC common stop time and vice versa
  // 
  Double_t t_stop = 3000.0; 
  
  Float_t corr_left = t_stop - 1500.00;
  Float_t corr_right = t_stop - (1500.00+difference);
  

  

  //Make Plots for S2
  TCanvas *c3[4],*c4[7],*c5[7];
  Int_t scin_count,pad_count;

  Int_t counter(16);
  Float_t pad[16],intercept[16],s2_speed[16],slope[16],slope_err[16];
  Float_t sum[16];

  
  

  for(Int_t i=0;i<4;i++){
    c3[i]= new TCanvas(Form("c3[%d]",i));
    c3[i]->Divide(2,4);
    pad_count=0;
    
    for(Int_t j=0;j<4;j++){
      scin_count = 4*i +j;
      pad_count++; c3[i]->cd(pad_count);
      hleft[scin_count]->Draw("col");
      
      pad_count++; c3[i]->cd(pad_count);
      hright[scin_count]->Draw("col");
    }
    
    if(i==0){
     c3[i]->Print(Form("plots/R_s2_timing_%i.pdf[",runno));
     c3[i]->Print(Form("plots/R_s2_timing_%i.pdf",runno));
    }
    else{
      c3[i]->Print(Form("plots/R_s2_timing_%i.pdf",runno));
    }
  }


  
  for(Int_t i=0;i<8;i++){
    c4[i]= new TCanvas(Form("c4[%d]",i));
    c4[i]->Divide(1,2);
    pad_count=0;
    for(Int_t j=0;j<2;j++){
      scin_count = 2*i + j;
      pad_count++; c4[i]->cd(pad_count);
    
        pdiff[scin_count]->Draw();
        pdiff[scin_count]->Fit("pol1","R","",-0.125,0.125);
        c4[i]->Update();
        
	intercept[scin_count] = pdiff[scin_count]->GetFunction("pol1")->GetParameter(0);
	temp_speed = pdiff[scin_count]->GetFunction("pol1")->GetParameter(1);
	s2_speed[scin_count] = abs((2./temp_speed)/ctos);
	 
      c4[i]->Modified();
    }
      c4[i]->Print(Form("plots/R_s2_timing_%i.pdf",runno));
  }




  
  // This loop calculates sum of S2 coeffecients

  for(Int_t i=0;i<nentries;i++){
    
    if(i%100000==0) cout << " events processed = " << i << endl;
    t->GetEntry(i);
   
    
    Int_t x=evtypebits;
    for(Int_t j=0;j<16;j++){
      if(R_tr_n==1 && R_cer_asum_c>2500&&R_s2_t_pads[0]==j&&R_s0_nthit==1&&(R_ps_e+R_sh_e)/R_tr_p[0]/1000>0.7){

	tof=(R_s2_trpath[0]-R_s0_trpath[0])/elec_v;
	TOF=tof/ctos;
	//	Double_t SUM=(R_s0_lt[0]+R_s0_rt[0])-(R_s2_rt[j]+R_s2_lt[j])+(L/c_speed/ctos)-(l/s2_speed/ctos)+3018.26-2*TOF;
	Double_t SUM=(R_s0_lt[0]+R_s0_rt[0])-(R_s2_rt[j]+R_s2_lt[j])+(L/c_speed/ctos)-(l/s2_speed[j]/ctos)+ (corr_left+corr_right) - 2*TOF;

	
	hsum[j]->Fill(SUM);

	
      }
    }
  }

  

  for(Int_t i=0;i<8;i++){
    c5[i]= new TCanvas(Form("c5[%d]",i));
    c5[i]->Divide(1,2);
    pad_count=0;
    for(Int_t j=0;j<2;j++){
      scin_count = 2*i + j;
      pad_count++;
      c5[i]->cd(pad_count);
      hsum[scin_count]->Draw();
      hsum[scin_count]->Fit("gaus");
      double xx=hsum[scin_count]->GetFunction("gaus")->GetParameter(2);
      double yy=hsum[scin_count]->GetFunction("gaus")->GetParameter(1);
      hsum[scin_count]->Fit("gaus","R","",yy-xx,yy+xx);
      c5[i]->Modified();
      sum[scin_count]=hsum[scin_count]->GetFunction("gaus")->GetParameter(1);
      cout<<"mean value="<<sum[scin_count]<<" haha"<<endl;
    }
    if(i==7){
      c5[i]->Print(Form("plots/R_s2_timing_%i.pdf",runno));
      c5[i]->Print(Form("plots/R_s2_timing_%i.pdf]",runno));
    }
    else c5[i]->Print(Form("plots/R_s2_timing_%i.pdf",runno));
  }

   cout<<"------------------------------------------------------\n";
    
   
   Double_t left[16],right[16];


   ofstream oofile(Form("DB/RHRS_TOF_timing_%i.txt",runno));
   oofile<<"RHRS Scintillator timing offsets: "<<endl<<endl;
   oofile<<"S2 Offsets: "<<endl;
   oofile<<"R.s2.L.off =  ";
    
   for (Int_t ii=0;ii<16;ii++)
     {
       //       left[ii]=(sum[ii]-intercept[ii])/2.0;
       left[ii]=t_stop- ((sum[ii]-intercept[ii])/2.0);
       cout<<"s2_left["<<ii+1<<"]="<<left[ii]<<";"<<endl;
       
       oofile<<left[ii]<<" ";
     }

   oofile<<endl;
   oofile<<"R.s2.R.off =  ";
   for (Int_t ii=0;ii<16;ii++)
     {
       //       right[ii]=(sum[ii]+intercept[ii])/2.0;
       right[ii]=t_stop - ((sum[ii]+intercept[ii])/2.0);       
       cout<<"s2_right["<<ii+1<<"]="<<right[ii]<<";"<<endl;
       
       oofile<<right[ii]<<" ";       
     }
   
   oofile <<endl<<endl;



    
    
  //Write some info to the screen
  cout<<"\n------------------------------------------------------\n";

  cout<<"Left - Right Intercept [channel]: "<<difference<<endl;
  
  cout<<"Left PMT Correction Factor: "<<corr_left<<endl;
  cout<<"Right PMT Correction Factor: "<<corr_right<<endl<<endl;

  oofile<<endl;
  oofile<<"S0_left = "<<corr_left<<endl;
  oofile<<"S0_Right = "<<corr_right<<endl;
  oofile.close();

  
  cout<<"Speed of Light in S0 = "<<c_speed<<" m/s !!!"<<endl;
  cout<<"Speed of Light in S0 = "<<temp_speed<<" m/chan (m/0.5ns) !!!"<<endl;
  cout<<"------------------------------------------------------\n";




 
  cout<<"Speed of Light in S2 = "<<s2_speed[3]<<" m/s !!!"<<endl;
  cout<<"------------------------------------------------------\n";


  // // calc means of left and right offsets

  Double_t left_off_mean = 0.0;
  Double_t right_off_mean = 0.0;

  for(Int_t i = 0; i < NS2Pad; i ++){

    left_off_mean += left[i];
    right_off_mean += right[i];

  }
  
  left_off_mean = left_off_mean/NS2Pad;
  right_off_mean = right_off_mean/NS2Pad;
  Double_t lr_off_mean = (left_off_mean+right_off_mean)/2.;


  cout << "left_off_mean = " << left_off_mean << endl;
  cout << "right_off_mean = " << right_off_mean << endl;
  cout << "lr_off_mean = " << lr_off_mean << endl;

  


  // save corrections to csv DB file

  ofstream oofile_csv(Form("DB/RHRS_TOF_timing_%i.csv",runno));

  oofile_csv<<"R.s2.L.off, R.s2.R.off" << endl;
    
   for (Int_t ii=0;ii<16;ii++)
     {
       oofile_csv<<left[ii]<<","<<right[ii]<<endl;       
     }

   oofile_csv.close();



   // plot results of corrections

   TH2F* hsum_all = new TH2F("hsum_all","Corrected S2 Left + Right",16,0,15,200,2490-lr_off_mean,2520-lr_off_mean);

  hsum_all->GetXaxis()->SetTitle("TDC Channel (sum)"); hsum_all->GetXaxis()->CenterTitle();

  TH2F* hsum_all_l = new TH2F("hsum_all_l","corrected S2 Left",16,0,15,30,2490-left_off_mean,2520-left_off_mean);

  hsum_all_l->GetXaxis()->SetTitle("TDC Channel (sum)"); hsum_all_l->GetXaxis()->CenterTitle();

  TH2F* hsum_all_r = new TH2F("hsum_all_r","corrected S2 right",16,0,15,25,2490-right_off_mean,2520-right_off_mean);

  hsum_all_r->GetXaxis()->SetTitle("TDC Channel (sum)"); hsum_all_r->GetXaxis()->CenterTitle();

  

  // uncorrected equivalent of plots
  
  TH2F* hsum_all_un = new TH2F("hsum_all_un","Uncorrected S2 Left + Right",16,0,15,200,2490,2520);

  hsum_all_un->GetXaxis()->SetTitle("TDC Channel (sum)"); hsum_all_un->GetXaxis()->CenterTitle();

  TH2F* hsum_all_un_l = new TH2F("hsum_all_un_l","Uncorrected S2 Left",16,0,15,30,2490,2520);

  hsum_all_un_l->GetXaxis()->SetTitle("TDC Channel (sum)"); hsum_all_un_l->GetXaxis()->CenterTitle();

  TH2F* hsum_all_un_r = new TH2F("hsum_all_un_r","Uncorrected S2 right",16,0,15,25,2490,2520);

  hsum_all_un_r->GetXaxis()->SetTitle("TDC Channel (sum)"); hsum_all_un_r->GetXaxis()->CenterTitle();


  Double_t lrtime_mean = 0.0;
  Double_t ltime_mean = 0.0;
  Double_t rtime_mean = 0.0;

  Int_t mean_count = 0;
  
  for(Int_t i=0;i<nentries;i++){
    
    if(i%100000==0) cout << " events processed = " << i << endl;
    t->GetEntry(i);
    //    cout << " events processed = " << i << endl;
    
    Int_t x=evtypebits;
    for(Int_t j=0;j<16;j++){
      if(R_tr_n==1 && R_cer_asum_c>1500 && (R_ps_e+R_sh_e)/(1000.*R_tr_p[0])>0.8 && R_s0_nthit==1){

	tof=(R_s2_trpath[0]-R_s0_trpath[0])/elec_v;
        TOF=tof/ctos;

	Double_t left_time = R_s2_lt[j];
	Double_t right_time = R_s2_rt[j];
	Double_t lr_time = (left_time+right_time)/2.;

	Double_t left_time_corr = R_s2_lt[j] - left[j];
	Double_t right_time_corr = R_s2_rt[j] - right[j];
	Double_t lr_time_corr = (left_time_corr+right_time_corr)/2.;

	
	hsum_all->Fill(j,lr_time_corr);
	hsum_all_l->Fill(j,left_time_corr);
	hsum_all_r->Fill(j,right_time_corr);


	lrtime_mean += lr_time_corr;
	ltime_mean += left_time_corr;
	rtime_mean += right_time_corr;
	
	hsum_all_un->Fill(j,lr_time);
	hsum_all_un_l->Fill(j,left_time);
	hsum_all_un_r->Fill(j,right_time);

	mean_count++;
      }
    }
  }



  TCanvas* c6 = new TCanvas("c6","c6");
  c6->Divide(3,2);

  c6->cd(1);
  hsum_all->Draw("colz");

  c6->cd(2);
  hsum_all_l->Draw("colz");

  c6->cd(3);
  hsum_all_r->Draw("colz");
  
  c6->cd(4);
  hsum_all_un->Draw("colz");

  c6->cd(5);
  hsum_all_un_l->Draw("colz");

  c6->cd(6);
  hsum_all_un_r->Draw("colz");

  
  lrtime_mean = lrtime_mean/mean_count;
  ltime_mean = ltime_mean/mean_count;
  rtime_mean = rtime_mean/mean_count;

  cout << "lrtime_mean = " << lrtime_mean << endl;
  cout << "ltime_mean = " << ltime_mean << endl;
  cout << "rtime_mean = " << rtime_mean << endl;

   
  
}

