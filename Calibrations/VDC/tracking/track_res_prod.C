/**************************************************************
  track_res_porod.C
  John Williamson    
  29th April, 2021

  Calculates sigma for timing offsets for clusters and resolution of local distance (for production runs)

  - timing offset: plot timing offset and get width of central gaussian (no contribution from accidentals)

  - distance: plot difference from global track and local track and take intercept

  
*************************************************************/


#include "Load_more_rootfiles.C"
#include "file_def.h"
#include "TTD_namespace.h"

#define NPLANE 4
#define MAX_ENTRIES 100000
#define MAX_HIT 1000

#define DEG_TO_RAD 0.017453278


void track_res_prod(const char *arm = "L", Int_t runno = 4668, const char* title = ""){


  const char plane[NPLANE][8] = {"u1", "u2", "v1", "v2"};

  const Int_t NAngle = 2;
  const Int_t NPROJ = 4;
  const char proj[NPROJ][8] = {"X12","X21","Y12","Y21"};
  const char projX[NPROJ][8] = {"X","X","Y","Y"};
  const char projFull[NPROJ][20] = {"X_{1} - P_{x,2} (m)","X_{2} - P_{x,1} (m)","Y_{1} - P_{y,2} (m)","Y_{2} - P_{y,1} (m)"};
  

  TString projROOTX[NPROJ] = {Form("%s.tr.UV12_X[0]",arm),Form("%s.tr.UV21_X[0]",arm),Form("%s.tr.UV12_Y[0]",arm),Form("%s.tr.UV21_Y[0]",arm)};
  TString projROOTPX[NPROJ] = {Form("%s.tr.UV12_PX[0]",arm),Form("%s.tr.UV21_PX[0]",arm),Form("%s.tr.UV12_PY[0]",arm),Form("%s.tr.UV21_PY[0]",arm)};

  

  Double_t XLim[NPROJ] = {1.5,1.5,0.06,0.06};
  Double_t DiffLim[NPROJ] = {0.05,0.05,0.03,0.03};
  
  // 1D Projection difference plots

  TH1F* hDistDiff[NPROJ];
  TH1F* hDistDiffCorr[NPROJ];

  TH1F* hDistDiffCorrNorm[NPROJ];
  
  
  // 2D projection vs plots

  TH2F* hDistDiffX[NPROJ];
  TH2F* hDistDiffXCorr[NPROJ];

  TH2F* hDistDiffPX[NPROJ];
  TH2F* hDistDiffPXCorr[NPROJ];


  // 1D angular difference plots

  TH1F* hXAngleDiff = new TH1F("hXAngleDiff","hXAngleDiff",100,-0.10,0.10);
  TH1F* hXAngleDiffTest = new TH1F("hXAngleDiffTest","hXAngleDiffTest",100,-0.10,0.10);
  TH1F* hYAngleDiff = new TH1F("hYAngleDiff","hYAngleDiff",100,-0.10,0.10);
  TH1F* hYAngleDiffTest = new TH1F("hYAngleDiffTest","hYAngleDiffTest",100,-0.10,0.10);


  for(Int_t i = 0; i<NPROJ; i++){

    hDistDiff[i] = new TH1F(Form("hDistDiff_%s",proj[i]),Form("%sHRS %s",arm,proj[i]),100, -0.05, 0.05);
    hDistDiff[i]->GetXaxis()->SetTitle( Form("%s",projFull[i]));

    hDistDiffCorr[i] = (TH1F*) hDistDiff[i]->Clone(Form("hDistDiffCorr_%s",proj[i]));
    hDistDiffCorr[i]->SetTitle(Form("Corrected Distance Difference %s",proj[i]));
    
    hDistDiffCorrNorm[i] = new TH1F(Form("hDistDiffCorrNorm_%s",proj[i]),Form("%sHRS %s",arm,proj[i]),100, -5, 5);
    hDistDiffCorrNorm[i]->GetXaxis()->SetTitle(Form("(%s%i - P%s%i)/#sigma_{%s%i}",projX[i],(i%2)+1,projX[i],(i+1)%2+1,projX[i],i%2+1));
        

    hDistDiffX[i] = new TH2F(Form("hDistDiffX_%s",proj[i]),Form("%sHRS %s vs %s",arm,proj[i],projX[i]),100,-XLim[i],XLim[i],100, -DiffLim[i], DiffLim[i]);   
    hDistDiffX[i]->GetXaxis()->SetTitle( Form("%s_{%i}",projX[i],(i%2)+1));
    hDistDiffX[i]->GetYaxis()->SetTitle( Form("%s",projFull[i]));

    hDistDiffXCorr[i] = (TH2F*) hDistDiffX[i]->Clone(Form("hDistDiffXCorr_%s",proj[i]));
	
    hDistDiffPX[i] = new TH2F(Form("hDistDiffPX_%s",proj[i]),Form("%sHRS %s vs P%s",arm,proj[i],projX[i]),100,-XLim[i],XLim[i],100, -DiffLim[i], DiffLim[i]);
    hDistDiffPX[i]->GetXaxis()->SetTitle( Form("P_{%s,%i}",projX[i],(i%2)+1));
    hDistDiffPX[i]->GetYaxis()->SetTitle( Form("%s",projFull[i]));

    hDistDiffPXCorr[i] = (TH2F*) hDistDiffPX[i]->Clone(Form("hDistDiffPXCorr_%s",proj[i]));
     
  }
  

  TH1F* hLoff[NPLANE];
  TF1* gaus_off[NPLANE];
  Double_t mean_off[NPLANE];
  Double_t rms_off[NPLANE];
  
  for(Int_t i = 0; i<NPLANE; i++){  

    hLoff[i] = new TH1F(Form("hLoff_%s",plane[i]),Form("%sHRS timing offset, %s",arm,plane[i]),100,-60,60);
    hLoff[i]->GetXaxis()->SetTitle("timing offset (ns)");

  }


  // 1D global - local angular differences

  TH1F* hmuDiff[NAngle];
  TH1F* hmvDiff[NAngle];

  TF1* fmuDiff[NAngle];
  Double_t mu_mean[NAngle];
  Double_t mu_rms[NAngle];  
  TF1* fmvDiff[NAngle];
  Double_t mv_mean[NAngle];
  Double_t mv_rms[NAngle];
   
  TH1F* hmThDiff[NAngle];
  TH1F* hmPhDiff[NAngle];

  TF1* fmThDiff[NAngle];
  Double_t mTh_mean[NAngle];
  Double_t mTh_rms[NAngle];
  TF1* fmPhDiff[NAngle];
  Double_t mPh_mean[NAngle];
  Double_t mPh_rms[NAngle];

  for(Int_t i = 0; i<NAngle; i++){

    hmuDiff[i] = new TH1F(Form("hmuDiff_%i",i),Form("mu local - global, VDC %i",i),100,-0.2,0.2);
    hmvDiff[i] = new TH1F(Form("hmvDiff_%i",i),Form("mv local - global, VDC %i",i),100,-0.2,0.2);
    hmThDiff[i] = new TH1F(Form("hmThDiff_%i",i),Form("#theta local - global, VDC %i",i),100,-0.2,0.2);
    hmPhDiff[i] = new TH1F(Form("hmPhDiff_%i",i),Form("#phi local - global, VDC %i",i),100,-0.2,0.2); 
  }

  

  TChain* T = Load_more_rootfiles(runno);
  TChain* TA  = Load_more_rootfiles(runno);

  
  TCut T6Cut = "DR.evtypebits&(1<<6)";
  

  Int_t coinc_bins = 300;
  Double_t coinc_start = 0;
  Double_t coinc_end = 200;
  
  TH1F *hCoinc = new TH1F("hCoinc"," S2-Time difference",coinc_bins,coinc_start,coinc_end);
  hCoinc->GetXaxis()->SetTitle("Coinc time (ns)");

  T->Draw("((L.s2.lt_c[L.s2.t_pads]+L.s2.rt_c[L.s2.t_pads])/2-(R.s2.lt_c[R.s2.t_pads]+R.s2.rt_c[R.s2.t_pads])/2)*1e9>>hCoinc",T6Cut,"0");

  // S2 cut

  Double_t CutCen = 163.5;
  Double_t CutWidth = 7.5;
  
  Double_t TLowCut = 156.0;
  Double_t THighCut = 171.0;

  TLine* TLowLine = new TLine(TLowCut,0,TLowCut,1.0*hCoinc->GetMaximum());
  TLowLine->SetLineColor(kMagenta);
  TLowLine->SetLineStyle(kDashed);

  TLine* THighLine = new TLine(THighCut,0,THighCut,1.0*hCoinc->GetMaximum());
  THighLine->SetLineColor(kMagenta);
  THighLine->SetLineStyle(kDashed);



  TCut S2Cut = Form("abs(((L.s2.lt_c[L.s2.t_pads]+L.s2.rt_c[L.s2.t_pads])/2-(R.s2.lt_c[R.s2.t_pads]+R.s2.rt_c[R.s2.t_pads])/2)*1e9-%f) < %f", CutCen, CutWidth);
  

  TH1F *hCoincC = new TH1F("hCoincC"," S2-Time difference with cut",coinc_bins,coinc_start,coinc_end);
  hCoincC->GetXaxis()->SetTitle("Coinc time (ns)");


  T->Draw("((L.s2.lt_c[L.s2.t_pads]+L.s2.rt_c[L.s2.t_pads])/2-(R.s2.lt_c[R.s2.t_pads]+R.s2.rt_c[R.s2.t_pads])/2)*1e9>>hCoincC",T6Cut + S2Cut,"0");    


  
  


  // obtain cut on difference between angles in upper and lower VDCs

  // for gaus fits
  Int_t max_bin;
  Double_t min, max, temp;
  
  TF1* gaus_Ang_X;
  TF1* gaus_Ang_Y;


  const Double_t seperation = 0.3348; // LHRS old = 0.3348 // 0.3327 (old RHRS) //0.3084 (new RHRS)

  T->Draw(Form("((%s-%s)+(%s-%s))/%f>>hXAngleDiff",projROOTPX[0].Data(),projROOTX[1].Data(),projROOTPX[1].Data(),projROOTX[0].Data(),seperation),"","0");
  
  hXAngleDiff->GetXaxis()->SetTitle("(#theta_{VDC,1} - #theta_{VDC,2})");
  
  hXAngleDiff->Fit("gaus","Q0","",-0.03,0.03);
  max_bin = hXAngleDiff->GetMaximumBin();
  min = hXAngleDiff->GetBinCenter(max_bin) - 0.025;
  max = hXAngleDiff->GetBinCenter(max_bin) + 0.025;
  hXAngleDiff->Fit("gaus","Q0","",min,max);
  
  gaus_Ang_X = hXAngleDiff->GetFunction("gaus");
  min = gaus_Ang_X->GetParameter(1) - 1.5*gaus_Ang_X->GetParameter(2);
  max = gaus_Ang_X->GetParameter(1) + 1.5*gaus_Ang_X->GetParameter(2);

  hXAngleDiff->Fit("gaus","Q0","",min,max);
  gaus_Ang_X = hXAngleDiff->GetFunction("gaus");

  
  TCut XAngCut = Form(" abs( ((%s-%s)+(%s-%s))/%f-(%f)) < %f ",projROOTPX[0].Data(),projROOTX[1].Data(),projROOTPX[1].Data(),projROOTX[0].Data(),seperation,gaus_Ang_X->GetParameter(1),1.5*gaus_Ang_X->GetParameter(2));


  T->Draw(Form("((%s-%s)+(%s-%s))>>hXAngleDiffTest",projROOTPX[0].Data(),projROOTX[1].Data(),projROOTPX[1].Data(),projROOTX[0].Data()),XAngCut,"0");
  
  cout << "XAngCut = " << XAngCut << endl;


    
  // equivalent fit for y
  
  T->Draw(Form("((%s-%s)+(%s-%s))/%f>>hYAngleDiff",projROOTPX[2].Data(),projROOTX[3].Data(),projROOTPX[2].Data(),projROOTX[3].Data(),seperation),"","");
  hYAngleDiff->GetXaxis()->SetTitle("(#phi_{VDC,1} - #phi_{VDC,2})");

  hYAngleDiff->Fit("gaus","Q0","",-0.03,0.03);
  max_bin = hYAngleDiff->GetMaximumBin();
  min = hYAngleDiff->GetBinCenter(max_bin) - 0.025;
  max = hYAngleDiff->GetBinCenter(max_bin) + 0.025;
  hYAngleDiff->Fit("gaus","Q0","",min,max);
  
  gaus_Ang_Y = hYAngleDiff->GetFunction("gaus");
  min = gaus_Ang_Y->GetParameter(1) - 1.5*gaus_Ang_Y->GetParameter(2);
  max = gaus_Ang_Y->GetParameter(1) + 1.5*gaus_Ang_Y->GetParameter(2);

  hYAngleDiff->Fit("gaus","Q0","",min,max);
  gaus_Ang_Y = hYAngleDiff->GetFunction("gaus");
  
    
  TCut YAngCut = Form(" abs( ((%s-%s)+(%s-%s))/%f-(%f)) < %f ",projROOTPX[2].Data(),projROOTX[3].Data(),projROOTPX[2].Data(),projROOTX[3].Data(),seperation,gaus_Ang_Y->GetParameter(1),1.5*gaus_Ang_Y->GetParameter(2));

  T->Draw(Form("(%s-%s)+(%s-%s)>>hYAngleDiffTest",projROOTPX[2].Data(),projROOTX[3].Data(),projROOTPX[2].Data(),projROOTX[3].Data()),YAngCut,"");
    
  cout << "YAngCut = " << YAngCut << endl;


  TCut AngCut = XAngCut + YAngCut; //+ T6Cut + S2Cut;




  TCut TimeCut = ""; // t0 cut for all planes
  
  // timing offset  
  
  for(Int_t i = 0; i<NPLANE; i++){
    
    // T->Draw(Form("(L.vdc.%s.trdist-L.vdc.%s.ltrdist)>>hDist_diff_%s",plane[i],plane[i],plane[i]));

    if(i==0){
      TimeCut = Form("%s.vdc.%s.t0!=0",arm,plane[i]);
    }
    else{
      TimeCut += Form("%s.vdc.%s.t0!=0",arm,plane[i]);
    }
    
    TA->Draw(Form("(%s.vdc.%s.t0)*1e9>>hLoff_%s",arm,plane[i],plane[i]),Form("%s.vdc.%s.t0!=0",arm,plane[i]) + AngCut,"0");

    max_bin = hLoff[i]->GetMaximumBin();
    min = hLoff[i]->GetBinCenter(max_bin) - 40;
    max = hLoff[i]->GetBinCenter(max_bin) + 40;
    hLoff[i]->Fit("gaus","Q0","",min,max);
    
    gaus_off[i] = hLoff[i]->GetFunction("gaus");
    min = gaus_off[i]->GetParameter(1) - 1.5*gaus_off[i]->GetParameter(2);
    max = gaus_off[i]->GetParameter(1) + 1.5*gaus_off[i]->GetParameter(2);
    
    hLoff[i]->Fit("gaus","Q0","",min,max);
    gaus_off[i] = hLoff[i]->GetFunction("gaus");
    
    mean_off[i] = gaus_off[i]->GetParameter(1);
    // rms_off[i] = gaus_off[i]->GetParameter(2);
    rms_off[i] = hLoff[i]->GetRMS();
  }


  
  // draw angular differences in lower and upper chambers

  for(Int_t i = 0; i<NAngle; i++){

    T->Draw(Form("%s.vdc.%s.lslope[0]-%s.vdc.%s.slope[0]>>hmuDiff_%i",arm,plane[2*i],arm,plane[2*i],i),AngCut,"0");

    max_bin = hmuDiff[i]->GetMaximumBin();
    min = hmuDiff[i]->GetBinCenter(max_bin) - 0.15;
    max = hmuDiff[i]->GetBinCenter(max_bin) + 0.15;
    hmuDiff[i]->Fit("gaus","Q0","",min,max);
    
    fmuDiff[i] = hmuDiff[i]->GetFunction("gaus");
    min = fmuDiff[i]->GetParameter(1) - 1.5*fmuDiff[i]->GetParameter(2);
    max = fmuDiff[i]->GetParameter(1) + 1.5*fmuDiff[i]->GetParameter(2);
    
    hmuDiff[i]->Fit("gaus","Q0","",min,max);
    fmuDiff[i] = hmuDiff[i]->GetFunction("gaus");
    
    mu_mean[i] = fmuDiff[i]->GetParameter(1);
    mu_rms[i] = fmuDiff[i]->GetParameter(2);
    

    T->Draw(Form("%s.vdc.%s.lslope[0]-%s.vdc.%s.slope[0]>>hmvDiff_%i",arm,plane[2*i+1],arm,plane[2*i+1],i),AngCut,"0");

    max_bin = hmvDiff[i]->GetMaximumBin();
    min = hmvDiff[i]->GetBinCenter(max_bin) - 0.15;
    max = hmvDiff[i]->GetBinCenter(max_bin) + 0.15;
    hmvDiff[i]->Fit("gaus","Q0","",min,max);
    
    fmvDiff[i] = hmvDiff[i]->GetFunction("gaus");
    min = fmvDiff[i]->GetParameter(1) - 1.5*fmvDiff[i]->GetParameter(2);
    max = fmvDiff[i]->GetParameter(1) + 1.5*fmvDiff[i]->GetParameter(2);
    
    hmvDiff[i]->Fit("gaus","Q0","",min,max);
    fmvDiff[i] = hmvDiff[i]->GetFunction("gaus");
    
    mv_mean[i] = fmvDiff[i]->GetParameter(1);
    mv_rms[i] = fmvDiff[i]->GetParameter(2);

  }

  T->Draw(Form("(%s-%s)/%f-%s.tr.d_th[0]>>hmThDiff_%i",projROOTPX[0].Data(),projROOTX[1].Data(),seperation,arm,0),AngCut,"0");

  max_bin = hmThDiff[0]->GetMaximumBin();
  min = hmThDiff[0]->GetBinCenter(max_bin) - 0.15;
  max = hmThDiff[0]->GetBinCenter(max_bin) + 0.15;
  hmThDiff[0]->Fit("gaus","Q0","",min,max);
  
  fmThDiff[0] = hmThDiff[0]->GetFunction("gaus");
  min = fmThDiff[0]->GetParameter(1) - 1.5*fmThDiff[0]->GetParameter(2);
  max = fmThDiff[0]->GetParameter(1) + 1.5*fmThDiff[0]->GetParameter(2);
    
  hmThDiff[0]->Fit("gaus","Q0","",min,max);
  fmThDiff[0] = hmThDiff[0]->GetFunction("gaus");
    
  mTh_mean[0] = fmThDiff[0]->GetParameter(1);
  mTh_rms[0] = fmThDiff[0]->GetParameter(2);


  T->Draw(Form("(%s-%s)/%f-%s.tr.d_th[0]>>hmThDiff_%i",projROOTX[0].Data(),projROOTPX[1].Data(),seperation,arm,1),AngCut,"0");

  max_bin = hmThDiff[1]->GetMaximumBin();
  min = hmThDiff[1]->GetBinCenter(max_bin) - 0.15;
  max = hmThDiff[1]->GetBinCenter(max_bin) + 0.15;
  hmThDiff[1]->Fit("gaus","Q0","",min,max);
  
  fmThDiff[1] = hmThDiff[1]->GetFunction("gaus");
  min = fmThDiff[1]->GetParameter(1) - 1.5*fmThDiff[1]->GetParameter(2);
  max = fmThDiff[1]->GetParameter(1) + 1.5*fmThDiff[1]->GetParameter(2);
    
  hmThDiff[1]->Fit("gaus","Q0","",min,max);
  fmThDiff[1] = hmThDiff[1]->GetFunction("gaus");
    
  mTh_mean[1] = fmThDiff[1]->GetParameter(1);
  mTh_rms[1] = fmThDiff[1]->GetParameter(2);


  T->Draw(Form("(%s-%s)/%f-%s.tr.d_ph[0]>>hmPhDiff_%i",projROOTPX[2].Data(),projROOTX[3].Data(),seperation,arm,0),AngCut,"0");

  max_bin = hmPhDiff[0]->GetMaximumBin();
  min = hmPhDiff[0]->GetBinCenter(max_bin) - 0.15;
  max = hmPhDiff[0]->GetBinCenter(max_bin) + 0.15;
  hmPhDiff[0]->Fit("gaus","Q0","",min,max);
  
  fmPhDiff[0] = hmPhDiff[0]->GetFunction("gaus");
  min = fmPhDiff[0]->GetParameter(1) - 1.5*fmPhDiff[0]->GetParameter(2);
  max = fmPhDiff[0]->GetParameter(1) + 1.5*fmPhDiff[0]->GetParameter(2);
    
  hmPhDiff[0]->Fit("gaus","Q0","",min,max);
  fmPhDiff[0] = hmPhDiff[0]->GetFunction("gaus");
    
  mPh_mean[0] = fmPhDiff[0]->GetParameter(1);
  mPh_rms[0] = fmPhDiff[0]->GetParameter(2);

  T->Draw(Form("(%s-%s)/%f-%s.tr.d_ph[0]>>hmPhDiff_%i",projROOTX[2].Data(),projROOTPX[3].Data(),seperation,arm,1),AngCut,"0");

  max_bin = hmPhDiff[1]->GetMaximumBin();
  min = hmPhDiff[1]->GetBinCenter(max_bin) - 0.15;
  max = hmPhDiff[1]->GetBinCenter(max_bin) + 0.15;
  hmPhDiff[1]->Fit("gaus","Q0","",min,max);
  
  fmPhDiff[1] = hmPhDiff[1]->GetFunction("gaus");
  min = fmPhDiff[1]->GetParameter(1) - 1.5*fmPhDiff[1]->GetParameter(2);
  max = fmPhDiff[1]->GetParameter(1) + 1.5*fmPhDiff[1]->GetParameter(2);
    
  hmPhDiff[1]->Fit("gaus","Q0","",min,max);
  fmPhDiff[1] = hmPhDiff[1]->GetFunction("gaus");
    
  mPh_mean[1] = fmPhDiff[1]->GetParameter(1);
  mPh_rms[1] = fmPhDiff[1]->GetParameter(2);

  
  
  

  // timing offset differences between planes

  const Int_t NComb = 6;
  
  TH1F* hToffDiff[NComb];

  // fit gaussian to normalised differences
  TF1* gaus_off_Norm[NComb];
  Double_t mean_off_Norm[NComb];
  Double_t rms_off_Norm[NComb];
  
  Int_t iToffDiff = 0; //iterator

  Double_t t0_res = 5.7e-9;

  TString ToffDiffSum = "0.5*(";
  
  for(Int_t i = 0; i<NPLANE-1; i++){
    for(Int_t j = i+1; j<NPLANE; j++){


      hToffDiff[iToffDiff] = new TH1F(Form("hToffDiff_%i",iToffDiff),Form("Timing Diff %s - %s",plane[i],plane[j]),100,-5,5);

      hToffDiff[iToffDiff]->GetXaxis()->SetTitle(Form("(t_{0,%s}-t_{0,%s})/#sigma_{t_{0}}",plane[i],plane[j]));
      
      if(iToffDiff == 0){
	ToffDiffSum += Form("(pow((%s.vdc.%s.t0-%s.vdc.%s.t0)/%.1e,2))",arm,plane[i],arm,plane[j],pow( pow(rms_off[i],2) + pow(rms_off[j],2), 0.5)*1e-9);
      }
      else{
	ToffDiffSum += Form(" + (pow((%s.vdc.%s.t0-%s.vdc.%s.t0)/%.1e,2))",arm,plane[i],arm,plane[j],pow( pow(rms_off[i],2) + pow(rms_off[j],2), 0.5)*1e-9);
      }

      TA->Draw(Form("(pow((%s.vdc.%s.t0-%s.vdc.%s.t0)/%.1e,1))>>hToffDiff_%i",arm,plane[i],arm,plane[j],pow( pow(rms_off[i],2) + pow(rms_off[j],2), 0.5)*1e-9,iToffDiff),Form("%s.vdc.%s.t0!=%s.vdc.%s.t0",arm,plane[i],arm,plane[j]) + AngCut,"0");

      
      //T->Draw(Form("(pow((%s.vdc.%s.t0-%s.vdc.%s.t0)/%.1e,1))>>hToffDiff_%i",arm,plane[i],arm,plane[j],(rms_off[i]+rms_off[j])*1e-9,iToffDiff),Form("%s.vdc.%s.t0!=%s.vdc.%s.t0",arm,plane[i],arm,plane[j]) + AngCut,"0");
      
      
      max_bin = hToffDiff[iToffDiff]->GetMaximumBin();
      min = hToffDiff[iToffDiff]->GetBinCenter(max_bin) - 3;
      max = hToffDiff[iToffDiff]->GetBinCenter(max_bin) + 3;
      hToffDiff[iToffDiff]->Fit("gaus","Q0","",min,max);
      
      gaus_off_Norm[iToffDiff] = hToffDiff[iToffDiff]->GetFunction("gaus");
      min = gaus_off_Norm[iToffDiff]->GetParameter(1) - 1.5*gaus_off_Norm[iToffDiff]->GetParameter(2);
      max = gaus_off_Norm[iToffDiff]->GetParameter(1) + 1.5*gaus_off_Norm[iToffDiff]->GetParameter(2);
      
      hToffDiff[iToffDiff]->Fit("gaus","Q0","",min,max);
      gaus_off_Norm[iToffDiff] = hToffDiff[iToffDiff]->GetFunction("gaus");
      
      mean_off_Norm[iToffDiff] = gaus_off_Norm[iToffDiff]->GetParameter(1);
      rms_off_Norm[iToffDiff] = gaus_off_Norm[iToffDiff]->GetParameter(2);
      
      iToffDiff++;

    }    
  }

  ToffDiffSum += ")";

  
  
  
  

  for(Int_t i = 0; i<NPROJ; i++){

    T->Draw(Form("(%s-%s)>>hDistDiff_%s",projROOTX[i].Data(),projROOTPX[i].Data(),proj[i]),AngCut,"0");

    T->Draw(Form("(%s-%s):%s>>hDistDiffX_%s",projROOTX[i].Data(),projROOTPX[i].Data(),projROOTX[i].Data(),proj[i]),AngCut,"0");

    T->Draw(Form("(%s-%s):%s>>hDistDiffPX_%s",projROOTX[i].Data(),projROOTPX[i].Data(),projROOTPX[i].Data(),proj[i]),AngCut,"0");
  }

  


  TF1* gaus_proj[NPLANE];
  Double_t mean_proj[NPLANE];
  Double_t rms_proj[NPLANE];
  
  for(Int_t i = 0; i<NPROJ; i++){
    // cDist->cd(i+1);
    // hDistDiff[i]->Draw();
    
    max_bin = hDistDiff[i]->GetMaximumBin();
    min = hDistDiff[i]->GetBinCenter(max_bin) - 0.02;
    max = hDistDiff[i]->GetBinCenter(max_bin) + 0.02;
    hDistDiff[i]->Fit("gaus","Q0","",min,max);

    gaus_proj[i] = hDistDiff[i]->GetFunction("gaus");
    min = gaus_proj[i]->GetParameter(1) - 1.5*gaus_proj[i]->GetParameter(2);
    max = gaus_proj[i]->GetParameter(1) + 1.5*gaus_proj[i]->GetParameter(2);
    
    hDistDiff[i]->Fit("gaus","Q0","",min,max);
    gaus_proj[i] = hDistDiff[i]->GetFunction("gaus");
    gaus_proj[i]->SetLineColor(kRed);
    //    gaus_proj[i]->Draw("same");
    
    mean_proj[i] = gaus_proj[i]->GetParameter(1);
    rms_proj[i] = gaus_proj[i]->GetParameter(2);

    cout << "For " << projFull[i] << ": " << endl;
    cout << "mean = " << mean_proj[i] << endl;
    cout << "rms = " << rms_proj[i] << endl << endl;    
    
  }
  
  


  // angle from VDC chamber vs 'true'/ global angle
  
  TH2F* hAngleVs[NAngle];

  const char Angle[NAngle] = {};

  Double_t AngleSign[NAngle] = {1.0,-1.0}; // sign flip in angle formula for projecting from chamber 1->2 compared to 2->1



  // line to show perfect correlation

  TLine* LXAngle = new TLine(0.7,0.7,1.4,1.4);
  LXAngle->SetLineColor(kRed);


  hAngleVs[0] = new TH2F("hAngleVs_0","VDC (1) Angle vs true",100,0.7,1.4,100,0.7,1.4);
  hAngleVs[0]->GetXaxis()->SetTitle("#theta_{g} ");
  hAngleVs[0]->GetYaxis()->SetTitle("#theta_{VDC,1} ");
  
  
  T->Draw(Form("(%s-%s)/%f:%s.tr.d_th[0]>>hAngleVs_%i",projROOTPX[0].Data(),projROOTX[1].Data(),seperation,arm,0),AngCut,"0");

  
  hAngleVs[1] = new TH2F("hAngleVs_1","VDC (2) Angle vs true",100,0.7,1.4,100,0.7,1.4);
  hAngleVs[1]->GetXaxis()->SetTitle("#theta_{g} ");
  hAngleVs[1]->GetYaxis()->SetTitle("#theta_{VDC,2} ");

  T->Draw(Form("(%s-%s)/%f:%s.tr.d_th[0]>>hAngleVs_%i",projROOTX[0].Data(),projROOTPX[1].Data(),seperation,arm,1),AngCut,"0");
  
  


  TH2F* hAngleDiffVs[NAngle];

  hAngleDiffVs[0] = new TH2F("hAngleDiffVs_0","(VDC (1) Angle - true) vs VDC(1) Angle",100,0.7,1.4,100,-0.3,0.3);
  hAngleDiffVs[0]->GetXaxis()->SetTitle("#theta_{VDC,1} ");
  hAngleDiffVs[0]->GetYaxis()->SetTitle("#theta_{VDC,1} - #theta_{g} ");
  
  T->Draw(Form("(%s-%s)/%f-%s.tr.d_th[0]:((%s-%s)/%f)>>hAngleDiffVs_%i",projROOTPX[0].Data(),projROOTX[1].Data(),seperation,arm,projROOTPX[0].Data(),projROOTX[1].Data(),seperation,0),AngCut,"0");
  // 30,80
  hAngleDiffVs[0]->FitSlicesY(0,20,70,10,"QN");
    
  
  hAngleDiffVs[1] = new TH2F("hAngleDiffVs_1","(VDC (2) Angle - true) vs VDC(1) Angle",100,0.7,1.4,100,-0.3,0.3);
  hAngleDiffVs[1]->GetXaxis()->SetTitle("#theta_{VDC,2} ");
  hAngleDiffVs[1]->GetYaxis()->SetTitle("#theta_{VDC,2} - #theta_{g} ");
  
  T->Draw(Form("(%s-%s)/%f-%s.tr.d_th[0]:((%s-%s)/%f)>>hAngleDiffVs_%i",projROOTX[0].Data(),projROOTPX[1].Data(),seperation,arm,projROOTX[0].Data(),projROOTPX[1].Data(),seperation,1),AngCut,"0");

  hAngleDiffVs[1]->FitSlicesY(0,20,70,10,"QN");
  


  Double_t XSigma_up = 2.0; // how many sigma to cut away from mean of 'real' distance distrib

  Double_t XSigma_down = 2.0; // how many sigma to cut away from mean of 'real' distance distrib

  // slice distributions
  TH1D *h_means[NAngle];
  TH1D *h_sigma[NAngle];
  TH1D *h_const[NAngle];

  // fit funtionc for mean
  TF1 *f_means[NAngle];
  
  // plot sigma limits above and below means of distribution
  TH1D *h_fit_up[NAngle];
  TH1D *h_fit_low[NAngle];

  
  // fit results
  Double_t c[NAngle] = {0};
  Double_t m[NAngle] = {0};
  Double_t m_1[NAngle] = {0};
  
  for(Int_t i = 0; i<NAngle; i++){
      
    
    h_const[i] = (TH1D*)gDirectory->Get(Form("hAngleDiffVs_%i_0",i));
        
    h_means[i] = (TH1D*)gDirectory->Get(Form("hAngleDiffVs_%i_1",i));

    h_means[i]->Fit("pol1","Q0");    
    f_means[i] = h_means[i]->GetFunction("pol1");


    c[i] = f_means[i]->GetParameter(0);
    m[i] = f_means[i]->GetParameter(1);
    m_1[i] = 1.0-m[i];


    c[i] = 0;
    m[i] = 1.0;
    m_1[i] = 1.0;
    
    cout << "For " << i << endl;
    cout << "c = " << c[i] << endl;
    cout << "m = " << m[i] << endl << endl;
    
    h_sigma[i] = (TH1D*)gDirectory->Get(Form("hAngleDiffVs_%i_2",i));   

    h_fit_up[i] = (TH1D*) h_means[i]->Clone();
    h_fit_up[i]->Add(h_sigma[i],XSigma_up);


    h_fit_low[i] = (TH1D*) h_means[i]->Clone();
    h_fit_low[i]->Add(h_sigma[i],-XSigma_down);


  }
  
  Double_t mX = (m[0]+m[1])/2.; // average of two slopes
  Double_t mX_1 = 1.0-mX; // average of two slopes
  mX_1 = 1.0;

  

  TH2F* hAngleVsCorr[NAngle];

  hAngleVsCorr[0] = new TH2F("hAngleVsCorr_0","VDC (1) Angle vs true (corrected)",100,0.7,1.4,100,0.7,1.4);
  hAngleVsCorr[0]->GetXaxis()->SetTitle("#theta_{g} ");
  hAngleVsCorr[0]->GetYaxis()->SetTitle("#theta_{VDC,1} ");
  
  
  T->Draw(Form("(%f*(%s-%s)-%f*%f)/%f:%s.tr.d_th[0]>>hAngleVsCorr_%i",mX_1,projROOTPX[0].Data(),projROOTX[1].Data(),c[0],seperation,seperation,arm,0),AngCut,"0");

  
  hAngleVsCorr[1] = new TH2F("hAngleVsCorr_1","VDC (2) Angle vs true (corrected)",100,0.7,1.4,100,0.7,1.4);
  hAngleVsCorr[1]->GetXaxis()->SetTitle("#theta_{g} ");
  hAngleVsCorr[1]->GetYaxis()->SetTitle("#theta_{VDC,2} ");

  T->Draw(Form("(%f*(%s-%s)-%f*%f)/%f:%s.tr.d_th[0]>>hAngleVsCorr_%i",mX_1,projROOTX[0].Data(),projROOTPX[1].Data(),c[1],seperation,seperation,arm,1),AngCut,"0");
  


  TH2F* hAngleDiffVsCorr[NAngle];
    
  hAngleDiffVsCorr[0] = new TH2F("hAngleDiffVsCorr_0","(VDC (1) Angle - true) vs VDC(1) Angle (corrected)",100,0.7,1.4,100,-0.3,0.3);
  hAngleDiffVsCorr[0]->GetXaxis()->SetTitle("#theta_{VDC,1} ");
  hAngleDiffVsCorr[0]->GetYaxis()->SetTitle("#theta_{VDC,1} - #theta_{g} ");
  
  T->Draw(Form("(%f*(%s-%s)-%f*%f)/%f-%s.tr.d_th[0]:((%s-%s)/%f)>>hAngleDiffVsCorr_%i",mX_1,projROOTPX[0].Data(),projROOTX[1].Data(),c[0],seperation,seperation,arm,projROOTPX[0].Data(),projROOTX[1].Data(),seperation,0),AngCut,"0");
    
  
  hAngleDiffVsCorr[1] = new TH2F("hAngleDiffVsCorr_1","(VDC (2) Angle - true) vs VDC(1) Angle (corrected)",100,0.7,1.4,100,-0.3,0.3);
  hAngleDiffVsCorr[1]->GetXaxis()->SetTitle("#theta_{VDC,2} ");
  hAngleDiffVsCorr[1]->GetYaxis()->SetTitle("#theta_{VDC,2} - #theta_{g} ");
  
  T->Draw(Form("(%f*(%s-%s)-%f*%f)/%f-%s.tr.d_th[0]:((%s-%s)/%f)>>hAngleDiffVsCorr_%i",mX_1,projROOTX[0].Data(),projROOTPX[1].Data(),c[1],seperation,seperation,arm,projROOTX[0].Data(),projROOTPX[1].Data(),seperation,1),AngCut,"0");



    
  TF1* gaus_X_Wid[NPLANE];
  Double_t meanX[NAngle];
  Double_t gausX[NAngle];
  
  TF1* gaus_X_Wid_Norm[NPLANE];
  Double_t meanX_Norm[NAngle];
  Double_t gausX_Norm[NAngle];

  TString X_Norm_sum = "("; // string to sum all x-px components
  
  for(Int_t i = 0; i<NAngle; i++){

    T->Draw(Form("%s-%s-%f*(%s-%s)+%f*%f*%f>>hDistDiffCorr_%s",projROOTX[i].Data(),projROOTX[(i+1)%2].Data(),mX_1,projROOTPX[(i)%2].Data(),projROOTX[(i+1)%2].Data(),AngleSign[i],c[i],seperation,proj[i]),AngCut,"0");

    cout << "For " << i << endl;
    
    
    //    cout << "c pre-correction = " << c[i] << endl;
    //    c[i] = -AngleSign[i]*hDistDiffCorr[i]->GetBinCenter(hDistDiffCorr[i]->GetMaximumBin());
    // c[i] -= AngleSign[i]*hDistDiffCorr[i]->GetBinCenter(hDistDiffCorr[i]->GetMaximumBin());
    //    cout << "c post-correction = " <<  c[i] << endl;
    

    //    T->Draw(Form("%s-%s-%f*(%s-%s)+%f*%f>>hDistDiffCorr_%s",projROOTX[i].Data(),projROOTX[(i+1)%2].Data(),mX_1,projROOTPX[(i)%2].Data(),projROOTX[(i+1)%2].Data(),c[i],AngleSign[i],proj[i]),AngCut,"0");

    max_bin = hDistDiffCorr[i]->GetMaximumBin();
    min = hDistDiffCorr[i]->GetBinCenter(max_bin) - 0.02;
    max = hDistDiffCorr[i]->GetBinCenter(max_bin) + 0.02;
    hDistDiffCorr[i]->Fit("gaus","Q0","",min,max);

    gaus_X_Wid[i] = hDistDiffCorr[i]->GetFunction("gaus");
    min = gaus_X_Wid[i]->GetParameter(1) - 1.5*gaus_X_Wid[i]->GetParameter(2);
    max = gaus_X_Wid[i]->GetParameter(1) + 1.5*gaus_X_Wid[i]->GetParameter(2);
    
    hDistDiffCorr[i]->Fit("gaus","Q0","",min,max);
    

    gaus_X_Wid[i] = hDistDiffCorr[i]->GetFunction("gaus");
    meanX[i] = gaus_X_Wid[i]->GetParameter(1);
    gausX[i] = gaus_X_Wid[i]->GetParameter(2);  
    

    T->Draw(Form("(%s-%s-%f*(%s-%s)+%f*%f*%f)/%f>>hDistDiffCorrNorm_%s",projROOTX[i].Data(),projROOTX[(i+1)%2].Data(),mX_1,projROOTPX[(i)%2].Data(),projROOTX[(i+1)%2].Data(),AngleSign[i],c[i],seperation,gausX[i],proj[i]),AngCut,"0");
 

    if(i==0){
      X_Norm_sum += Form("pow((%s-%s-%f*(%s-%s)+%f*%f*%f)/%f,2)",projROOTX[i].Data(),projROOTX[(i+1)%2].Data(),mX_1,projROOTPX[(i)%2].Data(),projROOTX[(i+1)%2].Data(),AngleSign[i],c[i],seperation,gausX[i]);
    }
    else{
      X_Norm_sum += Form(" + pow((%s-%s-%f*(%s-%s)+%f*%f*%f)/%f,2)",projROOTX[i].Data(),projROOTX[(i+1)%2].Data(),mX_1,projROOTPX[(i)%2].Data(),projROOTX[(i+1)%2].Data(),AngleSign[i],c[i],seperation,gausX[i]);

    }
    
    max_bin = hDistDiffCorrNorm[i]->GetMaximumBin();
    min = hDistDiffCorrNorm[i]->GetBinCenter(max_bin) - 3;
    max = hDistDiffCorrNorm[i]->GetBinCenter(max_bin) + 3;
    hDistDiffCorrNorm[i]->Fit("gaus","Q0","",min,max);

    gaus_X_Wid_Norm[i] = hDistDiffCorrNorm[i]->GetFunction("gaus");
    min = gaus_X_Wid_Norm[i]->GetParameter(1) - 1.5*gaus_X_Wid_Norm[i]->GetParameter(2);
    max = gaus_X_Wid_Norm[i]->GetParameter(1) + 1.5*gaus_X_Wid_Norm[i]->GetParameter(2);
    
    hDistDiffCorrNorm[i]->Fit("gaus","Q0","",min,max);
    

    gaus_X_Wid_Norm[i] = hDistDiffCorrNorm[i]->GetFunction("gaus");
    meanX_Norm[i] = gaus_X_Wid_Norm[i]->GetParameter(1);
    gausX_Norm[i] = gaus_X_Wid_Norm[i]->GetParameter(2);  
    
    
    
    /*
      T->Draw(Form("%s-%s-%f*(%s-%s)+%f*%f:%s>>hDistDiffXCorr_%s",projROOTX[i].Data(),projROOTX[(i+1)%2].Data(),mX_1,projROOTPX[(i)%2].Data(),projROOTX[(i+1)%2].Data(),c[i],AngleSign[i],projROOTX[i].Data(),proj[i]),AngCut,"0");
    
      T->Draw(Form("%s-%s-%f*(%s-%s)+%f*%f:%s>>hDistDiffPXCorr_%s",projROOTX[i].Data(),projROOTX[(i+1)%2].Data(),mX_1,projROOTPX[(i)%2].Data(),projROOTX[(i+1)%2].Data(),c[i],AngleSign[i],projROOTPX[i].Data(),proj[i]),AngCut,"0");    
    */
    
  }

  


  // Y equivalents of angle plots

  
  TLine* LYAngle = new TLine(-0.08,-0.08,0.08,0.08);
  LYAngle->SetLineColor(kRed);
  
  TH2F* hAngleYVs[NAngle];
  
  hAngleYVs[0] = new TH2F("hAngleYVs_0","VDC (1) Angle vs true",100,-0.08,0.08,100,-0.08,0.08);
  hAngleYVs[0]->GetXaxis()->SetTitle("#phi_{g} ");
  hAngleYVs[0]->GetYaxis()->SetTitle("#phi_{VDC,1} ");

  T->Draw(Form("(%s-%s)/%f:%s.tr.d_ph[0]>>hAngleYVs_%i",projROOTPX[2].Data(),projROOTX[3].Data(),seperation,arm,0),AngCut,"0");
    
  
  hAngleYVs[1] = new TH2F("hAngleYVs_1","VDC (2) Angle vs true",100,-0.08,0.08,100,-0.08,0.08);
  hAngleYVs[1]->GetXaxis()->SetTitle("#phi_{g} ");
  hAngleYVs[1]->GetYaxis()->SetTitle("#phi_{VDC,2} ");

  T->Draw(Form("(%s-%s)/%f:%s.tr.d_ph[0]>>hAngleYVs_%i",projROOTX[2].Data(),projROOTPX[3].Data(),seperation,arm,1),AngCut,"0");
  


  TH2F* hAngleYDiffVs[NAngle];

  hAngleYDiffVs[0] = new TH2F("hAngleYDiffVs_0","(VDC (1) Angle - true) vs VDC(1) Angle",100,-0.08,0.08,100,-0.1,0.1);
  hAngleYDiffVs[0]->GetXaxis()->SetTitle("#phi_{VDC,1} ");
  hAngleYDiffVs[0]->GetYaxis()->SetTitle("#phi_{VDC,1} - #phi_{g} ");
  
  T->Draw(Form("(%s-%s)/%f-%s.tr.d_ph[0]:((%s-%s)/%f)>>hAngleYDiffVs_%i",projROOTPX[2].Data(),projROOTX[3].Data(),seperation,arm,projROOTPX[2].Data(),projROOTX[3].Data(),seperation,0),AngCut,"0");

  hAngleYDiffVs[0]->FitSlicesY(0,30,55,10,"QN");
    
  
  hAngleYDiffVs[1] = new TH2F("hAngleYDiffVs_1","(VDC (2) Angle - true) vs VDC(2) Angle",100,-0.08,0.08,100,-0.1,0.1);
  hAngleYDiffVs[1]->GetXaxis()->SetTitle("#phi_{VDC,2} ");
  hAngleYDiffVs[1]->GetYaxis()->SetTitle("#phi_{VDC,2} - #phi_{g} ");
  
  T->Draw(Form("(%s-%s)/%f-%s.tr.d_ph[0]:((%s-%s)/%f)>>hAngleYDiffVs_%i",projROOTX[2].Data(),projROOTPX[3].Data(),seperation,arm,projROOTX[2].Data(),projROOTPX[3].Data(),seperation,1),AngCut,"0");

  hAngleYDiffVs[1]->FitSlicesY(0,30,55,10,"QN");
  

  
  // slice distributions
  TH1D *h_meansY[NAngle];
  TH1D *h_sigmaY[NAngle];
  TH1D *h_constY[NAngle];

  // fit funtionc for mean
  TF1 *f_meansY[NAngle];
  
  // plot sigma limits above and below means of distribution
  TH1D *h_fit_upY[NAngle];
  TH1D *h_fit_lowY[NAngle];

  
  // fit results
  Double_t cY[NAngle] = {0};
  Double_t mY_b[NAngle] = {0};
  Double_t mY_1_b[NAngle] = {0};
  
  for(Int_t i = 0; i<NAngle; i++){
    
    h_constY[i] = (TH1D*)gDirectory->Get(Form("hAngleYDiffVs_%i_0",i));    

    h_meansY[i] = (TH1D*)gDirectory->Get(Form("hAngleYDiffVs_%i_1",i));

    h_meansY[i]->Fit("pol1","Q");    
    f_meansY[i] = h_meansY[i]->GetFunction("pol1");

    cY[i] = f_meansY[i]->GetParameter(0);
    mY_b[i] = f_meansY[i]->GetParameter(1);
    mY_1_b[i] = 1.0-mY_b[i];

    
    cY[i] = 0;
    mY_b[i] = 1.0;
    mY_1_b[i] = 1.0;
    
    cout << "(Y) For " << i << endl;
    cout << "c = " << cY[i] << endl;
    cout << "m = " << mY_b[i] << endl << endl;
    
    h_sigmaY[i] = (TH1D*)gDirectory->Get(Form("hAngleYDiffVs_%i_2",i));   

    h_fit_upY[i] = (TH1D*) h_meansY[i]->Clone();
    h_fit_upY[i]->Add(h_sigmaY[i],XSigma_up);
    h_fit_upY[i]->SetLineColor(kRed);

    h_fit_lowY[i] = (TH1D*) h_meansY[i]->Clone();
    h_fit_lowY[i]->Add(h_sigmaY[i],-XSigma_down);
    h_fit_lowY[i]->SetLineColor(kRed);

    h_meansY[i]->SetLineColor(kBlack);

  }
  
  Double_t mY = (mY_b[0]+mY_b[1])/2.; // average of two slopes
  Double_t mY_1 = 1.0-mY; // average of two slopes
  mY_1 = 1.0;
 
  
  cout << "mY_1 = " << mY_1 << endl;

  TH2F* hAngleYVsCorr[NAngle];

  hAngleYVsCorr[0] = new TH2F("hAngleYVsCorr_0","VDC (1) Angle vs true (corrected)",100,-0.08,0.08,100,-0.08,0.08);
  hAngleYVsCorr[0]->GetXaxis()->SetTitle("#phi_{g} ");
  hAngleYVsCorr[0]->GetYaxis()->SetTitle("#phi_{VDC,1} ");
  
  
  T->Draw(Form("(%f*(%s-%s)-%f*%f)/%f:%s.tr.d_ph[0]>>hAngleYVsCorr_%i",mY_1,projROOTPX[2].Data(),projROOTX[3].Data(),cY[0],seperation,seperation,arm,0),AngCut,"0");

  
  hAngleYVsCorr[1] = new TH2F("hAngleYVsCorr_1","VDC (2) Angle vs true (corrected)",100,-0.08,0.08,100,-0.08,0.08);
  hAngleYVsCorr[1]->GetXaxis()->SetTitle("#phi_{g} ");
  hAngleYVsCorr[1]->GetYaxis()->SetTitle("#phi_{VDC,2} ");

  T->Draw(Form("(%f*(%s-%s)-%f*%f)/%f:%s.tr.d_ph[0]>>hAngleYVsCorr_%i",mY_1,projROOTX[2].Data(),projROOTPX[3].Data(),cY[1],seperation,seperation,arm,1),AngCut,"0");
  


  TH2F* hAngleYDiffVsCorr[NAngle];
    
  hAngleYDiffVsCorr[0] = new TH2F("hAngleYDiffVsCorr_0","(VDC (1) Angle - true) vs VDC(1) Angle (corrected)",100,-0.08,0.08,100,-0.1,0.1);
  hAngleYDiffVsCorr[0]->GetXaxis()->SetTitle("#phi_{VDC,1} ");
  hAngleYDiffVsCorr[0]->GetYaxis()->SetTitle("#phi_{VDC,1} - #phi_{g} ");
  
  T->Draw(Form("(%f*(%s-%s)-%f*%f)/%f-%s.tr.d_ph[0]:((%s-%s)/%f)>>hAngleYDiffVsCorr_%i",mY_1,projROOTPX[2].Data(),projROOTX[3].Data(),cY[0],seperation,seperation,arm,projROOTPX[2].Data(),projROOTX[3].Data(),seperation,0),AngCut,"0");
    
  
  hAngleYDiffVsCorr[1] = new TH2F("hAngleYDiffVsCorr_1","(VDC (2) Angle - true) vs VDC(1) Angle (corrected)",100,-0.08,0.08,100,-0.1,0.1);
  hAngleYDiffVsCorr[1]->GetXaxis()->SetTitle("#phi_{VDC,2} ");
  hAngleYDiffVsCorr[1]->GetYaxis()->SetTitle("#phi_{VDC,2} - #phi_{g} ");
  
  T->Draw(Form("(%f*(%s-%s)-%f*%f)/%f-%s.tr.d_ph[0]:((%s-%s)/%f)>>hAngleYDiffVsCorr_%i",mY_1,projROOTX[2].Data(),projROOTPX[3].Data(),cY[0],seperation,seperation,arm,projROOTX[2].Data(),projROOTPX[3].Data(),seperation,1),AngCut,"0");
  

  Double_t meanY[NAngle];
  Double_t gausY[NAngle];

  Double_t meanY_Norm[NAngle];
  Double_t gausY_Norm[NAngle];
  
  for(Int_t i = 0; i<NAngle; i++){

    T->Draw(Form("%s-%s-%f*(%s-%s)+%f*%f*%f>>hDistDiffCorr_%s",projROOTX[2+i].Data(),projROOTX[2+(i+1)%2].Data(),mY_1,projROOTPX[2+(i)%2].Data(),projROOTX[2+(i+1)%2].Data(),AngleSign[i],cY[i],seperation,proj[i+2]),AngCut,"0");


    X_Norm_sum += Form(" + pow((%s-%s-%f*(%s-%s)+%f*%f*%f)/%f,2)",projROOTX[i].Data(),projROOTX[(i+1)%2].Data(),mX_1,projROOTPX[(i)%2].Data(),projROOTX[(i+1)%2].Data(),AngleSign[i],c[i],seperation,gausX[i]);
    
    cout << "For " << i << endl;
    // cout << "Y: c pre-correction = " << cY[i] << endl;
    // cY[i] -= AngleSign[i]*hDistDiffCorr[i+2]->GetBinCenter(hDistDiffCorr[i+2]->GetMaximumBin());
    // cout << "c post-correction = " <<  cY[i] << endl;
    //    cY[i] += AngleSign[i]*hDistDiffCorr[i+2]->GetBinCenter(hDistDiffCorr[i+2]->GetMaximumBin());
        
    // T->Draw(Form("%s-%s-%f*(%s-%s)+%f*%f>>hDistDiffCorr_%s",projROOTX[2+i].Data(),projROOTX[2+(i+1)%2].Data(),mY_1,projROOTPX[2+(i)%2].Data(),projROOTX[2+(i+1)%2].Data(),cY[i],AngleSign[i],proj[i+2]),AngCut,"0");
    

    max_bin = hDistDiffCorr[i+2]->GetMaximumBin();
    min = hDistDiffCorr[i+2]->GetBinCenter(max_bin) - 0.02;
    max = hDistDiffCorr[i+2]->GetBinCenter(max_bin) + 0.02;
    hDistDiffCorr[i+2]->Fit("gaus","Q0","",min,max);

    gaus_X_Wid[i+2] = hDistDiffCorr[i+2]->GetFunction("gaus");
    min = gaus_X_Wid[i+2]->GetParameter(1) - 1.5*gaus_X_Wid[i+2]->GetParameter(2);
    max = gaus_X_Wid[i+2]->GetParameter(1) + 1.5*gaus_X_Wid[i+2]->GetParameter(2);
    
    hDistDiffCorr[i+2]->Fit("gaus","Q0","",min,max);
    gaus_X_Wid[i+2] = hDistDiffCorr[i+2]->GetFunction("gaus");

    meanY[i] = gaus_X_Wid[i+2]->GetParameter(1);
    gausY[i] = gaus_X_Wid[i+2]->GetParameter(2);


    T->Draw(Form("(%s-%s-%f*(%s-%s)+%f*%f*%f)/%f>>hDistDiffCorrNorm_%s",projROOTX[2+i].Data(),projROOTX[2+(i+1)%2].Data(),mY_1,projROOTPX[2+(i)%2].Data(),projROOTX[2+(i+1)%2].Data(),AngleSign[i],cY[i],seperation,gausY[i],proj[i+2]),AngCut,"0");
    
    max_bin = hDistDiffCorrNorm[i+2]->GetMaximumBin();
    min = hDistDiffCorrNorm[i+2]->GetBinCenter(max_bin) - 3;
    max = hDistDiffCorrNorm[i+2]->GetBinCenter(max_bin) + 3;
    hDistDiffCorrNorm[i+2]->Fit("gaus","Q0","",min,max);

    gaus_X_Wid_Norm[i+2] = hDistDiffCorrNorm[i+2]->GetFunction("gaus");
    min = gaus_X_Wid_Norm[i+2]->GetParameter(1) - 1.5*gaus_X_Wid_Norm[i+2]->GetParameter(2);
    max = gaus_X_Wid_Norm[i+2]->GetParameter(1) + 1.5*gaus_X_Wid_Norm[i+2]->GetParameter(2);
    
    hDistDiffCorrNorm[i+2]->Fit("gaus","Q0","",min,max);
    

    gaus_X_Wid_Norm[i+2] = hDistDiffCorrNorm[i+2]->GetFunction("gaus");
    meanY_Norm[i] = gaus_X_Wid_Norm[i+2]->GetParameter(1);
    gausY_Norm[i] = gaus_X_Wid_Norm[i+2]->GetParameter(2);  
    
    
        
    /*
      T->Draw(Form("%s-%s-%f*(%s-%s):%s>>hDistDiffXCorr_%s",projROOTX[2+i].Data(),projROOTX[2+(i+1)%2].Data(),mY_1,projROOTPX[2+(i)%2].Data(),projROOTX[2+(i+1)%2].Data(),projROOTX[2+i].Data(),proj[i+2]),AngCut,"0");

      T->Draw(Form("%s-%s-%f*(%s-%s)+%f*%f:%s>>hDistDiffPXCorr_%s",projROOTX[2+i].Data(),projROOTX[2+(i+1)%2].Data(),mY_1,projROOTPX[2+(i)%2].Data(),projROOTX[2+(i+1)%2].Data(),projROOTPX[2+i].Data(),proj[i+2]),AngCut,"0");          
    */
 
  }
  
  X_Norm_sum += ")";

  cout << "ToffDiffSum = " << ToffDiffSum << endl;
  cout << "X_Norm_sum = " << X_Norm_sum << endl;

  TString GOMtot = X_Norm_sum + " + " + ToffDiffSum;
  cout << "GOMtot = " << GOMtot << endl<<endl;

  cout << "TimeCut = " << TimeCut << endl;
  
  // Form overall 'goodness of match' measure from combining X and Y
  
  TH1F* hGOM = new TH1F("hGOM","Goodness of Match",100,0,100);
  hGOM->GetXaxis()->SetTitle("GOM");

  T->Draw(Form("%s>>hGOM",GOMtot.Data()),TimeCut,"0");


  TH1F* hGOMNorm = new TH1F("hGOMNorm","Goodness of Match (normalised histogram)",100,0,100);
  hGOMNorm->GetXaxis()->SetTitle("GOM");

  T->Draw(Form("%s>>hGOMNorm",GOMtot.Data()),TimeCut,"0");
  
  
  

  // calc 'alpha' values for x and y

  Double_t alp_X[NAngle] = {0};
  Double_t alp_Y[NAngle] = {0};

  for(Int_t i = 0; i<NAngle; i++){

    alp_X[i] = 1/(1-m[i]);
    alp_Y[i] = 1/(1-mY_b[i]);

    cout << "X_alpha_" << i+1 << " = " << alp_X[i] << endl;
    cout << "Y_alpha_" << i+1 << " = " << alp_Y[i] << endl;
  }


  
  // plots


  TCanvas* c_Coinc = new TCanvas("c_Coinc","c_Coinc", 1000, 1000);
  c_Coinc->Divide(2,1);

  c_Coinc->cd(1);
  hCoinc->Draw("hist");

  TLowLine->Draw("same");
  THighLine->Draw("same");

  c_Coinc->cd(2);
  hCoincC->Draw();

  TLowLine->Draw("same");
  THighLine->Draw("same");
  

  TCanvas* cAngleDiff = new TCanvas("cAngleDiff", "Angle Differences", 640, 480);
  cAngleDiff->Divide(2,1);


  cAngleDiff->cd(1);  
  hXAngleDiff->Draw();
  
  gaus_Ang_X->SetLineColor(kRed);

  gaus_Ang_X->Draw("same");

  
  cAngleDiff->cd(2);  
  hYAngleDiff->Draw();
  
  gaus_Ang_Y->SetLineColor(kRed);

  gaus_Ang_Y->Draw("same");

  
  
  TCanvas* cAngle = new TCanvas("cAngle", "X Angles", 640, 480);
  cAngle->Divide(2,2);
  
  for(Int_t i = 0; i<NAngle; i++){
    
    cAngle->cd(i+1);
    hAngleVs[i]->Draw("colz");    
    LXAngle->Draw();

    cAngle->cd(i+3);
    hAngleVsCorr[i]->Draw("colz");
    LXAngle->Draw();
  }


  TCanvas* cAngleDiffVs = new TCanvas("cAngleDiffVs", "X Angle Diff", 640, 480);
  cAngleDiffVs->Divide(2,2);
  
  for(Int_t i = 0; i<NAngle; i++){
    
    cAngleDiffVs->cd(i+1);
    hAngleDiffVs[i]->Draw("colz");

    cAngleDiffVs->cd(i+3);
    hAngleDiffVsCorr[i]->Draw("colz");

  }


  TCanvas* cAngleSlice[NAngle];

    
  for(Int_t i = 0; i<NAngle; i++){
    
    cAngleSlice[i]  = new TCanvas(Form("cAngleSlice_%i",i), Form("X Angle Slices %i",i), 640, 480);
  
    cAngleSlice[i]->Divide(2,2);

    
    cAngleSlice[i]->cd(1);
    hAngleDiffVs[i]->Draw("colz");
  
    cAngleSlice[i]->cd(2);
    h_const[i]->Draw();

    cAngleSlice[i]->cd(3);
    h_means[i]->Draw();
    f_means[i]->SetLineColor(kMagenta);

    
    cAngleSlice[i]->cd(4);
    h_sigma[i]->Draw();


    cAngleSlice[i]->cd(1);
    h_means[i]->SetLineColor(kBlack);
    h_means[i]->Draw("same");
    h_fit_up[i]->SetLineColor(kRed);
    h_fit_up[i]->Draw("hist  L same");
    h_fit_low[i]->SetLineColor(kRed);
    h_fit_low[i]->Draw("hist L same");
  }
  
  
  TCanvas* cDist = new TCanvas("cProj", "Position Projection Difference", 640, 480);
  cDist->Divide(2,2);

  TCanvas* cDistCorr = new TCanvas("cProjCorr", "Position Projection Difference (Corrected)", 640, 480);
  cDistCorr->Divide(2,2);

  TCanvas* cDistCorrNorm = new TCanvas("cProjCorrNorm", "Normalised Position Projection Difference (Corrected)", 640, 480);
  cDistCorrNorm->Divide(2,2);
  

  THStack* hDistDiffStack[NPROJ];
  
  TLegend *leg_DistDiff =  new TLegend(.65,.65,.9,.9,"Key");
  leg_DistDiff->SetFillColor(0);
  leg_DistDiff->SetTextSize(0.025);

  
  TLatex* texDistDiff[NPROJ];
  TLatex* texDistDiffNorm[NPROJ];
  
  for(Int_t i = 0; i<NPROJ; i++){
    
    cDistCorr->cd(i+1);
    hDistDiffCorr[i]->Draw();    

    gaus_X_Wid[i]->SetLineColor(kRed);
    gaus_X_Wid[i]->Draw("same");
    
    texDistDiff[i] = new TLatex( 0.65, 0.4, Form("#sigma = %.3f mm",1e3*gaus_X_Wid[i]->GetParameter(2)));
    texDistDiff[i]->SetNDC(1);
    texDistDiff[i]->SetTextFont(42);
    texDistDiff[i]->SetTextColor(1);
    texDistDiff[i]->SetTextSize(0.065);
    texDistDiff[i]->Draw();
    

    hDistDiffStack[i]  = new THStack(Form("hDistDiffStack_%i",i),Form("Dist Difference; %s",projFull[i]));
    
    cDist->cd(i+1);
    hDistDiff[i]->SetLineColor(kRed);


    hDistDiffStack[i]->Add(hDistDiffCorr[i]);
    hDistDiffStack[i]->Add(hDistDiff[i]);
    
    hDistDiffStack[i]->Draw("nostack");

    if(i==0){
      leg_DistDiff->AddEntry(hDistDiff[i],"uncorrected","l");
      leg_DistDiff->AddEntry(hDistDiffCorr[i],"corrected","l");
    }
    leg_DistDiff->Draw("same");
    

    // cDistXCorr->cd(i+1);
    // hDistDiffXCorr[i]->Draw("colz");
    
    // cDistPXCorr->cd(i+1);
    // hDistDiffPXCorr[i]->Draw("colz");

    

    cDistCorrNorm->cd(i+1);
    hDistDiffCorrNorm[i]->Draw();


    gaus_X_Wid_Norm[i]->SetLineColor(kRed);
    gaus_X_Wid_Norm[i]->Draw("same");
    
    texDistDiffNorm[i] = new TLatex( 0.65, 0.4, Form("#sigma_{norm} = %.3f",gaus_X_Wid_Norm[i]->GetParameter(2)));
    texDistDiffNorm[i]->SetNDC(1);
    texDistDiffNorm[i]->SetTextFont(42);
    texDistDiffNorm[i]->SetTextColor(1);
    texDistDiffNorm[i]->SetTextSize(0.065);
    texDistDiffNorm[i]->Draw();
    
    
    
  }

  

  
  TCanvas* cAngleY = new TCanvas("cAngleY", "Y Angles", 640, 480);
  cAngleY->Divide(2,2);

  for(Int_t i = 0; i<NAngle; i++){
    
    cAngleY->cd(i+1);
    hAngleYVs[i]->Draw("colz");
    LYAngle->Draw();

    cAngleY->cd(i+3);
    hAngleYVsCorr[i]->Draw("colz");
    LYAngle->Draw();
    
  }


  TCanvas* cAngleYDiffVs = new TCanvas("cAngleYDiffVs", "Y Angle Diff", 640, 480);
  cAngleYDiffVs->Divide(2,2);
  
  for(Int_t i = 0; i<NAngle; i++){
    
    cAngleYDiffVs->cd(i+1);
    hAngleYDiffVs[i]->Draw("colz");

    cAngleYDiffVs->cd(i+3);
    hAngleYDiffVsCorr[i]->Draw("colz");

  }


  TCanvas* cAngleYSlice[NAngle];

    
  for(Int_t i = 0; i<NAngle; i++){
    
    cAngleYSlice[i]  = new TCanvas(Form("cAngleYSlice_%i",i), Form("Y Angle Slices %i",i), 640, 480);
  
    cAngleYSlice[i]->Divide(2,2);

    
    cAngleYSlice[i]->cd(1);
    hAngleYDiffVs[i]->Draw("colz");
  
    cAngleYSlice[i]->cd(2);
    h_constY[i]->Draw();

    cAngleYSlice[i]->cd(3);
    h_meansY[i]->Draw();
    f_meansY[i]->SetLineColor(kMagenta);

    
    cAngleYSlice[i]->cd(4);
    h_sigmaY[i]->Draw();


    cAngleYSlice[i]->cd(1);
    h_meansY[i]->SetLineColor(kBlack);
    h_meansY[i]->Draw("same");
    h_fit_upY[i]->SetLineColor(kRed);
    h_fit_upY[i]->Draw("hist  L same");
    h_fit_lowY[i]->SetLineColor(kRed);
    h_fit_lowY[i]->Draw("hist L same");
  }


  TCanvas* cSlopeDiff[NAngle];

  TLatex* texmu[NAngle];
  TLatex* texmv[NAngle];
  TLatex* texmTh[NAngle];
  TLatex* texmPh[NAngle];
  
  for(Int_t i = 0; i<NAngle; i++){
    
    cSlopeDiff[i] = new TCanvas(Form("cSlopeDiff_%i",i),Form("Slope/ angle differences, chamber %i",i),640,480);

    cSlopeDiff[i]->Divide(2,2);

    cSlopeDiff[i]->cd(1);
    hmuDiff[i]->Draw();
    fmuDiff[i]->SetLineColor(kRed);
    fmuDiff[i]->Draw("same");

    texmu[i] = new TLatex( 0.65, 0.4, Form("#sigma = %.3f",fmuDiff[i]->GetParameter(2)));
    texmu[i]->SetNDC(1);
    texmu[i]->SetTextFont(42);
    texmu[i]->SetTextColor(1);
    texmu[i]->SetTextSize(0.065);
    texmu[i]->Draw();
    

    cSlopeDiff[i]->cd(2);
    hmvDiff[i]->Draw();
    fmvDiff[i]->SetLineColor(kRed);
    fmvDiff[i]->Draw("same");

    texmv[i] = new TLatex( 0.65, 0.4, Form("#sigma = %.3f",fmvDiff[i]->GetParameter(2)));
    texmv[i]->SetNDC(1);
    texmv[i]->SetTextFont(42);
    texmv[i]->SetTextColor(1);
    texmv[i]->SetTextSize(0.065);
    texmv[i]->Draw();

    cSlopeDiff[i]->cd(3);
    hmThDiff[i]->Draw();
    fmThDiff[i]->SetLineColor(kRed);
    fmThDiff[i]->Draw("same");

    texmTh[i] = new TLatex( 0.65, 0.4, Form("#sigma = %.3f",fmThDiff[i]->GetParameter(2)));
    texmTh[i]->SetNDC(1);
    texmTh[i]->SetTextFont(42);
    texmTh[i]->SetTextColor(1);
    texmTh[i]->SetTextSize(0.065);
    texmTh[i]->Draw();

    cSlopeDiff[i]->cd(4);
    hmPhDiff[i]->Draw();
    fmPhDiff[i]->SetLineColor(kRed);
    fmPhDiff[i]->Draw("same");

    texmPh[i] = new TLatex( 0.65, 0.4, Form("#sigma = %.3f",fmPhDiff[i]->GetParameter(2)));
    texmPh[i]->SetNDC(1);
    texmPh[i]->SetTextFont(42);
    texmPh[i]->SetTextColor(1);
    texmPh[i]->SetTextSize(0.065);
    texmPh[i]->Draw();
    
  }
  

  
  TCanvas* cT0 = new TCanvas("cT0", "timing offset", 640, 480);  
  cT0->Divide(2,2);
  
  TLatex* texToff[NPLANE];
    
  cout << "timingoffset" << endl;
  for(Int_t i = 0; i<NPLANE; i++){

    cT0->cd(i+1);
    hLoff[i]->Draw();
    gaus_off[i]->SetLineColor(kRed);
    gaus_off[i]->Draw("same");

    cout << plane[i] << ": mean = " << mean_off[i] << ", rms = " << rms_off[i] << endl;


    texToff[i] = new TLatex( 0.65, 0.4, Form("#sigma = %.3f ns",gaus_off[i]->GetParameter(2)));
    texToff[i]->SetNDC(1);
    texToff[i]->SetTextFont(42);
    texToff[i]->SetTextColor(1);
    texToff[i]->SetTextSize(0.065);
    texToff[i]->Draw();
    
  }


  TCanvas* cToffDiff = new TCanvas("cToffDiff", "timing offset differences", 640, 480);
  cToffDiff->Divide(2,3);

  TLatex* texToffDiff[NComb];
  
  for(Int_t i = 0; i<NComb; i++){

    cToffDiff->cd(i+1);

    hToffDiff[i]->Draw();

    gaus_off_Norm[i]->SetLineColor(kRed);
    gaus_off_Norm[i]->Draw("same");
    
    texToffDiff[i] = new TLatex( 0.65, 0.4, Form("#sigma_{norm} = %.3f",rms_off_Norm[i]));
    texToffDiff[i]->SetNDC(1);
    texToffDiff[i]->SetTextFont(42);
    texToffDiff[i]->SetTextColor(1);
    texToffDiff[i]->SetTextSize(0.065);
    texToffDiff[i]->Draw();

  }
  

  TCanvas* cGOM = new TCanvas("cGOM", "Goodness of Match", 640, 480);
  cGOM->Divide(2,1);
  
  cGOM->cd(1);
  hGOM->Draw();

  cGOM->cd(2);
  Double_t factor = 1;
  hGOMNorm->Scale(factor/hGOMNorm->GetEntries());
  //  hGOMNorm->Draw("hist");

  // calc quantile
  Int_t NQuant = 3; // number of quantils to be calcuated
  Double_t QuantLev[3] = {0.68,0.95,0.997};
  Double_t QuantVal[3];
  hGOMNorm->GetQuantiles(3,QuantVal,QuantLev);

  cout << "QuantVal = " << QuantVal[0] << ", " << QuantVal[1] << ", " << QuantVal[2] << endl;
  
  TF1* fChi2 = new TF1("fChi2","ROOT::Math::chisquared_pdf(x,7,0)",.01,100);
  //  fChi2->SetParameter(0,20);
  fChi2->SetLineColor(kRed);
  fChi2->Draw();

  
  hGOMNorm->Draw("hist same");
  
  
  /*

  


    TCanvas* cDistXCorr = new TCanvas("cProjXCorr", "Position Projection Difference vs Position (Corrected)", 640, 480);
    cDistXCorr->Divide(2,2);

    TCanvas* cDistPXCorr = new TCanvas("cProjPXCorr", "Position Projection Difference vs Projection (Corrected)", 640, 480);
    cDistPXCorr->Divide(2,2);


    THStack* hDistDiffStack[NPROJ];

    TLegend *leg_DistDiff =  new TLegend(.65,.65,.9,.9,"Key");
    leg_DistDiff->SetFillColor(0);
    leg_DistDiff->SetTextSize(0.025);

  
    TLatex* texDistDiff[NPROJ];
  
    for(Int_t i = 0; i<NPROJ; i++){

    cDistCorr->cd(i+1);
    hDistDiffCorr[i]->Draw();    

    gaus_X_Wid[i]->SetLineColor(kRed);
    gaus_X_Wid[i]->Draw("same");

    //    texDistDiff[i] = new TLatex( 0.5, 0.3, Form("#splitline{#color[4]{#bar{nhits} = %.2f}}{#color[2]{#bar{nhits} = %.2f}}",hLNHits_A->GetMean(),hLNHits_F->GetMean()));
    texDistDiff[i] = new TLatex( 0.65, 0.4, Form("#sigma = %.3f mm",1e3*gaus_X_Wid[i]->GetParameter(2)));
    texDistDiff[i]->SetNDC(1);
    texDistDiff[i]->SetTextFont(42);
    texDistDiff[i]->SetTextColor(1);
    texDistDiff[i]->SetTextSize(0.065);
    texDistDiff[i]->Draw();
    

    hDistDiffStack[i]  = new THStack(Form("hDistDiffStack_%i",i),Form("Dist Difference; %s",projFull[i]));
    
    cDist->cd(i+1);
    hDistDiff[i]->SetLineColor(kRed);


    hDistDiffStack[i]->Add(hDistDiffCorr[i]);
    hDistDiffStack[i]->Add(hDistDiff[i]);
    
    hDistDiffStack[i]->Draw("nostack");

    if(i==0){
    leg_DistDiff->AddEntry(hDistDiff[i],"uncorrected","l");
    leg_DistDiff->AddEntry(hDistDiffCorr[i],"corrected","l");
    }
    leg_DistDiff->Draw("same");
    
    //      hDistDiff[i]->Draw("same");
    
    
    //      hDistDiffStack[i]->Add(hDistDiffCorr[i]);

    //      leg_DistDiff->AddEntry(hDistDiff[i],"uncorrected","l");
    

    //      hDistDiff[i]->Draw();
      


    cDistXCorr->cd(i+1);
    hDistDiffXCorr[i]->Draw("colz");
    
    cDistPXCorr->cd(i+1);
    hDistDiffPXCorr[i]->Draw("colz");
    
    }



  

  
  
    TCanvas* cAngleCorr = new TCanvas("cAngleCorr", "AnglesCorr", 640, 480);
    cAngleCorr->Divide(2,1);  

    for(Int_t i = 0; i<NAngle; i++){
    
    cAngleCorr->cd(i+1);
    hAngleVsCorr[i]->Draw("colz");
      
    }


    // Y angle correlation plots
  
    TCanvas* cAngleYSlice[NAngle];

  
  
    for(Int_t i = 0; i<NAngle; i++){
    
    cAngleYSlice[i]  = new TCanvas(Form("cAngleYSlice_%i",i), Form("Angle Slices (Y) %i",i), 640, 480);
  
    cAngleYSlice[i]->Divide(2,2);
  
    
    cAngleYSlice[i]->cd(1);
    hAngleYVs[i]->Draw("colz");
  
  
    cAngleYSlice[i]->cd(2);
    h_constY[i]->Draw();

    cAngleYSlice[i]->cd(3);
    h_meansY[i]->Draw();
    f_meansY[i]->SetLineColor(kMagenta);

      
    cAngleYSlice[i]->cd(4);
    h_sigmaY[i]->Draw();


    cAngleYSlice[i]->cd(1);
    h_meansY[i]->SetLineColor(kBlack);
    h_meansY[i]->Draw("same");
    h_fit_upY[i]->SetLineColor(kRed);
    h_fit_upY[i]->Draw("hist  L same");
    h_fit_lowY[i]->SetLineColor(kRed);
    h_fit_lowY[i]->Draw("hist L same");
    }
  
  

  
    TCanvas* cDistX = new TCanvas("cProjX", "Position Projection Difference vs Position", 640, 480);
    cDistX->Divide(2,2);

    TCanvas* cDistPX = new TCanvas("cProjPX", "Position Projection Difference vs Projection", 640, 480);
    cDistPX->Divide(2,2);


 
    for(Int_t i = 0; i<NPROJ; i++){
    cDistX->cd(i+1);
    hDistDiffX[i]->Draw("colz");

    cDistPX->cd(i+1);
    hDistDiffPX[i]->Draw("colz");
    
    }

  */
    
  // save plots

  cAngleDiff->SaveAs(Form("plots/VDC_sep/%s_%s_%i_Angle_12_Diff.png",title,arm,runno));
  
  cAngle->SaveAs(Form("plots/VDC_sep/%s_%s_%i_XAngles.png",title,arm,runno));
  cAngleDiffVs->SaveAs(Form("plots/VDC_sep/%s_%s_%i_XAnglesDiff.png",title,arm,runno));

  cDist->SaveAs(Form("plots/VDC_sep/%s_%s_%i_bothDist.png",title,arm,runno));

  cDistCorr->SaveAs(Form("plots/VDC_sep/%s_%s_%i_bothDistCorr.png",title,arm,runno));

  cDistCorrNorm->SaveAs(Form("plots/VDC_sep/%s_%s_%i_bothDistCorrNorm.png",title,arm,runno));

  cAngleY->SaveAs(Form("plots/VDC_sep/%s_%s_%i_YAngles.png",title,arm,runno));
  cAngleYDiffVs->SaveAs(Form("plots/VDC_sep/%s_%s_%i_YAnglesDiff.png",title,arm,runno));

   for(Int_t i = 0; i<NAngle; i++){
     cSlopeDiff[i]->SaveAs(Form("plots/VDC_sep/%s_%s_%i_SlopeDiff_%i.png",title,arm,runno,i));
   }
  
  cT0->SaveAs(Form("plots/VDC_sep/%s_%s_%i_Toff.png",title,arm,runno));

  cToffDiff->SaveAs(Form("plots/VDC_sep/%s_%s_%i_ToffDiff.png",title,arm,runno));

  cGOM->SaveAs(Form("plots/VDC_sep/%s_%s_%i_GOM.png",title,arm,runno));
  
  /*
    cAngleDiff->SaveAs(Form("plots/VDC_sep/%s_%i_AngleDiff.png",arm,runno));
  
    cDist->SaveAs(Form("plots/VDC_sep/%s_%i_Dist_proj.png",arm,runno));

    cDistCorr->SaveAs(Form("plots/VDC_sep/%s_%i_Dist_proj_Corr.png",arm,runno));

    cDistX->SaveAs(Form("plots/VDC_sep/%s_%i_Dist_proj_v_Dist.png",arm,runno));

    cDistPX->SaveAs(Form("plots/VDC_sep/%s_%i_Dist_proj_v_proj.png",arm,runno));

    cDistXCorr->SaveAs(Form("plots/VDC_sep/%s_%i_Dist_proj_v_Dist_Corr.png",arm,runno));

    cDistPXCorr->SaveAs(Form("plots/VDC_sep/%s_%i_Dist_proj_v_proj_Corr.png",arm,runno));
  

    cT0->SaveAs(Form("plots/VDC_sep/%s_%i_toff.png",arm,runno));

    for(Int_t i = 0; i<NAngle; i++){
    cAngleSlice[i]->SaveAs(Form("plots/VDC_sep/%s_%i_X_Slice_%i.png",arm,runno,i));
    cAngleYSlice[i]->SaveAs(Form("plots/VDC_sep/%s_%i_X_Slice_%i.png",arm,runno,i));
    }
  
  */


  // save values to DB

  std::ofstream* outp = new std::ofstream;

  outp->open(Form("DB/VDC_sep/db_%s_Res_%s.vdc.%d.dat", arm, title, runno) );



  for(Int_t i = 0; i < NAngle; i++ ){
    *outp<<arm<<".vdc.uv"<<i+1<<".X.m ="<<endl;
    *outp<<1-m[i]<<endl<<endl;

    *outp<<arm<<".vdc.uv"<<i+1<<".X.c ="<<endl;
    *outp<<c[i]<<endl<<endl;
  }


  for(Int_t i = 0; i < NAngle; i++ ){
    *outp<<arm<<".vdc."<<i+1<<".X.res ="<<endl;
    *outp<<gausX[i]<<endl<<endl;
  }
    
  

  for(Int_t i = 0; i < NAngle; i++ ){
    *outp<<arm<<".vdc.uv"<<i+1<<".Y.m ="<<endl;
    *outp<<1-mY_b[i]<<endl<<endl;

    *outp<<arm<<".vdc.uv"<<i+1<<".Y.c ="<<endl;
    *outp<<cY[i]<<endl<<endl;
  }

  for(Int_t i = 0; i < NAngle; i++ ){
    *outp<<arm<<".vdc."<<i+1<<".Y.res ="<<endl;
    *outp<<gausY[i]<<endl<<endl;
  }

  

  for(Int_t i = 0; i<NPLANE; i++){

    *outp<<arm<<".vdc."<<plane[i]<<".t0_res ="<<endl;
    *outp<<rms_off[i]<<endl<<endl;
    
    *outp<<arm<<".vdc."<<plane[i]<<".t0_offset ="<<endl;
    *outp<<mean_off[i]<<endl<<endl;
  }



  outp->close();


  
}
