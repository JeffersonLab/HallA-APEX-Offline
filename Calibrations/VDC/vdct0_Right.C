R__LOAD_LIBRARY(GmpVDCt0Calib_C.so)

#include "GmpVDCt0Calib.h"
#include "Load_more_rootfiles.C"
#include "file_def.h"

  void vdct0_Right(int run) {

    //gStyle->SetOptStat(0);
    gStyle->SetTitleX(0.5);
    gStyle->SetTitleY(1.05);
    gStyle->SetTitleH(0.20);
    gStyle->SetTitleW(0.50);
    gStyle->SetLabelSize(0.04,"X");
    gStyle->SetLabelSize(0.04,"Y");
    gStyle->SetTitleSize(0.04,"X");
    gStyle->SetTitleSize(0.04,"Y");

    //    const TString ROOTfilePath = "/chafs1/work1/tritium/rootfiles_1.6_root6";

    //    TChain* T = new TChain("T");

    // std::ostringstream str;
    // str << ROOTfilePath << "/triton_R_" << run;

    // TString basename = str.str().c_str();
    // TString rootfile = basename + "-1.6.root";

    // Long_t i=0;
    // while ( !gSystem->AccessPathName(rootfile.Data()) ) {
    //     T->Add(rootfile.Data());
    //     cout << "ROOT file " << rootfile << " added to TChain." << endl;
    //     i++;
    //     //rootfile = basename + Form("_%d",i) + ".root";
    //     rootfile = basename + "_" + i + "-1.6.root";
    // }

    // if (i==0) {
    //     cerr << "The specified run does not exist or has not been replayed." << endl;
    //     return;
    // }

    TChain* T = Load_more_rootfiles(run);

    TString plotname = Form("plots/R_vdct0_%d.pdf",run);
    GmpVDCt0Calib* vdct0 = new GmpVDCt0Calib;
    //    vdct0->Calibrate(T,0,"rhrs.evtypebits==2",plotname);

    vdct0->Calibrate(T,0,"(DR.evtype==5)",plotname); 
    
    vdct0->SaveDatabase(Form("DB/db_R_slope.vdc.%d.dat",run));

    Double_t min = 10000, max = 0;
    for (UInt_t ii=0; ii<GmpVDCt0Calib::GetNumberOfPlanes(); ii++) {
        for (UInt_t jj=0; jj<GmpVDCt0Calib::GetNumberOfWiresPerPlane(); jj++) {
            Double_t offset = vdct0->GetOffset(ii,jj);
            if (offset<min) min=offset;
            if (offset>max) max=offset;
        }
    }

    TCanvas* c1 = new TCanvas;
    //TH1F* hFrame = c1->DrawFrame(0,1480,GmpVDCt0Calib::GetNumberOfWiresPerPlane(),1560,"VDC t0 calibration reuslt ( (DL.evtypebits&1<<3) == 1<<3 );Wire number;Time offset (TDC channel)");
    TH1F* hFrame = c1->DrawFrame(0,2900,GmpVDCt0Calib::GetNumberOfWiresPerPlane(),3000,"VDC t0 calibration reuslt;Wire number;Time offset (TDC channel)");
    hFrame->GetXaxis()->CenterTitle();
    hFrame->GetYaxis()->CenterTitle();
    hFrame->GetYaxis()->SetTitleOffset(1.05);

    TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    for (UInt_t ii=0; ii<GmpVDCt0Calib::GetNumberOfPlanes(); ii++) {
        TGraph* graph = new TGraph;
        graph->SetName(vdct0->GetPlaneName(ii));
        UInt_t npoint = 0;	
        for (UInt_t j=0; j<GmpVDCt0Calib::GetNumberOfWiresPerPlane(); j++) {
            Double_t offset = vdct0->GetOffset(ii,j);
            if (offset>0){
	      graph->SetPoint(npoint++,j,offset);
	    }
	    else{
	    // add default offsets	      
	      	     
	      if(j/(vdct0->GetNumberOfWiresPerGroup()) > 0 ){
		// case of not first group in plane and previous group offset is non-zero
		offset = vdct0->GetOffset(ii,j-16); // use previous  group offset
		graph->SetPoint(npoint++,j,offset);

	      }
	      else if( (j/vdct0->GetNumberOfWiresPerGroup()) <((GmpVDCt0Calib::GetNumberOfWiresPerPlane()/16)-1) ){
		// case of not last group in plane and next group offset is non-zero
		offset = vdct0->GetOffset(ii,j+16);
		graph->SetPoint(npoint++,j,offset);	
	      }
	      
	    }
	    
        }
        graph->SetMarkerColor(ii+1);
        graph->SetMarkerStyle(ii+20);
        graph->Draw("P");
        leg->AddEntry(vdct0->GetPlaneName(ii),vdct0->GetPlaneName(ii),"p");
    }

    leg->Draw();
    c1->Update();

    //    c1->SaveAs(Form("vdct0plot_%d.C",run));
}
