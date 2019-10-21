void projection(){

  //Macro takes graphical cuts and makes phi and theta projections from them and writes some results to the csv file containing all the information
  
  TFile* tcuts = new TFile("../Sieve/xfp_-10_10/apex_4647.root.FullCut.root","READ");
  TChain * t = new TChain("T");
  t->Add("/home/sean/Grad/Research/APEX/Rootfiles/apex_4647.root");
  
  TCut GeneralCut = "R.tr.n==1 && (R.cer.asum_c>500) && abs(R.tr.r_x) < 0.10";
  double phi_rms, th_rms;
  int  stat;

  ifstream cutcsv("../Sieve/xfp_-10_10/apex_4647.root.cuts_full.csv");
  ofstream cutcsvnew("apex_4647.root.EllipseCuts.Sieve.csv");
  string line;
  cutcsvnew<<fixed<<setprecision(1);
  getline(cutcsv,line);
  cutcsvnew<<line<<endl;
  getline(cutcsv,line);
  cutcsvnew<<line<<endl;
  
  for(int n_col = 0; n_col < 27; n_col++){
    for(int n_row = 0; n_row < 17; n_row++){
      TCutG* g = NULL;
      tcuts->GetObject(Form("hcut_R_1_%d_%d",n_col,n_row), g);
      getline(cutcsv,line);
      
      if (!g){
	cutcsvnew<<line<<",0,0,0"<<endl;
	continue;
      }
      
      TString name = Form("ID = %d:%d",n_col,n_row);
      TString name2 = Form("ID  = %d:%d",n_col,n_row);

       for(int i = 0; i<g->GetN();i++){
	  g->GetX()[i] /= 100;
	  g->GetY()[i] /= 100;
	}
       
      TCut id_cut = TCut(Form("hcut_R_1_%d_%d",n_col,n_row));

      
      t->Draw("R.tr.tg_ph*1000>>" + name, GeneralCut && id_cut,"");
      //t->Draw("Sieve.y*100>>" + name, GeneralCut && id_cut,"");
      

      
      TH1F *htemp = (TH1F*)gPad->GetPrimitive(name);   
      htemp->SetTitle("");
      htemp->GetXaxis()->SetTitle("Target #phi (rad)");
      htemp->GetYaxis()->SetTitle("Entries");
      stat = htemp->GetEntries();
      cout<<stat<<endl;
      TCanvas *c = new TCanvas("c","",800,600);
      htemp->Draw();
      phi_rms = htemp->GetRMS();
      //c->SaveAs(Form("./projections/ph_%d_%d.gif",n_col,n_row));
      
      t->Draw("R.tr.tg_th*1000>>" + name2, GeneralCut && id_cut,"");
      //t->Draw("Sieve.x*100>>" + name2, GeneralCut && id_cut,"");

      TH1F *htemp2 = (TH1F*)gPad->GetPrimitive(name2);   
      htemp2->SetTitle("");
      htemp2->GetXaxis()->SetTitle("Target #theta (rad)");
      htemp2->GetYaxis()->SetTitle("Entries");

      htemp2->Draw();
      th_rms = htemp2->GetRMS();
      //c->SaveAs(Form("./projections/th_%d_%d.gif",n_col,n_row));
      
      cutcsvnew<<line<<","<<phi_rms<<","<<th_rms<<","<<stat<<endl;
    }
  }
  
  
}
