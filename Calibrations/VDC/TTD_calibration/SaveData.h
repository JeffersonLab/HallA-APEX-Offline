// Free function for saving TTD DB

int SaveNewTTDData(vector<Float_t> table, Double_t NBins, TString arm, TString plane, Int_t runnumber){


  TString filename = Form("DB/lookup_tables/db_%s_%s_lookup_test_TTD.vdc.%d.dat",arm.Data(), plane.Data(), runnumber);
  
  std::ofstream* outp = new std::ofstream;

  outp->open(filename.Data());

  *outp << Form("%s.vdc.%s.ttd_table.nbins = ",arm.Data(),plane.Data()) << endl <<NBins << endl << endl;


  *outp << Form("%s.vdc.%s.ttd_table.table =",arm.Data(),plane.Data()) << endl;

  for(Int_t j=0; j<NBins; j++){
    if (j%10 == 0 && j>0){
      *outp << endl;
    }
    *outp << table[j] << " ";
  }
  *outp << endl;

  
  return 1;
}
