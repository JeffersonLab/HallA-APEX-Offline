void Distill(){

  //Macro to take full root file and get only needed variables
  
  TString run = "4653";

  TChain * t = new TChain("T");
  TString rootfile = "/lustre19/expphy/volatile/halla/apex/jeffas/apex_root/Rootfiles/";
  t->Add(rootfile + "apex_"+run+"*.root");

  int entries = t->GetEntries();

  t->SetBranchStatus("*",0);
  t->SetBranchStatus("R.tr.n",1);
  t->SetBranchStatus("R.tr.x",1);
  t->SetBranchStatus("R.tr.y",1);
  t->SetBranchStatus("R.tr.th",1);
  t->SetBranchStatus("R.tr.ph",1);
  t->SetBranchStatus("R.tr.p",1);
  t->SetBranchStatus("R.cer.asum_c",1);
  t->SetBranchStatus("R.s0.nthit",1);
  t->SetBranchStatus("R.tr.tg_y",1);
  t->SetBranchStatus("R.tr.tg_th",1);
  t->SetBranchStatus("R.tr.tg_ph",1);
  t->SetBranchStatus("R.tr.tg_dp",1);
  t->SetBranchStatus("Rrb.x",1);
  t->SetBranchStatus("Rrb.y",1);
  t->SetBranchStatus("Rurb.x",1);
  t->SetBranchStatus("Rurb.y",1);
  t->SetBranchStatus("R.tr.r_x",1);
  t->SetBranchStatus("R.tr.r_y",1);
  t->SetBranchStatus("R.tr.r_th",1);
  t->SetBranchStatus("R.tr.r_ph",1);
  
  t->SetBranchStatus("R.tr.vz",1);
  t->SetBranchStatus("Rrb.BPMA.x",1);
  t->SetBranchStatus("Rrb.BPMA.y",1);
  t->SetBranchStatus("Rrb.BPMB.x",1);
  t->SetBranchStatus("Rrb.BPMB.y",1);
  

  int max_n = 2000000;

  TFile* f = new TFile(rootfile + "/Distilled/apex_"+run+".root","recreate");

  TTree* tree = t->CloneTree(0);
  
  int file_stop = max_n;
  if(max_n > entries) file_stop = entries;

  for(int i=0;i < file_stop;i++){
    t->GetEntry(i);
    tree->Fill();
  }

  f->Write();
  f->Close();

  if(file_stop == max_n){
    TFile* f2 = new TFile(rootfile + "/Distilled/apex_"+run+"_1.root","recreate");
    TTree* tree2 = t->CloneTree(0);

    for(int i=file_stop;i < entries;i++){
      t->GetEntry(i);
      tree2->Fill();
    }

    f2->Write();
    f2->Close();

  }


}
