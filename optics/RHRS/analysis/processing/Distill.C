void Distill(){

  //Macro to take full root file and get only needed variables
  
  TChain * t = new TChain("T");
  t->Add("../root/old/apex_4647*.root");

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
  t->SetBranchStatus("R.tr.r_x",1);
  t->SetBranchStatus("R.tr.r_y",1);
  t->SetBranchStatus("R.tr.r_th",1);
  t->SetBranchStatus("R.tr.r_ph",1);
  
  t->SetBranchStatus("R.tr.vz",1);

  TFile* f = new TFile("../apex_4647.root","recreate");
  TTree* tree = t->CloneTree(0);
  
  tree->CopyEntries(t);

  f->Write();
  
}
