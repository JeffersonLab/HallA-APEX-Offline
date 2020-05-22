{

  gSystem->AddIncludePath("-I$ANALYZER/src/");
  gSystem->AddIncludePath("-I$ANALYZER/hana_decode");
  gInterpreter->AddIncludePath("$ANALYZER/src/");
  gInterpreter->AddIncludePath("$ANALYZER/hana_decode");

  gROOT->LoadMacro("ROpticsOpt.C+");
  gROOT->LoadMacro("replay.C");
  gROOT->LoadMacro("Distill.C");
  //gROOT->LoadMacro("proj.C");
  
  gStyle->SetPalette(1);
}
