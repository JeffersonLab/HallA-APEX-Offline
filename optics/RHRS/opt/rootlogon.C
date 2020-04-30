{
    gSystem->AddIncludePath("-I$ANALYZER/src/");
    gSystem->AddIncludePath("-I$ANALYZER/hana_decode");
    gInterpreter->AddIncludePath("$ANALYZER/src/");
    gInterpreter->AddIncludePath("$ANALYZER/hana_decode");

    gROOT->LoadMacro("ROpticsOpt.C+");
    gROOT->LoadMacro("ROpticsOptScript.C");
    gROOT->LoadMacro("cut_R.C");

    
    gStyle->SetPalette(1);
}
