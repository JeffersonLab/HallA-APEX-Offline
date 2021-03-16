{
    gSystem->AddIncludePath("-I$ANALYZER/src/");
    gSystem->AddIncludePath("-I$ANALYZER/hana_decode");
    gInterpreter->AddIncludePath("$ANALYZER/src/");
    gInterpreter->AddIncludePath("$ANALYZER/hana_decode");


    gROOT->LoadMacro("NoiseCut.C+");
    gROOT->LoadMacro("TTDTable.C+");
    
    gStyle->SetPalette(1);
}
