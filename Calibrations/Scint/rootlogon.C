{
    gSystem->AddIncludePath("-I$ANALYZER/src/");
    gSystem->AddIncludePath("-I$ANALYZER/hana_decode");
    gInterpreter->AddIncludePath("$ANALYZER/src/");
    gInterpreter->AddIncludePath("$ANALYZER/hana_decode");


    gROOT->LoadMacro("CsvParser.C+");
    
    gStyle->SetPalette(1);
}
