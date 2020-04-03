//////////////////////////////////////////////////////////////////////////
//
// rootlogon.C
//
// Load Lib, paths and possible scripts to the analyzer upon start
//
//////////////////////////////////////////////////////////////////////////
//
// Author : Jin Huang <mailto:jinhuang@jlab.org>    Mar 2008
//
//////////////////////////////////////////////////////////////////////////
{
// 	printf("\nrootlogon.C: Loading BigBite Library...");
// 	gSystem->Load("libBigBite.so");
// 	gSystem->Load("onlineGUI/online_C.so");
//	gSystem->Load("$ROOTSYS/lib/libThread.so");
     //	gSystem->Load("libRICH.so");
// 	gSystem->Load("libBigBite.so");

//     //Load more libs here, if necessary.
//     //Make sure it's in path of $LD_LIBRARY_PATH
// 	printf("\nrootlogon.C: Adding include directories...");
//	gSystem->AddIncludePath("-I$ANALYZER");
//	gInterpreter->AddIncludePath("$ANALYZER/");
//

    const char* macros[] =
    {   //list of scripts to load, end with "\0"
      //      "Cuts.C",
      //      "HallA_style.cxx",
      //      "RootFileLoader.C",
	//
	//"LOpticsCalib.C",
	//   "LOpticsCalib.C",
      "SaveCanvas.C",

      "LOpticsOpt.C+",
      // 	//	"LOpticsOptScriptVertex.C",
      // 	//	"cut_L.C",
      "LOpticsOptScript.C",
	//          "cut_R.C",
        "\0"
    };

    if (*macros[0]!=0)
        cout << "\nrootlogon.C: Loading macro's:" << endl;
    for(UInt_t i=0; *macros[i]!=0; i++) {
        cout << "\t " << macros[i] << endl;
        gROOT->LoadMacro(macros[i]);
    }


    TVirtualFitter::SetDefaultFitter("Minuit2");
    gROOT->ProcessLine(".x SetOKStyle.C");
    //    std::cout << "Execute command: " << std::endl;
    //    gSystem->Exec("cd $ROOTSYS/bin; source thisroot.csh; cd -");
    //    gSystem->Exec("cd $ROOTSYS/bin; bash thisroot.sh; cd -");
    
 
}
