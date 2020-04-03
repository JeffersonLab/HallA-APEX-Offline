///////////////////////////////////////////////////////////////////////////////
//
// minimiser_comp
//
// Run LOpticsOptScript.C with different minimisers and algorithms and save
// results (rms value for theta, phi, y) to csv file
//
//
///////////////////////////////////////////////////////////////////////////////


void minimiser_comp(){



  // create list of minimisers and associated algorithms
  
 std:vector<std::pair<char *, char*> > Min_list = {{(char*)"Minuit2",(char *)"Migrad"},{(char *)"Minuit2",(char *)"Simplex"},{(char *)"Minuit2",(char *)"Combined"},{(char *)"Minuit2",(char *)"Scan"},{(char *)"Minuit2",(char *)"Fumili"},{(char *)"GSLMultiMin",(char *)"ConjugateFR"},{(char *)"GSLMultiMin",(char *)"ConjugatePR"},{(char *)"GSLMultiMin",(char *)"BFGS"},{(char *)"GSLMultiMin",(char *)"BFGS2"},{(char *)"GSLMultiMin",(char *)"SteepestDescent"},{(char *)"GSLMultiMin",(char *)""},{(char *)"GSLSimAn",(char *)""}};

  // std:vector<std::pair<char *, char*> > Min_list = {{(char*)"Minuit2",(char *)"Migrad"},{(char *)"Minuit2",(char *)"Simplex"}};

  // std:vector<std::pair<char *, char*> > Min_list = {{(char *)"Minuit2",(char *)"Simplex"}};
  
  
  Int_t i = 0;
  for( auto const& pair : Min_list){
    cout << "For pair " << i << " : " << std::get<0>(pair) << " & " << std::get<1>(pair) << endl; 


    
    LOpticsOptScript("theta","test_j",(TString)std::get<0>(pair) + (char*)"_" + std::get<1>(pair) + "_T",std::get<0>(pair),std::get<1>(pair),1);


    LOpticsOptScript("phi",(TString)std::get<0>(pair) + (char*)"_" + std::get<1>(pair) + "_T",(TString)std::get<0>(pair) + (char*)"_" + std::get<1>(pair) + "_TP",std::get<0>(pair),std::get<1>(pair),1);

    LOpticsOptScript("y",(TString)std::get<0>(pair) + (char*)"_" + std::get<1>(pair) + "_TP",(TString)std::get<0>(pair) + (char*)"_" + std::get<1>(pair) + "_TPY",std::get<0>(pair),std::get<1>(pair),1);


    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl << endl << endl;

    i++;

  }

  

  // void LOpticsOptScript(TString select, TString SourceDataBase, TString DestDataBase, const char* min = min_def, const char* algo = al_def, Int_t save = 0)


   



}
