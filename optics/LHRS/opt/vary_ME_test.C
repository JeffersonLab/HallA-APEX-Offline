////////////////////////////////////////////////////
//   vary_ME_test          
//
//   Script designed to randomly vary input Matrix elements
//   for optimisation to see if they converge to the same
//   value for different initial inputs
//
//
//   John Williamson 22/11/2019
///////////////////////////////////////////////////

#include <TRandom.h>
 



Int_t LoadDataBase( TString DataBaseName, Double_t (&ME_matrix)[6][6][6][6], Char_t tg_var);



void vary_ME_test(Char_t var_tg){



  // for T martix elements (for calculating th_tg)


  //  TTree* out_T = new TTree("out_T","output tree");
 


  // randomise matrix elements and print DataBase


  TFile* newfile = TFile::Open(Form("ME_vary/rootfiles/L_%c_vary.root",var_tg),"recreate");

  TTree* out_T = new TTree("out_T","output tree");
  
  Double_t ME_initial[6][6][6][6] = {0};
  Double_t ME_final[6][6][6][6] = {0};
  


  

  //  out_T->();

  

  const Double_t up_lim = 1e+3;
  const Double_t low_lim = -1e+3;


  const  Int_t order = 5;


  TRandom rand;


  //  Double_t ME_initial[order][order][order][order] = {0};

  //.q  ME_initial[6][6][6][6] = {0};

  for(Int_t th_p = 0; th_p <= order; th_p++){
    for(Int_t y_p = 0; y_p <= order; y_p++){
      for(Int_t ph_p = 0; ph_p <= order; ph_p++){
	for(Int_t x_p = 0; x_p <= order; x_p++){

	  
	  //	  cout << " th_p = " << th_p << ", y_p = " << y_p << ", ph_p = " << ph_p << ", x_p = " << x_p  << ", with sum (th_p + y_p + ph_p + x_p) = " << (th_p + y_p + ph_p + x_p) << " (and order = " << order << ")" << endl;
	  if( (th_p + y_p + ph_p + x_p) <= order){
	    //	    cout << "Passed condition!" << endl;
	    ME_initial[th_p][y_p][ph_p][x_p] = rand.Uniform(low_lim,up_lim);
	  }
	  else{
	    //	    cout << "Failed condition :(" << endl;
	    ME_initial[th_p][y_p][ph_p][x_p] = 0;
	  }
	  

	}
      }
    }
  }
    

  
  // Save generated MEs into a DB (so that LOpticsOptScript can read it later)

  
  TString DataBaseName = "test_DB";
  
  FILE* file = fopen( DataBaseName,"w" );


  TDatime dt;
  fprintf(file,"# -------------------------------------------------------------");	fprintf(file,"\n");
  fprintf(file,"# Optimized by John Williamson @ %s",dt.AsString());fprintf(file,"\n");
  fprintf(file,"# Saved to %s",DataBaseName.Data());fprintf(file,"\n");


  fprintf(file,"[ L.global ]");fprintf(file,"\n");
  fprintf(file,"0.33l 27 1 0.0 270.2 0.0 -1.6e-03        VDC Angle, Plane Spacing, Gamma Coefficents");fprintf(file,"\n");
  fprintf(file,"L.vdc.matrixelem= ");fprintf(file,"\n");
  fprintf(file,"t 0 0 0  -1.001135e+00 -3.313373e-01 -4.290819e-02  4.470852e-03  0.000000e+00  0.000000e+00  0.000000e+00  0");fprintf(file,"\n");
  fprintf(file,"y 0 0 0  -8.060915e-03  1.071977e-03  9.019102e-04 -3.239615e-04  0.000000e+00  0.000000e+00  0.000000e+00  0");fprintf(file,"\n");
  fprintf(file,"p 0 0 0  -2.861912e-03 -2.469069e-03  8.427172e-03  2.274635e-03  0.000000e+00  0.000000e+00  0.000000e+00  0");fprintf(file,"\n");





  // print optics elements relevant to other elements 

  fprintf(file,"Y 0 0 0 1 0 0 0 0 0 0 0 0 \n");
  fprintf(file,"P 0 0 0 1 0 0 0 0 0 0 0 0 \n");
  fprintf(file,"D 0 0 0 1 0 0 0 0 0 0 0 0 \n");






  // test printing of random-prodcued matrix

 // for(Int_t th_p = 0; th_p <= order; th_p++){
 //    for(Int_t y_p = 0; y_p <= order; y_p++){
 //      for(Int_t ph_p = 0; ph_p <= order; ph_p++){
 // 	for(Int_t x_p = 0; x_p <= order; x_p++){

 // 	  if(x_p == 0){
	    
 // 	    //	    cout << " T " << th_p << " " << y_p << " " << ph_p << " ~~~~~~~ " << ME_initial[th_p][y_p][ph_p][x_p];
 // 	  }
 // 	  else{
 // 	    //	    cout << " " <<  ME_initial[th_p][y_p][ph_p][x_p];
 // 	  }	  
 // 	}
 // 	cout << endl;
 //      }
 //    }
 // }





   for(Int_t th_p = 0; th_p <= order; th_p++){
    for(Int_t y_p = 0; y_p <= order; y_p++){
      for(Int_t ph_p = 0; ph_p <= order; ph_p++){
	//	for(Int_t x_p = 0; x_p <order; x_p++){
	

	Bool_t newline = false;

	Int_t x_p = 0;

	
	while(  ME_initial[th_p][y_p][ph_p][x_p] != 0 && (th_p + y_p + ph_p + x_p) <= order){
	  
	  newline = true;

	  if( x_p == 0){
	    
	    fprintf(file,"T %i %i %i %13.6e", th_p, y_p, ph_p, ME_initial[th_p][y_p][ph_p][x_p]);
	  }
	  else{
	    
	    fprintf(file," %13.6e",ME_initial[th_p][y_p][ph_p][x_p]);
	  }
	  
	  
	  x_p++;
	  
	}


	Int_t No_of_zeros  = order+1 - (x_p) ;


	//	cout << "order = " << order << ", (x_p) = " << (x_p) << ",  No_of_zeros = " << No_of_zeros << endl;
	
	// for(Int_t No_zeros = 0; No_zeros < No_of_zeros; No_zeros++){

	//   fprintf(file," 0 ");
	  
	// }

	if (newline){

	  for(Int_t No_zeros = 0; No_zeros < No_of_zeros; No_zeros++){

	    fprintf(file," 0 ");
	  
	  }


	  fprintf(file," \n");


	}

	  
	  
	  // if(  ME_initial[th_p][y_p][ph_p][x_p] != 0){
	  //   if( x_p == 0){
	      
	  //     fprintf(file,"%i %i %i %13.6e", th_p, y_p, ph_p, ME_initial[th_p][y_p][ph_p][x_p]);
	  //   }
	  //   else if( x_p > 0){
	      

	  //   }
	  


      }
	  
    }
   }




   // next stage: use this DB with the optics opt script 

   fclose(file);

   
   LOpticsOptScript("theta","test_DB","test_DB_2");
   //   LOpticsOptScript("theta","opt_output/Vars_with_cuts_4179_1_7_2019_dp_opt.dat","test_DB_2");





   // after this: read output DB into matrix and tree as 'after opt'



   //   ME_final[6][6][6][6] = {1};

   LoadDataBase("test_DB_2", ME_final, 'T');

   // might want to save chi^2 (before and after somehow)




   // test of printing of 2nd DB

   
   // cout << "Printing final DB " << endl;
   // cout << " ~~~~~~~~~~~~~~ " << endl << endl;
   
   // for(Int_t th_p = 0; th_p <= order; th_p++){
   //   for(Int_t y_p = 0; y_p <= order; y_p++){
   //     for(Int_t ph_p = 0; ph_p <= order; ph_p++){
   // 	 for(Int_t x_p = 0; x_p <= order; x_p++){
	   
   // 	   if(x_p == 0){
	     
   // 	     cout << " T " << th_p << " " << y_p << " " << ph_p << " ~~~~~~~ " << ME_final[th_p][y_p][ph_p][x_p];
   // 	   }
   // 	   else{
   // 	     cout << " " <<  ME_final[th_p][y_p][ph_p][x_p];
   // 	   }	  
   // 	 }
   // 	 cout << endl;
   //     }
   //   }
   // }
   



   



   

}






Int_t LoadDataBase( TString DataBaseName_new, Double_t (&ME_matrix)[6][6][6][6], Char_t tg_var)
{


  cout << "Entered LoadDataBase" << endl;


  // Read VDC database

  TString OldComments;


  FILE* new_file = fopen( DataBaseName_new,"r" );

  if( !new_file ) {
    cout << "File not load" << endl;      
  }

  //load global VDC parameters

  const int LEN = 200;
  char buff[LEN];


  auto fPrefix = new char [1000];

  fPrefix = "L";

  TString tag(fPrefix);
  Ssiz_t tagpos = tag.Index(".");
  if (tagpos != kNPOS)
    tag = tag(0, tagpos + 1);
  else
        tag.Append(".");
    //    tag.Prepend("[");
  tag.Append("vdc.matrixelem=");
  TString line, tag2(tag);
  //  TString tag2(tag);
  tag.ToLower();

  bool found = false;
  while (!found && fgets(buff, LEN, new_file) != NULL) {

    //read in comments
    TString tmpline = buff;
    if ( tmpline.BeginsWith("#"))
      {
	OldComments += tmpline;
	OldComments += "\n";
      }
    
    char* buf = ::Compress(buff);  //strip blanks
    line = buf;
    delete [] buf;
    if( line.EndsWith("\n") ) line.Chop();
    
    line.ToLower();
    if ( tag == line )
      found = true;
      }
  if( !found ) {
    //      Error(Here(here), "Database entry %s not found!", tag2.Data() );
    fclose(new_file);
    assert(0);//
    return -1;
  }
  
    
    
    // We found the section, now read the data

    // read in some basic constants first
    //  fscanf(file, "%lf", &fSpacing);
    // fSpacing is calculated from the actual z-positions in Init()
  fgets(buff, LEN, new_file); // Skip rest of line	
  fgets(buff, LEN, new_file); // Skip comment line
  

  typedef vector<string>::size_type vsiz_t;
  map<string,vsiz_t> power;
  power["t"] = 3;  // transport to focal-plane tensors
  power["y"] = 3;
  power["p"] = 3;
  power["D"] = 3;  // focal-plane to target tensors
  power["T"] = 3;
  power["Y"] = 3; // (JW: Defined above references)
  power["YTA"] = 4;
  power["P"] = 3;
  power["PTA"] = 4;
  power["L"] = 4;  // pathlength from z=0 (target) to focal plane (meters) (JW: Defined above references)
  power["XF"] = 5; // forward: target to focal-plane (I think)
  power["TF"] = 5;
  power["PF"] = 5;
  power["YF"] = 5;



  // define integers to hold powers of MEs
  Int_t th_power = 0;
  Int_t y_power = 0;
  Int_t ph_power = 0;
  Int_t x_power = 0;

	


  // Read in as many of the matrix elements as there are.
  // Read in line-by-line, so as to be able to handle tensors of
  // different orders.
  while( fgets(buff, LEN, new_file) ) {



    string line(buff);
    // Erase trailing newline
    if( line.size() > 0 && line[line.size()-1] == '\n' ) {
      buff[line.size()-1] = 0;
      line.erase(line.size()-1,1);
    }
    // Split the line into whitespace-separated fields
    vector<string> line_spl = Split(line);
    



    // Stop if the line does not start with a string referring to
    // a known type of matrix element. In particular, this will
    // stop on a subsequent timestamp or configuration tag starting with "["
    if(line_spl.empty())
      continue; //ignore empty lines
    const char* w = line_spl[0].c_str();



    if( *w == tg_var){ // test for 'correct' variety of ME (ie t if theta is being optimised
 
      vsiz_t npow = power[w];
    
  


    // Looks like a good line, go parse it.
      


	// use this to select powes of theta, y and phi
      th_power = atoi(line_spl[1].c_str());
      y_power = atoi(line_spl[2].c_str());
      ph_power = atoi(line_spl[3].c_str());
      
      vsiz_t pos = 0;

      vsiz_t p_cnt = 4; // first line entry after element and power info ie T 0 0 0 A B C etc.. (5th entry in line is first ME)

      while( atoi(line_spl[p_cnt].c_str()) != 0 && p_cnt<line_spl.size()){
	

	  	  
	ME_matrix[th_power][y_power][ph_power][pos] = atof(line_spl[p_cnt].c_str());
	  
	  pos++;
	  p_cnt++;

      }



      // for ( p_cnt=4; pos<line_spl.size();pos++,p_cnt++ )
      // 	{
	  
      // 	  ME_matrix[th_power][y_power][ph_power][pos] = atof(line_spl[p_cnt].c_str());
      
    



    } // if loop for correct variable end
    
  
  }// while loop end


  
    // vsiz_t p_cnt;
    // for ( p_cnt=0; pos<line_spl.size() && p_cnt<kPORDER && pos<=npow+kPORDER;
    // 	  pos++,p_cnt++ )
    //   {
    // 	ME.poly[p_cnt] = atof(line_spl[pos].c_str());
    // 			if (ME.poly[p_cnt] != 0.0) {
    // 			  ME.iszero = false;
    // 			  ME.order = p_cnt+1;
    // 			}
    //   }
    // if (p_cnt < 1) {
    //   Error(Here(here), "Could not read in Matrix Element %s%d%d%d!",
    // 	    w, ME.pw[0], ME.pw[1], ME.pw[2]);
    //   Error(Here(here), "Line looks like: %s",line.c_str());
    //   fclose(new_file);
    //   return kInitError;
    // }
    
    
// 		// JW: Altered from old reading of OptOrder to setting above in this function
		
// 		// Olden way
// 		//order optimize to
// 		// ~~~~~~~~
// 		//		ME.OptOrder = atoi(line_spl[line_spl.size()-1].c_str());
// 		// ~~~~~~~~


// 		Int_t opt_order = power_opt[w];
		
// 		Int_t sum_powers = 0;

// 		for( Int_t pc = 0; pc < npow; pc++){ //pc acronym for power count (goes through i, j and k to get size of these)
// 		sum_powers  += ME.pw[pc];
// 		}

// 		ME.OptOrder = opt_order - sum_powers; 
// 		cout << "Type: " << w << " :  opt_order = " << opt_order << ", sum_powers = " << sum_powers << ", ME.OptOrder = " << ME.OptOrder << endl;




//     // Don't bother with all-zero matrix elements
// 		if( ME.iszero )
// 			continue;

//     // Add this matrix element to the appropriate array
// 		vector<THaMatrixElement> *mat = matrix_map[w];
// 		if (mat) {
//       // Special checks for focal plane matrix elements
// 			if( mat == &fFPMatrixElems ) {
// 				if( ME.pw[0] == 0 && ME.pw[1] == 0 && ME.pw[2] == 0 ) {
// 					THaMatrixElement& m = (*mat)[fp_map[w]];
// 					if( m.order > 0 ) {
// 						Warning(Here(here), "Duplicate definition of focal plane "
// 								"matrix element: %s. Using first definition.", buff);
// 					} else
// 						m = ME;
// 				} else
// 					Warning(Here(here), "Bad coefficients of focal plane matrix "
// 							"element %s", buff);
// 			}
// 			else {
// 	// All other matrix elements are just appended to the respective array
// 	// but ensure that they are defined only once!
// 				bool match = false;
// 				for( vector<THaMatrixElement>::iterator it = mat->begin();
// 								 it != mat->end() && !(match = it->match(ME)); it++ ) {}
// 				if( match ) {
// 					Warning(Here(here), "Duplicate definition of "
// 							"matrix element: %s. Using first definition.", buff);
// 				} else
// 					mat->push_back(ME);
// 			}
// 		}
// 		else if ( fDebug > 0 )
// 			Warning(Here(here), "Not storing matrix for: %s !", w);

// 	} //while(fgets)

// //   // Compute derived quantities and set some hardcoded parameters
// //   const Double_t degrad = TMath::Pi()/180.0;
// //   fTan_vdc  = fFPMatrixElems[T000].poly[0];
// //   fVDCAngle = TMath::ATan(fTan_vdc);
// //   fSin_vdc  = TMath::Sin(fVDCAngle);
// //   fCos_vdc  = TMath::Cos(fVDCAngle);
// 	//
// //   // Define the VDC coordinate axes in the "detector system". By definition,
// //   // the detector system is identical to the VDC origin in the Hall A HRS.
// //   DefineAxes(0.0*degrad);
// 	//
// //   fNumIter = 1;      // Number of iterations for FineTrack()
// //   fErrorCutoff = 1e100;
// 	//
// //   // figure out the track length from the origin to the s1 plane
// 	//
// //   // since we take the VDC to be the origin of the coordinate
// //   // space, this is actually pretty simple
// //   const THaDetector* s1 = GetApparatus()->GetDetector("s1");
// //   if(s1 == NULL)
// //     fCentralDist = 0;
// //   else
// //     fCentralDist = s1->GetOrigin().Z();

// 	CalcMatrix(1.,fLMatrixElems); // tensor without explicit polynomial in x_fp


// 	fIsInit = true;
// 	fclose(new_file);
// 	return kOK;
  return 0;
}
  
