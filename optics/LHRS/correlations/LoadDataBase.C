////////////////////////////////////////////////////
//          LoadDatabase
//
//   Copied from OpticOpt Script:
//   - designed to Load VDC DataBase matrix 
//    coeffecients into matrices
//    
//   
//////////////////////////////////////////////////




#include "LoadDataBase.h"



//using namespace std;
//using THaString::Split;



//void CalcMatrix( const Double_t x, vector<THaMatrixElement> matrix );


void LoadDataBase(TString DataBaseName){

  TString OldComments;


  auto Prefix = "L";
  auto fPrefix = new char [1000];
  sprintf(fPrefix, "%s", Prefix);


//   std::vector<THaMatrixElement> * fCurrentMatrixElems;
// // initial matrix elements
//   std::vector<THaMatrixElement> fTMatrixElems;
//   std::vector<THaMatrixElement> fDMatrixElems;
//   std::vector<THaMatrixElement> fPMatrixElems;
//   std::vector<THaMatrixElement> fPTAMatrixElems; // involves abs(theta_fp)
//   std::vector<THaMatrixElement> fYMatrixElems;
//   std::vector<THaMatrixElement> fYTAMatrixElems; // involves abs(theta_fp)
//   std::vector<THaMatrixElement> fFPMatrixElems;  // matrix elements used in
//   // focal plane transformations
//   // { T, Y, P }
  
//   std::vector<THaMatrixElement> fLMatrixElems;   // Path
//-length corrections (meters)
  
  
  
  OldComments = "";
  
  // Read VDC database
  
  FILE* file = fopen( DataBaseName,"r" );
  if( !file ) {
    Error("LoadDataBase","%s can not be opened", DataBaseName.Data());
									//		assert(0);//
		std::cout << "kerror" << endl;

  }
  else DEBUG_INFO("LoadDataBase","Parsing Database %s", DataBaseName.Data());
  
  // load global VDC parameters
  static const char* const here = "LoadDataBase";
  const int LEN = 200;
  char buff[LEN];
  
  //Look for the section [<prefix>.global] in the file, // e.g. [ R.global ]
  TString tag(fPrefix);
  	// Ssiz_t pos = tag.Index(".");
  	// if( pos != kNPOS )
  	// 	tag = tag(0,pos+1);
  	// else
  	// 	tag.Append(".");
  	// tag.Prepend("[");
  	// tag.Append("global]");
  	TString line, tag2(tag);
  	// tag.ToLower();

	//    TString tag(fPrefix);
    Ssiz_t tagpos = tag.Index(".");
    if (tagpos != kNPOS)
        tag = tag(0, tagpos + 1);
    else
        tag.Append(".");
    //    tag.Prepend("[");
    tag.Append("vdc.matrixelem=");
    // 	TString line, tag2(tag);
	//    TString tag2(tag);
    tag.ToLower();

	bool found = false;
	while (!found && fgets(buff, LEN, file) != NULL) {

		//read in comments
		TString tmpline = buff;
		if ( tmpline.BeginsWith("#"))
		{
			OldComments += tmpline;
//			OldComments += "\n";
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
	  Error("DB entry finding", "Database entry %s not found!", tag2.Data() );
	  fclose(file);
	  assert(0);//
	  //	  return kInitError;
	}
	
  // We found the section, now read the data

  // read in some basic constants first
  //  fscanf(file, "%lf", &fSpacing);
  // fSpacing is calculated from the actual z-positions in Init()
	fgets(buff, LEN, file); // Skip rest of line	
	fgets(buff, LEN, file); // Skip comment line

	fTMatrixElems.clear();
	fDMatrixElems.clear();
	fPMatrixElems.clear();
	fPTAMatrixElems.clear();
	fYMatrixElems.clear();
	fYTAMatrixElems.clear();
	fLMatrixElems.clear();

	fFPMatrixElems.clear();
	fFPMatrixElems.resize(3);

	typedef vector<string>::size_type vsiz_t;
	std::map<string,vsiz_t> power;
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




	// JW added to input order being optimised to
	// where order = i + j + k + l for optics elements
	std::map<string,Int_t> power_opt;
	power_opt["t"] = 4;  // transport to focal-plane tensors
	power_opt["y"] = 4;
	power_opt["p"] = 4;
	power_opt["D"] = 5;  // focal-plane to target tensors
	power_opt["T"] = 5;
	power_opt["Y"] = 5; 

	power_opt["YTA"] = 4;
	power_opt["P"] = 5;
	power_opt["PTA"] = 4;
	power_opt["L"] = 4;  

	power_opt["XF"] = 5; 
	power_opt["TF"] = 5;
	power_opt["PF"] = 5;
	power_opt["YF"] = 5;

	


	map<string,vector<THaMatrixElement>*> matrix_map;
	matrix_map["t"] = &fFPMatrixElems;
	matrix_map["y"] = &fFPMatrixElems;
	matrix_map["p"] = &fFPMatrixElems;
	matrix_map["D"] = &fDMatrixElems;
	matrix_map["T"] = &fTMatrixElems;
	matrix_map["Y"] = &fYMatrixElems;
	matrix_map["YTA"] = &fYTAMatrixElems;
	matrix_map["P"] = &fPMatrixElems;
	matrix_map["PTA"] = &fPTAMatrixElems;
	matrix_map["L"] = &fLMatrixElems;

	map <string,int> fp_map;
	fp_map["t"] = 0;
	fp_map["y"] = 1;
	fp_map["p"] = 2;

  // Read in as many of the matrix elements as there are.
  // Read in line-by-line, so as to be able to handle tensors of
  // different orders.
	while( fgets(buff, LEN, file) ) {
		string line(buff);
    // Erase trailing newline
		if( line.size() > 0 && line[line.size()-1] == '\n' ) {
			buff[line.size()-1] = 0;
			line.erase(line.size()-1,1);
		}
    // Split the line into whitespace-separated fields
		vector<string> line_spl = THaString::Split(line);

    // Stop if the line does not start with a string referring to
    // a known type of matrix element. In particular, this will
    // stop on a subsequent timestamp or configuration tag starting with "["
 if(line_spl.empty())
   continue; //ignore empty lines
 const char* w = line_spl[0].c_str();
 vsiz_t npow = power[w]; // JW line picks out 'L' or 'Y' and gets power as defined above (JW: defined above)
 if( npow == 0 )
   break;
 
 
#if DEBUG_LEVEL>=4
		cout<<"Matrix Line = ";
		for (pos=1; (UInt_t)pos<(UInt_t)line_spl.size(); pos++) {
			cout<< pos <<"("<<line_spl[pos].c_str()<<"), ";
		}
		cout<<endl;
#endif

    // Looks like a good line, go parse it.
		THaMatrixElement ME;
		ME.pw.resize(npow);
		ME.iszero = true;  ME.order = 0;
		vsiz_t pos;
		for (pos=1; pos<=npow && pos<line_spl.size(); pos++) {
			ME.pw[pos-1] = atoi(line_spl[pos].c_str());
		}
		vsiz_t p_cnt;
		for ( p_cnt=0; pos<line_spl.size() && p_cnt<kPORDER && pos<=npow+kPORDER;
					pos++,p_cnt++ )
		{
			ME.poly[p_cnt] = atof(line_spl[pos].c_str());
			if (ME.poly[p_cnt] != 0.0) {
				ME.iszero = false;
				ME.order = p_cnt+1;
			}
		}
		if (p_cnt < 1) {
		  //	Error(Here(here), "Could not read in Matrix Element %s%d%d%d!",
		  //		  w, ME.pw[0], ME.pw[1], ME.pw[2]);
		//	Error(Here(here), "Line looks like: %s",line.c_str());
			fclose(file);
			//			return kInitError;
			cout << "kiniterror" << endl;
			//			return -3;

		}


		// JW: Altered from old reading of OptOrder to setting above in this function
		
		// Olden way
		//order optimize to
		// ~~~~~~~~
		//		ME.OptOrder = atoi(line_spl[line_spl.size()-1].c_str());
		// ~~~~~~~~


		Int_t opt_order = power_opt[w];
		
		Int_t sum_powers = 0;

		for( UInt_t pc = 0; pc < npow; pc++){ //pc acronym for power count (goes through i, j and k to get size of these)
		sum_powers  += ME.pw[pc];
		}

		ME.OptOrder = opt_order - sum_powers; 
		//		cout << "opt_order = " << opt_order << ", sum_powers = " << sum_powers << ", ME.OptOrder = " << ME.OptOrder << endl;




    // Don't bother with all-zero matrix elements
		if( ME.iszero )
			continue;

    // Add this matrix element to the appropriate array
		vector<THaMatrixElement> *mat = matrix_map[w];
		if (mat) {
      // Special checks for focal plane matrix elements
			if( mat == &fFPMatrixElems ) {
				if( ME.pw[0] == 0 && ME.pw[1] == 0 && ME.pw[2] == 0 ) {
					THaMatrixElement& m = (*mat)[fp_map[w]];
					if( m.order > 0 ) {
					  //	Warning(Here(here), "Duplicate definition of focal plane ","matrix element: %s. Using first definition.", buff);
					} else
					  m = ME;
				} else
				  //				  	Warning(Here(here), "Bad coefficients of focal plane matrix ""element %s", buff);
				  cout << "warning: Bad coefficients of focal plane matrix" << endl;
				  }
			else {
	// All other matrix elements are just appended to the respective array
	// but ensure that they are defined only once!
				bool match = false;
				for( vector<THaMatrixElement>::iterator it = mat->begin();
								 it != mat->end() && !(match = it->match(ME)); it++ ) {}
				if( match ) {
				  //					Warning(Here(here), "Duplicate definition of "
				  //							"matrix element: %s. Using first definition.", buff);
				} else
					mat->push_back(ME);
			}
		}
		// else if ( fDebug > 0 )
		// 	Warning(Here(here), "Not storing matrix for: %s !", w);

	} //while(fgets)

//   // Compute derived quantities and set some hardcoded parameters
//   const Double_t degrad = TMath::Pi()/180.0;
//   fTan_vdc  = fFPMatrixElems[T000].poly[0];
//   fVDCAngle = TMath::ATan(fTan_vdc);
//   fSin_vdc  = TMath::Sin(fVDCAngle);
//   fCos_vdc  = TMath::Cos(fVDCAngle);
	//
//   // Define the VDC coordinate axes in the "detector system". By definition,
//   // the detector system is identical to the VDC origin in the Hall A HRS.
//   DefineAxes(0.0*degrad);
	//
//   fNumIter = 1;      // Number of iterations for FineTrack()
//   fErrorCutoff = 1e100;
	//
//   // figure out the track length from the origin to the s1 plane
	//
//   // since we take the VDC to be the origin of the coordinate
//   // space, this is actually pretty simple
//   const THaDetector* s1 = GetApparatus()->GetDetector("s1");
//   if(s1 == NULL)
//     fCentralDist = 0;
//   else
//     fCentralDist = s1->GetOrigin().Z();

	CalcMatrix(1.,fLMatrixElems); // tensor without explicit polynomial in x_fp


	//	fIsInit = true;
	fclose(file);
	//	return kOK;
}





void CalcMatrix( const Double_t x, vector<THaMatrixElement>& matrix )
{
  // calculates the values of the matrix elements for a given location
  // by evaluating a polynomial in x of order it->order with
  // coefficients given by it->poly


  //  cout << "value from first matrix" <<   matrix[0].v << endl;


  for( vector<THaMatrixElement>::iterator it=matrix.begin(); it!=matrix.end(); it++ ){


    it->v = 0.0;

    if(it->order > 0) {
      for(int i=it->order-1; i>=1; i--){
	// value = x * ( value + coeff)
	it->v = x * (it->v + it->poly[i]);
      }
      it->v += it->poly[0];
      //      cout << "element value = " << it->v << endl;
    }


    // 0: value = x * (0.0 + coeff(0))
    // 1: value = x * ( coeff(0)x + coeff(1) )
    // 2: value = x * ( (coeff(0)x^2 + coeff(1)x ) + coeff(2)) 
    // etc



  }

  // for(int i=0; i<=it->order-1; i++){
  //   cout << "it->poly[" << i << "] = " << it->poly[i] << endl;

    


  // }




}



void CalcMatrixElem(MatrixElems_Vals &MV)
{
  // calculates value for each event of each combination of polynomial ie for T A B C this calculates for theta^A, y^B and phi^C and calculates for all powers of X related to this as well


  // For:
  // T 2 1 0  -1.575120e+01 -5.741883e+02  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
  // values will be calculated for th^2*y^1*ph^0*x^0 and for th^2*y^1*ph^0*x^1

  MV.tg_v = 0;

  for(Int_t i = 0; i<MV.no_elements; i++){


    

    
    MV.values.at(i) =  powers[MV.DB_entries.at(i).x_p][0]*powers[MV.DB_entries.at(i).th_p][1]*powers[MV.DB_entries.at(i).y_p][2]*powers[MV.DB_entries.at(i).ph_p][3] * MV.DB_entries.at(i).DB_coeff; 


      //MV.DB_entries.at(i).x_p


    //    cout << " calculated value for element " << i << " = " << powers[MV.DB_entries.at(i).x_p][0]*powers[MV.DB_entries.at(i).th_p][1]*powers[MV.DB_entries.at(i).y_p][2]*powers[MV.DB_entries.at(i).ph_p][3] * MV.DB_entries.at(i).DB_coeff << endl;



    // cout << "testing: " << endl;
    // cout << "MV.DB_entries.at(" << i << ").x_p = " << MV.DB_entries.at(i).x_p << " with powers[MV.DB_entries.at(" << i << ").x_p][0] = " << powers[MV.DB_entries.at(i).x_p][0] <<  endl;
    // cout << "MV.DB_entries.at(" << i << ").th_p = " << MV.DB_entries.at(i).th_p << " with powers[MV.DB_entries.at(" << i << ").th_p][1] = " << powers[MV.DB_entries.at(i).th_p][1] <<  endl;
    // cout << "MV.DB_entries.at(" << i << ").y_p = " << MV.DB_entries.at(i).y_p << " with powers[MV.DB_entries.at(" << i << ").y_p][2] = " << powers[MV.DB_entries.at(i).y_p][2] <<  endl;
    // cout << "MV.DB_entries.at(" << i << ").ph_p = " << MV.DB_entries.at(i).ph_p << " with powers[MV.DB_entries.at(" << i << ").ph_p][3] = " << powers[MV.DB_entries.at(i).ph_p][3] <<  endl;


    // cout << "MV.DB_entries.at(" << i << ").DB_coeff = " << MV.DB_entries.at(i).DB_coeff << endl << endl;;

    MV.tg_v +=  MV.values.at(i);

  }
  


  

}

void MatrixElemExist(vector<THaMatrixElement>& matrix, MatrixElems_Vals &MV)
{
  // test if element exists in DB corresponding to (x^i)(theta^j)(y^k)(phi^l) 


  Int_t alt_no_entries = 0;
  
  // loop through all of one kind of element
  for( vector<THaMatrixElement>::iterator it=matrix.begin(); it!=matrix.end(); it++ ) {




    for(int i=0; i<=it->order-1; i++){
      // order here is power of X
      // it->order = (max power of x) +1
      // it->poly[i] is coeffecient of x^i
      // it->pw[j] is For j = 0 power of theta,
      //                  j = 1 power of y,
      //                  j = 2 power of phi,
      
      MV.Exist[i][it->pw[0]][it->pw[1]][it->pw[2]] = 1;


      if(MV.no_elements == 0){
	// to deal with first elements of vectors
	MV.values.at(0)  = 0;
      }
      else{
	
	MV.values.push_back(0);
      }


      DB_entry db_entry;

      db_entry.x_p = i;
      db_entry.th_p = it->pw[0];
      db_entry.y_p = it->pw[1];
      db_entry.ph_p = it->pw[2];
      db_entry.DB_coeff = it->poly[i];

      
      MV.DB_entries.push_back(db_entry);



      alt_no_entries++;

    }


  }
  
  MV.no_elements = alt_no_entries;

  


}
    



  // }












//_____________________________________________________________________________
bool THaMatrixElement::match(const THaMatrixElement& rhs) const
{
  // Compare coefficients of this matrix element to another

	if( pw.size() != rhs.pw.size() )
		return false;
	for( vector<int>::size_type i=0; i<pw.size(); i++ ) {
		if( pw[i] != rhs.pw[i] )
			return false;
	}
	return true;
}
//_____________________________________________________________________________
void THaMatrixElement::SkimPoly()
{
  //reduce order to highest non-zero poly

	if (iszero) return;

	while(!poly[order-1] && order >0)
	{
		poly.pop_back();
		order = order-1;
	}

	if (order==0) iszero = kTRUE;
}


//____________________________________
//
Double_t CalcTargetVar(const vector<THaMatrixElement>& matrix,
	      Double_t powers[][5])
{
  // calculates the value of a variable at the target
  // the x-dependence is already in the matrix, so only 1-3 (or np) used
  Double_t retval=0.0;
  Double_t v=0;
  for( vector<THaMatrixElement>::const_iterator it=matrix.begin();
       it!=matrix.end(); it++ )
    if(it->v != 0.0) {
      v = it->v;
      unsigned int np = it->pw.size(); // generalize for extra matrix elems.
      for (unsigned int i=0; i<np; i++){
	v *= powers[it->pw[i]][i+1];
	//	cout << " v = " << v << endl;
      }
      
      retval += v;
      
      //      retval += it->v * powers[it->pw[0]][1]
      //	              * powers[it->pw[1]][2]
      //	              * powers[it->pw[2]][3];
    }
    else{
      cout << "it->v == 0" << endl;
    }
  
  return retval;
}




void Read_Mat_DB(MatrixElems_Vals matrix, vector< Int_t> &DB_powers, vector< Double_t> &DB_coeffs){


   for(std::vector<DB_entry>::iterator it = matrix.DB_entries.begin(); it != matrix.DB_entries.end(); ++it){

    DB_coeffs.push_back(it->DB_coeff);
    DB_powers.push_back(it->x_p);
    DB_powers.push_back(it->th_p);
    DB_powers.push_back(it->y_p);
    DB_powers.push_back(it->ph_p);
  }
  


}


void write_Mat_DB(MatrixElems_Vals matrix, std::ofstream &DB){


   for(std::vector<DB_entry>::iterator it = matrix.DB_entries.begin(); it != matrix.DB_entries.end(); ++it){


     DB << it->x_p << " " << it->th_p << " " << it->y_p << " " << it->ph_p << " " << it->DB_coeff << endl;
     
  }
  


}
