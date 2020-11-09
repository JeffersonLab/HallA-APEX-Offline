///////////////////////////////////////////////////////////////////////////////
//
// Script designed to compare and plots two optics matrix BDs
//
//

// Author: John Williamson <jwilli@jlab.org>
// 
//  29th June 2020
//
///////////////////////////////////////////////////////////////////////////////



#include "LOpticsOpt.h"
#include "compare_DBs.h"
#include "TGraph.h"
#include "TH1F.h"

#include <typeinfo>

using namespace std;


void compare_DBs(TString DB_name_1, TString DB_name_2, Int_t EL /* element T,P,Y,D*/){


  //  TString element_string[4]; // array of strings for T,P,Y and D elements

  element_string[T] = "T";
  element_string[P] = "P";
  element_string[Y] = "Y";
  element_string[D] = "D";


  // load in first database
  LOpticsOpt* DB_1 = new LOpticsOpt();

  DB_1->LoadDataBase("DB/" + DB_name_1);

  // load in second database
  LOpticsOpt* DB_2 = new LOpticsOpt();

  DB_2->LoadDataBase("DB/" + DB_name_2);


  // read first database
  DB_MEs EL_1 = ReadDB(DB_1, EL);

  ME_Ordered coeff_names_1;
  extract_elems(EL_1,coeff_names_1);

  
  // read second database
  DB_MEs EL_2 = ReadDB(DB_2, EL); 

  ME_Ordered coeff_names_2;
  extract_elems(EL_2,coeff_names_2);

		 
  // retreive difference of MEs between DBs (and remove elements of both which are not present in the other)
  ME_Ordered diff_elems = compare_elems(coeff_names_1,coeff_names_2,0);

  // get proportional differences beween DBs
  ME_Ordered diff_elems_prop = compare_elems(coeff_names_1,coeff_names_2,1);
  

  // vector to hold graph of values for first DB
  vector<TGraph*> ME_graphs_1;		 		 
  ME_graphs_1.resize(EL_1.max_order + 1); // size to number of orders

  // function fill graphs from first DB/MEs
  fill_graphs(ME_graphs_1, coeff_names_1);


  
  // vector to hold graph of values for second DB
  vector<TGraph*> ME_graphs_2;
  ME_graphs_2.resize(EL_2.max_order + 1);
    
  fill_graphs(ME_graphs_2, coeff_names_2);
  

  // vector to hold differences between two DBs
  vector<TGraph*> ME_graphs_diff;
  ME_graphs_diff.resize(diff_elems.size());

  fill_graphs(ME_graphs_diff, diff_elems);


  // vector to hold differences between two DBs
  vector<TGraph*> ME_graphs_diff_prop;
  ME_graphs_diff_prop.resize(diff_elems_prop.size());

  fill_graphs(ME_graphs_diff_prop, diff_elems_prop);

  
  Int_t No_orders = diff_elems.size();

  // declare as many canvases are there are orders (present in both DBs)
  vector<TCanvas*> C;
  C.resize(No_orders);

  // create similar vector for multigraphs (used to draw both DD/MEs on same pad)
  vector<TMultiGraph*> multi_graphs;
  multi_graphs.resize(No_orders);

  // create legends for all orers
  vector<TLegend*> plot_legs;
  plot_legs.resize(No_orders);

  
  // loop through canvasses and fill
  for (Int_t i = 0; i<C.size() ; i++){
    cout << "in the canvas loop " << i << endl;

    C[i] = new TCanvas(Form("Order %i",i),Form("Order %i",i));
        
    C[i]->Divide(1,3);
    C[i]->cd(1);
    gPad->SetBottomMargin(0.25); 
    

    // maximum of first graph
    Double_t max_1 = TMath::MaxElement(ME_graphs_1[i]->GetN(),ME_graphs_1[i]->GetY());
    // maximum of second graph 
    Double_t max_2 = TMath::MaxElement(ME_graphs_2[i]->GetN(),ME_graphs_2[i]->GetY());
    // maximum of both graphs
    Double_t max_both  = std::max(max_1,max_2);

    // minimum of first graph
    Double_t min_1 = TMath::MinElement(ME_graphs_1[i]->GetN(),ME_graphs_1[i]->GetY());
    // minimum of second graph 
    Double_t min_2 = TMath::MinElement(ME_graphs_2[i]->GetN(),ME_graphs_2[i]->GetY());
    // minimum of both graphs
    Double_t min_both  = std::min(min_1,min_2);


    // parameter for setting margins at top and bottom of graph
    Double_t plot_margin = std::max(std::abs(max_both),std::abs(min_both));


    
    ME_graphs_1[i]->Draw("AP");
    // make sure y-axis encompasses all data points
    ME_graphs_1[i]->GetYaxis()->SetRangeUser(min_both-0.15*plot_margin, max_both+0.15*plot_margin);
    ME_graphs_1[i]->SetMarkerStyle(20);
    ME_graphs_1[i]->SetMarkerColor(kBlue);    


    ME_graphs_1[i]->GetYaxis()->SetTitle("Value of element");
    ME_graphs_1[i]->GetYaxis()->SetTitleSize(0.05);
    ME_graphs_1[i]->GetYaxis()->SetTitleOffset(0.40);
    
    // plot second DB on same subplot
    ME_graphs_2[i]->Draw("same P ");
    ME_graphs_2[i]->SetMarkerStyle(3);
    ME_graphs_2[i]->SetMarkerColor(kRed);


    // add legend to distinguish both DBs
    plot_legs[i] = new TLegend(.85,.65,0.99,.9,"Key");
    plot_legs[i]->SetFillColor(0);
    plot_legs[i]->AddEntry(ME_graphs_1[i],DB_name_1,"p");
    plot_legs[i]->AddEntry(ME_graphs_2[i],DB_name_2,"p");
    plot_legs[i]->Draw("same");

    // add gridlines
    gPad->SetGridx();
    gPad->SetGridy();

    
    C[i]->cd(2);

    // find maximum element
    Double_t max_diff_prop = TMath::MaxElement(ME_graphs_diff_prop[i]->GetN(),ME_graphs_diff_prop[i]->GetY());
    // find lowest power of 10 greater than  maximum element
    // used as upper limit as plot is in log scale
    Double_t diff_upper = powf(10.0f, ceilf(log10f(max_diff_prop)));

    
    gPad->SetBottomMargin(0.25);
    gPad->SetLogy();
    ME_graphs_diff_prop[i]->SetMaximum(diff_upper);
    ME_graphs_diff_prop[i]->Draw("AP");
    ME_graphs_diff_prop[i]->SetMarkerStyle(23);
    ME_graphs_diff_prop[i]->SetMarkerColor(kGreen+3);

    ME_graphs_diff_prop[i]->GetYaxis()->SetTitle("Proportional Difference");
    ME_graphs_diff_prop[i]->GetYaxis()->SetTitleSize(0.05);
    ME_graphs_diff_prop[i]->GetYaxis()->SetTitleOffset(0.40);

    // sets grid to be only powers of 10
    ME_graphs_diff_prop[i]->GetYaxis()->SetNdivisions(10);    
    
    gPad->SetGrid(1,1);


    C[i]->cd(3);

    Double_t max_diff_raw= TMath::MaxElement(ME_graphs_diff[i]->GetN(),ME_graphs_diff[i]->GetY());
    diff_upper = powf(10.0f, ceilf(log10f(max_diff_raw)));
    
    
    gPad->SetBottomMargin(0.25);
    gPad->SetLogy();
    ME_graphs_diff[i]->SetMaximum(diff_upper);
    ME_graphs_diff[i]->Draw("AP");
    ME_graphs_diff[i]->SetMarkerStyle(22);
    ME_graphs_diff[i]->SetMarkerColor(kMagenta+3);

    ME_graphs_diff[i]->GetYaxis()->SetTitle("Absolute Difference");
    ME_graphs_diff[i]->GetYaxis()->SetTitleSize(0.05);
    ME_graphs_diff[i]->GetYaxis()->SetTitleOffset(0.40);

    ME_graphs_diff[i]->GetYaxis()->SetNdivisions(10);
    
    gPad->SetGrid(1,1);
    
  }
  

  // loop through canvases and print
  for (Int_t i = 0; i<C.size() ; i++){

    C[i]->Print(Form("DB_comp_plots/%s_%s_%i_order_%s.pdf",DB_name_1.Data(),DB_name_2.Data(),i,element_string[EL].Data()),"pdf");

  }
  
  


  delete DB_1;
  delete DB_2;
  
}




// function designed to read in matrix elements to new struct ME for chosen DB and element type (T,P,Y or D)
DB_MEs ReadDB(LOpticsOpt* DB, Int_t variable){

  Int_t element_type; // int to hold element type
  
 
  switch(variable){
  case T:
    cout << "T elements" << endl;
    element_type = T;
    DB->fCurrentMatrixElems = &(DB->fTMatrixElems);
    break;
  case P:
    cout << "P elements" << endl;
    element_type = P;
    DB->fCurrentMatrixElems = &(DB->fPMatrixElems);
    break;
  case Y:
    cout << "Y elements" << endl;
    element_type = Y;
    DB->fCurrentMatrixElems = &(DB->fYMatrixElems);
    break;
  case D:
    cout << "D elements" << endl;
    element_type = D;
    DB->fCurrentMatrixElems = &(DB->fDMatrixElems);
    break;
  }
  

  vector<ME> ME_vect;

  
  // cout << typeid(DB->fCurrentMatrixElems).name() << endl;

  // cout << typeid(*(DB->fCurrentMatrixElems)).name() << endl;
  
  // cout << typeid((*(DB->fCurrentMatrixElems))[0]).name() << endl;

  Int_t max_order = 0; // keep track of largest order element

  for (Int_t i = 0; i < DB->fCurrentMatrixElems->size(); i++) {
    const THaMatrixElement& m = (*(DB->fCurrentMatrixElems))[i];

    // create ints to store powers of theta, Y and phi
    Int_t th_pw = m.pw[0];
    Int_t y_pw = m.pw[1];
    Int_t ph_pw = m.pw[2];
    Int_t non_x_order = th_pw + y_pw + ph_pw;
   
    
    for (Int_t j = 0; j < m.pw.size(); j++) {
      
      printf("  %2d", m.pw[j]);
    }
    for (int j = 0; j < m.order; j++) {
      ME new_el; // create new element
      new_el.type = element_type;
      
      Int_t x_pw = j;
      new_el.Order = non_x_order + x_pw;

      if(new_el.Order > max_order){
	max_order = new_el.Order;
      }
      new_el.x_pow = x_pw;
      new_el.th_pow = th_pw;
      new_el.y_pow = y_pw;
      new_el.ph_pow = ph_pw;
      new_el.Coeff = m.poly[j];
      printf("\t%g", m.poly[j]);
      ME_vect.push_back(new_el);
      
    }
    printf(" : Opt -> %d", m.OptOrder);
    if ((UInt_t) m.order != m.OptOrder) printf(" != Matrix Order !!");
    printf("\n");
  }
  

  DB_MEs DB_ME_vect;

  DB_ME_vect.MEs = ME_vect;
  DB_ME_vect.max_order = max_order;
  return DB_ME_vect;

  
}


// function that extracts vectors of co-effecients and names of elements for various orders
// produces 0th, 1st, 2nd etc order vectors

// vector of pairs should hold co-effecients and elements names for one order
// vector of vector of pairs hold the vectors for different orders
void extract_elems(DB_MEs &ME_vect, vector<vector<pair<Double_t, TString>>> &Coeff_names){

  // make as many vectors as there are orders of elements (e.g 3 if max element if 3rd order)
  // Coeff_names.reserve(ME_vect.max_order);

  Coeff_names.resize(ME_vect.max_order + 1);

  
  for (auto elem : ME_vect.MEs){
            
    for(Int_t i = 0; i<ME_vect.max_order + 1; i++){
      
      if(i == elem.Order){
	// correct order of event
	
	// inserts coeffecient of element and name (eg T1001)	
	Coeff_names[i].push_back(make_pair(elem.Coeff, element_string[elem.type] + Form("%i%i%i%i",elem.x_pow,elem.th_pow,elem.y_pow,elem.ph_pow)));
	
      }
      
    }

  }
  
}



// function designed to compare MEs from two DBs and remove elements which are not common to both
ME_Ordered compare_elems(ME_Ordered &Coeff_names_1, ME_Ordered &Coeff_names_2, Int_t flag = 0){


  if (flag !=0 && flag != 1){
    cout << "Incorrect flag for compare_elems_order function (must be 0 or 1) " << endl;
  }

 

  // loop should go through each order of matrix elements and use second function to compare MEs within orders
  ME_Ordered::iterator i = Coeff_names_1.begin();
  ME_Ordered::iterator j = Coeff_names_2.begin();

  
  // create ME_Ordered object to store differences between two DBs

  ME_Ordered diff;
  diff.resize(std::max(Coeff_names_1.size(),Coeff_names_2.size()));

  Int_t count = 0;
  
  for(i,j; i != Coeff_names_1.end() && j != Coeff_names_2.end(); ++i, ++j){

    // create new ME_names object to store differences between two DBs for one order
    ME_Names diff_order;
    
    // function that compares MEs for particular order    
    diff_order = compare_elems_order(*i,*j,flag);

    diff[count] = diff_order;

	       
    count++;
  }
    
  return diff;
}


// function to compare Matrix elements of two vectors and remove ones which are not in both
ME_Names compare_elems_order(ME_Names &Order_1, ME_Names &Order_2, Int_t flag = 0){

  if (flag !=0 && flag != 1){
    cout << "Incorrect flag for compare_elems_order function (must be 0 or 1) " << endl;
    //    exit();
  }

  
  // array to keep track of which elements in first vector have equivalents in second
  vector<Bool_t> Exist_1;
  Exist_1.resize(Order_1.size());
  std::fill(Exist_1.begin(), Exist_1.end(), false);
  Int_t count_1 = 0;
  
  vector<Bool_t> Exist_2;
  Exist_2.resize(Order_2.size());
  std::fill(Exist_2.begin(), Exist_2.end(), false);
  Int_t count_2 = 0;

  // create new ME_Names object to store differences between two DBs (for one order)
  ME_Names diff_order;

    
  for (auto elem_1 : Order_1){
    count_2 = 0;

    for (auto elem_2 : Order_2){

      // compare names of elements to check if they are the same
      if(elem_1.second == elem_2.second){
	Exist_1[count_1] = true;
	Exist_2[count_2] = true;
	if(flag == 0){
	  diff_order.push_back(make_pair(std::abs(elem_1.first-elem_2.first),elem_1.second));
	}
	else if(flag == 1){
	  diff_order.push_back(make_pair(std::abs((elem_1.first-elem_2.first)/(std::min(elem_1.first,elem_2.first))),elem_1.second));
	}
      }

      count_2++;
    }    
    count_1++;
  }

  
  for(Int_t i = 0; i<Order_1.size(); i++){
    if(!Exist_1[i]){
      // remove MEs which do not exist in other
      Order_1.erase(Order_1.begin() + i);
    }
  }

  
  for(Int_t i = 0; i<Order_2.size(); i++){
    if(!Exist_2[i]){
      // remove MEs which do not exist in other
      Order_2.erase(Order_2.begin() + i);
    }
  }


  return diff_order; 
  
}



// function used to fill graphs for each order from ME_Ordered
void fill_graphs(vector<TGraph*> &graphs, ME_Ordered coeff_names){

  
  
  // loop over all orders and create TGraphs
  for(Int_t i = 0; i < graphs.size(); i++){

    ME_Names order_names = coeff_names[i];

    graphs[i] = new TGraph(order_names.size());
    graphs[i]->GetXaxis()->SetTitle("Matrix element (x^{i}_{FP}#theta^{j}_{FP}y^{k}_{FP}#phi^{l}_{FP})");
    graphs[i]->GetXaxis()->SetTitleOffset(2.25); //1.75
    graphs[i]->GetXaxis()->SetTitleSize(0.05);
    graphs[i]->GetXaxis()->SetLabelSize(0.055);
    graphs[i]->SetTitle("");
   
    for (Int_t j = 0; j < graphs[i]->GetN(); j++){
      graphs[i]->SetPoint(j,j+1.,coeff_names[i][j].first);
    }

    for (Int_t j = 0; j < graphs[i]->GetN(); j++){
      graphs[i]->GetXaxis()->SetBinLabel(graphs[i]->GetXaxis()->FindBin(j + 1. ),coeff_names[i][j].second);
    }

  }

}


