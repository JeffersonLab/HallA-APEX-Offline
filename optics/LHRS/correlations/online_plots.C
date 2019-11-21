#include "corr_select.C"



void online_plots(TString DB_name, Int_t runnum){


  gSystem->Exec("mkdir /home/johnw/public_html/correlation_plots/" + DB_name);



  
  corr_select(DB_name, runnum,-1,-1);

  //  columns
  for(Int_t i = 0; i < 27; i++){
    corr_select(DB_name, runnum,i,-1);
  }

  // rows
  for(Int_t i = 0; i < 17; i++){
    corr_select(DB_name, runnum,-1,i);
  }


  // holes
  for(Int_t i = 0; i < 17; i++){
    for(Int_t j = 0; j < 27; j++){
    
      corr_select(DB_name, runnum,j,i);
      
    }
  }


  gSystem->Exec("mkdir /home/johnw/public_html/correlation_plots/" + DB_name);





  gSystem->Exec("mkdir /home/johnw/public_html/correlation_plots/" + DB_name + "/hole_cuts_used");


  TString optics_dir = "/home/johnw/HallA_scripts/APEX_scripts/optics";

  gSystem->Exec("cp " + optics_dir + "/Hole_selection_th_ph_4179.png " + optics_dir + "/Hole_selection_FP_4179.png " + optics_dir + "/Hole_selection_th_ph_4181.png " + optics_dir + "/Hole_selection_FP_4181.png" + " /home/johnw/public_html/correlation_plots/" + DB_name + "/hole_cuts_used");


  gSystem->Exec("mkdir /home/johnw/public_html/correlation_plots/" + DB_name + "/opt_output");


  TString opt_output = "/home/johnw/HallA_scripts/APEX_scripts/optics/opt_output/";

		gSystem->Exec("cp " + opt_output + DB_name + ".Sieve.Opt.png " + opt_output + DB_name + ".Vertex.Opt.png " + opt_output + DB_name + ".VertexRes.Opt.png" + " /home/johnw/public_html/correlation_plots/" + DB_name + "/opt_output/");


  





}