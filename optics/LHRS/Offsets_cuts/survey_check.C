#include "InputAPEXL.h"

using namespace std;


struct hole_info
{
  Int_t col;
  Int_t row;
  Double_t surv_x;
  Double_t surv_y;
  Double_t surv_z;  
};

 
void survey_check(){

  // for LHRS sieve
  // top hole is col 14, row 8
  // bottom hole is col 6, row 12

  vector<hole_info> holes;

  holes.push_back({14,8,81.0/1000,1.4/1000,791.8/1000});
  holes.push_back({6,12,61.9/1000,-22.1/1000,793.6/1000});

  // number of holes used in check
  Int_t no_holes = holes.size();
  
  
  vector<Double_t> x_offsets;
  vector<Double_t> y_offsets;
  vector<Double_t> z_offsets;
  
  
  
  Int_t hole_count = 0;
  

  for(auto&& hole : holes)
    {

      cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl << "        For hole      " << hole_count << endl << endl;
      

      Int_t Col = hole.col;
      Int_t Row = hole.row;

      Double_t exp_x, exp_y, exp_z;

      cout<<"Units in mm"<<endl;

      TVector3 SieveHoleTCS(SieveXbyRow[Row], SieveYbyCol[Col], 0);


      SieveHoleTCS.RotateX(-(yaw - HRSAngle));
      SieveHoleTCS.RotateY(pitch - TMath::Pi()/2);
      SieveHoleTCS.SetZ(SieveHoleTCS.Z() + ZPos);


      exp_x = SieveHoleTCS.X();
      exp_y = SieveHoleTCS.Y();
      exp_z = SieveHoleTCS.Z();
  
      /// Expected positions in TCS ///
      cout<<"Exp TCS "<<SieveHoleTCS.X()*1000<<" "<<SieveHoleTCS.Y()*1000<<" "<<SieveHoleTCS.Z()*1000<<endl;
      
      SieveHoleTCS.RotateX(-HRSAngle);
      SieveHoleTCS.SetXYZ(SieveHoleTCS.Y(),-1*SieveHoleTCS.X(),SieveHoleTCS.Z());
  
      /// Expected position in HCS /////
      cout<<"Exp HCS "<<SieveHoleTCS.X()*1000<<" "<<SieveHoleTCS.Y()*1000<<" "<<SieveHoleTCS.Z()*1000<<endl;
  
      /// Put in Survey info for this hole ////
      TVector3 Survey(hole.surv_x,hole.surv_y,hole.surv_z);

      cout<<"Survey HCS "<<Survey.X()*1000<<" "<<Survey.Y()*1000<<" "<<Survey.Z()*1000<<endl;
  
      Survey.SetXYZ(-Survey.Y(),Survey.X(),Survey.Z());
      Survey.RotateX(HRSAngle);

      /// Survey position in TCS ////
      cout<<"Survey TCS "<<Survey.X()*1000<<" "<<Survey.Y()*1000<<" "<<Survey.Z()*1000<<endl;


      x_offsets.push_back(Survey.X()*1000 - exp_x*1000);
      y_offsets.push_back(Survey.Y()*1000 - exp_y*1000);
      z_offsets.push_back(Survey.Z()*1000 - exp_z*1000);
           

      cout<<"TCS SieveOffX: "<<x_offsets[hole_count]<<endl;
      cout<<"TCS SieveOffY: "<<y_offsets[hole_count]<<endl;
      cout<<"TCS SieveOffZ: "<<z_offsets[hole_count]<<endl;
      
      cout  << endl << endl;
      
      hole_count++;      
    }
  

  // get average offsets

  auto xsize = x_offsets.size();
  Double_t x_offset_av = 0;
  
  if(xsize != 0)
    {
      x_offset_av = accumulate( x_offsets.begin(), x_offsets.end(),0.0)/ xsize;
    }
  cout << "TCS SieveOffX average = " << x_offset_av << endl;
  
  auto ysize = y_offsets.size();
  Double_t y_offset_av = 0;
  
  if(ysize != 0)
    {
      y_offset_av = accumulate( y_offsets.begin(), y_offsets.end(),0.0)/ ysize;
    }
  cout << "TCS SieveOffY average = " << y_offset_av << endl;

  auto zsize = z_offsets.size();
  Double_t z_offset_av = 0;
  
  if(zsize != 0)
    {
      z_offset_av = accumulate( z_offsets.begin(), z_offsets.end(),0.0)/ zsize;
    }
  cout << "TCS SieveOffZ average = " << z_offset_av << endl;
  
  
}
