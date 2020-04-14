void survey_check(){

  int Col = 14;
  int Row = 0;

  double exp_x, exp_y, exp_z;

  cout<<"Units in mm"<<endl;
  
  TVector3 SieveHoleTCS(SieveOffX + SieveXbyRow[Row], SieveOffY + SieveYbyCol[Col], 0);

  SieveHoleTCS.RotateX(yaw - fabs(HRSAngle));
  SieveHoleTCS.RotateY(pitch - TMath::Pi()/2);
  SieveHoleTCS.SetZ(SieveHoleTCS.Z() + ZPos);

  exp_x = SieveHoleTCS.X();
  exp_y = SieveHoleTCS.Y();
  exp_z = SieveHoleTCS.Z();
  
  /// Expected positions in TCS ///
  cout<<"Exp TCS "<<SieveHoleTCS.X()*1000<<" "<<SieveHoleTCS.Y()*1000<<" "<<SieveHoleTCS.Z()*1000<<endl;

  SieveHoleTCS.RotateX(fabs(HRSAngle));
  SieveHoleTCS.SetXYZ(SieveHoleTCS.Y(),-1*SieveHoleTCS.X(),SieveHoleTCS.Z());

  /// Expected position in HCS /////
  cout<<"Exp HCS "<<SieveHoleTCS.X()*1000<<" "<<SieveHoleTCS.Y()*1000<<" "<<SieveHoleTCS.Z()*1000<<endl;

  /// Put in Survey info for this hole ////
  TVector3 Survey(-78.5/1000,48.4/1000,790.9/1000);

  cout<<"Survey HCS "<<Survey.X()*1000<<" "<<Survey.Y()*1000<<" "<<Survey.Z()*1000<<endl;
  
  Survey.SetXYZ(-Survey.Y(),Survey.X(),Survey.Z());
  Survey.RotateX(-fabs(HRSAngle));

  /// Survey position in TCS ////
  cout<<"Survey TCS "<<Survey.X()*1000<<" "<<Survey.Y()*1000<<" "<<Survey.Z()*1000<<endl;


  cout<<"TCS SieveOffX: "<<Survey.X()*1000 - exp_x*1000<<endl;
  cout<<"TCS SieveOffY: "<<Survey.Y()*1000 - exp_y*1000<<endl;
  cout<<"TCS SieveOffZ: "<<Survey.Z()*1000 - exp_z*1000<<endl;
  
}
