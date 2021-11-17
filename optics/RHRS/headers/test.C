void test(){

  double angles[2];
  double xy[2];

  for(int i = 0; i<27;i++){

    Sieve_hole_pos(1,i,5,angles,xy);

    cout<<i<<" "<<angles[0]*1000<<endl;


  }


}
