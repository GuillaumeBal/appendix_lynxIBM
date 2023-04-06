cppFunction(
  'dispersalGB( // DataFrame, NumericVector
    int sMaxPs, // dispersal
    DataFrame disperser,
    IntegerVector trick
) {
  
  // trick def
  //int int_0 = 0;
  //int int_1 = 1;
  
  IntegerVector vec_ref_size = disperser["who"];
  int nDisp=disperser.size();
  IntegerVector stepsDisp = sample(sMaxPs, nDisp, true); // number of samples is argument 2
  int maxDisp = max(stepsDisp);

  // for(int step = 0; step<maxDisp; step++){
  //   IntegerVector still_disperser;
  //   if(disperser.nrow()>0){
  //     //int dp = 0;
  //     for(int d = 0; d<disperser.nrow(); d++){
  //       if(step<stepsDisp(d)){
  //         still_disperser(d) = 1;
  //       }else{
  //         still_disperser(d) = 0;
  //       }
  //     }
  //   }
  //   // return(step);
  //   // return still_disperser;
  //   // if(still_disperser.length()>0){
  //   //  stop("going within loop");
  //   // }
  // }
  
  return nDisp;
  
}
')'
