cppFunction(
  'IntegerVector dispersalGB(// DataFrame, NumericVector
      int sMaxPs, // dispersal
      DataFrame disperser,
      IntegerVector trick
  ){
    
    // trick def
    //int int_0 = 0;
    //int int_1 = 1;
    
    int nDisp = disperser.nrow();
    IntegerVector stepsDisp = sample(sMaxPs, nDisp, true); // number of samples is argument 2
    int maxDisp = max(stepsDisp);
    
    for(int step = 0; step<maxDisp; step++){
      IntegerVector still_disperser = nDisp;
      if(disperser.nrow()>0){
        int dp = 0;
        for(int d = 0; d<disperser.nrow(); d++){
          if(step<stepsDisp(d)){
            still_disperser(dp) = d;
            dp++;
          }
        }
        // return(step);
        // return still_disperser;
        // if(still_disperser.length()>0){
        //  stop("going within loop");
      }
      int nDispLeft = still_disperser.size(); // i.e some steps left
      if(nDispLeft>0){
        IntegerMatrix Cells_Disp(nDispLeft * 9, 3);
        IntegerVector lastDispX = disperser["lastDispX"]; // because cannot do lastDispX(still_disperser(i))
        IntegerVector lastDispY = disperser["lastDispY"];
        for(int i = 0; i<nDispLeft; i++){
          for(int l2 = 0; l2<3; l2++){
            for(int l1 = 0; l1<3; l1++){
              Cells_Disp(i * 9 + l1 + l2 * 3, 0) = lastDispX(still_disperser(i)) - 1 + l1;// -1 + l1; // because want square around indiv
              Cells_Disp(i * 9 + l1 + l2 * 3, 1) = lastDispY(still_disperser(i)) - 1 + l1;// -1 + l1; 
              Cells_Disp(i * 9 + l1 + l2 * 3, 2) = i;
            }
          }
        }
        return Cells_Disp;
      }
      
    }// end step loop
    //return disperser;
  }// end of function
 
')
  