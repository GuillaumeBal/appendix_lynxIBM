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
      nDispLeft <- still_disperser.size // i.e some steps left
      if(nDispLeft>0){
        IntegerMatrix Cells_Dip(nDispLeft*9, 3);
        for(int i = 0; i<nDispLeft, i++){
          for(int l = 0, l<3, l++){
            for(int c = 0, c<3, c++){
              IntegerMatrix(l + 0 + i * 9, 1) = disperser["lastDispX"](still_disperser(i)) -1 + l
              IntegerMatrix(l + 1 + i * 9, 2) = disperser["lastDispY"](still_disperser(i)) -1 + l
              IntegerMatrix(l + 2 + i * 9, 3) = i
            }
          }
        }
      }
      
      return still_disperser;
      
    }// end step loop
    
  }

')
  