cppFunction(
  'IntegerVector dispersalGB(// DataFrame, NumericVector
      DataFrame lynx,
      int sMaxPs, // dispersal
      IntegerMatrix HabitatMap,
      double pMat
  ){
    
    ///////////////////////////////////////////////////////////////////////////////////////
    // TRICKS AND FUNCTIONS FOR CODE
    
    // some integer def because IntegerVector v(1) = 1 does not work
    int int_0 = 0;
    int int_1 = 1;
    
    // function check values CharacterVector equal to a criteria
    int EqualString(
        CharacterVector ToCheck, 
        char Crit
    ){
      int n_equal = 0;
      for(in i=0; i<CharacterVec.size();i++){
        if (ToCheck(i) == Crit){
          n_equal++
        }
      }
      return n_equal;
    }
      //     return true;
      //   else
      //     return false;
      // }
      
      ///////////////////////////////////////////////////////////////////////////////////////
      // LYNX DATA PROCESSING
      
      // attach all lynx data
      NumericVector lynx_xcor = lynx["xcor"], lynx_ycor = lynx["ycor"], lynx_who = lynx["who"], lynx_heading = lynx["heading"], lynx_prevX = lynx["prevX"], lynx_prevY = lynx["prevY"], lynx_age = lynx["age"], lynx_lastDispX = lynx["lastDispX"], lynx_lastDispY = lynx["lastDispY"], lynx_nMat = lynx["nMat"], lynx_maleID = lynx["maleID"], lynx_nFem = lynx["nFem"], lynx_rdMortTerr = lynx["rdMortTerr"]; 
    CharacterVector lynx_breed = lynx["breed"], lynx_color = lynx["color"], lynx_pop = lynx["pop"], lynx_sex = lynx["sex"], lynx_status = lynx["status"];
    
    // create disperser and non disperser var
    int nLynx = lynx_xcor.size();
    int nDisp = lynx_status;
    IntegerVector Dispersers_xcor;
    IntegerVector Dispersers_ycor;
    for(int i = 0, i_disp = 0, i_ndisp = 0 ; i<lynx_status.size(); i++){
      if(lynx_status(i) == "disp"){
        Dispersers_xcor(i_disp) = lynx_xcor(i), Dispersers_ycor(i_disp) = lynx_ycor(i); //Dispersers_who(i_disp) = lynx_who(i), Dispersers_heading(i_disp) = lynx_heading(i), Dispersers_prevX(i_disp) = lynx_prevX(i), Dispersers_prevY(i_disp) = lynx_prevY(i), Dispersers_breed(i_disp) = lynx_breed(i), Dispersers_color(i_disp) = lynx_color(i), Dispersers_pop(i_disp) = lynx_pop(i), Dispersers_sex(i_disp) = lynx_sex(i), Dispersers_age(i_disp) = lynx_age(i), Dispersers_status(i_disp) = lynx_status(i), Dispersers_lastDispX(i_disp) = lynx_lastDispX(i), Dispersers_lastDispY(i_disp) = lynx_lastDispY(i), Dispersers_nMat(i_disp) = lynx_nMat(i), Dispersers_maleID(i_disp) = lynx_maleID(i), Dispersers_nFem(i_disp) = lynx_nFem(i), Dispersers_rdMortTerr(i_disp) = lynx_rdMortTerr(i);
        i_disp++;
        //nDisp = i_disp;
      }
      else{
        i_ndisp ++;
      }
    }
    return Dispersers_xcor;
    
    
    // attach all data.Frame vectors within cpp because dataframe not supported for computations
    //NumericVector disperser_xcor = disperser["xcor"], disperser_ycor = disperser["ycor"], disperser_who = disperser["who"], disperser_heading = disperser["heading"], disperser_prevX = disperser["prevX"], disperser_prevY = disperser["prevY"], disperser_age = disperser["age"], disperser_lastDispX = disperser["lastDispX"], disperser_lastDispY = disperser["lastDispY"], disperser_nMat = disperser["nMat"], disperser_maleID = disperser["maleID"], disperser_nFem = disperser["nFem"], disperser_rdMortTerr = disperser["rdMortTerr"]; 
    //CharacterVector disperser_breed = disperser["breed"], disperser_color = disperser["color"], disperser_pop = disperser["pop"], disperser_sex = disperser["sex"], disperser_status = disperser["status"];
    
    // // trick def
    // int int_0 = 0;
    // int int_1 = 1;
    // 
    // int nDisp = disperser.nrow();
    // IntegerVector stepsDisp = sample(sMaxPs, nDisp, true); // number of samples is argument 2
    // int maxDisp = max(stepsDisp);
    // 
    // for(int step = 0; step<maxDisp; step++){
    //   IntegerVector still_disperser = nDisp;
    //   if(disperser.nrow()>0){
    //     //int dp = 0;
    //     for(int d = 0, dp = 0; d<disperser.nrow(); d++){
    //       if(step<stepsDisp(d)){
    //         still_disperser(dp) = d;
    //         dp++;
    //       }
    //     }
    //   }
    //   //return still_disperser; //works up to here without crashing after a while
    //   int nDispLeft = still_disperser.size(); // i.e are there some steps left for the indiv ?
    //   if(nDispLeft>0){
    //     // mat with X / Y / id ind / Habitat cell
    //     IntegerMatrix Cells_Disp(nDispLeft * 9, 4);
    //     IntegerVector lastDispX = disperser["lastDispX"]; // because cannot do lastDispX(still_disperser(i))
    //     IntegerVector lastDispY = disperser["lastDispY"];
    //     // hab freq indiv
    //     IntegerMatrix HabFreq_Disp(nDispLeft, 3);
    //     // vector recording whether next step is in matrix
    //     IntegerVector Mat_chosen = nDispLeft;
    //     for(int i = 0; i<nDispLeft; i++){
    //       for(int l2 = 0; l2<3; l2++){
    //         for(int l1 = 0; l1<3; l1++){
    //           Cells_Disp(i * 9 + l1 + l2 * 3, 0) = lastDispX(still_disperser(i)) - 1 + l1;// because want square around indiv
    //           Cells_Disp(i * 9 + l1 + l2 * 3, 1) = lastDispY(still_disperser(i)) - 1 + l1;
    //           Cells_Disp(i * 9 + l1 + l2 * 3, 2) = i;
    //           Cells_Disp(i * 9 + l1 + l2 * 3, 3) = HabitatMap(Cells_Disp(i * 9 + l1 + l2 * 3, 0), Cells_Disp(i * 9 + l1 + l2 * 3, 1));
    //           // here count habitat occurences
    //           if(Cells_Disp(i * 9 + l1 + l2 * 3, 3) == 0){ // barrier
    //             HabFreq_Disp(i, 0)++;
    //           }
    //           if(Cells_Disp(i * 9 + l1 + l2 * 3, 3) == 2){ // matrix
    //             HabFreq_Disp(i, 1)++;
    //           }
    //           if(Cells_Disp(i * 9 + l1 + l2 * 3, 3) == 3){ // dispersing
    //             HabFreq_Disp(i, 2)++;
    //           }
    //           if(Cells_Disp(i * 9 + l1 + l2 * 3, 3) == 4){ // breeding
    //             HabFreq_Disp(i, 2)++;
    //           }
    //         }
    //       }
    //       // does indiv pick matrix type habitat for next move (pMat * freq hab 2)
    //       //Mat_chosen(i) = R::rbinom(int_1, (pMat / 100) * HabFreq_Disp(i, 1)); // number of samples is argument 2  
    //       
    //     }// end i loop indiv
    
    //return Cells_Disp; //nextCellsType;
    
    //     // new ind loop for matrix dip kept because integrating with previous was messy
    //     // for(int i =0, p = 0; i<nDispLeft; i++){
    //     //   // cell left to disperse
    //     //   IntegerMatrix nextCellsType;
    //     //   for(int nci = 0; nci<9; nci++){// nci number cell per indiv, fill final dispersal matrix base on choice
    //     //     for(int c = 0; c<3; c++){ // ncolumn in final table to fill
    //     //       if((Mat_chosen[i] == 1) & (Cells_Disp(nci + i * nci, 3) == 2)){
    //     //         //nextCellsType(p, c) = 1; // Cells_Disp(nci + i * nci, c);
    //     //         p++;
    //     //       }else{
    //     //         if((Mat_chosen[i] == 0) & (Cells_Disp(nci + i * nci, 3) == 3)){
    //     //           //stop("got within");
    //     //           //nextCellsType(p, c) = 1;// Cells_Disp(nci + i * nci, c);
    //     //           p++;
    //     //         }else{
    //     //           if((Mat_chosen[i] == 0) & (Cells_Disp(nci + i * nci, 3) == 4)){
    //     //             //   //stop("got within");
    //     //             //   //nextCellsType(p, c) = 1;// Cells_Disp(nci + i * nci, c);
    //     //             p++;
    //     //           }
    //     //         }
    //     //       }
    //     //     }
    //     //   }   
    //     // } // end second ind loop
    //     return HabFreq_Disp;
    // }// end if still some dips
    //}// end step loop
  }// end of function

')
  