require(inline)

#agents.is.agentmatrix <- inherits(agents, "agentMatrix") * 1
agents.is.agentmatrix = TRUE
#agents2.is.agentmatrix <- inherits(agents2, "agentMatrix") * 1
agents2.is.agentmatrix = TRUE
Matrice1 <- matrix(1:4, ncol = 2)
Matrice2 <- matrix(5:8, ncol = 2)

toward.cpp <-
  cfunction(
    signature(
      agents= "matrix",
      agents2= "matrix",
      agents_is_agentmatrix= "integer",
      agents2_is_agentmatrix= "integer"),
    body = "
#include <R.h>
//#include <Rmath.h>
//#include <Rinternals.h>
   
 int n=ncols(agents2);
 agents_is_agentmatrix=coerceVector(agents_is_agentmatrix,INTSXP);
 agents2_is_agentmatrix=coerceVector(agents_is_agentmatrix,INTSXP);
 //SEXP out = PROTECT(allocInterger(INTSXP));
 SEXP out = PROTECT(allocVector(INTSXP,1));
 if(agents_is_agentmatrix==agents2_is_agentmatrix){
   INTEGER(out)[0]=50;
 }
 else{
   INTEGER(out)[0]=100;
 }
 UNPROTECT(1);
 return out;
"
  )
toward.cpp(agents = Matrice1, agents2 = Matrice2, 
           agents_is_agentmatrix = 1, agents2_is_agentmatrix = 2)

#require(inline)
seqn = cfunction(signature(n_="integer", start_="integer"), 
                 body="
 int i, n = asInteger(n_), start=asInteger(start_);
 SEXP out = PROTECT(allocVector(INTSXP, n));
 for(i=0; i<n; ++i){
   INTEGER(out)[i]=i+start;
 }
 UNPROTECT(1);
 return out;
")
seqn(8, 3)










