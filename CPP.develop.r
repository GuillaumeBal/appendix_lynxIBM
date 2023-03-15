require(inline)

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
