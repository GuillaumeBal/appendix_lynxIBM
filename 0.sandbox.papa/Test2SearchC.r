require(inline)
require(readr)

Vsexe<-c("F","M","F","F","M")
Vvivant<-c(1,0,1,0,0)
Vrace=c("Linx","Loup","Mme","Elle","Lui")
sourceCpp("Search2.cpp")

matrix_out <- f_crea_matrice( Vsexe,Vvivant, Vrace )
 
str(matrix_out)
 
matrix_out
