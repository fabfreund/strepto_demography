#' Replicated from B. Eldon's work:
## Copyright (2014):  Bjarki Eldon 
## Distributed under GNU GPL v3 or later
## compile  .c file in terminal as: #c file in this folder, same copyright
## R CMD SHLIB -O3 CREBiepg.c -lm -lgsl -lgslcblas 
## code based on Polanski and Kimmel (2003)
## and Polanski, Bobrowski, Kimmel (2003)
## Nl is number leaves = sample size
## bp is the growth parameter (b)
## returns expected branch lengths E_b[B_i] where B_i
## is the random length of branches subtending i leaves

"REBiepg" <- function( Nl, bp )
{
  if( is.loaded( "CREBiepg.so" ) ){
    dyn.unload( "CREBiepg.so")
    dyn.load( "CREBiepg.so") }
  else
    dyn.load("CREBiepg.so")
  
  vx <- as.double( (1:(Nl)) * 0 )
  
  out =  .C("EBi_W", Nleaves = as.integer(Nl), betaparam = as.double(bp),   x = vx )
  
  return( out$x[-1] )
  
}
