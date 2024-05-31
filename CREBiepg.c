#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h> 
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_elementary.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_expint.h>


double d_aj( double beta, int j )
{
  // return a_j computed using double
  double bjb = gsl_sf_choose( j,2)/beta ;

  if( bjb < 50 ){
    return ( gsl_sf_exp( bjb ) * gsl_sf_expint_E1( bjb) / ( (double)beta ) ) ; }
  else{
    // use approximation to E_1( bjb )
    int n;
    double s = 0.0 ;
    for( n = 0 ; n < 30 ; n++ ){
      s =  s  +  ( gsl_sf_fact( n ) / gsl_sf_pow_int( -bjb, n) ) ; }
    return( s/( bjb * beta ) ) ; }
}




void EBi_W( int * Nleaves,  double * betaparam, double * x )
{

  double beta = betaparam[0] ;
  int N = Nleaves[0] ;
  
  gsl_vector * Wj = gsl_vector_calloc( N + 2 ) ;
  int i, j ;
  double dN = (double)N ;
  double di, dj ;
  for( i = 1 ; i < N ; i++ ){
    di = (double)i ;

    gsl_vector_set_zero( Wj ) ;

    gsl_vector_set( Wj, 2, 6.0/(dN + 1) ) ;
    gsl_vector_set( Wj, 3, 30.0*(dN -  2*di)/( (dN+1)*(dN+2) ) ) ;

    for( j = 2; j <= N-2 ; j++ ){
      dj = (double)j ;
      gsl_vector_set(  Wj,j+2, (3+  2*dj)*(dN -  2*di)*gsl_vector_get(Wj,j+1)/( dj*(dN + dj+1) )  -   ( (1+dj)*(3+  2*dj)*(dN -  dj)*gsl_vector_get(Wj,j)/( dj*(2*dj  - 1)*(dN +  dj+1) ) ) ) ; }

    for( j = 2 ; j <= N ; j++ ){
      // double d_aj( double beta, int j )
      x[i] =   x[i]   +   (gsl_vector_get(Wj,j) * d_aj( beta, j) ) ; } }

  gsl_vector_free ( Wj ) ;

  //for( j = 1 ; j < N ; j++ ){
  //  printf("%.10g\n", gsl_vector_get( x,j) ) ; }
}


