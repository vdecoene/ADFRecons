#include "port_i.h"
#include <stdio.h>
#define N_F_DIM 10

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void a_func( int& N, SCALAR* X, int& NF, SCALAR* F, int* UIPARM, SCALAR* URPARM, void* UFPARM )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Function to minimise is: F( x1, ..., xN ) = sum{ ( xi - i )^2, i }
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  F[ 0 ] = 0.;
  for ( int i = 0; i < N_F_DIM; i++ ) 
    F[ 0 ]+= ( X[ i ] - (SCALAR)i )*( X[ i ] - (SCALAR)i );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void a_grad( int& N, SCALAR* X, int& NF, SCALAR* G, int* UIPARM, SCALAR* URPARM, void* UFPARM )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Analytical Gradient of the function to minimise.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{   
  for ( int i = 0; i < N_F_DIM; i++ ) 
    G[ i ]= 2.*( X[ i ] - (SCALAR)i );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void a_hess( int& N, SCALAR* X, int& NF, SCALAR* G, SCALAR* H, int* UIPARM, SCALAR* URPARM, void* UFPARM )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Analytical Gradient and Hessian of the function to minimise.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  for ( int i = 0; i < N_F_DIM; i++ ) 
    G[ i ]= 2.*( X[ i ] - (SCALAR)i );
  
  int k = 0;
  for ( int i = 0; i < N_F_DIM; i++ ) 
    for ( int j = 0; j <= i; j++ ) 
    {
      if ( i == j ) 
        H[ k ] = 2.;
      else 
        H[ k ] = 0.;
      k++; 
    }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void print_X( std::string str, PORT_i* d )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  std::cout << str << ":"<< std::endl;
  for ( int i = 1; i <= N_F_DIM; i++ ) 
    printf( "X(%2d) = %10.3lf\n", i, d->X( i ) );
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main()
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  //                      Part I: unbounded fit
  //
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // Get a DMN toolbox
  //===
  PORT_i* x_PORT_i = new PORT_i( N_F_DIM );
  
  // Steer DMN algorithms 'a la FORTRAN'
  //===
  x_PORT_i->IV( 17 ) = 1000;  // Max func. evaluations
  x_PORT_i->IV( 18 ) =  100;  // Max iterations
  x_PORT_i->IV( 21 ) =    0;  // Mute stdio printout
  
  // Function to minimize and its Gradient, Hessian
  //===
  x_PORT_i->f   = &a_func;
  x_PORT_i->g   = &a_grad;  // Gradient
  x_PORT_i->gh  = &a_hess;  // Gradient and Hessian
  
  // Initial guess
  //===
  std::vector<SCALAR> x0( N_F_DIM, 0.0 );
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Solve with DMNF
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  
  // Initial Guess
  //===
  for( int i = 1; i <= N_F_DIM; i++ ) 
    x_PORT_i->X( i ) = 10.0;
  
  // Minimisation
  //===
  if ( x_PORT_i->MNF() ) 
    print_X( "DMNF SOLUTION ", x_PORT_i );
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Solve with DMNG: using analytical Gradient
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
  // Initial Guess
  //===
  for( int i = 1; i <= N_F_DIM; i++ ) 
    x_PORT_i->X( i ) = 10.0;
  
  // Minimisation
  //===
  if ( x_PORT_i->MNG() ) 
    print_X( "DMNG SOLUTION ", x_PORT_i );
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Solve with DMNH: using analytical Hessian
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  // Initial Guess
  //===
  for( int i = 1; i <= N_F_DIM; i++ ) 
    x_PORT_i->X( i ) = 10.0;
  
  // Minimisation using std::vector I/O
  //===
  if ( x_PORT_i->MNH() ) 
    print_X( "DMNH SOLUTION ", x_PORT_i );
   
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  //                      Part II: bounded fit
  //
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
  // Steer using higher level toolbox options 
  //===
  x_PORT_i->max_func_evals() = 1000; // Max func. evaluations
  x_PORT_i->max_iterations() =  100; // Max iterations
  x_PORT_i->printing_off();          // Mute stdio printout
  
  // Set some boundings
  //===
  for( int i = 1; i <= N_F_DIM; i++ ) 
    x_PORT_i->B( 1, i ) = 2.5; // xi min value
  for( int i = 1; i <= N_F_DIM; i++ ) 
    x_PORT_i->B( 2, i ) = 5.5; // xi max value
   
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  // Solve with DMNFB: bounded minimisation
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  
  // Initial Guess
  //===
  for( int i = 1; i <= N_F_DIM; i++ ) 
    x_PORT_i->X( i ) = 10.0;
  
  // Minimisation
  //===
  if ( x_PORT_i->MNFB() ) 
    print_X( "DMNFB SOLUTION", x_PORT_i );
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Solve with DMNGB: bounded and using analytical Gradient
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  // Initial Guess
  //===
  for( int i = 1; i <= N_F_DIM; i++ ) 
    x_PORT_i->X( i ) = 10.0;
  
  // Minimise
  //===
  if ( x_PORT_i->MNGB() ) 
    print_X( "DMNGB SOLUTION", x_PORT_i );
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Solve with DMNHB: bounded and using analytical Hessian
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  // Initial Guess
  //===
  for( int i = 1; i <= N_F_DIM; i++ ) 
    x_PORT_i->X( i ) = 10.0;
  
  // Minimisation using std::vector I/O
  //===
  if ( x_PORT_i->MNHB() ) 
    print_X( "DMNHB SOLUTION", x_PORT_i );
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Release the toolbox
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  delete x_PORT_i;
  
  return 0;
}
