#include "port_i.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PORT_i::PORT_i( int n ):
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
_dummy_int( 0 ),
_dummy_real( 0. ),
_verbosity( 1 ),
_LIV( 59+3*PORT_N_MAX_DIM ),
_LV( 78+PORT_N_MAX_DIM*(PORT_N_MAX_DIM+32) )
{
  reset( n );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PORT_i::~PORT_i()
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PORT_i::reset( int n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  _N = n;
 
  init();
  f  = NULL;
  g  = NULL;
  gh = NULL;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PORT_i::init( int alg )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define MAX_R 1e10 
{
  if ( _N == 0 ) return;  
  
  for ( unsigned  i = 0; i < _N; i++ ) 
  { 
    _D[ i ] = 1.; 
    _X[ i ] = 0.;
    B( 1, i+1 ) = -MAX_R;
    B( 2, i+1 ) =  MAX_R; 
  }
  
  _IV[ 0 ] = 0;
  CONCAT( PORT_PREFIX, ivset_( alg, _IV, _LIV, _LV, _V ) );
}
#undef MAX_R

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SCALAR& PORT_i::X( int n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( n <= 0 || n > _N ) 
  {
    warning( "access to X is out of range" );
    _dummy_real = 0.;
    return _dummy_real;
  }
  return _X[ n - 1 ];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SCALAR& PORT_i::D( int n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( n <= 0 || n > _N ) 
  {
    warning( "access to D is out of range" );
    _dummy_real = 0.;
    return _dummy_real;
  }
  return _D[ n - 1 ];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SCALAR& PORT_i::B( int i, int n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( ( n <= 0 || n > _N ) || ( i <= 0 || i > 2 ) ) 
  {
    warning( "access to B is out of range" );
    _dummy_real = 0.;
    return _dummy_real;
  }
  return _B[ 2*( n - 1 ) + i - 1 ];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SCALAR& PORT_i::V( int n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( n <= 0 || n > _LV ) 
  {
    warning( "access to V is out of range" );
    _dummy_real = 0.;
    return _dummy_real;
  }
  return _V[ n - 1 ];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int& PORT_i::IV( int n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( n <= 0 || n > _LIV ) 
  {
    warning( "access to IV is out of range" );
    _dummy_int = 0;
    return _dummy_int;
  }
  return _IV[ n - 1 ];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PORT_i::printing_on()
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  _IV[ 20 ] = 6;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PORT_i::printing_off()
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  _IV[ 20 ] = 0;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool PORT_i::MNF( int* UIPARM, SCALAR* URPARM, void (*UFPARM)(), int Neff, bool fresh_start )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{ 
  if ( f == NULL ) 
  {
    warning( "no function to minimise. Aborting" );
    return false;
  }
  if ( Neff == -1 ) 
    Neff = _N;
    
  if ( Neff > _N ) 
  {
    warning( "Space dimension to minimise is greater than the parameter space. Aborting" );
    return false;
  }
  
  if ( fresh_start ) 
    _IV[ 0 ] = 12;
    
  CONCAT( PORT_PREFIX, mnf_( Neff, _D, _X, f, _IV, _LIV, _LV, _V, UIPARM, URPARM, UFPARM ) );
  
  return true; 
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool PORT_i::MNG( int* UIPARM, SCALAR* URPARM, void (*UFPARM)(), int Neff, bool fresh_start )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{ 
  if ( f == NULL ) 
  {
    warning( "no function to minimise. Aborting" );
    return false;
  }
  
  if ( g == NULL )
  {
    warning( "no analytical Gradient provided. Aborting" );
    return false;
  }
  
  if ( Neff == -1 ) 
    Neff = _N;
    
  if ( Neff > _N ) 
  {
    warning( "Space dimension to minimise is greater than the parameter space. Aborting" );
    return false;
  }

  if ( fresh_start ) 
    _IV[ 0 ] = 12;
    
  CONCAT( PORT_PREFIX, mng_( Neff, _D, _X, f, g, _IV, _LIV, _LV, _V, UIPARM, URPARM, UFPARM ) );
  
  return true;  
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool PORT_i::MNH( int* UIPARM, SCALAR* URPARM, void (*UFPARM)(), int Neff, bool fresh_start )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( f == NULL ) {
    warning( "no F function to minimise. Aborting" );
    return false;
  }
  
  if ( gh == NULL )
  {
    warning( "no analytical Gradient and Hessian provided. Aborting" );
    return false;
  }
  
  if ( Neff == -1 ) 
    Neff = _N;
    
  if ( Neff > _N ) 
  {
    warning( "Space dimension to minimise is greater than the parameter space. Aborting" );
    return false;
  }

  if ( fresh_start ) 
    _IV[ 0 ] = 12;
    
  CONCAT( PORT_PREFIX, mnh_( Neff, _D, _X, f, gh, _IV, _LIV, _LV, _V, UIPARM, URPARM, UFPARM ) );
  
  return true;  
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool PORT_i::MNFB( int* UIPARM, SCALAR* URPARM, void (*UFPARM)(), int Neff, bool fresh_start )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{ 
  if ( f == NULL )
  {
    warning( "no function to minimise. Aborting" );
    return false;
  }
  
  if ( Neff == -1 )
    Neff = _N;
    
  if ( Neff > _N ) 
  {
    warning( "Space dimension to minimise is greater than the parameter space. Aborting" );
    return false;
  }
  
  if ( fresh_start )
    _IV[ 0 ] = 12;
    
  CONCAT( PORT_PREFIX, mnfb_( Neff, _D, _X, (SCALAR*)_B, f, _IV, _LIV, _LV, _V, UIPARM, URPARM, UFPARM ) );
  
  return true;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool PORT_i::MNGB( int* UIPARM, SCALAR* URPARM, void (*UFPARM)(), int Neff, bool fresh_start )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{ 
  if ( f == NULL ) 
  {
    warning( "no F function to minimise. Aborting" );
    return false;
  }
  if ( g == NULL ) 
  {
    warning( "no analytical Gradient provided. Aborting" );
    return false;
  }
  if ( Neff == -1 ) 
    Neff = _N;
  if ( Neff > _N ) 
  {
    warning( "Space dimension to minimise is greater than the parameter space. Aborting" );
    return false;
  }

  if ( fresh_start ) 
    _IV[ 0 ] = 12;
    
  CONCAT( PORT_PREFIX, mngb_( Neff, _D, _X, (SCALAR*)_B, f, g, _IV, _LIV, _LV, _V, UIPARM, URPARM, UFPARM ) );
  
  return true;  
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool PORT_i::MNHB( int* UIPARM, SCALAR* URPARM, void (*UFPARM)(), int Neff, bool fresh_start )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( f == NULL ) 
  {
    warning( "no function to minimise. Aborting" );
    return false;
  }
  
  if ( gh == NULL ) 
  {
    warning( "no analytical Gradient and Hessian provided. Aborting" );
    return false;
  }
  
  if ( Neff == -1 ) 
    Neff = _N;
    
  if ( Neff > _N ) 
  {
    warning( "Space dimension to minimise is greater than the parameter space. Aborting" );
    return false;
  }

  if ( fresh_start ) 
    _IV[ 0 ] = 12;
    
  CONCAT( PORT_PREFIX, mnhb_( Neff, _D, _X, (SCALAR*)_B, f, gh, _IV, _LIV, _LV, _V, UIPARM, URPARM, UFPARM ) );
  
  return true;  
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SCALAR& PORT_i::value_g( int n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( n <= 0 || n > _N ) {
    warning( "wrong index value while acessing to gradient G." );
    return _dummy_real;
  }

  return( _V[ _IV[ 27 ]-1 + ( n - 1 ) ] );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PORT_i::warning( const char* strwar )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( _verbosity < 1 ) return;
  
  std::cout << "PORT_i WARNING: " << strwar << std::endl;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PORT_i::call_f( int* UIPARM, SCALAR* URPARM, void (*UFPARM)() )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  SCALAR* f0 = &_V[ 9 ];
  int   nf  = 0;
    
  f( _N, _X, nf, f0, UIPARM, URPARM, &UFPARM );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PORT_i::call_g( int* UIPARM, SCALAR* URPARM, void (*UFPARM)() )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  SCALAR* g0 = &_V[ _IV[ 27 ]-1 ];
  int   nf  = 0;
    
  g( _N, _X, nf, g0, UIPARM, URPARM, &UFPARM );
}
