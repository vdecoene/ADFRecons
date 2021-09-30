//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//  C++ interface to FORTRAN PORT functions
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifndef PORT_I_H
#define PORT_I_H 1

#include <iostream>
#include <vector>

#include "port_i.def"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extern "C" { 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Prototyping of some FORTRAN functions
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  void CONCAT( PORT_PREFIX, mnf_(  int& N, SCALAR* D, SCALAR* X, 
               void (*CALCF)(int&, SCALAR*, int&, SCALAR*, int*, SCALAR*, void*),
	       int* IV, int& LIV, int& LV, SCALAR* V, int* UIPARM, SCALAR* URPARM, 
	       void (*UFPARM)() ) );
 
  void CONCAT( PORT_PREFIX, mng_(  int& N, SCALAR* D, SCALAR* X, 
               void (*CALCF)(int&, SCALAR*, int&, SCALAR*, int*, SCALAR*, void*), 
	       void (*CALCG)(int&, SCALAR*, int&, SCALAR*, int*, SCALAR*, void*), 
	       int* IV, int& LIV, int& LV, SCALAR* V, int* UIPARM, SCALAR* URPARM, 
	       void (*UFPARM)() ) );
	      
  void CONCAT( PORT_PREFIX, mnh_(  int& N, SCALAR* D, SCALAR* X, 
               void (*CALCF)(int&, SCALAR*, int&, SCALAR*, int*, SCALAR*, void*), 
	       void (*CALCGH)(int&, SCALAR*, int&, SCALAR*, SCALAR*, int*, SCALAR*, void*), 
	       int* IV, int& LIV, int& LV, SCALAR* V, int* UIPARM, SCALAR* URPARM, 
	       void (*UFPARM)() ) );
	      
  void CONCAT( PORT_PREFIX, mnfb_( int& N, SCALAR* D, SCALAR* X, SCALAR* B,
               void (*CALCF)(int&, SCALAR*, int&, SCALAR*, int*, SCALAR*, void*),
	       int* IV, int& LIV, int& LV, SCALAR* V, int* UIPARM, SCALAR* URPARM, 
	       void (*UFPARM)() ) );
 
  void CONCAT( PORT_PREFIX, mngb_( int& N, SCALAR* D, SCALAR* X, SCALAR* B, 
               void (*CALCF)(int&, SCALAR*, int&, SCALAR*, int*, SCALAR*, void*), 
	       void (*CALCG)(int&, SCALAR*, int&, SCALAR*, int*, SCALAR*, void*), 
	       int* IV, int& LIV, int& LV, SCALAR* V, int* UIPARM, SCALAR* URPARM, 
	       void (*UFPARM)() ) );
	      
  void CONCAT( PORT_PREFIX, mnhb_( int& N, SCALAR* D, SCALAR* X, SCALAR* B, 
               void (*CALCF)(int&, SCALAR*, int&, SCALAR*, int*, SCALAR*, void*), 
	       void (*CALCGH)(int&, SCALAR*, int&, SCALAR*, SCALAR*, int*, SCALAR*, void*), 
	       int* IV, int& LIV, int& LV, SCALAR* V, int* UIPARM, SCALAR* URPARM, 
	       void (*UFPARM)() ) );
	      
  void CONCAT( PORT_PREFIX, ivset_( int& ALG, int* IV, int& LIV, int& LV, SCALAR* V ) );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class PORT_i {
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Interface class to PORT library
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
public:

  // Constructor/Destructor
  //===
  			PORT_i( int n=0 );
  			~PORT_i();

  // Pointers to function f to minimise, Gradient g and Hessian h
  //===			
  void 			(*f)( int&, SCALAR*, int&, SCALAR*, int*, SCALAR*, void*);
  void 			(*g)( int&, SCALAR*, int&, SCALAR*, int*, SCALAR*, void*);
  void 			(*gh)(int&, SCALAR*, int&, SCALAR*, SCALAR*, int*, SCALAR*, void*);

  // Public initialisation
  //===
  void			reset( int n=0 );
  void			init( int algorithm=2 );			

  //  Accessors to data members
  //===
  int			N()   { return _N;   } // set via init only
  int			LIV() { return _LIV; } //
  int			LV()  { return _LV;  } //

  SCALAR&		X(  int n );
  SCALAR&		D(  int n );
  SCALAR&		B(  int i, int n );
  int&			IV( int n );
  SCALAR&		V(  int n );

  // Minimisation algorithms: DMNF, DMNG and DMNH and their bounded versions
  //===
  bool			MNF(  int* UIPARM=NULL, SCALAR* URPARM=NULL, void (*UFPARM)()=NULL, int Neff=-1, bool b=true );    
  bool			MNG(  int* UIPARM=NULL, SCALAR* URPARM=NULL, void (*UFPARM)()=NULL, int Neff=-1, bool b=true );
  bool			MNH(  int* UIPARM=NULL, SCALAR* URPARM=NULL, void (*UFPARM)()=NULL, int Neff=-1, bool b=true );
  bool			MNFB( int* UIPARM=NULL, SCALAR* URPARM=NULL, void (*UFPARM)()=NULL, int Neff=-1, bool b=true );    
  bool			MNGB( int* UIPARM=NULL, SCALAR* URPARM=NULL, void (*UFPARM)()=NULL, int Neff=-1, bool b=true );
  bool			MNHB( int* UIPARM=NULL, SCALAR* URPARM=NULL, void (*UFPARM)()=NULL, int Neff=-1, bool b=true );
  
  //  Minimiser I/O
  //===
  void			call_f( int* UIPARM=NULL, SCALAR* URPARM=NULL, void (*UFPARM)()=NULL );
  void			call_g( int* UIPARM=NULL, SCALAR* URPARM=NULL, void (*UFPARM)()=NULL );
  
  SCALAR&		value_f()		{ return _V[ 9 ]; }
  SCALAR&		value_g( int n );
  
  int			iteration()		{ return _IV[ 30 ]; }
  int			f_eval_counter()	{ return _IV[  5 ]; }
  int			g_eval_counter()	{ return _IV[ 29 ]; }
  
  int&			status()		{ return _IV[ 0 ];  }
  			  
  int&			max_func_evals()	{ return _IV[ 16 ]; }
  int&			max_iterations()	{ return _IV[ 17 ]; }
  
  SCALAR&		tolerance_f()		{ return _V[ 31 ];  }
  SCALAR&		tolerance_x()		{ return _V[ 32 ];  }
  SCALAR&		tolerance_false_conv()	{ return _V[ 33 ];  }
  SCALAR&		threshold_false_conv()	{ return _V[ 26 ];  }
  
  int&			printing_unit()		{ return _IV[ 20 ]; }
  void			printing_on();
  void			printing_off();

  // Warnings
  //===
  void			warningsOn()  { _verbosity = 1; }
  void			warningsOff() { _verbosity = 0; }
  void 			warning( const char* );	

private:
  int			_N;
  SCALAR		_D[ PORT_N_MAX_DIM ];
  SCALAR		_X[ PORT_N_MAX_DIM ];
  SCALAR		_B[ 2*PORT_N_MAX_DIM ];
  int			_IV[ 59+3*PORT_N_MAX_DIM ];
  int			_LIV;
  int			_LV;
  SCALAR		_V[ 78+PORT_N_MAX_DIM*(PORT_N_MAX_DIM+32) ];
  
  int			_dummy_int;
  SCALAR   		_dummy_real;
  
  char			_verbosity;
  
  			PORT_i( const PORT_i& ); // copy is protected
			
};

#endif
