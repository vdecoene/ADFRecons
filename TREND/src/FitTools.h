#ifndef FITTOOLS_H
#define FITTOOLS_H 1

#include "../port_i/src/port_i.h"
#include <math.h>
#include <stdio.h>
#include <setjmp.h>


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// DEFINE
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define MAX_ANTENNA 10000
                         // Note that the memory space for input antennas
                         // data is 'static'. We allocate more space than
                         // we require. This is convenient since the
			 // minimiser is written in FORTRAN.

#define N_ANTENNA_DATA 8
                         // Number of antenna related data: position x,y,z,
                         // signal arrival time t, signal amplitude A and
			 //+ X_xmax, Y_xmax, Z_xmax(Valentin Decoene 9/11/2018)

#define ADD_USER_SPACE 3
                         // Additional user space for fit data. e.i. to
                         // store the cascade direction in the amplitude fit.

#define MAX_PARAMETERS 10
                         // Maximum number of parameters for any model.


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ENUMERATIONS
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
enum fitType
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  PLAN_WAVE		=	0,
  SPHERE_WAVE = 1,
  ADF = 2
};


extern "C" {
  double gradsol_( double*, const int&, const double&, const int& );
};


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Toolbox class for fitting (FitTools)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class FitTools
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
public:
  			FitTools();
			~FitTools();

  // References to Plane Wave parameters
  //===
  double&		theta() { return _theta; }  // wave vector theta angle
  double&		phi()   { return _phi;   }  // wave vector phi angle
  // References to Sphere Wave parameters (in addition of thoses from Plane Wave)
  //===
  double&		x_xmax() { return _x_xmax; }  // sphere wave
  double&		y_xmax() { return _y_xmax; }  // sphere wave
  double&		z_xmax() { return _z_xmax; }  // sphere wave
  double&   r_xmax() {return _r_xmax;  }  // sphere wave
  double&		ts() { return _ts; }          // sphere wave
  // References to Sphere Wave parameters (in addition of thoses from Plane Wave)
  //===
  double& cerenkov_angle() { return _cerenkov_angle; } //ADF
  double& cerenkov_width() { return _cerenkov_width; } //ADF
  double& amp_corr() { return _amp_corr; } //ADF
  double& asym_corr() { return _asym_corr; } //ADF
  // Reference to wave speed
  //===
  double&		cr()    { return _cr;    }  // wave normalised speed

  // References to antennas
  //===
  int&			Na() { return _Na; }  // Number of antennas for this minimisation
  double&		xa( const int& n );   // nth antenna x position
  double&		ya( const int& n );   // nth antenna y position
  double&		za( const int& n );   // nth antenna z position
  double&		ta( const int& n );   // nth antenna signal arrival time
  double&   amplitude( const int& n ); // nth antenna amplitude

  // Reference to plane wave parameters as input parameter
  //===
  double&		theta_input() { return _theta_input; }  // input theta angle
  double&		phi_input()   { return _phi_input;   }  // input phi angle

  // Reference to Xmax position as input parameter
  //===
  double&   x_Xmax( const int& n ); // nth xmax position componant
  double&   y_Xmax( const int& n ); // nth xmax position componant
  double&   z_Xmax( const int& n ); // nth xmax position componant

  // References to error terms (one must call the propagate() method 1st to compute these terms)
  //===
  double&		sigma_t()         { return _sigma_t;      }  // Error on time measurements (in m).
  double		error( const int& i, const int& j ); // Error propagation coefficient on parameter i from the time measurement on antenna j.
  double		error( const int& i );               // Error coefficient on parameter i for identical and independently distributed
                                                             // errors on inputs with unit sigma. Multiply by the common sigma value of input.
							     // measurements to get the true parameters error.

  double		theta_error()       { return ( _errorII[ 0 ]*_sigma_t ); }  // Typical error on theta.
  double		phi_error()         { return ( _errorII[ 1 ]*_sigma_t ); }  // Typical error on phi.
  double		chi2_significance() { return _chi2_significance;         }  // Chi2 significance.

  int&			gchi2_algorithm() { return _gchi2_algorithm; } // algorithm used for the generalised chi2 computation.

  // I/O to fit options
  //===
  fitType&		fitModel()        { return _fitModel;      } // Model to fit
  bool&			fixedSpeed()      { return _fixedSpeed;    } // Is speed fixed or is it a free parameter in the fit?
  bool&			fixedSource()     { return _fixedSource;   } // Is source position fixed or is it a free parameter? (spherical fit only)
  bool&			fixedZs()         { return _fixedZs;       } // Is source z coordinate fixed or is it a free parameter? (spherical fit only)

  bool&			computeErrors()   { return _computeErrors; } // To compute errors or not to ...

  // Fit methods
  //===
  double		fit( int na=-1 );                     // Single fit
  double		scan( int na=-1 );                    // Multiple fit with a scan of parameter space for initial conditions
  void			errorPropagation( double chi2=-1.0 ); // Compute errors propagation and, if provided, the chi2 significance (only PLAN WAVE fit for now).

  double		randn();                                                       // Gaussian normal centred random variable.
  double		erfinv( double );                                              // Inverse error function.
  void 			eigenv( double** A, double* V, const int& );                   // Compute the eigen values V of matrix A. Overwrites eigen vectors in A.
  double		gchi2cdf( double*, double*, int*, const int&, const double& ); // Generalised chi2 cdf.

private:
  PORT_i		_xPORT; // C++ interface to the FORTRAN minimiser

  fitType		_fitModel;
  bool			_fixedSpeed;
  bool			_fixedSource;
  bool			_fixedZs;
  bool			_computeErrors;

  double		_Xa[ N_ANTENNA_DATA*MAX_ANTENNA + ADD_USER_SPACE ];
  double    _Xmax[ 3 ];
  int			  _Na;
  double		_errorM[ MAX_PARAMETERS ][ MAX_ANTENNA ];
  double		_errorII[ MAX_PARAMETERS ];
  double		_sigma_t;
  double		_chi2_significance;


  double		_dummyReal;
  bool			_dummyBool;

  void			_initFit();
  void			_copyFitParameters();
  int			  _nParameters;

  double		_scan_TDoA();

  void			_tred2( double**, const int&, double*, double* );
  void 			_tqli( double*, double*, const int&, double** );
  double		_pythag ( const double&, const double& );

  double 		_imhof( double*, int*, double*, const int&, const double&, double&,
  			  const double&, const double&, const double&, int& );
  double		_imhofbd( double*, int*, double*, const int&, const double&,
  			  const double&, int& );
  double 		_imhofint( double*, int*, double*, const int&, const double&, const double&,
  			  const double&, int& );
  double 		_imhoffn( double*, int*, double*, const int&, const double&, const double& );
  double 		_ruben( double*, int*, double*, const int&, const double&, const double&, int& );
  double 		_qf( double*, double*, int*, const int&, const double&, const double&,
  			  const int&, const double&, double*, int* );
  void 			_qf_counter();
  double 		_qf_exp1( double );
  double 		_qf_log1( const double&, const bool& );
  void 			_qf_order();
  double 		_qf_errbd( double, double* );
  double 		_qf_ctff( double, double* );
  double 		_qf_truncation( double, double );
  void 			_qf_findu( double*, double );
  void 			_qf_integrate( int, double, double, bool );
  double 		_qf_cfe( double );

  double _EffectiveRefractionIndex_fittools(double x0,double y0,double z0,double ns, double kr, double groundz,double xant,double yant,int stepsize);  //Compute effective refraction index as in Aires simulation (cf Matias Tueros for details)

  double 		_qf_sigsq, _qf_lmax, _qf_lmin, _qf_mean, _qf_c;
  double 		_qf_intl, _qf_ersm;
  int 			_qf_count, _qf_r, _qf_lim;
  bool 			_qf_ndtsrt, _qf_fail;
  int 			*_qf_n,*_qf_th;
  double 		*_qf_lb,*_qf_nc;
  jmp_buf 		_qf_env;

  bool			_randn_parity;
  double  		_randn_d1;
  double  		_randn_d2;
  double  		_randn_sq;

  int			_gchi2_algorithm;
  int			_mult[  MAX_ANTENNA ];
  double		_lambda[ MAX_ANTENNA ];
  double		_delta[ MAX_ANTENNA ];
  double		_lambda_153[ MAX_ANTENNA+1 ];
  double*		_p_lambda_153;
  double		_trace_155[ 7 ];

  double		_ts;

  double		_theta;
  double    _theta_input;
  double		_phi;
  double    _phi_input;
  double		_cr;

  double    _x_xmax;
  double    _y_xmax;
  double    _z_xmax;
  double    _r_xmax;

  double _cerenkov_angle;
  double _cerenkov_width;
  double _amp_corr;
  double _asym_corr;
};

// Objective functions to minimise and their gradients
//===
void PWF_function_F( int& N, double* X, int& NF, double* F, int* Na, double* Xa, void* UFPARM ); // Plane wave
void PWF_function_G( int& N, double* X, int& NF, double* G, int* Na, double* Xa, void* UFPARM );
void SWF_function_F( int& N, double* X, int& NF, double* F, int* Na, double* Xa, void* UFPARM ); // Sphere Wave
void SWF_function_G( int& N, double* X, int& NF, double* G, int* Na, double* Xa, void* UFPARM );
void ADF_function_F( int& N, double* X, int& NF, double* F, int* Na, double* Xa, void* UFPARM); // Cerenkov pattern
void ADF_function_G( int& N, double* X, int& NF, double* G, int* Na, double* Xa, void* UFPARM);

double _EffectiveRefractionIndex(double x0,double y0,double z0,double ns, double kr, double groundz,double xant,double yant,int stepsize); //Compute effective refraction index as in Aires simulation (cf Matias Tueros for details)
double _ZHSEffectiveRefractionIndex(double x0,double y0,double z0, double xa,double ya, double za, double ns, double kr); //Compute effective refraction index as in Aires simulation (cf Matias Tueros for details)
double GetRefractionIndexAtXmax(double x0,double y0,double z0,double ns, double kr); //Compute refraction index at Xmax location as in Aires simulation (cf Matias Tueros for details)
double master_equation(double w_, double X_, double Delta_, double alpha_, double n0_, double n1_); //Function of the Cerenkov angle equation
void cross_product(double &vecx, double &vecy, double &vecz, double u1x, double u1y, double u1z, double u2x, double u2y, double u2z); // Procedure to compute the cross product of two vectors
void compute_observer_position(double &x, double &y, double &z, double w, double kx, double ky, double kz, double ux, double uy, double uz, double x_Xmax, double y_Xmax, double z_Xmax, double GroundAltitude); //Procedure to compute the observer position given a specific direction and angle and Xmax position
double Cerenkov_dichotomie(double w_start, double eta, double kx, double ky, double kz, double XmaxDist, double x_Xmax, double y_Xmax, double z_Xmax, double delta, double GroundAltitude); // Solve via dichotomie the Cerenkov angle
#endif
