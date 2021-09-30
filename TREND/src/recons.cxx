#include "FitTools.h"
#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <math.h>

#define RAD2DEG 180./M_PI

int main(int argc, char* argv[])
{
  int recons_type = atoi(argv[1]);
  std::cout << "### Reconstuction ###" << std::endl;

  char* text_dir = argv[2];

  // Load antennas positions
  long int antenna;
  long int init_ant;
  double x;
  double y;
  double z;
  std::vector<double> coord(10000000,0);
  std::ostringstream file_antennapos;
  file_antennapos << text_dir << "/coord_antennas.txt";
  FILE* fid_antennapos  = fopen( file_antennapos.str().c_str(),   "r" );
  std::cout << "Loading antennas positions file " << file_antennapos.str().c_str() << std::endl;
  int i = 0;
  while ( fscanf( fid_antennapos, "%ld  %lf %lf %lf ", &antenna,&x,&y,&z )>0 )  {
    if (i==0)  {  // first line
      init_ant = antenna;
    }
    coord.at(i)=x;
    coord.at(i+1)=y;
    coord.at(i+2)=z;
    i=i+3;
  }
  // Load data files (I/Os)
  std::ostringstream file_in;
  file_in << text_dir << "/Rec_coinctable.txt";
  FILE* fid_in  = fopen( file_in.str().c_str(),	"r" );
  std::cout<<"Looking for input file..."<<std::endl;
  if (fid_in==0)  {
    std::cout << "Could not find file " << file_in.str().c_str() << std::endl;
    return 0;
  }
  else  {
    std::cout << "File " << file_in.str().c_str() << " successfully opened." <<  std::endl;
  }
  FILE* fid_out ;
  FitTools fBox;    // The fit toolbox.

  fid_in  = fopen( file_in.str().c_str(),	"r" );
  std::ostringstream file_out;

  // Fit settings
  if (recons_type==0) {
    file_out << text_dir << "/Rec_plane_wave_recons.txt";
    std::cout << "Now writing plane wave recons results to file " << file_out.str() << std::endl;
    fBox.fitModel()  = PLAN_WAVE;
  }
  else if (recons_type==1) {
    file_out << text_dir << "/Rec_sphere_wave_recons.txt";
    std::cout << "Now writing Spherical wave recons results to file " << file_out.str() << std::endl;
    fBox.fitModel()  = SPHERE_WAVE;
  }
  else if (recons_type==2) {
    file_out << text_dir << "/Rec_adf_recons.txt";
    std::cout << "Now writing Angular Distribution Function recons results to file " << file_out.str() << std::endl;
    fBox.fitModel()  = ADF;
  }
  else {
    std::cout<<"No analysis technic found"<<std::endl;
  }

  fBox.fixedSpeed() = true;	 // true: wave speed is a fixed parameter in the fit / false: wave speed is an additional free parameter in the fit
  fBox.cr() = 1;        // Speed value to use if fixed. By default it is set to 1.0 at initialisation.
  fBox.sigma_t()     = 3.0; // error on times, in m.
  fid_out = fopen( file_out.str().c_str(), "w" );

  // Initialize variables
  int j = 1;
  int k =1;
  int l =1;
  long int dummy_int;
  double dummy;
  long int iCoinc=0;
  long int iCoinc_prev=0;
  long int iCoinc_rec=0;
  long int Id_antenna=0;
  long int index=0;
  double trig=0;
  double amp=0;
  double x_xmax=0;
  double y_xmax=0;
  double z_xmax=0;
  double theta_rec=0;
  double phi_rec=0;

  std::ostringstream file_input_angles; //input file from previous reconsructions
  std::ostringstream file_input_xmax; //input file from previous reconsructions
  FILE* fid_input_angles;
  FILE* fid_input_xmax;

  double GroundAltitude = 1086.;//2900.;
  std::cout<<"Ground Altitude set at "<<GroundAltitude<<std::endl;

  if (recons_type==1 || recons_type==2){
    file_input_angles << text_dir << "/Rec_plane_wave_recons.txt";
    fid_input_angles = fopen( file_input_angles.str().c_str(), "r" );
  }
  if (recons_type==2){
    file_input_xmax << text_dir << "/Rec_sphere_wave_recons.txt";
    fid_input_xmax = fopen( file_input_xmax.str().c_str(), "r" );
    fscanf( fid_input_xmax, "%ld %ld %lf %lf %lf %lf %lf %lf",&iCoinc_rec, &dummy_int, &dummy, &dummy, &x_xmax ,&y_xmax, &z_xmax, &dummy); // read recons of Xmax from spherical recons for the first coinc !!!
  }

  while ( fscanf( fid_in, "%ld %ld %lf %lf",&Id_antenna,&iCoinc,&trig,&amp)>0 )  {  // New coinctable.txt format

      if (iCoinc_prev != iCoinc || i==0 )  {  // New coinc

        std::cout << "Now writting recons from coinc " << iCoinc_prev << " to file." << std::endl;

        // Reconstruction for previous coinc
        fBox.Na() = j;
        if (fBox.Na()>3){

          if (recons_type==1 || recons_type==2){
            fscanf( fid_input_angles, "%ld %ld %lf %lf %lf %lf %lf %lf \n", &iCoinc_rec, &dummy_int, &theta_rec, &dummy, &phi_rec, &dummy, &dummy, &dummy); // read recons of Xmax from spherical recons
            fBox.theta_input() = theta_rec;
            fBox.phi_input() = phi_rec;
            if (iCoinc_rec != iCoinc_prev){std::cout<<"Warning wrong input coinc event"<<std::endl;}
            std::cout<<"theta_rec = "<<theta_rec<<" phi_rec = "<<phi_rec<<std::endl;
          }
          if (recons_type==2){
            //fscanf( fid_input_xmax, "%ld %ld %lf %lf %lf %lf %lf %lf",&iCoinc_rec, &dummy_int, &dummy, &dummy, &x_xmax ,&y_xmax, &z_xmax, &dummy); // read recons of Xmax from spherical recons
            //if (iCoinc_rec != iCoinc_prev){std::cout<<"Warning wrong input coinc event"<<std::endl;}
            std::cout<<"x_xmax = "<<x_xmax<<" y_xmax = "<<y_xmax<<" z_xmax = "<<z_xmax<<std::endl; //print the Xmax position for the current coinc !!
          }

            std::cout<<"Now calling scan" <<std::endl;
            double chi2 = fBox.scan();  // Perform a 3D scan for the fBox.Na() first antennas
            std::cout<<"Done" <<std::endl;

          if (recons_type==0) {	// Plan wave
            fprintf( fid_out, "%ld %3.0d %12.5le %12.5le %12.5le %12.8le %12.5le %12.5le \n", iCoinc_prev, fBox.Na(), fBox.theta()*RAD2DEG, fBox.theta_error()*RAD2DEG, fBox.phi()*RAD2DEG, fBox.phi_error()*RAD2DEG, chi2,fBox.chi2_significance()); // plan wave parameters
          }
          else if (recons_type==1){ // Spherical wave
              fprintf( fid_out, "%ld %3.0d %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le\n", iCoinc_prev, fBox.Na(), chi2, fBox.chi2_significance(), - fBox.r_xmax()*(sin(fBox.theta())*cos(fBox.phi())), - fBox.r_xmax()*(sin(fBox.theta())*sin(fBox.phi())), GroundAltitude - fBox.r_xmax()*cos(fBox.theta()), fBox.ts()); // Physical spherical wave parameters
              std::cout<<"x = "<<- fBox.r_xmax()*(sin(fBox.theta())*cos(fBox.phi()))<<" y = "<<- fBox.r_xmax()*(sin(fBox.theta())*sin(fBox.phi()))<<" z = "<<- fBox.r_xmax()*cos(fBox.theta())<<std::endl;
            }
            else if (recons_type==2){ //ADF
                fprintf( fid_out, "%ld %3.0d %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le \n", iCoinc_prev, fBox.Na(), fBox.theta()*RAD2DEG, fBox.theta_error()*RAD2DEG, fBox.phi()*RAD2DEG, fBox.phi_error()*RAD2DEG, chi2,fBox.chi2_significance(), fBox.cerenkov_width(), fBox.amp_corr()); // ADF parameters
              }
            //------
            if (recons_type==2){
              fscanf( fid_input_xmax, "%ld %ld %lf %lf %lf %lf %lf %lf",&iCoinc_rec, &dummy_int, &dummy, &dummy, &x_xmax ,&y_xmax, &z_xmax, &dummy); // read recons of Xmax from spherical recons for the next coinc !!
              if (iCoinc_rec != iCoinc_prev){std::cout<<"Warning wrong input coinc event"<<std::endl;}
              //std::cout<<"x_xmax = "<<x_xmax<<" y_xmax = "<<y_xmax<<" z_xmax = "<<z_xmax<<std::endl;
            }
            //-----
          }

        // Re-initialize variables
        j=1;
        k=0;
        iCoinc_prev = iCoinc;
        l++;
      }
      else  {   // Same coinc
        j++;
      }
      index = Id_antenna-init_ant;  // Note:this means that ALL antennas have to be given in table starting from init_ant.
      // Now read antenna positions
      fBox.xa(j) = coord.at(index*3);
      fBox.ya(j) = coord.at(index*3+1);
      fBox.za(j) = coord.at(index*3+2);

      fBox.ta(j) = trig*2.997924580e8;
      fBox.amplitude(j) = amp;

      fBox.x_Xmax(j) = x_xmax;
      fBox.y_Xmax(j) = y_xmax;
      fBox.z_Xmax(j) = z_xmax;
      //std::cout<<"fBox.x_Xmax(j) = "<<fBox.x_Xmax(j)<<" fBox.y_Xmax(j) = "<<fBox.y_Xmax(j)<<" fBox.z_Xmax(j) = "<<fBox.z_Xmax(j)<<std::endl;

    }  // while

  // Reconstruction for last coinc
  fBox.Na() = j;
  if ( fBox.Na()>3  && (fBox.ta(2)+fBox.ta(3)>0 ) ) {  // Timing performed for intercorrelation
    double chi2 = fBox.scan();  // Perform the source reconstruction for the fBox.Na() first antennas
    if (recons_type==0)  {   // Plan wave
      fprintf( fid_out, "%ld %3.0d %12.5le %12.5le %12.5le %12.8le %12.5le %12.5le \n", iCoinc_prev, fBox.Na(), fBox.theta()*RAD2DEG, fBox.theta_error()*RAD2DEG, fBox.phi()*RAD2DEG, fBox.phi_error()*RAD2DEG, chi2,fBox.chi2_significance() ); // plan wave parameters
    }
    else if (recons_type==1){ //Spherical wave
      fprintf( fid_out, "%ld %3.0d %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le\n", iCoinc_prev, fBox.Na(), chi2, fBox.chi2_significance(), - fBox.r_xmax()*(sin(fBox.theta())*cos(fBox.phi())), - fBox.r_xmax()*(sin(fBox.theta())*sin(fBox.phi())), GroundAltitude - fBox.r_xmax()*cos(fBox.theta()), fBox.ts()); // Physical Hyperbolic wave parameters
    }
    else if (recons_type==2){ //ADF
        fprintf( fid_out, "%ld %3.0d %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le\n", iCoinc_prev, fBox.Na(), fBox.theta()*RAD2DEG, fBox.theta_error()*RAD2DEG, fBox.phi()*RAD2DEG, fBox.phi_error()*RAD2DEG, chi2,fBox.chi2_significance(), fBox.cerenkov_width(), fBox.amp_corr()); // ADF parameters
      }
  }
  // Close I/Os
  fclose( fid_out );
  fclose( fid_in  );

  return 0;
}
