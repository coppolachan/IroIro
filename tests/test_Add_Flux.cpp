/*!
  @file test_Add_Flux.cpp
  @brief File for testing the Staples and Mapper classes
 */
#include "test_Add_Flux.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include "Tools/sunMat.hpp"
#include "Tools/fieldUtils.hpp"
#include "Tools/sunMatUtils.hpp"
#include <stdio.h>

using namespace std;

#define NO_APPEND_MODE false
#define APPEND_MODE    true

enum {Lambda3 = 3, Lambda8 = 8};


SU3mat SetLambda(int Lambda, double B, double Area){
  
  CCIO::cout << "Phase : "<< B*2*PI/Area << "\n";

  double exp_re = cos(B*2*PI/Area);
  double exp_im = sin(B*2*PI/Area);
  
  SU3mat Abelian;
  Abelian.zero();

  switch(Lambda){
  case Lambda3:
    Abelian.set(0, 0, exp_re, exp_im);
    Abelian.set(1, 1, exp_re, -exp_im);
    Abelian.set(2, 2, 1.0, 0.0);   
    break;
  case Lambda8:
    Abelian.set(0, 0, exp_re, exp_im);
    Abelian.set(1, 1, exp_re, exp_im);
    exp_re = cos(-B*4*PI/Area);
    exp_im = sin(-B*4*PI/Area);
    Abelian.set(2, 2, exp_re, exp_im);
    break;
  default:
    Abelian.unity();
  }

  return Abelian;
}


int Test_Add_Flux::run() {
  CCIO::header("Test Add Flux");

  Mapping::init_shiftField();

  int res1 = plaquette();


  // Add a U(1) Flux 
  double BField = 1.0; // Magnetic field strength default
  std::string OutputFile = "output_flux.bin";  //default 
  int Lambda;
  XML::node Flux_node = Gauge_node;

  XML::read(Flux_node,"MagneticField",BField, MANDATORY);
  XML::read(Flux_node,"Lambda",Lambda, MANDATORY);
  XML::read(Flux_node,"Output",OutputFile);

  // U(1) matrix
  SU3mat U1;

  // Field on z axis
  int Lx = CommonPrms::instance()->Lx();
  int Ly = CommonPrms::instance()->Ly();

  int Nx = CommonPrms::instance()->Nx();
  int Ny = CommonPrms::instance()->Ny();
  int Nz = CommonPrms::instance()->Nz();
  int Nt = CommonPrms::instance()->Nt();

  double B;

  for (int x = 0; x < Nx; x++){
    for (int y = 0; y < Ny; y++){
      for (int mu = 0; mu < 2; mu++){

	// x field
	if (mu == 0){
	  if (SiteIndex::instance()->global_x(x) == (Lx-1))
	    B = - BField*(SiteIndex::instance()->global_y(y))*Lx; 
	  else
	    B = 0.0;
	}
	// y field
	if (mu == 1)
	  B = BField*(SiteIndex::instance()->global_x(x));

	

	CCIO::cout << " x = "<< SiteIndex::instance()->global_x(x) 
		   << " y = "<< SiteIndex::instance()->global_y(y)
		   << " mu = "<< mu << "\n";

	U1 = SetLambda(Lambda, B, (double)(Lx*Ly));
	
	for (int t = 0; t < Nt; t++){
	  for (int z = 0; z < Nz; z++){
	    
	    //Get Field matrix 
	    int site = SiteIndex::instance()->site(x,y,z,t);
	    SU3mat Gmat = FieldUtils::mat(conf_, site, mu);
	    Gmat *= U1;
	    if (!Gmat.is_unitary())
	      CCIO::cout << "error matrix not unitary\n";
	    FieldUtils::SetMat(conf_, Gmat, site, mu);
	    
	  }
	}
      }
    }
  }

  // Plaquette after
  CCIO::cout << "After applying the abelian flux\n";
  plaquette();

  
  // Save configuration
  CCIO::SaveOnDisk< Format::Format_G >(conf_.data, OutputFile.c_str(), NO_APPEND_MODE);  
  

  return (res1);
}

int Test_Add_Flux::plaquette(){
  using namespace FieldUtils;
  using namespace SUNmatUtils;
  Staples wl;
  CCIO::cout<<" Plaquette = "<<  wl.plaquette(conf_) << endl;
  CCIO::cout<<" Plaquette (xy) = "<<  wl.plaq_mu_nu(conf_, 0,1) << endl;
  CCIO::cout<<" Plaquette (xz) = "<<  wl.plaq_mu_nu(conf_, 0,2) << endl;
  CCIO::cout<<" Plaquette (yz) = "<<  wl.plaq_mu_nu(conf_, 1,2) << endl;
  CCIO::cout<<" Plaquette (yt) = "<<  wl.plaq_mu_nu(conf_, 1,3) << endl;
  CCIO::cout<<" Plaquette (xt) = "<<  wl.plaq_mu_nu(conf_, 0,3) << endl;
  CCIO::cout<<" Plaquette (zt) = "<<  wl.plaq_mu_nu(conf_, 2,3) << endl;


  // DEBUG
  /*
  int Nx = CommonPrms::instance()->Nx();
  int Ny = CommonPrms::instance()->Ny();
  int Nz = CommonPrms::instance()->Nz();
  int Nt = CommonPrms::instance()->Nt();
  
  GaugeField1D stpl = wl.lower(conf_,0,1);
  for (int x = 0; x < Nx; x++){
    for (int y = 0; y < Ny; y++){
      for (int z = 0; z < Nz; z++){
	for (int t = 0; t < Nt; t++){
	  int site = SiteIndex::instance()->site(x,y,z,t);
	  double plaq = ReTr(mat(conf_,site,0)*mat_dag(stpl,site))/3.0;  // P_ij   
	  std::cout << Communicator::instance()->id()<< " " << SiteIndex::instance()->global_x(x) 
		    << " " << SiteIndex::instance()->global_y(y)
		    << " " << SiteIndex::instance()->global_z(z)
		    << " " << SiteIndex::instance()->global_t(t) 
		    << " - Plaquette : "<< plaq << "\n";
	}
      }
    }
  }
  */

  return 0;
}
