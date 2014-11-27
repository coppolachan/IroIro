/*!
 * @file test_smear.cpp
 * @brief Tests for the propagators 
 */
#include "test_smear.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include "Measurements/GaugeM/polyakovLoop.hpp"
#include "Smearing/stoutSmear.hpp"
#include "Smearing/smearingFactories.hpp"
#include <stdio.h>

using namespace std;
using namespace Format;

#define NO_APPEND_MODE false
#define APPEND_MODE    true

void determinant_sub2x2 (double* det,SU3mat& matrix, int i , int j){
  det[0] = matrix.r(1,i)*matrix.r(2,j)-matrix.i(1,i)*matrix.i(2,j);
  det[1] = matrix.i(1,i)*matrix.r(2,j)+matrix.r(1,i)*matrix.i(2,j);

  det[0] -= matrix.r(1,j)*matrix.r(2,i)-matrix.i(1,j)*matrix.i(2,i);
  det[1] -= matrix.i(1,j)*matrix.r(2,i)+matrix.r(1,j)*matrix.i(2,i);
}


void determinant3x3 (double* det, SU3mat& matrix){
  double subdet[2];
  
  determinant_sub2x2(subdet, matrix, 1,2);

  det[0] = matrix.r(0,0)*subdet[0] - matrix.i(0,0)*subdet[1];
  det[1] = matrix.i(0,0)*subdet[0] + matrix.r(0,0)*subdet[1];
  determinant_sub2x2(subdet, matrix, 0,2);
  det[0] -= matrix.r(0,1)*subdet[0] - matrix.i(0,1)*subdet[1];
  det[1] -= matrix.i(0,1)*subdet[0] + matrix.r(0,1)*subdet[1];
  determinant_sub2x2(subdet, matrix, 0,1);
  det[0] += matrix.r(0,2)*subdet[0] - matrix.i(0,2)*subdet[1];
  det[1] += matrix.i(0,2)*subdet[0] + matrix.r(0,2)*subdet[1];
}

enum {Lambda3 = 3, Lambda8 = 8};

int plaquette(){
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


void Add_Flux(GaugeField& conf_){
  // Add a U(1) Flux 
  double BField = 8.0; // Magnetic field strength default
  std::string OutputFile = "output_flux.bin";  //default 
  int Lambda= Lambda3; //default
  //XML::node Flux_node = Gauge_node;

  //XML::read(Flux_node,"MagneticField",BField, MANDATORY);
  //XML::read(Flux_node,"Lambda",Lambda, MANDATORY);
  //XML::read(Flux_node,"Output",OutputFile);

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


}





extern  int write_eig(int, int, int, int, double*, double*);

int Test_Smear::run(){
  using namespace FieldUtils;
  Smear* SmearingObj;
  std::string OutputFile = "output_PolyakovGaugeFixed.bin";  //default

  // Smearing 
  XML::node SmearObjNode = Smear_node_; 
  /* because descend updates the node object and we want to generate 
     propagator too (it lives on the same level). 
     So we just copy Smear_node_ into a new object */

  XML::descend(SmearObjNode, "Smearing");
  SmearingFactory* Sm_Factory = 
    Smearings::createSmearingFactory(SmearObjNode);
  // Create smearing objects
  SmearingObj = Sm_Factory->getSmearing();
  int Nsmear = 10;
  XML::read(Smear_node_,"Nsmearing", Nsmear); 
  XML::read(Smear_node_,"Output", OutputFile); 

  /////////////////////////////////////////////
  // Just a check on configuration
  Staples Staple;
  PolyakovLoop PLoop(TDIR);//time dir
  CCIO::cout<< " - Plaquette : "<< Staple.plaquette(conf_)<< std::endl;
  //////////////////////////////////////
  smeared_u_ = conf_; // Copy thin links to the initial smearing

  std::complex<double> pl_avg;
  double pl_avg_adj;
   
  // Smearing and quark propagator
  for(int smear_step = 0; smear_step <= Nsmear; ++smear_step){
    CCIO::cout << "Smearing step #"<<smear_step<<"\n";
    
    if(smear_step>0){
      previous_u_ = smeared_u_;
      SmearingObj->smear(smeared_u_, previous_u_);
    }
    pl_avg = PLoop.calc_SUN(smeared_u_);
    pl_avg_adj = PLoop.calc_SUNadj(smeared_u_);
    CCIO::cout<< smear_step<<" - Plaquette (s/t): "<< Staple.plaq_s(smeared_u_) 
	      << " "<< Staple.plaq_t(smeared_u_) << std::endl;
    CCIO::cout<< smear_step<<" - Polyakov  : "<< pl_avg.real() << " "
	      <<  pl_avg.imag() << std::endl;
    CCIO::cout<< smear_step<<" - ADJ_Polyakov  : "<< pl_avg_adj << std::endl;
  }


  
  GaugeField1D PL_field = PLoop.get_PLField(smeared_u_, true);
  GaugeField1D GaugeFix;

    //    for (int s = 0; s< (SiteIndex::instance()->slsize(0,3)); s++){
    for (int s = 0; s< CommonPrms::instance()->Lvol(); s++){
      int global_site = SiteIndex::instance()->get_gsite(s);
      int gx = SiteIndex::instance()->g_x(global_site);
      int gy = SiteIndex::instance()->g_y(global_site);
      int gz = SiteIndex::instance()->g_z(global_site);
      int gt = SiteIndex::instance()->g_t(global_site);
      
      SUNmat PL_mat = mat(PL_field,s);
  
 
      PL_mat.print();
      varray_double matrix = PL_mat.getva();
      varray_double gauge_fix(matrix.size());
      write_eig(gx, gy, gz, gt, &matrix[0],&gauge_fix[0] );

      SUNmat gfixmat = SUNmat(gauge_fix);
      gfixmat.print();
      
      double det[2];  
      determinant3x3(det, gfixmat);
      double phase = atan(det[1]/det[0]);
      CCIO::cout << "Det GF: " << det[0] << " "<<det[1] << "  phase: " << phase <<"\n";
      double sign = ((det[0] > 0) - (det[0] < 0));
      double fact_r = sign*cos(-phase/3.0);
      double fact_i = sign*sin(-phase/3.0);
      for (int p = 0; p < 9; p++) {
	double re = gfixmat.r(p)*fact_r - gfixmat.i(p)*fact_i;
	double im = gfixmat.i(p)*fact_r + gfixmat.r(p)*fact_i;
	gfixmat.setr(p, re);
	gfixmat.seti(p, im);

      }
      determinant3x3(det, gfixmat);
      CCIO::cout << "New Det GF: " << det[0] << " "<<det[1] << "\n";
      
      CCIO::cout << "Is unitary:  "<< gfixmat.is_unitary() << "\n";
      SUNmat gfix = gfixmat;
      SUNmat fixed = gfixmat.dag();
      fixed *= PL_mat;
      fixed *= gfix;
      fixed.print();
      
      SetMat(GaugeFix, gfix, s);

    }

    // Apply Gauge Fix
    for (int mu = 0; mu < CommonPrms::instance()->Ndim(); mu++){
      GaugeField1D GfixMu = Mapping::shiftField(GaugeFix,mu,Mapping::Forward());
      for (int s = 0; s< CommonPrms::instance()->Lvol(); s++){
	SUNmat Umu_up = mat(conf_,s,mu);
	double det[2];
	determinant3x3(det, Umu_up);
	CCIO::cout << "Det Umu: " << det[0] << " "<<det[1] << "\n";

	CCIO::cout << "U Is unitary       :  "<< Umu_up.is_unitary() << "\n";

	SUNmat Omega = mat(GaugeFix, s).dag();
	determinant3x3(det, Omega);
	CCIO::cout << "Det Omega: " << det[0] << " "<<det[1] << "\n";
	CCIO::cout << "Omega Is unitary   :  "<< Omega.is_unitary() << "\n";
	SUNmat Omega_up = mat(GfixMu, s);
	determinant3x3(det, Omega_up);
	CCIO::cout << "Det Omega_up: " << det[0] << " "<<det[1] << "\n";
	CCIO::cout << "Omega up Is unitary:  "<< Omega_up.is_unitary() << "\n";
	
	Omega *= Umu_up;
	Omega *= Omega_up; // new U = Omega^dag * U * Omega_up 
	CCIO::cout << "Res Is unitary     :  "<< Omega.is_unitary() << "\n";

	determinant3x3(det, Omega);
	CCIO::cout << "Det Final Omega: " << det[0] << " "<<det[1] << "\n";


	
	// Test diagonal product
	double detr=1.0;
	double deti=0.0;
	for (int p=0; p< 3; p++){
	  double oldr = detr;
	  double oldi = deti;
	  detr = oldr*Omega.r(p,p) - oldi*Omega.i(p,p);
	  deti = oldi*Omega.r(p,p) + oldr*Omega.i(p,p);
	}
	printf("Determinant: %f   %f \n", detr, deti);



	// Abelian projection
	// Take the phases of the diagonal part of each gauge link
	double phase[3], new_phase[3];
	double defect=0.0;
	for (int p=0; p< 3; p++){
	  double norm = sqrt(Omega.r(p,p)*Omega.r(p,p)+Omega.i(p,p)*Omega.i(p,p));
	  //normalize
	  double real = Omega.r(p,p)/norm;
	  double imag = Omega.i(p,p)/norm;

	  phase[p] = atan2(imag, real);
	 
	  printf("norm: %f \tphase[%d]: %f\n", norm, p,phase[p]);
	  defect += phase[p];
	}
	
	defect =  defect - round((defect)/(2*PI))*2*PI;
	SU3mat abelian_mat;
	double phase_sum = 0.0;
	for (int p=0; p< 3; p++){
	  new_phase[p] = phase[p] - defect/3;
	  abelian_mat.set(p,p, cos(new_phase[p]), sin(new_phase[p]));
	  phase_sum += new_phase[p];
	}
	printf("s: %6d\t mu: %2d   defect: %f   newPhaseSum: %f\n", s, mu, defect, phase_sum);

	abelian_mat.print();
	
	SU3mat offdiag = Omega;
	SU3mat abelian = abelian_mat;
	offdiag *= abelian.dag();

	offdiag.print();

	// Apply flux here


	Omega.print();

	SetMat(conf_, abelian_mat, s, mu);



      }
    }  
    
    //Check correctness of the Polyakov gauge fix
    pl_avg = PLoop.calc_SUN(conf_);
    pl_avg_adj = PLoop.calc_SUNadj(conf_);
    CCIO::cout<< " - Check Plaquette : "<< Staple.plaquette(conf_)<< std::endl;
    CCIO::cout<< " - Plaquette (s/t): "<< Staple.plaq_s(conf_) 
	      << " "<< Staple.plaq_t(conf_) << std::endl;
    CCIO::cout<< " - Polyakov  : "<< pl_avg.real() << " "
	      <<  pl_avg.imag() << std::endl;
    CCIO::cout<< " - ADJ_Polyakov  : "<< pl_avg_adj << std::endl;


 







    // Save gauge fixed configuration 
    CCIO::SaveOnDisk< Format::Format_G >(conf_.data, OutputFile.c_str(), NO_APPEND_MODE);


  


  
  
  return 0;
}

