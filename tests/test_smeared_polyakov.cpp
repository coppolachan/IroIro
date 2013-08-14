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

extern  int write_eig(int, int, int, int, double*);

int Test_Smear::run(){
  using namespace FieldUtils;
  Smear* SmearingObj;
  
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

  /////////////////////////////////////////////
  // Just a check on configuration
  Staples Staple;
  PolyakovLoop PLoop(TDIR);//time dir
  CCIO::cout<< "Plaquette : "<< Staple.plaquette(conf_)<< std::endl;
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
    CCIO::cout<< smear_step<<"-Plaquette (s/t): "<< Staple.plaq_s(smeared_u_) 
	      << " "<< Staple.plaq_t(smeared_u_) << std::endl;
    CCIO::cout<< smear_step<<"-Polyakov  : "<< pl_avg.real() << " "
	      <<  pl_avg.imag() << std::endl;
    CCIO::cout<< smear_step<<"-ADJ_Polyakov  : "<< pl_avg_adj << std::endl;
  }


  
  GaugeField1D PL_field = PLoop.get_PLField(smeared_u_);
  for (int s = 0; s< (SiteIndex::instance()->slsize(0,3)); s++){
    int global_site = SiteIndex::instance()->get_gsite(s);
    int gx = SiteIndex::instance()->g_x(global_site);
    int gy = SiteIndex::instance()->g_y(global_site);
    int gz = SiteIndex::instance()->g_z(global_site);
    int gt = SiteIndex::instance()->g_t(global_site);

      SUNmat PL_mat = mat(PL_field,s);
      //PL_mat.print();
      varray_double matrix = PL_mat.getva();
      write_eig(gx, gy, gz, gt, &matrix[0]);
    }
  
  return 0;
}

