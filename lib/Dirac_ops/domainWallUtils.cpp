#include <stdlib.h>
#include <stdio.h>
#include <cassert>
#include <math.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_elljac.h>
#include "domainWallUtils.hpp"

using namespace std;

// Constructors for DomainWall Parameters classes
//======================================================================
DomainWallParams::DomainWallParams(XML::node dw_node,DWFType Type){
  /* temporal hack*/
  XML::node mynode = dw_node;
  XML::descend(mynode,"BaseKernel", MANDATORY);
  XML::read(mynode, "mass", M0_,MANDATORY);

  XML::read(dw_node, "N5d", N5_,MANDATORY);
  XML::read(dw_node, "b",   b_, MANDATORY);
  XML::read(dw_node, "c",   c_, MANDATORY);
  XML::read(dw_node, "mass",mq_,MANDATORY);

  if(Type == PauliVillars) mq_= 1.0;
  
  XML::node ApproxNode = dw_node.child("approximation");
  if(ApproxNode !=NULL) {
    const char* Approx_name = ApproxNode.attribute("name").value();
    if(!strcmp(Approx_name, "Zolotarev")){
      double lambda_min, lambda_max;
      XML::read(ApproxNode, "lambda_min", lambda_min); 
      XML::read(ApproxNode, "lambda_max", lambda_max); 
      omega_= DWF::getOmega(N5_,lambda_min,lambda_max);
    }
    if (!strcmp(Approx_name, "Tanh"))  
      for(int s=0; s<N5_; ++s) omega_.push_back(1.0);
  }else{
    CCIO::cout << "Error: missing [approximation] node or wrong entry\n";
    abort();
  }
  set_arrays();   // setup of the member arrays
}

DomainWallParams::DomainWallParams(const DomainWallParams& prms,DWFType Type)
  :omega_(prms.omega_),N5_(omega_.size()),
   b_(prms.b_),c_(prms.c_),M0_(prms.M0_),mq_(prms.mq_){
  if(Type == PauliVillars) mq_= 1.0;
  set_arrays();
}

DomainWallParams::DomainWallParams(double b,double c,double M0,double mq,
				   const std::vector<double>& omega)
  :omega_(omega),N5_(omega.size()),
   b_(b),c_(c),M0_(M0),mq_(mq){set_arrays();}

void DomainWallParams::set_arrays(){
  
  for(int s=0; s<N5_; ++s){
    bs_.push_back(0.5*(b_*omega_[s] +c_));
    cs_.push_back(0.5*(b_*omega_[s] -c_));
    dp_.push_back(bs_[s]*(4.0 +M0_)+1.0);
    dm_.push_back(1.0 -cs_[s]*(4.0 +M0_));
  }
  es_.push_back(dm_[N5_-1]/dp_[0]);
  fs_.push_back(dm_[0]);
  for (int s=1; s<N5_-1; ++s) {
    es_.push_back(es_[s-1]*(dm_[s-1]/dp_[s]));
    fs_.push_back(fs_[s-1]*(dm_[s]/dp_[s-1]));
  }
}

/* @brief Namespace definining useful functions for DWF*/
namespace DWF {

  inline double set_vs( int is, int ns, double kprime ){
    double ekprime = gsl_sf_ellint_Kcomp( kprime , 0 ); 
    double vs = is * ekprime / ns; 
    return vs;
  }

  const vector<double> getOmega(int Ns,double lambda_min,double lambda_max){
    double u, m;
    double sn, cn, dn; 
    double kprime = sqrt(1.0-(lambda_min/lambda_max)*(lambda_min/lambda_max));
    
    vector<double> omegas(Ns);

    for(int ii=0; ii<Ns; ++ii){
      int is = 2*ii + 1;
      m = kprime*kprime;
      double vs = set_vs( is, Ns*2, kprime );
      gsl_sf_elljac_e( vs, m, &sn, &cn, &dn );
      double sn2 = sn * sn;
      double kappaprime2 = kprime * kprime;
      omegas[ii] = (1.0/lambda_min)* sqrt(1.0-kappaprime2*sn2);
    }
#if VERBOSITY>2
    for( int ii=0; ii<Ns; ++ii) printf("%24.16E\n", omegas[ii] );
#endif
    return omegas;
  }
}
