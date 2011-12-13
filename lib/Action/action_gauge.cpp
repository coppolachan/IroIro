/*!
  @file action_gauge.cc

  @brief Implementation of the Action_gauge class
*/


#include "action_gauge.hpp"
#include "Tools/randNum_MP.h"
#include "include/format_A.h"
#include "Communicator/comm_io.hpp"

using namespace std;
using namespace SUNmat_utils;

inline SUNmat ActionGaugeWilson::u(const Field& g,int site,int dir){
  return SUNmat(g[gf_.cslice(0,site,dir)]);
}
inline SUNmat ActionGaugeWilson::u_dag(const Field& g,int site,int dir){
  return SUNmat(g[gf_.cslice(0,site,dir)]).dag();
}

inline SUNmat ActionGaugeWilson::u(const Field& g,int site){
  return SUNmat(g[sf_->cslice(0,site)]);
}
inline SUNmat ActionGaugeWilson::u_dag(const Field& g,int site){
  return SUNmat(g[sf_->cslice(0,site)]).dag();
}

double ActionGaugeWilson::calc_H(){
  //CCIO::cout<<"ActionGaugeWilson::calc_H"<<endl;
  double plaq = stpl_->plaquette(*u_);
  int NP = CommonPrms::instance()->NP();

  CCIO::cout<<"Plaq = "<<plaq<<endl;
 
  double Hgauge = Params.beta*Nvol_*NP*Ndim_*(Ndim_-1)/2*(1-plaq);
  if(nodeid_==0) printf ("H_gauge= %.15e\n",Hgauge);
  return Hgauge;
}

Field ActionGaugeWilson::md_force(const void*){
  using namespace SUNmat_utils;

  Field force(gf_.size());
  Field tmp(sf_->size());

  for(int m = 0; m < Ndim_; ++m){
    tmp=0.0;
    for(int n=0; n< Ndim_; ++n){
      if(n != m){
	tmp += stpl_->upper(*u_,m,n);
	tmp += stpl_->lower(*u_,m,n);
      }
    }
    for(int site=0; site<Nvol_; ++site){

      SUNmat staple_sum = u_dag(tmp,site);
      SUNmat link = u(*u_,site,m);
      SUNmat pl = link*staple_sum;
      
      force.set(gf_.cslice(0,site,m), anti_hermite(pl));
    }
  }
  force *= Params.beta/Nc_/2.0;
  return force;
}

void ActionGaugeWilson::
init(Field& P,const RandNum& rand,const void*){
  CCIO::cout<<"ActionGaugeWilson::init"<<endl;
  md_mom(P,rand);
}

void ActionGaugeWilson::md_mom(Field& P,const RandNum& rand){
  if(CommonPrms::instance()->Nc() == 3){
    md_mom_su3(P,rand);
  }else{
    throw "MD momentum P is not implemented for Nc_!=3.";
  }
}

void ActionGaugeWilson::md_mom_su3(Field& P,const RandNum& rand){

  static const double sq3i = 1.0/sqrt(3.0);
  
  Format::Format_A fmt(CommonPrms::instance()->Nvol());
  valarray<double> pj(fmt.size());

  MPrand::mp_get_gauss(pj,rand,fmt);
  pj *= sqrt(2.0); // to get Px,mu^a with the distribution exp(-1/2*P*P)
  pj /= 2.0;       // to get i*Px,mu =i*Px,mu^a*lambda^a/2 below
  for(int d=0; d<Ndim_; ++d){
    for(int site=0; site<Nvol_; ++site){

      double pjp[] = {0.0,   
		      pj[fmt.index(2,site,d)]+sq3i*pj[fmt.index(7,site,d)],
		      pj[fmt.index(1,site,d)], 
		      pj[fmt.index(0,site,d)],          
		      pj[fmt.index(4,site,d)], 
		      pj[fmt.index(3,site,d)],          
		     -pj[fmt.index(1,site,d)], 
		      pj[fmt.index(0,site,d)],           
		      0.0,  
		     -pj[fmt.index(2,site,d)]+sq3i*pj[fmt.index(7,site,d)],
		      pj[fmt.index(6,site,d)], 
		      pj[fmt.index(5,site,d)],          
		     -pj[fmt.index(4,site,d)], 
		      pj[fmt.index(3,site,d)],           
		     -pj[fmt.index(6,site,d)], 
		      pj[fmt.index(5,site,d)],           
		      0.0,  
		      -2*sq3i*pj[fmt.index(7,site,d)]};
      
      P.set(gf_.cslice(0,site,d),valarray<double>(pjp,18));
    }
  }
}
