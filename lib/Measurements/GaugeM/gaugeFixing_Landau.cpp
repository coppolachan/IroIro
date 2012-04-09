/*!@file  gaugeFixing_Landau.cpp 
   @brief implementation of the Landau gauge fixing
   Original varaion is written by H. Matsufuru
*/
#include "gaugeFixing_Landau.hpp"
#include "Main/Geometry/mapping.hpp"
#include "Communicator/communicator.h"
#include "include/macros.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/fieldUtils.hpp"

using namespace std;
using namespace FieldUtils;
using namespace SUNmatUtils;
using namespace Mapping;

const GaugeField GaugeFixing_Landau::fix(const GaugeField& Uin)const{
  assert(Uin.Nvol()==Nvh_);

  GaugeField Ue = get_even(Uin);
  GaugeField Uo = get_odd(Uin);
  int Niter= 0;

  for(int it=0; it<max_iter_; ++it){

    if(it%Nmeas_==0){
      double Fval = calc_F(Ue,Uo);
      double sg = calc_SG(Ue,Uo);
      Niter = it;
      CCIO::cout<<" iter= "<<Niter<<" sg = "<<sg<<" Fval = "<<Fval<<endl;
      if(sg < prec_) break;
    }

    double wp = (it%Nreset_<Niter_) ? wp_: 1.0;
    gstep_->gfix_step(Ue,Uo,wp);

    if((it%Nreset_)== 0 && it>0) random_gtr(Ue,Uo);
  }
  CCIO::cout<<"converged at iter = "<<Niter <<endl;;
  return combine_eo(Ue,Uo);
}

void GaugeFixing_Landau::random_gtr(GaugeField& Ue, GaugeField& Uo)const{
  CCIO::cout<<"  random gauge transformation performed."<<std::endl;
  GaugeField1D Gtr(2*Nvh_);
  GaugeConf_rand gtr(Gtr.format,rnd_);
  gtr.init_conf(Gtr.data);
  gstep_->gauge_tr_even(Ue,Uo,get_even(Gtr));
  gstep_->gauge_tr_odd( Ue,Uo,get_odd(Gtr));
}

double GaugeFixing_Landau::
calc_F(const GaugeField& Ue,const GaugeField& Uo)const{
  
  double Fval = 0.0;
  for(int mu=0; mu<NDIM_; ++mu){
    for(int site=0; site<Nvh_; ++site){
      Fval += ReTr(mat(Ue,site,mu));
      Fval += ReTr(mat(Uo,site,mu));
    }
  }
  Communicator::instance()->reduce_sum(Fval);
  return Fval/NDIM_/CommonPrms::instance()->Lvol();
}

double GaugeFixing_Landau::
calc_SG(const GaugeField& Ue,const GaugeField& Uo)const{

  // on even sites
  GaugeField1D Delta(Nvh_);
  for(int mu=0; mu<NDIM_; ++mu){
    Delta -= DirSlice(Ue,mu);
    Delta += shiftField_eo(DirSlice(Uo,mu),mu,Forward());
  }
  GaugeField1D ta = TracelessAntihermite(Delta);
  double sg = ta.data*ta.data;

  // on odd sites  
  Delta = 0.0;
  for(int mu=0; mu<NDIM_; ++mu){
    Delta -= DirSlice(Uo,mu);
    Delta += shiftField_oe(DirSlice(Ue,mu),mu,Forward());
  }
  ta = TracelessAntihermite(Delta);
  sg += ta.data*ta.data;
  
  return 4.0*sg/NDIM_/NC_/CommonPrms::instance()->Lvol();
}

