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

const GaugeField GaugeFixing_Landau::do_fix(const GaugeField& Uin)const{
  assert(Uin.Nvol()==2*Nvh_);

  GaugeField Ue = get_even(Uin);
  GaugeField Uo = get_odd(Uin);
  int Nconv = -1;

  CCIO::cout<<setw(8) <<" iter"<<setw(25)<<setiosflags(ios_base::left )
	    <<"     residual"<<"        F_value"<<endl;

  for(int it=0; it<Niter_; ++it){
    if(it%Nmeas_==0){
      double gc = gauge_cond(Ue,Uo);
      double Fval = calc_F(Ue,Uo);

      CCIO::cout<< setiosflags(ios_base::scientific);
      CCIO::cout<<setw(8) <<setiosflags(ios_base::right)<<it;
      CCIO::cout<<setw(25)<<setiosflags(ios_base::left )<<gc;
      CCIO::cout<<"  ";
      CCIO::cout<<setw(25)<<setiosflags(ios_base::right)<<Fval<<endl;
      CCIO::cout<< resetiosflags(ios_base::scientific);

      if(gc < prec_){
	Nconv = it;
	break;
      }
    }
    double wp = (it%Nreset_<Nor_) ? 1.0: wp_;
    gstep_->gfix_step(Ue,Uo,wp);

    if((it%Nreset_)== 0 && it>0) random_gtr(Ue,Uo);
  }
  if(Nconv<0){
    CCIO::cout<<"did not converge."<<endl;
    abort();
  }else{
    CCIO::cout<<"converged at iter = "<<Nconv <<endl;
  }
  return combine_eo(Ue,Uo);
}

double GaugeFixing_Landau::calc_F(const GaugeField& Ue,
				  const GaugeField& Uo)const{
  double Fval = 0.0;
  for(int mu=0; mu<NDIM_; ++mu){
    for(int site=0; site<Ue.Nvol(); ++site){
      Fval += ReTr(mat(Ue,site,mu));
      Fval += ReTr(mat(Uo,site,mu));
    }
  }
  Communicator::instance()->reduce_sum(Fval);
  return Fval/NDIM_/CommonPrms::instance()->Lvol();
}

double GaugeFixing_Landau::gauge_cond(const GaugeField& Ue,
				      const GaugeField& Uo)const{
  GaugeField1D Dele(Nvh_),Delo(Nvh_);

  for(int mu=0; mu<NDIM_; ++mu){
    Dele -= DirSlice(Ue,mu);
    Dele += shiftField_eo(DirSlice(Uo,mu),mu,Backward());

    Delo -= DirSlice(Uo,mu);
    Delo += shiftField_oe(DirSlice(Ue,mu),mu,Backward());
  }
  GaugeField1D ta = TracelessAntihermite(Dele);
  double gc = ta.data*ta.data;
  ta = TracelessAntihermite(Delo);
  gc += ta.data*ta.data;
  
  return 4.0*gc/NDIM_/NC_/CommonPrms::instance()->Lvol();
}

void GaugeFixing_Landau::random_gtr(GaugeField& Ue,GaugeField& Uo)const{
  CCIO::cout<<"  random gauge transformation performed."<<std::endl;

  GaugeField1D Gtr(2*Nvh_);
  GaugeConf_rand gtr(Gtr.format,rng_);
  gtr.init_conf(Gtr.data);
  gstep_->gauge_tr_even(Ue,Uo,get_even(Gtr));
  gstep_->gauge_tr_odd( Ue,Uo,get_odd( Gtr));
}

