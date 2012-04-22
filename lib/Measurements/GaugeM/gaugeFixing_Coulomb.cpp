/*!@file  gaugeFixing_Coulomb.cpp 
   @brief implementation of the Coulomb gauge fixing
   Original varaion is written by H. Matsufuru
*/
#include "gaugeFixing_Coulomb.hpp"
#include "Main/Geometry/mapping.hpp"
#include "Communicator/communicator.h"
#include "include/macros.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/fieldUtils.hpp"
//#include "staples.hpp"
#include <numeric>

using namespace std;
using namespace FieldUtils;
using namespace SUNmatUtils;
using namespace Mapping;
using namespace SiteMap;

const GaugeField GaugeFixing_Coulomb::do_fix(const GaugeField& Uin)const{
  assert(Uin.Nvol()==2*Nvh_);
  
  //Staples stpl;
  GaugeField Ue = get_even(Uin);
  GaugeField Uo = get_odd(Uin);
  int Nconv = -1;

  CCIO::cout<<setw(8) <<" iter"<<setw(25)<<setiosflags(ios_base::left )
	    <<"     residual"<<"        F_value"<<endl;

  for(int it=0; it<Niter_; ++it){
    if(it%Nmeas_==0){
      //CCIO::cout<< "Plaquette (GaugeFixing_Coulomb)= "
      //	<< stpl.plaquette(combine_eo(Ue,Uo)) << std::endl;

      vector<double> gc = gauge_cond(Ue,Uo);
      double gc_sum = accumulate(gc.begin(),gc.end(),0.0);
      vector<double> Fval = calc_F(Ue,Uo);
      double Fval_sum = accumulate(Fval.begin(),Fval.end(),0.0);
      
      CCIO::cout<< setiosflags(ios_base::scientific);
      CCIO::cout<<setw(8) <<setiosflags(ios_base::right)<<it;
      CCIO::cout<<setw(25)<<setiosflags(ios_base::left )<<gc_sum;
      CCIO::cout<<"  ";
      CCIO::cout<<setw(25)<<setiosflags(ios_base::right)<<Fval_sum<<endl;
      CCIO::cout<< resetiosflags(ios_base::scientific);

      if(gc_sum < prec_){
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

const vector<double> GaugeFixing_Coulomb::calc_F(const GaugeField& Ue,
						 const GaugeField& Uo)const{
  vector<double> Fval(CommonPrms::instance()->Lt());

  for(int t=0; t<CommonPrms::instance()->Nt(); ++t){
    int gt = SiteIndex::instance()->global_t(t);

    for(int mu=0; mu<NDIM_-1; ++mu){
      for(int n=0; n<shiftSite_eo.slice_size(t,TDIR); ++n)
	Fval[gt] += ReTr(mat(Ue,shiftSite_eo.xslice(t,n,TDIR),mu));
      for(int n=0; n<shiftSite_oe.slice_size(t,TDIR); ++n)
	Fval[gt] += ReTr(mat(Uo,shiftSite_oe.xslice(t,n,TDIR),mu));
    }
  }
  for(int gt=0; gt<Fval.size();++gt){
    Communicator::instance()->reduce_sum(Fval[gt]);
    Fval[gt] /= (NDIM_-1)*CommonPrms::instance()->Lvol();
  }
  return Fval;
}

const vector<double> GaugeFixing_Coulomb::gauge_cond(const GaugeField& Ue,
						     const GaugeField& Uo)const{

  GaugeField1D Dele(Nvh_),Delo(Nvh_);

  for(int mu=0; mu<NDIM_-1; ++mu){
    Dele -= DirSlice(Ue,mu);
    Dele += shiftField_eo(DirSlice(Uo,mu),mu,Backward());

    Delo -= DirSlice(Uo,mu);
    Delo += shiftField_oe(DirSlice(Ue,mu),mu,Backward());
  }
  GaugeField1D tae = TracelessAntihermite(Dele);
  GaugeField1D tao = TracelessAntihermite(Delo);

  vector<double> gc(CommonPrms::instance()->Lt(),0.0);

  for(int t=0; t<CommonPrms::instance()->Nt(); ++t){
    int gt = SiteIndex::instance()->global_t(t);
    Field fe(tae.data[tae.get_sub(shiftSite_eo.xslice_map(t,TDIR))]);
    gc[gt] += fe*fe;
    Field fo(tao.data[tao.get_sub(shiftSite_oe.xslice_map(t,TDIR))]);
    gc[gt] += fo*fo;
  }

  for(int gt=0; gt<gc.size(); ++gt)
    gc[gt] *= 4.0/(NDIM_-1)/NC_/CommonPrms::instance()->Lvol();

  return gc;
}
void GaugeFixing_Coulomb::random_gtr(GaugeField& Ue,GaugeField& Uo)const{
  CCIO::cout<<"  random gauge transformation performed."<<std::endl;
  
  GaugeField1D Gtr(2*Nvh_);
  GaugeConf_rand gtr(Gtr.format,rng_);
  gtr.init_conf(Gtr.data);

  SUNmat uni = unity();
  vector<double> gc = gauge_cond(Ue,Uo);
  for(int t=0; t<CommonPrms::instance()->Nt(); ++t){
    int gt = SiteIndex::instance()->global_t(t);
    if(gc[gt]<prec_){ 
      for(int n=0; n<shiftSite.slice_size(t,TDIR); ++n)
	SetMat(Gtr,uni,shiftSite.xslice(t,n,TDIR));
    }
  }
  gstep_->gauge_tr_even(Ue,Uo,get_even(Gtr));
  gstep_->gauge_tr_odd( Ue,Uo,get_odd( Gtr));
}
