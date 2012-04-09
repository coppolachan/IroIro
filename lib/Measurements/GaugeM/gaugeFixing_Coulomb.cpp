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
#include <numeric>

using namespace std;
using namespace FieldUtils;
using namespace SUNmatUtils;
using namespace Mapping;
using namespace SiteMap;

const GaugeField GaugeFixing_Coulomb::fix(const GaugeField& Uin)const{
  assert(Uin.Nvol()==Nvh_);

  GaugeField Ue = get_even(Uin);
  GaugeField Uo = get_odd(Uin);
  int Niter= 0;

  for(int it=0; it<max_iter_; ++it){

    if(it%Nmeas_==0){
      vector<double> sg = calc_SG(Ue,Uo);
      double sg_sum = accumulate(sg.begin(),sg.end(),0);

      vector<double> Fval = calc_F(Ue,Uo);
      double Fval_sum = accumulate(Fval.begin(),Fval.end(),0);

      Niter = it;
      CCIO::cout<<" iter= "<< Niter
		<<" sg = "<< sg_sum <<" Fval = "<<Fval_sum <<endl;
      if(sg_sum < prec_) break;
    }

    double wp = (it%Nreset_<Niter_) ? wp_: 1.0;
    gstep_->gfix_step(Ue,Uo,wp);

    if((it%Nreset_)== 0 && it>0) random_gtr(Ue,Uo);
  }
  CCIO::cout<<"converged at iter = "<<Niter <<endl;;
  return combine_eo(Ue,Uo);
}

void GaugeFixing_Coulomb::random_gtr(GaugeField& Ue, GaugeField& Uo)const{
  CCIO::cout<<"  random gauge transformation performed."<<std::endl;
  
  GaugeField1D Gtr(2*Nvh_);
  GaugeConf_rand gtr(Gtr.format,rnd_);
  gtr.init_conf(Gtr.data);

  SUNmat uni = unity();
  vector<double> sg = calc_SG(Ue,Uo);
  for(int t=0; t<CommonPrms::instance()->Nt(); ++t){
    int gt = SiteIndex::instance()->global_t(t);
    if(sg[gt]<prec_){ 
      for(int n=0; n<shiftSite.slice_size(t,TDIR); ++n)
	SetMat(Gtr,uni,shiftSite.xslice(t,n,TDIR));
    }
  }
  gstep_->gauge_tr_even(Ue,Uo,get_even(Gtr));
  gstep_->gauge_tr_odd( Ue,Uo,get_odd(Gtr));
}

const vector<double> GaugeFixing_Coulomb::
calc_F(const GaugeField& Ue,const GaugeField& Uo)const{
  
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

const vector<double> GaugeFixing_Coulomb::
calc_SG(const GaugeField& Ue,const GaugeField& Uo)const{

  vector<double> sg(CommonPrms::instance()->Lt());
  GaugeField1D Dele(Nvh_),Delo(Nvh_);

  for(int mu=0; mu<NDIM_-1; ++mu){
    Dele -= DirSlice(Ue,mu);
    Dele += shiftField_eo(DirSlice(Uo,mu),mu,Forward());

    Delo -= DirSlice(Uo,mu);
    Delo += shiftField_oe(DirSlice(Ue,mu),mu,Forward());
  }
  GaugeField1D tae = TracelessAntihermite(Dele);
  GaugeField1D tao = TracelessAntihermite(Delo);

  for(int t=0; t<CommonPrms::instance()->Nt(); ++t){
    int gt = SiteIndex::instance()->global_t(t);
    Field fe(tae.data[tae.get_sub(shiftSite_eo.xslice_map(t,TDIR))]);
    sg[gt] += fe*fe;
    Field fo(tao.data[tao.get_sub(shiftSite_eo.xslice_map(t,TDIR))]);
    sg[gt] += fo*fo;
  }
  for(int gt=0; gt<sg.size(); ++gt)
    sg[gt] *= 4.0/(NDIM_-1)/NC_/CommonPrms::instance()->Lvol();
  
  return sg;
}

