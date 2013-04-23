/*!@file  gaugeFixing_Coulomb.cpp 
   @brief implementation of the Coulomb gauge fixing
   Original varaion is written by H. Matsufuru
*/
#include "gaugeFixing_Coulomb.hpp"
#include "Main/Geometry/mapping.hpp"
#include "Main/gaugeConf.hpp"
#include "Communicator/communicator.hpp"
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
    double gc_sum;
    if(it%Nmeas_==0){
      //double plaq =stpl.plaquette(combine_eo(Ue,Uo)); 
      //CCIO::cout<< "Plaquette = "<< plaq << std::endl;
      
      vector<double> gc = gauge_cond(Ue,Uo);
      gc_sum = accumulate(gc.begin(),gc.end(),0.0);
      vector<double> Fval = calc_F(Ue,Uo);
      
      //for(int t=0;t<gc.size();++t)
      //CCIO::cout<<"gc["<<t<<"]="<<gc[t]<<" Fval["<<t<<"]="<<Fval[t]<<"\n";

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
    // iteration step 
    int iter = it%Nreset_;
    if(iter==0 && it>0){
      random_gtr(Ue,Uo);
    }else{
      if(gc_sum < esdm_){
	gstep_->step_sdm(Ue,Uo);
      }else{ 
	if(iter>=Nor_) gstep_->step_SU2(Ue,Uo,ORon);
	else           gstep_->step_SU2(Ue,Uo,ORoff);
      }
    }
  }
  if(Nconv<0){
    CCIO::cout<<"did not converge."<<endl;
    abort();
  }else{
    CCIO::cout<<"converged at iter = "<<Nconv <<endl;
  }
  return combine_eo(Ue,Uo);
}

const vector<double> 
GaugeFixing_Coulomb::calc_F(const GaugeField& Ue,
			    const GaugeField& Uo)const{

  vector<double> Fval(CommonPrms::instance()->Lt(),0.0);

  for(int t=0; t<CommonPrms::instance()->Nt(); ++t){
    int gt = SiteIndex::instance()->global_t(t);

    for(int mu=0; mu<NDIM_-1; ++mu){
      for(int n=0; n<shiftSite_eo.slice_size(t,TDIR); ++n)
	Fval[gt] += ReTr(mat(Ue,shiftSite_eo.xslice(t,n,TDIR),mu));
      for(int n=0; n<shiftSite_oe.slice_size(t,TDIR); ++n)
      	Fval[gt] += ReTr(mat(Uo,shiftSite_oe.xslice(t,n,TDIR),mu));
    }
  }
  double nrm = 1.0/(NDIM_-1)/CommonPrms::instance()->Lvol();

  for(int gt=0; gt<Fval.size(); ++gt){
    double fval = Communicator::instance()->reduce_sum(Fval[gt]);
    Fval[gt] = fval*nrm;
    //CCIO::cout<<"Fval["<<gt<<"]="<<Fval[gt]<<endl;
  }
  return Fval;
}

const vector<double> 
GaugeFixing_Coulomb::gauge_cond(const GaugeField& Ue,
				const GaugeField& Uo)const{

  GaugeField1D dle = TracelessAntihermite(gstep_->umu_even(Ue,Uo));
  GaugeField1D dlo = TracelessAntihermite(gstep_->umu_odd( Ue,Uo));

  vector<double> gc(CommonPrms::instance()->Lt(),0.0);


  for(int t=0; t<CommonPrms::instance()->Nt(); ++t){
    int gt = SiteIndex::instance()->global_t(t);
    valarray<double> ve= dle.data[dle.get_sub(shiftSite_eo.xslice_map(t,TDIR))];
    gc[gt] = (ve*ve).sum();
    valarray<double> vo= dlo.data[dlo.get_sub(shiftSite_oe.xslice_map(t,TDIR))];
    gc[gt] += (vo*vo).sum();
  }
  double nrm = 4.0/(NDIM_-1)/NC_/CommonPrms::instance()->Lvol();
  
  for(int gt=0; gt<gc.size(); ++gt){
    double gcval = Communicator::instance()->reduce_sum(gc[gt]);
    //CCIO::cout<<"gcval["<<gt<<"]="<<gcval<<"\n";
    gc[gt] = gcval*nrm;
  }
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

/*
void GaugeFixing_Coulomb::re_overrelax(GaugeField& Ue,GaugeField& Uo)const{
  CCIO::cout<<" overrelaxation on ill timeslices performed."<<std::endl;

  GaugeField1D Gte(Nvh_),Gto(Nvh_);
  gstep_->step_ovrlx(Gte,Gto,Ue,Uo);

  SUNmat uni = unity();
  vector<double> gc = gauge_cond(Ue,Uo);
  for(int t=0; t<CommonPrms::instance()->Nt(); ++t){
    int gt = SiteIndex::instance()->global_t(t);
    if(gc[gt]<prec_){ 
      for(int n=0; n<shiftSite_eo.slice_size(t,TDIR); ++n)
	SetMat(Gte,uni,shiftSite_eo.xslice(t,n,TDIR));

      for(int n=0; n<shiftSite_oe.slice_size(t,TDIR); ++n)
	SetMat(Gto,uni,shiftSite_oe.xslice(t,n,TDIR));
    }
  }
  gstep_->gauge_tr_even(Ue,Uo,Gte);
  gstep_->gauge_tr_odd( Ue,Uo,Gto);
}
*/
