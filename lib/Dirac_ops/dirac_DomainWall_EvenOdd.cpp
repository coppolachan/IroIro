/*!--------------------------------------------------------------------------
 * @file dirac_DomainWall_EvenOdd.cpp
 *
 * @brief Definition of class methods for Dirac_optimalDomainWall_EvenOdd (5d op.)
 *-------------------------------------------------------------------------*/
#include "dirac_DomainWall_EvenOdd.hpp"
#include "Communicator/comm_io.hpp"
#include<stdlib.h>
#include<stdio.h>
#include<cassert>
#include<math.h>

using namespace std;


//-----------------------------------------------------------------------------
const Field Dirac_optimalDomainWall_EvenOdd::mult_ee(const Field& f)const{
  return Deo_.mult_hop5(f);
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_oo(const Field& f)const{
  return Deo_.mult_hop5(f);
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_ee_inv(const Field& f)const{
  return Deo_.mult_hop5_inv(f);
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_ee_dinv(const Field& f)const{
  return Deo_.mult_hop5_dinv(f);
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_oo_inv(const Field& f)const{
  return Deo_.mult_hop5_inv(f);
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_oo_dinv(const Field& f)const{
  return Deo_.mult_hop5_dinv(f);
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_eo(const Field& f)const{
  return Deo_.mult_hop5_inv(Deo_.mult(f));
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_oe(const Field& f)const{
  return Doe_.mult_hop5_inv(Doe_.mult(f));
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_eo_dag(const Field& f)const{
  return Doe_.mult_dag(Deo_.mult_hop5_dinv(f));
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_oe_dag(const Field& f)const{
  return Deo_.mult_dag(Doe_.mult_hop5_dinv(f));
}

const Field Dirac_optimalDomainWall_EvenOdd::mult(const Field& f) const{
  timeval start_, end_;
  gettimeofday(&start_,NULL);
  
  Field w(f);
  w -= mult_eo(mult_oe(f));
  return w;
  

  /*
  // just slightly faster (but only BGQ)
  Field res(Deo_.fsize());
  Doe_.mult_hop(res, f);
  //BGQ_EO_mult(res,f);
  return res;
  */

  gettimeofday(&end_,NULL);
  mult_timer += (end_.tv_sec - start_.tv_sec)*1000.0;
  mult_timer += (end_.tv_usec - start_.tv_usec) / 1000.0;   // us to ms
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_dag(const Field& f) const{
  timeval start_, end_;
  gettimeofday(&start_,NULL);

  Field w(f);
  w -= mult_oe_dag(mult_eo_dag(f));

  gettimeofday(&end_,NULL);
  multdag_timer += (end_.tv_sec - start_.tv_sec)*1000.0;
  multdag_timer += (end_.tv_usec - start_.tv_usec) / 1000.0;   // us to ms

  return w;
}

void Dirac_optimalDomainWall_EvenOdd::
md_force_eo(Field& fce, const Field& eta,const Field& zeta) const{
  Deo_.md_force_p(fce,eta,mult_ee_dinv(zeta));
  Doe_.md_force_m(fce,eta,mult_ee_dinv(zeta));
}

void Dirac_optimalDomainWall_EvenOdd::
md_force_oe(Field& fce, const Field& eta,const Field& zeta) const{
  Doe_.md_force_p(fce,eta,mult_oo_dinv(zeta));
  Deo_.md_force_m(fce,eta,mult_oo_dinv(zeta));
}

const Field Dirac_optimalDomainWall_EvenOdd::
md_force(const Field& eta,const Field& zeta) const{
  Field fce(gsize());
  md_force_eo(fce,mult_oe(eta),zeta);
  md_force_oe(fce,eta,mult_eo_dag(zeta));
  fce *= 0.5;
  return fce;
}

///////////////////////////////////////////////////////////////////////
// OPTIMIZED LIBRARIES
///////////////////////////////////////////////////////////////////////



/*
#ifdef IBM_BGQ_WILSON
void Dirac_optimalDomainWall_EvenOdd::BGQ_EO_mult(Field& w5, const Field&f5 ) const {
  timeval start_, end_;
  int f4size_ = f4size();
  int fsize_ = fsize();

  Field lpf(f4size_);
  Field lmf(f4size_), w(f4size_), v(f4size_), ey(f4size_), fy(f4size_);
  Field temp(fsize_);

  double* w_ptr    = w.getaddr(0);
  double* v_ptr    = v.getaddr(0);
  double* lpf_ptr  = lpf.getaddr(0);
  double* lmf_ptr  = lmf.getaddr(0);
  double* ey_ptr   = ey.getaddr(0);
  double* fy_ptr   = fy.getaddr(0);
  double* temp_ptr = temp.getaddr(0);

  double mass_fact   = 4.0+M0_;
  double minus_kappa = -Dw_.getKappa();

  double* pU         = const_cast<Field *>(u_)->getaddr(0);
  double* f5_ptr     = const_cast<Field&>(f5).getaddr(0);

  register int Nvol = CommonPrms::instance()->Nvol()/2;

  Format::Format_F Fformat = Dw_.get_fermionFormat();

  ///////////////////////////////  doe

  // s = 0
  BGWilsonLA_Proj_P(lpf_ptr,const_cast<Field&>(f5).getaddr(Fformat.index(0,0,N5_-1)),Nvol);
  BGWilsonLA_Proj_M(lmf_ptr,const_cast<Field&>(f5).getaddr(Fformat.index(0,0,1)),Nvol);
  BGWilsonLA_MultScalar_Add(lpf_ptr,lmf_ptr,-mq_,Nvol);
  BGWilsonLA_AXPBY(v_ptr, f5_ptr, lpf_ptr, Params.bs_[0],Params.cs_[0],Nvol);
  BGWilson_MultEO(w_ptr, pU, v_ptr, minus_kappa , 2, BGWILSON_DIRAC);
  BGWilsonLA_MultScalar(temp_ptr, w_ptr, mass_fact, Nvol);

  for(int s=1; s<N5_-1; ++s) {
    f5_ptr   += f4size_;
    temp_ptr += f4size_;

    BGWilsonLA_Proj_P(lpf_ptr,const_cast<Field&>(f5).getaddr(Fformat.index(0,0,s-1)),Nvol);
    BGWilsonLA_Proj_M(lmf_ptr,const_cast<Field&>(f5).getaddr(Fformat.index(0,0,s+1)),Nvol);
    BGWilsonLA_Add(lpf_ptr,lmf_ptr,Nvol);
    BGWilsonLA_AXPBY(v_ptr, f5_ptr, lpf_ptr, Params.bs_[s],Params.cs_[s],Nvol);
    
    gettimeofday(&start_,NULL);
    BGWilson_MultEO(w_ptr, pU, v_ptr, minus_kappa , 2, BGWILSON_DIRAC);
    gettimeofday(&end_,NULL);
    wilson_mult_timer += (end_.tv_sec - start_.tv_sec)*1000.0;
    wilson_mult_timer += (end_.tv_usec - start_.tv_usec) / 1000.0;   // us to ms


    BGWilsonLA_MultScalar(temp_ptr, w_ptr, mass_fact, Nvol);
  }

  // s = N5-1
  f5_ptr   += f4size_;
  temp_ptr += f4size_;
  BGWilsonLA_Proj_P(lpf_ptr,const_cast<Field&>(f5).getaddr(Fformat.index(0,0,N5_-2)),Nvol);
  BGWilsonLA_Proj_M(lmf_ptr,const_cast<Field&>(f5).getaddr(Fformat.index(0,0,0)),Nvol);
  BGWilsonLA_MultAddScalar(lpf_ptr,lmf_ptr,-mq_,Nvol);  
  BGWilsonLA_AXPBY(v_ptr, f5_ptr, lpf_ptr, Params.bs_[N5_-1],Params.cs_[N5_-1],Nvol);
  BGWilson_MultEO(w_ptr, pU, v_ptr, minus_kappa , 2, BGWILSON_DIRAC);
  BGWilsonLA_MultScalar(temp_ptr, w_ptr, mass_fact, Nvol);

  /// 5d hopping term
  double* temp_ptr_bdry   = temp.getaddr((N5_-1)*f4size_);
  BGWilsonLA_Proj_M(ey_ptr,temp.getaddr(Fformat.index(0,0,0)),Nvol);
  BGWilsonLA_MultAddScalar(temp_ptr_bdry,     ey_ptr,-mq_* Params.es_[0],Nvol);

  for (int s=1; s<N5_-1; ++s) {
    double* temp_ptr   = temp.getaddr(s*f4size_);
    double fact_lpf = (Params.dm_[s]/Params.dp_[s-1]);
    double fact_ey =  mq_*Params.es_[s];

    BGWilsonLA_Proj_P(lpf_ptr,temp.getaddr(Fformat.index(0,0,s-1)),Nvol);
    BGWilsonLA_Proj_M(ey_ptr,temp.getaddr(Fformat.index(0,0,s)),Nvol);
 
    BGWilsonLA_MultAddScalar(temp_ptr,     lpf_ptr,fact_lpf,Nvol);
    BGWilsonLA_MultAddScalar(temp_ptr_bdry,ey_ptr,-fact_ey,Nvol);

  }
  BGWilsonLA_Proj_P(lpf_ptr,temp.getaddr(Fformat.index(0,0,N5_-2)),Nvol);
  BGWilsonLA_MultAddScalar(temp_ptr_bdry,     lpf_ptr,(Params.dm_[N5_-1]/Params.dp_[N5_-2]),Nvol);

  double fact= 1.0/(Params.dp_[N5_-1] +mq_*Params.dm_[N5_-2]*Params.es_[N5_-2]);
  BGWilsonLA_MultScalar(temp_ptr_bdry, temp_ptr_bdry, fact, Nvol);


  for(int s=N5_-2; s>=0; --s) {
    double* temp_ptr   = temp.getaddr(s*f4size_);
    BGWilsonLA_Proj_M(lmf_ptr,temp.getaddr(Fformat.index(0,0,s+1)),Nvol);
    BGWilsonLA_Proj_P(fy_ptr,temp.getaddr(Fformat.index(0,0,N5_-1)),Nvol);
    BGWilsonLA_MultAddScalar(temp_ptr,     lmf_ptr,Params.dm_[s],Nvol);
    BGWilsonLA_MultAddScalar(temp_ptr,     fy_ptr,-mq_*Params.fs_[s],Nvol);
    BGWilsonLA_MultScalar(temp_ptr, temp_ptr, 1.0/ Params.dp_[s], Nvol);
  }

  ////////////////deo
  for(int s=0; s<N5_; ++s) {
    double* temp_ptr = temp.getaddr(s*f4size_);
    double* w5_ptr   = w5.getaddr(s*f4size_);
    BGWilsonLA_Proj_P(lpf_ptr,temp.getaddr(Fformat.index(0,0,(s+N5_-1)%N5_)),Nvol);
    BGWilsonLA_Proj_M(lmf_ptr,temp.getaddr(Fformat.index(0,0,(s+1)%N5_)),Nvol);

    if(s==0){
      BGWilsonLA_MultScalar_Add(lpf_ptr,lmf_ptr,-mq_,Nvol);
    }
    else if(s==N5_-1){
      BGWilsonLA_MultAddScalar(lpf_ptr,lmf_ptr,-mq_,Nvol);
    }
    else{
      BGWilsonLA_Add(lpf_ptr,lmf_ptr,Nvol);
    }
    BGWilsonLA_AXPBY(v_ptr, temp_ptr, lpf_ptr,
		     Params.bs_[s],Params.cs_[s],Nvol);
    
    gettimeofday(&start_,NULL);
    BGWilson_MultEO(w_ptr, pU, v_ptr, minus_kappa , 1, BGWILSON_DIRAC);
    gettimeofday(&end_,NULL);
    wilson_mult_timer += (end_.tv_sec - start_.tv_sec)*1000.0;
    wilson_mult_timer += (end_.tv_usec - start_.tv_usec) / 1000.0;   // us to ms

    BGWilsonLA_MultScalar(w5_ptr, w_ptr, mass_fact, Nvol);
  }


  double* w5_ptr_bdry   = w5.getaddr((N5_-1)*f4size_);
  BGWilsonLA_Proj_M(ey_ptr,w5.getaddr(Fformat.index(0,0,0)),Nvol);
  BGWilsonLA_MultAddScalar(w5_ptr_bdry,     ey_ptr,-mq_* Params.es_[0],Nvol);

  for (int s=1; s<N5_-1; ++s) {
    double* w5_ptr   = w5.getaddr(s*f4size_);
    double fact_lpf = (Params.dm_[s]/Params.dp_[s-1]);
    double fact_ey =  mq_*Params.es_[s];

    BGWilsonLA_Proj_P(lpf_ptr,w5.getaddr(Fformat.index(0,0,s-1)),Nvol);
    BGWilsonLA_Proj_M(ey_ptr,w5.getaddr(Fformat.index(0,0,s)),Nvol);
 
    BGWilsonLA_MultAddScalar(w5_ptr,     lpf_ptr,fact_lpf,Nvol);
    BGWilsonLA_MultAddScalar(w5_ptr_bdry,ey_ptr,-fact_ey,Nvol);

  }
  BGWilsonLA_Proj_P(lpf_ptr,w5.getaddr(Fformat.index(0,0,N5_-2)),Nvol);
  BGWilsonLA_MultAddScalar(w5_ptr_bdry,     lpf_ptr,(Params.dm_[N5_-1]/Params.dp_[N5_-2]),Nvol);
  BGWilsonLA_MultScalar(w5_ptr_bdry, w5_ptr_bdry, fact, Nvol);


  for(int s=N5_-2; s>=0; --s) {
    double* w5_ptr   = w5.getaddr(s*f4size_);
    BGWilsonLA_Proj_M(lmf_ptr,w5.getaddr(Fformat.index(0,0,s+1)),Nvol);
    BGWilsonLA_Proj_P(fy_ptr,w5.getaddr(Fformat.index(0,0,N5_-1)),Nvol);
    BGWilsonLA_MultAddScalar(w5_ptr,     lmf_ptr,Params.dm_[s],Nvol);
    BGWilsonLA_MultAddScalar(w5_ptr,     fy_ptr,-mq_*Params.fs_[s],Nvol);
    BGWilsonLA_MultScalar(w5_ptr, w5_ptr, 1.0/ Params.dp_[s], Nvol);
  }

  BGWilsonLA_MultScalar_Add(w5.getaddr(0),const_cast<Field&>(f5).getaddr(0), -1.0, Nvol*N5_);

}
#endif
*/
