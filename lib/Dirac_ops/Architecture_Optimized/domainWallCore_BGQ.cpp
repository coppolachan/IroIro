/*!--------------------------------------------------------------------------
 * @file domainWallCore_BGQ.cpp
 *
 * @brief Definition of some BGQ optimized methods for domainWallCore (5d operator)
 * Time-stamp: <2013-12-10 18:29:35 noaki>

 *-------------------------------------------------------------------------*/
#include "domainWallCore_BGQ.hpp"
#include "Tools/utils_BGQ.hpp"
#include "Fields/field_expressions.hpp"

const Field DomainWallCore_BGQ::mult(const Field& f5) const{
  Field w5(f5.size());
  (this->*mult_core)(w5,f5);
  return w5;
}

const Field DomainWallCore_BGQ::mult_dag(const Field& f5) const{
  Field w5(f5.size());
  (this->*mult_dag_core)(w5,f5);
  return w5;
}

/*! @brief total MD-force */
const Field DomainWallCore_BGQ::md_force(const Field& phi,const Field& psi)const{
  Field fce(Dw_->gsize());
  md_force_p(fce,phi,psi);
  md_force_m(fce,phi,psi);
  fce *= -0.5;
  return fce;
}

const Field DomainWallCore_BGQ::mult_hop5(const Field& f5) const{
  Field w5(f5.size()), v(f4size_);
  double* v_ptr = v.getaddr(0);

  for(int s=0; s<N5_; ++s) {
    double* w5_ptr = w5.getaddr(ff_.index(0,0,s));
    double* f5_ptr = const_cast<Field&>(f5).getaddr(ff_.index(0,0,s));
    BGWilsonLA_MultScalar(w5_ptr,f5_ptr,prms_.dp_[s],Nvol_);
  }
  for(int s=0; s<N5_-1; ++s) {
    double* w5_ptr = w5.getaddr(ff_.index(0,0,s));
    BGWilsonLA_Proj_M(v_ptr,const_cast<Field&>(f5).getaddr(ff_.index(0,0,s+1)),Nvol_);
    BGWilsonLA_MultAddScalar(w5_ptr,v_ptr,-prms_.dm_[s],Nvol_);
  }
  for(int s=1; s<N5_; ++s) {
    double* w5_ptr = w5.getaddr(ff_.index(0,0,s));
    BGWilsonLA_Proj_P(v_ptr,const_cast<Field&>(f5).getaddr(ff_.index(0,0,s-1)),Nvol_);
    BGWilsonLA_MultAddScalar(w5_ptr,v_ptr,-prms_.dm_[s],Nvol_);
  }
  double* w5_ptr = w5.getaddr(ff_.index(0,0,0));
  BGWilsonLA_Proj_P(v_ptr,const_cast<Field&>(f5).getaddr(ff_.index(0,0,N5_-1)),Nvol_);
  BGWilsonLA_MultAddScalar(w5_ptr,v_ptr,prms_.mq_*prms_.dm_[0],Nvol_);

  w5_ptr = w5.getaddr(ff_.index(0,0,N5_-1));
  BGWilsonLA_Proj_M(v_ptr,const_cast<Field&>(f5).getaddr(ff_.index(0,0,0)),Nvol_);
  BGWilsonLA_MultAddScalar(w5_ptr, v_ptr,prms_.mq_*prms_.dm_[N5_-1],Nvol_);

  return w5;
}

const Field DomainWallCore_BGQ::mult_hop5_dag(const Field& f5) const{
  Field w5(f5.size()), v(f4size_);
  double* v_ptr = v.getaddr(0);

  for(int s=0; s<N5_; ++s){
    double* w5_ptr = w5.getaddr(ff_.index(0,0,s));
    double* f5_ptr = const_cast<Field&>(f5).getaddr(ff_.index(0,0,s));
    BGWilsonLA_MultScalar(w5_ptr,f5_ptr,prms_.dp_[s],Nvol_);
  }
  for(int s=0; s<N5_-1; ++s){
    double* w5_ptr = w5.getaddr(ff_.index(0,0,s));
    BGWilsonLA_Proj_P(v_ptr,const_cast<Field&>(f5).getaddr(ff_.index(0,0,s+1)),Nvol_);
    BGWilsonLA_MultAddScalar(w5_ptr,v_ptr,-prms_.dm_[s+1],Nvol_);
  }
  for(int s=1; s<N5_; ++s){
    double* w5_ptr = w5.getaddr(ff_.index(0,0,s));
    BGWilsonLA_Proj_M(v_ptr,const_cast<Field&>(f5).getaddr(ff_.index(0,0,s-1)),Nvol_);
    BGWilsonLA_MultAddScalar(w5_ptr,v_ptr,-prms_.dm_[s-1],Nvol_);
  }
  double* w5_ptr = w5.getaddr(ff_.index(0,0,0));
  BGWilsonLA_Proj_M(v_ptr,const_cast<Field&>(f5).getaddr(ff_.index(0,0,N5_-1)),Nvol_);
  BGWilsonLA_MultAddScalar(w5_ptr,v_ptr,prms_.mq_*prms_.dm_[N5_-1],Nvol_);

  w5_ptr = w5.getaddr(ff_.index(0,0,N5_-1));
  BGWilsonLA_Proj_P(v_ptr,const_cast<Field&>(f5).getaddr(ff_.index(0,0,0)),Nvol_);
  BGWilsonLA_MultAddScalar(w5_ptr,v_ptr,prms_.mq_*prms_.dm_[0],Nvol_);
  
  return w5;
}

const Field DomainWallCore_BGQ::mult_hop5_inv(const Field& f5) const{

  Field w5(f5),lpf(f4size_),ey(f4size_),lmf(f4size_),fy(f4size_);

  double* lpf_ptr = lpf.getaddr(0);
  double* ey_ptr = ey.getaddr(0);
  double* fy_ptr = fy.getaddr(0);
  double* lmf_ptr = lmf.getaddr(0);

  double* w5_ptr_bdry   = w5.getaddr((N5_-1)*f4size_);
  BGWilsonLA_Proj_M(ey_ptr,const_cast<Field&>(w5).getaddr(ff_.index(0,0,0)),Nvol_);
  BGWilsonLA_MultAddScalar(w5_ptr_bdry,ey_ptr,-prms_.mq_* prms_.es_[0],Nvol_);

  for(int s=1; s<N5_-1; ++s){
    double* w5_ptr  = w5.getaddr(s*f4size_);
    double fact_lpf = (prms_.dm_[s]/prms_.dp_[s-1]);
    double fact_ey =  prms_.mq_*prms_.es_[s];

    BGWilsonLA_Proj_P(lpf_ptr,const_cast<Field&>(w5).getaddr(ff_.index(0,0,s-1)),Nvol_);
    BGWilsonLA_Proj_M( ey_ptr,const_cast<Field&>(w5).getaddr(ff_.index(0,0,s)),Nvol_);
 
    BGWilsonLA_MultAddScalar(w5_ptr,    lpf_ptr,fact_lpf,Nvol_);
    BGWilsonLA_MultAddScalar(w5_ptr_bdry,ey_ptr,-fact_ey,Nvol_);
  }
  BGWilsonLA_Proj_P(lpf_ptr,const_cast<Field&>(w5).getaddr(ff_.index(0,0,N5_-2)),Nvol_);
  BGWilsonLA_MultAddScalar(w5_ptr_bdry,lpf_ptr,
			   (prms_.dm_[N5_-1]/prms_.dp_[N5_-2]),Nvol_);

  double fact= 1.0/(prms_.dp_[N5_-1] +prms_.mq_*prms_.dm_[N5_-2]*prms_.es_[N5_-2]);
  BGWilsonLA_MultScalar(w5_ptr_bdry, w5_ptr_bdry, fact, Nvol_);

  for(int s=N5_-2; s>=0; --s) {
    double* w5_ptr   = w5.getaddr(s*f4size_);
    BGWilsonLA_Proj_M(lmf_ptr,const_cast<Field&>(w5).getaddr(ff_.index(0,0,s+1)),Nvol_);
    BGWilsonLA_Proj_P(fy_ptr,const_cast<Field&>(w5).getaddr(ff_.index(0,0,N5_-1)),Nvol_);
    BGWilsonLA_MultAddScalar(w5_ptr,     lmf_ptr,prms_.dm_[s],Nvol_);
    BGWilsonLA_MultAddScalar(w5_ptr,     fy_ptr,-prms_.mq_*prms_.fs_[s],Nvol_);
    BGWilsonLA_MultScalar(w5_ptr, w5_ptr, 1.0/ prms_.dp_[s], Nvol_);
  }
  return w5;
}

const Field DomainWallCore_BGQ::mult_hop5_dinv(const Field& f5) const{
  
  Field w5(f5),lpf(f4size_),ey(f4size_),lmf(f4size_),v(f4size_);
  
  Spinor* w5_ptr  = (Spinor*)w5.getaddr(0);
  Spinor* v_ptr   = (Spinor*)v.getaddr(0);
  Spinor* lpf_ptr = (Spinor*)lpf.getaddr(0);
  Spinor* lmf_ptr = (Spinor*)lmf.getaddr(0);
  Spinor* ey_ptr  = (Spinor*)ey.getaddr(0);
  
#pragma omp parallel 
  {
    int nid = omp_get_num_threads();
    int is = omp_get_thread_num()*Nvol_/nid;
    int ns = Nvol_/nid;
    
    BGWilsonLA_MultScalar(w5_ptr+is,w5_ptr+is,1.0/prms_.dp_[0],ns);
    
    for(int s=1; s<N5_-1; ++s){
      BGWilsonLA_Proj_M(lmf_ptr+is,(Spinor*)(w5.getaddr((s-1)*f4size_))+is,ns);
      BGWilsonLA_MultAddScalar(w5_ptr+is+s*Nvol_,lmf_ptr+is,prms_.dm_[s-1],ns);
      BGWilsonLA_MultScalar(w5_ptr+is+s*Nvol_,w5_ptr+is+s*Nvol_,1.0/prms_.dp_[s], ns);   
    }
    
    BGWilsonLA_Proj_M(v_ptr+is,(Spinor*)(w5.getaddr((N5_-2)*f4size_))+is,ns);   
    BGWilsonLA_MultAddScalar(w5_ptr+is+(N5_-1)*Nvol_,v_ptr+is,prms_.dm_[N5_-2],ns);
    for(int s=0; s<N5_-1; ++s) {
      BGWilsonLA_Proj_P(ey_ptr+is,(Spinor*)(w5.getaddr(s*f4size_))+is,ns); 
      BGWilsonLA_MultAddScalar(w5_ptr+is+(N5_-1)*Nvol_,
			       ey_ptr+is,-prms_.mq_*prms_.fs_[s],ns);
    }
    BGWilsonLA_MultScalar(w5_ptr+is+(N5_-1)*Nvol_,
			  w5_ptr+is+(N5_-1)*Nvol_,
			  1.0/(prms_.dp_[N5_-1] +prms_.mq_*prms_.dm_[N5_-2]*prms_.es_[N5_-2]),ns);
    
    for(int s=N5_-2; s>=0; --s){ 
      double fact_lpf = (prms_.dm_[s+1]/prms_.dp_[s]);
      
      BGWilsonLA_Proj_P(lpf_ptr+is,(Spinor*)(w5.getaddr((s+1)*f4size_  ))+is,ns); 
      BGWilsonLA_Proj_M(ey_ptr +is,(Spinor*)(w5.getaddr((N5_-1)*f4size_))+is,ns);
      BGWilsonLA_AXPBYPZ(w5_ptr  +is+s*Nvol_,
			 lpf_ptr +is,
			 ey_ptr  +is,
			 w5_ptr  +is+s*Nvol_,
			 fact_lpf,-prms_.mq_*prms_.es_[s],ns);
    }
  }
  return w5;
}

/*! @brief definitions of D_dwf */
void DomainWallCore_BGQ::mult_full(Field& w5, const Field& f5) const{ 

#pragma disjoint
  Field lpf(f4size_), lmf(f4size_), v(f4size_),w(f4size_);
  double* v_ptr   = v.getaddr(0);
  double* lpf_ptr = lpf.getaddr(0);
  double* lmf_ptr = lmf.getaddr(0);
  double* w_ptr   = w.getaddr(0);

  for(int s=0; s<N5_; ++s){
    double* f5_ptr = const_cast<Field&>(f5).getaddr(s*f4size_);
    double* w5_ptr   = w5.getaddr(s*f4size_);

    BGWilsonLA_Proj_P(lpf_ptr,
		      const_cast<Field&>(f5).getaddr(ff_.index(0,0,(s+N5_-1)%N5_)),
		      Nvol_);
    BGWilsonLA_Proj_M(lmf_ptr,
		      const_cast<Field&>(f5).getaddr(ff_.index(0,0,(s+1)%N5_)),
		      Nvol_);
    
    if(     s==0)     BGWilsonLA_MultScalar_Add(lpf_ptr,lmf_ptr,-prms_.mq_,Nvol_);
    else if(s==N5_-1) BGWilsonLA_MultAddScalar(lpf_ptr,lmf_ptr,-prms_.mq_,Nvol_);
    else              BGWilsonLA_Add(lpf_ptr,lmf_ptr,Nvol_);

    BGWilsonLA_AXPBY(v_ptr,f5_ptr,lpf_ptr,prms_.bs_[s],prms_.cs_[s],Nvol_);
    
    Dw_->mult_ptr(w_ptr, v_ptr);  
    BGWilsonLA_AXPBYPZ(w5_ptr,w_ptr,lpf_ptr,f5_ptr,4.0+prms_.M0_,-1.0,Nvol_);
  }
}

void DomainWallCore_BGQ::mult_dag_full(Field& w5,const Field& f5) const{
#pragma disjoint
  assert(w5.size()==f5.size());

  Field v5(f5.size());
  Field f4(f4size_), w(f4size_);
  Field lpf(f4size_), lmf(f4size_);
  
  int spin_idx;
  double* w_ptr = w.getaddr(0);
  
  for(int s=0; s<N5_; ++s){
    spin_idx = s*f4size_;
    double* f5_ptr = const_cast<Field&>(f5).getaddr(spin_idx);
    double* w5_ptr = w5.getaddr(spin_idx);
    double* v5_ptr = v5.getaddr(spin_idx);
    
    Dw_->mult_dag_ptr(w_ptr, f5_ptr);
    BGWilsonLA_AXPY(w5_ptr,w_ptr,f5_ptr,(4.0+prms_.M0_)*prms_.bs_[s],Nvol_);
    BGWilsonLA_AXMY(v5_ptr,w_ptr,f5_ptr,(4.0+prms_.M0_)*prms_.cs_[s],Nvol_);
  }
  
  for(int s=0; s<N5_; ++s){
    spin_idx = s*f4size_;
    double* w5_ptr = w5.getaddr(spin_idx);
    double* v5_ptr = v5.getaddr(spin_idx);
    
    BGWilsonLA_Proj_P(lpf.getaddr(0),
		      const_cast<Field&>(v5).getaddr(ff_.index(0,0,(s+1)%N5_)),Nvol_);
    BGWilsonLA_Proj_M(lmf.getaddr(0),
		      const_cast<Field&>(v5).getaddr(ff_.index(0,0,(s+N5_-1)%N5_)),Nvol_);
    if(s==N5_-1)  
      BGWilsonLA_AXPBYPZ(w5_ptr,lpf.getaddr(0),lmf.getaddr(0),w5_ptr,-prms_.mq_,1.0,Nvol_);
    else if(s==0) 
      BGWilsonLA_AXPBYPZ(w5_ptr,lpf.getaddr(0),lmf.getaddr(0),w5_ptr,1.0,-prms_.mq_,Nvol_);
    else 
      BGWilsonLA_AXPBYPZ(w5_ptr, lpf.getaddr(0),lmf.getaddr(0),w5_ptr,1.0,1.0,Nvol_);
  }
}

void DomainWallCore_BGQ::mult_offdiag(Field& w5,const Field& f5) const{ 
  
  Field lpf(f4size_), lmf(f4size_), w(f4size_), v(f4size_);
  double* w_ptr = w.getaddr(0);
  double* v_ptr = v.getaddr(0);
  double* lpf_ptr = lpf.getaddr(0);
  double* lmf_ptr = lmf.getaddr(0);
  double mass_fact= 4.0+prms_.M0_;
  double* f5_ptr  = const_cast<Field&>(f5).getaddr(0);
  double* w5_ptr  = w5.getaddr(0);
  
  /* here, Nvol_ is half-size */
  for(int s=0; s<N5_; ++s) {
    BGWilsonLA_Proj_P(lpf_ptr,
		      const_cast<Field&>(f5).getaddr(ff_.index(0,0,(s+N5_-1)%N5_)),
		      Nvol_);
    BGWilsonLA_Proj_M(lmf_ptr,
		      const_cast<Field&>(f5).getaddr(ff_.index(0,0,(s+1)%N5_)),
		      Nvol_);

    if(     s==0)     BGWilsonLA_MultScalar_Add(lpf_ptr,lmf_ptr,-prms_.mq_,Nvol_);
    else if(s==N5_-1) BGWilsonLA_MultAddScalar(lpf_ptr,lmf_ptr,-prms_.mq_,Nvol_);
    else              BGWilsonLA_Add(lpf_ptr,lmf_ptr,Nvol_);

    BGWilsonLA_AXPBY(v_ptr,f5_ptr,lpf_ptr,prms_.bs_[s],prms_.cs_[s],Nvol_);
    Dw_->mult_ptr_EO(w_ptr,v_ptr);

    BGWilsonLA_MultScalar(w5_ptr,w_ptr,mass_fact,Nvol_);
    
    f5_ptr += f4size_;
    w5_ptr += f4size_;
  }
}

void DomainWallCore_BGQ::mult_dag_offdiag(Field& w5,const Field& f5) const{
  Field v5(f5.size()),lpf(f4size_), lmf(f4size_),w(f4size_);
  double* w_ptr = w.getaddr(0);

  /* here, Nvol_ is half-size */
  for(int s=0; s<N5_; ++s){
    int spin_idx = s*f4size_;
    double* f5_ptr = const_cast<Field&>(f5).getaddr(spin_idx);
    double* w5_ptr = w5.getaddr(spin_idx);
    double* v5_ptr = v5.getaddr(spin_idx);
    
    double bs = (4.0+prms_.M0_)*prms_.bs_[s];
    double cs = (4.0+prms_.M0_)*prms_.cs_[s];

    Dw_->mult_dag_ptr_EO(w_ptr,f5_ptr);
 
    BGWilsonLA_MultScalar(w5_ptr,w_ptr,bs,Nvol_);
    BGWilsonLA_MultScalar(v5_ptr,w_ptr,cs,Nvol_);
  }
  for(int s=0; s<N5_; ++s){
    int spin_idx = s*f4size_;
    double* w5_ptr = w5.getaddr(spin_idx);
    double* v5_ptr = v5.getaddr(spin_idx);

    BGWilsonLA_Proj_P(lpf.getaddr(0),
		      const_cast<Field&>(v5).getaddr(ff_.index(0,0,(s+1)%N5_)),Nvol_);
    BGWilsonLA_Proj_M(lmf.getaddr(0),
		      const_cast<Field&>(v5).getaddr(ff_.index(0,0,(s+N5_-1)%N5_)),Nvol_);
    
    if(s == N5_-1) 
      BGWilsonLA_AXPBYPZ(w5_ptr,lpf.getaddr(0),lmf.getaddr(0),w5_ptr,-prms_.mq_,1.0,Nvol_);
    else if(s==0) 
      BGWilsonLA_AXPBYPZ(w5_ptr,lpf.getaddr(0),lmf.getaddr(0),w5_ptr,1.0,-prms_.mq_,Nvol_);
    else 
      BGWilsonLA_AXPBYPZ(w5_ptr,lpf.getaddr(0),lmf.getaddr(0),w5_ptr,1.0,1.0,Nvol_);
  }
}

/*! @brief contribution to the MD-force from forward difference */
void DomainWallCore_BGQ::
md_force_p(Field& fce,const Field& phi,const Field& psi)const{
  using namespace FieldExpression;
  register int Nvol = CommonPrms::instance()->Nvol()/2;

  Field lpf(f4size_), lmf(f4size_);
  Field w(phi.size());
  
  Spinor* lpf_ptr = (Spinor*)lpf.getaddr(0);
  Spinor* lmf_ptr = (Spinor*)lmf.getaddr(0);
  Spinor* w_ptr = (Spinor*)w.getaddr(0);
  Spinor* phi_ptr = (Spinor*)const_cast<Field&>(phi).getaddr(0);

  Field xie(f4size_);
  double* xie_ptr = xie.getaddr(0);
  double* fce_ptr = fce.getaddr(0);
  double* pU      = const_cast<Field*>(Dw_->getGaugeField_ptr())->getaddr(0);
  Spinor* psi_ptr = (Spinor*)const_cast<Field&>(psi).getaddr(0);
  double* zeta_ptr;
  std::vector<int> global_sites;
  if (EOtag_ == 1)
    global_sites = SiteIndex_EvenOdd::instance()->esec();
  else if (EOtag_ == 2)
    global_sites = SiteIndex_EvenOdd::instance()->osec();

 #pragma omp parallel 
  {
    int nid = omp_get_num_threads();
    int is = omp_get_thread_num()*Nvol/nid;
    int ns = Nvol/nid;

    for(int s=0; s<N5_; ++s){
      
      BGWilsonLA_Proj_P(lpf_ptr+is,phi_ptr+((s+N5_-1)%N5_)*Nvol+is,ns); 
      //  lpf *= -prms_.mq_;
      if(s == 0)      BGWilsonLA_MultScalar(lpf_ptr+is,lpf_ptr+is,-prms_.mq_,ns);
      BGWilsonLA_Proj_M(lmf_ptr+is,phi_ptr +((s+1)%N5_)*Nvol+is,ns); 

      //lmf *= -prms_.mq_;
      if(s == N5_-1)  BGWilsonLA_MultScalar(lmf_ptr+is, lmf_ptr+is,-prms_.mq_,ns);
      
      BGWilsonLA_Add(lpf_ptr+is,lmf_ptr+is,ns);
      BGWilsonLA_MultScalar(w_ptr+is, phi_ptr+s*Nvol+is, prms_.bs_[s],ns);
      BGWilsonLA_MultAddScalar(w_ptr+is, lpf_ptr+is, prms_.cs_[s],ns);

      //Dw_->md_force_p(fce,w,get4d(psi,s));
      
      zeta_ptr = (double*)(psi_ptr+s*Nvol);
      for(int mu=0; mu<NDIM_; ++mu){
	BGWilson_MultEO_Dir(xie_ptr,pU,w_ptr,1.0,EOtag_,BGWILSON_DIRAC,mu,BGWILSON_FORWARD);
	
	for(int site=is; site<is+ns; ++site){
	  unsigned int index = 2*NC_*ND_*site;
	  //unsigned int g_idx = gf_.index(0,global_sites[site],mu); 
	  unsigned int g_idx = 2*NC_*NC_*(global_sites[site]+2*Nvol*mu); 
	  for(int a=0; a<NC_; ++a){
	    for(int b=0; b<NC_; ++b){
	      unsigned int fce_idx = g_idx+2*(NC_*a+b);
	      for(int s=0; s<ND_; ++s){
		unsigned int ra = index+ 2*(a+NC_*s);
		unsigned int rb = index+ 2*(b+NC_*s);
		
		fce_ptr[fce_idx  ]+=zeta_ptr[rb]*xie_ptr[ra  ]+zeta_ptr[rb+1]*xie_ptr[ra+1];
		fce_ptr[fce_idx+1]+=zeta_ptr[rb]*xie_ptr[ra+1]-zeta_ptr[rb+1]*xie_ptr[ra  ];
	      }//spin
	    }//b
	  }//a 
	} //site
      }//mu
      BGQThread_Barrier(0,nid);
    }
  }
} 

/*! @brief contribution to the MD-force from backward difference */
void DomainWallCore_BGQ::
md_force_m(Field& fce,const Field& phi,const Field& psi)const{

  using namespace FieldExpression;
  register int Nvol = CommonPrms::instance()->Nvol()/2;

  Field lpf(f4size_), lmf(f4size_);
  Field w(phi.size());
  Spinor* lpf_ptr = (Spinor*)lpf.getaddr(0);
  Spinor* lmf_ptr = (Spinor*)lmf.getaddr(0);
  Spinor* w_ptr   = (Spinor*)w.getaddr(0);
  Spinor* phi_ptr = (Spinor*)const_cast<Field&>(phi).getaddr(0);
 
  Field xie(f4size_);
  double* xie_ptr = xie.getaddr(0);
  double* fce_ptr = fce.getaddr(0);
  double* pU      = const_cast<Field*>(Dw_->getGaugeField_ptr())->getaddr(0);
  Spinor* psi_ptr = (Spinor*)const_cast<Field&>(psi).getaddr(0);
  double* zt5_ptr = lpf.getaddr(0);
  double* et5_ptr = lmf.getaddr(0);
  double* zeta_ptr;
  std::vector<int> global_sites;
  if (EOtag_== 1)
    global_sites = SiteIndex_EvenOdd::instance()->esec();
  else if (EOtag_== 2)
    global_sites = SiteIndex_EvenOdd::instance()->osec();
  
#pragma omp parallel 
  {
    int nid = omp_get_num_threads();
    int ns = Nvol/nid;
    int is = omp_get_thread_num()*ns; 
      
    for(int s=0; s<N5_; ++s){
      BGWilsonLA_Proj_P(lpf_ptr+is,phi_ptr+((s+N5_-1)%N5_)*Nvol+is,ns); 
      // lpf *= -prms_.mq_;
      if(s == 0)     BGWilsonLA_MultScalar(lpf_ptr+is, lpf_ptr+is,-prms_.mq_,ns);
      BGWilsonLA_Proj_M(lmf_ptr+is,phi_ptr +((s+1)%N5_)*Nvol+is,ns); 
      //lmf *= -prms_.mq_;
      if(s == N5_-1) BGWilsonLA_MultScalar(lmf_ptr+is, lmf_ptr+is,-prms_.mq_,ns);
      
      BGWilsonLA_Add(lpf_ptr+is,lmf_ptr+is,ns);
      BGWilsonLA_MultScalar(w_ptr+is, phi_ptr+s*Nvol+is, prms_.bs_[s],ns);
      BGWilsonLA_MultAddScalar(w_ptr+is, lpf_ptr+is, prms_.cs_[s],ns);
	
      //Dw_->md_force_m(fce,w,get4d(psi,s));
	
      BGWilsonLA_MultGamma5((Spinor*)(zt5_ptr)+is, psi_ptr+s*Nvol+is, ns);
      BGWilsonLA_MultGamma5((Spinor*)(et5_ptr)+is, w_ptr+is, ns);
	
      for(int mu=0; mu<NDIM_; ++mu){
	BGWilson_MultEO_Dir(xie_ptr,pU,zt5_ptr,1.0,EOtag_,BGWILSON_DIRAC,mu,BGWILSON_FORWARD);
	
	for(int site=is; site<is+ns; ++site){
	  unsigned int index = 2*NC_*ND_*site;
	  unsigned int g_idx = 2*NC_*NC_*(global_sites[site]+2*Nvol*mu); 
	  for(int a=0; a<NC_; ++a){
	    for(int b=0; b<NC_; ++b){
	      unsigned int fce_idx = g_idx+2*(NC_*a+b);
	      for(int s=0; s<ND_; ++s){
		unsigned int ra = index+ 2*(a+NC_*s);
		unsigned int rb = index+ 2*(b+NC_*s);
		
		fce_ptr[fce_idx  ]-= xie_ptr[rb]*et5_ptr[ra  ] +xie_ptr[rb+1]*et5_ptr[ra+1];
		fce_ptr[fce_idx+1]-= xie_ptr[rb]*et5_ptr[ra+1] -xie_ptr[rb+1]*et5_ptr[ra  ];
	      }//spin
	    }//b
	  }//a 
	} //site
      }//mu
    }
  }
}  


