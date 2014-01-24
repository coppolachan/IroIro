#include "domainWallCore.hpp"

#include "domainWallCore.hpp"
#include "Communicator/comm_io.hpp"
#include "Fields/field_expressions.hpp"

const Field DomainWallCore::mult(const Field& f5) const{
  Field w5(f5.size());
  (this->*mult_core)(w5,f5);
  return w5;
}

const Field DomainWallCore::mult_dag(const Field& f5) const{
  Field w5(f5.size());
  (this->*mult_dag_core)(w5,f5);
  return w5;
}

const Field DomainWallCore::md_force(const Field& phi,const Field& psi)const{
  Field fce(Dw_->gsize());
  md_force_p(fce,phi,psi);
  md_force_m(fce,phi,psi);
  fce *= -0.5;
  return fce;
}

/*
void DomainWallCore::Dminus(Field& w5,const Field& f5) const{
  //1-c_s * D_w(-M)
  w5 = f5;
  for(int s=0; s<N5_; ++s) {
    Field tmp =  Dw_->mult(get4d(f5,s));
    tmp *= -prms_.cs_[s]; // = [-c_s * D_w(-M)]f5
    add5d(w5,tmp,s); //= [1-c_s * D_w(-M)]f5
  }
}
*/
const Field DomainWallCore::get4d(const Field& f5,int s) const{
  Field w(f4size_);
  for(int i=0; i<f4size_; ++i) w.set(i, f5[s*f4size_+i]);
  return w;
}
void DomainWallCore::add5d(Field& f5,const Field& f4_1,const Field& f4_2,int s) const{
  for(int i=0; i<f4size_; ++i) f5.add(s*f4size_+i,f4_1[i]+f4_2[i]);
}
void DomainWallCore::add5d(Field& f5,const Field& f4,int s) const{
  f5.add(std::slice(s*f4size_,f4size_,1),f4.getva());
}
void DomainWallCore::set5d(Field& f5,const Field& f4,int s) const{
  for(int i=0; i<f4size_; ++i) f5.set(s*f4size_+i,f4[i]);
}
void DomainWallCore::mul5d(Field& f5,double fac,int s) const{
  for(int i=0; i<f4size_; ++i) f5.set(s*f4size_+i,fac*f5[s*f4size_+i]);
}

/*! @brief mult without 4D-hopping parameters */

const Field DomainWallCore::mult_hop5(const Field& f5) const{
  Field w5(f5.size());
  Field v(f4size_);

  for(int s=0; s<N5_; ++s){
    Field v = get4d(f5,s);
    v *= prms_.dp_[s];
    set5d(w5,v,s);
  }
  for(int s=0; s<N5_-1; ++s){
    projM_(v,f5,s+1);
    v *= -prms_.dm_[s];
    add5d(w5,v,s);
  }
  for(int s=1; s<N5_; ++s){
    projP_(v,f5,s-1);
    v *= -prms_.dm_[s];
    add5d(w5,v,s);
  }

  projP_(v,f5,N5_-1);
  v *= prms_.mq_*prms_.dm_[0];
  add5d(w5,v,0);

  projM_(v,f5,0);
  v *= prms_.mq_*prms_.dm_[N5_-1];
  add5d(w5,v,N5_-1);

  return w5;
}

const Field DomainWallCore::mult_hop5_dag(const Field& f5) const{
  Field w5(f5.size());
  Field v(f4size_);

  for(int s=0; s<N5_; ++s){
    v = get4d(f5,s);
    v *= prms_.dp_[s];
    set5d(w5,v,s);
  }
  for(int s=0; s<N5_-1; ++s){
    projP_(v,f5,s+1);
    v *= -prms_.dm_[s+1];
    add5d(w5,v,s);
  }
  for(int s=1; s<N5_; ++s){
    projM_(v,f5,s-1);
    v *= -prms_.dm_[s-1];
    add5d(w5,v,s);
  }
  projM_(v,f5,N5_-1);
  v *= prms_.mq_*prms_.dm_[N5_-1];
  add5d(w5,v,0);

  projP_(v,f5,0);
  v *= prms_.mq_*prms_.dm_[0];
  add5d(w5,v,N5_-1);

  return w5;
}

const Field DomainWallCore::mult_hop5_inv(const Field& f5) const{
  Field w5(f5);
  Field lf(f4size_),fy(f4size_);

  projM_(fy,w5,0);
  fy *= -prms_.mq_* prms_.es_[0];
  add5d(w5,fy,N5_-1);

  for(int s=1; s<N5_-1; ++s){
    projP_(lf,w5,s-1);
    projM_(fy,w5,s);

    lf *= prms_.dm_[s]/prms_.dp_[s-1];
    add5d(w5,lf,s);

    fy *= -prms_.mq_*prms_.es_[s];
    add5d(w5,fy,N5_-1);
  }
  projP_(lf,w5,N5_-2);
  lf *= (prms_.dm_[N5_-1]/prms_.dp_[N5_-2]);
  add5d(w5,lf,N5_-1);

  mul5d(w5,1.0/(prms_.dp_[N5_-1] +prms_.mq_*prms_.dm_[N5_-2]*prms_.es_[N5_-2]),
	N5_-1);  

  for(int s=N5_-2; s>=0; --s){
    projM_(lf,w5,s+1);
    lf *= prms_.dm_[s];
    add5d(w5,lf,s);
    projP_(fy,w5,N5_-1);
    fy *= -prms_.mq_*prms_.fs_[s];
    add5d(w5,fy,s);
    Field v = get4d(w5,s);
    v*= 1.0/ prms_.dp_[s];
    set5d(w5,v,s);
  }
  return w5;
}

const Field DomainWallCore::mult_hop5_dinv(const Field& f5) const{
  Field w5(f5);
  Field lpf(f4size_),ey(f4size_),lmf(f4size_),fy(f4size_),v(f4size_);
  double* lpf_ptr = lpf.getaddr(0);
  double* ey_ptr = ey.getaddr(0);

  v = get4d(w5,0);
  v *= 1.0/prms_.dp_[0];
  set5d(w5,v,0);

  for(int s=1; s<N5_-1; ++s){
    projM_(lmf,w5,s-1);                                         
    lmf *= prms_.dm_[s-1];                                                    

    add5d(w5,lmf,s);                                                            
    v = get4d(w5,s);                                                            
    v*= 1.0/prms_.dp_[s];                                                    
 
    set5d(w5,v,s);                                                              
  }                                                                             
  projM_(v,w5,N5_-2);                                                  
  v*= prms_.dm_[N5_-2];                                                      
 
  add5d(w5,v,N5_-1);                                                            
  for(int s=0; s<N5_-1; ++s) {                                                  
    projP_(fy,w5,s);                                             
    fy *= -prms_.mq_*prms_.fs_[s];                                            
    add5d(w5,fy,N5_-1);                                                         
  }                                                                             
  v = get4d(w5,N5_-1);                                                          
  v*= 1.0/(prms_.dp_[N5_-1] +prms_.mq_*prms_.dm_[N5_-2]*prms_.es_[N5_-2]);    
  set5d(w5,v,N5_-1);                                                            
                                                                                
  for(int s=N5_-2; s>=0; --s){ 
    double* w5_ptr   = w5.getaddr(s*f4size_);
    double fact_lpf = (prms_.dm_[s+1]/prms_.dp_[s]);
    double fact_ey =  prms_.mq_*prms_.es_[s];

    projP_(lpf,w5,  s+1);                                          
    projM_( ey,w5,N5_-1); 

    for(int i=0; i<f4size_; ++i)
      w5_ptr[i] += fact_lpf*lpf_ptr[i]- fact_ey*ey_ptr[i];
  }
  return w5;
}

/*! @brief definitions of D_dwf */
void DomainWallCore::mult_full(Field& w5, const Field& f5) const{ 

  Field v(f4size_),lpf(f4size_), lmf(f4size_),w(f4size_);
  double* v_ptr   = v.getaddr(0);
  double* lpf_ptr = lpf.getaddr(0);
  double* lmf_ptr = lmf.getaddr(0);
  double* w_ptr = w.getaddr(0);
  double mass_fact= 4.0+prms_.M0_;

  for(int s=0; s<N5_; ++s){
    double* f5_ptr = const_cast<Field&>(f5).getaddr(s*f4size_);
    double* w5_ptr = w5.getaddr(s*f4size_);

    projP_(lpf,f5,(s+N5_-1)%N5_);
    if(s==0)     lpf *= -prms_.mq_;
    projM_(lmf,f5,(s+1)%N5_);
    if(s==N5_-1) lmf *= -prms_.mq_;

    for(int i=0; i<f4size_; ++i){
      lpf_ptr[i] += lmf_ptr[i];
      v_ptr[i] = prms_.bs_[s]*f5_ptr[i]+prms_.cs_[s]*lpf_ptr[i];
    }
    w = Dw_->mult(v);
    for(int i=0; i<f4size_; ++i)
      w5_ptr[i] = mass_fact*w_ptr[i]+ f5_ptr[i] -lpf_ptr[i];
  }
}

void DomainWallCore::mult_dag_full(Field& w5,const Field& f5) const{
  Field v5(f5.size());
  Field lpf(f4size_),lmf(f4size_),w(f4size_);

  for(int s=0; s< N5_; ++s){
    int spin_idx = s*f4size_;
    double* f5_ptr = const_cast<Field&>(f5).getaddr(spin_idx);
    double* w5_ptr = w5.getaddr(spin_idx);
    double* v5_ptr = v5.getaddr(spin_idx);
    
    double bs = (4.0+prms_.M0_)*prms_.bs_[s];
    double cs = (4.0+prms_.M0_)*prms_.cs_[s];

    Field w = Dw_->mult_dag(get4d(f5,s));
    double* w_ptr = w.getaddr(0);
    
    for(int i=0; i<f4size_; ++i){
      w5_ptr[i] = bs*w_ptr[i];
      v5_ptr[i] = cs*w_ptr[i];
    }
    // do not change this
    for(int i=0; i<f4size_; ++i){
      w5_ptr[i] += f5_ptr[i];
      v5_ptr[i] -= f5_ptr[i];
    }
  }
  for(int s = 0; s < N5_; ++s){
    projP_(lpf,v5,(s+1)%N5_);
    if(s == N5_-1) lpf *= -prms_.mq_;
    projM_(lmf,v5,(s+N5_-1)%N5_);
    if(s == 0)     lmf *= -prms_.mq_;
    add5d(w5,lpf,lmf,s);
  }
}

void DomainWallCore::mult_offdiag(Field& w5,const Field& f5) const{ 
  using namespace FieldExpression;
  Field lpf(f4size_),lmf(f4size_),w(f4size_),v(f4size_);
  double mass_fact= 4.0+prms_.M0_;

  for(int s=0; s<N5_; ++s) {
    projP_(lpf,f5,(s+N5_-1)%N5_);
    if(s==0)     lpf *= -prms_.mq_;
    projM_(lmf,f5,(s+1)%N5_);
    if(s==N5_-1) lmf *= -prms_.mq_;

    lpf += lmf;
    v = get4d(f5,s);
    v *= prms_.bs_[s];
    v += prms_.cs_[s]*lpf;
    w = Dw_->mult(v);
    w *= mass_fact;
    set5d(w5,w,s);
  }
}

void DomainWallCore::mult_dag_offdiag(Field& w5,const Field& f5) const{
  assert(w5.size()==f5.size());
  Field v5(f5.size()),lpf(f4size_),lmf(f4size_);

  for(int s=0; s<N5_; ++s){
    int spin_idx = s*f4size_;
    double* f5_ptr = const_cast<Field&>(f5).getaddr(spin_idx);
    double* w5_ptr = w5.getaddr(spin_idx);
    double* v5_ptr = v5.getaddr(spin_idx);
    
    double bs = (4.0+prms_.M0_)*prms_.bs_[s];
    double cs = (4.0+prms_.M0_)*prms_.cs_[s];
    
    Field w = Dw_->mult_dag(get4d(f5,s));
    double* w_ptr = w.getaddr(0);
    
    for (int i=0; i<f4size_; i++){
      w5_ptr[i] = bs*w_ptr[i];
      v5_ptr[i] = cs*w_ptr[i];
    }
  }
  for(int s=0; s<N5_; ++s){
    projP_(lpf,v5,(s+1)%N5_);
    if(s == N5_-1) lpf *= -prms_.mq_;
    projM_(lmf,v5,(s+N5_-1)%N5_);
    if(s == 0)     lmf *= -prms_.mq_;
    add5d(w5,lpf,lmf,s);
  }
}

/*! @brief contribution to the MD-force from forward difference */
void DomainWallCore::
md_force_p(Field& fce,const Field& phi,const Field& psi)const{
  using namespace FieldExpression;

  Field lpf(f4size_), lmf(f4size_);
  Field xie(f4size_);
  Field w(f4size_);
  double* lpf_ptr  = lpf.getaddr(0);
  double* lmf_ptr  = lmf.getaddr(0);
  double* xie_ptr  = xie.getaddr(0);
  
  for(int s=0; s<N5_; ++s){
    projP_(lpf,phi,(s+N5_-1)%N5_);
    if(s == 0)     lpf *= -prms_.mq_;
    projM_(lmf,phi,(s+1    )%N5_);
    if(s == N5_-1) lmf *= -prms_.mq_;

    Field w = get4d(phi,s); 
    
    w *= prms_.bs_[s];
    w += prms_.cs_[s]*(lpf +lmf);
    
    Dw_->md_force_p(fce,w,get4d(psi,s));  
  }
}  

/*! @brief contribution to the MD-force from backward difference */
void DomainWallCore::
md_force_m(Field& fce,const Field& phi,const Field& psi)const{
  using namespace FieldExpression;
  Field lpf(f4size_), lmf(f4size_);

  for(int s=0; s<N5_; ++s){
    projP_(lpf,phi,(s+N5_-1)%N5_);
    if(s == 0)     lpf *= -prms_.mq_;
    projM_(lmf,phi,(s+1    )%N5_);
    if(s == N5_-1) lmf *= -prms_.mq_;

    Field w = get4d(phi,s); 

    w *= prms_.bs_[s];
    w += prms_.cs_[s]*(lpf +lmf);

    Dw_->md_force_m(fce,w,get4d(psi,s));
  }
}  








