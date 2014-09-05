/*! @file dirac_wilson_FiniteDensity.cpp
 * @brief memberfuncs of Dirac_Wilson_FiniteDensity class
 Time-stamp: <2014-08-29 20:52:48 noaki>
*/

#include "dirac_wilson_FiniteDensity.hpp"
using namespace std;

/// initialization of the generalized slices for the absolute boundaries
void Dirac_Wilson_FiniteDensity::init_Tbdry(){
  ffmt_t ff(CommonPrms::Nvol());
  
  if(Communicator::instance()->ipe(TDIR) == CommonPrms::NPEt()-1){
    vector<int> sl;
    int Nt = CommonPrms::Nt();

    for(int z=0; z<CommonPrms::Nz(); ++z)
      for(int y=0; y<CommonPrms::Ny(); ++y)
	for(int x=0; x<CommonPrms::Nx(); ++x)
	  sl.push_back(SiteIndex::instance()->site(x,y,z,Nt-1));
    sTmax_.resize(ff.size()/Nt);    
    sTmax_= ff.get_sub(sl);
  }

  if(Communicator::instance()->ipe(TDIR) == 0){
    vector<int> sl;
    int Nt = CommonPrms::Nt();
    
    for(int z=0; z<CommonPrms::Nz(); ++z)
      for(int y=0; y<CommonPrms::Ny(); ++y)
	for(int x=0; x<CommonPrms::Nx(); ++x)
	  sl.push_back(SiteIndex::instance()->site(x,y,z,0));
    sTmin_.resize(ff.size()/Nt);    
    sTmin_= ff.get_sub(sl);
  }
}

const Field Dirac_Wilson_FiniteDensity::mult(const Field& f)const{
  Field Df = Dw_.mult(f);

  Field fb(fsize());
  fb.set(sTmin_,f[sTmin_]);         
  fb *= (1.0-fgcty_)*Dw_.getKappa();
  (Dw_.*(Dw_.mult_p[TDIR]))(Df,fb); /// t = Lt-1 -> 0

  fb = 0.0;
  fb.set(sTmax_,f[sTmax_]);        
  fb *= (1.0-1.0/fgcty_)*Dw_.getKappa();
  (Dw_.*(Dw_.mult_m[TDIR]))(Df,fb); /// t = 0 -> Lt-1
  
  return Df;
}

const Field Dirac_Wilson_FiniteDensity::mult_Ds(const Field& f)const{
  /*
  Field Df = Dw_.mult(f);
  Df -= f;
  Df /= -Dw_.getKappa();
  (Dw_.*(Dw_.mult_p[TDIR]))(Df,f); 
  (Dw_.*(Dw_.mult_m[TDIR]))(Df,f);
  */
  Field Df(fsize());
  for(int d=0; d<CommonPrms::Ndim()-1;++d){
    (Dw_.*(Dw_.mult_p[d]))(Df,f); 
    (Dw_.*(Dw_.mult_m[d]))(Df,f);
  }
  
  Df *= -Dw_.getKappa();
  return Df;
}

const Field Dirac_Wilson_FiniteDensity::mult_Dtp(const Field& f)const{
  Field Df(fsize());
  Field fb(fsize());
  fb.set(sTmin_,f[sTmin_]);         
  fb *= (-1.0+fgcty_);
  fb += f; 
  (Dw_.*(Dw_.mult_p[TDIR]))(Df,fb); 

  fb = 0.0;
  fb.set(sTmax_,f[sTmax_]);        
  fb *= (-1.0+1.0/fgcty_);
  fb += f;
  (Dw_.*(Dw_.mult_m[TDIR]))(Df,fb);

  Df *= -Dw_.getKappa();
  return Df;
}

const Field Dirac_Wilson_FiniteDensity::mult_Dtm(const Field& f)const{
  Field Df(fsize());
  Field fb(fsize());
  fb.set(sTmin_,f[sTmin_]);         
  fb *= (-1.0+fgcty_);
  fb += f; 
  (Dw_.*(Dw_.mult_p[TDIR]))(Df,fb); 

  fb = 0.0;
  fb.set(sTmax_,f[sTmax_]);        
  fb *= (1.0-1.0/fgcty_);
  fb -= f;
  (Dw_.*(Dw_.mult_m[TDIR]))(Df,fb);
  
  Df *= -Dw_.getKappa();
  return Df;
}

const Field Dirac_Wilson_FiniteDensity::mult_dag(const Field& f)const{
  Field g5f = Dw_.gamma5(f);
  Field Df = Dw_.mult(g5f);

  Field fb(fsize());
  fb.set(sTmin_,g5f[sTmin_]);         
  fb *= (1.0-1.0/fgcty_)*Dw_.getKappa();
  (Dw_.*(Dw_.mult_p[TDIR]))(Df,fb); /// t = Lt-1 -> 0

  fb = 0.0;
  fb.set(sTmax_,g5f[sTmax_]);        
  fb *= (1.0-fgcty_)*Dw_.getKappa();
  (Dw_.*(Dw_.mult_m[TDIR]))(Df,fb); /// t = 0 -> Lt-1

  return Dw_.gamma5(Df);
}

const Field Dirac_Wilson_FiniteDensity::
md_force(const Field& eta,const Field& zeta)const{

  Field fce = Dw_.md_force(eta,zeta);

  Field eb(fsize()), zb(fsize());
  eb.set(sTmin_, eta[sTmin_]);        
  zb.set(sTmax_,zeta[sTmax_]);       
  eb *= (1.0 -fgcty_)*Dw_.getKappa();
  Dw_.mkfrc(fce,eb,zb,TDIR);

  eb = 0.0; zb = 0.0;
  eb.set(sTmax_,gamma5( eta)[sTmax_]); 
  zb.set(sTmin_,gamma5(zeta)[sTmin_]); 
  eb *= (1.0 -1.0/fgcty_)*Dw_.getKappa();
  Dw_.mkfrc(fce,zb,eb,TDIR);

  return fce;
}
  
#ifdef IBM_BGQ_WILSON
void Dirac_Wilson_FiniteDensity::
mult_ptr(double* w,double* const f)const{

  Dw_.mult_ptr(w,f);

  valarray<double> fb(fsize());
  for(int i=0; i<sTmin_.size(); ++i) fb[sTmin_[i]] = f[sTmin_[i]];

  BGWilson_MultAdd_Dir(w,w,const_cast<Field*>(u_)->getaddr(0),&fb[0],
		       (1.0-fgcty_)*Dw_.getKappa(),
		       BGWILSON_DIRAC,BGWILSON_T,BGWILSON_FORWARD);
  fb = 0.0;
  for(int i=0; i<sTmax_.size(); ++i) fb[sTmax_[i]] = f[sTmax_[i]];

  BGWilson_MultAdd_Dir(w,w,const_cast<Field*>(u_)->getaddr(0),&fb[0],
		       (1.0-1.0/fgcty_)*Dw_.getKappa(),
		       BGWILSON_DIRAC,BGWILSON_T,BGWILSON_BACKWARD);
}

void Dirac_Wilson_FiniteDensity::
mult_dag_ptr(double* w,double* const f)const{

  int Nvol= CommonPrms::Nvol();
  valarray<double> g5f(fsize());

  BGWilsonLA_MultGamma5(&g5f[0],f,Nvol);
  
  Dw_.mult_ptr(w,&g5f[0]);

  valarray<double> fb(fsize());
  for(int i=0; i<sTmin_.size(); ++i) fb[sTmin_[i]] = g5f[sTmin_[i]];

  BGWilson_MultAdd_Dir(w,w,const_cast<Field*>(u_)->getaddr(0),&fb[0],
		       (1.0-1.0/fgcty_)*Dw_.getKappa(),
		       BGWILSON_DIRAC,BGWILSON_T,BGWILSON_FORWARD);
  fb = 0.0;
  for(int i=0; i<sTmax_.size(); ++i) fb[sTmax_[i]] = g5f[sTmax_[i]];

  BGWilson_MultAdd_Dir(w,w,const_cast<Field*>(u_)->getaddr(0),&fb[0],
		       (1.0-fgcty_)*Dw_.getKappa(),
		       BGWILSON_DIRAC,BGWILSON_T,BGWILSON_BACKWARD);

  BGWilsonLA_MultGamma5(&g5f[0],w,Nvol);
  BGWilsonLA_Equate(w,&g5f[0],Nvol);
}

#endif
