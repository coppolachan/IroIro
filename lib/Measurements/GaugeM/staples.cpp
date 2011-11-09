//----------------------------------------------------------------------
// staples.cpp
//----------------------------------------------------------------------
#include "staples.h"

using namespace SUNmat_utils;
using namespace Format;
using namespace std;

// for variables with the direction unfixed 
inline SUNmat Staples::u(const Field& g,int site,int dir)const{
  return SUNmat(g[gf_.cslice(0,site,dir)]);
}
inline SUNmat Staples::u_dag(const Field& g,int site,int dir)const{
  return SUNmat(g[gf_.cslice(0,site,dir)]).dag();
}

// for variables with a specific direction
inline SUNmat Staples::u(const Field& g,int site)const{
  return SUNmat(g[sf_->cslice(0,site)]);
}
inline SUNmat Staples::u(const std::valarray<double>& vu,int site)const{
  return SUNmat(vu[sf_->cslice(0,site)]);
}
inline SUNmat Staples::u(const ShiftField& su,int site)const{
  return SUNmat(su.cv(0,site));
}
inline SUNmat Staples::u_dag(const Field& g,int site)const{
  return SUNmat(g[sf_->cslice(0,site)]).dag();
}
inline SUNmat Staples::u_dag(const std::valarray<double>& vu,int site)const{
  return SUNmat(vu[sf_->cslice(0,site)]).dag();
}
inline SUNmat Staples::u_dag(const ShiftField& su,int site)const{
  return SUNmat(su.cv(0,site)).dag();
}

double Staples::plaquette(const Field& g)const{ 
  return (plaq_s(g) +plaq_t(g))/2;
}
double Staples::plaquette(const ShiftField& gs)const{
  Field g = gs.field();
  return (plaq_s(g) +plaq_t(g))/2;
}

double Staples::plaq_s(const Field& g) const{
  double plaq = 0.0;
  Field stpl(sf_->size());

  for(int i=0;i<Ndim_-1;++i){
    int j = (i+1)%(Ndim_-1);
    
    stpl = lower(g,i,j);
    //cout<<"["<<j<<"]  stpl(0)="<<stpl[0]<<" stpl(1)="<<stpl[1]<<endl;

    for(int site=0; site<Nvol_; ++site)
      plaq += ReTr(u(g,site,i)*u_dag(stpl,site));  // P_ij
  }
  plaq = com_->reduce_sum(plaq);
  return plaq/(Lvol_*Nc_*3.0);
}

double Staples::plaq_s(const ShiftField& gs) const{
  Field g = gs.field();
  plaq_s(g);
}

double Staples::plaq_t(const Field& g)const{
  double plaq = 0.0;
  Field stpl(sf_->size());

  for(int nu=0; nu < Ndim_-1; ++nu){
    stpl = lower(g,3,nu);
    for(int site=0; site<Nvol_; ++site)
      plaq += ReTr(u(g,site,3)*u_dag(stpl,site));  // P_zx
  }
  plaq = com_->reduce_sum(plaq);
  return plaq/(Lvol_*Nc_*3.0);
}

double Staples::plaq_t(const ShiftField& gs)const{
  Field g = gs.field();
  plaq_t(g);
}

void Staples::staple(Field& W, const Field& g, int mu) const{
  W = 0.0;
  for(int nu = 0; nu < Ndim_; nu++){
    if(nu != mu){
      W += upper(g,mu,nu);
      W += lower(g,mu,nu);
    }
  }
}

void Staples::staple(Field& W, const ShiftField& gs, int mu)const {
  Field g = gs.field();
  staple(W,g,mu);
}

Field Staples::upper(const Field& g, int mu, int nu) const{
  //       mu,v                               
  //      +-->--+                                                    
  // nu,w |     |w_dag(site+mu,nu)
  //  site+     +                                                             

  valarray<double> w = g[gf_.dir_slice(nu)];
  valarray<double> v = g[gf_.dir_slice(mu)];
  ShiftField_up<Format_G> um(w,sf_,mu);
  ShiftField_up<Format_G> un(v,sf_,nu);

  valarray<double> c(sf_->size());
  for(int site=0; site<Nvol_; ++site){
    /*
    cout<<" nodeid="<< Communicator::instance()->nodeid()
	<<" w[0,"<<site <<"]="<<u(w, site).r(0)
	<<" un[0,"<<site<<"]="<<u(un,site).r(0)
      	<<" um[0,"<<site<<"]="<<u(um,site).r(0)
	<<endl;
    */
    /*
    SUNmat tmp1 = u(w,site)*u(un,site);
    cout<<" nodeid="<< Communicator::instance()->nodeid()
	<<" site="<<site
	<<" tmp[0]="<<(tmp1.getva())[0]
	<<" tmp[1]="<<(tmp1.getva())[1]
      	<<" tmp[2]="<<(tmp1.getva())[2]
	<<endl;
    */
    /*
    valarray<double> tmp=(tmp1*u_dag(um,site)).getva();
    cout<<" nodeid="<< Communicator::instance()->nodeid()
	<<" site="<<site
	<<" tmp[0]="<<tmp[0]
	<<" tmp[1]="<<tmp[1]
      	<<" tmp[2]="<<tmp[2]
	<<endl;
    */
    c[sf_->cslice(0,site)] = (u(w,site)*u(un,site)*u_dag(um,site)).getva();
  }
  return Field(c);
}
		 
Field Staples::upper(const ShiftField& gs, int mu,int nu) const{
  Field g = gs.field();
  return upper(g,mu,nu);
}

Field Staples::lower(const Field& g, int mu, int nu) const{
  //         +     +
  // nu,w_dag|     |w(site+mu,nu) 
  //     site+-->--+ 
  //           mu,v              

  valarray<double> v = g[gf_.dir_slice(mu)];
  valarray<double> w = g[gf_.dir_slice(nu)];

  ShiftField_up<Format_G> um(w,sf_,mu);

  valarray<double> c(sf_->size());
  for(int site=0; site<Nvol_; ++site)
    c[sf_->cslice(0,site)] = (u_dag(w,site)*u(v,site)*u(um,site)).getva();

  ShiftField_dn<Format_G> un(c,sf_,nu);

  for(int site=0; site<Nvol_; ++site)
    v[sf_->cslice(0,site)] = u(un,site).getva();
  return Field(v);
}

Field Staples::lower(const ShiftField& gs, int mu,int nu) const{
  Field g = gs.field();
  return lower(g,mu,nu);
}

