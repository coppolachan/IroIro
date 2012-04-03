/*! @file qprop_mom.cpp
  @brief implementation of QpropMom class */

#include "Main/Geometry/siteIndex.hpp"
#include "include/numerical_const.hpp"
#include "include/macros.hpp"
#include "qprop_mom.hpp"

void QpropMom::init_mom(){
  std::vector<int> mom(NDIM_);

  int Mx = Lx_/2/PI;
  int My = Ly_/2/PI;
  int Mz = Lz_/2/PI;
  int Mt = Lt_/2/PI;

  for(int px=-Mx; px<=Mx; ++px){
    mom[XDIR] = px;
    for(int py=-My; py<=My; ++py){
      mom[YDIR] = py;
      for(int pz=-Mz; pz<=Mz; ++pz){
	mom[ZDIR] = pz;
	for(int pt=-Mt; pt<=Mt; ++pt){
	  mom[TDIR] = pt;
	  platt_.push_back(mom);
	}
      }
    }
  }
}

void QpropMom::fourier_tr(std::vector<double>& Sp_re,
			  std::vector<double>& Sp_im, 
			  const std::vector<int>& p) const{

  double p_x = 2.0*PI/Lx_*p[XDIR];
  double p_y = 2.0*PI/Ly_*p[YDIR];
  double p_z = 2.0*PI/Lz_*p[ZDIR];
  double p_t = 2.0*PI/Lt_*p[TDIR];

  CCIO::cout<<"fourier_tr  p_latt= "
	    <<p_x<<" "
	    <<p_y<<" "
	    <<p_z<<" "
	    <<p_t<<std::endl;

  for(int in=0; in<NC_*ND_; ++in){
    for(int s=0; s<ND_; ++s){
      for(int c=0; c<NC_; ++c){

	double re =0.0;
	double im =0.0;

	for(int site=0; site<fmt_.Nvol(); ++site){
	  int gsite = SiteIndex::instance()->get_gsite(site);
	  double pdotx 
	    = p_x*SiteIndex::instance()->g_x(gsite)
	     +p_y*SiteIndex::instance()->g_y(gsite)
 	     +p_z*SiteIndex::instance()->g_z(gsite)
	     +p_t*SiteIndex::instance()->g_t(gsite);

	  //	  CCIO::cout<<"pdotx="<<pdotx<<std::endl;

	  re += cos(pdotx)*Sq_[in][fmt_.index_r(c,s,site)]
       	       +sin(pdotx)*Sq_[in][fmt_.index_i(c,s,site)];
	  im +=-sin(pdotx)*Sq_[in][fmt_.index_r(c,s,site)]
	       +cos(pdotx)*Sq_[in][fmt_.index_i(c,s,site)];
	}
	Sp_re.push_back(Communicator::instance()->reduce_sum(re));
	Sp_im.push_back(Communicator::instance()->reduce_sum(im));
      }
    }
  }
}
  
void QpropMom::output()const{

  CCIO::cout<<"S(p_latt)[s1 c1; s2 c2], p_latt=(nx,ny,nz,nt)"<<std::endl;

  for(int m=0; m<platt_.size(); ++m){

    std::vector<double> Sp_re;
    std::vector<double> Sp_im;

    fourier_tr(Sp_re,Sp_im,platt_[m]);

    CCIO::cout<<"p_latt= "
	      <<platt_[m][XDIR]<<" "
	      <<platt_[m][YDIR]<<" "
	      <<platt_[m][ZDIR]<<" "
	      <<platt_[m][TDIR]<<std::endl;
    
    for(int s1=0; s1<ND_; ++s1){
      for(int c1=0; c1<NC_; ++c1){
	for(int s2=0; s2<ND_; ++s2){
	  for(int c2=0; c2<NC_; ++c2){

	    int sc_idx = (s1*NC_+c1)*NC_*ND_+s2*NC_+c2;

	    CCIO::cout<<s1<<" "<<c1<<" "<<s2<<" "<<c2<<"   "
		      <<Sp_re[sc_idx]<<"  "<<Sp_im[sc_idx]<<std::endl;
	  }
	}
      }
    }
  }
}
