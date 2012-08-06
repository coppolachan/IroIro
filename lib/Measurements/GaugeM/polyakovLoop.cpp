/*! @file polyakovLoop.cpp
 *  @brief implementation of polyakovLoop.cpp
 */
#include "polyakovLoop.hpp"
#include "Tools/sunAdjMatUtils.hpp"

using namespace std;

complex<double> PolyakovLoop::calc_SUN(const GaugeField& G)const{

  valarray<double> plp(0.0,SUNmat::size()*slsize_);
  calc<SUNmat>(plp,G);

  int Nin = SUNmat::size();
  double pl_re = 0.0;
  double pl_im = 0.0;

  for(int s=0; s<slsize_; ++s){
    std::slice ms(Nin*s,Nin,1);
    pl_re += SUNmatUtils::ReTr(SUNmat(plp[ms]));
    pl_im += SUNmatUtils::ImTr(SUNmat(plp[ms]));
  }
  pl_re/= double(slsize_*NC_);
  pl_im/= double(slsize_*NC_);

  int NP = CommonPrms::instance()->NP();
  return std::complex<double>(com_->reduce_sum(pl_re)/double(NP),
			      com_->reduce_sum(pl_im)/double(NP));
  /*
  // possible output (if the propagator of Ploop is wanted)
  valarray<double> plx(2*slsize_);
  for(int s=0; s<slsize_; ++s){ 
    slice ms = slice(Nin*s,Nin,1);
    plx[2*s  ] = ReTr(MAT(plp[ms]));
    plx[2*s+1] = ImTr(MAT(plp[ms]));

    int site = SiteMap:shiftSite.xslice(0,s,mu_dir_);             
    int gsite = SiteIndex::instance()->get_gsite(site); 
    int lx = SiteIndex::instance()->g_x(site);
    int ly = SiteIndex::instance()->g_y(site);
    int lz = SiteIndex::instance()->g_z(site);
    int lt = SiteIndex::instance()->g_t(site);
  }
  */
}

double PolyakovLoop::calc_SUNadj(const GaugeField& G)const{
  
  valarray<double> plp(0.0,SUNadjMat::size()*slsize_);
  calc<SUNadjMat>(plp,G);

  int Nin = SUNadjMat::size();
  double pla = 0.0;
  
  for(int s=0; s<slsize_; ++s){
    std::slice ms(Nin*s,Nin,1);
    pla += SUNadjMatUtils::Tr(SUNadjMat(plp[ms]));
  }
  pla/= double(slsize_*(NC_*NC_-1));

  int NP = CommonPrms::instance()->NP();
  return com_->reduce_sum(pla)/double(NP);
}
