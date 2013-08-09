/*! @file polyakovLoop.cpp
 *  @brief implementation of polyakovLoop.cpp
 */
#include "polyakovLoop.hpp"
#include "Tools/sunAdjMatUtils.hpp"

// For eigenvalues
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

using namespace std;

GaugeField1D PolyakovLoop::get_PLField(const GaugeField& G) const {
  using namespace FieldUtils;
  GaugeField1D PL_field;
  calc<SUNmat>(PL_field,G);
  return PL_field;
}



complex<double> PolyakovLoop::calc_SUN(const GaugeField& G)const{
  using namespace FieldUtils;
  GaugeField1D PL_field;
  calc<SUNmat>(PL_field,G);
  
  int Nin = SUNmat::size();
  double pl_re = 0.0;
  double pl_im = 0.0;

  for(int s=0; s<slsize_; ++s){
    SUNmat PL = mat(PL_field,s);
    pl_re += SUNmatUtils::ReTr(PL);
    pl_im += SUNmatUtils::ImTr(PL);
  }
  // Normalization
  pl_re/= double(slsize_*NC_);
  pl_im/= double(slsize_*NC_);

  // Global sum
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
  GaugeField1D PL_field;
  calc<SUNadjMat>(PL_field,G);  

  int Nin = SUNadjMat::size();
  double pla = 0.0;
  
  for(int s=0; s<slsize_; ++s){
    std::slice ms(Nin*s,Nin,1);
    pla += SUNadjMatUtils::Tr(SUNadjMat(PL_field.data[ms]));
  }
  pla/= double(slsize_*(NC_*NC_-1));

  int NP = CommonPrms::instance()->NP();
  return com_->reduce_sum(pla)/double(NP);
}
