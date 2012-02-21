//---------------------------------------------------------------------
/*! @file sunMatUtils.hpp
  @brief \f$SU(N)\f$ Matrices linear algebra Utilities

  Class declarations
*/ 
//---------------------------------------------------------------------

#ifndef SUNMAT_UTILS_H_
#define SUNMAT_UTILS_H_

#include "sunMat.h"
#include "sunVec.h"
#include "include/common_fields.hpp"
#include "lib/Main/Geometry/shiftField.h"

namespace SUNmat_utils{
  SUNmat unity();
  SUNmat zero();
  double ReTr(const SUNmat&);  
  double ImTr(const SUNmat&);
  const SUNmat dag(const SUNmat&);  
  const SUNmat xI(const SUNmat&);  
  const SUNmat operator+(const SUNmat&, const SUNmat&);
  const SUNmat operator-(const SUNmat&, const SUNmat&);
  const SUNmat operator*(const SUNmat&, const SUNmat&);
  const SUNmat reunit(const SUNmat&);
  const SUNmat outer_prod(const SUNvec&, const SUNvec&);
  const std::valarray<double> trace_less(const SUNmat&);
  const std::valarray<double> anti_hermite(const SUNmat&);
  const std::valarray<double> anti_hermite_traceless(const SUNmat&);

//////////////////////////////////////////////////////////////////////////////

  // for variables with the direction unfixed 
  inline SUNmat u(const Field& g,const Format::Format_G& gf_, 
		  int site,int dir){
    return SUNmat(g[gf_.cslice(0,site,dir)]);
  }
  inline SUNmat u(const GaugeField& g, int site,int dir){
    return SUNmat(g.U[g.Format.cslice(0,site,dir)]);
  }
  inline SUNmat u_dag(const Field& g,const Format::Format_G& gf_,
		      int site,int dir){
    return SUNmat(g[gf_.cslice(0,site,dir)]).dag();
  }
  inline SUNmat u_dag(const GaugeField& g, int site,int dir){
    return SUNmat(g.U[g.Format.cslice(0,site,dir)]).dag();
  }

  // for variables with a specific direction
  inline SUNmat u(const Field& g,const Format::Format_G& sf_,int site){
    return SUNmat(g[sf_.cslice(0,site)]);
  }
  inline SUNmat u(const GaugeField1D& g,int site){
    return SUNmat(g.U[g.Format.cslice(0,site)]);
  }

  inline SUNmat u(const std::valarray<double>& vu,
		  const Format::Format_G& sf_,int site){
    return SUNmat(vu[sf_.cslice(0,site)]);
  }
  inline SUNmat u(const ShiftField& su,int site){
    return SUNmat(su.cv(0,site));
  }
  inline SUNmat u_dag(const Field& g,const Format::Format_G& sf_,int site){
    return SUNmat(g[sf_.cslice(0,site)]).dag();
  }
  inline SUNmat u_dag(const GaugeField1D& g,int site){
    return SUNmat(g.U[g.Format.cslice(0,site)]).dag();
  }
  inline SUNmat u_dag(const std::valarray<double>& vu,
		      const Format::Format_G& sf_,int site){
    return SUNmat(vu[sf_.cslice(0,site)]).dag();
  }
  inline SUNmat u_dag(const ShiftField& su,int site){
    return SUNmat(su.cv(0,site)).dag();
  }

}//endof namespace SUNmat_utils

#endif
