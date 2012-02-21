/*!
 * @file common_fields.cpp
 * @brief function definitions in FieldUtils
 */
#include "include/common_fields.hpp"
#include "Tools/sunMatUtils.hpp"

namespace FieldUtils{
  const Field TracelessAntihermite(const GaugeField& fce){
    using namespace SUNmat_utils;
    
    int Ndim = CommonPrms::instance()->Ndim();
    int Nvol = CommonPrms::instance()->Nvol();
    
    Field fp(fce.Format.size());
    for(int mu=0; mu< Ndim; ++mu){
      for(int site=0; site<Nvol; ++site){
	SUNmat a_h = u(fce,site,mu);
	fp.set(fce.Format.cslice(0,site,mu),
	       anti_hermite_traceless(a_h));
      }
    }
    return fp;
  }
  
  const Field field_oprod(const Field& f1,const Field& f2){
    using namespace SUNmat_utils;

    int Nd = CommonPrms::instance()->Nd();
    int Nvol = CommonPrms::instance()->Nvol();

    Format::Format_F ff(Nvol);
    Format::Format_G gf(Nvol,1);
    Field f(gf.size());
    
    for(int site=0; site<Nvol; ++site){
      SUNmat mat(0.0);
      for(int s=0; s<Nd; ++s){
	mat += outer_prod(SUNvec(f1[ff.cslice(s,site)]),
			  SUNvec(f2[ff.cslice(s,site)]));
      }
      f.set(gf.islice(site),mat.getva());
    }
    return f;
  }
}
