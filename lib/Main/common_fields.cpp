/*!
 * @file common_fields.cpp
 * @brief function definitions in FieldUtils
 */
#include "include/common_fields.hpp"
#include "Tools/sunMatUtils.hpp"

namespace FieldUtils{
  const GaugeField1D field_oprod(const FermionField& f1,
				 const FermionField& f2){
    using namespace SUNmatUtils;

    GaugeField1D f;
    SUNmat mat;
    int Nd = CommonPrms::instance()->Nd();
    int Nvol = CommonPrms::instance()->Nvol();

    for(int site=0; site<Nvol; ++site){
      mat.zero();
      for(int s=0; s<Nd; ++s){
	mat += outer_prod(vect(f1, s, site), vect(f2, s, site));
      }
      SetMatrix(f, mat, site);
    }
    return f;
  }

  const GaugeField TracelessAntihermite(const GaugeField& G){
    using namespace SUNmatUtils;
    int Ndim = CommonPrms::instance()->Ndim();
    int Nvol = CommonPrms::instance()->Nvol();

    GaugeField TAField;
    for(int mu=0; mu< Ndim; ++mu){
      for(int site=0; site<Nvol; ++site){
	SetMatrix(TAField, anti_hermite_traceless(matrix(G,site,mu)),site,mu);
      }
    }
    return TAField;
  }

  GaugeField1D DirSlice(const GaugeField& F, int dir){
    #ifdef IBM_WATSON
    GaugeField1D out;
    out.data = F.data[F.format.dir_slice(dir)];
    return out;
    #else
    return GaugeField1D(Field(F.data[F.format.dir_slice(dir)]));
    #endif
  }

  void SetSlice(GaugeField& G, const GaugeField1D& Gslice, int dir){
    G.data.set(G.format.dir_slice(dir), Gslice.data.getva());
  }

  void AddSlice(GaugeField& G, const GaugeField1D& Gslice, int dir){
    G.data.add(G.format.dir_slice(dir), Gslice.data.getva());
  }

  void SetMatrix(GaugeField& F, const SUNmat& mat, int site, int dir){
    F.data.set(F.format.cslice(0,site,dir), mat.getva());
  }
  void SetMatrix(GaugeField1D& F, const SUNmat& mat, int site){
    F.data.set(F.format.cslice(0,site), mat.getva());
  }

  void AddMatrix(GaugeField& F, const SUNmat& mat, int site, int dir){
    F.data.add(F.format.cslice(0,site,dir), mat.getva());
  }
  void AddMatrix(GaugeField1D& F, const SUNmat& mat, int site){
    F.data.add(F.format.cslice(0,site), mat.getva());
  }

  void SetVector(FermionField& F, const SUNvec& vec, int spin, int site){
    F.data.set(F.format.cslice(spin, site), vec.getva());
  }
  void AddVector(FermionField& F, const SUNvec& vec, int spin, int site){
    F.data.add(F.format.cslice(spin, site), vec.getva());
  }

}
