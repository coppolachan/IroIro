/*!
 * @file common_fields.cpp
 * @brief function definitions in FieldUtils
 */
#include "include/common_fields.hpp"
#include "Tools/sunMatUtils.hpp"

namespace FieldUtils{
  const Field field_oprod(const Field& f1,const Field& f2){
    using namespace SUNmatUtils;

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
    return GaugeField1D(Field(F.data[F.format.dir_slice(dir)]));
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
