/*!
 * @file common_fields.cpp
 * @brief function definitions in FieldUtils
 */
#include "include/common_fields.hpp"
#include "Tools/sunMatUtils.hpp"

namespace FieldUtils{
  const GaugeField TracelessAntihermite(const GaugeField& G){
    using namespace SUNmatUtils;
    int Ndim = CommonPrms::instance()->Ndim();
    int Nvol = CommonPrms::instance()->Nvol();

    GaugeField TAField;
    for(int mu=0; mu< Ndim; ++mu){
      for(int site=0; site<Nvol; ++site){
	SetMatrix(TAField, anti_hermite(matrix(G,site,mu)),site,mu);
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

  void SetMatrix(GaugeField& F, SUNmat mat, int site, int dir){
    F.data.set(F.format.cslice(0,site,dir), mat.getva());
  }
  void SetMatrix(GaugeField1D& F, SUNmat mat, int site){
    F.data.set(F.format.cslice(0,site), mat.getva());
  }

  void AddMatrix(GaugeField& F, SUNmat mat, int site, int dir){
    F.data.add(F.format.cslice(0,site,dir), mat.getva());
  }
  void AddMatrix(GaugeField1D& F, SUNmat mat, int site){
    F.data.add(F.format.cslice(0,site), mat.getva());
  }

}
