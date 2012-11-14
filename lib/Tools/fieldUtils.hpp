/*! @file fieldUtils.hpp
    @brief header file for utility functions for GeneralField
 */
#ifndef FIELDUTILS_INCLUDED
#define FIELDUTILS_INCLUDED

#include "sunVec.hpp"
#include "include/common_fields.hpp"
#include "Main/Geometry/siteIndex_EvenOdd.hpp"

class RandNum;

namespace FieldUtils{
  const GaugeField TracelessAntihermite(const GaugeField&);
  const GaugeField1D TracelessAntihermite(const GaugeField1D&);
  const GaugeField ReUnit(const GaugeField&);
  const GaugeField1D ReUnit(const GaugeField1D&);
  const GaugeField Exponentiate(const GaugeField&, const double, const int);

  const GaugeField1D field_oprod(const FermionField&,const FermionField&);

  void RandomGauge(GaugeField1D& G,const RandNum& rand);
  void RandomGauge(GaugeField& G,const RandNum& rand);

  // Field type-type transformations
  GaugeField1D DirSlice(const GaugeField& F, int dir);
#ifdef IBM_BGQ_WILSON
  void DirSliceBGQ(GaugeField1D &G, const GaugeField& F, int dir);
#endif

  void SetSlice(GaugeField&, const GaugeField1D&, int dir);
  void AddSlice(GaugeField&, const GaugeField1D&, int dir);

  void SetMat(GaugeField& F, const SUNmat& mat, int site, int dir);
  void SetMat(GaugeField1D& F, const SUNmat& mat, int site);
  
  void AddMat(GaugeField& F, const SUNmat& mat, int site, int dir);
  void AddMat(GaugeField1D& F, const SUNmat& mat, int site);

  void SetVec(FermionField&, const SUNvec&, int spin, int site); 
  void AddVec(FermionField&, const SUNvec&, int spin, int site); 

  // Inline functions
  inline SUNmat mat(const GaugeField& F,int site,int dir){
    return SUNmat(F.data[F.format.islice(site,dir)]);  }

  inline SUNmat mat(const GaugeField1D& F,int site){
    return SUNmat(F.data[F.format.islice(site)]);  }

  inline SUNmat mat_dag(const GaugeField& F,int site,int dir){
    return SUNmat(F.data[F.format.islice(site,dir)]).dag();  }

  inline SUNmat mat_dag(const GaugeField1D& F,int site){
    return SUNmat(F.data[F.format.islice(site)]).dag();  }

  inline SUNvec vec(const FermionField& F,int spin,int site){
    return SUNvec(F.data[F.format.cslice(spin,site)]);  }

  inline SUNvec vec(const FermionField1sp& F,int site){
    return SUNvec(F.data[F.format.islice(site)]);  }

  // function templates
  template<typename GF>
  GF get_even(const GF& Fin){
    GF Fout(Fin.Nvol()/2);
    Fout.data = Fin.data[Fin.get_sub(SiteIndex_EvenOdd::instance()->esec())];
    return Fout;
  }

  template<typename GF>
  GF get_odd(const GF& Fin){
    GF Fout(Fin.Nvol()/2);
    Fout.data = Fin.data[Fin.get_sub(SiteIndex_EvenOdd::instance()->osec())];
    return Fout;
  }

  template<typename GF>
  GF combine_eo(const GF& Fe,const GF& Fo){
    GF Fout(Fe.Nvol()+Fo.Nvol());
    Fout.data.set(Fout.get_sub(SiteIndex_EvenOdd::instance()->esec()),Fe.getva());
    Fout.data.set(Fout.get_sub(SiteIndex_EvenOdd::instance()->osec()),Fo.getva());
    return Fout;
  }
}

#endif
