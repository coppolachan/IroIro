/*!
  @file siteMap.hpp
  @brief Declares IndexOp, IndexOp_eo and IndexOp_oe classes
*/
#ifndef SITEMAP_H_
#define SITEMAP_H_

#include "siteIndex.hpp"
#include "siteIndex_EvenOdd.hpp"

enum TopBtm {Top,Btm};

namespace SiteMap{
  class IndexOp{
    typedef int (SiteIndex::*SiteOp)(int hs) const;
    typedef int (SiteIndex::*SliceOp)(int x,int n) const;

    static const SiteOp xc[];
    static const SiteOp xp[];
    static const SiteOp xm[];
    static const SiteOp xg[];
    static const SliceOp sl[];
  public:
    int x_c(int x,int dir)const{return (SiteIndex::instance()->*xc[dir])(x);}
    int x_p(int x,int dir)const{return (SiteIndex::instance()->*xp[dir])(x);}
    int x_m(int x,int dir)const{return (SiteIndex::instance()->*xm[dir])(x);}
    int x_g(int x,int dir)const{return (SiteIndex::instance()->*xg[dir])(x);}
    int slice_size(int x,int dir)const{
      return SiteIndex::instance()->slsize(x,dir);}
    int xslice(int x,int n,int dir)const{
      return (SiteIndex::instance()->*sl[dir])(x,n);}

    const std::vector<int> xslice_map(int x,int dir)const;
    const std::vector<int> bdry_map(int dir,TopBtm tb)const;
    const std::vector<int> bulk_map(int dir,TopBtm tb)const;
  };
  
  class IndexOp_eo{
    typedef int (SiteIndex_EvenOdd::*SiteOp)(int)const;
    typedef int (SiteIndex_EvenOdd::*SliceOp)(int x,int n) const;

    static const SiteOp xc[];
    static const SiteOp xp[];
    static const SiteOp xm[];
    static const SiteOp xg[];
    static const SiteOp slsize[];
    static const SliceOp sl[];
  public:
    int x_c(int x,int dir)const {
      return (SiteIndex_EvenOdd::instance()->*xc[dir])(x);}
    int x_p(int x,int dir)const {
      return (SiteIndex_EvenOdd::instance()->*xp[dir])(x);}
    int x_m(int x,int dir)const {
      return (SiteIndex_EvenOdd::instance()->*xm[dir])(x);}
    int x_g(int x,int dir)const {
      return (SiteIndex_EvenOdd::instance()->*xg[dir])(x);}
    int slice_size(int x,int dir)const {
      return (SiteIndex_EvenOdd::instance()->*slsize[dir])(x);}
    int xslice(int x,int n,int dir)const {
      return (SiteIndex_EvenOdd::instance()->*sl[dir])(x,n);}

    const std::vector<int> xslice_map(int x,int dir)const;
    const std::vector<int> bdry_map(int dir,TopBtm tb)const;
    const std::vector<int> bulk_map(int dir,TopBtm tb)const;
  };

  class IndexOp_oe{
    typedef int (SiteIndex_EvenOdd::*SiteOp)(int)const;
    typedef int (SiteIndex_EvenOdd::*SliceOp)(int x,int n) const;

    static const SiteOp xc[];
    static const SiteOp xp[];
    static const SiteOp xm[];
    static const SiteOp xg[];
    static const SiteOp slsize[];
    static const SliceOp sl[];
  public:
    int x_c(int x,int dir)const {
      return (SiteIndex_EvenOdd::instance()->*xc[dir])(x);}
    int x_p(int x,int dir)const {
      return (SiteIndex_EvenOdd::instance()->*xp[dir])(x);}
    int x_m(int x,int dir)const {
      return (SiteIndex_EvenOdd::instance()->*xm[dir])(x);}
    int x_g(int x,int dir)const {
      return (SiteIndex_EvenOdd::instance()->*xg[dir])(x);}
    int slice_size(int x,int dir)const {
      return (SiteIndex_EvenOdd::instance()->*slsize[dir])(x);}
    int xslice(int x,int n,int dir)const {
      return (SiteIndex_EvenOdd::instance()->*sl[dir])(x,n);}

    const std::vector<int> xslice_map(int x,int dir)const;
    const std::vector<int> bdry_map(int dir,TopBtm tb)const;
    const std::vector<int> bulk_map(int dir,TopBtm tb)const;
  };

  // declaration of global objects
  extern IndexOp    shiftSite;
  extern IndexOp_eo shiftSite_eo;
  extern IndexOp_oe shiftSite_oe;
}
#endif //MANUALMAP_H_
