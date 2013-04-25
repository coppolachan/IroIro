/*!
  @file siteMap.cpp
  @brief Definition of IndexOp,IndexOp_eo and IndexOop_oe classes
*/
#include "siteMap.hpp"
#include <algorithm>

namespace SiteMap{
// arrays of the function pointers 
// IndexOp
  const IndexOp::SiteOp IndexOp::xc[] ={
    &SiteIndex::c_x, &SiteIndex::c_y, 
    &SiteIndex::c_z, &SiteIndex::c_t,};

  const IndexOp::SiteOp IndexOp::xg[] ={
    &SiteIndex::g_x, &SiteIndex::g_y, 
    &SiteIndex::g_z, &SiteIndex::g_t,};
  
  const IndexOp::SiteOp IndexOp::xp[] ={
    &SiteIndex::p_x, &SiteIndex::p_y, 
    &SiteIndex::p_z, &SiteIndex::p_t,};

  const IndexOp::SiteOp IndexOp::xm[] ={
    &SiteIndex::m_x, &SiteIndex::m_y, 
    &SiteIndex::m_z, &SiteIndex::m_t,};

  const IndexOp::SliceOp IndexOp::sl[] ={
    &SiteIndex::slice_x, &SiteIndex::slice_y, 
    &SiteIndex::slice_z, &SiteIndex::slice_t,};

  const std::vector<int> IndexOp::xslice_map(int x,int dir)const{
    std::vector<int> map;
    for(int n=0; n<slice_size(x,dir);++n) 
      map.push_back(xslice(x,n,dir));
    return map;
  }

  const std::vector<int> IndexOp::bdry_map(int dir,TopBtm tb)const{
    int bd = (tb==Top)? SiteIndex::instance()->Bdir(dir): 0;
    return xslice_map(bd,dir);
  }

  const std::vector<int> IndexOp::bulk_map(int dir,TopBtm tb)const{
    int ini = 0;
    int fin = SiteIndex::instance()->Bdir(dir);
    if(tb==Btm){
      ini++;
      fin++;
    }
    std::vector<int> bulk;
    for(int x=ini; x<fin; ++x)
      for(int n=0; n<slice_size(x,dir);++n) 
	bulk.push_back(xslice(x,n,dir));
    //sort(bulk.begin(),bulk.end());
    return bulk;
  }

  // IndexOp_oe
  const IndexOp_oe::SiteOp IndexOp_oe::xc[] ={
    &SiteIndex_EvenOdd::c_xo,&SiteIndex_EvenOdd::c_y, 
    &SiteIndex_EvenOdd::c_z, &SiteIndex_EvenOdd::c_t,};

  const IndexOp_oe::SiteOp IndexOp_oe::xp[] ={
    &SiteIndex_EvenOdd::p_xoe,&SiteIndex_EvenOdd::p_y, 
    &SiteIndex_EvenOdd::p_z,  &SiteIndex_EvenOdd::p_t,};

  const IndexOp_oe::SiteOp IndexOp_oe::xm[] ={
    &SiteIndex_EvenOdd::m_xoe,&SiteIndex_EvenOdd::m_y, 
    &SiteIndex_EvenOdd::m_z,  &SiteIndex_EvenOdd::m_t,};
  
  const IndexOp_oe::SiteOp IndexOp_oe::slsize[] ={
    &SiteIndex_EvenOdd::slsize_xo,&SiteIndex_EvenOdd::slsize_y,
    &SiteIndex_EvenOdd::slsize_z, &SiteIndex_EvenOdd::slsize_t,};

  const IndexOp_oe::SliceOp IndexOp_oe::sl[] ={
    &SiteIndex_EvenOdd::slice_xo,&SiteIndex_EvenOdd::slice_y, 
    &SiteIndex_EvenOdd::slice_z, &SiteIndex_EvenOdd::slice_t,};

  const std::vector<int> IndexOp_oe::xslice_map(int x,int dir)const{
    std::vector<int> map;
    for(int n=0; n<slice_size(x,dir);++n) 
      map.push_back(xslice(x,n,dir));
    return map;
  }

  const std::vector<int> IndexOp_oe::bdry_map(int dir,TopBtm tb)const{
    int bd = (tb==Top)? SiteIndex_EvenOdd::instance()->Bdir(dir): 0;
    return xslice_map(bd,dir);
  }

  const std::vector<int> IndexOp_oe::bulk_map(int dir,TopBtm tb)const{
    int ini = 0;
    int fin = SiteIndex_EvenOdd::instance()->Bdir(dir);
    if(tb==Btm){
      ini++;
      fin++;
    }
    std::vector<int> bulk;
    for(int x=ini; x<fin; ++x)
      for(int n=0; n<slice_size(x,dir);++n) 
	bulk.push_back(xslice(x,n,dir));
    //sort(bulk.begin(),bulk.end());
    return bulk;
  }

  // IndexOp_eo
  const IndexOp_eo::SiteOp IndexOp_eo::xc[] ={
    &SiteIndex_EvenOdd::c_xe,&SiteIndex_EvenOdd::c_y, 
    &SiteIndex_EvenOdd::c_z, &SiteIndex_EvenOdd::c_t,};

  const IndexOp_eo::SiteOp IndexOp_eo::xp[] ={
    &SiteIndex_EvenOdd::p_xeo,&SiteIndex_EvenOdd::p_y, 
    &SiteIndex_EvenOdd::p_z,  &SiteIndex_EvenOdd::p_t,};

  const IndexOp_eo::SiteOp IndexOp_eo::xm[] ={
    &SiteIndex_EvenOdd::m_xeo,&SiteIndex_EvenOdd::m_y, 
    &SiteIndex_EvenOdd::m_z,  &SiteIndex_EvenOdd::m_t,};

  const IndexOp_eo::SiteOp IndexOp_eo::slsize[] ={
    &SiteIndex_EvenOdd::slsize_xe,&SiteIndex_EvenOdd::slsize_y,
    &SiteIndex_EvenOdd::slsize_z, &SiteIndex_EvenOdd::slsize_t,};

  const IndexOp_eo::SliceOp IndexOp_eo::sl[] ={
    &SiteIndex_EvenOdd::slice_xe,&SiteIndex_EvenOdd::slice_y, 
    &SiteIndex_EvenOdd::slice_z, &SiteIndex_EvenOdd::slice_t,};

  const std::vector<int> IndexOp_eo::xslice_map(int x,int dir)const{
    std::vector<int> map;
    for(int n=0; n<slice_size(x,dir);++n) 
      map.push_back(xslice(x,n,dir));
    return map;
  }

  const std::vector<int> IndexOp_eo::bdry_map(int dir,TopBtm tb)const{
    int bd = (tb==Top)? SiteIndex_EvenOdd::instance()->Bdir(dir): 0;
    return xslice_map(bd,dir);
  }

  const std::vector<int> IndexOp_eo::bulk_map(int dir,TopBtm tb)const{
    int ini = 0;
    int fin = SiteIndex_EvenOdd::instance()->Bdir(dir);
    if(tb==Btm){
      ini++;
      fin++;
    }
    std::vector<int> bulk;
    for(int x=ini; x<fin; ++x)
      for(int n=0; n<slice_size(x,dir);++n) 
	bulk.push_back(xslice(x,n,dir));
    //sort(bulk.begin(),bulk.end());
    return bulk;
  }

  // global instances
  IndexOp shiftSite;
  IndexOp_eo shiftSite_eo;
  IndexOp_oe shiftSite_oe;
}
