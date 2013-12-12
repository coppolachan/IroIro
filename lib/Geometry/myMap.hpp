/*!
  @file myMap.hpp
  @brief Declares MyMap class
*/
#ifndef MYMAP_H_
#define MYMAP_H_

#include "siteIndex.hpp"
#include "siteIndex3d.hpp"

namespace SiteMap{

  template<typename IDX>
  class MyMap{
    IDX idx_;
    typedef int (IDX::*Xsl)(int x,int n)const;
    Xsl xsl[NDIM_];
  public:    
    MyMap(){
      xsl[XDIR] = &IDX::slice_x;
      xsl[YDIR] = &IDX::slice_y;
      xsl[ZDIR] = &IDX::slice_z;
      xsl[TDIR] = &IDX::slice_t;
    }
  public:
    const std::vector<int> bdry_map(int dir,TopBtm tb)const{
      int bd = (tb==Top)? idx_.Bdir(dir): 0;
      int slsize = idx_.slsize(dir);

      std::vector<int> map;
      for(int n=0; n<slsize; ++n) map.push_back((idx_.*xsl[dir])(bd,n));
      return map;
    }
    const std::vector<int> bulk_map(int dir,TopBtm tb)const{
      int ini = 0;
      int fin = idx_.Bdir(dir);

      if(tb==Btm){
	ini++;
	fin++;
      }
      std::vector<int> bulk;
      int slsize = idx_.slsize(dir);
      for(int x=ini; x<fin; ++x)
	for(int n=0; n<slsize; ++n) bulk.push_back((idx_.*xsl[dir])(x,n));
      return bulk;
    }
  };
  
  typedef MyMap<SiteIndex3d> Map3d;
}
#endif
