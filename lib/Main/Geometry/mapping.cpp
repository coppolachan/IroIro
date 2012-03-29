/*!
  @file mapping.cpp
  @brief Definition of ShiftField class
*/
#include "mapping.hpp"
#include "include/macros.hpp"

namespace Mapping{
  void ShiftField::init_maps(){
    if(maps_.size()==0){
      for(int dir=0; dir<NDIM_; ++dir)	
	maps_.push_back(AutoMap(dir));
    }
  }

  // global instance
  ShiftField shiftField;
  void init_shiftField(){ shiftField.init_maps();}
}

