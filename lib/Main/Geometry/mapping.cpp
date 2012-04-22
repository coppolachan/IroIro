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

  void ShiftField_eo::init_maps(){
    if(maps_.size()==0){
      for(int dir=0; dir<NDIM_; ++dir)	
	maps_.push_back(AutoMap_EvenOdd(dir,EOtag()));
    }
  }

  void ShiftField_oe::init_maps(){
    if(maps_.size()==0){
      for(int dir=0; dir<NDIM_; ++dir)	
	maps_.push_back(AutoMap_EvenOdd(dir,OEtag()));
    }
  }

  // global instance
  ShiftField shiftField;
  void init_shiftField(){ shiftField.init_maps();}

  ShiftField_eo shiftField_eo;
  ShiftField_oe shiftField_oe;
  void init_shiftField_EvenOdd(){ 
    shiftField_eo.init_maps();
    shiftField_oe.init_maps();
  }
}

