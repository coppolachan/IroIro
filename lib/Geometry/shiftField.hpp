/*!
  @file shiftField.hpp
  @brief Declares ShiftField class for the shift
*/
#ifndef SHIFTFIELD_INCLUDED
#define SHIFTFIELD_INCLUDED

#include "mapping.hpp"

namespace Mapping{
  
  class ShiftField{
    std::vector<AutoMap> maps_;
  public:
    template<typename DATA,typename FB>
    DATA operator()(const DATA& F,int dir,FB fb)const{return maps_[dir](F,fb);}

    template <typename FB>
    void operator()(GaugeField1D& Fout,const double* Fin,int dir,FB fb)const{
      maps_[dir](Fout,Fin,fb); }

    void init_maps();
  };

  class ShiftField_eo{
    std::vector<AutoMap_EvenOdd> maps_;
  public:
    template<typename DATA,typename FB>
    DATA operator()(const DATA& F,int dir,FB fb)const{return maps_[dir](F,fb);}

    template <typename FB>
    void operator()(GaugeField1D& Fout,const double* Fin,int dir,FB fb)const{
      maps_[dir](Fout,Fin,fb); }

    void init_maps();
  };

  class ShiftField_oe{
    std::vector<AutoMap_EvenOdd> maps_;
  public:
    template<typename DATA,typename FB>
    DATA operator()(const DATA& F,int dir,FB fb)const{return maps_[dir](F,fb);}

    template <typename FB>
    void operator()(GaugeField1D& Fout,const double* Fin,int dir,FB fb)const{
      maps_[dir](Fout,Fin,fb); }

    void init_maps();
  };

  // declaration of a global object
  extern ShiftField shiftField;
  void init_shiftField();
  
  extern ShiftField_eo shiftField_eo;
  extern ShiftField_oe shiftField_oe;
  void init_shiftField_EvenOdd();
}



#endif
