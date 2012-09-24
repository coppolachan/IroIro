/*!
 * @file common_fields.hpp
 * @brief Code for general fields initializations
 */
#ifndef COMMON_FIELDS_H_
#define COMMON_FIELDS_H_

#include "include/macros.hpp"
#include "include/geometry.hpp"
#include "include/field.h"
#include "include/format_A.h"
#include "include/format_G.h"
#include "include/format_F.h"
#include "include/format_S.h"
#include<cassert>

struct OneDimTag {};
struct FiveDim{
  int fifth;
  int LocalVol;
};

/*!@brief A Class to handle gauge fields
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a> 
 * 
 * This class handles the Format and Field class in order to bind them 
 * into a single self contained object that knows the details of storage
 */
template <typename DATA,typename FORMAT,typename TAG = NullType> 
class GeneralField {
public:
  FORMAT format;
  DATA data;

  explicit GeneralField(int LocalVol = CommonPrms::instance()->Nvol())
    :format(LocalVol),data(format.size()){}

  explicit GeneralField(const Field& Fin,
			int LocalVol = CommonPrms::instance()->Nvol())
    :format(LocalVol),data(Fin){ assert(Fin.size()== format.size());}

  GeneralField(const GeneralField& Fin):format(Fin.format),data(Fin.data){}

  GeneralField& operator=(const GeneralField& rhs){
    assert(format.size() == rhs.format.size());
    data = rhs.data;
    return *this;
  }
  GeneralField& operator=(const Field& rhs){
    assert(format.size() == rhs.size());
    data = rhs;
    return *this;
  }
  GeneralField& operator=(double rhs){
    data = rhs;
    return *this;
  }

  GeneralField& operator+=(const GeneralField& rhs){
    assert(format.size() == rhs.format.size());
    data += rhs.data;
    return *this;
  }
  GeneralField& operator+=(const Field& rhs){
    assert(format.size() == rhs.size());
    data += rhs;
    return *this;
  }

  GeneralField& operator-=(const GeneralField& rhs){
    assert(format.size() == rhs.format.size());
    data -= rhs.data;
    return *this;
  }
  GeneralField& operator-=(const Field& rhs){
    assert(format.size() == rhs.size());
    data -= rhs;
    return *this;
  }

  GeneralField& operator*=(double rhs){
    data *= rhs;
    return *this;
  }

  GeneralField& operator/=(double rhs){
    data /= rhs;
    return *this;
  }

  int size()const { return data.size(); }
  int Nin()const {return format.Nin();}
  int Nvol()const {return format.Nvol();}
  int Nex()const {return format.Nex();}

  double norm(){ return data.norm();}

  const FORMAT get_Format() const{ return format;}
  const double operator[](int i)const {return data[i];}
  
  std::valarray<size_t> get_sub(const std::vector<int>& map,int ex)const{
    return format.get_sub(map,ex);}
  std::valarray<size_t> get_sub(const std::vector<int>& map)const{
    return format.get_sub(map);}
  const std::valarray<double> getva()const{return data.getva();}
};

//// list of typedefs ////
typedef GeneralField<Field,Format::Format_A>              AdjGaugeField;
typedef GeneralField<Field,Format::Format_G>              GaugeField;
typedef GeneralField<Field,Format::Format_G,OneDimTag>    GaugeField1D;
typedef GeneralField<Field,Format::Format_F,OneDimTag>    FermionField;
typedef GeneralField<Field,Format::Format_S,OneDimTag>    FermionField1sp;
typedef GeneralField<std::vector<Field>,Format::Format_F> PropagatorField;

template <> 
inline GaugeField1D::GeneralField(int LocalVol)
  :format(LocalVol,1),data(format.size()){}

/*! Since Nex is fixed to be 1 in this case, Nvol can be determined 
  as Fin.size()/Format::Format_G::Nin(), then e/o-field is also possible */
template <> 
inline GaugeField1D::GeneralField(const Field& Fin,int LocalVol)
  :format(Fin.size()/Format::Format_G::Nin(),1),data(Fin){}
template <> 
inline FermionField::GeneralField(const Field& Fin,int LocalVol)
  :format(Fin.size()/Format::Format_F::Nin(),1),data(Fin){}
template <> 
inline FermionField1sp::GeneralField(const Field& Fin,int LocalVol)
  :format(Fin.size()/Format::Format_S::Nin(),1),data(Fin){}

#endif //COMMON_FIELDS_H_
