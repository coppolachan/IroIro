/*!
 * @file common_fields.hpp
 * @brief Code for general fields initializations
 */
#ifndef COMMON_FIELDS_H_
#define COMMON_FIELDS_H_

#include "Main/gaugeConf.hpp"
#include "include/geometry.hpp"
#include "include/field.h"
#include "include/format_A.h"
#include "include/format_G.h"
#include "include/format_F.h"
#include "lib/Tools/sunMat.hpp"
#include "lib/Tools/sunVec.hpp"
#include "include/macros.hpp"

struct OneDimTag {};
//struct HalfSiteTag {};
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
template <class DATA,class FORMAT,typename TAG = NullType> 
class GeneralField {
public:
  FORMAT format;
  DATA data;
  GeneralField();
  GeneralField(int);
  
  const FORMAT get_Format() const{ return format;}

  explicit GeneralField(const DATA&);

  GeneralField& operator=(const GeneralField&);
  GeneralField& operator=(const Field&);
  GeneralField& operator=(const std::valarray<double>&);
  GeneralField& operator=(double);

  GeneralField& operator+=(const GeneralField&);
  GeneralField& operator+=(const Field&);
  GeneralField& operator+=(const std::valarray<double>&);

  GeneralField& operator-=(const GeneralField&);
  GeneralField& operator-=(const Field&);
  GeneralField& operator-=(const std::valarray<double>&);

  GeneralField& operator*=(double);
  GeneralField& operator/=(double);
  
  double norm();
  int size()const { return data.size(); }
  int Nin()const {return format.Nin();}
  int Nvol()const {return format.Nvol();}
  int Nex()const {return format.Nex();}
  const double operator[](int i)const {return data[i];}
  
  std::valarray<size_t> get_sub(const std::vector<int>& map)const{
    return format.get_sub(map);}
  const std::valarray<double> getva()const{return data.getva();}
};

// Constructor with no argument ///////////////////////
template <class DATA,class FORMAT,typename TAG> 
GeneralField<DATA,FORMAT,TAG>:: 
GeneralField():format(FORMAT(CommonPrms::instance()->Nvol())),
	       data(format.size()){}
/*
// Specialization for even/odd field
template <>
inline GeneralField<Field,Format::Format_F,HalfSiteTag>:: 
GeneralField():format(CommonPrms::instance()->Nvol()/2),
	       data(format.size()){}
template <>
inline GeneralField<Field,Format::Format_G,HalfSiteTag>:: 
GeneralField():format(CommonPrms::instance()->Nvol()/2),
	       data(format.size()){}
*/
// Specialization for one dimension gauge field
template <> 
inline GeneralField<Field,Format::Format_G,OneDimTag>::
GeneralField():format(Format::Format_G(CommonPrms::instance()->Nvol(),1)),
	       data(format.size()){}

// Constructor with 1 argument //////////////////////
template <class DATA,class FORMAT,typename TAG> 
GeneralField<DATA,FORMAT,TAG>::GeneralField(int LocalVol)
  :format(FORMAT(LocalVol)),
   data(format.size()){}

// Specialization for one dimension gauge field
template <> 
inline GeneralField<Field,Format::Format_G,OneDimTag>::
GeneralField(int LocalVol):format(Format::Format_G(LocalVol,1)),
			   data(format.size()){}

// Copy constructor ////////////////////////////////
template <class DATA,class FORMAT,typename TAG> 
GeneralField<DATA,FORMAT,TAG>::GeneralField(const DATA& Fin)
  :format(CommonPrms::instance()->Nvol()), // cannot accommodate e/o-field
   data(Fin){  assert(Fin.size()== format.size());}
/*
// Specialization for even/odd field
template <>
inline GeneralField<Field,Format::Format_F,HalfSiteTag>::
GeneralField(const Field& Fin):format(CommonPrms::instance()->Nvol()/2),
			       data(Fin){ assert(Fin.size()== format.size());}
template <>
inline GeneralField<Field,Format::Format_G,HalfSiteTag>::
GeneralField(const Field& Fin):format(CommonPrms::instance()->Nvol()/2),
			       data(Fin){ assert(Fin.size()== format.size());}
*/

// Specialization for one dimension gauge field
/*! Since Nex is fixed to be 1 in this case, Nvol can be determined 
  as Fin.size()/Format::Format_G::Nin(), then e/o-field is also possible */
template <> 
inline GeneralField<Field,Format::Format_G,OneDimTag>::
GeneralField(const Field& Fin)
  :format(Fin.size()/Format::Format_G::Nin(),1),
   data(Fin){  assert(Fin.size()== format.size());}

// Assignment operators ////////////////////////////
template <class DATA,class FORMAT,typename TAG> 
inline GeneralField<DATA,FORMAT,TAG>& 
GeneralField<DATA,FORMAT,TAG>::operator=(const GeneralField& rhs){
  assert(format.size() == rhs.format.size());
  data = rhs.data;
  return *this;
}

template <class DATA,class FORMAT,typename TAG> 
inline GeneralField<DATA,FORMAT,TAG>& 
GeneralField<DATA,FORMAT,TAG>::operator=(const Field& rhs){
  assert(format.size() == rhs.size());
  data = rhs;
  return *this;
}

template <class DATA,class FORMAT,typename TAG> 
inline GeneralField<DATA,FORMAT,TAG>& 
GeneralField<DATA,FORMAT,TAG>::operator=(const std::valarray<double>& rhs){
  assert(format.size() == rhs.size());
  data = rhs;
  return *this;
}

template <class DATA,class FORMAT,typename TAG> 
inline GeneralField<DATA,FORMAT,TAG>&
GeneralField<DATA,FORMAT,TAG>::operator=(double rhs){
  data = rhs;
  return *this;
}

template <class DATA,class FORMAT,typename TAG>
inline GeneralField<DATA,FORMAT,TAG>&
GeneralField<DATA,FORMAT,TAG>::operator+=(const GeneralField& rhs){
  assert(format.size() == rhs.format.size());
  data += rhs.data;
  return *this;
}

template <class DATA,class FORMAT,typename TAG>
inline GeneralField<DATA,FORMAT,TAG>&
GeneralField<DATA,FORMAT,TAG>::operator+=(const Field& rhs){
  assert(format.size() == rhs.size());
  data += rhs;
  return *this;
}

template <class DATA,class FORMAT,typename TAG>
inline GeneralField<DATA,FORMAT,TAG>&
GeneralField<DATA,FORMAT,TAG>::operator+=(const std::valarray<double>& rhs){
  assert(format.size() == rhs.size());
  data += rhs;
  return *this;
}

template <class DATA,class FORMAT,typename TAG> 
inline GeneralField<DATA,FORMAT,TAG>&
GeneralField<DATA,FORMAT,TAG>::operator-=(const GeneralField& rhs){
  assert(format.size() == rhs.format.size());
  data -= rhs.data;
  return *this;
}

template <class DATA,class FORMAT,typename TAG> 
inline GeneralField<DATA,FORMAT,TAG>&
GeneralField<DATA,FORMAT,TAG>::operator-=(const Field& rhs){
  assert(format.size() == rhs.size());
  data -= rhs;
  return *this;
}

template <class DATA,class FORMAT,typename TAG> 
inline GeneralField<DATA,FORMAT,TAG>&
GeneralField<DATA,FORMAT,TAG>::operator-=(const std::valarray<double>& rhs){
  assert(format.size() == rhs.size());
  data -= rhs;
  return *this;
}

template <class DATA,class FORMAT,typename TAG> 
inline GeneralField<DATA,FORMAT,TAG>&
GeneralField<DATA,FORMAT,TAG>::operator*=(double rhs){
  data *= rhs;
  return *this;
}

template <class DATA,class FORMAT,typename TAG> 
inline GeneralField<DATA,FORMAT,TAG>&
GeneralField<DATA,FORMAT,TAG>::operator/=(double rhs){
  data /= rhs;
  return *this;
}

template <class DATA,class FORMAT,typename TAG> 
double GeneralField<DATA,FORMAT,TAG>::norm(){ return data.norm();}

//////////////////////////////////////////////////////////////////////////////
// Notes
// GaugeField shoud contain also methods to create from geometry and initialize
typedef GeneralField<Field,Format::Format_A>              AdjGaugeField;
typedef GeneralField<Field,Format::Format_G>              GaugeField;
typedef GeneralField<Field,Format::Format_G,OneDimTag>    GaugeField1D;
typedef GeneralField<Field,Format::Format_F>              FermionField;
//typedef GeneralField<Field,Format::Format_F,HalfSiteTag>  FermionEOField;
typedef GeneralField<std::vector<Field>,Format::Format_F> PropagatorField;

class GaugeGlobal: public GaugeField{
public:
  GaugeGlobal(Geometry geom){
    format = Format::Format_G( geom.prms_->Nvol());
    data = Field(format.size());
  }

  GaugeField& operator=(GaugeGlobal& Source){
    assert (format.size() == Source.format.size());
    data = Source.data;
    return *this;
  }

  /*! Initializes the gauge field with random number */
  void initializeRand(const RandNum& rand){
    GaugeConf_rand gconf(format,rand);
    gconf.init_conf(data);
  }

  /*! Initializes the gauge field with an \n
   * external configuration in file <Filename>
   *
   * Configuration type in XML: TextFile
   * @param Filename String containing the filename
   */
  void initializeTxt(const std::string &Filename){
    GaugeConf_txt gconf(format,Filename);
    gconf.init_conf(data);
  }

  /*! Initializes the gauge field with an \n
   * external configuration in binary format contained in file <Filename>
   *
   * Configuration type in XML: Binary
   * @param Filename String containing the filename
   */
  void initializeBin(const std::string &Filename){
    GaugeConf_bin gconf(format,Filename);
    gconf.init_conf(data);
  }
  void initializeJLQCDlegacy(const std::string &Filename) {
    GaugeConf_JLQCDLegacy gconf(format,Filename);
    gconf.init_conf(data);
  }

 /*! @brief Initializes the gauge field with \n unit matrices
   *
   * Configuration type in XML: Unit
   */
  void initializeUnit(){
    GaugeConf_unit gconf(format);
    gconf.init_conf(data);
  }

  int initialize(XML::node node) {
    try {
      XML::descend(node, "Configuration");
      if (!XML::attribute_compare(node,"Type","TextFile")){
	std::string filename(node.child_value());
	initializeTxt(filename);
	return 0;
      }
      if (!XML::attribute_compare(node,"Type","Unit")){
	initializeUnit();
	return 0;
      }
      if (!XML::attribute_compare(node,"Type","Binary")){
	std::string filename(node.child_value());
	initializeBin(filename);
	return 0;
      }
      if (!XML::attribute_compare(node,"Type","JLQCDlegacy")){
	std::string filename(node.child_value());
	initializeJLQCDlegacy(filename);
	return 0;
      }
      std::cout << "Configuration type unknown\n";
      exit(1);
    } catch(...) {
      std::cout << "Error in initialization of gauge field "<< std::endl;
      abort();

    }
    return 0;
  }

  int initialize(XML::node node,std::string filename){
    try {
      XML::descend(node,"Configuration");
      if (!XML::attribute_compare(node,"Type","TextFile")){
	initializeTxt(filename);
	return 0;
      }
      if (!XML::attribute_compare(node,"Type","Binary")){
	initializeBin(filename);
	return 0;
      }
      if (!XML::attribute_compare(node,"Type","JLQCDlegacy")){
	initializeJLQCDlegacy(filename);
	return 0;
      }
      std::cout << "Configuration type unknown\n";
      exit(1);
    } catch(...) {
      std::cout << "Error in initialization of gauge field "<< std::endl;
    }
    return 0;
  }

private:
  GaugeGlobal(const GaugeGlobal&);//hide copy constructor
};

#endif //COMMON_FIELDS_H_
