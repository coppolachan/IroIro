/*!
 * @file common_fields.hpp
 *
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
#include "lib/Tools/sunMat.h"
#include "include/macros.hpp"

typedef Format::Format_G GaugeFieldFormat;/**< Format of gauge field
					     at compilation time */


/*!
 * @brief A Class to handle gauge fields
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a> 
 * 
 * This class handles the Format 
 * and Field class in order to bind them 
 * into a single self contained object
 * that knows the details of storage
 *
 */
class GaugeField {
public:
  GaugeFieldFormat Format;/**< Format specifier */
  Field U; /**< Configuration container */
  
  /*! 
   * Default constructor\n 
   * Takes Geometry class as input
   * @param geom Defines the lattice geometry
   */
  GaugeField(Geometry geom): 
    Format(GaugeFieldFormat( geom.parameters->Nvol() )),
    U(Field(Format.size())){}

  /*! 
   * Default constructor\n 
   * No parameters, automatically creates a 4D field 
   * of local dimension Nvol
   */
  GaugeField(): 
    Format(GaugeFieldFormat( CommonPrms::instance()->Nvol() )),
    U(Field(Format.size())){}

  /*! 
   * Constructor\n 
   * to contain given Field data
   */
  explicit GaugeField(Field& Uin): 
    Format(GaugeFieldFormat( CommonPrms::instance()->Nvol())),
    U(Uin){ 
    assert(U.size()==Format.size());
  }

  /*!
   * Initializes the gauge field with an \n
   * external configuration in file <Filename>
   *
   * Configuration type in XML: TextFile
   * @param Filename String containing the filename
   */
  void initializeTxt(const std::string &Filename) {
    GaugeConf_txt gconf(Format,Filename);
    gconf.init_conf(U);
  }

  /*!
   * Initializes the gauge field with an \n
   * external configuration in binary format
   * contained in file <Filename>
   *
   * Configuration type in XML: Binary
   * @param Filename String containing the filename
   */
  void initializeBin(const std::string &Filename) {
    GaugeConf_bin gconf(Format,Filename);
    gconf.init_conf(U);
  }

  void initializeJLQCDlegacy(const std::string &Filename) {
    GaugeConf_JLQCDLegacy gconf(Format,Filename);
    gconf.init_conf(U);
  }

 /*!
   * Initializes the gauge field with \n
   * unit matrices
   *
   * Configuration type in XML: Unit
   */
  void initializeUnit(){
    GaugeConf_unit gconf(Format);
    gconf.init_conf(U);
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
    } catch(...) {
      std::cout << "Error in initialization of gauge field "<< std::endl;
    }

    return 0;
  }

  inline GaugeField& operator=(const GaugeField& rhs){
    U = rhs.U;
    return *this;
  }

  inline SUNmat matrix(int site, int ex) const {
    return SUNmat(U[Format.cslice(0,site,ex)]);    
  }
};

class GaugeField1D {

public:
  const GaugeFieldFormat Format;/**< Format specifier */
  Field U; /**< Field container */
  
  /*! 
   * Default constructor\n
   *
   * Contains links only in one direction
   */
  GaugeField1D(): 
    Format(GaugeFieldFormat(CommonPrms::instance()->Nvol(),1)),
    U(Field(Format.size())){}

  GaugeField1D(Field& f): 
    Format(GaugeFieldFormat(CommonPrms::instance()->Nvol(),1)),
    U(f){}

  GaugeField1D(std::valarray<double> f): 
    Format(GaugeFieldFormat(CommonPrms::instance()->Nvol(),1)),
    U(Field(f)){}

  void set_matrix(int site, SUNmat matrix) {
    U.set(Format.cslice(0,site),matrix.getva());
  }

  inline GaugeField1D& operator=(const GaugeField1D& rhs){
    U = rhs.U;
    return *this;

  }

  inline SUNmat matrix(int site, int ex){
    return SUNmat(U[Format.cslice(0,site,ex)]);    
  }
};

namespace FieldUtils{
  const Field TracelessAntihermite(const GaugeField&);
}










//////////////////////////////////////////////////////////////////////
struct ExtraDimTag {};
struct OneDimTag {};


template < class DATA, class FORMAT, typename TAG = NullType> 
class GeneralField {
  FORMAT format;
  DATA data;
  GeneralField(const GeneralField& rhs); //hide copy constructor
public:
  GeneralField();
  GeneralField(int);

  /*! 
   * Constructor\n 
   * to store given data
   */
  explicit GeneralField(DATA&);


  GeneralField& operator=(const GeneralField& rhs);
  GeneralField& operator+=(const GeneralField& rhs);
  GeneralField& operator-=(const GeneralField& rhs);



};

template < class DATA, class FORMAT, typename TAG> 
GeneralField<DATA,FORMAT,TAG>:: 
GeneralField():format(FORMAT( CommonPrms::instance()->Nvol())),
	       data(format.size()){}

// Specialization for one dimension gauge field
template <> 
inline GeneralField< Field, Format::Format_G, OneDimTag>:: 
GeneralField():format(Format::Format_G( CommonPrms::instance()->Nvol(),1)),
	       data(format.size()){}

template < class DATA, class FORMAT, typename TAG> 
GeneralField<DATA,FORMAT,TAG>:: 
GeneralField(int LocalVol):format(FORMAT( LocalVol)),
			   data(format.size()){}

// Specialization for one dimension gauge field
template <> 
inline GeneralField< Field, Format::Format_G, OneDimTag>:: 
GeneralField(int LocalVol):format(Format::Format_G( LocalVol,1)),
			   data(format.size()){}


template < class DATA, class FORMAT, typename TAG> 
GeneralField<DATA,FORMAT,TAG>:: 
GeneralField(DATA& FieldIn):format(FORMAT( CommonPrms::instance()->Nvol()))
{
  assert(FieldIn.size()== format.size());
  data = FieldIn;
}

// Specialization for one dimension gauge field
template <> 
inline GeneralField< Field, Format::Format_G, OneDimTag>::
GeneralField(Field& FieldIn):format(Format::Format_G( CommonPrms::instance()->Nvol(),1))
{
  assert(FieldIn.size()== format.size());
  data = FieldIn;
}

template < class DATA, class FORMAT, typename TAG> 
inline GeneralField<DATA,FORMAT,TAG>& 
GeneralField<DATA,FORMAT,TAG>::operator=(const GeneralField& rhs)
{
  data = rhs.data;
  return *this;
}

template < class DATA, class FORMAT, typename TAG>
inline GeneralField<DATA,FORMAT,TAG>&
GeneralField<DATA,FORMAT,TAG>::operator+=(const GeneralField& rhs) 
{
  data += rhs.data;
  return *this;
}

template < class DATA, class FORMAT, typename TAG> 
inline GeneralField<DATA,FORMAT,TAG>&
GeneralField<DATA,FORMAT,TAG>::operator-=(const GeneralField& rhs) 
{
  data -= rhs.data;
  return *this;
}
////////////////////////////////////////////////////////////////
// Specialization for extradimension

template < class DATA, class FORMAT> 
class GeneralField<DATA, FORMAT, ExtraDimTag> {
  FORMAT format;
  DATA data;
  GeneralField(){}; //hide default constructor
  GeneralField(const GeneralField& rhs); //hide copy constructor
public:
  GeneralField(int Ls):format(FORMAT( CommonPrms::instance()->Nvol()*Ls)),
		       data(format.size()){}
  
  GeneralField(int Ls, int LocalVol):format(FORMAT( LocalVol*Ls)),
				     data(format.size()){}

  /*! 
   * Constructor\n 
   * to store given data
   */
  explicit GeneralField(DATA&);
  GeneralField& operator=(const GeneralField& rhs);
  GeneralField& operator+=(const GeneralField& rhs);
  GeneralField& operator-=(const GeneralField& rhs);

};

//////////////////////////////////////////////////////////////////////////////
typedef GeneralField< Field, Format::Format_A >              AdjGaugeField;
typedef GeneralField< Field, Format::Format_G >              GaugeFieldType;
typedef GeneralField< Field, Format::Format_G, OneDimTag >   GaugeField1DType;
typedef GeneralField< Field, Format::Format_F >              FermionField;
typedef GeneralField< Field, Format::Format_F, ExtraDimTag > FermionFieldExtraDim;
typedef GeneralField< std::vector<Field>, Format::Format_F > PropagatorField;


#endif //COMMON_FIELDS_H_
