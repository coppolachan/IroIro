/*!
 * @file Smearing_Factory.hpp 
 *
 * @brief Declaration of Smearing Operators factories
 */

#ifndef SMEARING_FACT_H_
#define SMEARING_FACT_H_

/*!
 * @brief Abstract base class for creating smearing operators
 */
class SmearingOperatorFactory {
public:
  virtual Smear* getSmearingOperator() = 0;
}


/////////////////////////////////////////////////////////
/*!
 * @brief Concrete class for creating APE smearing operators
 */
  class APESmearingFactory: public SmearingOperatorFactory {
    const XML::node APE_node;

  public:
    APESmearingFactory(const XML::node node):APE_node(node){}

    Smear* getSmearingOperator(){
      return new Smear_APE(node);
    }


  }





#endif
