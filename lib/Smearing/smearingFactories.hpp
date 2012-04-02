/*!
 * @file smearing_Factories.hpp 
 * @brief Declaration of Smearing Operators factories
 */
#ifndef SMEARING_FACT_H_
#define SMEARING_FACT_H_

#include "Smearing/APEsmear.hpp"
#include "Smearing/stoutSmear.hpp"

#include "Tools/RAIIFactory.hpp"
#include "include/pugi_interface.h"
#include "Communicator/comm_io.hpp"
/*!
 * @brief Abstract base class for creating smearing operators
 */
class SmearingOperatorFactory {
public:
  virtual Smear* getSmearingOperator() = 0;
};

namespace SmearingOperators {
  SmearingOperatorFactory* createSmearingOperatorFactory(const XML::node);
}

/*! @brief Concrete class for creating APE smearing operators */
class APESmearingFactory: public SmearingOperatorFactory {
  std::vector<double> rho;
public:
  APESmearingFactory(const XML::node node){
    int dims = CommonPrms::instance()->Ndim();
    XML::read_array(node, "rho", rho, MANDATORY);
    if (rho.size() == 1){
      double rho_val = rho[0];
      rho.resize(dims*dims);
      for(int i = 0; i < dims*dims; ++i) rho[i] = rho_val;
      for(int mu = 0; mu < dims; ++mu) 	rho[mu + mu*dims] = 0.0;
    }
    if((rho.size() !=1) && (rho.size() != dims*dims)){
      CCIO::cout <<"[APESmearingFactory] Error in rho size\n";
      abort();
    }
  }
  
  Smear_APE* getSmearingOperator(){ return new Smear_APE(rho);}
};

/*! @brief Concrete class for creating Stout smearing operators */
class StoutSmearingFactory: public SmearingOperatorFactory {
  RaiiFactoryObj<SmearingOperatorFactory> BaseSmearingObj;

public:
  StoutSmearingFactory(XML::node node){
    XML::descend(node, "Base", MANDATORY);
    BaseSmearingObj.save(SmearingOperators::createSmearingOperatorFactory(node));
  }
  
  Smear_Stout* getSmearingOperator(){
    return new Smear_Stout(BaseSmearingObj.get()->getSmearingOperator()); }
};

#endif
