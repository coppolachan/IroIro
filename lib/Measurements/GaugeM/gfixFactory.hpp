/*! @file gfixFactory.hpp
    @brief Declaration of the GFixFactory class 
*/
#ifndef GFIX_FACTORY_INCLUDED
#define GFIX_FACTORY_INCLUDED

#include "gaugeFixing_Free.hpp"
#include "gaugeFixing_Landau.hpp"
#include "gaugeFixing_Coulomb.hpp"
#include "include/pugi_interface.h"

class RandNum;

class GFixFactory{
public:
  virtual ~GFixFactory(){}
  virtual GaugeFixing* getGaugeFixing(const RandNum&) = 0;
};

class GFixFactory_Free: public GFixFactory {
  GaugeFixing* getGaugeFixing(const RandNum&){ 
    return new GaugeFixing_Free;  
  }
};

class GFixFactory_Landau: public GFixFactory {
  const XML::node gfix_node_;
public:
  GFixFactory_Landau(XML::node node):gfix_node_(node){}
  GaugeFixing* getGaugeFixing(const RandNum& rng){
    return new GaugeFixing_Landau(rng,gfix_node_);
  }
};

class GFixFactory_Coulomb: public GFixFactory {
  const XML::node gfix_node_;
public:
  GFixFactory_Coulomb(XML::node node):gfix_node_(node){}
  GaugeFixing* getGaugeFixing(const RandNum& rng){ 
    return new GaugeFixing_Coulomb(rng,gfix_node_);
  }
};

namespace GaugeFix{
  GFixFactory* createGaugeFixingFactory(XML::node&);
}
#endif
