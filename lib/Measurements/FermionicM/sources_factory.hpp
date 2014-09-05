/*!
 * @file sources_factory.hpp 
 * @brief Propagator Source factories
 */
#ifndef SOURCES_FACT_
#define SOURCES_FACT_

#include <string>
#include "include/pugi_interface.h"
#include "include/errors.hpp"

#include "Measurements/FermionicM/source_types.hpp"
#include "Tools/RAIIFactory.hpp"
#include "Tools/randNum_Factory.hpp"

/*!
 * @brief Abstract base class for creating Source
 */
class SourceFactory {
public:
  virtual Source* getSource(int Local_Dim = CommonPrms::instance()->Nvol()) = 0;
  virtual ~SourceFactory(){}
};
/////////////////////////////////////////////////
/*!
 @brief Concrete class for creating Local Source operator
*/
template <typename Format>
class LocalSourceFactory: public SourceFactory {
  const XML::node srcNode_;

public:
  LocalSourceFactory(XML::node node):srcNode_(node){}

  Source* getSource(int Local_Dim = CommonPrms::instance()->Nvol()){
    std::vector<int> position;
    XML::read_array(srcNode_,"position",position,MANDATORY);
    return new Source_local<Format>(position,Local_Dim);
  }
};

template <typename Format>
class ExpSourceFactory: public SourceFactory {
  const XML::node srcNode_;

public:
  ExpSourceFactory(XML::node node):srcNode_(node){}

  Source* getSource(int Local_Dim = CommonPrms::instance()->Nvol()){
    std::vector<int> position;
    double alpha;
    XML::read_array(srcNode_,"position", position, MANDATORY);
    XML::read(srcNode_,"alpha", alpha, MANDATORY);
    return new Source_exp<Format>(position,alpha,Local_Dim);
  }
};

template <typename Format>
class GaussSourceFactory: public SourceFactory {
  const XML::node srcNode_;

public:
  GaussSourceFactory(XML::node node):srcNode_(node){}

  Source* getSource(int Local_Dim = CommonPrms::instance()->Nvol()){
    std::vector<int> position;
    double alpha;
    XML::read_array(srcNode_,"position", position, MANDATORY);
    XML::read(srcNode_,"alpha", alpha, MANDATORY);
    return new Source_Gauss<Format>(position,alpha,Local_Dim);
  }
};

template <typename Format>
class WallSourceFactory: public SourceFactory {
  const XML::node srcNode_;
public:
  WallSourceFactory(XML::node node):srcNode_(node){}

  Source* getSource(int Local_Dim = CommonPrms::instance()->Nvol()){
    int position;
    XML::read(srcNode_, "position", position, MANDATORY);
    return new Source_wall<Format>(position,Local_Dim);
  }
};

template <typename Format>
class WhiteNoiseSourceFactory: public SourceFactory {
  const XML::node srcNode_;
  RaiiFactoryObj<RandNum> rng_;
public:
  WhiteNoiseSourceFactory(XML::node node):srcNode_(node){}

  Source* getSource(int Local_Dim = CommonPrms::instance()->Nvol()){
    rng_.save(RNG_Env::RandNumG::instance().getRNG());
    if (rng_.get() == NULL) {
      Errors::XMLerr("Random number generator not defined\nPlease provide a correct definition in the XML\n");
    }  
    return new Source_wnoise<Format>(*rng_.get(),Local_Dim);
  }
};

template <typename Format>
class Z2noiseSourceFactory: public SourceFactory {
  const XML::node srcNode_;
  RaiiFactoryObj<RandNum> rng_;
public:
  Z2noiseSourceFactory(XML::node node):srcNode_(node){}

  Source* getSource(int Local_Dim = CommonPrms::instance()->Nvol()){
    Z2Type NoiseType;
    std::string Z2Type_name;
    XML::read(srcNode_,"NoiseType", Z2Type_name, MANDATORY);
    if (!EnumString<Z2Type>::To( NoiseType, Z2Type_name )){
      CCIO::cerr<<"Error: string ["<< Z2Type_name <<"] not valid"<<std::endl;
      CCIO::cerr<<"Request from node ["<< srcNode_.name()<<"]\n";
      abort();
    } else {
      CCIO::cout<<"Choosing Z2Type type: "<< Z2Type_name << std::endl;
    }
    rng_.save(RNG_Env::RandNumG::instance().getRNG());
    if (rng_.get() == NULL) {
      Errors::XMLerr("Random number generator not defined\nPlease provide a correct definition in the XML\n");
    }
    return new Source_Z2noise<Format>(*rng_.get(),Local_Dim,NoiseType);
  }
};

///////////////////////////////////////////////////
/*!
 * @brief Namespace for Sources factory caller
 */
namespace Sources{
  template <typename Format>
  SourceFactory* createSourceFactory(const XML::node node){
    
    if (node !=NULL) {
      const char* Source_name = node.attribute("type").value();
      
      if (!strcmp(Source_name, "Local")) { 
        return new LocalSourceFactory<Format>(node);
      }
      if (!strcmp(Source_name, "Exp")) { 
        return new ExpSourceFactory<Format>(node);
      }
      if (!strcmp(Source_name, "Gauss")) { 
        return new GaussSourceFactory<Format>(node);
      }
      if (!strcmp(Source_name, "Z2noise")) { 
        return new Z2noiseSourceFactory<Format>(node);
      }
      if (!strcmp(Source_name, "Wall")) { 
        return new WallSourceFactory<Format>(node);
      }
      if (!strcmp(Source_name, "WhiteNoise")) { 
        return new WhiteNoiseSourceFactory<Format>(node);
      }
      ErrorString msg;
      msg << "No Source available with name [" << Source_name << "]\n";
      Errors::XMLerr(msg);
    } else {
      std::cout << "Requested node is missing in input file "
		<< "(Source Object)\n" 
		<< "Request by " << node.parent().name() << std::endl;
      abort();
    }
  }  
};

#endif 
