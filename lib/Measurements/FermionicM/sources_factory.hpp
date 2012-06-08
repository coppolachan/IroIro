/*!
 * @file sources_factory.hpp 
 * @brief Propagator Source factories
 */
#ifndef SOURCES_FACT_
#define SOURCES_FACT_

#include <string>
#include "include/pugi_interface.h"
#include "Measurements/FermionicM/source_types.hpp"
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
  const XML::node Source_node;

public:
  LocalSourceFactory(XML::node node):Source_node(node){}

  Source* getSource(int Local_Dim = CommonPrms::instance()->Nvol()){
    std::vector<int> position;
    XML::read_array(Source_node, "position", position, MANDATORY);
    return new Source_local<Format>(position,
				    Local_Dim);
  }

};

template <typename Format>
class ExpSourceFactory: public SourceFactory {
  const XML::node Source_node;

public:
  ExpSourceFactory(XML::node node):Source_node(node){}

  Source* getSource(int Local_Dim = CommonPrms::instance()->Nvol()){
    std::vector<int> position;
    double alpha;
    XML::read_array(Source_node, "position", position, MANDATORY);
    XML::read(Source_node, "alpha", alpha, MANDATORY);
    return new Source_exp<Format>(position,
				  alpha,
				  Local_Dim);
  }
};

template <typename Format>
class GaussSourceFactory: public SourceFactory {
  const XML::node Source_node;

public:
  GaussSourceFactory(XML::node node):Source_node(node){}

  Source* getSource(int Local_Dim = CommonPrms::instance()->Nvol()){
    std::vector<int> position;
    double alpha;
    XML::read_array(Source_node, "position", position, MANDATORY);
    XML::read(Source_node, "alpha", alpha, MANDATORY);
    return new Source_Gauss<Format>(position,
				    alpha,
				    Local_Dim);
  }
};

template <typename Format>
class WallSourceFactory: public SourceFactory {
  const XML::node Source_node;

public:
  WallSourceFactory(XML::node node):Source_node(node){}

  Source* getSource(int Local_Dim = CommonPrms::instance()->Nvol()){
    int position;
    XML::read(Source_node, "position", position, MANDATORY);
    return new Source_wall<Format>(position,
				   Local_Dim);
  }
};

template <typename Index,typename Format>
class WhiteNoiseSourceFactory: public SourceFactory {
  const XML::node Source_node;

public:
  WhiteNoiseSourceFactory(XML::node node):Source_node(node){}

  Source* getSource(int Local_Dim = CommonPrms::instance()->Nvol()){
    return new Source_wnoise<Index,Format>(*RNG_Env::RNG->getRandomNumGenerator(),
					   Local_Dim);
  }
};

template <typename Index,typename Format>
class Z2noiseSourceFactory: public SourceFactory {
  const XML::node Source_node;

public:
  Z2noiseSourceFactory(XML::node node):Source_node(node){}

  Source* getSource(int Local_Dim = CommonPrms::instance()->Nvol()){
    Z2Type NoiseType;
    std::string Z2Type_name;
    XML::read(Source_node, "NoiseType", Z2Type_name, MANDATORY);
    if (!EnumString<Z2Type>::To( NoiseType, Z2Type_name )){
      CCIO::cerr<<"Error: string ["<< Z2Type_name <<"] not valid"<<std::endl;
      CCIO::cerr<<"Request from node ["<< Source_node.name()<<"]\n";
      abort();
    } else {
      CCIO::cout<<"Choosing Z2Type type: "<< Z2Type_name << std::endl;
    }
    return new 
      Source_Z2noise<Index,Format>(*RNG_Env::RNG->getRandomNumGenerator(),
				   Local_Dim,
				   NoiseType);
  }
};

///////////////////////////////////////////////////
/*!
 * @brief Namespace for Sources factory caller
 */
namespace Sources{
  template <typename Index,typename Format>
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
        return new Z2noiseSourceFactory<Index,Format>(node);
      }
      if (!strcmp(Source_name, "Wall")) { 
        return new WallSourceFactory<Format>(node);
      }
      if (!strcmp(Source_name, "WhiteNoise")) { 
        return new WhiteNoiseSourceFactory<Index,Format>(node);
      }
      if (!strcmp(Source_name, "Z2noise")) { 
        return new Z2noiseSourceFactory<Index,Format>(node);
      }
      std::cerr << "No Source available with name ["
		<< Source_name << "]" << std::endl;
      abort();
    } else {
      std::cout << "Requested node is missing in input file "
		<< "(Source Object)\n" 
		<< "Request by " << node.parent().name() << std::endl;
      abort();
    }
  }  
};

#endif 
