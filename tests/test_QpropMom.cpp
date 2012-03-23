/*!
 * @file test_QpropMom.cpp
 * @brief Definition of classes for calculating Sq(p)
 */
#include "test_QpropMom.hpp"
#include "include/factories.hpp"
#include "Measurements/FermionicM/quark_prop_meas_factory.hpp"
#include "Measurements/FermionicM/qprop_mom.hpp"
#include "Measurements/GaugeM/staples.hpp"

int Test_QpropMom::run(){
  //  Staples Staple;
  //  CCIO::cout<< "Plaquette (thin): "<< Staple.plaquette(conf_) <<"\n";

  //// smearing ////
  XML::node smr_node = node_;  // copy of the root node
  XML::descend(smr_node,"Smearing"); 

  int Nsmear;                                    
  XML::read(smr_node,"Nsmear",Nsmear,MANDATORY);  

  CCIO::cout<<"SmrFactory initialized"<<std::endl;
  SmearingOperatorFactory* SmrFactory = 
    SmearingOperators::createSmearingOperatorFactory(smr_node);
  
  Smear* SmearingObj = SmrFactory->getSmearingOperator();
  CCIO::cout<<"SmearingObj obtained"<<std::endl;

  smeared_u_= conf_; // Copy original configuration to smeared_u_ 
  for(int i=0; i<Nsmear; i++){ // Do the actual smearing 
    GaugeField previous_u_= smeared_u_;
    SmearingObj->smear(smeared_u_,previous_u_);
  }
  //CCIO::cout<< "Plaquette (smeared): "<< Staple.plaquette(smeared_u_)
  //	    << std::endl;
  //// Quark Propagator ////
  XML::descend(node_,"QuarkProp");
  QuarkPropagatorFactory* 
    qpfact = QuarkPropagators::createQuarkPropagatorFactory(node_);
  QuarkPropagator* qprop = qpfact->getQuarkProp(smeared_u_);

  //// source creation ////
  XML::next_sibling(node_,"Source");
  SourceFactory* SrcFactory 
    = Sources::createSourceFactory<SiteIndex,Format::Format_F>(node_);
  Source* src = SrcFactory->getSource();

  prop_t sq;  //Defines a vector of fields
  CCIO::cout << " ---- Calculating propagator\n";
  qprop->calc(sq,*src);

  //  CCIO::ReadFromDisk < Format::Format_F >(sq, "propagator.bin", 12);
  //  CCIO::SaveOnDisk < Format::Format_F >(sq, "propagator.bin");

  //// Fourier transformation & save ////
  QpropMom qprop_mom(sq);
  CCIO::cout<<"output starts"<<std::endl;
  qprop_mom.output();

  return 0;
}
