/*!
 * @file test_MesonSpectrum.cpp
 * @brief Definition of classes for calculating meson correlators
 */
#include "test_MesonSpectrum.hpp"
#include "include/factories.hpp"
#include "Measurements/FermionicM/quark_prop_meas_factory.hpp"
#include "Measurements/FermionicM/qprop_mom.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include "Measurements/GaugeM/staples.hpp"

using namespace std;

int Test_MesonSpectrum::run(){
  Staples Staple;
  CCIO::cout<< "Plaquette (thin): "<< Staple.plaquette(conf_) <<"\n";


  //// smearing ////
  XML::node smr_node = node_; 
  XML::descend(smr_node,"Smearing"); 

  int Nsmear;                                    
  XML::read(smr_node,"Nsmear",Nsmear,MANDATORY);  

  SmearingOperatorFactory* SmrFactory = 
    SmearingOperators::createSmearingOperatorFactory(smr_node);
  Smear* SmearingObj = SmrFactory->getSmearingOperator();

  smeared_u_= conf_; // Copy original configuration to smeared_u_ 
  for(int i=0; i<Nsmear; ++i){ // Do the actual smearing 
    GaugeField previous_u_= smeared_u_;
    SmearingObj->smear(smeared_u_,previous_u_);
  }
  CCIO::cout<< "Plaquette (smeared): "<< Staple.plaquette(smeared_u_)
    	    << std::endl;

  //// random number generator ////
  XML::node rng_node = node_;  // copy of the root node
  RNG_Env::RNG = RNG_Env::createRNGfactory(rng_node);

  //// gauge fixing ////
  XML::node gfix_node = node_;
  GFixFactory* gffctry = GaugeFix::createGaugeFixingFactory(gfix_node);
  GaugeFixing* gfix 
    = gffctry->getGaugeFixing(*(RNG_Env::RNG->getRandomNumGenerator()));

  fixed_u_= gfix->do_fix(smeared_u_);

  CCIO::cout<< "Plaquette (gauge fixed): "<< Staple.plaquette(fixed_u_)
	    << std::endl;
  
  //// Quark Propagator ////
  XML::descend(node_,"QuarkProp");
  QuarkPropagatorFactory* 
    qpfact = QuarkPropagators::createQuarkPropagatorFactory(node_);
  QuarkPropagator* qprop = qpfact->getQuarkProp(fixed_u_);

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

  //// Meson correlators ////
  //GammaMatrices::Unit Gamma; //pion
  //MesonCorrelator meson(Gamma,Gamma);
  
  
  {// pion correlation function 
    MesonCorrelator meson(Pion);
    
    std::vector<double> mcorr = meson.calculate<Format::Format_F>(sq,sq);  
    CCIO::cout<<"---pp meson correlator---"<<endl;
    for(int t=0; t<mcorr.size(); ++t){
      CCIO::cout<< setiosflags(  ios_base::scientific);
      CCIO::cout<< setw(3) <<setiosflags(ios_base::right)<< t;
      CCIO::cout<< setw(25)<<setiosflags(ios_base::left )<< mcorr[t]<<endl;
      CCIO::cout<< resetiosflags(ios_base::scientific);
    }
  }
  {// vector1 correlation function 
    MesonCorrelator meson(Vector1);

    std::vector<double> mcorr = meson.calculate<Format::Format_F>(sq,sq);  
    CCIO::cout<<"---v1v1 meson correlator---"<<endl;
    for(int t=0; t<mcorr.size(); ++t){
      CCIO::cout<< setiosflags(  ios_base::scientific);
      CCIO::cout<< setw(3) <<setiosflags(ios_base::right)<< t;
      CCIO::cout<< setw(25)<<setiosflags(ios_base::left )<< mcorr[t]<<endl;
      CCIO::cout<< resetiosflags(ios_base::scientific);
    }
  }
  return 0;
}
