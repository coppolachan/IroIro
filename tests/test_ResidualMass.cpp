/*!
 * @file test_ResidualMass.cpp
 * @brief Definition of classes for calculating the Residual Mass
 */
#include "test_ResidualMass.hpp"
#include "include/factories.hpp"
#include "Measurements/FermionicM/qprop_DomainWall.hpp"
#include "Measurements/FermionicM/utils_DWF4d.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include <ctime>

int Test_ResMass::run(XML::node node) {
//create the factory for Rand Numbers Gen
  RNG_Env::RNG = RNG_Env::createRNGfactory(node); 

  Staples Staple(conf_.Format);
  CCIO::cout<< "Plaquette : "<< Staple.plaquette(conf_.U)<< std::endl;
  
  XML::descend(node, "QuarkDWFProp");
  QPropDWFFactory QP_DomainWallFact(node);//specific factory (as a test program)
  QpropDWF* QuarkPropDW = QP_DomainWallFact.getQuarkPropDW(conf_);

  XML::next_sibling(node, "Source");
  SourceFactory* Source_Factory 
    = Sources::createSourceFactory<SiteIndex,Format::Format_F>(node);
  Source* SourceObj = Source_Factory->getSource();

  prop_t sq;  //Defines a vector of fields
  CCIO::cout << "Calculating propagator\n";
  QuarkPropDW->calc(sq,*SourceObj);
  //CCIO::ReadFromDisk <Format::Format_F>(sq, "propagator.bin", 12);
  //CCIO::SaveOnDisk <Format::Format_F>(sq, "propagator.bin");

  //Pion Correlator
  MesonCorrelator Meson(Pion);
  std::vector<double> corr =  Meson.calculate <Format::Format_F>(sq,sq);
  for(int i=0; i<corr.size(); ++i) 
    CCIO::cout << i << "  " << corr[i] << std::endl;

  // Cycle among Dirac and color indexes and contract
  // D^-1 * Delta * D^-1
  double mres_numerator = 0;
  double im_check = 0;
  double mres_denominator = 0;
  Field Delta, Denom;

  int Nc = CommonPrms::instance()->Nc();
  int Nd = CommonPrms::instance()->Nd();

  for (int s=0; s<Nd; ++s){
    for (int c=0; c<Nc; ++c){
      //(Delta * D^-1)*source
      Delta = Utils_DWF4D::delta(*(QuarkPropDW->getKernel()),sq[c+Nc*s]); 
      //Contracting 
      mres_numerator += sq[c+Nc*s]*Delta;   // Re(sq[],Delta) sq[]=D^-1*source
      im_check       += sq[c+Nc*s].im_prod(Delta);//should be always zero (test)
      
      //Denominator
      Denom = sq[c+Nc*s];
      Denom -= SourceObj->mksrc(s,c); // (D^-1 - 1)*src
      Denom /= (1.0 - (QuarkPropDW->getKernel()->getMass()) );
      mres_denominator += Denom*Denom;
    }
  }

  CCIO::cout << "Numerator = ("<<mres_numerator<<","<<im_check<<")\n";
  CCIO::cout << "Denominator = " << mres_denominator << std::endl;
  CCIO::cout << "Residual mass = " << mres_numerator/mres_denominator 
	     << std::endl;
  return 0;
}
