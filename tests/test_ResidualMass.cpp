/*!
 * @file test_ResidualMass.cpp
 * @brief Definition of classes for calculating the Residual Mass
 */
#include "test_ResidualMass.hpp"
#include "include/factories.hpp"
#include "Measurements/FermionicM/qprop_DomainWall.hpp"
#include "Measurements/FermionicM/utils_DWF4d.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include "EigenModes/eigenModes_IRL.hpp"
#include "Fields/field_expressions.hpp"
#include "EigenModes/sortEigen.h"
#include <vector>

using namespace FieldExpression;

int Test_ResMass::run() {
//create the factory for Rand Numbers Gen
  RNG_Env::RNG = RNG_Env::createRNGfactory(node_); 

// Prints plaquette (thin) link
  Staples Staple(conf_.Format);
  CCIO::cout << "Plaquette (thin): " << Staple.plaquette(conf_.U) << std::endl;
  ///////////////////////////////////////////////////////////////////////////////////////
  // Smearing objects
  Smear* SmearingObj;         // Empty pointer
  int Nsmear;                 // Number of smearing steps
  XML::node SmearObjNode = node_;  // Copy the node ("descend" function updates it)
                                   // and we want to use again for QuarkPropagator
  XML::descend(SmearObjNode, "Smearing");   // SmearObjNode now points to <Smearing> node
  XML::read(SmearObjNode, "Nsmear", Nsmear, MANDATORY);  // Reads in <Nsmear>

  // Create smearing factory from node information
  SmearingOperatorFactory* Sm_Factory = 
    SmearingOperators::createSmearingOperatorFactory(SmearObjNode);
  // Create smearing objects from the factory
  SmearingObj = Sm_Factory->getSmearingOperator();


  // Copy original configuration to smeared_u_ 
  // smeared_u_ will be passed to the operators
  smeared_u_ = conf_;


  // Do the actual smearing 
  for (int i = 0; i < Nsmear; i++) {
    previous_u_ = smeared_u_;
    SmearingObj->smear(smeared_u_.U, previous_u_.U);
  }
  CCIO::cout << "Plaquette (smeared): " << Staple.plaquette(smeared_u_.U) << std::endl;
  //////////////////////////////////////////////////////////////////////////////////////
  // Eigenvalue calculation
  Format::Format_F ff(CommonPrms::instance()->Nvol());
  double mq  = -1.6;
  Fopr_H Hw(new Dirac_Wilson(mq, &(smeared_u_.U)));
  SortEigen_low sort;
  int    Nk = 20;
  int    Np = 50;
  double enorm = 1.e-22;
  double vthrs = 0.15;
  int    Niter = 500;
  EigenModes_IRL eigen(&Hw,&sort,Nk,Np,enorm,vthrs,Niter);
  Field b(ff.size());

  int Nmm = 100;
  int Nsbt = -1;
  int Nconv = -100;

  std::vector<double> lmd(Nmm);
  std::vector<Field>  evec(Nmm);
  
  for(int k = 0; k < Nmm; ++k) evec[k].resize(ff.size());
  eigen.calc(lmd,evec,b,Nsbt,Nconv);

  CCIO::cout << " --- Eigenvalues of H_W \n " << std::endl;
  Field v(ff.size());
  for(int i = 0; i< Nsbt + 1; ++i){
    v       = Hw.mult(evec[i]);
    v      -= lmd[i]*evec[i];
    lmd[i] *= (4+mq);                // "re"normalize
    CCIO::cout << "["<<i<<"] "<< lmd[i] << "  "<< v*v << "\n";
  }
  CCIO::cout << "\n";
  //////////////////////////////////////////////////////////////////////////////////////

  // Quark Propagator and source creation 
  XML::descend(node_, "QuarkDWFProp");
  QPropDWFFactory QP_DomainWallFact(node_);
  QpropDWF* QuarkPropDW = static_cast<QpropDWF*>(QP_DomainWallFact.getQuarkProp(smeared_u_));

  XML::next_sibling(node_, "Source");
  SourceFactory* Source_Factory 
    = Sources::createSourceFactory<SiteIndex,Format::Format_F>(node_);
  Source* SourceObj =  Source_Factory->getSource();

  prop_t sq;  //Defines a vector of fields
  CCIO::cout << " ---- Calculating propagator\n";
  QuarkPropDW->calc(sq,*SourceObj);
  //  CCIO::ReadFromDisk < Format::Format_F >(sq, "propagator.bin", 12);
  //  CCIO::SaveOnDisk < Format::Format_F >(sq, "propagator.bin");

  //////////////////////////////////////////////////////////////////////////////////////
  //Pion Correlator
  MesonCorrelator Meson(Pion);
  std::vector<double> corr =  Meson.calculate <Format::Format_F>(sq,sq);
  for(int i=0; i<corr.size(); ++i) 
    CCIO::cout << i << "  " << corr[i] << std::endl;

  //////////////////////////////////////////////////////////////////////////////////////
  // Residual mass calculation from quark propagator data

  // Cycle among Dirac and color indexes and contract
  // D^-1 * Delta * D^-1
  double mres_numerator = 0.0;
  double mres_denominator = 0.0;
  double im_check = 0.0;
  int Nc = CommonPrms::instance()->Nc();
  int Nd = CommonPrms::instance()->Nd();

  for (int s=0; s<Nd; ++s){
    for (int c=0; c<Nc; ++c){
      CCIO::cout << "s= " << s << ",  c= " << c << std::endl;
      //(Delta * D^-1)*source
      Field Delta = Utils_DWF4D::delta(*(QuarkPropDW->getKernel()),sq[c+Nc*s]); 
      //Contracting 
      mres_numerator += sq[c+Nc*s]*Delta;   // Re(sq[],Delta) sq[]=D^-1*source
      im_check       += sq[c+Nc*s].im_prod(Delta);//should be always zero (test)
      //Denominator
      Field Denom = sq[c+Nc*s];
      Denom -= SourceObj->mksrc(s,c); // (D^-1 - 1)*src
      Denom /= (1.0 - (QuarkPropDW->getKernel()->getMass()) );
      mres_denominator += Denom*Denom;
    }
  }
  CCIO::cout<<"---------------------------------------------------------"
	    <<std::endl;
  CCIO::cout<<"Numerator = ("   << mres_numerator<<","<<im_check     <<")\n";
  CCIO::cout<<"Denominator = "  << mres_denominator           	     <<"\n";
  CCIO::cout<<"Residual mass = "<< mres_numerator/mres_denominator   <<"\n";
  CCIO::cout<<"---------------------------------------------------------\n";
  ////////////////////////////////////////////////////////////////////////////////////////
  return 0;
}
