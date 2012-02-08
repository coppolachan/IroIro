/*!
 * @file test_ResidualMass.cpp
 *
 * @brief Definition of classes for calculating the Residual Mass
 *
 */
#include "test_ResidualMass.hpp"
#include "include/factories.hpp"
#include "Measurements/FermionicM/qprop_DomainWall.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include "Smearing/StoutSmear.hpp"
#include "Smearing/Smearing_Factories.hpp"
#include "EigenModes/eigenModes_IRL.hpp"
#include "Fields/field_expressions.hpp"
#include "EigenModes/sortEigen.h"
#include <vector>
#include <ctime>

using namespace FieldExpression;


Field Test_ResMass::delta(const Dirac_optimalDomainWall_4D* DWF, const Field& phi){
  //Delta function = 1/4 * (1- sign^2(Hw))
  Field sign = DWF->signKernel(phi);
  Field delta = DWF->signKernel(sign); //sign^2(Hw)
  delta -= phi;  //sign^2(Hw) -1
  delta *= -0.25; // 1/4*(1-sign^2(Hw))
  return delta;
}

int Test_ResMass::run() {
  RNG_Env::RNG = RNG_Env::createRNGfactory(ResMassNode); //create the factory for Rand Numbers Gen

  // First with a thin link
  Staples Staple(conf_.Format);
  CCIO::cout << "Plaquette (thin): " << Staple.plaquette(conf_.U) << std::endl;

  // Smearing 
  Smear* SmearingObj;
  XML::node SmearObjNode = ResMassNode; 
  XML::descend(SmearObjNode, "Smearing");
  SmearingOperatorFactory* Sm_Factory = 
    SmearingOperators::createSmearingOperatorFactory(SmearObjNode);
  // Create smearing objects
  SmearingObj = Sm_Factory->getSmearingOperator();

  smeared_u_ = conf_;

  int Nsmear = 2;
  for (int i = 0; i < Nsmear; i++) {
    previous_u_ = smeared_u_;
    SmearingObj->smear(smeared_u_.U, previous_u_.U);
  }

  /*  
  XML::descend(Smear_node_, "QuarkPropagator");
  QuarkPropagatorFactory* QP_Factory = 
    QuarkPropagators::createQuarkPropagatorFactory(Smear_node_);
  */
  CCIO::cout << "Plaquette (smeared): " << Staple.plaquette(smeared_u_.U) << std::endl;
  
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
  for(int k=0; k<Nmm; ++k) evec[k].resize(ff.size());
  eigen.calc(lmd,evec,b,Nsbt,Nconv);

  CCIO::cout << "Eigenvalues of H_W" << std::endl;
  Field v(ff.size());
  for(int i=0; i<Nsbt+1; ++i){
    v = Hw.mult(evec[i]);
    v -= lmd[i]*evec[i];
    lmd[i] *= (4+mq); // "re"normalize
    CCIO::cout << "["<<i<<"] "<< lmd[i] << "  "<<v*v<<std::endl;
  }
  CCIO::cout << std::endl;

  XML::descend(ResMassNode, "QuarkDWFProp");
  QPropDWFFactory QP_DomainWallFact(ResMassNode);
  //  QpropDWF* QuarkPropDW = static_cast<QpropDWF*>(QP_DomainWallFact.getQuarkProp(conf_));
  QpropDWF* QuarkPropDW = static_cast<QpropDWF*>(QP_DomainWallFact.getQuarkProp(smeared_u_));

  XML::next_sibling(ResMassNode, "Source");
  SourceFactory* Source_Factory = Sources::createSourceFactory<SiteIndex,Format::Format_F>(ResMassNode);
  Source* SourceObj =  Source_Factory->getSource();

  prop_t sq;  //Defines a vector of fields
  CCIO::cout << "Calculating propagator\n";
  QuarkPropDW->calc(sq,*SourceObj);
  //  CCIO::ReadFromDisk < Format::Format_F >(sq, "propagator.bin", 12);
  //  CCIO::SaveOnDisk < Format::Format_F >(sq, "propagator.bin");

  //Pion Correlator
  MesonCorrelator Meson(Pion);
  std::vector<double> corr =  Meson.calculate < Format::Format_F >(sq,sq);
  for(int i = 0; i < corr.size(); ++i) {
    CCIO::cout << i << "  " << corr[i] << std::endl;
  }
 
  // Cycle among Dirac and color indexes and contract
  // D^-1 * Delta * D^-1
  double mres_numerator = 0;
  double im_check = 0;
  double mres_denominator = 0;
  Field Delta(sq[0].size()), Denom(sq[0].size());
  
  for (int s = 0; s < 4; ++s) {
    for (int c = 0; c < 3; ++c) {
      CCIO::cout << "s= " << s << ",  c= " << c << std::endl;
      Delta = delta(QuarkPropDW->getKernel(),sq[c+3*s]); // (Delta * D^-1)*source
      //Contracting 
      mres_numerator += sq[c+3*s]*Delta;   
      // Re(sq[],Delta)    sq[]=D^-1*source
      im_check       += sq[c+3*s].im_prod(Delta);
      //should be always zero (just a check)
      
      //Denominator
      Denom = sq[c+3*s];
      Denom -= SourceObj->mksrc(s,c); // (D^-1 - 1)*src
      Denom /= (1.0 - (QuarkPropDW->getKernel()->getMass()) );
      mres_denominator += Denom*Denom;
    }
  }

  CCIO::cout << "Numerator = ("<<mres_numerator<<","<<im_check<<")\n";
  CCIO::cout << "Denominator = " << mres_denominator << std::endl;
  CCIO::cout << "Residual mass = " << mres_numerator/mres_denominator << std::endl;

  return 0;
}
