/*!
 * @file test_ResidualMass.cpp
 * @brief Definition of classes for calculating the Residual Mass
 */
#include "test_ResidualMass.hpp"
#include "include/factories.hpp"
#include "Measurements/FermionicM/qprop_DomainWall.hpp"
#include "Measurements/FermionicM/quark_prop_meas_factory.hpp"
#include "Measurements/FermionicM/utils_DWF4d.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include "EigenModes/eigenModesSolver_IRL.hpp"
#include "Fields/field_expressions.hpp"
#include "EigenModes/eigenSorter.hpp"
#include <vector>

using namespace FieldExpression;

int Test_ResMass::run() {
  //create the factory for Rand Numbers Gen

  //RNG_Env::RNG = RNG_Env::createRNGfactory(node_); 
  RNG_Env::initialize(node_);

  // Prints plaquette (thin) link
  Staples Staple;
  CCIO::cout << "Plaquette (thin): " << Staple.plaquette(conf_) << std::endl;

  ///////////////////////////////////////////////////////////////////////////////
  // Smearing objects
  Smear* SmearingObj;         // Empty pointer
  int Nsmear;                 // Number of smearing steps
  XML::node SmearObjNode = node_;// Copy the node ("descend" function updates it)
                                 // and we want to use again for QuarkPropagator
  XML::descend(SmearObjNode, "Smearing");//SmearObjNode now points to <Smearing> node
  XML::read(SmearObjNode, "Nsmear", Nsmear, MANDATORY);  // Reads in <Nsmear>

  // Create smearing factory from node information
  SmearingFactory* Sm_Factory = 
    Smearings::createSmearingFactory(SmearObjNode);
  // Create smearing objects from the factory
  SmearingObj = Sm_Factory->getSmearing();

  // Copy original configuration to smeared_u_ 
  // smeared_u_ will be passed to the operators
  smeared_u_ = conf_;

  // Do the actual smearing 
  for(int i=0; i<Nsmear; i++) {
    previous_u_= smeared_u_;
    SmearingObj->smear(smeared_u_,previous_u_);
  }
  CCIO::cout<<"Plaquette (smeared): "<< Staple.plaquette(smeared_u_)<<std::endl;
  //////////////////////////////////////////////////////////////////////////////
  // Eigenvalue calculation
  XML::node eigen_node= node_;
  XML::descend(eigen_node, "EigenModes");
  const char* eigen= eigen_node.attribute("exec").value();
  if(!strcmp(eigen,"DoCalc")){

    double mq,vthrs,enorm;
    int Nk,Np,Niter;

    // names are obscure
    XML::read(eigen_node,"mq",mq,MANDATORY);
    XML::read(eigen_node,"vthrs",vthrs,MANDATORY);
    XML::read(eigen_node,"enorm",enorm,MANDATORY);
    XML::read(eigen_node,"Nk",Nk,MANDATORY);
    XML::read(eigen_node,"Np",Np,MANDATORY);
    XML::read(eigen_node,"Niter",Niter,MANDATORY);

    Format::Format_F ff(CommonPrms::instance()->Nvol());
    Fopr_H Hw(new Dirac_Wilson(mq, &(smeared_u_.data)));

    CCIO::cout << "Calculating eigenvalues of H_W" << std::endl;
    CCIO::cout << " -- M0= " << mq << std::endl;
    CCIO::cout << " -- vthrs= " << vthrs << std::endl;
    EigenSorter_low esort(vthrs);
    EigenModesSolver_IRL eigen(&Hw,&esort,Nk,Np,enorm,Niter);
    
    int Nmm = 100;
    int Nsbt = -1;
    
    std::vector<double> lmd(Nmm);
    std::vector<Field>  evec(Nmm);
    
    for(int k = 0; k < Nmm; ++k) evec[k].resize(ff.size());
    eigen.calc(lmd,evec,Nsbt);
    
    CCIO::cout << " --- Eigenvalues of H_W" << std::endl;
    Field v(ff.size());
    for(int i = 0; i< Nsbt + 1; ++i){
      v       = Hw.mult(evec[i]);
      v      -= lmd[i]*evec[i];
      lmd[i] *= (4+mq);                // "re"normalize
      CCIO::cout << "["<<i<<"] "<< lmd[i] << "  "<< v*v << "\n";
    }
    CCIO::cout << std::endl;
  }else if(!strcmp(eigen,"NoCalc")){ }

  /////////////////////////////////////////////////////////////////////////////

  // Quark Propagator and source creation 
  XML::descend(node_, "QuarkDWFProp");
  QPropDWFFactory QP_DomainWallFact(node_);

  InputConfig input(&smeared_u_);
  QpropDWF* QuarkPropDW 
    = static_cast<QpropDWF*>(QP_DomainWallFact.getQuarkProp(input));

  XML::next_sibling(node_, "Source");
  SourceFactory* Source_Factory 
    = Sources::createSourceFactory<SiteIndex,Format::Format_F>(node_);
  Source* SourceObj =  Source_Factory->getSource();

  prop_t sq;  //Defines a vector of fields
  CCIO::cout << " ---- Calculating propagator\n";
  QuarkPropDW->calc(sq,*SourceObj);
  //  CCIO::ReadFromDisk < Format::Format_F >(sq, "propagator.bin", 12);
  //  CCIO::SaveOnDisk < Format::Format_F >(sq, "propagator.bin");

  /////////////////////////////////////////////////////////////////////////////
  //Pion Correlator
  MesonCorrelator Meson(Pion);
  std::vector<double> corr =  Meson.calculate <Format::Format_F>(sq,sq);
  for(int i=0; i<corr.size(); ++i) 
    CCIO::cout << i << "  " << corr[i] << std::endl;

  //////////////////////////////////////////////////////////////////////////////
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
  CCIO::cout<<"Numerator = ("   <<mres_numerator<<","<<im_check<<")"<<std::endl;
  CCIO::cout<<"Denominator = "  <<mres_denominator           	    <<std::endl;
  CCIO::cout<<"Residual mass = "<<mres_numerator/mres_denominator   <<std::endl;
  CCIO::cout<<"---------------------------------------------------------"
	    <<std::endl;
  ////////////////////////////////////////////////////////////////////////////////
  return 0;
}
