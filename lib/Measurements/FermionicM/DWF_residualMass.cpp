/*!
 * @file DWF_residualMass.cpp
 * @brief Definition of class to calculate the Residual Mass
 */
#include "DWF_residualMass.hpp"

double DWFresidualMass::calc() {

  XML::descend(node, "DiracOperator");
  // operator
  // here using a specific factory since we are testing the DWF-4d operator
  DiracDWF4DfullFactory DWF_4d_Factory(node);
  Dirac_DomainWall_4D* DiracDWF_4d = DWF_4d_Factory.getDiracOperator(&(conf_.U));

  //Propagator
  QpropDWF QuarkPropagator(*DiracDWF_4d);

  vector<int> spos(4,0); 
  //Source generator
  Source_local<Format_F> src(spos,CommonPrms::instance()->Nvol());

  prop_t sq;  //Defines a vector of fields
  CCIO::cout << "Calculating propagator\n";
  QuarkPropagator.calc(sq,src);
 
  // Cycle among Dirac and color indexes and contract
  // D^-1 * Delta * D^-1
  double mres_numerator = 0;
  double im_check = 0;
  double mres_denominator = 0;
  Field Delta, Denom;

  int Nd = CommonPrms::instance()->Nd();
  int Nc = CommonPrms::instance()->Nc();

  for(int s =0; s<Nd; ++s){
    for(int c=0; c<Nc; ++c){
      Delta = delta(DiracDWF_4d,sq[c+Nc*s]); // (Delta * D^-1)*source
      //Contracting 
      mres_numerator += sq[c+Nc*s]*Delta;         //Re(sq[],Delta)    sq[]=D^-1*source
      im_check       += sq[c+Nc*s].im_prod(Delta);//should be always zero (just a check)
      CCIO::cout<< "Numerator = ("<<mres_numerator<<","<<im_check<<")\n";
      
      //Denominator
      Denom = sq[c+Nc*s];
      Denom -= src.mksrc(s,c); // (D^-1 - 1)*src
      Denom /= (1.0 -DiracDWF_4d->getMass());
      mres_denominator += Denom*Denom;
      CCIO::cout << "Denominator = " << mres_denominator << endl;
      CCIO::cout << "Residual mass = " << mres_numerator/mres_denominator << endl;
    }
  }
  return 0;
}
