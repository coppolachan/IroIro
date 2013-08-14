/*!
  @file smartConf_base.cpp
  @brief Defines the SmartConf class member base functions 
*/
#include "smartConf.hpp"
#include "include/messages_macros.hpp"

typedef std::complex<double> dcomplex;

//====================================================================
/*! @brief Returns top level configuration */
GaugeField* SmartConf::get_current_conf() const{
  return const_cast<GaugeField*>(&(SmearedSet[smearingLevels-1]));
}
//====================================================================
/*! @brief Returns smeared configuration at level 'Level' */
const GaugeField& SmartConf::get_smeared_conf(int Level) const{
  return SmearedSet[Level];
}
//====================================================================
void SmartConf::fill_smearedSet(){
  GaugeField previous_u;
 
  _Message(DEBUG_VERB_LEVEL, "[SmartConf] Filling SmearedSet\n");
 
  previous_u = *ThinLinks;
  for(int smearLvl = 0; smearLvl < smearingLevels; ++smearLvl){
    StoutSmearing.smear(SmearedSet[smearLvl],previous_u);
    previous_u = SmearedSet[smearLvl];
  }
}
//====================================================================
void SmartConf::set_fj(dcomplex& f0, dcomplex& f1, dcomplex& f2,
		       const double& u, const double& w)const{
  double xi0 = func_xi0(w);
  double u2 = u*u;
  double w2 = w*w;
  double cosw = cos(w);

  dcomplex ixi0 = dcomplex(0.0,xi0);
  dcomplex emiu = dcomplex(cos(u),-sin(u));
  dcomplex e2iu = dcomplex(cos(2.0*u),sin(2.0*u));

  dcomplex h0 = e2iu * dcomplex(u2-w2,0.0) 
    + emiu*(dcomplex(8.0*u2*cosw,0.0) +dcomplex(2.0*u*(3.0*u2 + w2),0.0)*ixi0);
  
  dcomplex h1 = dcomplex(2*u,0.0) * e2iu 
    - emiu*(dcomplex(2.0*u*cosw,0.0) -dcomplex(3.0*u2-w2,0.0)*ixi0);

  dcomplex h2 = e2iu -emiu*(dcomplex(cosw,0.0) +dcomplex(3.0*u,0.0)*ixi0);

  double fden = 1.0/(9.0*u2 -w2);

  f0 = h0 * dcomplex(fden,0.0);
  f1 = h1 * dcomplex(fden,0.0);
  f2 = h2 * dcomplex(fden,0.0);
}

//====================================================================
double SmartConf::func_xi0(double w)const{

  double xi0 = sin(w)/w;
  if(w<0.0001) CCIO::cout<<"w too small: "<<w<<"\n";
  return xi0;
}

//====================================================================
double SmartConf::func_xi1(double w)const{

  double xi1 = cos(w)/(w*w) - sin(w)/(w*w*w);
  if(w<0.0001) CCIO::cout<<"w too small: "<<w<<"\n";

  return xi1;
}
