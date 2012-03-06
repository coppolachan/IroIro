/*!
  @file SmartConf.cpp

  @brief Defines the SmartConf class member functions
*/
#include "SmartConf.hpp"
#include "Main/Geometry/mapper.hpp"
#include "Tools/sunMatUtils.hpp"
#include "include/messages_macros.hpp"

typedef complex<double> dcomplex;

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
void SmartConf::smeared_force(GaugeField& SigmaTilde)const {
  using namespace FieldUtils;
  using namespace SUNmatUtils;

  GaugeField force;
  SUNmat ut;

  int Nvol = CommonPrms::Nvol();
  
  force = SigmaTilde;//actually = U*SigmaTilde

  for(int mu = 0; mu < NDIM_; ++mu){
    for (int site = 0; site < Nvol; ++site){
      SetMatrix(force,
		matrix_dag(SmearedSet[smearingLevels-1],site,mu) * matrix(force,site,mu),
		site,
		mu);
    }
  } 
  
  for(int ismr = smearingLevels - 1; ismr > 0; --ismr){
    force = AnalyticSmearedForce(force,get_smeared_conf(ismr));
  }
  
  force = AnalyticSmearedForce(force,*ThinLinks);
  
  for(int mu = 0; mu < NDIM_; ++mu){
    for (int site = 0; site < Nvol; ++site){
      SetMatrix(SigmaTilde, 
		matrix(*ThinLinks,site,mu) * matrix(force,site,mu),
		site,
		mu);
    }
  }  
}
//====================================================================
GaugeField SmartConf::AnalyticSmearedForce(const GaugeField& SigmaKPrime,
					   const GaugeField& GaugeK) const{
  using namespace SUNmatUtils;
  using namespace FieldUtils;

  int Nvol = CommonPrms::Nvol();
  
  GaugeField iQ, e_iQ;
  GaugeField C, iLambda;
  GaugeField SigmaK;
  
  StoutSmearing.BaseSmear(C,GaugeK);

  for(int mu = 0; mu < NDIM_; ++mu){
    for(int site = 0; site < Nvol; ++site){
      SetMatrix(iQ, 
		anti_hermite_traceless(matrix(C,site,mu) * matrix_dag(GaugeK,site,mu)), 
		site, 
		mu);
    }
  }// created iQ
  
  set_iLambda(iLambda, e_iQ, iQ, SigmaKPrime, GaugeK);
  
  for(int mu = 0; mu < NDIM_; ++mu){
    for (int site = 0; site < Nvol; ++site){
      SetMatrix(SigmaK,matrix(SigmaKPrime,site,mu) * matrix(e_iQ,site,mu), site, mu);
      AddMatrix(SigmaK,matrix_dag(C,site,mu) * matrix(iLambda,site,mu)   , site, mu);
    }
  }

  StoutSmearing.derivative(SigmaK, iLambda, GaugeK);
  
  return SigmaK;
}
//====================================================================
void SmartConf::set_iLambda(GaugeField& iLambda, 
			    GaugeField& e_iQ,
			    const GaugeField& iQ,
			    const GaugeField& Sigmap, 
			    const GaugeField& GaugeK)const{
  using namespace SUNmatUtils;
  using namespace FieldUtils;

  int Nvol = CommonPrms::Nvol();

  SUNmat iQ0, iQ1, iQ2, iQ3;
  SUNmat B1, B2;
  SUNmat USigmap,iQUS,iUSQ,iGamma;
  
  iQ0 = unity();

  double u_, w, u2, w2, cosw, xi0, xi1, fden;
  dcomplex f0, f1, f2, h0, h1, h2, e2iu, emiu, qt;
  dcomplex r01, r11, r21, r02, r12, r22, tr1, tr2;
  dcomplex b10, b11, b12, b20, b21, b22;

  for(int mu = 0; mu < NDIM_; ++mu){
    for(int site = 0; site < Nvol; ++site){

      iQ1 = matrix(iQ,site,mu);
      iQ2 = iQ1 * iQ1;
      iQ3 = iQ1 * iQ2;

      set_uw(u_,w,iQ2,iQ3);
      set_fj(f0,f1,f2,u_,w);

      for(int cc = 0; cc < NC_*NC_; ++cc){
	qt =  f0 * dcomplex(iQ0.r(cc), iQ0.i(cc))
          + f1 * dcomplex(iQ1.i(cc),-iQ1.r(cc))
          - f2 * dcomplex(iQ2.r(cc), iQ2.i(cc));
	e_iQ.data.set(e_iQ.format.index(2*cc,site,mu),  qt.real());
        e_iQ.data.set(e_iQ.format.index(2*cc+1,site,mu),qt.imag());
      }

      xi0 = func_xi0(w);
      xi1 = func_xi1(w);
      u2 = u_ * u_;
      w2 = w * w;
      cosw = cos(w);

      emiu = dcomplex(cos(u_),-sin(u_));
      e2iu = dcomplex(cos(2.0*u_),sin(2.0*u_));

      r01 = dcomplex(2.0*u_,2.0*(u2-w2)) * e2iu
	+ emiu * dcomplex(16.0*u_*cosw + 2.0*u_*(3.0*u2+w2)*xi0,
			  -8.0*u2*cosw + 2.0*(9.0*u2+w2)*xi0);

      r11 = dcomplex(2.0,4.0*u_) * e2iu
	+ emiu * dcomplex(-2.0*cosw + (3.0*u2-w2)*xi0,
			  2.0*u_*cosw + 6.0*u_*xi0);

      r21 = dcomplex(0.0,2.0) * e2iu
	+ emiu * dcomplex(-3.0*u_*xi0, cosw - 3.0*xi0);

      r02 = dcomplex(-2.0,0.0) * e2iu
	+ emiu * dcomplex(-8.0*u2*xi0,
			  2.0*u_*(cosw + xi0 + 3.0*u2*xi1));

      r12 = emiu * dcomplex(2.0*u_*xi0,
			    -cosw - xi0 + 3.0*u2*xi1);

      r22 = emiu * dcomplex(xi0, -3.0*u_*xi1);

      fden = 1.0/(2*(9.0*u2-w2)*(9.0*u2-w2));

      b10 =  dcomplex(2.0*u_, 0.0)*r01 + dcomplex(3.0*u2-w2, 0.0)*r02
	- dcomplex(30.0*u2+2.0*w2, 0.0)*f0;

      b11 =  dcomplex(2.0*u_, 0.0)*r11 + dcomplex(3.0*u2-w2, 0.0)*r12
	- dcomplex(30.0*u2+2.0*w2, 0.0)*f1;

      b12 =  dcomplex(2.0*u_, 0.0)*r21 + dcomplex(3.0*u2-w2, 0.0)*r22
	- dcomplex(30.0*u2+2.0*w2, 0.0)*f2;

      b20 = r01 - dcomplex(3.0*u_, 0.0)*r02 - dcomplex(24.0*u_, 0.0)*f0;

      b21 = r11 - dcomplex(3.0*u_, 0.0)*r12 - dcomplex(24.0*u_, 0.0)*f1;

      b22 = r21 - dcomplex(3.0*u_, 0.0)*r22 - dcomplex(24.0*u_, 0.0)*f2;

      b10 *= dcomplex(fden,0.0);
      b11 *= dcomplex(fden,0.0);
      b12 *= dcomplex(fden,0.0);
      b20 *= dcomplex(fden,0.0);
      b21 *= dcomplex(fden,0.0);
      b22 *= dcomplex(fden,0.0);

      for(int cc = 0; cc < NC_*NC_; ++cc){
	qt =  b10 * dcomplex(iQ0.r(cc), iQ0.i(cc))
          + b11 * dcomplex(iQ1.i(cc),-iQ1.r(cc))
          - b12 * dcomplex(iQ2.r(cc), iQ2.i(cc));
	B1.set(cc,qt.real(),qt.imag());
	qt =  b20 * dcomplex(iQ0.r(cc), iQ0.i(cc))
          + b21 * dcomplex(iQ1.i(cc),-iQ1.r(cc))
          - b22 * dcomplex(iQ2.r(cc), iQ2.i(cc));
	B2.set(cc,qt.real(),qt.imag());
      }

      USigmap = matrix(GaugeK,site,mu) * matrix(Sigmap,site,mu);

      tr1 = dcomplex(ReTr(USigmap*B1),ImTr(USigmap*B1));
      tr2 = dcomplex(ReTr(USigmap*B2),ImTr(USigmap*B2));

      iQUS = iQ1 * USigmap;
      iUSQ = USigmap * iQ1;

      for(int cc = 0; cc < NC_*NC_; ++cc){
	qt =  tr1 * dcomplex(iQ1.i(cc),-iQ1.r(cc))
          - tr2 * dcomplex(iQ2.r(cc), iQ2.i(cc))
          + f1  * dcomplex(USigmap.r(cc), USigmap.i(cc))
          + f2  * dcomplex(iQUS.i(cc),-iQUS.r(cc))
          + f2  * dcomplex(iUSQ.i(cc),-iUSQ.r(cc));
	iGamma.set(cc,-qt.imag(),qt.real());
      }
      SetMatrix(iLambda,anti_hermite_traceless(iGamma), site, mu);
    }
  }

}

//====================================================================
void SmartConf::set_fj(dcomplex& f0, dcomplex& f1, dcomplex& f2,
		       const double& u, const double& w)const{

  dcomplex h0, h1, h2, e2iu, emiu, ixi0, qt;
  double c0, c1, xi0, c0max, u2, w2, cosw, fden;

  xi0 = func_xi0(w);
  u2 = u * u;
  w2 = w * w;
  cosw = cos(w);

  ixi0 = dcomplex(0.0,xi0);
  emiu = dcomplex(cos(u),-sin(u));
  e2iu = dcomplex(cos(2.0*u),sin(2.0*u));

  h0 = e2iu * dcomplex(u2-w2,0.0)
    + emiu * (dcomplex(8.0*u2*cosw,0.0)
	      + dcomplex(2.0*u*(3.0*u2 + w2),0.0)*ixi0);
  h1 = dcomplex(2*u,0.0) * e2iu
    - emiu*(dcomplex(2.0*u*cosw,0.0)
	    - dcomplex(3.0*u2-w2,0.0) * ixi0);
  h2 = e2iu - emiu * (dcomplex(cosw,0.0)
		      + dcomplex(3.0*u,0.0)*ixi0);

  fden = 1.0/(9.0*u2 - w2);
  f0 = h0 * dcomplex(fden,0.0);
  f1 = h1 * dcomplex(fden,0.0);
  f2 = h2 * dcomplex(fden,0.0);

}

//====================================================================
void SmartConf::set_uw(double& u, double& w,
		       const SUNmat& iQ2,
		       const SUNmat& iQ3)const{

  double c0, c1, c0max, theta;

  c0 = 0.0;
  c1 = 0.0;
  for(int cc = 0; cc < NC_; ++cc){
    c0 += iQ3.i(cc,cc);
    c1 += iQ2.r(cc,cc);
  }
  c0 = -c0/3.0;
  c1 = -c1/2.0;
  c0max = 2.0 * pow(c1/3.0,1.5);

  theta = acos(c0/c0max);
  u = sqrt(c1/3.0) * cos(theta/3.0);
  w = sqrt(c1) * sin(theta/3.0);

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
//====================================================================
//============================================================END=====
