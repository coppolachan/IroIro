/*!
 * @file baryonCorrelator.hpp
 * @brief Class for calculation of generic baryon correlator
 * It uses the GammaMatrices namespace
 */
#ifndef BARYON_CORRELATOR_HPP_
#define BARYON_CORRELATOR_HPP_

#include "include/format_F.h"
#include "include/commonPrms.h"
#include "Tools/colorEpsilon.hpp"
#include "Tools/gammaMatrices.hpp"
#include "Communicator/communicator.h"
#include "Main/Geometry/siteIndex.hpp"
#include <vector>
#include <valarray>
#include <complex>

class Field;
typedef std::vector<Field> prop_t;
typedef std::vector<std::complex<double> > correl_t;

enum UpDn{UP=0,DN=1};

class BaryonCorrelator{
private:
  const int Nc_;
  const int Nd_; 
  const int Nep_;
  const int Nvol_;
  const int Nt_;

  Format::Format_F fmt_;
  ColorUtils::ColorEpsilon epsilon_;
  Communicator* comm_;

  const prop_t& Sl_;
  const prop_t& Sh_;
 
  std::complex<double> spTrTr(const prop_t& S1,const prop_t& S2,const prop_t& S3,
			      const GammaMatrices::Gamma& G,int site,UpDn)const;

  std::complex<double> spTr(const prop_t& S1,const prop_t& S2,const prop_t& S3,
			    const GammaMatrices::Gamma& G,int site,UpDn)const;
  
  correl_t global_correl(const std::vector<std::complex<double> >& local);

  const correl_t octet(const prop_t& S1,const prop_t& S2,UpDn);
  const correl_t decuplet_del(const prop_t& S,site_dir,UpDn);
  const correl_t decuplet_sgm(const prop_t& S1,const prop_t& S2,site_dir,UpDn);
public:
  BaryonCorrelator(const prop_t& S_light,const prop_t& S_heavy)
    :Nc_(CommonPrms::instance()->Nc()),
     Nd_(CommonPrms::instance()->Nd()),
     Nep_(Nc_*(Nc_-1)),
     Nvol_(CommonPrms::instance()->Nvol()),
     Nt_(CommonPrms::instance()->Nt()),
     fmt_(Nvol_),
     comm_(Communicator::instance()),
     Sl_(S_light),Sh_(S_heavy){}

  // octet
  const correl_t nucleon(UpDn);
  const correl_t sigma8(UpDn);
  const correl_t xi8(UpDn);
  const correl_t lambda(UpDn);
  // decuplet
  const correl_t delta(site_dir k,UpDn);
  const correl_t omega(site_dir k,UpDn);
  const correl_t sigma10(site_dir k,UpDn);
  const correl_t xi10(site_dir k,UpDn);
  
};

#endif
