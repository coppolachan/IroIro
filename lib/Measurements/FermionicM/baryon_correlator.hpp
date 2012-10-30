/*!
 * @file baryon_correlator.hpp
 * @brief Class for calculation of generic baryon correlator
 * It uses the GammaMatrices namespace
 */
#ifndef BARYON_CORRELATOR_HPP_
#define BARYON_CORRELATOR_HPP_

#include <vector>

class Field;
typedef std::vector<Field> prop_t;

class BaryonCorrelator{
private:
  const int Nc_;
  const int Nd_; 
  const int Nep_;
  const int Nvol_;
  const int Nt_;

  Format::Format_F fmt_;
  ColorEpsilon epsilon_;
  Communicator* comm_;

  const prop_t& Sl_;
  const prop_t& Sh_;
 
  std::complex<double> spTrTr(const prop_t& S1,
			      const prop_t& S2,
			      const prop_t& S3,
			      const GammaMatrices::Gamma& G,
			      int site,int i,int j) const;

  std::complex<double> spTr(const prop_t& S1,
			    const prop_t& S2,
			    const prop_t& S3,
			    const GammaMatrices::Gamma& G,
			    int site,int i,int j) const;

  void global_correl(std::vector<complex<double> >& global, 
		     const std::vector<complex<double> >& local);

  const std::vector<std::complex<double> > octet(const prop_t& S1,
						 const prop_t& S2);
  const std::vector<std::complex<double> > decuplet_del(const prop_t& S,int k);
  const std::vector<std::complex<double> > decuplet_sgm(const prop_t& S1,
							const prop_t& S2,int k);
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
  const std::vector<std::complex<double> > nucleon();
  const std::vector<std::complex<double> > sigma8();
  const std::vector<std::complex<double> > xi8();
  const std::vector<std::complex<double> > lambda();
  // decuplet
  const std::vector<std::complex<double> > delta(int k);
  const std::vector<std::complex<double> > omega(int k);
  const std::vector<std::complex<double> > sigma10(int k);
  const std::vector<std::complex<double> > xi10(int k);
  
};


