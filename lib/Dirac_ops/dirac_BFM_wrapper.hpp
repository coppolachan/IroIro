/*!
 * @file dirac_BFM_wrapper.hpp
 * @brief Defines the wrapper classs for P. Boyle Bagel/BFM libs
 * Time-stamp: <2013-07-03 10:21:24 cossu>
 */
#ifndef DIRAC_BFM_WRAPPER_
#define DIRAC_BFM_WRAPPER_

/*
 *! @class Dirac_BFM_Wrapper
 * @brief Declaration of abstract base class for Dirac operators
 */
class Dirac_BFM_Wrapper: public DiracWilsonLike_EvenOdd {
public:
  Dirac_BFM_Wrapper();
  ~Dirac_BFM_Wrapper(){}
  size_t fsize() const;
  size_t gsize() const;

  const Field mult    (const Field&)const = 0;
  const Field mult_dag(const Field&)const = 0;

  const Field md_force(const Field& eta,const Field& zeta)const{};
  void update_internal_state(){};
};
