/*!--------------------------------------------------------------------------
 * @file dirac_DomainWall_EvenOdd.cpp
 *
 * @brief Definition of class methods for Dirac_optimalDomainWall_EvenOdd (5d op.)
 *-------------------------------------------------------------------------*/
#include "dirac_DomainWall_EvenOdd.hpp"
#include "Communicator/comm_io.hpp"
#include<stdlib.h>
#include<stdio.h>
#include<cassert>
#include<math.h>

using namespace std;

//-----------------------------------------------------------------------------
const Field Dirac_optimalDomainWall_EvenOdd::mult_ee(const Field& f)const{
  return Deo_.mult_hop5(f);
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_oo(const Field& f)const{
  return Deo_.mult_hop5(f);
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_ee_inv(const Field& f)const{
  return Deo_.mult_hop5_inv(f);
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_ee_dinv(const Field& f)const{
  return Deo_.mult_hop5_dinv(f);
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_oo_inv(const Field& f)const{
  return Deo_.mult_hop5_inv(f);
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_oo_dinv(const Field& f)const{
  return Deo_.mult_hop5_dinv(f);
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_eo(const Field& f)const{
  return Deo_.mult_hop5_inv(Deo_.mult(f));
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_oe(const Field& f)const{
  return Doe_.mult_hop5_inv(Doe_.mult(f));
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_eo_dag(const Field& f)const{
  return Doe_.mult_dag(Deo_.mult_hop5_dinv(f));
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_oe_dag(const Field& f)const{
  return Deo_.mult_dag(Doe_.mult_hop5_dinv(f));
}

const Field Dirac_optimalDomainWall_EvenOdd::mult(const Field& f) const{
  
  Field w(f);
  w -= mult_eo(mult_oe(f));
  return w;
  
  
  /*
  // just slightly faster (only BGQ)
  //  Field w5(Deo_.fsize());
  Field res(Deo_.fsize());
  Doe_.mult_hop(res,f);
  //Deo_.mult_hop(res,w5);
  //res -= f;

  return res;
  */
}
const Field Dirac_optimalDomainWall_EvenOdd::mult_dag(const Field& f) const{
  Field w(f);
  w -= mult_oe_dag(mult_eo_dag(f));
  return w;
}

void Dirac_optimalDomainWall_EvenOdd::
md_force_eo(Field& fce, const Field& eta,const Field& zeta) const{
  Deo_.md_force_p(fce,eta,mult_ee_dinv(zeta));
  Doe_.md_force_m(fce,eta,mult_ee_dinv(zeta));
}

void Dirac_optimalDomainWall_EvenOdd::
md_force_oe(Field& fce, const Field& eta,const Field& zeta) const{
  Doe_.md_force_p(fce,eta,mult_oo_dinv(zeta));
  Deo_.md_force_m(fce,eta,mult_oo_dinv(zeta));
}

const Field Dirac_optimalDomainWall_EvenOdd::
md_force(const Field& eta,const Field& zeta) const{
  Field fce(gsize());
  md_force_eo(fce,mult_oe(eta),zeta);
  md_force_oe(fce,eta,mult_eo_dag(zeta));
  fce *= 0.5;
  return fce;
}
