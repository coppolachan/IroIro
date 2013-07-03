/*! @file dirac_DWoverlap.cpp
 *  @brief methods of the Dirac_DomainWall_4D class
 * Time-stamp: <2013-07-03 10:42:24 noaki>
 */
#include "dirac_DWoverlap.hpp"
#include "subSpaceProjector.hpp"
#include "Fields/field_expressions.hpp"

using namespace FieldExpression;

const Field Dirac_DWoverlap::mult(const Field& f)const{
  Field fp;
  SubSpace::project_real(fp,f,ems_->evecs_,sgnEv_); // sgn(lmd)*f_low

  Field w = gamma5(fp);          
  w *= 0.5*(1.0-mq_);

  SubSpace::project(fp,f,ems_->evecs_);    
  w += 0.5*(1.0+mq_)*fp;                              
  fp -= f;
  w -= Ddw4d_->mult(fp);                        // D4*f_high

  return w;
}

const Field Dirac_DWoverlap::mult_dag(const Field& f)const{ 
  return gamma5(mult(gamma5(f)));}


