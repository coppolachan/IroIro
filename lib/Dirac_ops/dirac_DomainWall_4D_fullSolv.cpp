/*!
 * @file dirac_DomainWall_4D_fullSolv.cpp
 * @brief Declaration of Dirac_DomainWall_4D_fullSolv class 
 Time-stamp: <2013-12-05 10:03:35 noaki>
 */
#include "dirac_DomainWall_4D_fullSolv.hpp"
#include "Fields/field_expressions.hpp"

using namespace std;
using namespace FieldExpression;

void Dirac_DomainWall_4D_fullSolv::mult_std(Field& w,const Field& f)const{
  // Dpv_^-1
  Field src = Dpv_->mult_dag(Dodw_->mult(Bproj_dag(f)));
  Field sol5(Dodw_->fsize());
  SolverOutput monitor = slv_pv_->solve(sol5,src);
#if VERBOSITY > BASE_VERB_LEVEL
  monitor.print();
#endif
  w = Bproj(sol5);
}

void Dirac_DomainWall_4D_fullSolv::mult_inv_std(Field& w,const Field& f)const{
  // D_dw^-1
  Field src = Dodw_->mult_dag(Dpv_->mult(Bproj_dag(f)));
  Field sol5(Dodw_->fsize());
  SolverOutput monitor = slv_odw_->solve(sol5,src);
#if VERBOSITY > BASE_VERB_LEVEL
  monitor.print();
#endif
  w = Bproj(sol5);
}

void Dirac_DomainWall_4D_fullSolv::mult_LU(Field& w,const Field& f)const{
  // Dpv_^-1
  Field src = Dpv_->mult_hop5_dinv(Dpv_->mult_dag(Dodw_->mult(Bproj_dag(f))));
  Field sol5(Dodw_->fsize());
  SolverOutput monitor = slv_pv_->solve(sol5,src);
#if VERBOSITY > BASE_VERB_LEVEL
  monitor.print();
#endif
  w = Bproj(Dpv_->mult_hop5_inv(sol5));
}

void Dirac_DomainWall_4D_fullSolv::mult_inv_LU(Field& w,const Field& f)const{
  // D_dw^-1
  Field src = Dodw_->mult_hop5_dinv(Dodw_->mult_dag(Dpv_->mult(Bproj_dag(f))));
  Field sol5(Dodw_->fsize());
  SolverOutput monitor = slv_odw_->solve(sol5,src);
#if VERBOSITY > BASE_VERB_LEVEL
  monitor.print();
#endif
  w = Bproj(Dodw_->mult_hop5_inv(sol5));
}

const Field Dirac_DomainWall_4D_fullSolv::mult(const Field& f)const{
  Field w(fsize_);
  (this->*mult_core)(w,f);
  return w;
}

const Field Dirac_DomainWall_4D_fullSolv::mult_inv(const Field& f)const{
  Field w(fsize_);
  (this->*mult_inv_core)(w,f);
  return w;
}

const Field Dirac_DomainWall_4D_fullSolv::mult_dag(const Field& f)const{
  return gamma5(mult(gamma5(f)));}

const Field Dirac_DomainWall_4D_fullSolv::mult_dag_inv(const Field& f)const{
  return gamma5(mult_inv(gamma5(f)));}

const Field Dirac_DomainWall_4D_fullSolv::gamma5(const Field& f)const{ 
  Field w(Dodw_->f4size());
  Dodw_->gamma5_4d(w,f);
  return w;
}

const Field Dirac_DomainWall_4D_fullSolv::Bproj(const Field& f5) const{ 
  Field f4(fsize_);
  int N5 = Dodw_->fsize()/fsize_;
  int Nvol = fsize_/ffmt_t::Nin();
  ffmt_t ff(Nvol);
  for(int site=0; site<Nvol; ++site)
    d5_.BprojCore(f4.getaddr(ff.index(0,site)),
		  const_cast<Field&>(f5).getaddr(ff.index(0,site,0)),
		  const_cast<Field&>(f5).getaddr(ff.index(0,site,N5-1)));
  return f4;
}

const Field Dirac_DomainWall_4D_fullSolv::Bproj_dag(const Field& f4) const{
  int N5 = Dodw_->fsize()/fsize_;
  int Nvol = fsize_/ffmt_t::Nin();
  ffmt_t ff(Nvol);
  Field f5(Dodw_->fsize());
  for(int site=0; site<Nvol; ++site)
    d5_.BprojCore_dag(f5.getaddr(ff.index(0,site,0)),
		      f5.getaddr(ff.index(0,site,N5-1)),
		      const_cast<Field&>(f4).getaddr(ff.index(0,site)));
  return f5;
}
