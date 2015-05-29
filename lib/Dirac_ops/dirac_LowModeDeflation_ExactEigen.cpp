/*!@file  dirac_LowModeDeflation_ExactEigen.cpp
 * @brief methods for low-mode preconditioning with exact low-modes
 *Time-stamp: <2014-10-05 11:14:14 noaki>
 */
#include "dirac_LowModeDeflation_ExactEigen.hpp"
#include "EigenModes/subSpaceProjector.hpp"

const Field Dirac_LowModeDeflation_ExactEigen::mult(const Field& f)const{
  return Dw_->mult(mult_left(f));
}

const Field Dirac_LowModeDeflation_ExactEigen::mult_dag(const Field& f)const{
  return gamma5(mult(gamma5(f)));
}

const Field Dirac_LowModeDeflation_ExactEigen::mult_little(const Field& f)const{
  Field Dfl;
  SubSpace::project_complex(Dfl,f,*evecs_,*evals_);
  return Dfl;
}

const Field Dirac_LowModeDeflation_ExactEigen::mult_little_dag(const Field& f)const{
  Field Dfl;
  SubSpace::project_complex(Dfl,f,*evecs_,evals_conj_);
  return Dfl;
}

const Field Dirac_LowModeDeflation_ExactEigen::mult_little_inv(const Field& f)const{
  Field Difl;
  SubSpace::project_complex(Difl,f,*evecs_,ievals_);
  return Difl;
}

const Field Dirac_LowModeDeflation_ExactEigen::mult_little_inv_dag(const Field& f)const{
  Field Difl;
  SubSpace::project_complex(Difl,f,*evecs_,ievals_conj_);
  return Difl;
}

const Field Dirac_LowModeDeflation_ExactEigen::mult_left(const Field& f)const{
  Field fh;
  SubSpace::projectOut(fh,f,*evecs_);
  return fh;
}

void Dirac_LowModeDeflation_ExactEigen::setEigenArrays(){
  int Neigen = evecs_->size();
  assert(Neigen == evals_->size()/2);

  evals_conj_= *evals_;

  for(int i=0; i<Neigen; ++i){
    double Nev = (*evals_)[2*i]*(*evals_)[2*i] +(*evals_)[2*i+1]*(*evals_)[2*i+1];
    evals_conj_[2*i+1] *= -1.0;

    ievals_[2*i]   =  (*evals_)[2*i]/Nev;
    ievals_[2*i+1] = -(*evals_)[2*i+1]/Nev;
    
    ievals_conj_[2*i]   = (*evals_)[2*i]/Nev;
    ievals_conj_[2*i+1] = (*evals_)[2*i+1]/Nev;
  }
}
