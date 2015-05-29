/*!@file  dirac_LowModeDeflation_Approx.cpp
 * @brief methods for low-mode preconditioning with given subspace
 *Time-stamp: <2014-10-05 11:14:46 noaki>
 */
#include "dirac_LowModeDeflation_Approx.hpp"
#include "EigenModes/subSpaceProjector.hpp"
#include "Tools/complexMatrix.hpp"
//#include "Fields/field_expressions.hpp"
//using namespace FieldExpression;

using namespace std;

const Field Dirac_LowModeDeflation_Approx::mult(const Field& f)const{
  Field w; 
  rightPrec(w,f);
  return mult_left(Dw_->mult(w));
}

const Field Dirac_LowModeDeflation_Approx::mult_dag(const Field& f)const{
  return gamma5(mult(gamma5(f)));
}

const Field Dirac_LowModeDeflation_Approx::mult_little(const Field& f)const{
  Field DPDf(f.size());
  lowModeProj(DPDf,f,Dsub_);
  return DPDf;
}


const Field Dirac_LowModeDeflation_Approx::mult_little_dag(const Field& f)const{
  Field DPDf(f.size());
  lowModeProj(DPDf,f,Dsub_dag_);
  return DPDf;
}

const Field Dirac_LowModeDeflation_Approx::mult_little_inv(const Field& f)const{
  Field DPDf(f.size());
  lowModeProj(DPDf,f,Dsub_inv_);
  return DPDf;
}

const Field Dirac_LowModeDeflation_Approx::mult_little_inv_dag(const Field& f)const{
  Field DPDf(f.size());
  lowModeProj(DPDf,f,Dsub_dag_inv_);
  return DPDf;
}

const Field Dirac_LowModeDeflation_Approx::mult_left(const Field& f)const{
  Field w,v;
  SubSpace::project(w,f,*subVecs_);
  SubSpace::project(v,mult_little_inv(w),*subVecs_);  
  w = f;
  w -= Dw_->mult(v);
  return w;
}

const Field Dirac_LowModeDeflation_Approx::mult_right(const Field& f)const{
  Field w,v;
  SubSpace::project(w,Dw_->mult(f),*subVecs_);
  SubSpace::project(v,mult_little_inv(w),*subVecs_);  
  w = f;
  w -= v;

  rightPrec(v,w);
  return v;
}

void Dirac_LowModeDeflation_Approx::rightPrec(Field& w,const Field& f)const{
  w=f; 
}


void Dirac_LowModeDeflation_Approx::calcLittleOp(){
  for(int i=0; i<Nsub_; ++i){
    for(int j=0; j<Nsub_; ++j){
      Dsub_[2*(Nsub_*i +j)  ] = (*subVecs_)[i]*Dw_->mult((*subVecs_)[j]);
      Dsub_[2*(Nsub_*i +j)+1] = (*subVecs_)[i].im_prod(Dw_->mult((*subVecs_)[j]));
    }
  }
  ComplexMatrix::invert(Dsub_inv_,Dsub_);
  ComplexMatrix::hermite(Dsub_dag_,Dsub_);
  ComplexMatrix::hermite(Dsub_dag_inv_,Dsub_inv_);
}

void Dirac_LowModeDeflation_Approx::
lowModeProj(Field& w,const Field& f,const valarray<double>& Dkl)const{
  slice re(0,f.size()/2,2);
  slice im(1,f.size()/2,2);

  w = 0.0;
  for(int i=0; i<Nsub_; ++i){
    double DPfr = 0.0;
    double DPfi = 0.0;

    for(int j=0; j<Nsub_; ++j){
      double vfr = (*subVecs_)[j]*f;
      double vfi = (*subVecs_)[j].im_prod(f);

      DPfr += Dkl[2*(Nsub_*i+j)  ]*vfr -Dkl[2*(Nsub_*i+j)+1]*vfi;
      DPfi += Dkl[2*(Nsub_*i+j)+1]*vfr +Dkl[2*(Nsub_*i+j)  ]*vfi;
    }      
    w.add(re,DPfr*((*subVecs_)[i][re] -DPfi*(*subVecs_)[i][im]));
    w.add(im,DPfi*((*subVecs_)[i][re] +DPfr*(*subVecs_)[i][im]));
  }
}
