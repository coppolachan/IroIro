//---------------------------------------------------------------------------
// fopr.h
//---------------------------------------------------------------------------
#ifndef FOPR_INCLUDED
#define FOPR_INCLUDED

#ifndef FIELD_INCLUDED
#include "field.h"
#endif

#include "Dirac_ops/dirac.h"
#include "Solver/solver.h"

class Fopr{
public:
  virtual ~Fopr(){}
  virtual const Field mult(const Field&) const =0;
  virtual const Field mult_dag(const Field&) const =0;
  virtual size_t fsize() const =0;
};

class Fopr_D:public Fopr{
private:
  const Dirac* D_;
public:
  Fopr_D(const Dirac* D):D_(D){}

  const Field mult(    const Field& f)const{return D_->mult(f);}
  const Field mult_dag(const Field& f)const{return D_->mult_dag(f);}
  size_t fsize()const {return D_->fsize();}
};

class Fopr_Ddag:public Fopr{
private:
  const Dirac* D_;
public:
  Fopr_Ddag(const Dirac* D):D_(D){}

  const Field mult(const Field& f) const {return D_->mult_dag(f);}
  const Field mult_dag(const Field& f) const {return D_->mult(f);}
  size_t fsize()const {return D_->fsize();}
};

class Fopr_H : public Fopr {
private:
  const DiracWilsonLike* D_;
public:
  Fopr_H(const DiracWilsonLike* D):D_(D){}
  
  const Field mult(const Field& f) const{ return D_->gamma5(D_->mult(f));}
  const Field mult_dag(const Field& f) const{return mult(f);}
  const Field gamma5(const Field&f) const { return D_->gamma5(f); }

  size_t fsize()const {return D_->fsize();}
};


////////////// Hermitian operators

class Fopr_Herm : public Fopr {
public:
  const Field mult_dag(const Field& f) const{return mult(f);}
};

class Fopr_DdagD : public Fopr_Herm {
private:
  const Dirac* D_;
public:
  Fopr_DdagD(const Dirac* D):D_(D){}
  
  const Field mult(const Field& f) const{return D_->mult_dag(D_->mult(f));}
  size_t fsize()const {return D_->fsize();}
};

class Fopr_DDdag: public Fopr_Herm {
private:
  const Dirac* D_;
public:
  Fopr_DDdag(const Dirac* D):D_(D){}

  const Field mult(const Field& f) const{return D_->mult(D_->mult_dag(f));}
  size_t fsize()const {return D_->fsize();}
};

#endif

