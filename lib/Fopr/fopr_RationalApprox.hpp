#ifndef FOPR_RATIONALAPPROX_INCLUDED
#define FOPR_RATIONALAPPROX_INCLUDED

#include "fopr.h"
#include "lib/Tools/RationalApprox/rationalapprox.hpp"
#include "lib/Fields/field_expressions.hpp"
#include "lib/Solver/multiShiftSolver.hpp"

class Fopr_RationalApprox :public Fopr_Herm{
  const Fopr_Herm* Op_;
  const MultiShiftSolver* mslv_;

  double r0_;
  std::vector<double> res_;
  std::vector<double> pole_;
public:
  Fopr_RationalApprox(double r0,
		      const std::vector<double>& res, 
		      const std::vector<double>& pole, 
		      const Fopr_Herm* Op,
		      const MultiShiftSolver* mslv)
    :r0_(r0),res_(res),pole_(pole),Op_(Op),mslv_(mslv){}
  
  const Field mult(const Field& f)const{
    using namespace FieldExpression;

    std::vector<Field> sol(res_.size());
    double r;
    int Nconv;

    try{
      SolverOutput so = mslv_->solve(sol,f,pole_,r,Nconv);
    }catch(const char* error){
      CCIO::cerr<<error<<std::endl;
      return EXIT_FAILURE;
    }

    Field w(f);
    w *= r0_;
    for(int i=0; i<res_.size(); ++i) w += res_[i]*sol[i];
    return w;
  }

  double func(double x)const{ 
    double w = r0_;
    for(int i=0; i<res_.size(); ++i)
      w += res_[i]/(Op_->func(x) +pole_[i]);
    return w;
  }

  size_t fsize()const{return Op_->fsize();}
};

#endif
