/*!
 * @file solver_CG.hpp
 *
 * @brief Declaration of Solver_CG class
 *
 * Time-stamp: <2013-07-04 16:09:36 cossu>
 */
#ifndef SOLVER_CG_INCLUDED
#define SOLVER_CG_INCLUDED

#include <typeinfo>

#include "Dirac_ops/dirac_Preconditioners.hpp"
#include "include/fopr.h"
#include "solver.hpp"
#include "solver_CG_params.hpp"

/*!
 * @brief Solves \f$Dx = b\f$ using 
 * <a href="http://en.wikipedia.org/wiki/Conjugate_gradient_method">Conjugate Gradient method</a>
 *
 * An hermitian operator is assumed 
 * 
 * WARNING: Never use with (non trivially) preconditioned operators
 */
class Solver_CG: public Solver{
private:
  const Fopr_Herm* opr_;/*!< @brief Hermitian input operator */
  const Solver_CG_Prms Params;/*!< @brief Inputs container */
  const int nodeid_;

  void solve_step(Field&, Field&, Field&, double&, double&) const;

public:
  Solver_CG(const double prec,const int MaxIterations,const Fopr_Herm* fopr)
    :opr_(fopr),
     nodeid_(Communicator::instance()->nodeid()),
     Params(Solver_CG_Prms(prec, MaxIterations)){}

  Solver_CG(const XML::node Solver_node,const Fopr_Herm* fopr)
    :opr_(fopr),
     Params(Solver_CG_Prms(Solver_node)),
     nodeid_(Communicator::instance()->nodeid()){}

  ~Solver_CG(){}

  SolverOutput solve(Field& solution, const Field& source) const;

  bool check_DdagD() const {
    return (typeid(*opr_) == typeid(Fopr_DdagD));
  }
};

/*!
 * @brief Solves \f$Dx = b\f$ using 
 * <a href="http://en.wikipedia.org/wiki/Conjugate_gradient_method">Conjugate Gradient method</a>
 * 
 * An hermitian preconditioned operator is assumed 
 *
 */
class Solver_CG_Precondition : public Solver {
private:
  const Fopr_Herm_Precondition* opr_;/*!< @brief Hermitian Preconditioned input operator */
  const Solver_CG_Prms Params;/*!< @brief Inputs container */
  const int nodeid_;

  void solve_step(Field&, Field&, Field&, double&) const;

public:
  Solver_CG_Precondition(const double prec, 
			 const int MaxIterations,
			 const Fopr_Herm_Precondition* fopr)
    :opr_(fopr),
     nodeid_(Communicator::instance()->nodeid()),
     Params(Solver_CG_Prms(prec, MaxIterations)){}
  
  Solver_CG_Precondition(const XML::node Solver_node,
			 const Fopr_Herm_Precondition* fopr)
    :opr_(fopr),
     Params(Solver_CG_Prms(Solver_node)),
     nodeid_(Communicator::instance()->nodeid()){}
  
  ~Solver_CG_Precondition(){}
  
  SolverOutput solve(Field& solution, const Field& source) const;
  
  bool check_DdagD() const  {
    return  (typeid(*opr_) == typeid(Fopr_DdagD_Precondition));
  }
};


class Solver_CG_Precondition_New : public Solver {
private:
  const Fopr_Herm_Precondition* opr_;/*!< @brief Hermitian Preconditioned input operator */
  const Solver_CG_Prms Params;/*!< @brief Inputs container */
  const int nodeid_;
  const Preconditioner& prec;


  void solve_step(Field&, Field&, Field&, double&) const;

  // notes: take the preconditioner directly from the fopr
  // give public access to preconditioner 
  // function inf fopr_herm_precondition to access the preconditioner

  // set mode to solve D or DdagD (the source preconditioning is affected)
  // at creation time (a general solver should't know about this
  // so no functions beside solve() must be created)

public:
  Solver_CG_Precondition_New(const double prec, 
			     const int MaxIterations,
			     const Fopr_Herm_Precondition* fopr)
    :opr_(fopr),
     nodeid_(Communicator::instance()->nodeid()),
     Params(Solver_CG_Prms(prec, MaxIterations)),
     prec(*(new DefaultPreconditioner())){}

  Solver_CG_Precondition_New(const double prec, 
			     const int MaxIterations,
			     const Fopr_Herm_Precondition* fopr,
			     const Preconditioner& SetPrec)
    :opr_(fopr),
     nodeid_(Communicator::instance()->nodeid()),
     Params(Solver_CG_Prms(prec, MaxIterations)),
     prec(SetPrec){}
  
  Solver_CG_Precondition_New(const XML::node Solver_node,
			     const Fopr_Herm_Precondition* fopr)
    :opr_(fopr),
     Params(Solver_CG_Prms(Solver_node)),
     nodeid_(Communicator::instance()->nodeid()),
     prec(*(new DefaultPreconditioner())){}

  Solver_CG_Precondition_New(const XML::node Solver_node,
			     const Fopr_Herm_Precondition* fopr,
			     const Preconditioner& SetPrec)
    :opr_(fopr),
     Params(Solver_CG_Prms(Solver_node)),
     nodeid_(Communicator::instance()->nodeid()),
     prec(SetPrec){}
  
  ~Solver_CG_Precondition_New(){}
  
  SolverOutput solve(Field& solution, const Field& source) const;
  
  bool check_DdagD() const  {
    return  (typeid(*opr_) == typeid(Fopr_DdagD));
  }
};



#endif    
