/*!
 * @file solver_CG.hpp
 * @brief Declaration of Solver_CG class
 * Time-stamp: <2014-07-02 21:31:27 noaki>
 */
#ifndef SOLVER_CG_INCLUDED
#define SOLVER_CG_INCLUDED

#include <typeinfo>
#include "Fopr/fopr.h"
#include "solver.hpp"
#include "solver_CG_params.hpp"

/*!
 * @brief Solves \f$Dx = b\f$ using 
 * <a href="http://en.wikipedia.org/wiki/Conjugate_gradient_method">Conjugate Gradient method</a>
 * An hermitian operator is assumed 
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


#endif    
