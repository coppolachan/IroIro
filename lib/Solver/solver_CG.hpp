/*!
 * @file solver_CG.hpp
 * @brief Declaration of Solver_CG class
 * Time-stamp: <2013-07-05 18:41:12 noaki>
 */
#ifndef SOLVER_CG_INCLUDED
#define SOLVER_CG_INCLUDED

#include <typeinfo>

#include "include/pugi_interface.h"
#include "include/fopr.h"
#include "solver.hpp"

/*!
 * @brief Structure containing parameters for the Solver_CG class
 *
 */
struct Solver_CG_Prms{
  int MaxIter;/*!< Maximum number of iteration for the solver */
  double GoalPrecision; /*!< Threshold for the final residual */
  
  Solver_CG_Prms(const XML::node node){
    XML::read(node, "MaxIter", MaxIter);
    XML::read(node, "Precision", GoalPrecision);
  }
  
  Solver_CG_Prms(const double prec_, const double MaxIter_){
    MaxIter       = MaxIter_;
    GoalPrecision = prec_;
  }
};

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

#ifdef IBM_BGQ_WILSON
#include "Dirac_ops/dirac_DomainWall_EvenOdd.hpp"
/*!
 * @brief Solves \f$Dx = b\f$ using 
 * <a href="http://en.wikipedia.org/wiki/Conjugate_gradient_method">Conjugate Gradient method</a>
 * Domain Wall fermions optimized operator for BGQ Wrapper class
 */
class Solver_CG_DWF_Optimized: public Solver{
private:
  const Dirac_optimalDomainWall_EvenOdd* opr_;
  const Solver_CG_Prms Params;/*!< @brief Inputs container */
  const int nodeid_;

public:
  Solver_CG_DWF_Optimized(const double prec,const int MaxIterations,
			  const Dirac_optimalDomainWall_EvenOdd* DWFopr)
    :opr_(DWFopr),
     nodeid_(Communicator::instance()->nodeid()),
     Params(Solver_CG_Prms(prec, MaxIterations)){}

  Solver_CG_DWF_Optimized(const XML::node Solver_node,
			  const Dirac_optimalDomainWall_EvenOdd* DWFopr)
    :opr_(DWFopr),
     Params(Solver_CG_Prms(Solver_node)),
     nodeid_(Communicator::instance()->nodeid()){}

  ~Solver_CG_DWF_Optimized(){}

  SolverOutput solve(Field& solution, const Field& source) const;

  bool check_DdagD() const {
    //empty in this case
    return 1;
  }
};

#endif

#endif    
