/*!
 * @file solver_CG.h
 *
 * @brief Declaration of Solver_CG class
 *
 */
#ifndef SOLVER_CG_INCLUDED
#define SOLVER_CG_INCLUDED

#include <typeinfo>
#include "include/pugi_interface.h"
#include "include/fopr.h"
#include "solver.h"

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
  };
  
  Solver_CG_Prms(const double prec_, const double MaxIter_){
    MaxIter       = MaxIter_;
    GoalPrecision = prec_;
  };

};

/*!
 * @brief Solves \f$Dx = b\f$ using 
 * <a href="http://en.wikipedia.org/wiki/Conjugate_gradient_method">Conjugate Gradient method</a>
 *
 *
 * An hermitian operator is assumed 
 */
class Solver_CG: public Solver{
protected:
  const Fopr_Herm* opr_;/*!< @brief Hermitian input operator */
  const Solver_CG_Prms Params;/*!< @brief Inputs container */
  const int nodeid_;


  void solve_step(Field&, Field&, Field&, double&) const;

public:
  Solver_CG(const double prec, 
	    const int MaxIterations,
	    const Fopr_Herm* fopr)
    :opr_(fopr),
     nodeid_(Communicator::instance()->nodeid()),
     Params(Solver_CG_Prms(prec, MaxIterations)){};

 Solver_CG(const XML::node Solver_node,
	   const Fopr_Herm* fopr)
    :opr_(fopr),
     Params(Solver_CG_Prms(Solver_node)),
     nodeid_(Communicator::instance()->nodeid()){};

  ~Solver_CG(){}

  void solve(Field& solution, const Field& source, 
	     double& diff, int& Nconv) const;

  bool check_DdagD() const  {
    return  (typeid(*opr_) == typeid(Fopr_DdagD));
  }

};

class Solver_CG_Preconditioned : public Solver_CG {
};


#endif    
