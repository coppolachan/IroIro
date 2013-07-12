/*!
 * @file solver_CG_params.hpp
 *
 * @brief Declaration of Solver_CG_Prms structure
 *
 * Time-stamp: <2013-07-04 16:09:52 cossu>
 */
#ifndef SOLVER_CG_PARAMS_INCLUDED
#define SOLVER_CG_PARAMS_INCLUDED

#include "include/pugi_interface.h"
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

#endif
