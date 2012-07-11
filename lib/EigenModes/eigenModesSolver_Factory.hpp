/*! @file eigenModsSolver_Factory.hpp
 *  @brief factory of EigenmodesSolver
 */
#ifndef EIGENMODESSOLVER_FACTORY_INCLUDED
#define EIGENMODESSOLVER_FACTORY_INCLUDED

#include "include/common_fields.hpp"
#include "eigenModesSolver.hpp"
#include "eigenModesSolver_IRL.hpp"
#include "eigenSorter_Factory.hpp"
#include "foprHermFactory.hpp"
#include "include/pugi_interface.h"

class EigenSolverFactory{
public:
  virtual EigenModesSolver* getEigenSolver(const Fopr_Herm*)const =0;
  virtual ~EigenSolverFactory(){}
};

class EigenSolverFactory_IRL: public EigenSolverFactory{
  const XML::node eigslv_node_;
public:
  EigenSolverFactory_IRL(XML::node node):eigslv_node_(node){}
  
  EigenModesSolver* getEigenSolver(const Fopr_Herm* opr)const{
    return new EigenModesSolver_IRL(opr,eigslv_node_);
  }
};

namespace EigenModes{
  EigenSolverFactory* createEigenSolverFactory(XML::node);
}
#endif
