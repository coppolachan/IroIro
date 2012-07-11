/* @file test_EigenModesSolver.cpp
 * @brief implementation of the Test_EigenModesSolver class
 */
#include "test_EigenModesSolver.hpp"
#include "EigenModes/eigenModesSolver_Factory.hpp"
#include "Fields/field_expressions.hpp"
#include "include/field.h"
#include <memory>

using namespace std;
using namespace FieldExpression;

int Test_EigenModesSolver::run(){

  //// EigenModesSolver Factory ////
  XML::descend(node_,"EigenModesSolver");
  auto_ptr<EigenSolverFactory> eslv_fact(EigenModes::createEigenSolverFactory(node_));

  //// HermitianOperator Factory ////
  XML::next_sibling(node_,"HermiteOpr");
  auto_ptr<FoprHermFactory> opfact(FoprHerm::createFoprHermFactory(node_));
  XML::next_sibling(node_,"AccelOpr");
  auto_ptr<FoprHermFactory> aofact(FoprHerm::createFoprHermFactory(node_));

  //// WilsonLikeDiracOperator Factory ////
  XML::next_sibling(node_,"WilsonLikeDirac");
  auto_ptr<DiracWilsonLikeOperatorFactory> 
    dfact(DiracOperators::createGeneralDiracWilsonLikeOperatorFactory(node_));

  //// creating Objects ////  
  auto_ptr<DiracWilsonLike>  dirac(dfact->getDiracOperatorWL(&(conf_.data)));
  auto_ptr<Fopr_Herm>        opr(aofact->getFoprHerm(dirac.get()));
  auto_ptr<EigenModesSolver> emslv(eslv_fact->getEigenSolver(opr.get()));

  //// actual calculation ////
  vector<double> eval;  
  vector<Field> evec;
  int Neig;
  emslv->calc(eval,evec,Neig);  /*!< @brief eval and evec are resized inside */

  //// true_eigenvalue & residual ////  
  opr = auto_ptr<Fopr_Herm>(opfact->getFoprHerm(dirac.get()));

  CCIO::cout<< setiosflags(ios_base::scientific);
  for(int i=0; i<=Neig; ++i){
    Field Av = opr->mult(evec[i]);
    double vv = evec[i]*evec[i];
    eval[i] = evec[i]*Av;    
    eval[i] /= vv;    
    Av -= eval[i]*evec[i]; 
    double res = Av*Av; // residual 
    
    CCIO::cout<<" ["<<setw( 3)<<setiosflags(ios_base::right)<< i<<"] ";
    CCIO::cout<<      setw(25)<<setiosflags(ios_base::left )<< eval[i];
    CCIO::cout<<"  "<<setw(25)<<setiosflags(ios_base::right)<< res;
    CCIO::cout<<"  "<<setw(25)<<setiosflags(ios_base::right)<< vv-1.0 <<endl;
  }
  CCIO::cout<< resetiosflags(ios_base::scientific);

  ///// output here if necessary /////
  
  return 0;
}
