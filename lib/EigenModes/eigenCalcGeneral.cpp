/*!@file eigenCalcGeneral.cpp
 * @brief generalizes the eigenmodes calculation of fermion operators
 */
#include "eigenCalcGeneral.hpp"
#include "foprHermFactory_ChebyshevDdagDLin.hpp"
#include "eigenSorter_Factory.hpp"
#include "include/field.h"
#include "Fields/field_expressions.hpp"
#include <cassert>

using namespace std;

EigenCalcGeneral::EigenCalcGeneral(const XML::node& node){

  XML::node diracNode = node;
  XML::descend(diracNode,"WilsonLikeDirac");
  diracFact_.reset(DiracOperators::createDiracWilsonLikeOperatorFactory(diracNode));

  XML::node oprNode = node;
  XML::descend(oprNode,"Operator");
  opOrigFact_.reset(createFoprHermFactory(oprNode));
  opAccelFact_.reset(createAccelOpFactory(oprNode));
  esortFact_.reset(createEigenSorterFactory(oprNode));
  
  XML::node eslvNode = node;
  XML::descend(eslvNode,"EigenModesSolver");
  eslvFact_.reset(EigenModes::createEigenSolverFactory(eslvNode));  
}

void EigenCalcGeneral::do_calc(Field* const conf){
  
  auto_ptr<DiracWilsonLike> dirac(diracFact_->getDiracOperatorWL(conf));
  auto_ptr<Fopr_Herm>       aopr(opAccelFact_->getFoprHerm(dirac.get()));
  auto_ptr<EigenSorter>     sorter(esortFact_->getEigenSorter(aopr.get()));
  auto_ptr<EigenModesSolver> emslv(eslvFact_->getEigenSolver(aopr.get(),sorter.get()));
  
  emslv->calc(eval_,evec_,Neig_); /*!< @brief solving eigenproblem of aopr
				    eval and evec are resized inside */

  auto_ptr<Fopr_Herm> oopr(opOrigFact_->getFoprHerm(dirac.get()));
  get_eval(oopr.get()); // eigenvalues of oopr 
}

void EigenCalcGeneral::get_eval(const Fopr_Herm* opr){
  using namespace FieldExpression;
  assert(Neig_>0);

  CCIO::cout<< setiosflags(ios_base::scientific);
  for(int i=0; i<=Neig_; ++i){
    Field Av = opr->mult(evec_[i]);
    double vv = evec_[i]*evec_[i];
    eval_[i] = evec_[i]*Av;    
    eval_[i] /= vv;    
    Av -= eval_[i]*evec_[i]; 
    double res = Av*Av; // residual 
    
    CCIO::cout<<" ["<<setw( 3)<<setiosflags(ios_base::right)<< i<<"] ";
    CCIO::cout<<      setw(25)<<setiosflags(ios_base::left )<< eval_[i];
    CCIO::cout<<"  "<<setw(25)<<setiosflags(ios_base::right)<< res;
    CCIO::cout<<"  "<<setw(25)<<setiosflags(ios_base::right)<< vv-1.0 <<endl;
  }
  CCIO::cout<< resetiosflags(ios_base::scientific);
}

void EigenCalcGeneral::output(ofstream& writer){

  if(Communicator::instance()->primaryNode()){
    CCIO::cout<<"starting output\n";
    for(int i=0; i<Neig_; ++i){
      writer<< setw(2) <<setiosflags(ios_base::right)<< i;
      writer<< setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << eval_[i]<<endl;

      for(int k=0; k<evec_[i].size()/2; ++k){
	writer<< setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	      << evec_[i][2*k];
	writer<< setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	      << evec_[i][2*k+1]
	      << endl;
      }
    }
    CCIO::cout<<"output finished\n";
  }
}

FoprHermFactory* EigenCalcGeneral::createAccelOpFactory(const XML::node& node)const{
  XML::node ac_node = node;
  XML::descend(ac_node,"Acceleration");
  const char* ac_name = ac_node.attribute("name").value();
  
  if(!strcmp(ac_name,"None"))    
    return createFoprHermFactory(node);
  if(!strcmp(ac_name,"Chebyshev")) 
    return new FoprHermFactory_ChebyshevDdagDLin(node);
  CCIO::cout<<ac_name<<" is not compatible with current implementation.\n";
  abort();
}

FoprHermFactory* EigenCalcGeneral::createFoprHermFactory(const XML::node& node)const{
  const char* name = node.attribute("name").value();
  
  if(!strcmp(name,"Hx"))    return new FoprHermFactory_H;
  if(!strcmp(name,"DdagD")) return new FoprHermFactory_DdagD;
  if(!strcmp(name,"DDdag")) return new FoprHermFactory_DDdag;
  CCIO::cout<<name<<" is not compatible with current implementation.\n";
  abort();
}

EigenSorterFactory* EigenCalcGeneral::createEigenSorterFactory(const XML::node& node)const{

  const char* st_name = node.attribute("sorting").value();
  double thrs;
  XML::read(node,"threshold",thrs);

  XML::node ac_node = node;
  XML::descend(ac_node,"Acceleration");
  const char* ac_name = ac_node.attribute("name").value();

  if(!strcmp(ac_name,"None")){
    if(     !strcmp(st_name,"Lowest") ) return new EigenSorterFactory_low(thrs);
    else if(!strcmp(st_name,"Highest")) return new EigenSorterFactory_high(thrs);
    else CCIO::cout<<st_name<<" is not compatible with current implementation.\n";

  }else if(!strcmp(ac_name,"Chebyshev")){
    return new EigenSorterFactory_high(thrs);

  }else CCIO::cout<<ac_name<<" is not compatible with current implementation.\n";
  abort();
}

