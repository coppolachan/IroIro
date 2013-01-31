/*!@file eigenCalcGeneral.cpp
 * @brief generalizes the eigenmodes calculation of fermion operators
 */
#include "eigenCalcGeneral.hpp"
#include "foprHermFactory_ChebyshevDdagDLin.hpp"
#include "eigenSorter_Factory.hpp"
#include "include/field.h"
#include "Communicator/fields_io.hpp"
#include "Fields/field_expressions.hpp"
#include <cassert>

using namespace std;

EigenCalcGeneral::EigenCalcGeneral(const XML::node& node){

  XML::node diracNode = node;
  XML::descend(diracNode,"WilsonLikeDirac");
  diracFptr_.reset(DiracOperators::createDiracWilsonLikeOperatorFactory(diracNode));

  XML::node oprNode = node;
  XML::descend(oprNode,"Operator");
  opOrigFptr_.reset(createFoprHermFactory(oprNode));
  opAccelFptr_.reset(createAccelOpFactory(oprNode));
  esortFptr_.reset(createEigenSorterFactory(oprNode));
  
  XML::node eslvNode = node;
  XML::descend(eslvNode,"EigenModesSolver");
  eslvFptr_.reset(EigenModes::createEigenSolverFactory(eslvNode));  
}

void EigenCalcGeneral::do_calc(Field* const conf){

  const auto_ptr<DiracWilsonLike> diracPtr(diracFptr_->getDiracOperatorWL(conf));
  const auto_ptr<Fopr_Herm>       aoprPtr(opAccelFptr_->getFoprHerm(diracPtr.get()));
  const auto_ptr<EigenSorter>     sorterPtr(esortFptr_->getEigenSorter(aoprPtr.get()));
  const auto_ptr<EigenModesSolver> emslvPtr(eslvFptr_->getEigenSolver(aoprPtr.get(),sorterPtr.get()));
  emslvPtr->calc(evals_,evecs_,Neig_); /*!< @brief solving eigenproblem of aopr
					 eval and evec are resized inside */

  const auto_ptr<Fopr_Herm> ooprPtr(opOrigFptr_->getFoprHerm(diracPtr.get()));

  if(Neig_> 0){
    CCIO::cout<<"Calcuration successfully finished. Eigenvalues are:\n";
    get_eval(ooprPtr.get()); // eigenvalues of oopr 
    
  }else if(Neig_== 0){
    CCIO::cout<<"NO converged eigenmode.\n";
    throw "Calcuration did not successfully finished.";

  }else if(Neig_< 0){/*!<@brief it means emslvPtr->calc() ended abnormally. */
    Neig_*= -1;           
    CCIO::cout<<Neig_<<" eigenvalues are obtained:\n";

    get_eval(ooprPtr.get()); // eigenvalues of oopr 
    throw "Calcuration abnormally finished.\n";
  }
}

void EigenCalcGeneral::get_eval(const Fopr_Herm* opr){
  using namespace FieldExpression;
  assert(Neig_>= 0);

  CCIO::cout<< setiosflags(ios_base::scientific);
  for(int i=0; i<Neig_; ++i){
    Field Av = opr->mult(evecs_[i]);
    double vv = evecs_[i]*evecs_[i];

    evals_[i] = evecs_[i]*Av;    
    evals_[i] /= vv;    
    Av -= evals_[i]*evecs_[i]; 
    double res = Av*Av; // residual 
    
    CCIO::cout<<" ["<<setw( 3)<<setiosflags(ios_base::right)<< i<<"] ";
    CCIO::cout<<      setw(25)<<setiosflags(ios_base::left )<< evals_[i];
    CCIO::cout<<"  "<<setw(25)<<setiosflags(ios_base::right)<< res;
    CCIO::cout<<"  "<<setw(25)<<setiosflags(ios_base::right)<< vv-1.0 <<endl;
  }

  CCIO::cout<< resetiosflags(ios_base::scientific);
}

void EigenCalcGeneral::output_txt(const string& output)const{
  if(Communicator::instance()->primaryNode()){
    ofstream writer(output.c_str()); 
    CCIO::cout<<"starting output\n";
    for(int i=0; i<Neig_; ++i){
      writer<< setw(2) <<setiosflags(ios_base::right)<< i;
      writer<< setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << evals_[i]<<endl;

      for(int k=0; k<evecs_[i].size()/2; ++k){
	writer<< setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	      << evecs_[i][2*k];
	writer<< setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	      << evecs_[i][2*k+1]
	      << endl;
      }
    }
    CCIO::cout<<"output finished\n";
  }
}

void EigenCalcGeneral::output_bin(const string& output)const{
  for(int i=0; i<Neig_; ++i)
    CCIO::SaveOnDisk<Format::Format_F>(evecs_[i],output.c_str(),true);
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

