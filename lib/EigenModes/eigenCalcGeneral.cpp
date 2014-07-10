/*!@file eigenCalcGeneral.cpp
 * @brief generalizes the eigenmodes calculation of fermion operators
 */
#include "PugiXML/xmlUtilities.hpp"
#include "eigenCalcGeneral.hpp"
#include "eigenSorter_Factory.hpp"
#include "chebyshevAccelFuncFactory.hpp"
#include "Fields/field_expressions.hpp"
#include "Fopr/foprHermFactoryCreator.hpp"
#include "inputConfig.hpp"
#include "field.h"
#include <cassert>

using namespace std;

EigenCalcGeneral::EigenCalcGeneral(const XML::node& node){
  XML::node opNode = node;
  XML::descend(opNode,"HermitianOperator", MANDATORY);
  opOrigFptr_.reset(HermiteOp::createFoprHermFactory(opNode));

  XML::node setupNode = node;
  XML::descend(setupNode,"Setup", MANDATORY);
  createEigenSorterFactory(setupNode);  /// initializing esortFptr_
  
  XML::node acNode = setupNode;
  XML::descend(acNode,"Acceleration");
  const char* ac_name = acNode.attribute("name").value();
  if(!strcmp(ac_name,"Chebyshev"))
    opAccelFptr_.reset(new ChebyshevAccelFuncFactory(setupNode));

  XML::node eslvNode = node;
  XML::descend(eslvNode,"EigenModesSolver", MANDATORY);
  eslvFptr_.reset(EigenSolver::createEigenSolverFactory(eslvNode));  
}

EigenSorterFactory* EigenCalcGeneral::createEigenSorterFactory(const XML::node& node){
  const char* st_name = node.attribute("sorting").value();
  double thrs;
  XML::read(node,"threshold",thrs);

  XML::node ac_node = node;
  XML::descend(ac_node,"Acceleration");
  const char* ac_name = ac_node.attribute("name").value();

  if(!strcmp(ac_name,"None")){
    if(     !strcmp(st_name,"Lowest") )  esortFptr_.reset(new EigenSorterFactory_low( thrs));
    else if(!strcmp(st_name,"Highest"))  esortFptr_.reset(new EigenSorterFactory_high(thrs));
    else XML::stopMsg(ac_node,st_name);
  }else if(!strcmp(ac_name,"Chebyshev")){esortFptr_.reset(new EigenSorterFactory_high(thrs));
  }else{ XML::stopMsg(ac_node,ac_name);}
}

void EigenCalcGeneral::do_calc(InputConfig& input){

  auto_ptr<Fopr_Herm> opOrigPtr(opOrigFptr_->getFoprHerm(input));

  auto_ptr<Fopr_Herm> opAccelPtr;
  if(opAccelFptr_.get()) opAccelPtr.reset(opAccelFptr_->getOp(opOrigPtr.get()));
  else                   opAccelPtr.reset(opOrigPtr.get());

  auto_ptr<EigenSorter>      sorterPtr(esortFptr_->getEigenSorter(opAccelPtr.get()));
  auto_ptr<EigenModesSolver> emslvPtr(eslvFptr_->getEigenSolver(opAccelPtr.get(),
								sorterPtr.get()));
  emslvPtr->calc(evals_,evecs_,Neig_); 
  if(Neig_> 0){
    CCIO::cout<<"Calculation successfully finished. Eigenvalues are:\n";
    get_eval(opOrigPtr.get()); // eigenvalues of oopr 
    
  }else if(Neig_== 0){
    CCIO::cout<<"NO desired eigenmode.\n";
    throw "Calculation finished with error";

  }else if(Neig_< 0){/*!<@brief it means emslvPtr->calc() ended abnormally. */
    Neig_*= -1;           
    CCIO::cout<<Neig_<<" eigenvalues are obtained:\n";

    get_eval(opOrigPtr.get()); // eigenvalues of oopr 
    throw "Calculation finished with error\n";
  }
}

void EigenCalcGeneral::get_eval(const Fopr_Herm* opr){
  using namespace FieldExpression;
  assert(Neig_>= 0);

  CCIO::cout<< setiosflags(ios_base::scientific);
  for(int i=0; i<Neig_; ++i){
    Field Av = opr->mult(evecs_[i]);
    evals_[i] = evecs_[i]*Av;    
    Av -= evals_[i]*evecs_[i]; 
    double res = Av.norm(); // residual 
    double vndif = evecs_[i].norm() -1.0;

    CCIO::cout<<" ["<<setw( 3)<<setiosflags(ios_base::right)<< i<<"] ";
    CCIO::cout<<      setw(25)<<setiosflags(ios_base::left )<< evals_[i];
    CCIO::cout<<"  "<<setw(25)<<setiosflags(ios_base::right)<< res;
    CCIO::cout<<"  "<<setw(25)<<setiosflags(ios_base::right)<< vndif <<"\n";
  }
  CCIO::cout<< resetiosflags(ios_base::scientific);
}

void EigenCalcGeneral::output_txt(const string& output)const{
  if(Communicator::instance()->primaryNode()){
    ofstream writer(output.c_str()); 
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
  }
}

