/*!@file eigenCalcGeneral.cpp
 * @brief generalizes the eigenmodes calculation of fermion operators
 */
#include "eigenCalcGeneral.hpp"
#include "eigenSorter_Factory.hpp"
#include "chebyshevAccelFunc.hpp"
#include "IO/fields_io.hpp"
#include "Fields/field_expressions.hpp"
#include "inputConfig.hpp"
#include "field.h"

#include <cassert>

using namespace std;

EigenCalcGeneral::EigenCalcGeneral(const XML::node& node){

  XML::node targetOpNode = node;
  XML::descend(targetOpNode,"HermitianOperator");
  opOrigFptr_.reset(createFoprHermFactory(targetOpNode));

  XML::node setupNode = node;
  XML::descend(setupNode,"Setup");
  esortFptr_.reset(createEigenSorterFactory(setupNode));
  XML::descend(setupNode,"Acceleration");
  opAccelFptr_.reset(createAccelOpFunc(setupNode));

  XML::node eslvNode = node;
  XML::descend(eslvNode,"EigenModesSolver");
  eslvFptr_.reset(EigenSolver::createEigenSolverFactory(eslvNode));  
}

FoprHermFunc* EigenCalcGeneral::createAccelOpFunc(const XML::node& node)const{
  XML::node ac_node = node;
  XML::descend(ac_node,"Acceleration");
  const char* ac_name = ac_node.attribute("name").value();
  
  if(!strcmp(ac_name,"None"))      return new FoprAsIsFunc();
  if(!strcmp(ac_name,"Chebyshev")) return new ChebyshevAccelFunc(node);
  CCIO::cout<<ac_name<<" is not compatible with current implementation.\n";
  abort();
}

FoprHermFactory* EigenCalcGeneral::createFoprHermFactory(const XML::node& node)const{
  const char* hf_name = node.attribute("name").value();
  
  if(!strcmp(hf_name,"g5D"))   return new FoprHermFactory_H(node);
  if(!strcmp(hf_name,"DdagD")) return new FoprHermFactory_DdagD(node);
  if(!strcmp(hf_name,"DDdag")) return new FoprHermFactory_DDdag(node);
  if(!strcmp(hf_name,"Laplacian")) return new FoprHermFactory_Scalar(node);
  CCIO::cout<<hf_name<<" is not compatible with current implementation.\n";
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

void EigenCalcGeneral::do_calc(InputConfig& input){

  const auto_ptr<Fopr_Herm> opOrigPtr(opOrigFptr_->getFoprHerm(input));
  const auto_ptr<Fopr_Herm> opAccelPtr(opAccelFptr_->getFoprHerm(opOrigPtr.get()));
  const auto_ptr<EigenSorter> sorterPtr(esortFptr_->getEigenSorter(opAccelPtr.get()));
  const auto_ptr<EigenModesSolver> emslvPtr(eslvFptr_->getEigenSolver(opAccelPtr.get(),
								      sorterPtr.get()));

  emslvPtr->calc(evals_,evecs_,Neig_); /*!< @brief solving eigenproblem of aopr
					 eval and evec are resized inside */
  if(Neig_> 0){
    CCIO::cout<<"Calculation successfully finished. Eigenvalues are:\n";
    get_eval(opOrigPtr.get()); // eigenvalues of oopr 
    
  }else if(Neig_== 0){
    CCIO::cout<<"NO converged eigenmode.\n";
    throw "Calculation did not successfully finished.";

  }else if(Neig_< 0){/*!<@brief it means emslvPtr->calc() ended abnormally. */
    Neig_*= -1;           
    CCIO::cout<<Neig_<<" eigenvalues are obtained:\n";

    get_eval(opOrigPtr.get()); // eigenvalues of oopr 
    throw "Calculation abnormally finished.\n";
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
  std::string output_evals = output + "_evals.txt";
  ofstream writer(output_evals.c_str());
  for(int i=0; i<Neig_; ++i){
    CCIO::SaveOnDisk<Format::Format_F>(evecs_[i],output.c_str(),true);
    writer<< setw(2) <<setiosflags(ios_base::right)<< i;
    writer<< setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	  << evals_[i]<<endl;
  }
}
