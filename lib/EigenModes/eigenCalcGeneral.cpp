/*!@file eigenCalcGeneral.cpp
 * @brief generalizes the eigenmodes calculation of fermion operators
 */
#include "PugiXML/xmlUtilities.hpp"
#include "eigenCalcGeneral.hpp"
#include "chebyshevAccelFuncFactory.hpp"
#include "Fields/field_expressions.hpp"
#include "Fopr/foprHermFactoryCreator.hpp"
#include "inputConfig.hpp"
#include "field.h"
#include <cassert>

using namespace std;

EigenCalcGeneral::EigenCalcGeneral(XML::node node):eslvNode_(node){
  XML::node opNode = node;
  XML::descend(opNode,"HermitianOperator", MANDATORY);
  opOrigFptr_.reset(HermiteOp::createFoprHermFactory(opNode));

  XML::node setupNode = node;
  XML::descend(setupNode,"Setup", MANDATORY);
  esortPtr_.reset(eigenSorterFactory(setupNode));

  XML::node acNode = setupNode;
  XML::descend(acNode,"Acceleration");
  const char* ac_name = acNode.attribute("name").value();
  if(strcmp(ac_name,"None"))
    opAccelFptr_.reset(new ChebyshevAccelFuncFactory(setupNode));

  XML::descend(eslvNode_,"EigenModesSolver", MANDATORY);
  const char* eslv_name = eslvNode_.attribute("name").value();
  assert(!strcmp(eslv_name,"ImplicitRestartedLanczos"));
}

EigenSorter* EigenCalcGeneral::eigenSorterFactory(const XML::node& node)const{
  const char* st_name = node.attribute("sorting").value();
  double thrs;
  XML::read(node,"threshold",thrs);

  XML::node ac_node = node;
  XML::descend(ac_node,"Acceleration");
  const char* ac_name = ac_node.attribute("name").value();

  if(!strcmp(ac_name,"None")){
    if(     !strcmp(st_name,"Lowest") )  return new EigenSorter_low( thrs);
    else if(!strcmp(st_name,"Highest"))  return new EigenSorter_high(thrs);
    else XML::stopMsg(ac_node,st_name);

  }else if(!strcmp(ac_name,"Chebyshev")){ return new EigenSorter_high(thrs);
  }else{ XML::stopMsg(ac_node,ac_name);}
}

EigenModesSolver* EigenCalcGeneral::eigSlvFactory(const Fopr_Herm* op,
						  const EigenSorter* esort){
  int Nk,Np,Niter,Nthrs,max_iter;
  double prec;
  XML::read(eslvNode_,"Nthreshold",Nthrs,MANDATORY);
  XML::read(eslvNode_,"Nk",Nk,MANDATORY);
  XML::read(eslvNode_,"Np",Np,MANDATORY);
  XML::read(eslvNode_,"precision",prec,MANDATORY);
  XML::read(eslvNode_,"max_iter",Niter,MANDATORY);

  if(opAccelFptr_.get()){
    CCIO::cout<<"Acceleration is ON\n";
    double thrs = op->func(esort->thrs());
    CCIO::cout<<"threshold = "<<thrs<<"\n";
  
    double scale = 0.5/prec*fabs(op->func(esort->thrs() +prec)
				 -op->func(esort->thrs() -prec));
    assert(scale > 1.0);
    CCIO::cout<<"rescaling factor = "<<scale<<"\n";
    prec*= scale;
    CCIO::cout<<"precision = "<<prec<<"\n";
  }else{
    CCIO::cout<<"Acceleration is OFF\n";
  }
  return new EigenModesSolver_IRL(op,esort,Nk,Np,prec,Niter,Nthrs);
}

void EigenCalcGeneral::do_calc(InputConfig& input){

  auto_ptr<Fopr_Herm> opPtr(opOrigFptr_->getFoprHerm(input));
  if(opAccelFptr_.get())
    opPtr.reset(opAccelFptr_->getOp(opOrigFptr_->getFoprHerm(input)));
  auto_ptr<EigenModesSolver> 
    emslvPtr(eigSlvFactory(opPtr.get(),esortPtr_.get()));

  emslvPtr->calc(evals_,evecs_,Neig_); 
  if(Neig_> 0){
    CCIO::cout<<"Calculation successfully finished. Eigenvalues are:\n";
    get_eval(opOrigFptr_->getFoprHerm(input));  // eigenvalues of oopr 

  }else if(Neig_== 0){
    CCIO::cout<<"NO desired eigenmode.\n";
    throw "Calculation finished with error";

  }else if(Neig_< 0){/*!<@brief it means emslvPtr->calc() ended abnormally. */
    Neig_*= -1;           
    CCIO::cout<<Neig_<<" eigenvalues are obtained:\n";
    get_eval(opOrigFptr_->getFoprHerm(input)); // eigenvalues of oopr 
    throw "Calculation finished with error\n";
  }
}

void EigenCalcGeneral::get_eval(const Fopr_Herm* opr){
  using namespace FieldExpression;
  assert(Neig_>= 0);

  CCIO::cout<<" index ";
  CCIO::cout<<      setw(20)<<setiosflags(ios_base::left)<<"eigenvalues";
  CCIO::cout<<"  "<<setw(22)<<setiosflags(ios_base::left)<<"residual";
  CCIO::cout<<"  "<<setw(24)<<setiosflags(ios_base::left)<<"|vec|-1.0\n";
  CCIO::cout<<"---------------------------------------------";
  CCIO::cout<<"---------------------------------------------\n";
  CCIO::cout<< setiosflags(ios_base::scientific);

  for(int i=0; i<Neig_; ++i){
    Field Av = opr->mult(evecs_[i]);
    evals_[i] = evecs_[i]*Av;    
    Av -= evals_[i]*evecs_[i]; 
    double res = Av.norm(); // residual 
    double vndif = evecs_[i].norm() -1.0;

    CCIO::cout<<" ["<<setw( 3)<<setiosflags(ios_base::left)<< i<<"] ";
    CCIO::cout<<      setw(25)<<setiosflags(ios_base::left)<< evals_[i];
    CCIO::cout<<"  "<<setw(25)<<setiosflags(ios_base::left)<< res;
    CCIO::cout<<"  "<<setw(25)<<setiosflags(ios_base::left)<< vndif <<"\n";
  }

  CCIO::cout<< resetiosflags(ios_base::scientific);
  CCIO::cout<<"---------------------------------------------";
  CCIO::cout<<"---------------------------------------------\n";
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

