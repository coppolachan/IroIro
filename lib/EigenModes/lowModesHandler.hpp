/*!@file lowmodeHandler.hpp
 * @brief declaration of the low modes handling classes
 */
#ifndef LOWMODESHANDLER_INCLUDED
#define LOWMODESHANDLER_INCLUDED

#include "include/pugi_interface.h"
#include <string>
#include <string.h>
#include "include/field.h"
#include "Communicator/comm_io.hpp"
#include <vector>
#include <cmath>

class Field;

class LowModesHandler{
private:
  int Neig_;
  double threshold_;
  std::vector<double> evals_;
  std::vector<Field> evecs_;
  void eval_input_txt(const std::string&);
  void evec_input_bin(const std::string&);
  std::vector<double> applied_evals_;

public:
  LowModesHandler(const XML::node eigen_node){
    std::string evecfile, evalfile;
    XML::node enode = eigen_node;

    if(enode !=NULL){
      CCIO::cout<<"LowModesHandler called.\n";
      const char* lmp_type = enode.attribute("type").value();
      double thre = 100.0; // default value should be greater than the theoretical maximum.
      XML::read(enode,"threshold", thre);// not mandatory: if not given, determined from the evalfile.
      threshold_= thre;

      XML::read(enode,"eval_file",evalfile,MANDATORY);
      eval_input_txt(evalfile);            // reading eval file. (Neig_ is determined here.)

      XML::read(enode,"evec_file",evecfile,MANDATORY);
      evec_input_bin(evecfile);            // reading evec file

      if(!strcmp(lmp_type, "inv")){        // list of 1/evals
	for(int i=0; i<Neig_; ++i) 
	  applied_evals_.push_back(1.0/evals_[i]);
      }else if (!strcmp(lmp_type, "sgn")){ // list of sgn(evals)
	for(int i=0; i<Neig_; ++i) 
	  applied_evals_.push_back(evals_[i]/fabs(evals_[i]));
      }else{
	CCIO::cout<< "LowModesHandler :type "<< lmp_type <<" not valid.\n";
	abort();
      }
    }
  }

  LowModesHandler(int Neigen,
		 const std::string& evecfile,
		 const std::string& evalfile,
		 const char* lmp_type="none"):Neig_(Neigen){
    if(strcmp(lmp_type, "none")){
      eval_input_txt(evalfile);
      evec_input_bin(evecfile);

      if(!strcmp(     lmp_type, "inv")){  // list of 1/evals
	for(int i=0; i<Neig_; ++i)
	  applied_evals_.push_back(1.0/evals_[i]);
      }else if(!strcmp(lmp_type, "sgn")){ // list of sgn(evals)
	for(int i=0; i<Neig_; ++i)
	  applied_evals_.push_back(evals_[i]/fabs(evals_[i]));
      }else{
	CCIO::cout<< "LowModesHandler: type "<< lmp_type <<"not valid.\n";
	abort();
      }
    }
  }

  const Field proj_high(const Field&) const;
  const Field proj_appliedLow(const Field&) const;

  const Field get_evec(int i){return evecs_[i];}
  double get_eval(int i){return evals_[i];}
  double get_appliedEval(int i){return applied_evals_[i];}
};

#endif
