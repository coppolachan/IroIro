/*!@file eigenCalcGeneral.hpp
 * @brief generalizes the eigenmodes calculation of fermion operators
 */
#ifndef EIGENCALCGENERAL_INCLUDED
#define EIGENCALCGENERAL_INCLUDED

#include "Dirac_ops/dirac_Operator_FactoryCreator.hpp"
#include "eigenModesSolver_Factory.hpp"
#include "eigenSorter_Factory.hpp"
#include "Communicator/comm_io.hpp"
#include "IO/fields_io.hpp"
#include "foprHermFunc.hpp"
#include "include/foprHermFactory.hpp"
#include <memory>
#include <vector>

class Field;
class Fopr_Herm;

class EigenCalcGeneral{
  std::auto_ptr<FoprHermFactory> opOrigFptr_; 
  std::auto_ptr<FoprHermFunc> opAccelFptr_; 
  std::auto_ptr<EigenSorterFactory> esortFptr_; 
  std::auto_ptr<EigenSolverFactory> eslvFptr_;
    
  std::vector<double> evals_;
  std::vector<Field> evecs_;
  int Neig_;

  FoprHermFactory* createFoprHermFactory(const XML::node&)const;
  FoprHermFunc* createAccelOpFunc(const XML::node&)const;
  EigenSorterFactory* createEigenSorterFactory(const XML::node&)const;
  void get_eval(const Fopr_Herm*);

public:
  EigenCalcGeneral(const XML::node& node);
  ~EigenCalcGeneral(){}

  void do_calc(InputConfig& input);
  void output_txt(const std::string&)const;

  template<typename FMT>
  void output_bin(const std::string&)const;

  template<typename FMT>
  void output_bin3D(const std::string&)const;
};

template<typename FMT>
void EigenCalcGeneral::output_bin(const std::string& output)const{
  std::string output_evals = output + "_evals.txt";
  std::ofstream writer(output_evals.c_str());

  for(int i=0; i<Neig_; ++i){
    CCIO::SaveOnDisk<FMT>(evecs_[i],output.c_str(),true);
    writer<< std::setw(2) <<setiosflags(std::ios_base::right)<< i;
    writer<< std::setw(25)<<std::setprecision(16)
	  <<setiosflags(std::ios_base::left )
	  << evals_[i]<<std::endl;
  }
}

template<typename FMT>
void EigenCalcGeneral::output_bin3D(const std::string& output)const{
  std::string output_evals = output + "_evals.txt";
  std::ofstream writer(output_evals.c_str());

  for(int i=0; i<Neig_; ++i){
    CCIO::SaveOnDisk3D<FMT>(evecs_[i],output.c_str(),true);
    writer<< std::setw(2) <<setiosflags(std::ios_base::right)<< i;
    writer<< std::setw(25)<<std::setprecision(16)
	  <<setiosflags(std::ios_base::left )
	  << evals_[i]<<std::endl;
  }
}

#endif
