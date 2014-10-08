/*!@file eigenCalcGeneral.hpp
 * @brief generalizes the eigenmodes calculation of fermion operators
 */
#ifndef EIGENCALCGENERAL_INCLUDED
#define EIGENCALCGENERAL_INCLUDED

#include "Dirac_ops/dirac_Operator_FactoryCreator.hpp"
#include "eigenModesSolver_Factory.hpp"
#include "eigenSorter.hpp"
#include "Communicator/comm_io.hpp"
#include "IO/fields_io.hpp"
#include "Fopr/foprHermFuncFactory.hpp"
#include "Fopr/foprHermFactory.hpp"
#include <memory>
#include <vector>

class Field;
class Fopr_Herm;

class OpFunc{
public:
  virtual Fopr_Herm* getOp(const Fopr_Herm*)const = 0;
  virtual ~OpFunc(){}
};

class EigenCalcGeneral{
  std::auto_ptr<EigenSorter>     esortPtr_;
  std::auto_ptr<OpFunc>          opAccelFptr_; 
  std::auto_ptr<FoprHermFactory> opOrigFptr_; 
  std::auto_ptr<EigenSolverFactory> eslvFptr_;
  
  std::vector<double> evals_;
  std::vector<Field> evecs_;
  int Neig_;

  bool useExModes_;
  int esId_;
  int Nex_;

  EigenSorter* eigenSorterFactory(const XML::node&)const;
  EigenModesSolver* eigSlvFactory(const Fopr_Herm*,const EigenSorter*);

  void get_eval(const Fopr_Herm*);
  void test_exModes(const Fopr_Herm*,const EigenModes* em);

  XML::node eslvNode_;
public:
  EigenCalcGeneral(XML::node node);

  void do_calc(const InputConfig& input);
  void output_txt(const std::string&,bool append = false)const;

  template<typename FMT>
  void output_bin(const std::string&,bool append = false)const;

  template<typename FMT>
  void output_bin3D(const std::string&,bool append = false)const;
};

template<typename FMT>
void EigenCalcGeneral::output_bin(const std::string& output,bool append)const{
  std::string output_evals = output +"_evals.txt";

  std::ofstream writer(output_evals.c_str());
  if(append) writer.open(output_evals.c_str(),std::ios_base::app);

  for(int i=0; i<Neig_; ++i){
    CCIO::SaveOnDisk<FMT>(evecs_[i],output.c_str(),append);
    writer<< std::setw(2) <<setiosflags(std::ios_base::right)<< Nex_+i;
    writer<< std::setw(25)<<std::setprecision(16)
	  << setiosflags(std::ios_base::left )
	  << evals_[i]<<std::endl;
  }
}

template<typename FMT>
void EigenCalcGeneral::output_bin3D(const std::string& output,bool append)const{
  std::string output_evals = output +"_evals.txt";

  std::ofstream writer(output_evals.c_str());
  if(append) writer.open(output_evals.c_str(),std::ios_base::app);

  for(int i=0; i<Neig_; ++i){
    CCIO::SaveOnDisk3D<FMT>(evecs_[i],output.c_str(),append);
    writer<< std::setw(2) <<setiosflags(std::ios_base::right)<< Nex_+i;
    writer<< std::setw(25)<<std::setprecision(16)
	  << setiosflags(std::ios_base::left )
	  << evals_[i]<<std::endl;
  }
}

#endif
