/*!@file eigenCalcGeneral.hpp
 * @brief generalizes the eigenmodes calculation of fermion operators
 */
#include "Dirac_ops/dirac_Operator_FactoryCreator.hpp"
#include "eigenModesSolver_Factory.hpp"
#include "eigenSorter_Factory.hpp"
#include "Communicator/comm_io.hpp"
#include "include/foprHermFactory.hpp"
#include <memory>
#include <vector>

class Field;
class Fopr_Herm;

class EigenCalcGeneral{
  std::auto_ptr<DiracWilsonLikeOperatorFactory> diracFptr_;
  std::auto_ptr<FoprHermFactory> opOrigFptr_; 
  std::auto_ptr<FoprHermFactory> opAccelFptr_; 
  std::auto_ptr<EigenSorterFactory> esortFptr_; 
  std::auto_ptr<EigenSolverFactory> eslvFptr_;
    
  std::vector<double> evals_;
  std::vector<Field> evecs_;
  int Neig_;

  FoprHermFactory* createAccelOpFactory(const XML::node&)const;
  FoprHermFactory* createFoprHermFactory(const XML::node&)const;
  EigenSorterFactory* createEigenSorterFactory(const XML::node&)const;
  void get_eval(const Fopr_Herm*);

public:
  EigenCalcGeneral(const XML::node& node);
  ~EigenCalcGeneral(){}

  void do_calc(Field* const conf);
  void output_txt(const std::string&)const;
  void output_bin(const std::string&)const;
};

