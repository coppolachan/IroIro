/*!@file eigenCalcGeneral.hpp
 * @brief generalizes the eigenmodes calculation of fermion operators
 */
#include "Dirac_ops/dirac_Operator_Factory.hpp"
#include "eigenModesSolver_Factory.hpp"
#include "eigenSorter_Factory.hpp"
#include "Communicator/comm_io.hpp"
#include "include/foprHermFactory.hpp"
#include <memory>
#include <vector>

class Field;
class Fopr_Herm;

class EigenCalcGeneral{
  std::auto_ptr<DiracWilsonLikeOperatorFactory> diracFact_;
  std::auto_ptr<FoprHermFactory> opOrigFact_; 
  std::auto_ptr<FoprHermFactory> opAccelFact_; 
  std::auto_ptr<EigenSorterFactory> esortFact_; 
  std::auto_ptr<EigenSolverFactory> eslvFact_;
    
  std::vector<double> eval_;
  std::vector<Field> evec_;
  int Neig_;

  FoprHermFactory* createAccelOpFactory(const XML::node&)const;
  FoprHermFactory* createFoprHermFactory(const XML::node&)const;
  EigenSorterFactory* createEigenSorterFactory(const XML::node&)const;
  void get_eval(const Fopr_Herm*);

public:
  EigenCalcGeneral(const XML::node& node);

  void do_calc(Field* const conf);
  void output(std::ofstream&);
};

