/*!@file eigenModes.hpp
 * @brief declaration of the classes to hold eigenmodes from files
 */
#ifndef EIGENMODES_INCLUDED
#define EIGENMODES_INCLUDED

#include "field.h"
#include "IO/fields_io.hpp"
#include <vector>
#include <string>
#include <string.h>

struct EigenModes{
  double thold_;
public:
  std::vector<double> evals_;
  std::vector<Field> evecs_;

  EigenModes(XML::node node):thold_(0){
    XML::read(node,"threshold",thold_); }
  
  EigenModes(double threshold = 100.0):thold_(threshold){}

  template<typename FMT>  
  void initialize(const std::string& evalfile,
		  const std::string& evecfile){
    evals_.clear();
    CCIO::cout<<"Reading eigenvalues from "<< evalfile<<".\n";
    std::ifstream reader(evalfile.c_str()); 
    int idummy, i=0;
    double eval=0.0;

    while(reader >> idummy && eval < thold_){
      reader >> eval;
      evals_.push_back(eval);
      i++;
    }

    evecs_.clear();
    CCIO::cout<<"Reading eigenvectors from "<<evecfile<<"\n";
    int Neig;
    CCIO::ReadFromDisk<FMT>(evecs_,evecfile.c_str(),Neig);
    CCIO::cout<< Neig << "eigenmodes are loaded."<<"\n";
  }
};

#endif
