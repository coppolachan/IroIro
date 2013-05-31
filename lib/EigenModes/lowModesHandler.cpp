/*!@file lowModesHandler.cpp
 * @brief implementation of the low-mode handling class
*/

#include "lowModesHandler.hpp"
#include "Communicator/fields_io.hpp"
#include "Fields/field_expressions.hpp"

using namespace std;

void LowModesHandler::eval_input_txt(const string& input){
  ifstream reader(input.c_str()); 
  int idummy;
  double eval=0.0;
  CCIO::cout<<"Reading eigenvalues from "<< input << ".\n";
  int i=0;
  while(reader >> idummy && eval < threshold_){
    reader >> eval;
    evals_.push_back(eval);
    //CCIO::cout << idummy << "  " << evals_[i]<<endl;
   i++;
  }
  CCIO::cout << i << " eigenvalues are loaded." <<endl;
  Neig_=i;
}

void LowModesHandler::evec_input_bin(const string& input){
  CCIO::cout<<"LowModesHandler : Reading eigenvectors from "<< input << endl;
  CCIO::ReadFromDisk<Format::Format_F>(evecs_,input.c_str(),Neig_);
  CCIO::cout<< Neig_ << "eigenmodes are loaded."<< endl;
}

const Field LowModesHandler::proj_high(const Field& f) const{
  Field f_h = f;
  size_t size = f.size();
  assert(size%2 ==0);

  std::slice re(0,size/2,2);
  std::slice im(1,size/2,2);

  vector<double> sr(Neig_);
  vector<double> si(Neig_);

  for(int i=0; i<Neig_; ++i){
    sr[i]=evecs_[i]*f;
    si[i]=evecs_[i].im_prod(f);
  }
  for(int j=0; j<Neig_; ++j){
    f_h.add(re, -sr[j]*evecs_[j][re]+si[j]*evecs_[j][im]);
    f_h.add(im, -sr[j]*evecs_[j][im]-si[j]*evecs_[j][re]);
  }
  return f_h;
}

const Field LowModesHandler::proj_appliedLow(const Field& f) const{
  size_t size = f.size();
  assert(size%2 ==0);
  Field f_l(size,0.0);

  std::slice re(0,size/2,2);
  std::slice im(1,size/2,2);

  vector<double> sr(Neig_);
  vector<double> si(Neig_);  

  for(int i=0; i<Neig_; ++i){
    sr[i]=evecs_[i]*f;
    si[i]=evecs_[i].im_prod(f);
  }
  for(int i=0; i<Neig_; ++i){
    f_l.add(re, applied_evals_[i]*(sr[i]*evecs_[i][re]-si[i]*evecs_[i][im]));
    f_l.add(im, applied_evals_[i]*(sr[i]*evecs_[i][im]+si[i]*evecs_[i][re]));
  }
  return f_l;
}

