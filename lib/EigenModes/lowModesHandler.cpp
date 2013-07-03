/*!@file lowModesHandler.cpp
 * @brief implementation of the low-mode handling class
*/

#include "lowModesHandler.hpp"
#include "Communicator/fields_io.hpp"
#include "Fields/field_expressions.hpp"

using namespace std;

namespace LowModes{
  
  const LowModesHandler* createHandler(XML::node node,
				       const EigenModes* const ems){

    if(     !XML::attribute_compare(node,"InnerLowModes","Inverse")) 
      return new LowModesHandler(ems,LowModes::Inverse());
    else if(!XML::attribute_compare(node,"InnerLowModes","Sign")) 
      return new LowModesHandler(ems,LowModes::Sign());
    else if(!XML::attribute_compare(node,"InnerLowModes","Off")) 
      return NULL;
    else abort();
  }
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
    sr[i]=eigs_->evec(i)*f;
    si[i]=eigs_->evec(i).im_prod(f);
  }
  for(int i=0; i<Neig_; ++i){
    f_h.add(re, -sr[i]*eigs_->evec(i,re)+si[i]*eigs_->evec(i,im));
    f_h.add(im, -sr[i]*eigs_->evec(i,im)-si[i]*eigs_->evec(i,re));
  }
  return f_h;
}

const Field LowModesHandler::proj_low(const Field& f) const{
  size_t size = f.size();
  assert(size%2 ==0);
  Field f_l(size,0.0);

  std::slice re(0,size/2,2);
  std::slice im(1,size/2,2);

  vector<double> sr(Neig_);
  vector<double> si(Neig_);  

  for(int i=0; i<Neig_; ++i){
    sr[i]=eigs_->evec(i)*f;
    si[i]=eigs_->evec(i).im_prod(f);
  }
  for(int i=0; i<Neig_; ++i){
    f_l.add(re,evals_[i]*(sr[i]*eigs_->evec(i,re)-si[i]*eigs_->evec(i,im)));
    f_l.add(im,evals_[i]*(sr[i]*eigs_->evec(i,im)+si[i]*eigs_->evec(i,re)));
  }
  return f_l;
}

