/*!@file lowmodeHandler.hpp
 * @brief declaration of the low modes handling classes
 */
#ifndef LOWMODESHANDLER_INCLUDED
#define LOWMODESHANDLER_INCLUDED

#include "pugi_interface.h"
#include "field.h"
#include "Communicator/comm_io.hpp"
#include "EigenModes/eigenModes.hpp"
#include <vector>
#include <cmath>

namespace LowModes{
  struct Inverse{};
  struct Sign{};
  const LowModesHandler* createHandler(XML::node,const EigenModes* const);
}

class LowModesHandler{
private:
  const EigenModes* const eigs_;
  int Neig_;
  std::vector<double> evals_;
public:
  LowModesHandler(const EigenModes* const ems,LowModes::Inv)
    :eigs_(ems),Neig_(ems->size()){
    for(int i=0; i<Neig_; ++i) evals_.push_back(1.0/eigs_->eval(i));
  }
  LowModesHandler(const EigenModes* const ems,LowModes::Sgn)
    :eigs_(ems),Neig_(ems->size()){
    for(int i=0; i<Neig_; ++i) 
      evals_.push_back(eigs_->eval(i)/fabs(eigs_->eval(i)));
  }

  const Field proj_high(const Field&) const;
  const Field proj_low(const Field&) const;
};

#endif
