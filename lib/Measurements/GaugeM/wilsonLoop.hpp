/*! @file wilsonLoop.hpp
 *  @brief calculation of the wilson loops with all possible areas
 */
#include "include/commonPrms.h"
#include "include/common_fields.hpp"
#include "Communicator/communicator.h"
#include "Main/Geometry/mapping.hpp"

class WilsonLoop{
  int mu_dir_;
  int mu_size_;
  Communicator* com_;
  
public:
  WilsonLoop(int mu)
    :mu_dir_(mu),
     mu_size_(CommonPrms::instance()->global_size(mu)),
     com_(Communicator::instance()){
    Mapping::init_shiftField();
  }
  
  std::vector<double> calc(const GaugeField&) const;
};


