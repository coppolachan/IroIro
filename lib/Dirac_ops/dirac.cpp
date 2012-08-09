#include "dirac.hpp"
#include "include/pugi_interface.h"

namespace Dop{
  double read_mass(const XML::node& node){
    double mass;
    XML::read(node, "mass", mass);
    return mass;
  }
}
