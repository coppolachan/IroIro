//-------------------------------------------------------------
/*!@file   staggered_Utils.hpp
 * @brief  utilities for staggered fermions
 */
//-------------------------------------------------------------
#ifndef STAGGERED_UTILS_INCLUDED
#define STAGGERED_UTILS_INCLUDED

#include "include/format_G.h"

typedef Format::Format_G gfmt_t;     // link valiables (fundamental rep.)

namespace Dstagg{
  enum Dtype{DdagDee=0,DdagDoo=1,Dfull=2};

  void set_ksphase(std::valarray<double>& ev,std::valarray<double>& od,int Nv);
}
#endif
