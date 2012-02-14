/*!

  @file mapper.hpp

  @brief Declares Mapping class for shifts

*/
#ifndef MAPPER_H_
#define MAPPER_H_

#include "include/commonPrms.h"
#include "include/common_fields.hpp"
#include "siteIndex.h"
#include "Communicator/comm_io.hpp"

#ifndef HAVE_MPI
#include "map_scalar.hpp"
#else
#include "map_parallelscalar.hpp"
#endif


class Mapper{
  std::vector<Map> ShiftsMap;

public:
  Mapper();

  template<class DATA, class FORMAT, class TAG>
  GeneralField<DATA,FORMAT,TAG>
  operator()(const GeneralField<DATA,FORMAT,TAG>& F, int dir, int sign){
    return ShiftsMap[(sign+1)>>1+2*dir](F);
  }


};


#endif //MAPPER_H_
