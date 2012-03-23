/*! @file mapper.hpp
    @brief Declares Mapping class for shifts
*/
#ifndef MAPPER_H_
#define MAPPER_H_

#include "include/common_fields.hpp"

#ifndef HAVE_MPI
#include "map_scalar.hpp"
#else
#include "map_parallelscalar.hpp"
#endif

class Mapper{
  std::vector<Map> ShiftsMap;

public:
  Mapper(){};
  void fill();

  template<class DATA, class FORMAT, class TAG>
  GeneralField<DATA,FORMAT,TAG>
  operator()(const GeneralField<DATA,FORMAT,TAG>& F, int dir, int sign)const {
    CCIO::cout<<"Mapper::operator()"<<std::endl;
    return ShiftsMap[((sign+1)>>1)+2*dir](F);
  }
};

namespace MapsEnv{
  extern Mapper shift;
  void initialize_mapper();
}

#endif //MAPPER_H_
