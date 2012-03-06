/*!

  @file mapper.cpp

  @brief Defined Mapping class functions for shifts

*/

#include "mapper.hpp"
#include "include/macros.hpp"


namespace MapsEnv{
  Mapper shift;

  void initialize_mapper() {
    shift.fill();
  }
}

void Mapper::fill(){
  for (int dir = 0; dir < NDIM_ ; ++dir ) {
    ShiftsMap.push_back(Map(dir, Backward));
    ShiftsMap.push_back(Map(dir, Forward));
  }
}
