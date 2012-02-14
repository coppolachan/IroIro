/*!

  @file mapper.cpp

  @brief Defined Mapping class functions for shifts

*/

#include "include/macros.hpp"
#include "mapper.hpp"

Mapper::Mapper(){
  for (int dir = 0; dir < NDIM_ ; ++dir ) {
    ShiftsMap.push_back(Map(dir, Backward));
    ShiftsMap.push_back(Map(dir, Forward));
  }
}
