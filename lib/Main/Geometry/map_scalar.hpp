#ifndef MAP_SCALAR_H_
#define MAP_SCALAR_H_

#include <vector>
#include "include/common_fields.hpp"

enum ShiftSign {Forward = 1, Backward = -1};

class Map{
  // Field offsets
  std::vector<int> site_map;
  int sites;

  Map(){}; // hide default constructor
public:
  Map(int, ShiftSign);

  template<class DATA, class FORMAT, class TAG>
  GeneralField<DATA,FORMAT,TAG> 
  operator()(const GeneralField<DATA,FORMAT,TAG>& InField) const{
    GeneralField<DATA,FORMAT,TAG> OutField;
    register int Nex = InField.get_Format().Nex();

    for (int index = 0; index < sites; ++index) {
      for (int outer_idx = 0; outer_idx < Nex; ++outer_idx){
	//very general without any assumption on ordering
	OutField.data[OutField.get_Format().islice(          index,outer_idx)] = 
	  InField.data[InField.get_Format().islice(site_map[index],outer_idx)];
      }
    }
    return OutField;
  }
  
};

#endif //MAP_SCALAR_H_
