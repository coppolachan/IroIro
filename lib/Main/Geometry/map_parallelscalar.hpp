#ifndef MAP_PARALLELSCALAR_H_
#define MAP_PARALLELSCALAR_H_

#include <vector>
#include "include/common_fields.hpp"

enum ShiftSign {Forward = 1, Backward = -1};

typedef double basic_type;

class Map{
  // Field offsets
  std::vector<int> node_site_map;
  std::vector<int> bdry_site_map;
  std::vector<bool> bdry_site;

  int sites;
  int bdry_size;
  int direction;

  typedef void (Communicator::*comm_transfer_func)(basic_type*, basic_type*, int, int);
  comm_transfer_func comm_transfer;
  void transfer(basic_type*, basic_type*, size_t);

  Map(){}; // default constructor empty
public:
  Map(const int, const ShiftSign);
  
  template<class DATA, class FORMAT, class TAG>
  GeneralField<DATA,FORMAT,TAG> 
  operator()(const GeneralField<DATA,FORMAT,TAG>& InField){
    
    GeneralField<DATA,FORMAT,TAG> OutField;
    register int Nex = InField.get_Format().Nex();
    register int Nin = InField.get_Format().Nin();
    
    FORMAT bdry_form(bdry_size, Nex);
  
    basic_type* send_bdry = new basic_type[bdry_form.size()];
    basic_type* recv_bdry = new basic_type[bdry_form.size()];
 
    // Gather boundary
    for (int idx = 0, sidx = 0; idx < bdry_size; ++idx) {
      //CCIO::cout << "boundary site map  "<< bdry_site_map[idx] << "\n";
      for (int outer_idx = 0; outer_idx < Nex; ++outer_idx){
	for (int inner_idx = 0; inner_idx < Nin; ++inner_idx){
	  send_bdry[sidx++] = InField.data[InField.get_Format().index(inner_idx,
								      bdry_site_map[idx],
								      outer_idx)];
	}
      }
    }

    // Send
    transfer(recv_bdry,send_bdry,bdry_form.size());
    
    for (int idx = 0, ridx = 0; idx < sites; ++idx) {
      if (bdry_site[idx]) {
	// site on boundary
        for (int outer_idx = 0; outer_idx < Nex; ++outer_idx){
	  for (int inner_idx = 0; inner_idx < Nin; ++inner_idx){ 
	    // the ordering of sites in the two boundaries is the same
	    OutField.data.set(OutField.get_Format().index(inner_idx,
							  idx,
							  outer_idx), recv_bdry[ridx++]);
	    
	  }
	}
	
      }
      else{
	// site on bulk
	for (int outer_idx = 0; outer_idx < Nex; ++outer_idx){
	  OutField.data[OutField.get_Format().islice(               idx,outer_idx)] =
	    InField.data[InField.get_Format().islice(node_site_map[idx],outer_idx)];
	}
      }
    }
    
    //check
    /*
    for (int idx = 0; idx < sites; ++idx) {
      CCIO::cout << "Site: "<< idx<<"\n";
      for (int inner_idx = 0; inner_idx < Nin; ++inner_idx) {
	CCIO::cout << "Out["<<idx<<"]: "<<OutField.data[OutField.get_Format().index(inner_idx,idx,0)]
		   << " In["<<idx<<"]: "<< InField.data[OutField.get_Format().index(inner_idx,idx,0)]
		   << "\n";
	
      }
    }
    */
    return OutField;
  }
  
};

#endif //MAP_PARALLELSCALAR_H_
