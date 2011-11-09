/*!
 * @file geometry.h
 * @brief Defines the Geometry class
 */

#ifndef INIT_H_
#define INIT_H_

#include <iostream>
#include <vector>
#include <stdlib.h>

#include "include/commonPrms.h"
#include "Main/Geometry/siteIndex.h"
#include "include/pugi_interface.h"
#include "Communicator/comm_io.hpp"

/*!
 * @class Geometry
 * @brief Initializes several geometry related objects
 * 
 * From input files initializes the CommonPrms singleton\n
 * and the SiteIndex singleton classes.
 */
class Geometry {
  /*!
   * @brief initializer function
   */
  void initialize(pugi::xml_node node) {
    std::vector<int> Dim;
    std::vector<int> Node;

    Dim.resize(4);
    Node.resize(4);

    XML::read_array(node,"Lattice",Dim);
    XML::read_array(node,"Nodes",Node);
      
    Lattice t_latt = {0,
		      Dim[0],
		      Dim[1],
		      Dim[2], 
		      Dim[3],
		      Node[0],
		      Node[1],
		      Node[2],
		      Node[3]};
		      
    latt = t_latt;
    #ifndef HAVE_MPI
    // Check on number of cores - should be 1,1,1,1
    // for single core compilation
    if (latt.NPEx!=1 || latt.NPEy!=1  || latt.NPEz!=1  || latt.NPEt!=1 ) {
      std::cerr << "Number of nodes uncorrect for a single core" << std::endl;
      exit(1);
    }
     #endif

    parameters = CommonPrms::instance(latt);
    
    idx = SiteIndex::instance();
    
  }

 public:
  Lattice latt; /*!< Lattice structure containing the lattice dimensions 
		 * and MPI nodes distribution on 4 dimensions */ 
  SiteIndex* idx;  /*!< Singleton containing the geometry (nearest neighbors...) */
  CommonPrms* parameters; /*!< Singleton handling the global lattice parameters */

  /*!
    fdfsdjkflsdjkl
  */
  Geometry(pugi::xml_node node) {
    pugi::xml_node geom_node = node.child("Geometry");
    initialize(geom_node);
  };
  
  
};


#endif
