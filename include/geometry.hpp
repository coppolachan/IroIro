/*!
 * @file geometry.hpp
 * @brief Defines the Geometry class
 */
#ifndef GEOMETORY_H_
#define GEOMETORY_H_

#include "include/commonPrms.h"
#include "include/pugi_interface.h"
#include "Main/Geometry/siteIndex.hpp"
#include <iostream>
#include <vector>
#include "macros.hpp"

#ifdef IBM_BGQ_WILSON
#include "bgqwilson.h"
#endif

/*!
 * @class Geometry
 * @brief Initializes several geometry related objects
 * 
 * From input files initializes the CommonPrms singleton\n
 * and the SiteIndex singleton classes.
 */
class Geometry {
  /*! @brief Initializer function */
  void initialize(pugi::xml_node node){
    std::vector<int> Dim(NDIM_);
    std::vector<int> Node(NDIM_);

    XML::read_array(node,"Lattice",Dim,MANDATORY);
    XML::read_array(node,"Nodes",Node,MANDATORY);
      
    Lattice latt = {0,Dim[XDIR],Dim[YDIR],Dim[ZDIR],Dim[TDIR],
		    Node[XDIR],Node[YDIR],Node[ZDIR],Node[TDIR]};
#ifndef HAVE_MPI
    // Check on number of cores - should be 1,1,1,1
    // for single core compilation
    if(latt.NPEx!=1 || latt.NPEy!=1 || latt.NPEz!=1 || latt.NPEt!=1){
      std::cerr << "Number of nodes uncorrect for a single core" << std::endl;
      exit(1);
    }
#endif
    prms_= CommonPrms::instance(latt);
    //idx_= SiteIndex::instance();

    #ifdef IBM_BGQ_WILSON
    BGWilson_Init(Dim[XDIR],Dim[YDIR],Dim[ZDIR],Dim[TDIR],
    		  Node[XDIR],Node[YDIR],Node[ZDIR],Node[TDIR]);
    #endif

  }

public:
  //Lattice latt_; /*!< Lattice structure containing the lattice dimensions 
  //  		  * and MPI nodes distribution on 4 dimensions */ 
  CommonPrms* prms_; /*!< Singleton handling the general parameters */
  //SiteIndex* idx_;  /*!< Singleton containing the local geometry */
  
  /*! @brief Constructor - Initialized geometry object */
  Geometry(XML::node node){ 
    initialize(node.child("Geometry")); }  


};

#endif
