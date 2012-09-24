/*! @file measGeneral.cpp
 *  @brief definition of the MeasGeneral class
 */
#ifndef MEASGENERAL_INCLUDED
#define MEASGENERAL_INCLUDED

#include "include/pugi_interface.h"
#include "include/common_fields.hpp"
#include "Tools/randNum_Factory.h"
#include "Main/gaugeGlobal.hpp"
#include <string>
#include <vector>
#include <sstream>

class RandNum;

class MeasGeneral{
private:
  const XML::node node_;
  XML::node meas_node_;
  std::vector<std::string> config_list_;
  GaugeGlobal Uin_;
  
  int starting_id_;
  int increment_;
  std::string output_prefix_;
  std::string gauge_prefix_;
  std::string seed_prefix_;

  bool gauge_output_;
  bool seed_output_;

  void list_config();
  void pre_process(GaugeField&,const RandNum&,int);
  void post_process(GaugeField&,const RandNum&,int);

public:
  MeasGeneral(XML::node node)
    :node_(node),meas_node_(node),Uin_(Geometry(node)),
     starting_id_(0),increment_(1),
     output_prefix_("output_"),
     gauge_output_(false),seed_output_(false){

    CCIO::cout<<"constructing MeasGeneral"<<std::endl;

    list_config();       //initialization of config_list_
    RNG_Env::RNG = RNG_Env::createRNGfactory(node);// RNG is static value.

    XML::descend(meas_node_,"Measurement");

    XML::node onode = node;
    XML::descend(onode,"Output");
    if(XML::read(onode,"starting_idx",starting_id_))
      CCIO::cout<<"[default] starting_id_= "<<starting_id_<<std::endl;
    if(XML::read(onode,"idx_increment",increment_))
      CCIO::cout<<"[default] inclement_= "  <<increment_  <<std::endl;

    // output file for the measurement
    if(XML::read(onode,"output_prefix",output_prefix_))
      CCIO::cout<<"[default] output_prefix_= "<<output_prefix_<<std::endl;
    else CCIO::cout<<       "output_prefix_= "<<output_prefix_<<std::endl;
    
    // output of the pre-processed gauge config.
    if(XML::read(onode,"gauge_prefix",gauge_prefix_)){
      CCIO::cout<<"#precond gauge is NOT saved"<<std::endl;
    }else{
      gauge_output_= true;
      CCIO::cout<<"#precond gauge is saved in "<<output_prefix_<<std::endl;
    }
    // output of the rng seed
    if(XML::read(onode,"seed_prefix",seed_prefix_)){
      CCIO::cout<<"#rng seed is NOT saved"<<std::endl;
    }else{
      seed_output_= true;
      CCIO::cout<<"#rng seed is saved in "<<seed_prefix_<<std::endl;
    }
  }

  template <typename MeasObj> void do_meas();
};

template<typename MeasObj> void MeasGeneral::do_meas(){

  const RandNum* rng = RNG_Env::RNG->getRandomNumGenerator();
  
  int id = starting_id_;
  
  for(int c=0; c<config_list_.size(); ++c){

    Uin_.initialize(node_,config_list_[c]);
    CCIO::cout<<"Uin_ initialized from "<<config_list_[c]<<std::endl;
    
    pre_process(Uin_,*rng,id);  // monitoring + smr + gfix
    //------- actual measurement -------
    std::stringstream outfile;
    outfile << output_prefix_<< id;

    MeasObj meas(meas_node_,Uin_,*rng,outfile.str()); 

    CCIO::cout<<"Starting measurement"<<std::endl;
    meas.run();
    //----------------------------------
    post_process(Uin_,*rng,id); // seed saving
    CCIO::cout<<"\n";

    id += increment_;
  }
}
#endif
