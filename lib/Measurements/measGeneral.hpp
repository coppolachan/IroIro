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
  GaugeGlobal Uin_;

  int meas_num_;
  int starting_id_;
  int increment_;
  std::vector<std::string> config_list_;
  std::string output_prefix_;
  std::string gauge_prefix_;
  std::string seed_prefix_;

  bool file_list_;
  bool gauge_output_;
  bool seed_output_;

  void pre_process(GaugeField&,const RandNum&,int);
  void post_process(GaugeField&,const RandNum&,int);

public:
  MeasGeneral(XML::node node)
    :node_(node),meas_node_(node),Uin_(Geometry(node)),
     meas_num_(1),starting_id_(0),increment_(1),
     config_list_(0),
     output_prefix_("output_"),
     file_list_(false),gauge_output_(false),seed_output_(false){
    
    CCIO::cout<<"constructing MeasGeneral"<<std::endl;

    XML::node inode = node_;
    XML::descend(inode,"Configuration");
  
    /*! @brief RegularStep deals with config files with incremental numbers,
      while FileList accommodate a list of files with any name.
    */
    if(!XML::attribute_compare(inode,"Input","RegularStep")){
      std::string input_prefix;
      XML::read(inode,"input_prefix",input_prefix,MANDATORY);
      config_list_.push_back(input_prefix);   
      /*!< input_prefix is stored in the first element of config_list_*/
      
      if(XML::read(inode,"total_num",meas_num_))
	CCIO::cout<<"[default] meas_num_= "<<meas_num_<<std::endl;
      
      if(XML::read(inode,"starting_idx",starting_id_))
	CCIO::cout<<"[default] starting_id_= "<<starting_id_<<std::endl;
      
      if(XML::read(inode,"idx_increment",increment_))
	CCIO::cout<<"[default] inclement_= "  <<increment_  <<std::endl;
      
    }else if(!XML::attribute_compare(inode,"Input","FileList")){
      file_list_= true;
      for(XML::iterator it=inode.begin(); it!=inode.end();++it)
	config_list_.push_back(it->child_value());
    }else {
      abort();
    }
    
    RNG_Env::RNG = RNG_Env::createRNGfactory(node);// RNG is static value.
    
    XML::descend(meas_node_,"Measurement");
    
    /*! @brief When RegularStep is set, the numbering of the output files is 
      synchronized with that of config files. On the other hand, if FileList
      is set, the numbering rule can be specified in the Output section.
    */
    XML::node onode = node_;
    XML::descend(onode,"Output");
    if(file_list_){ 
      meas_num_= config_list_.size();
      if(XML::read(onode,"starting_idx",starting_id_))
	CCIO::cout<<"[default] starting_id_= "<<starting_id_<<std::endl;
    
      if(XML::read(onode,"idx_increment",increment_))
	CCIO::cout<<"[default] inclement_= "  <<increment_  <<std::endl;
    }
    
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
  for(int c=0; c<meas_num_; ++c){
    
    std::stringstream infile;
    if(file_list_) infile << config_list_[c];
    else           infile << config_list_[0]<< id;
    
    Uin_.initialize(node_,infile.str());
    CCIO::cout<<"Uin_ initialized from "<<infile.str()<<std::endl;
    
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
