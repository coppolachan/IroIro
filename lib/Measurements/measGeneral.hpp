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
  std::string output_prefix_;
  std::string gauge_prefix_;
  std::string seed_prefix_;

  bool file_list_;
  bool gauge_output_;
  bool seed_output_;

  std::vector<std::string> config_list_;
  std::vector<int> number_list_;

  void pre_process(GaugeField&,const RandNum&,int);
  void post_process(GaugeField&,const RandNum&,int);

  MeasGeneral(const MeasGeneral&);
  MeasGeneral& operator=(const MeasGeneral&);
  
public:

  //MeasGeneral::MeasGeneral(XML::node node,Geometry geom)
  MeasGeneral(XML::node node,Geometry geom)
  :node_(node),meas_node_(node),Uin_(geom),
   meas_num_(1),output_prefix_("output_"),
   file_list_(false),gauge_output_(false),seed_output_(false){

  RNG_Env::RNG = RNG_Env::createRNGfactory(node);// RNG is static value.

  CCIO::cout<<"constructing MeasGeneral"<<std::endl;
  
  XML::node inode = node_;
  XML::descend(inode,"Configuration");
  
  if(!XML::attribute_compare(inode,"Input","RegularStep")){
    /*!@brief RegularStep deals with config files with incremental numbers */ 

    std::string input_prefix;
    XML::read(inode,"input_prefix",input_prefix,MANDATORY);
    config_list_.push_back(input_prefix);   
    /*!< input_prefix is stored in the 1st element of config_list_*/

    std::string input_postfix ="";
    XML::read(inode,"input_postfix",input_postfix);
    config_list_.push_back(input_postfix);   
    /*!<input_postfix is optionally stored in the 2nd element of config_list_*/
      
    if(XML::read(inode,"total_num",meas_num_))
      CCIO::cout<<"[default] meas_num_= "<<meas_num_<<"\n";

    int starting=0;
    if(XML::read(inode,"starting_idx",starting))
      CCIO::cout<<"[default] starting_num= "<<starting<<"\n";
    int increment=1;
    if(XML::read(inode,"idx_increment",increment))
      CCIO::cout<<"[default] inclement_num= "<<increment<<"\n";

    for(int c=0; c<meas_num_; ++c) 
      number_list_.push_back(starting+increment*c);
    
  }else if(!XML::attribute_compare(inode,"Input","NumberList")){
    /*!@brief NumberList deals with a list of config number following prefix.*/
    std::string input_prefix;
    XML::read(inode,"input_prefix",input_prefix,MANDATORY);
    config_list_.push_back(input_prefix);   

    std::string input_postfix ="";
    XML::read(inode,"input_postfix",input_postfix);
    config_list_.push_back(input_postfix);   
    
    XML::read_array(inode,"numbers",number_list_,MANDATORY);
    meas_num_= number_list_.size();      

  }else if(!XML::attribute_compare(inode,"Input","FileList")){
    /*!@brief FileList deals with a list of files with any name. */

    for(XML::iterator it=inode.begin(); it!=inode.end();++it)
      config_list_.push_back(it->child_value());
    file_list_= true;
    meas_num_= config_list_.size();

    /*! @brief In this case, the number_list_ is only used for output.*/
    int starting=0;
    if(XML::read(inode,"starting_idx",starting))
      CCIO::cout<<"[default] output: starting_num = "<<starting<<"\n";
    int increment=1;
    if(XML::read(inode,"idx_increment",increment))
      CCIO::cout<<"[default] output: inclement_num= "<<increment<<"\n";

    for(int c=0; c<meas_num_; ++c) 
      number_list_.push_back(starting+increment*c);
  }else {
    abort();
  }

  ///
  XML::descend(meas_node_,"Measurement");
    
  XML::node onode = node_;
  XML::descend(onode,"Output");

  // output file for the measurement
  if(XML::read(onode,"output_prefix",output_prefix_))
    CCIO::cout<<"[default] output_prefix_= "<<output_prefix_<<std::endl;
  else CCIO::cout<<       "output_prefix_= "<<output_prefix_<<std::endl;
  
  // output of the pre-processed gauge config.
  if(XML::read(onode,"gauge_prefix",gauge_prefix_)){
    CCIO::cout<<"#gauge fixed config is NOT saved"<<std::endl;
  }else{
    gauge_output_= true;
    CCIO::cout<<"#gauge fixed config is saved in "<<gauge_prefix_<<std::endl;
  }
  // output of the rng seed
  if(XML::read(onode,"seed_prefix",seed_prefix_)){
    CCIO::cout<<"#rng seed is NOT saved"<<std::endl;
  }else{
    seed_output_= true;
    CCIO::cout<<"#rng seed is saved in "<<seed_prefix_<<std::endl;
  }
}
  //  MeasGeneral(XML::node node,Geometry geom);
  template <typename MeasObj> void do_meas();
};

template<typename MeasObj> void MeasGeneral::do_meas(){

  const RandNum* rng = RNG_Env::RNG->getRandomNumGenerator();

  for(int c=0; c<meas_num_; ++c){
    std::stringstream infile;
    if(file_list_) infile << config_list_[c];
    else           infile << config_list_[0]<< number_list_[c]
			  << config_list_[1];

    Uin_.initialize(node_,infile.str());
    CCIO::cout<<"Uin_ initialized from "<<infile.str()<<std::endl;
    
    pre_process(Uin_,*rng,number_list_[c]);// monitoring + smr + gfix
    //------- actual measurement -------
    std::stringstream outfile;
    outfile << output_prefix_<< number_list_[c];
    
    MeasObj meas(meas_node_,Uin_,*rng,outfile.str()); 

    CCIO::cout<<"Starting measurement"<<std::endl;
    meas.run();
    //----------------------------------
    post_process(Uin_,*rng,number_list_[c]); // seed saving
    CCIO::cout<<"\n";
  }
}
#endif
