/*! @file measGeneral.cpp
 *  @brief definition of the MeasGeneral class
 */
#ifndef MEASGENERAL_INCLUDED
#define MEASGENERAL_INCLUDED

#include "pugi_interface.h"
#include "common_fields.hpp"
#include "inputConfig.hpp"
#include "include/messages_macros.hpp"
#include "Tools/randNum_Factory.hpp"
#include "Main/gaugeGlobal.hpp"
#include "EigenModes/eigenModes.hpp"
#include "Dirac_ops/BoundaryConditions/boundaryCond.hpp"

#include <string>
#include <vector>
#include <sstream>

class RandNum;

namespace Measurements{
  struct Input{
    XML::node node;
    GaugeField* gconf;
    EigenModes* eigen;
    RandNum* rng;
    std::string output;
    
    Input(XML::node inode):node(inode),gconf(NULL),eigen(NULL),rng(NULL){
      XML::descend(node,"Measurement");
      _Message(DEBUG_VERB_LEVEL, "initialising RNG\n");
      RNG_Env::initialize(inode);
      rng = RNG_Env::RandNumG::instance().getRNG();
    }
    ~Input(){if(eigen) delete eigen;}
    InputConfig getConfig()const{ return InputConfig(gconf,eigen); }
  };
}

class MeasGeneral{
private:
  const XML::node node_;
  GaugeGlobal Uin_;

  Measurements::Input input_;
  
  int meas_num_;
  std::string output_prefix_;
  std::string gauge_prefix_;
  std::string seed_prefix_;

  bool file_list_;
  bool has_eigen_;
  bool gauge_output_;
  bool seed_output_;

  std::vector<std::string> config_list_;
  std::vector<std::string> eigen_list_;
  std::vector<int> number_list_;

  void pre_process(GaugeField&,const RandNum&,int)const;
  void post_process(GaugeField&,const RandNum&,int)const;
  
  void setup(XML::node);
  void input_RegularStep(XML::node);
  void input_NumberList(XML::node);
  void input_FileList(XML::node);

  MeasGeneral(const MeasGeneral&);   //prohibiting copy constructor
  MeasGeneral& operator=(const MeasGeneral&); //prohibiting substitution
public:
  MeasGeneral(XML::node node,Geometry geom,bool check_config=false)
    :node_(node),input_(node),Uin_(geom,check_config),
     meas_num_(1),output_prefix_("output_"),
     file_list_(false),
     gauge_output_(false),
     seed_output_(false),
     has_eigen_(false){ 

    _Message(DEBUG_VERB_LEVEL, "Constructing MeasGeneral object\n");
    setup(node);
    _Message(DEBUG_VERB_LEVEL, "MeasGeneral object constructed\n");
  }
  template <typename MeasObj> void do_meas();
};

template<typename MeasObj> void MeasGeneral::do_meas(){
  _Message(DEBUG_VERB_LEVEL, "Starting MeasGeneral::do_meas()\n");
  
  for(int c=0; c<meas_num_; ++c){

    /* setting  gauge configuration */
    std::stringstream infile;
    if(file_list_) infile<< config_list_[c];
    else           infile<< config_list_[0]<< number_list_[c]
			 << config_list_[1];

    Uin_.initialize(node_,infile.str());
    input_.gconf = &Uin_;
    CCIO::cout<<"Gauge configuration loaded from "<<infile.str()<<std::endl;

    /* setting corresponding eigenmodes */
    if(has_eigen_){ 
      std::stringstream evalfile;
      if(file_list_) evalfile<< eigen_list_[c];
      else           evalfile<< eigen_list_[0]<< number_list_[c]
			     << eigen_list_[2];
      std::stringstream evecfile;
      if(file_list_) evecfile<< eigen_list_[meas_num_+c];
      else           evecfile<< eigen_list_[1]<< number_list_[c]
			     << eigen_list_[3];
      
      input_.eigen->initialize<Format::Format_F>(evalfile.str(),
						 evecfile.str());
    }
    CCIO::cout<<"pre_process begins\n";
    pre_process(Uin_,*(input_.rng),number_list_[c]);// monitoring + smr + gfix
    CCIO::cout<<"pre_process ends\n";

    std::stringstream outfile;
    outfile << output_prefix_<< number_list_[c];
    input_.output = outfile.str();

    //------- actual measurement -------    
    CCIO::cout<<"Creating measurement"<<std::endl;
    MeasObj meas(input_); 
    CCIO::cout<<"Starting measurement"<<std::endl;
    meas.run();
    //----------------------------------
    post_process(Uin_,*(input_.rng),number_list_[c]); // seed saving
    CCIO::cout<<"\n";
  }
}
#endif
