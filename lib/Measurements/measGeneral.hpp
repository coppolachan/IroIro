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
    InputConfig config;
    RandNum* rng;
    std::string output;
    bool app_out;

    Input(XML::node inode):node(inode),rng(NULL),app_out(false){
      XML::descend(node,"Measurement");
      _Message(DEBUG_VERB_LEVEL, "initialising RNG\n");
      RNG_Env::initialize(inode);
      rng = RNG_Env::RandNumG::instance().getRNG();
    }
    Field* getGconf()const{ return config.getGconf();}
    EigenModes* getEmodes(int e=0)const{return config.emodes[e];}
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
  bool gauge_output_;
  bool seed_output_;

  std::vector<std::string> config_list_;
  std::vector<std::string> eigen_list_;
  std::vector<std::string> eval_list_;
  std::vector<int> number_list_;

  std::vector<Eigen::Predic*> eig_pred_;

  void pre_process(GaugeField&,const RandNum&,int)const;
  void post_process(GaugeField&,const RandNum&,int)const;
  void post_process_last(GaugeField&,const RandNum&)const;

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
     seed_output_(false){

    _Message(DEBUG_VERB_LEVEL, "Constructing MeasGeneral object\n");
    setup(node);
    _Message(DEBUG_VERB_LEVEL, "MeasGeneral object constructed\n");
  }
  ~MeasGeneral(){
    for(int e=0; e<eig_pred_.size(); ++e) delete input_.config.emodes[e];
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
    
    /* Initialize the gauge configuration from the given filename */
    Uin_.initialize(node_,infile.str());
    input_.config.gconf = &Uin_;   
  
    CCIO::cout<<"Gauge configuration loaded from "<<infile.str()<<std::endl;

    /* setting corresponding eigenmodes */
    input_.config.emodes.clear();
    for(int e=0; e<eig_pred_.size(); ++e){
      std::stringstream eigfile, evalfile;
      if(file_list_) eigfile<< eigen_list_[e*meas_num_+c];
      else           eigfile<< eigen_list_[e]<< number_list_[c];

      if(!eval_list_.empty()){
	if (!file_list_) evalfile << eval_list_[e] << number_list_[c];
	input_.config.emodes.push_back(Eigen::initFromFile<Format::Format_F>(eig_pred_[e],evalfile.str(),eigfile.str()));
      } else {
	input_.config.emodes.push_back(Eigen::initFromFile<Format::Format_F>(eig_pred_[e],eigfile.str()));
      }
    }
    
    CCIO::cout<<"pre_process begins\n";
    pre_process(Uin_,*(input_.rng),number_list_[c]);// monitoring +smr +gfix
    CCIO::cout<<"pre_process ends\n";
    
    std::stringstream outfile;
    if(file_list_) {
      // Strip the directory name from the config_list_ entry - only when in FileList mode
      std::string filename = config_list_[c].substr(config_list_[c].find_last_of("/")+1,
						    config_list_[c].size());
      outfile<< output_prefix_<< filename;
    }
    else           outfile<< output_prefix_<< number_list_[c];
    

    input_.output = outfile.str();

    //------- actual measurement -------    
    CCIO::cout << "Creating measurement" << std::endl;
    MeasObj* meas = new MeasObj(input_); 
    CCIO::cout << "Starting measurement" << std::endl;
    meas->run();
    delete meas;
    //----------------------------------

    post_process(Uin_,*(input_.rng),number_list_[c]); // seed saving
    if(c == meas_num_-1)
      post_process_last(Uin_,*(input_.rng)); // seed saving for the last
    CCIO::cout<<"\n";



  }
}
#endif
