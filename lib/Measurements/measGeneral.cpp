/*! @file measGeneral.cpp
 *  @brief implementing member functions of the MeasGeneral class
 *
 * Time-stamp: <2014-11-12 17:50:22 cossu>
 */

#include "measGeneral.hpp"
#include "GaugeM/staples.hpp"
#include "GaugeM/gfixFactory.hpp"
#include "Smearing/smearingFactories.hpp"
#include "IO/fields_io.hpp"
#include <iostream>
#include <string.h>

using namespace std;

void MeasGeneral::setup(XML::node inode){
  // input files 
  XML::descend(inode,"Configuration");
  if(     !XML::attribute_compare(inode,"Input","RegularStep")) input_RegularStep(inode);
  else if(!XML::attribute_compare(inode,"Input","NumberList" )) input_NumberList(inode);
  else if(!XML::attribute_compare(inode,"Input","FileList"   )) input_FileList(inode);
  else{
    CCIO::cerr<<"No correct Configuration node is found\n";
     abort();
  }

  CCIO::cout<<"MeasGeneral::setup - Initializing Output\n";
  
  /// Sets output file for the measurements    
  XML::node onode = node_;
  XML::descend(onode,"Output");
  XML::read(onode,"append_mode",input_.app_out);

  if(XML::read(onode,"output_prefix",output_prefix_))
    CCIO::cout<<"[default] output_prefix_= "<<output_prefix_<<"\n";
  else CCIO::cout<<       "output_prefix_= "<<output_prefix_<<"\n";
  
  // output of the pre-processed gauge config.
  if(XML::read(onode,"gauge_prefix",gauge_prefix_)){
    CCIO::cout<<"Warning: gauge fixed configuration is NOT saved\n";
  }else{
    gauge_output_= true;
    CCIO::cout<<"gauge fixed config is saved in "<<gauge_prefix_<<"\n";
  }
  // output of the rng seed
  if(XML::read(onode,"seed_prefix",seed_prefix_)){
    CCIO::cout<<"Warning: random number seed is NOT saved\n";
  }else{
    seed_output_= true;
    CCIO::cout<<"random number seed is saved in "<<seed_prefix_<<"\n";
  }
}

/*--------------------- input with RegularStep ----------------------*/
void MeasGeneral::input_RegularStep(XML::node inode){
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
    CCIO::cout<<"[default] increment_num= "<<increment<<"\n";
  
  for(int c=0; c<meas_num_; ++c) number_list_.push_back(starting+increment*c);

  /// EigenModes section (optional input)
  for(XML::node eig_node = inode.child("EigenModes");
      eig_node;
      eig_node = eig_node.next_sibling("EigenModes")){

    eig_pred_.push_back(Eigen::predFactory(eig_node.child("ReadCondition")));

    XML::read(eig_node,"eigen_prefix",input_prefix,MANDATORY);
    eigen_list_.push_back(input_prefix);

    // Optional - declare explicitly the eigenvalues file list
    if (!XML::read(eig_node,"eval_prefix",input_prefix))
      eval_list_.push_back(input_prefix);


  }
}

/*--------------------- input with NumberList ----------------------*/
void MeasGeneral::input_NumberList(XML::node inode){
  /*!@brief NumberList deals with a list of config number following prefix.*/
  std::string input_prefix;
  XML::read(inode,"input_prefix",input_prefix,MANDATORY);
  config_list_.push_back(input_prefix);   
  
  std::string input_postfix ="";
  XML::read(inode,"input_postfix",input_postfix);
  config_list_.push_back(input_postfix);   
  
  XML::read_array(inode,"numbers",number_list_,MANDATORY);
  meas_num_= number_list_.size();      

  /// EigenModes section (optional input)  
  for(XML::node eig_node = inode.child("EigenModes");
      eig_node;
      eig_node = eig_node.next_sibling("EigenModes")){

    eig_pred_.push_back(Eigen::predFactory(eig_node.child("ReadCondition")));

    XML::read(eig_node,"eigen_prefix",input_prefix,MANDATORY);
    eigen_list_.push_back(input_prefix);
  }
}

/*--------------------- input with FileList ----------------------*/
void MeasGeneral::input_FileList(XML::node inode){
  /*!@brief FileList deals with a list of files with any name. 
    In this case, the content of number_list_ is not used.*/
  int starting=0;
  if(XML::read(inode,"starting_idx",starting))
    CCIO::cout<<"[default] output: starting_num = "<<starting<<"\n";

  int increment=1;
  if(XML::read(inode,"idx_increment",increment))
    CCIO::cout<<"[default] output: increment_num= "<<increment<<"\n";

  XML::node conf_node = inode;
  XML::descend(conf_node,"GaugeConfigs",MANDATORY);
  for(XML::iterator it=conf_node.begin(); it!=conf_node.end();++it)
    config_list_.push_back(it->child_value());
  file_list_= true;
  meas_num_= config_list_.size();

  for(int c=0; c<meas_num_; ++c)
    number_list_.push_back(starting+increment*c);

  /// EigenModes section (optional input)  
  for(XML::node eig_node = inode.child("EigenModes");
      eig_node;
      eig_node = eig_node.next_sibling("EigenModes")){

    eig_pred_.push_back(Eigen::predFactory(eig_node.child("ReadCondition")));

    XML::node ev_node = eig_node;
    XML::descend(ev_node,"EvecFiles",MANDATORY);
    for(XML::iterator it=ev_node.begin(); it!=ev_node.end();++it)
      eigen_list_.push_back(it->child_value());

    assert(eigen_list_.size() == meas_num_);
  }
}

void MeasGeneral::pre_process(GaugeField& U,const RandNum& rng,int id)const{
  Staples Staple;
  BoundaryCond* BC;
  CCIO::cout<< "Plaquette (thin)      : "<< Staple.plaquette(U) <<"\n";

  GaugeField Ubuf = U;

  //// gauge fixing ////
  XML::node gfix_node = node_;
  XML::descend(gfix_node,"GaugeFixing");
  GFixFactory* gffct = GaugeFix::createGaugeFixingFactory(gfix_node);
  GaugeFixing* gfix = gffct->getGaugeFixing(rng);
  
  U = gfix->do_fix(Ubuf);

  CCIO::cout<<"Plaquette (gauge fixed): "<<Staple.plaquette(U)<<endl;

  if(gauge_output_){
    std::stringstream gout;
    gout << gauge_prefix_<< id;
    if(CCIO::SaveOnDisk<Format::Format_G> (U.data,gout.str().c_str()))
      CCIO::cout << "Some error occurred in saving the file "<<gout.str() <<"\n";
  }

  //// smearing ////
  XML::node smr_node = node_; 
  XML::descend(smr_node,"Smearing"); 
  
  if(XML::attribute_compare(smr_node,"type","Off")){
    int Nsmear;                                    
    XML::read(smr_node,"Nsmear",Nsmear,MANDATORY);  

    SmearingFactory* SmrFactory = 
      Smearings::createSmearingFactory(smr_node);
    Smear* SmearingObj = SmrFactory->getSmearing();
    
    for(int i=0; i<Nsmear; ++i){ // Do the actual smearing 
      Ubuf = U;
      SmearingObj->smear(U,Ubuf);
    }
    CCIO::cout<<"Plaquette (smeared): "<<Staple.plaquette(U)<<endl;
  }

  /* Check if boundary conditions are requested */
  
  bool AntiPeriodicBC = false; // default
  XML::read(input_.node, "AntiPeriodicBC", AntiPeriodicBC);

  /* Apply boundary conditions AFTER THE PREPROCESS, if applicable */
  if (AntiPeriodicBC){
    CCIO::cout << "Applying antiperiodic Boundary conditions on direction T\n";
    BoundaryCond_antiPeriodic apbc(TDIR);
    apbc.apply_bc(*input_.config.gconf);
  }
}

void MeasGeneral::post_process(GaugeField&,const RandNum& rng,int id)const{
  if(seed_output_){
    std::stringstream seed;
    seed << seed_prefix_<< id;
    rng.saveSeed(seed.str());
  }
}

void MeasGeneral::post_process_last(GaugeField&,const RandNum& rng)const{
  if(seed_output_){
    std::stringstream seed;
    seed << seed_prefix_<<"last";
    rng.saveSeed(seed.str());
  }
}
