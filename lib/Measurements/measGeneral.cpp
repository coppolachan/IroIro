/*! @file measGeneral.cpp
 *  @brief implementing member functions of the MeasGeneral class
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
   /// input files 
  XML::descend(inode,"Configuration");
  if(     !XML::attribute_compare(inode,"Input","RegularStep")) input_RegularStep(inode);
  else if(!XML::attribute_compare(inode,"Input","NumberList" )) input_NumberList(inode);
  else if(!XML::attribute_compare(inode,"Input","FileList"   )) input_FileList(inode);
  else{
     cerr<<"no correct Configuration node is found\n";
     abort();
  }

  CCIO::cout<<"initializing Output\n";
  /// output file for the measurement    
  XML::node onode = node_;
  XML::descend(onode,"Output");

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
    CCIO::cout<<"[default] inclement_num= "<<increment<<"\n";
  
  for(int c=0; c<meas_num_; ++c) number_list_.push_back(starting+increment*c);

  XML::node eig_node= inode.child("EigenModes");
  if(eig_node!= NULL){
    XML::node option_node = eig_node;
    XML::descend(option_node,"Options");
    input_.eigen = new EigenModes(option_node);
    has_eigen_= true;

    XML::read(eig_node,"eval_prefix",input_prefix,MANDATORY);
    eigen_list_.push_back(input_prefix);
    XML::read(eig_node,"evec_prefix",input_prefix,MANDATORY);
    eigen_list_.push_back(input_prefix);
    
    XML::read(eig_node,"eval_postfix",input_postfix);
    eigen_list_.push_back(input_postfix);
    XML::read(eig_node,"evec_postfix",input_postfix);
    eigen_list_.push_back(input_postfix);
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
  
  XML::node eig_node= inode.child("EigenModes");
  if(eig_node!= NULL){
    XML::node option_node = eig_node;
    XML::descend(option_node,"Options");
    input_.eigen = new EigenModes(option_node);
    has_eigen_= true;

    XML::read(eig_node,"eval_prefix",input_prefix,MANDATORY);
    eigen_list_.push_back(input_prefix);
    XML::read(eig_node,"evec_prefix",input_prefix,MANDATORY);
    eigen_list_.push_back(input_prefix);
    
    XML::read(eig_node,"eval_postfix",input_postfix);
    eigen_list_.push_back(input_postfix);
    XML::read(eig_node,"evec_postfix",input_postfix);
    eigen_list_.push_back(input_postfix);
  }
}

/*--------------------- input with FileList ----------------------*/
void MeasGeneral::input_FileList(XML::node inode){
  /*!@brief FileList deals with a list of files with any name. 
    In this case, the number_list_ is only used for output.*/
  int starting=0;
  if(XML::read(inode,"starting_idx",starting))
    CCIO::cout<<"[default] output: starting_num = "<<starting<<"\n";
  int increment=1;
  if(XML::read(inode,"idx_increment",increment))
    CCIO::cout<<"[default] output: inclement_num= "<<increment<<"\n";
  
  for(int c=0; c<meas_num_; ++c) number_list_.push_back(starting+increment*c);

  XML::node conf_node = inode;
  XML::descend(conf_node,"GaugeConfigs",MANDATORY);
  for(XML::iterator it=conf_node.begin(); it!=conf_node.end();++it)
    config_list_.push_back(it->child_value());
  file_list_= true;
  meas_num_= config_list_.size();
  
  XML::node eig_node= inode.child("EigenModes");
  if(eig_node!= NULL){
    XML::node ev_node = eig_node.child("Options");
    input_.eigen = new EigenModes(ev_node);
    has_eigen_= true;

    ev_node = eig_node;
    XML::descend(ev_node,"EvalFiles",MANDATORY);
    for(XML::iterator it=ev_node.begin(); it!=ev_node.end();++it)
      eigen_list_.push_back(it->child_value());

    ev_node = eig_node;
    XML::descend(ev_node,"EvecFiles",MANDATORY);
    for(XML::iterator it=ev_node.begin(); it!=ev_node.end();++it)
      eigen_list_.push_back(it->child_value());

    assert(eigen_list_.size()/2 == meas_num_);
  }
}

void MeasGeneral::pre_process(GaugeField& U,const RandNum& rng,int id)const{
  Staples Staple;
  CCIO::cout<< "Plaquette (thin): "<< Staple.plaquette(U) <<"\n";

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
      CCIO::cout << "Some error occurred in saving file\n";
  }

  //// smearing ////
  XML::node smr_node = node_; 
  XML::descend(smr_node,"Smearing"); 
  
  if(!XML::attribute_compare(smr_node,"type","Off")){
    return;
  }else{
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
}

void MeasGeneral::post_process(GaugeField&,const RandNum& rng,int id)const{
  if(seed_output_){
    std::stringstream seed;
    seed << seed_prefix_<< id;
    rng.saveSeed(seed.str());
  }
}
