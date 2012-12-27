/*! @file measGeneral.cpp
 *  @brief implementing member functions of the MeasGeneral class
 */
#include "measGeneral.hpp"
#include "GaugeM/staples.hpp"
#include "GaugeM/gfixFactory.hpp"
#include "Smearing/smearingFactories.hpp"
#include "Communicator/fields_io.hpp"
#include <iostream>
#include <string.h>
using namespace std;

void MeasGeneral::pre_process(GaugeField& U,const RandNum& rng,int id){

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

    SmearingOperatorFactory* SmrFactory = 
      SmearingOperators::createSmearingOperatorFactory(smr_node);
    Smear* SmearingObj = SmrFactory->getSmearingOperator();
    
    for(int i=0; i<Nsmear; ++i){ // Do the actual smearing 
      Ubuf = U;
      SmearingObj->smear(U,Ubuf);
    }
    CCIO::cout<<"Plaquette (smeared): "<<Staple.plaquette(U)<<endl;
  }
}

void MeasGeneral::post_process(GaugeField&,const RandNum& rng,int id){
  if(seed_output_){
    std::stringstream seed;
    seed << seed_prefix_<< id;
    rng.saveSeed(seed.str());
  }
}
