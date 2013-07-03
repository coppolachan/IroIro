//--------------------------------------------------------------------
/*! @file heatbath.cpp
 *
 * @brief Definitions of functions for HeatBath update
 *
 * Time-stamp: <2013-06-04 15:24:20 neo>
 */
//--------------------------------------------------------------------
#include "heatbath.hpp"
#include "Communicator/comm_io.hpp"
#include "IO/fields_io.hpp"

#include <string>
#include <sstream>

void HeatBathGeneral::update(GaugeField& Uin)const{
  std::string seedfile_base = "seed_file_";
  std::string seedfile;
  std::stringstream seedfile_string;
  
  CCIO::cout << "Evolving using heatbath updates...\n";

  Action_->init(*rand_);

 // Thermalizations
  Action_->calc_H();
  for(int iter=1; iter <= Params_.ThermalizationSteps; ++iter){
    CCIO::cout << "-- # Thermalization step = "<< iter <<  "\n";

    double timer;
    TIMING_START;
    Action_->hb_update(*rand_);
    TIMING_END(timer);
    Action_->calc_H();
    CCIO::cout<< "Time for sweep (s) : "<< timer/1000.0 << "\n";
    //    Uin = md_->get_U();  //accept every time
  }

  // Actual updates
  for(int iter=Params_.StartingConfig; 
      iter < Params_.Nsweeps+Params_.StartingConfig; ++iter){
    CCIO::cout << "-- # Sweep = "<< iter <<  "\n";
    CCIO::cout << "---------------------------\n";
    double timer;
    TIMING_START;
    Action_->hb_update(*rand_);
    TIMING_END(timer);
    Action_->calc_H();
    CCIO::cout<< "Time for sweep (s) : "<< timer/1000.0 << "\n";

    if (Params_.SaveInterval!=0 && iter%Params_.SaveInterval == 0) {
      std::stringstream file;
      file << Params_.Filename_prefix << iter;
      CCIO::cout << "Saving configuration on: "<< file << "\n";
      if(CCIO::SaveOnDisk<Format::Format_G> (Uin.data, file.str().c_str())){
        CCIO::cout << "Some error occurred in saving file\n";
      }
      seedfile_string.str("");//clears up the string
      seedfile_string << seedfile_base << iter;
      CCIO::cout << "Saving seed file\n";
      rand_->saveSeed(seedfile_string.str());
    }
  }
  seedfile_string.str("");//clears up the string
  seedfile_string << seedfile_base << "last";
  rand_->saveSeed(seedfile_string.str());


}

