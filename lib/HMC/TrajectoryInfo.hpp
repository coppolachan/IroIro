//--------------------------------------------------------------------
/*! @file TrajectoryInfo.hpp
 * @brief Declaration of structure to hold some info on the trajectory
 */
//--------------------------------------------------------------------
#ifndef TRAJINFO_INCLUDED
#define TRAJINFO_INCLUDED

struct TrajInfo{
  int SaveInterval;
  int ThermalizationSteps;
  int StartingConfig;
  std::string Filename_prefix;
  
  // Default constructor
  TrajInfo():SaveInterval(1),
	     ThermalizationSteps(0),
	     StartingConfig(1),
	     Filename_prefix("Conf_"){}; 

};

#endif
