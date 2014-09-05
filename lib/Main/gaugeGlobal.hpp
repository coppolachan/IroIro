/*!
 * @file gaugeGlobal.hpp
 * @brief definition of the GaugeGlobal class
 *
 * Time-stamp: <2014-08-11 16:54:12 noaki>
 */
#ifndef GAUGEGLOBAL_INCLUDED
#define GAUGEGLOBAL_INCLUDED

#include "include/geometry.hpp"
#include "include/common_fields.hpp"
#include "HMC/TrajectoryInfo.hpp"

class RandNum;

class GaugeGlobal: public GaugeField{
  bool do_check_;
  GaugeGlobal(const GaugeGlobal&);//hide copy constructor
public:
  GaugeGlobal(const Geometry& geom,bool do_check=true)
  :GaugeField(geom.prms_->Nvol()),do_check_(do_check){}
  
  /*! Intialization with "Unit" */
  void initializeUnit();
  /*! Initializes the gauge field with random number */
  void initializeRand(const RandNum& rand);
  /*! Intialization with "TextFile" */
  void initializeTxt(const std::string &Filename);
  /*! Intialization with "Binary" */
  void initializeBin(const std::string &Filename);
  /*! Intialization with "CSDTbinary" */
  void initializeCSDTbin(const std::string &Filename);
  /*! Intialization with "JLQCDlegacy" */
  void initializeJLQCDlegacy(const std::string &Filename);
  /*! Intialization with "NERSC" */
  void initializeNERSC(const std::string &Filename);
 /*! Intialization with "ILDG" */
  void initializeILDG(const std::string &Filename);
 /*! Intialization with "ILDG_FORTRAN" */
  void initializeILDG_FORTRAN(const std::string &Filename);
  
  TrajInfo initialize(XML::node node);
  int initialize(XML::node node,std::string filename);
};

#endif 
