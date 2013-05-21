/*!
 * @file gaugeGlobal.hpp
 * @brief definition of the GaugeGlobal class
 *
 * Time-stamp: <2013-04-26 17:02:23 neo>
 */
#ifndef GAUGEGLOBAL_INCLUDED
#define GAUGEGLOBAL_INCLUDED

#include "include/geometry.hpp"
#include "include/common_fields.hpp"

class RandNum;

class GaugeGlobal: public GaugeField{
  GaugeGlobal(const GaugeGlobal&);//hide copy constructor
public:
  GaugeGlobal(const Geometry& geom):GaugeField(geom.prms_->Nvol()){}
  
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
  
  int initialize(XML::node node);
  int initialize(XML::node node,std::string filename);
};

#endif 
