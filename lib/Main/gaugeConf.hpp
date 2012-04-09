/*! 
  @file gaugeConf.hpp
  @brief Declaration of classes to manage gauge configurations
*/
#ifndef GAUGECONF_INCLUDED
#define GAUGECONF_INCLUDED

#include "Main/Geometry/siteIndex.hpp"
#include "include/format_G.h"
#include "Tools/randNum_MP.h"

class Field;
class SiteIndex;

class GaugeConf{
protected:
  const Format::Format_G& fg_;
  int Nx_, Ny_, Nz_, Nt_;
  int Nc_;
  int Ndim_;
  int getid(int x,int y,int z,int t);
public:
  GaugeConf(const Format::Format_G& fg)
    :fg_(fg),
     Nx_(CommonPrms::instance()->Nx()),
     Ny_(CommonPrms::instance()->Ny()),
     Nz_(CommonPrms::instance()->Nz()),
     Nt_(CommonPrms::instance()->Nt()),
     Nc_(CommonPrms::instance()->Nc()),
     Ndim_(CommonPrms::instance()->Ndim()){}
  virtual ~GaugeConf(){}
  virtual void init_conf(Field&) = 0;
};  

class GaugeConf_bin:public GaugeConf{
  std::string file;
public:
  GaugeConf_bin(const Format::Format_G& fg, const std::string& fname)
    :GaugeConf(fg), file(fname){}
  void init_conf(Field&);
};

class GaugeConf_JLQCDLegacy:public GaugeConf{
  std::string file;
public:
  GaugeConf_JLQCDLegacy(const Format::Format_G& fg, const std::string& fname)
    :GaugeConf(fg), file(fname){}
  void init_conf(Field&);
};

class GaugeConf_txt:public GaugeConf{
  std::string file;
public:
  GaugeConf_txt(const Format::Format_G& fg, const std::string& fname)
    :GaugeConf(fg),file(fname){}
  void init_conf(Field&);
};

class GaugeConf_unit:public GaugeConf{
public:
  GaugeConf_unit(const Format::Format_G& fg):GaugeConf(fg){}
  void init_conf(Field&);
};

class GaugeConf_rand:public GaugeConf{
  const RandNum& rand_;
public:
  GaugeConf_rand(const Format::Format_G& fg,const RandNum& rand)
    :GaugeConf(fg),rand_(rand){}
  void init_conf(Field&);
};
#endif
