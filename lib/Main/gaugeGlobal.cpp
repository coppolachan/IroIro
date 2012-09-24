#include "gaugeGlobal.hpp"
#include "gaugeConf.hpp" 
#include "include/geometry.hpp"

void GaugeGlobal::initializeUnit(){
  GaugeConf_unit gconf(format);
  gconf.init_conf(data);
}

void GaugeGlobal::initializeRand(const RandNum& rand){
  GaugeConf_rand gconf(format,rand);
  gconf.init_conf(data);
}

void GaugeGlobal::initializeTxt(const std::string &Filename){
  GaugeConf_txt gconf(format,Filename);
  gconf.init_conf(data);
}

void GaugeGlobal::initializeBin(const std::string &Filename){
  GaugeConf_bin gconf(format,Filename);
  gconf.init_conf(data);
}

void GaugeGlobal::initializeCSDTbin(const std::string &Filename){
  GaugeConf_csdt gconf(format,Filename);
  gconf.init_conf(data);
}

void GaugeGlobal::initializeJLQCDlegacy(const std::string &Filename) {
  GaugeConf_JLQCDLegacy gconf(format,Filename);
  gconf.init_conf(data);
}

int GaugeGlobal::initialize(XML::node node){
  try{
    XML::descend(node, "Configuration");
    if(!XML::attribute_compare(node,"Type","TextFile")){
      std::string filename(node.child_value());
      initializeTxt(filename);
      return 0;
    }
    if(!XML::attribute_compare(node,"Type","Unit")){
      initializeUnit();
      return 0;
    }
    if(!XML::attribute_compare(node,"Type","Binary")){
      std::string filename(node.child_value());
      initializeBin(filename);
      return 0;
    }
    if(!XML::attribute_compare(node,"Type","CSDTbinary")){
      std::string filename(node.child_value());
      initializeCSDTbin(filename);
      return 0;
    }
    if(!XML::attribute_compare(node,"Type","JLQCDlegacy")){
      std::string filename(node.child_value());
      initializeJLQCDlegacy(filename);
      return 0;
    }
    std::cout << "Configuration type unknown\n";
    exit(1);
  }catch(...) {
    std::cout << "Error in initialization of gauge field "<< std::endl;
    abort();
  }
  return 0;
}

int GaugeGlobal::initialize(XML::node node,std::string filename){
  try{
    XML::descend(node,"Configuration");
    if(!XML::attribute_compare(node,"Type","TextFile")){
      initializeTxt(filename);
      return 0;
    }
    if(!XML::attribute_compare(node,"Type","Binary")){
      initializeBin(filename);
      return 0;
    }
    if(!XML::attribute_compare(node,"Type","CSDTbinary")){
      initializeCSDTbin(filename);
      return 0;
    }
    if(!XML::attribute_compare(node,"Type","JLQCDlegacy")){
      initializeJLQCDlegacy(filename);
      return 0;
    }
    std::cout << "Configuration type unknown\n";
    exit(1);
  }catch(...) {
    std::cout << "Error in initialization of gauge field "<< std::endl;
  }
  return 0;
}

