/*!
 * @file gaugeGlobal.cpp
 * @brief Declaration of the GaugeGlobal class
 *
 * Time-stamp: <2013-09-25 13:33:28 cossu>
 */

#include "gaugeGlobal.hpp"
#include "gaugeConf.hpp" 
#include "include/errors.hpp"

#include <dirent.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

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

void GaugeGlobal::initializeNERSC(const std::string &Filename) {
  GaugeConf_NERSC gconf(format,Filename);
  gconf.init_conf(data);
}

void GaugeGlobal::initializeILDG(const std::string &Filename) {
  GaugeConf_ILDG gconf(format,Filename);
  gconf.init_conf(data);
}


TrajInfo GaugeGlobal::initialize(XML::node node){
  using namespace std;

  TrajInfo Info;
  XML::node top_node = node;
  XML::descend(node, "Configuration");
  string filename(node.child_value());
  
  
  if(!XML::attribute_compare(node,"Class","Trajectory")){
    // Search for a file starting with the node content
    // and appending the trajectory number at the end
    
    // Separate the filename into directory and prepend
    const char* str = filename.c_str();
    const char* token;
    string directory, prepend ;
    string filepath, fname;

    if ((token = strrchr(str, '/'))!=NULL){
      directory = filename.substr(0,token-str+1);
      prepend = filename.substr(token-str+1, filename.length());
    } else {
      directory = "./";
      prepend = filename;
    }
 
    string latest_conf = prepend;

    Info.Filename_prefix = filename;

    // Looks for the latest created file starting 
    // with filename

    DIR *dir;
    struct dirent *ent;
    time_t latest_access = 0;
    if ((dir = opendir (directory.c_str())) != NULL) {
      // saves all the filenames in the dir 
      struct stat filestat;
      while (ent = readdir (dir)){
	fname = ent->d_name;
	filepath = directory.c_str()  + fname;
	stat( filepath.c_str(), &filestat );
	if (S_ISREG(filestat.st_mode) && fname.length() >= prepend.length()){
	  //starts with [filename]
	  if (0 == fname.compare(0, prepend.length(), prepend))
	    if (filestat.st_mtime > latest_access){
	      latest_access = filestat.st_mtime; 
	      latest_conf = fname;
	      Info.StartingConfig = 
		atoi(latest_conf.substr(prepend.length(),latest_conf.length()).c_str())+1;
	    }
	}
      }
      closedir (dir);
    } else {
      // Directory not found
      ostringstream msg;
      msg << "The directory "<< directory << " was not found in your path.\n";
      Errors::IOErr(msg);
    }
    filename = directory.c_str() + latest_conf;
  } 
  
  int status = initialize(top_node, filename);

  return Info;
}

int GaugeGlobal::initialize(XML::node node,std::string filename){
  try{
    XML::descend(node,"Configuration");
   if(!XML::attribute_compare(node,"Type","Unit")){
      initializeUnit();
      return 0;
    }
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
    if(!XML::attribute_compare(node,"Type","NERSC")){
      initializeNERSC(filename);
      return 0;
    }
    if(!XML::attribute_compare(node,"Type","ILDG")){
      initializeILDG(filename);
      return 0;
    }
    
    std::ostringstream msg;
    msg << "Configuration type unknown\n";
    Errors::BaseErr("GaugeGlobal::initialize",msg);
    
  }catch(...) {
    std::ostringstream msg;
    msg << "Error in initialization of the Gauge field\n ";
    Errors::BaseWarning("GaugeGlobal::initialize",msg);
    return -1;
  }
  return 0;
}

