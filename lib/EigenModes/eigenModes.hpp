/*!@file eigenModes.hpp
 * @brief declaration of the classes to hold eigenmodes from files
 */
#ifndef EIGENMODES_INCLUDED
#define EIGENMODES_INCLUDED

#include "field.h"
#include "IO/fields_io.hpp"
#include <vector>
#include <string>
#include <string.h>


/// eigenmodes data structure ///
struct EigenModes{
  std::vector<double> evals;
  std::vector<Field> evecs;
  int size()const{ return evals.size();}
};

/// predicates to determine how to read the data ///
namespace Eigen{

  static std::string eval_ext = "_evals.txt";

  class Predic{
  public:
    virtual bool operator()(double)const{return true;}
    virtual ~Predic(){}
  };
  
  class BelowEq :public Predic{
    double thold_;
  public:
    BelowEq(XML::node node){ XML::read(node,"threshold",thold_);}
    BelowEq(double thold):thold_(thold){}
    bool operator()(double a)const{return a <= thold_;}
  };

  class Beyond :public Predic{
    double thold_;
  public:
    Beyond(XML::node node){ XML::read(node,"threshold",thold_);}
    Beyond(double thold):thold_(thold){}
    bool operator()(double a)const{return a > thold_;}
  };

  class SmallerEq :public Predic{
    double thold_;
  public:
    SmallerEq(XML::node node){ XML::read(node,"threshold",thold_);}
    SmallerEq(double thold):thold_(thold){}
    bool operator()(double a)const{return fabs(a) <= thold_;}
  };

  class LargerThan :public Predic{
    double thold_;
  public:
    LargerThan(XML::node node){ XML::read(node,"threshold",thold_);}
    LargerThan(double thold):thold_(thold){}
    bool operator()(double a)const{return fabs(a) > thold_;}
  };

  Predic* predFactory(XML::node node);


  /// functions to load eigenmodes from files ///

  template<typename FMT>  
  void initFromFile(EigenModes& emodes,const Predic* filter,
		    const std::string& evalfile,
		    const std::string& evecfile){

    std::ifstream evf(evalfile.c_str()); 
    if(evf.good()){
      CCIO::cout<<"Reading eigenvalues from "<< evalfile<<"\n";
    }else{
      CCIO::cout<<"Reading from "<< evalfile<<" failed.\n";
      abort();
    }
    emodes.evals.clear();
    int id;
    double ev;

    evf>>id; evf>>ev;
    while(!evf.eof()){
      if((*filter)(ev)) emodes.evals.push_back(ev);
      evf>>id; evf>>ev;
    }
    int Neig = emodes.evals.size();
    emodes.evecs.clear();
    CCIO::cout<<"Reading eigenvectors from "<< evecfile<<"\n";
    CCIO::ReadFromDisk<FMT>(emodes.evecs,evecfile.c_str(),Neig);
    CCIO::cout<<Neig<<" eigenmodes are loaded.\n";
  }

  template<typename FMT>  
  EigenModes* initFromFile(const Predic* filter,
			   const std::string& evalfile,
			   const std::string& evecfile){
    EigenModes* emode = new EigenModes;
    initFromFile<FMT>(*emode,filter,evalfile,evecfile);
    return emode;
  }

  template<typename FMT>  
  EigenModes* initFromFile(const Predic* filter,const std::string& eigfile){
    std::string evalfile = eigfile+eval_ext;
    return initFromFile<FMT>(filter,evalfile,eigfile);
  }

  template<typename FMT>  
  int pickUpFromFile(double& eval,Field& evec,int ireq,
		     const Predic* filter,
		     const std::string& evalfile,
		     const std::string& evecfile){
    FMT fmt(CommonPrms::instance()->Nvol());

    std::ifstream evf(evalfile.c_str()); 
    if(evf.good()){
      CCIO::cout<<"Reading eigenvalue"<< ireq <<" from "<< evalfile<<"\n";
    }else{
      CCIO::cout<<"Reading from "<< evalfile<<" failed.\n";
      abort();
    }
    
    int id;
    double ev;
    while(!evf.eof()){
      evf>>id; evf>>ev;
      if(id == ireq && (*filter)(ev)){
	eval = ev;
	CCIO::cout<<"Eigenvalue: "<<eval<<"\n";
	break;
      }
    }
    if(evf.eof()){
      CCIO::cout<<"Required eigenmode was not found.\n";
      return 1;
    }

    uint64_t offset 
      = sizeof(double)*ireq*fmt.size()*CommonPrms::instance()->NP();
    
    CCIO::cout<<"Reading eigenvector"<< ireq <<" from "<<evecfile<<"\n";
    CCIO::ReadFromDisk<FMT>(evec,evecfile.c_str(),offset);
    CCIO::cout<<"Associated eigenvector loaded.\n";
    return 0;
  }

}

#endif
