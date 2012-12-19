/*!--------------------------------------------------------------------------
 * @file boundaryCond.cpp
 * @brief implementation of the BoundaryCond class
 *-------------------------------------------------------------------------*/
#include "boundaryCond.hpp"
#include "Main/Geometry/siteMap.hpp"
#include "include/numerical_const.hpp"
#include "Tools/fieldUtils.hpp"
#include "Tools/sunMatUtils.hpp"
#include <cstring>

using namespace std;
using namespace FieldUtils;
using namespace SUNmatUtils;

////// BC creator ///////
BoundaryCond* createBC(const XML::node& node){
  if(node !=NULL){
    const char* bcname = node.attribute("name").value();

    if (!strcmp(bcname,"")) {
      std::cerr<< "No name provided for boundary condition. Request by <"
	       << node.name() << ">\n";
      abort();
    }
    if(!strcmp(bcname,"periodic")) 
      return new BoundaryCond_periodic(node);
    if(!strcmp(bcname,"anti_periodic")) 
      return new BoundaryCond_antiPeriodic(node);
    if(!strcmp(bcname,"U1phase")) 
      return new BoundaryCond_U1phase(node);
    if(!strcmp(bcname,"SUNdiag")) 
      return new BoundaryCond_SUNmat(node);
  }
}

////////////////// member functions ///////////////////////
////// anti-periodic BC //////
BoundaryCond_antiPeriodic::BoundaryCond_antiPeriodic(const XML::node& bcnode){
  const char* dir_name = bcnode.attribute("dir").value();
  if(     !strcmp(dir_name,"X")) dir_= XDIR;
  else if(!strcmp(dir_name,"Y")) dir_= YDIR;
  else if(!strcmp(dir_name,"Z")) dir_= ZDIR;
  else if(!strcmp(dir_name,"T")) dir_= TDIR;
  else {
    CCIO::cout<<"No valid direction available"<<endl;
    abort();
  }
}

void BoundaryCond_antiPeriodic::apply_bc(GaugeField& u)const{
  if(Communicator::ipe(dir_)==CommonPrms::NPE(dir_)-1){

    int Nbd = SiteIndex::instance()->Bdir(dir_);
    int slsize = SiteIndex::instance()->slsize(Nbd,dir_);

    for(int n=0; n<slsize; ++n){
      int site = SiteMap::shiftSite.xslice(Nbd,n,dir_);
      u.data.set(u.format.islice(site,dir_),-mat(u,site,dir_).getva());
    }
  }
}

void BoundaryCond_antiPeriodic::apply_bc(GaugeField& ue,GaugeField& uo)const{
  if(Communicator::ipe(dir_)==CommonPrms::NPE(dir_)-1){

    int Nbd = SiteIndex_EvenOdd::instance()->Bdir(dir_);

    int slsize = SiteMap::shiftSite_eo.slice_size(Nbd,dir_);
    for(int n=0; n<slsize; ++n){
      int hs = SiteMap::shiftSite_eo.xslice(Nbd,n,dir_);
      ue.data.set(ue.format.islice(hs,dir_),-mat(ue,hs,dir_).getva());
    }
    slsize = SiteMap::shiftSite_oe.slice_size(Nbd,dir_);
    for(int n=0; n<slsize; ++n){
      int hs = SiteMap::shiftSite_oe.xslice(Nbd,n,dir_);
      uo.data.set(uo.format.islice(hs,dir_),-mat(uo,hs,dir_).getva());
    }
  }
}

/////// U(1) phase BC /////////
BoundaryCond_U1phase::BoundaryCond_U1phase(const XML::node& bcnode)
  :bc_(1.0,0.0){
  const char* dir_name = bcnode.attribute("dir").value();
  if(     !strcmp(dir_name,"X")) dir_= XDIR;
  else if(!strcmp(dir_name,"Y")) dir_= YDIR;
  else if(!strcmp(dir_name,"Z")) dir_= ZDIR;
  else if(!strcmp(dir_name,"T")) dir_= TDIR;
  else {
    CCIO::cout<<"No valid direction specified"<<endl;
    abort();
  }
  int Nc = CommonPrms::instance()->Nc();
  double theta; /*! @brief phase */
  XML::read(bcnode,"phase",theta,MANDATORY);
  bc_= complex<double>(cos(theta),sin(theta));
}

void BoundaryCond_U1phase::apply_bc(GaugeField& u)const{
  if(Communicator::ipe(dir_)==CommonPrms::NPE(dir_)-1){

    int Nbd = SiteIndex::instance()->Bdir(dir_);
    int slsize = SiteIndex::instance()->slsize(Nbd,dir_);

    for(int n=0; n<slsize; ++n){
      int site = SiteMap::shiftSite.xslice(Nbd,n,dir_);
      SetMat(u,mat(u,site,dir_)*bc_,site,dir_);
    }
  }
}

void BoundaryCond_U1phase::apply_bc(GaugeField& ue,GaugeField& uo)const{
  if(Communicator::ipe(dir_)==CommonPrms::NPE(dir_)-1){

    int Nbd = SiteIndex_EvenOdd::instance()->Bdir(dir_);

    int slsize = SiteMap::shiftSite_eo.slice_size(Nbd,dir_);
    for(int n=0; n<slsize; ++n){
      int hs = SiteMap::shiftSite_eo.xslice(Nbd,n,dir_);
      SetMat(ue,mat(ue,hs,dir_)*bc_,hs,dir_);
    }
    slsize = SiteMap::shiftSite_oe.slice_size(Nbd,dir_);
    for(int n=0; n<slsize; ++n){
      int hs = SiteMap::shiftSite_oe.xslice(Nbd,n,dir_);
      SetMat(uo,mat(uo,hs,dir_)*bc_,hs,dir_);
    }
  }
}

/////// SU(N) matrix BC /////////
/*!@brief this constructor is compatible only with diagonal matrix*/
BoundaryCond_SUNmat::BoundaryCond_SUNmat(const XML::node& bcnode)
 :bm_(SUNmatUtils::unity()){
  const char* dir_name = bcnode.attribute("dir").value();
  if(     !strcmp(dir_name,"X")) dir_= XDIR;
  else if(!strcmp(dir_name,"Y")) dir_= YDIR;
  else if(!strcmp(dir_name,"Z")) dir_= ZDIR;
  else if(!strcmp(dir_name,"T")) dir_= TDIR;
  else {
    CCIO::cout<<"No valid direction specified"<<endl;
    abort();
  }
  int Nc = CommonPrms::instance()->Nc();
  vector<double> theta(Nc-1); /*! @brief phase */
  
  XML::read_array(bcnode,"phase",theta,MANDATORY);
  theta.push_back(2*PI);      /*! @brief now theta.size() = Nc */
  for(int c=0; c<Nc-1; ++c) theta[Nc-1] -= theta[c];
  for(int c=0; c<Nc; ++c) bm_.set(c,c,cos(theta[c]),sin(theta[c]));
}

void BoundaryCond_SUNmat::apply_bc(GaugeField& u)const{
  if(Communicator::ipe(dir_)==CommonPrms::NPE(dir_)-1){

    int Nbd = SiteIndex::instance()->Bdir(dir_);
    int slsize = SiteIndex::instance()->slsize(Nbd,dir_);

    for(int n=0; n<slsize; ++n){
      int site = SiteMap::shiftSite.xslice(Nbd,n,dir_);
      SetMat(u,bm_*mat(u,site,dir_),site,dir_);
    }
  }
}

void BoundaryCond_SUNmat::apply_bc(GaugeField& ue,GaugeField& uo)const{
  if(Communicator::ipe(dir_)==CommonPrms::NPE(dir_)-1){

    int Nbd = SiteIndex_EvenOdd::instance()->Bdir(dir_);

    int slsize = SiteMap::shiftSite_eo.slice_size(Nbd,dir_);
    for(int n=0; n<slsize; ++n){
      int hs = SiteMap::shiftSite_eo.xslice(Nbd,n,dir_);
      SetMat(ue,bm_*mat(ue,hs,dir_),hs,dir_);
    }
    slsize = SiteMap::shiftSite_oe.slice_size(Nbd,dir_);
    for(int n=0; n<slsize; ++n){
      int hs = SiteMap::shiftSite_oe.xslice(Nbd,n,dir_);
      SetMat(uo,bm_*mat(uo,hs,dir_),hs,dir_);
    }
  }
}

