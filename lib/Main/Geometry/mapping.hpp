/*!
  @file mapping.hpp
  @brief Declares AutoMap and ShiftField class for the shift
*/
#ifndef MAPPING_H_
#define MAPPING_H_

#include "include/common_fields.hpp"
#include "Communicator/communicator.h"
#include "siteMap.hpp"
#include <vector>
#include <valarray>
//   ____________   ->-    ____________
//  |*|          |   mu   |          |*|
//  |*|  bulk_b_ |        | bulk_t_  |*|
//  |*|__________|        |__________|*|
//   ^ bdry_b_                bdry_t_ ^

namespace Mapping{

  struct EOtag{};
  struct OEtag{};

  struct Forward{};
  struct Backward{};

  class AutoMap{
    int dir_;
    std::vector<int> bdry_t_;
    std::vector<int> bulk_t_;
    std::vector<int> bdry_b_;
    std::vector<int> bulk_b_;

    AutoMap(){}//default constructor

  public:
    // using the reguler site indexing 
    AutoMap(int dir)
      :bdry_t_(SiteMap::shiftSite.bdry_map(dir,Top)),
       bulk_t_(SiteMap::shiftSite.bulk_map(dir,Top)),
       bdry_b_(SiteMap::shiftSite.bdry_map(dir,Btm)),
       bulk_b_(SiteMap::shiftSite.bulk_map(dir,Btm)),
       dir_(dir){}

    template<class DATA,class FORMAT,class TAG>
    GeneralField<DATA,FORMAT,TAG> 
    operator()(const GeneralField<DATA,FORMAT,TAG>& Fin,Forward)const{
      GeneralField<DATA,FORMAT,TAG> Fout(Fin.Nvol());
      std::valarray<double> recv_bdry(bdry_t_.size()*Fin.Nin()*Fin.Nex());
      Communicator::instance()->transfer_fw(recv_bdry,
					    Fin.data[Fin.get_sub(bdry_b_)],
					    dir_);
      Fout.data.set(Fin.get_sub(bdry_t_),recv_bdry);
      Fout.data.set(Fin.get_sub(bulk_t_),Fin.data[Fin.get_sub(bulk_b_)]);
      return Fout;
    }

    template<class DATA,class FORMAT,class TAG>
    GeneralField<DATA,FORMAT,TAG> 
    operator()(const GeneralField<DATA,FORMAT,TAG>& Fin,Backward) const{
      GeneralField<DATA,FORMAT,TAG> Fout(Fin.Nvol());
      
      std::valarray<double> recv_bdry(bdry_b_.size()*Fin.Nin()*Fin.Nex());
      Communicator::instance()->transfer_bk(recv_bdry,
					    Fin.data[Fin.get_sub(bdry_t_)],
					    dir_);
      Fout.data.set(Fin.get_sub(bdry_b_),recv_bdry);
      Fout.data.set(Fin.get_sub(bulk_b_),Fin.data[Fin.get_sub(bulk_t_)]);
      return Fout;
    }
  };

  class AutoMap_EvenOdd{
    int dir_;
    std::vector<int> recv_bdry_t_;
    std::vector<int> recv_bulk_t_;
    std::vector<int> recv_bdry_b_;
    std::vector<int> recv_bulk_b_;

    std::vector<int> send_bdry_t_;
    std::vector<int> send_bulk_t_;
    std::vector<int> send_bdry_b_;
    std::vector<int> send_bulk_b_;

    AutoMap_EvenOdd(){}//default constructor

  public:
    // using the e/o site indexing (e<-o)
    AutoMap_EvenOdd(int dir,EOtag)
      :send_bdry_t_(SiteMap::shiftSite_oe.bdry_map(dir,Top)),
       send_bulk_t_(SiteMap::shiftSite_oe.bulk_map(dir,Top)),
       send_bdry_b_(SiteMap::shiftSite_oe.bdry_map(dir,Btm)),
       send_bulk_b_(SiteMap::shiftSite_oe.bulk_map(dir,Btm)),
       recv_bdry_t_(SiteMap::shiftSite_eo.bdry_map(dir,Top)),
       recv_bulk_t_(SiteMap::shiftSite_eo.bulk_map(dir,Top)),
       recv_bdry_b_(SiteMap::shiftSite_eo.bdry_map(dir,Btm)),
       recv_bulk_b_(SiteMap::shiftSite_eo.bulk_map(dir,Btm)),
       dir_(dir){}

    // using the e/o site indexing (o<-e)
    AutoMap_EvenOdd(int dir,OEtag)
      :send_bdry_t_(SiteMap::shiftSite_eo.bdry_map(dir,Top)),
       send_bulk_t_(SiteMap::shiftSite_eo.bulk_map(dir,Top)),
       send_bdry_b_(SiteMap::shiftSite_eo.bdry_map(dir,Btm)),
       send_bulk_b_(SiteMap::shiftSite_eo.bulk_map(dir,Btm)),
       recv_bdry_t_(SiteMap::shiftSite_oe.bdry_map(dir,Top)),
       recv_bulk_t_(SiteMap::shiftSite_oe.bulk_map(dir,Top)),
       recv_bdry_b_(SiteMap::shiftSite_oe.bdry_map(dir,Btm)),
       recv_bulk_b_(SiteMap::shiftSite_oe.bulk_map(dir,Btm)),
       dir_(dir){}
    
    template<class DATA,class FORMAT,class TAG>
    GeneralField<DATA,FORMAT,TAG> 
    operator()(const GeneralField<DATA,FORMAT,TAG>& Fin,Forward)const{
      GeneralField<DATA,FORMAT,TAG> Fout(Fin.Nvol());
      std::valarray<double> recv_bdry(recv_bdry_t_.size()*Fin.Nin()*Fin.Nex());
      Communicator::instance()->transfer_fw(recv_bdry,
					    Fin.data[Fin.get_sub(send_bdry_b_)],
					    dir_);
      Fout.data.set(Fin.get_sub(recv_bdry_t_),recv_bdry);
      Fout.data.set(Fin.get_sub(recv_bulk_t_),
		    Fin.data[Fin.get_sub(send_bulk_b_)]);
      return Fout;
    }

    template<class DATA,class FORMAT,class TAG>
    GeneralField<DATA,FORMAT,TAG> 
    operator()(const GeneralField<DATA,FORMAT,TAG>& Fin,Backward) const{
      GeneralField<DATA,FORMAT,TAG> Fout(Fin.Nvol());
      
      std::valarray<double> recv_bdry(recv_bdry_b_.size()*Fin.Nin()*Fin.Nex());
      Communicator::instance()->transfer_bk(recv_bdry,
					    Fin.data[Fin.get_sub(send_bdry_t_)],
					    dir_);
      Fout.data.set(Fin.get_sub(recv_bdry_b_),recv_bdry);
      Fout.data.set(Fin.get_sub(recv_bulk_b_),
		    Fin.data[Fin.get_sub(send_bulk_t_)]);
      return Fout;
    }
  };

  class ShiftField{
    std::vector<AutoMap> maps_;
  public:
    template<typename DATA,typename FORMAT,typename TAG,typename FB>
    GeneralField<DATA,FORMAT,TAG>
    operator()(const GeneralField<DATA,FORMAT,TAG>& F,int dir, FB fb)const{
      return maps_[dir](F,fb);}
    void init_maps();
  };

  class ShiftField_eo{
    std::vector<AutoMap_EvenOdd> maps_;
  public:
    template<typename DATA,typename FORMAT,typename TAG,typename FB>
    GeneralField<DATA,FORMAT,TAG>
    operator()(const GeneralField<DATA,FORMAT,TAG>& F,int dir, FB fb)const{
      return maps_[dir](F,fb);}
    void init_maps();
  };

  class ShiftField_oe{
    std::vector<AutoMap_EvenOdd> maps_;
  public:
    template<typename DATA,typename FORMAT,typename TAG,typename FB>
    GeneralField<DATA,FORMAT,TAG>
    operator()(const GeneralField<DATA,FORMAT,TAG>& F,int dir, FB fb)const{
      return maps_[dir](F,fb);}
    void init_maps();
  };

  // declaration of a global object
  extern ShiftField shiftField;
  void init_shiftField();
  
  extern ShiftField_eo shiftField_eo;
  extern ShiftField_oe shiftField_oe;
  void init_shiftField_EvenOdd();
}

#endif
