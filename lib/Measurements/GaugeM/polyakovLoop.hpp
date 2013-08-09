/*!@file polyakovLoop.hpp
 * @brief calculation of the Polyakov loop
 */
#ifndef POLYAKOVLOOP_INCLUDED
#define POLYAKOVLOOP_INCLUDED

#include "include/field.h"
#include "Geometry/mapping.hpp"
#include "Tools/fieldUtils.hpp"
#include "Tools/sunMatUtils.hpp"
#include <complex>

class PolyakovLoop {
private:
  Communicator* com_;
  int mu_dir_,Nmu_,NPmu_,slsize_;

  template<typename MAT> 
  void calc(GaugeField1D& PL_matrices, const GaugeField& G)const;

public:
  PolyakovLoop(site_dir mu)
    :com_(Communicator::instance()), 
     mu_dir_(mu),
     Nmu_(CommonPrms::instance()->local_size(mu)),
     NPmu_(CommonPrms::instance()->node_num(mu)),
     slsize_(SiteIndex::instance()->slsize(0,mu)){}

  GaugeField1D get_PLField(const GaugeField&)const;

  std::complex<double> calc_SUN(const GaugeField&) const;
  double calc_SUNadj(const GaugeField&) const;
}; 

template<typename MAT>
void PolyakovLoop::calc(GaugeField1D& PL_matrix_field,const GaugeField& G)const{
  using namespace SUNmatUtils;
  using namespace FieldUtils;
  
  int Nin = MAT::size();

  for(int s=0; s<slsize_; ++s){
    MAT pt = MAT(mat(G,SiteMap::shiftSite.xslice(0,s,mu_dir_),mu_dir_));

    for(int xmu=1; xmu<Nmu_; ++xmu)
      pt *= MAT(mat(G,SiteMap::shiftSite.xslice(xmu,s,mu_dir_),mu_dir_));

    
    PL_matrix_field.data.set(std::slice(Nin*s,Nin,1), pt.getva());
  }

  GaugeField1D pc = PL_matrix_field;
  varray_double pfw(pc.size());

  for(int pmu=1; pmu<NPmu_; ++pmu){
    com_->transfer_fw(pfw,pc.data.getva(),mu_dir_);

    for(int s=0; s<slsize_; ++s){
      std::slice ms(Nin*s,Nin,1);
      MAT pt = MAT(PL_matrix_field.data[ms]);

      pt *= MAT(pfw[ms]);
      PL_matrix_field.data.set(ms,pt.getva());
    }
    pc.data = pfw;
  }
}

#endif
