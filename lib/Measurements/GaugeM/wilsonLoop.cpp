/*! @file wilsonLoop.cpp
 *  @brief implementation of wilsonLoop.cpp
 */
#include "wilsonLoop.hpp"

using namespace std;

/*
vector<double> WilsonLoop::calc(const GaugeField& u)const{
  GaugeField1D ut = DirSlice(u,mu_dir_);
  for(int nu=0; nu<NDIM_;++nu){
    if(nu != mu_dir_){
      GaugeField1D v = DirSlice(u,)
      uarm = shift(,mu,Forward());
}
*/

void set_buf(vector<GaugeField1D>& vg,const GaugeField& u, int dir){
  
  for(int mu=0; mu<NDIM_; ++mu){
    if(mu != dir){

      int rlen = CommonPrms::instance()->global_size(mu);
      vector<GaugeField1D> vg(rlen);
      GaugeField1D gt = DirSlice(u,mu);
      vg[0] = gt;
	
      for(int r=0; r<rlen; ++r){
	GaugeField1D sg = shiftField(vg[r],mu,Forward());

	for(int site=0; site<Nvol_; ++site)
	  vg[r+1].data[vg[r+1].format.islice(site)] 
	    = (mat(gt,site)*mat(sg,site)).getva();
      }
      vg.push_back(vg);
    }
  }
}

void shift_fw(GaugeField1D& Gto,const GaugeField1D& Gfrom,int N,int dir){

  Gto = shiftField(Gfrom,dir,Forward());
  for(int s=0; s<N-1; ++s){
    Gbuf_= shiftField(Gto,dir,Forward());
    Gto = Gbuf_;
  }
}

void shift_bk(GaugeField1D& Gto,const GaugeField1D& Gfrom,int N,int dir){

  Gto = shiftField(Gfrom,dir,Backward());
  for(int s=0; s<N-1; ++s){
    Gbuf_= shiftField(Gto,dir,Backward());
    Gto = Gbuf_;
  }
}

