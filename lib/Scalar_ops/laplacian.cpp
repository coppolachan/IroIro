#include "laplacian.hpp"
#include "Tools/fieldUtils.hpp"
#include "Tools/sunMat.hpp"

/*! 
  Calculates the laplacian operator on a scalar object (no spin d.o.f)
            --                          +
  sol (x) = >   [ u_i(x)*psi(x+i) + u_i(x-i)*psi(x-i) ] -6* psi(x) 
            -- 
             i=(xdir,ydir,zdir)

  That is rewritten in blocks
            --                              
  sol (x) = >   [ chi_1(x) + chi_2(x-i) ] -6* psi(x) 
            -- 
             i=(xdir,ydir,zdir)
  where                                         +
  chi_1 (x) = u_i(x)*psi(x+i),   chi_2 (x) = u_i(x)*psi(x)

  So that we need only two shifts
 */

void Laplacian::setup(){
  /* this class works only with single node in t-direction */
  assert(CommonPrms::instance()->NPEt() == 1); 
  tsl_size_= SiteIndex::instance()->slsize(t_,TDIR);
  
  /* initialization of shiftField3d */
  shiftField3d_.init_maps(map_);
}

const Field Laplacian::mult(const Field& f)const{
  using namespace Mapping;
  using namespace SUNvecUtils;
  using namespace SUNmatUtils;
  using namespace FieldUtils;

  FermionField1sp F(f,sfmt_.Nvol());
  FermionField1sp tmp(sfmt_.Nvol());
  FermionField1sp chi1(sfmt_.Nvol());
  FermionField1sp chi2(sfmt_.Nvol());

  for(int i=0; i<NDIM_-1; ++i){ 
    FermionField1sp psi(shiftField3d_(F,i,Forward()));

    for(int ts=0; ts<tsl_size_; ++ts){
      int site = SiteIndex::instance()->slice_t(t_,ts);
      SUNmat Ui((*u_)[gfmt_.islice(site,i)]);
      FieldUtils::AddVec(chi1,Ui*vec(psi,ts),ts);
      FieldUtils::SetVec(tmp,dag(Ui)*vec(F,ts),ts);
    }   
    chi2 += shiftField3d_(tmp,i,Backward());
  }

  F *= -6.0;
  F += chi1;
  F += chi2;
  
  return F.data;    
}

