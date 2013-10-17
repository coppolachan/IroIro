#include "laplacian.hpp"
#include "Tools/fieldUtils.hpp"
#include "include/commonPrms.hpp"
#include "Tools/sunMat.hpp"
#include "Geometry/shiftField.hpp"

Laplacian::Laplacian(GaugeField *GF):u(GF){
  Mapping::init_shiftField();
};



/*! 
  Calculates the laplacian operator on a scalar object (no spin d.o.f)
            -- 3                             +
  sol (x) = >      [ u_i (x) * psi(x+i) + u_i (x-i) psi(x-i) ] - 6* psi(x) 
            -- i=2 

  That is rewritten in blocks

            -- 3                             
  sol (x) = >      [ chi_1(x) + chi_2(x-i) ] - 6* psi(x) 
            -- i=2 
  
  where
 
  chi_1 (x) = u_i (x) * psi(x+i)

                 +
  chi_2 (x) = u_i (x) psi(x)

  So that we need only two shifts

 */
FermionField1sp Laplacian::apply(const FermionField1sp& F){
  using namespace Mapping;
  using namespace SUNvecUtils;
  using namespace FieldUtils;

  int Nvol = CommonPrms::instance()->Nvol();
  int internal = 2*NC_;
  Communicator* comm = Communicator::instance();
  //int slice_t_size = SiteIndex::instance()->slsize(t, TDIR);

  FermionField1sp output1;
  FermionField1sp chi_1;

  FermionField1sp result;

  SUNvec null;
  null.zero();

  CCIO::cout << " result size:  "<< result.size() << "\n";
  
  for (int mu = 0; mu < NDIM_-1; mu++){ 
    FermionField1sp psi =  shiftField(F, mu, Forward());
   
    for (int site = 0 ; site < Nvol; site++){
      SU3mat loc_m = mat(*u, site, mu);
      SUNvec loc_v = vec(psi, site);
      SUNvec prod = loc_m * loc_v;
      FieldUtils::SetVec(chi_1, prod, site);

      loc_m.dag();
      loc_v = vec(F, site);
      prod = loc_m * loc_v;
      FieldUtils::SetVec(output1, prod, site);
    }   
        
    FermionField1sp chi_2 = shiftField(output1, mu, Backward());
    
    /*
    CCIO::cout << " ---- mu = "<< mu << "\n";
    for (int i = 0; i < F.size(); i++){
      CCIO::cout <<"psi     ["<<i<<"] = "<< psi.data[i] << "\n";
      CCIO::cout <<"chi_2   ["<<i<<"] = "<< chi_2.data[i] << "\n";
      CCIO::cout <<"output2 ["<<i<<"] = "<< output2.data[i] << "\n";
    }
    comm->instance()->sync();
    */
    
    for (int site = 0 ; site < Nvol; site++) {
      for (int i = 0; i < internal; i++) {
	int idx = site*internal+i;
	result.data.add(idx, chi_2.data[idx] + chi_1.data[idx] - 2*F.data[idx]);
      } 
    }
    
    
  }// mu 
  
  return result;
  
};
