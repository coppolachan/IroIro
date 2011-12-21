//--------------------------------------------------------------------
// hmcGeneral.h
//--------------------------------------------------------------------
#ifndef HMCEXEC_INCLUDED
#define HMCEXEC_INCLUDED

#ifndef HMCPRMS_INCLUDED
//#include "hmcPrms.h"
#endif
#ifndef COMMUNICATOR_INCLUDED
#include "Communicator/communicator.h"
#endif
#include "Tools/randNum_Factory.h"
#include "HMC/mdExec_Factory.hpp"
#include <vector>

#ifndef PUGI_INTERFACE_H_
#include "include/pugi_interface.h"
#endif

class Field;
class MDexec;
class RandNum;

struct HMCGeneralParams {
  int Nsweeps;
  int ThermalizationSteps;

  HMCGeneralParams(pugi::xml_node node) {
    XML::read(node, "Nsweeps", Nsweeps, MANDATORY);
    XML::read(node, "Thermalization", ThermalizationSteps, MANDATORY);
  }
};

class HMCgeneral{
private:
  MDexec* md_;
  const RandNum* rand_;
  const HMCGeneralParams Params;
  int delete_me_;
  bool metropolis_test(const double Hdiff)const;
  double evolve_step(Field&)const;

public:
  HMCgeneral(pugi::xml_node node, 
	     MDexec& md_exec, 
	     const RandNum& rand_num)
    :md_(&md_exec),
     rand_(&rand_num),
     Params(HMCGeneralParams(node)),
     delete_me_(0){}
  
  HMCgeneral(pugi::xml_node node, 
	     MDexec& md_exec)
    :md_(&md_exec),
     rand_(RNG_Env::RNG->getRandomNumGenerator()),
     Params(HMCGeneralParams(node)),
     delete_me_(1){}
  
  HMCgeneral(pugi::xml_node node)
    :md_(Integrators::Integr->getMDIntegrator()),
     rand_(RNG_Env::RNG->getRandomNumGenerator()),
     Params(HMCGeneralParams(node)),
     delete_me_(2){}
  
  ~HMCgeneral(){
    if (delete_me_ >= 1) delete rand_;
    if (delete_me_ == 2) delete md_;
  }
  void evolve(Field& U)const;
};

#endif
