/*
 * @file mdExec_leapfrog.h
 *
 * @brief Declarations of MDexec_leapfrog class and Parameters 
 */
#ifndef MD_LEAPFROG_INCLUDED
#define MD_LEAPFROG_INCLUDED

#ifndef MDEXEC_INCLUDED
#include "mdExec.h"
#endif

#ifndef ACTION_INCLUDED
#include "Action/action.h"
#endif

#include "Measurements/GaugeM/staples.h"
#include "include/pugi_interface.h"

#include<vector>

class Field;
class SUNmat;
class RandNum;
namespace Format{
  class Format_G;
}

typedef std::vector<Action*> ActionLevel;
typedef std::vector<ActionLevel> ActionSet;

struct MDexec_leapfrogParams{
  int Nexp;
  int MDsteps;
  double step_size;

  MDexec_leapfrogParams(XML::node node){
    XML::read(node, "MDsteps", MDsteps);
    XML::read(node, "step_size" , step_size);
    XML::read(node, "exp_approx", Nexp);
  }

  MDexec_leapfrogParams(int Nexp_,int MDsteps_,double step)
    :Nexp(Nexp_),MDsteps(MDsteps_),step_size(step){}
};

class MDexec_leapfrog : public MDexec {
private:
  const MDexec_leapfrogParams Params;
  const std::vector<int> Nrel_;
  const ActionSet as_;
  const Format::Format_G& gf_;
  const Staples* stpl_;
  Field* const CommonField;

  SUNmat u(const Field& g,int site,int dir) const;
  void update_P(Field& P,int lv,double ep) const;
  void update_U(const Field& P,double ep) const;
  void integrator_step(Field& P,int level,std::vector<int>& clock) const;
  
public:
  MDexec_leapfrog(int Nexp, int MDiter, double step,
		  const ActionSet as,
		  const std::vector<int> multipliers,
		  const Format::Format_G& gf,
		  Field* const CommonF)
    :as_(as),gf_(gf),stpl_(new Staples(gf)),
     Params(MDexec_leapfrogParams(Nexp,MDiter,step)),
     Nrel_(multipliers),
     CommonField(CommonF){}
  
  MDexec_leapfrog(XML::node node,
		  const ActionSet as,
		  const std::vector<int> multipliers,
		  const Staples* stpl,
		  const Format::Format_G& gf,
		  Field* const CommonF)
    :as_(as),gf_(gf),stpl_(stpl),
     Params(MDexec_leapfrogParams(node)),
     Nrel_(multipliers),
     CommonField(CommonF){}
  
  void init(std::vector<int>& clock,Field& P,const Field& U,
	    const RandNum& rand)const;
  
  void integrator(Field& P,int level,std::vector<int>& clock) const;
  double calc_H(const Field& P)const;
};

#endif  //MD_LEAPFROG_INCLUDED
