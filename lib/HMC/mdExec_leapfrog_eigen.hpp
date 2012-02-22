//--------------------------------------------------------------------
// mdExec_leapfrog_eigen.hpp
//--------------------------------------------------------------------
#ifndef MD_LEAPFROG_EIGEN_INCLUDED
#define MD_LEAPFROG_EIGEN_INCLUDED

#ifndef COMMONPRMS_INCLUDED
#include "include/commonPrms.h"
#endif

#ifndef MDEXEC_INCLUDED
#include "mdExec.h"
#endif

#ifndef ACTION_INCLUDED
#include "Action/action.h"
#endif

#include "hmcPrms.h"
#include "include/typeDefs.h"

#include<vector>

class Field;
class Staples;
class SUNmat;
namespace Format{
  class Format_G;
}

typedef std::vector<Action*> ActionLevel;
typedef std::vector<ActionLevel> ActionSet;
typedef std::vector<EigenProc_Zolotarev*> EigenLevel;
typedef std::vector<EigenLevel> EigenSet;

//typedef Fopr_signH_Zolotarevn> Fopr_signH;
//typedef Fopr_signH::EigenData EigDat;

struct MDexec_leapfrog_EigenParams {
  int Nexp;
  int MDsteps;
  double step_size;

  MDexec_leapfrog_EigenParams(XML::node node) {
    XML::read(node, "MDsteps", MDsteps);
    XML::read(node, "step_size" , step_size);
    XML::read(node, "exp_approx", Nexp);
  }

  MDexec_leapfrog_EigenParams(int Nexp_, 
			      int MDsteps_,
			      double step):
    Nexp(Nexp_),MDsteps(MDsteps_),step_size(step){};

};


class MDexec_leapfrog_eigen : public MDexec {
private:
  const MDexec_leapfrog_EigenParams Params;
  const ActionSet& as_;
  const EigenSet& es_; 
  const Format::Format_G& gf_;
  const Staples* stpl_;         
  Field* const CommonField;

  std::vector<int> Nrel_;

  SUNmat u(const Field& g,int site,int dir) const;

  void update_P(Field& P,int lv,double ep) const;
  void update_U(const Field& P,double ep) const;
  
public:
  MDexec_leapfrog_eigen(const ActionSet& as,
  			const EigenSet& es,
  			const Format::Format_G& gf,
  			const Staples& stpl)
    :as_(as),es_(es),gf_(gf),stpl_(&stpl),
     CommonField(new Field(gf.size())),
     Params(MDexec_leapfrog_EigenParams(HMCprms::instance()->Nexp(),
  					HMCprms::instance()->MDsteps(),
  					HMCprms::instance()->epsilon())),
     Nrel_(   HMCprms::instance()->Nrel()){};

  MDexec_leapfrog_eigen(XML::node node,
			const ActionSet& as,
			const EigenSet& es,
			const std::vector<int> multipliers,
			const Format::Format_G& gf,
			Field* const CommonF)
    :as_(as),es_(es),gf_(gf),stpl_(new Staples(gf)),
     Params(MDexec_leapfrog_EigenParams(node)),
     Nrel_(multipliers),
     CommonField(CommonF){};

  void init(std::vector<int>& clock,
	    Field& P,
	    const Field& U,
	    const RandNum& rand) const;

  void integrator(Field& P,
		  int level,
		  std::vector<int>& clock) const;
  double calc_H()const;

};
#endif
