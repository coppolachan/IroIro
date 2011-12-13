//------------------------------------------------------------------------
// typeDefs.h
//-----------------------------------------------------------------------
#ifndef TYPEDEFS_INCLUDED
#define TYPEDEFS_INCLUDED

class Action;

#ifndef SOLVER_CG_INCLUDED
#include "Solver/solver_CG.h"
#endif
#ifndef ACTION_GAUGE_INCLUDED
#include "Action/action_gauge.h"
#endif
#ifndef DIRAC_WILSON_INCLUDED
#include "Dirac_ops/dirac_wilson.h"
#endif
#ifndef ACTION_NF2_INCLUDED
#include "Action/action_Nf2.h"
#endif
#ifndef ACTION_NF2_RATIO_INCLUDED
#include "Action/action_Nf2_ratio.h"
#endif
#ifndef FOPR_SIGNH_ZOLOTAREV_INCLUDED
#include "include/fopr_signH_Zolotarev.h"
#endif
#ifndef EIGENPROCESS_INCLUDED
#include "EigenModes/eigenProc_Zolotarev.h"
#endif
#include <vector>


typedef ActionGaugeWilson Sg;

typedef std::vector<Action*> ActionLevel;
typedef std::vector<ActionLevel> ActionSet;
typedef std::vector<EigenProc_Zolotarev*> EigenLevel;
typedef std::vector<EigenLevel> EigenSet;
/*
typedef Fopr_signH_Zolotarev<Dirac_Wilson> Fopr_signH;
typedef Fopr_signH::EigenData EigDat;
typedef Fopr_signH::EigenPrms EigPms;
*/
#endif
