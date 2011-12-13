//----------------------------------------------------------------------
// mdExec.h
//----------------------------------------------------------------------
#ifndef MDEXEC_INCLUDED
#define MDEXEC_INCLUDED

#include<vector>

class Action;
class Field;
class RandNum;

typedef std::vector<Action*> ActionLevel;
typedef std::vector<ActionLevel> ActionSet;

class MDexec{
private:
  virtual void update_P(Field& P,int lv,double ep)const=0;
  virtual void update_U(const Field& P,double ep)const=0;
public:
  virtual ~MDexec(){};
  virtual void init(std::vector<int>& clock,Field& P,const Field& U,
		    const RandNum& rand)const=0;
  virtual double calc_H(const Field& P)const =0;
  virtual void integrator(Field& P,int level,std::vector<int>& clock)const =0;
};
#endif//MDEXEC_INCLUDED
