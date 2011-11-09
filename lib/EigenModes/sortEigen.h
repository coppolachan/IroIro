//--------------------------------------------------------------------
// sortEigen.h 
//--------------------------------------------------------------------
#ifndef SORTEIGEN_INCLUDED
#define SORTEIGEN_INCLUDED
#include <vector>
#include <utility>

class Field;

class SortEigen{
public:
  virtual void push(std::vector<double>& lmd,int N) const = 0;
  virtual void push(std::vector<double>& lmd,
		    std::vector<Field>& evec,int N) const = 0;
  virtual bool saturated(double lmd, double thrs) const =0;
  virtual ~SortEigen(){}
};


class SortEigen_high :public SortEigen{
private:
  static bool greater_lmd(double left, double right);
  static bool greater_pair(std::pair<double,Field>& left, 
			   std::pair<double,Field>& right);
public:
  void push(std::vector<double>& lmd,int N) const;
  void push(std::vector<double>& lmd,std::vector<Field>& evec,int N) const;
  bool saturated(double lmd, double thrs) const;
};

class SortEigen_low :public SortEigen{
private:
  static bool less_lmd(double left, double right);
  static bool less_pair(std::pair<double,Field>& left, 
			std::pair<double,Field>& right);
public:
  void push(std::vector<double>& lmd,int N) const;
  void push(std::vector<double>& lmd,std::vector<Field>& evec,int N) const;
  bool saturated(double lmd, double thrs) const;
};

#endif
