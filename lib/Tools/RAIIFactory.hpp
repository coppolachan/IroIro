/*!
 * @file RAIIFactory.hpp
 *
 * @brief Implements Resource Acquisition Is Initialization behaviour
 *
 */
#ifndef RAII_FACTORY_H_
#define RAII_FACTORY_H_

#include <vector>


/*!
 * @brief Class for RAII(Resource Acquisition Is Initialization)
 *
 */
template <typename T>
class RaiiFactoryObj
{
public:
  RaiiFactoryObj(){};

  ~RaiiFactoryObj()
  {
    delete pointer;
  }
  
  T* get(){
    return pointer;
  }

  T* save (T* p)
  {
    pointer = p;//save the object pointer
    return p;
  }

private:
  T* pointer;

  // non-copyable

  RaiiFactoryObj (const RaiiFactoryObj&);
  RaiiFactoryObj& operator= (const RaiiFactoryObj&);
};



/*!
 * @brief Class for RAII(Resource Acquisition Is Initialization) vector
 *
 */
template <typename T>
class RaiiFactoryObjVec
{
public:
  RaiiFactoryObjVec(){};

  ~RaiiFactoryObjVec()
  {
    while (!pointer.empty())
    {
      T* p = pointer.back();
      pointer.pop_back();
      delete p;
    }
  }
  
  T* get(const int i){
    return pointer[i];
  }

  T* save (T* p)
  {
    pointer.push_back(p);//save the object pointer
    return p;
  }

private:
  std::vector<T*> pointer;

  // non-copyable

  RaiiFactoryObjVec (const RaiiFactoryObjVec&);
  RaiiFactoryObjVec& operator= (const RaiiFactoryObjVec&);
};


#endif
