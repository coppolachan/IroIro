#ifndef MESSAGES_MACRO_HPP_
#define MESSAGES_MACRO_HPP_

#include "Communicator/comm_io.hpp"
#include "include/common_fields.hpp"

#define STRINGIFY(tokens) /****************/                            \
  static_cast<std::ostringstream&>(                                     \
  std::ostringstream().flush() << std::setprecision(15) << tokens ).str()

template< bool Condition, typename T >
class IF_VERB {
public:
  static inline void Msg(std::string M){
    CCIO::cout << M;
  }
  static inline void Monitor(T* obj, GaugeField& f, std::string M){
    obj->monitor_force(f,M);
  }
};

template<typename T>
class IF_VERB< false, T > {
public:
  static inline void Msg(std::string M){}
  static inline void Monitor(T* Obj, GaugeField& f, std::string M){};
};


#define _Message(VerbLevel,Message_)                    \
  IF_VERB<(VERBOSITY>=VerbLevel),NullType>::Msg(STRINGIFY(Message_))

#define _MonitorMsg(VerbLevel, Class_, Field_, String_)                \
  IF_VERB<(VERBOSITY>=VerbLevel), Class_ >::Monitor(this, Field_, String_)


#endif
