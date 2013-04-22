/*!
  @file dcmt_wrapper.cpp
  @brief Class to handle the dcmt C library by Matsumoto and Nishimura

 Reference:
 Makoto Matsumoto and Takuji Nishimura,
 "Dynamic Creation of Pseudorandom Number Generators",
 Monte Carlo and Quasi-Monte Carlo Methods 1998,
 Springer, 2000, pp 56--69. 

 Open code here:
 http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/DC/dc.html
*/

#include <stdint.h>

#include "dcmt_wrapper.hpp"
#include "include/messages_macros.hpp"
#include "include/errors.hpp"
#define BITS_FACTOR32   2.32830643708079737543e-10
#define BITS_FACTOR32_1 2.3283064365386962890625e-10
#define BITS_FACTOR53   1.11022302462515654042e-16
#define BITS_FACTOR64   5.42101086242752217033113759205528043e-20

int RandNum_DCMT::init(int ID, int genSeed, int seed){
  // genSeed 16 bit integer, seed 32 bit
  // Get the independent MT
  _Message(DEBUG_VERB_LEVEL, "Initialization of DC Mersenne Twister\n");
  mts = get_mt_parameter_id_st(w, p, ID, genSeed);
  if (mts == NULL) {
    Errors::AllocationErr("RandNum_DCMT::init - mts structure not allocated");
  }
  // Initialize RNG mts with 32-bit seed;
  sgenrand_mt(seed, mts);
  return 0;
}

mt_struct *RandNum_DCMT::alloc_mt_struct(int n){
  mt_struct *mts;
  
  mts = (mt_struct*)malloc(sizeof(mt_struct));
  if (NULL == mts) return NULL;
  mts->state = (uint32_t*)malloc(n*sizeof(uint32_t));
  if (NULL == mts->state) {
    free(mts);
    return NULL;
  }
  
  return mts;
}


//generates a random number on [0,0x7fffffff]
long RandNum_DCMT::rand_int31()const{ return (long)(genrand_mt(mts)>>1);}

//generates a random number on [0,1] double number
double RandNum_DCMT::rand_double1()const{
  return static_cast<double>(genrand_mt(mts)*BITS_FACTOR32);
}

//generates a random number on [0,1) double number
double RandNum_DCMT::rand_double2()const{
  return static_cast<double>(genrand_mt(mts)*BITS_FACTOR32_1);
}

// generates a random number on (0,1) double number
double RandNum_DCMT::rand_double3()const{
  return static_cast<double>((genrand_mt(mts)+0.5)*BITS_FACTOR32_1);
}

// generates a random number on [0,1) with 53-bit resolution
inline double RandNum_DCMT::rand_res53()const{
  return ((genrand_mt(mts)>>5)*67108864.0+(genrand_mt(mts)>>6))*BITS_FACTOR53;
}

// generates a 64 bit random number on [0,1]
inline double RandNum_DCMT::rand_64bit()const{
  register uint64_t low  = genrand_mt(mts);
  register uint64_t high = genrand_mt(mts);
  register uint64_t random = low | (high << 32);
  return random * BITS_FACTOR64;

}

double RandNum_DCMT::do_rand()const {
  return rand_res53();
}

double RandNum_DCMT::do_rand_closed()const {
  return rand_64bit();
}

RandNum_DCMT::~RandNum_DCMT(){
  if (mts != NULL)
    free_mt_struct(mts);
}


void RandNum_DCMT::saveSeed(const std::string& file) const{
  CCIO::cout << "Saving DC Mersenne Twister in file ["
	     << file << "]\n";
  int state_size = (p/w+1);
  std::ofstream writer;
  for (int n = 0; n < Communicator::instance()->size(); n++) {
    if (n == Communicator::instance()->id()) {
      if (n == 0){
	writer.open(file.c_str()); //open new file
	//Write header with some useful info
	writer << CommonPrms::instance()->node_num(0) << " "
	       << CommonPrms::instance()->node_num(1) << " "
	       << CommonPrms::instance()->node_num(2) << " "
	       << CommonPrms::instance()->node_num(3) << std::endl;
      } else {
	//append data 
	writer.open(file.c_str(),std::ofstream::app); 
      }
      
      writer << mts->aaa << std::endl;
      writer << mts->mm << " " << mts->nn << " "
	     << mts->rr << " "<< mts->ww << std::endl;
      writer << mts->wmask << " " << mts->umask << " "
	     << mts->lmask << std::endl;     
      writer << mts->shift0 << " " << mts->shift1 << " "
	     << mts->shiftB << " "<< mts->shiftC << std::endl;      
      writer << mts->maskB << " " << mts->maskC << std::endl; 
      writer << mts->i << std::endl;
      for (int i = 0; i < state_size; i++) {
	writer << mts->state[i] << std::endl;
      }
      writer << "\n";
      writer.close();
    }
    Communicator::instance()->sync();
  }
};

void RandNum_DCMT::loadSeed(const std::string& file){
  CCIO::cout << "Loading DC Mersenne Twister structures from file ["
	     << file << "]\n";
  int state_size = (p/w+1);
  int np0, np1, np2, np3;
  std::ifstream reader(file.c_str());
  //check that the processes match
  //Read header with some useful info
  reader >> np0 >> np1 >> np2 >> np3;
  if (Communicator::instance()->primaryNode()){
    if ((np0 == CommonPrms::instance()->node_num(0)) &&
	(np1 == CommonPrms::instance()->node_num(1)) &&
	(np2 == CommonPrms::instance()->node_num(2)) &&
	(np3 == CommonPrms::instance()->node_num(3))){
      std::cout << "Number and distribution of nodes matches\n";
    } else {
      std::cout << "[Error] The Number and distribution of nodes does not match.";
      exit(1);
    }
  }
  Communicator::instance()->sync();
  //allocate space for the structure.
  mts = alloc_mt_struct(state_size);
  // load and store data
  // each process reads the same data but stores only the ones referring to its own process.
  mt_struct *temp_mt;
  temp_mt = alloc_mt_struct(state_size);
  for (int n = 0; n < Communicator::instance()->size(); n++) {
    reader >> temp_mt->aaa ;
    reader >> temp_mt->mm >> temp_mt->nn >> temp_mt->rr >> temp_mt->ww;
    reader >> temp_mt->wmask >> temp_mt->umask >> temp_mt->lmask;
    reader >> temp_mt->shift0 >> temp_mt->shift1 >> temp_mt->shiftB >> temp_mt->shiftC;
    reader >> temp_mt->maskB >> temp_mt->maskC;
    reader >> temp_mt->i ;
    for (int i = 0; i < state_size; i++) {
      reader >> temp_mt->state[i];
    }
    if (n == Communicator::instance()->id()){
      // deep copy to the real mts structure
      mts->aaa    = temp_mt->aaa;
      mts->mm     = temp_mt->mm;
      mts->nn     = temp_mt->nn;
      mts->rr     = temp_mt->rr;
      mts->ww     = temp_mt->ww;
      mts->wmask  = temp_mt->wmask;
      mts->umask  = temp_mt->umask;
      mts->lmask  = temp_mt->lmask;
      mts->shift0 = temp_mt->shift0;
      mts->shift1 = temp_mt->shift1;
      mts->shiftB = temp_mt->shiftB;
      mts->shiftC = temp_mt->shiftC;
      mts->maskB  = temp_mt->maskB;
      mts->maskC  = temp_mt->maskC;
      mts->i      = temp_mt->i;
      for (int i = 0; i < state_size; i++) {
	mts->state[i] = temp_mt->state[i];
      }
      break;
    }
  }
  Communicator::instance()->sync();
  free_mt_struct(temp_mt);
};
