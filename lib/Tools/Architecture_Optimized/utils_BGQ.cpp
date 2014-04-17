/*!--------------------------------------------------------------------------
 * @file utils_BGQ.cpp

 * Time-stamp: <2014-04-17 12:50:14 neo>
 *-------------------------------------------------------------------------*/

#include "utils_BGQ.hpp"
#include "Communicator/communicator.hpp"
#include <spi/include/kernel/memory.h> // memory tools

double omp_norm(Spinor* pointer,int is,int ns,int nid,int tid){
  double tSum;
  BGWilsonLA_Norm(&tSum,pointer+is,ns);
  tSum = BGQThread_GatherDouble(tSum,0,tid,nid);
  if(tid == 0) tSum = Communicator::instance()->reduce_sum(tSum);
  tSum = BGQThread_ScatterDouble(tSum,0,tid,nid);
  return sqrt(tSum);
}

void MonitorBGQMemory(){
  uint64_t shared, persist, heapavail, estheapavail, stackavail, stack, heap, guard, mmap;
  Communicator::instance()->sync();
  CCIO::cout << "Memory check summary after free of BFM_Wrapper\n";

  Kernel_GetMemorySize(KERNEL_MEMSIZE_SHARED, &shared);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_PERSIST, &persist);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPAVAIL, &heapavail);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_ESTHEAPAVAIL, &estheapavail);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACKAVAIL, &stackavail);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACK, &stack);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_GUARD, &guard);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_MMAP, &mmap);
  
  if (Communicator::instance()->primaryNode()){
    printf("Allocated heap: %.2f MB, avail. heap: %.2f MB, estimated avail. heap: %.2f MB\n", 
	   (double)heap/(1024*1024),(double)heapavail/(1024*1024),(double)estheapavail/(1024*1024));
    printf("Allocated stack: %.2f MB, avail. stack: %.2f MB\n", 
	   (double)stack/(1024*1024), (double)stackavail/(1024*1024));
    printf("Memory: shared: %.2f MB, persist: %.2f MB, guard: %.2f MB, mmap: %.2f MB\n", 
	   (double)shared/(1024*1024), (double)persist/(1024*1024), (double)guard/(1024*1024), (double)mmap/(1024*1024));
  }
  Communicator::instance()->sync();

}
