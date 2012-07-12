// General threaded allocator for BGQ
// Based on Peter Boyle bfm code

#include <spi/include/l2/atomic.h>
#include <spi/include/kernel/location.h>
#include <pthread.h>

#define L1D_CACHE_LINE_SIZE 64

typedef struct {
  volatile __attribute__((aligned(L1D_CACHE_LINE_SIZE)))
    uint64_t start;  /*!< Thread count at start of current round. */
  volatile __attribute__((aligned(L1D_CACHE_LINE_SIZE)))
    uint64_t count;  /*!< Current thread count. */
} L2_Barrier_t;

L2_Barrier_t init_b = {0,0};

/*Locate in atomic op space*/
L2_Barrier_t  b ;

int L2_BarrierWithTicket(L2_Barrier_t *b, int numthreads);

int TLBmapped=0;

void SpiL2AtomicInit(void)
{
    Kernel_L2AtomicsAllocate(&b,sizeof(L2_Barrier_t));
}

int L2_BarrierWithTicket(L2_Barrier_t *b, int numthreads)
{
  if( !TLBmapped)  SpiL2AtomicInit();
  TLBmapped=1;

  uint64_t start = b->start;
  ppc_msync();  // make sure we pick up the correct start for this round
  uint64_t count = L2_AtomicLoadIncrement(&b->count);

  uint64_t target = start + numthreads;
  uint64_t current = count + 1;

  if (current == target) {
    b->start = current;  // advance to next round
  } else {
    int count = 0;
    while (b->start < current) {
      count++;
      if ( count > 10*1024 ) {
        count =0; 
        pthread_yield();
      }
    }  // wait for advance to next round
  }

  return (int) (count - start);
}

int thread_barrier(void) 
{
  int me ;
  me = L2_BarrierWithTicket(&b,nthread);

  // if ( PhysicalThreadMapping ) me = bgq_hwid();

  //  Delay(me*20);

  return me;
}

void  thread_bcast(int me,void * val) 
{ 

  if (me == 0) bcast_ptr = val;
  thread_barrier();
  val = bcast_ptr;
  thread_barrier();
  return val; 
};


////////////////////////////////////////////////////////////


void * threadedAlloc (int block_size)
{
  void* ret;
  int me = thread_barrier();
  if ( me == 0 ) {
    ret = this->allocBlock(block_size);
  }
  ret = thread_bcast(me,ret);
  thread_barrier();
  return ret;
}

void threadedFree (void* handle)
{
  int me = thread_barrier();
  if ( me == 0 ) { 
     freeBlock(ptr);
  }
  thread_barrier();
}

void* allocBlock (int mem_size)
{
  void *ptr = (double *)malloc(mem_size);
  if ( ptr == NULL ) {
    std::cout << "bad bfm_alloc\n";
    exit(-1);
  }
  return ptr;
}
void freeBlock    (void* handle)
{
  free(handle);
}
