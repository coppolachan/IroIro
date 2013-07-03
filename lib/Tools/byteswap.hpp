/*! 
 * @file byteswap.hpp 
 * @brief Definition of function for trasforming a vector of big-endian doubles to litte-endian format
 *
 * Portable solution, performance is not necessary
 */

/*!
 * @brief Tool to trasform an array of doubles from big-endian to little-endian format
 */

typedef unsigned int uint32_t;



inline void byte_swap_double(void *ptr, size_t nmemb)
{
  unsigned int j;
  char char_in[8];              /* characters used in byte swapping */
                                                                                
  char *in_ptr;
  double *double_ptr;           /* Pointer used in the double routines */
  
  for(j = 0, double_ptr = (double *) ptr;
      j < nmemb;
      j++, double_ptr++){
    
    
    in_ptr = (char *) double_ptr; /* Set the character pointer to
                                       point to the start of the double */
    
    /*
     *  Assign all the byte variables to a character
     */
    
    char_in[0] = in_ptr[0];
    char_in[1] = in_ptr[1];
    char_in[2] = in_ptr[2];
    char_in[3] = in_ptr[3];
    char_in[4] = in_ptr[4];
    char_in[5] = in_ptr[5];
    char_in[6] = in_ptr[6];
    char_in[7] = in_ptr[7];
    /*
     *  Now just swap the order
     */
    
    in_ptr[0] = char_in[7];
    in_ptr[1] = char_in[6];
    in_ptr[2] = char_in[5];
    in_ptr[3] = char_in[4];
    in_ptr[4] = char_in[3];
    in_ptr[5] = char_in[2];
    in_ptr[6] = char_in[1];
    in_ptr[7] = char_in[0];
  }
}

inline void byte_swap_float(void *ptr, size_t nmemb)
{
  union {
    uint32_t uint_wrd;
    char c[sizeof(uint32_t)];
  } wrd;
  register char chr;
  int s;
  uint32_t* _u = (uint32_t*) ptr;
  for (s=0;s<nmemb;s++) {
    wrd.uint_wrd = _u[s];
    chr = wrd.c[0];
    wrd.c[0] = wrd.c[3];
    wrd.c[3] = chr;
    chr = wrd.c[2];
    wrd.c[2] = wrd.c[1];
    wrd.c[1] = chr;
    _u[s] = wrd.uint_wrd;
  }
}

//Defines overloaded function
inline void byte_swap(double *ptr, size_t nmemb){
  byte_swap_double(ptr, nmemb);
}

inline void byte_swap(float *ptr, size_t nmemb){
  byte_swap_float(ptr, nmemb);
}
