/*! 
 * @file byteswap.hpp 
 * @brief Definition of function for trasforming a vector of big-endian doubles to litte-endian format
 *
 * Portable solution, performance is not necessary
 */

/*!
 * @brief Tool to trasform an array of doubles from big-endian to little-endian format
 */

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
  unsigned int j;
  char char_in[8];              /* characters used in byte swapping */
                                                                                
  char *in_ptr;
  float *float_ptr;           /* Pointer used in the float routines */
  
  for(j = 0, float_ptr = (float *) ptr;
      j < nmemb;
      j++, float_ptr++){
    
    
    in_ptr = (char *) float_ptr; /* Set the character pointer to
                                       point to the start of the float */
    
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
