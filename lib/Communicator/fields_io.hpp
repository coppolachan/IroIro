/*! 
 * @file fields_io.hpp 
 * @brief Definition of MPI safe input/output routines for fields
 *
 * These are the parallel versions of the STL cout,cerr,cin.
 * The user can control which nodes are involved in output/input.
 */

#include "include/field.h"

namespace CCIO
{
  int SaveOnDisk(const std::vector<Field>& f, const char* filename);
}
