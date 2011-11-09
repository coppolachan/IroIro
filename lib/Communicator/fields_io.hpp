/*! 
 * @file fields_io.hpp 
 * @brief Definition of MPI safe input/output routines for fields
 *
 * These are the parallel versions of the STL cout,cerr,cin.
 * The user can control which nodes are involved in output/input.
 */

#include "include/common_code.hpp"

namespace CCIO
{
  class BinaryFileReader{
    //here the layout parameters are needed
    //for mpi compatibility
    int NumProcesses_;
    int NodeID_;
    bool isPrimary_;
    std::vector<int> local_lattice_;
    std::vector<int> global_lattice_;


  public:
    BinaryFileReader();

    ~BinaryFileReader();
    
    BinaryFileReader(const std::string& filename);
   
    //Specialized fuctions
    void read(Field& f);

    void read(std::vector<Field>& field_vec);


  };

  class BinaryFileWriter{
    //here the layout parameters are needed
    //for mpi compatibility
  public:
    BinaryFileWriter();

    ~BinaryFileWriter();
    
    BinaryFileWriter(const std::string& filename);
   
    //Specialized fuctions
    void write(Field& f);

    void write(std::vector<Field>& field_vec);


  };


  int SaveOnDisk(std::vector<Field>& f, const char* filename);

}
