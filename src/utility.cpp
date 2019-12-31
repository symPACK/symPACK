#include "sympack/Environment.hpp"
#include "sympack/utility.hpp"

#include <typeinfo>




using namespace std;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::cerr;
namespace symPACK{


// *********************************************************************
// IO functions
// TODO Move this to utility.hpp and make them inline functions
// *********************************************************************
//---------------------------------------------------------
Int SeparateRead(std::string name, std::istringstream& is)
{
  MPI_Barrier(MPI_COMM_WORLD);
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  //
  char filename[100];
  sprintf(filename, "%s_%d_%d", name.c_str(), mpirank, mpisize);  //cerr<<filename<<std::endl;
  std::ifstream fin(filename);
	if( !fin.good() ){
		throw std::logic_error( "File cannot be openeded!" );
	}
 
 	is.str( std::string(std::istreambuf_iterator<char>(fin), std::istreambuf_iterator<char>()) );
  fin.close();
  //
  MPI_Barrier(MPI_COMM_WORLD);
  return 0;
}

//---------------------------------------------------------
Int SeparateWrite(std::string name, std::ostringstream& os)
{
   MPI_Barrier(MPI_COMM_WORLD);
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  //
  char filename[100];
  sprintf(filename, "%s_%d_%d", name.c_str(), mpirank, mpisize);
  ofstream fout(filename);
	if( !fout.good() ){
		throw std::logic_error( "File cannot be openeded!" );
	}
  fout<<os.str();
  fout.close();
  //
  MPI_Barrier(MPI_COMM_WORLD);
  return 0;
}


//---------------------------------------------------------
Int SharedWrite(std::string name, std::ostringstream& os)
{
  MPI_Barrier(MPI_COMM_WORLD);
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  //
  if(mpirank==0) {
    ofstream fout(name.c_str());
		if( !fout.good() ){
			throw std::logic_error( "File cannot be openeded!" );
		}
    fout<<os.str();
    fout.close();
  }
  MPI_Barrier(MPI_COMM_WORLD);
  return 0;
}


//---------------------------------------------------------
Int SeparateWriteAscii(std::string name, std::ostringstream& os)
{
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  return 0;
}


//---------------------------------------------------------



void SetValue( std::vector<char>& vec, bool val ){
  fill(vec.begin(),vec.end(),val);
}




  void PtrSum( void *in, void *inout, int *len, MPI_Datatype *dptr ) 
  { 
    int i; 

    Ptr * pinout = (Ptr*)inout;
    Ptr * pin = (Ptr*)in;
#pragma unroll
    for (i=0; i< *len; ++i) { 
      pinout[i] += pin[i];
    } 
  } 


  void PtrMax( void *in, void *inout, int *len, MPI_Datatype *dptr ) 
  { 
    int i; 

    Ptr * pinout = (Ptr*)inout;
    Ptr * pin = (Ptr*)in;
#pragma unroll
    for (i=0; i< *len; ++i) { 
      pinout[i] = std::max(pinout[i], pin[i]);
    } 
  } 




  bool is_complex_type<std::complex<float> >::value = true;
  bool is_complex_type<std::complex<double> >::value = true;





}  // namespace SYMPACK
