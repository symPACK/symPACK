/// @file utility.cpp
/// @brief Implementation of the non-templated utility subroutines.
/// @author Lin Lin
/// @date 2012-09-20
#include "sympack/Environment.hpp"
#include "sympack/utility.hpp"

#include <typeinfo>

using namespace std;
using std::ifstream;
using std::ofstream;
using SYMPACK::vector;
using std::cerr;

namespace SYMPACK{


Int iam,np;

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
  sprintf(filename, "%s_%d_%d", name.c_str(), mpirank, mpisize);  //cerr<<filename<<endl;
  ifstream fin(filename);
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
/////  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
/////  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
/////  //
/////  char filename[100];
/////  sprintf(filename, "%s_%d_%d", name.c_str(), mpirank, mpisize);
/////  ofstream fout(filename, ios::trunc);
/////	if( !fout.good() ){
/////		throw std::logic_error( "File cannot be opened!" );
/////	}
/////  fout<<os;
/////  fout.close();
  //
  MPI_Barrier(MPI_COMM_WORLD);
  return 0;
}


//---------------------------------------------------------
void ReadDistSparseMatrix ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm )
{
	// Get the processor information within the current communicator
  MPI_Barrier( comm );
  Int mpirank;  MPI_Comm_rank(comm, &mpirank);
  Int mpisize;  MPI_Comm_size(comm, &mpisize);
	MPI_Status mpistat;
	std::ifstream fin;

  // Read basic information
	if( mpirank == 0 ){
		fin.open(filename);
		if( !fin.good() ){
			throw std::logic_error( "File cannot be openeded!" );
		}
		fin.read((char*)&pspmat.size, sizeof(Int));
		fin.read((char*)&pspmat.nnz,  sizeof(Int));
	}
	
	MPI_Bcast(&pspmat.size, 1, MPI_INT, 0, comm);
	MPI_Bcast(&pspmat.nnz,  1, MPI_INT, 0, comm);

	// Read colptr
	SYMPACK::vector<Int>  colptr(pspmat.size+1);
	if( mpirank == 0 ){
		Int tmp;
		fin.read((char*)&tmp, sizeof(Int));  

		if( tmp != pspmat.size+1 ){
			throw std::logic_error( "colptr is not of the right size." );
		}

		fin.read((char*)&colptr[0], sizeof(Int)*tmp);

	}

	MPI_Bcast(&colptr[0], pspmat.size+1, MPI_INT, 0, comm);

	// Compute the number of columns on each processor
	SYMPACK::vector<Int> numColLocalVec(mpisize);
	Int numColLocal, numColFirst;
	numColFirst = std::max(1,pspmat.size / mpisize);
  SetValue( numColLocalVec, numColFirst );
  numColLocalVec[mpisize-1] = pspmat.size - numColFirst * (mpisize-1);  // Modify the last entry	
	numColLocal = numColLocalVec[mpirank];

	pspmat.Local_.colptr.resize( numColLocal + 1 );
	for( Int i = 0; i < numColLocal + 1; i++ ){
		pspmat.Local_.colptr[i] = colptr[mpirank * numColFirst+i] - colptr[mpirank * numColFirst] + 1;
	}

	// Calculate nnz_loc on each processor
	pspmat.Local_.nnz = pspmat.Local_.colptr[numColLocal] - pspmat.Local_.colptr[0];

  pspmat.Local_.rowind.resize( pspmat.Local_.nnz );
	pspmat.nzvalLocal.resize ( pspmat.Local_.nnz );

	// Read and distribute the row indices
	if( mpirank == 0 ){
		Int tmp;
		fin.read((char*)&tmp, sizeof(Int));  

		if( tmp != pspmat.nnz ){
			std::ostringstream msg;
			msg 
				<< "The number of nonzeros in row indices do not match." << std::endl
				<< "nnz = " << pspmat.nnz << std::endl
				<< "size of row indices = " << tmp << std::endl;
			throw std::logic_error( msg.str().c_str() );
		}
		SYMPACK::vector<Idx> buf;
		Ptr numRead;
		for( Int ip = 0; ip < mpisize; ip++ ){
			numRead = colptr[ip*numColFirst + numColLocalVec[ip]] - 
				colptr[ip*numColFirst];
			buf.resize(numRead);
			fin.read( (char*)&buf[0], numRead*sizeof(Idx) );

			if( ip > 0 ){
				MPI_Send(&numRead, sizeof(numRead), MPI_BYTE, ip, 0, comm);
				MPI_Send(&buf[0], numRead*sizeof(Idx), MPI_BYTE, ip, 1, comm);
			}
			else{
        pspmat.Local_.rowind = buf;
			}
		}
	}
	else{
		Int numRead;
		MPI_Recv(&numRead, sizeof(numRead), MPI_BYTE, 0, 0, comm, &mpistat);
		if( numRead != pspmat.Local_.nnz ){
			std::ostringstream msg;
			msg << "The number of columns in row indices do not match." << std::endl
				<< "numRead  = " << numRead << std::endl
				<< "nnzLocal = " << pspmat.Local_.nnz << std::endl;
			throw std::logic_error( msg.str().c_str() );
		}

    pspmat.Local_.rowind.resize( numRead );
		MPI_Recv( &pspmat.Local_.rowind[0], numRead*sizeof(Idx), MPI_BYTE, 0, 1, comm, &mpistat );
	}
		
	// Read and distribute the nonzero values
	if( mpirank == 0 ){
		Int tmp;
		fin.read((char*)&tmp, sizeof(Int));  

		if( tmp != pspmat.nnz ){
			std::ostringstream msg;
			msg 
				<< "The number of nonzeros in values do not match." << std::endl
				<< "nnz = " << pspmat.nnz << std::endl
				<< "size of values = " << tmp << std::endl;
			throw std::logic_error( msg.str().c_str() );
		}
		SYMPACK::vector<Real> buf;
		Int numRead;
		for( Int ip = 0; ip < mpisize; ip++ ){
			numRead = colptr[ip*numColFirst + numColLocalVec[ip]] - 
				colptr[ip*numColFirst];
			buf.resize(numRead);
			fin.read( (char*)&buf[0], numRead*sizeof(Real) );

			if( ip > 0 ){
				MPI_Send(&numRead, 1, MPI_INT, ip, 0, comm);
				MPI_Send(&buf[0], numRead, MPI_DOUBLE, ip, 1, comm);
			}
			else{
        pspmat.nzvalLocal = buf;
			}
		}
	}
	else{
		Int numRead;
		MPI_Recv(&numRead, 1, MPI_INT, 0, 0, comm, &mpistat);
		if( numRead != pspmat.Local_.nnz ){
			std::ostringstream msg;
			msg << "The number of columns in values do not match." << std::endl
				<< "numRead  = " << numRead << std::endl
				<< "nnzLocal = " << pspmat.Local_.nnz << std::endl;
			throw std::logic_error( msg.str().c_str() );
		}

    pspmat.nzvalLocal.resize( numRead );
		MPI_Recv( &pspmat.nzvalLocal[0], numRead, MPI_DOUBLE, 0, 1, comm, &mpistat );
	}

	// Close the file
	if( mpirank == 0 ){
    fin.close();
	}



  MPI_Barrier( comm );


	return ;
}		// -----  end of function ReadDistSparseMatrix  ----- 


void ParaWriteDistSparseMatrix ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm )
{
  // Get the processor information within the current communicator
  MPI_Barrier( comm );
  Int mpirank;  MPI_Comm_rank(comm, &mpirank);
  Int mpisize;  MPI_Comm_size(comm, &mpisize);
  MPI_Status mpistat;
  Int err = 0;



  int filemode = MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_UNIQUE_OPEN;

  MPI_File fout;
  MPI_Status status;



  err = MPI_File_open(comm,(char*) filename, filemode, MPI_INFO_NULL,  &fout);

  if (err != MPI_SUCCESS) {
    throw std::logic_error( "File cannot be openeded!" );
  }

  // Write header
  if( mpirank == 0 ){
    err = MPI_File_write_at(fout, 0,(char*)&pspmat.size, 1, MPI_INT, &status);
    err = MPI_File_write_at(fout, sizeof(Int),(char*)&pspmat.nnz, 1, MPI_INT, &status);
  }


  // Compute the number of columns on each processor
  Int numColLocal = pspmat.Local_.colptr.size()-1;
  Int numColFirst = pspmat.size / mpisize;
  SYMPACK::vector<Int>  colptrChunk(numColLocal+1);

  Int prev_nz = 0;
  MPI_Exscan(&pspmat.Local_.nnz, &prev_nz, 1, MPI_INT, MPI_SUM, comm);

  for( Int i = 0; i < numColLocal + 1; i++ ){
    colptrChunk[i] = pspmat.Local_.colptr[i] + prev_nz;
  }


  MPI_Datatype memtype, filetype;
  MPI_Aint disps[6];
  int blklens[6];
  MPI_Datatype types[6] = {MPI_INT,MPI_INT, MPI_INT,MPI_INT, MPI_INT,MPI_DOUBLE};

  /* set block lengths (same for both types) */
  blklens[0] = (mpirank==0)?1:0;
  blklens[1] = numColLocal+1;
  blklens[2] = (mpirank==0)?1:0;
  blklens[3] = pspmat.Local_.nnz;
  blklens[4] = (mpirank==0)?1:0;
  blklens[5] = pspmat.Local_.nnz;




  //Calculate offsets
  MPI_Offset myColPtrOffset, myRowIdxOffset, myNzValOffset;
  myColPtrOffset = 3*sizeof(int) + (mpirank*numColFirst)*sizeof(Int);
  myRowIdxOffset = 3*sizeof(int) + (pspmat.size +1  +  prev_nz)*sizeof(Int);
  myNzValOffset = 4*sizeof(int) + (pspmat.size +1 +  pspmat.nnz)*sizeof(Int)+ prev_nz*sizeof(Real);
  disps[0] = 2*sizeof(int);
  disps[1] = myColPtrOffset;
  disps[2] = myRowIdxOffset;
  disps[3] = sizeof(int)+myRowIdxOffset;
  disps[4] = myNzValOffset;
  disps[5] = sizeof(int)+myNzValOffset;



#if ( _DEBUGlevel_ >= 1 )
  char msg[200];
  char * tmp = msg;
  tmp += sprintf(tmp,"P%d ",mpirank);
  for(int i = 0; i<6; ++i){
    if(i==5)
      tmp += sprintf(tmp, "%d [%d - %d] | ",i,disps[i],disps[i]+blklens[i]*sizeof(double));
    else
      tmp += sprintf(tmp, "%d [%d - %d] | ",i,disps[i],disps[i]+blklens[i]*sizeof(int));
  }
  tmp += sprintf(tmp,"\n");
  printf("%s",msg);
#endif




  MPI_Type_create_struct(6, blklens, disps, types, &filetype);
  MPI_Type_commit(&filetype);

  /* create memory type */
  Int np1 = pspmat.size+1;
  MPI_Address( (void *)&np1,  &disps[0]);
  MPI_Address(&colptrChunk[0], &disps[1]);
  MPI_Address( (void *)&pspmat.nnz,  &disps[2]);
  MPI_Address((void *)&pspmat.Local_.rowind[0],  &disps[3]);
  MPI_Address( (void *)&pspmat.nnz,  &disps[4]);
  MPI_Address((void *)&pspmat.nzvalLocal[0],   &disps[5]);

  MPI_Type_create_struct(6, blklens, disps, types, &memtype);
  MPI_Type_commit(&memtype);



  /* set file view */
  err = MPI_File_set_view(fout, 0, MPI_BYTE, filetype, "native",MPI_INFO_NULL);

  /* everyone writes their own row offsets, columns, and 
   * data with one big noncontiguous write (in memory and 
   * file)
   */
  err = MPI_File_write_all(fout, MPI_BOTTOM, 1, memtype, &status);

  MPI_Type_free(&filetype);
  MPI_Type_free(&memtype);





  MPI_Barrier( comm );

  MPI_File_close(&fout);

  return ;
}		// -----  end of function ParaWriteDistSparseMatrix  ----- 


void ParaReadDistSparseMatrix ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm )
{





  // Get the processor information within the current communicator
  MPI_Barrier( comm );
  Int mpirank;  MPI_Comm_rank(comm, &mpirank);
  Int mpisize;  MPI_Comm_size(comm, &mpisize);
  MPI_Status mpistat;
  MPI_Datatype type;
  int lens[3];
  MPI_Aint disps[3];
  MPI_Datatype types[3];





  Int err = 0;

  int filemode = MPI_MODE_RDONLY | MPI_MODE_UNIQUE_OPEN;

  MPI_File fin;
  MPI_Status status;

  // FIXME Maybe change to MPI_Comm_Dup.
  pspmat.comm = comm;

  err = MPI_File_open(comm,(char*) filename, filemode, MPI_INFO_NULL,  &fin);

  if (err != MPI_SUCCESS) {
    #ifdef USE_ABORT
abort();
#endif
throw std::logic_error( "File cannot be opened!" );
  }

  // FIXME Note that nnz uses the Int data type for consistency of writing / reading
  // Read header
  if( mpirank == 0 ){
    err = MPI_File_read_at(fin, 0,(char*)&pspmat.size, 1, MPI_INT, &status);
    err = MPI_File_read_at(fin, sizeof(Int),(char*)&pspmat.nnz, 1, MPI_INT, &status);
  }


  /* define a struct that describes all our data */
  lens[0] = sizeof(pspmat.size);
  lens[1] = sizeof(pspmat.nnz);
  MPI_Address(&pspmat.size, &disps[0]);
  MPI_Address(&pspmat.nnz, &disps[1]);
  types[0] = MPI_BYTE;
  types[1] = MPI_BYTE;
  MPI_Type_struct(2, lens, disps, types, &type);
  MPI_Type_commit(&type);


  /* broadcast the header data to everyone */
  MPI_Bcast(MPI_BOTTOM, 1, type, 0, comm);

  MPI_Type_free(&type);

  pspmat.Local_.size = pspmat.size;
  pspmat.Global_.size = pspmat.size;
  pspmat.Global_.nnz = pspmat.nnz;
 

  // Compute the number of columns on each processor
  SYMPACK::vector<Int> numColLocalVec(mpisize);
  Int numColLocal, numColFirst;
  numColFirst = pspmat.size / mpisize;
  SetValue( numColLocalVec, numColFirst );
  numColLocalVec[mpisize-1] = pspmat.size - numColFirst * (mpisize-1);  // Modify the last entry  
  numColLocal = numColLocalVec[mpirank];
  pspmat.Local_.colptr.resize( numColLocal + 1 );



  MPI_Offset myColPtrOffset = (2 + ((mpirank==0)?0:1) )*sizeof(int) + (mpirank*numColFirst)*sizeof(Int);

  Int np1 = 0;
  lens[0] = ((mpirank==0)?1:0)*sizeof(int);
  lens[1] = (numColLocal + 1)*sizeof(int);

  MPI_Address(&np1, &disps[0]);
  MPI_Address(&pspmat.Local_.colptr[0], &disps[1]);

  MPI_Type_hindexed(2, lens, disps, MPI_BYTE, &type);
  MPI_Type_commit(&type);

  err= MPI_File_read_at_all(fin, myColPtrOffset, MPI_BOTTOM, 1, type, &status);

  if (err != MPI_SUCCESS) {
#ifdef USE_ABORT
    abort();
#endif
    throw std::logic_error( "error reading colptr" );
  }
  MPI_Type_free(&type);

  if(typeid(int)!=typeid(Ptr)){
    np1 = (Int)*((int*)&np1);
    int * cptr = (int*)&pspmat.Local_.colptr[0];
    for(int64_t i = pspmat.Local_.colptr.size()-1 ; i>=0; i--){
      pspmat.Local_.colptr[i] = (Ptr)cptr[i];
    } 
  }


  // Calculate nnz_loc on each processor
  pspmat.Local_.nnz = pspmat.Local_.colptr[numColLocal] - pspmat.Local_.colptr[0];


  pspmat.Local_.rowind.resize( pspmat.Local_.nnz );
  pspmat.nzvalLocal.resize ( pspmat.Local_.nnz );

  //read rowIdx
  MPI_Offset myRowIdxOffset = (3 + ((mpirank==0)?0:1) )*sizeof(int) + (pspmat.size+1 + (pspmat.Local_.colptr[0]-1))*sizeof(int);

  lens[0] = (mpirank==0)?1:0;
  lens[1] = pspmat.Local_.nnz;

  MPI_Address(&np1, &disps[0]);
  MPI_Address(&pspmat.Local_.rowind[0], &disps[1]);

  MPI_Type_hindexed(2, lens, disps, MPI_INT, &type);
  MPI_Type_commit(&type);

  err= MPI_File_read_at_all(fin, myRowIdxOffset, MPI_BOTTOM, 1, type,&status);

  if (err != MPI_SUCCESS) {
    #ifdef USE_ABORT
abort();
#endif
throw std::logic_error( "error reading rowind" );
  }
  MPI_Type_free(&type);


  if(typeid(int)!=typeid(Idx)){
    np1 = (Int)*((int*)&np1);
    int * rptr = (int*)&pspmat.Local_.rowind[0];
    for(int64_t i = pspmat.Local_.rowind.size()-1 ; i>=0; i--){
      pspmat.Local_.rowind[i] = (Idx)rptr[i];
    } 
  }





  //read nzval
  MPI_Offset myNzValOffset = (4 + ((mpirank==0)?0:1) )*sizeof(int) + (pspmat.size+1 + pspmat.nnz)*sizeof(Int) + (pspmat.Local_.colptr[0]-1)*sizeof(double);

  lens[0] = (mpirank==0)?1:0;
  lens[1] = pspmat.Local_.nnz;

  MPI_Address(&np1, &disps[0]);
  MPI_Address(&pspmat.nzvalLocal[0], &disps[1]);

  types[0] = MPI_INT;
  types[1] = MPI_DOUBLE;

  MPI_Type_create_struct(2, lens, disps, types, &type);
  MPI_Type_commit(&type);

  err = MPI_File_read_at_all(fin, myNzValOffset, MPI_BOTTOM, 1, type,&status);

  if (err != MPI_SUCCESS) {
    #ifdef USE_ABORT
abort();
#endif
throw std::logic_error( "error reading nzval" );
  }

  MPI_Type_free(&type);

  if(typeid(double)!=typeid(Real)){
    np1 = (Int)*((int*)&np1);
    double * rptr = (double*)&pspmat.nzvalLocal[0];
    for(int64_t i = pspmat.nzvalLocal.size()-1 ; i>=0; i--){
      pspmat.nzvalLocal[i] = (Idx)rptr[i];
    } 
  }




  //convert to local references
  for( Int i = 1; i < numColLocal + 1; i++ ){
    pspmat.Local_.colptr[i] = pspmat.Local_.colptr[i] -  pspmat.Local_.colptr[0] + 1;
  }
  pspmat.Local_.colptr[0]=1;

  MPI_Barrier( comm );

  MPI_File_close(&fin);

  return ;





//#ifndef _RELEASE_
//  PushCallStack("ParaReadDistSparseMatrix");
//#endif
//  // Get the processor information within the current communicator
//  MPI_Barrier( comm );
//  Int mpirank;  MPI_Comm_rank(comm, &mpirank);
//  Int mpisize;  MPI_Comm_size(comm, &mpisize);
//  MPI_Status mpistat;
//  MPI_Datatype type;
//  int lens[3];
//  MPI_Aint disps[3];
//  MPI_Datatype types[3];
//  Int err = 0;
//
//
//
//  int filemode = MPI_MODE_RDONLY | MPI_MODE_UNIQUE_OPEN;
//
//  MPI_File fin;
//  MPI_Status status;
//
//
//  err = MPI_File_open(comm,(char*) filename, filemode, MPI_INFO_NULL,  &fin);
//
//  if (err != MPI_SUCCESS) {
//    throw std::logic_error( "File cannot be openeded!" );
//  }
//
//  // Read header
//  if( mpirank == 0 ){
//    err = MPI_File_read_at(fin, 0,(char*)&pspmat.size, 1, MPI_INT, &status);
//    err = MPI_File_read_at(fin, sizeof(Int),(char*)&pspmat.nnz, 1, MPI_INT, &status);
//  }
//
//  pspmat.comm = comm;
//
//  /* define a struct that describes all our data */
//  lens[0] = 1;
//  lens[1] = 1;
//  MPI_Address(&pspmat.size, &disps[0]);
//  MPI_Address(&pspmat.nnz, &disps[1]);
//  types[0] = MPI_INT;
//  types[1] = MPI_INT;
//  MPI_Type_struct(2, lens, disps, types, &type);
//  MPI_Type_commit(&type);
//
//  /* broadcast the header data to everyone */
//  MPI_Bcast(MPI_BOTTOM, 1, type, 0, comm);
//
//  MPI_Type_free(&type);
//
//  // Compute the number of columns on each processor
//  SYMPACK::vector<Int> numColLocalVec(mpisize);
//  Int numColLocal, numColFirst;
//  numColFirst = pspmat.size / mpisize;
//  SetValue( numColLocalVec, numColFirst );
//  numColLocalVec[mpisize-1] = pspmat.size - numColFirst * (mpisize-1);  // Modify the last entry	
//  numColLocal = numColLocalVec[mpirank];
//  pspmat.Local_.colptr.resize( numColLocal + 1 );
//
//
//
//  MPI_Offset myColPtrOffset = (2 + ((mpirank==0)?0:1) )*sizeof(int) + (mpirank*numColFirst)*sizeof(Int);
//
//  Int np1 = 0;
//  lens[0] = (mpirank==0)?1:0;
//  lens[1] = numColLocal + 1;
//
//  MPI_Address(&np1, &disps[0]);
//  MPI_Address(pspmat.Local_.colptr.Data(), &disps[1]);
//
//  MPI_Type_hindexed(2, lens, disps, MPI_INT, &type);
//  MPI_Type_commit(&type);
//
//  err= MPI_File_read_at_all(fin, myColPtrOffset, MPI_BOTTOM, 1, type, &status);
//
//  if (err != MPI_SUCCESS) {
//    throw std::logic_error( "error reading colptr" );
//  }
//  MPI_Type_free(&type);
//
//  // Calculate nnz_loc on each processor
//  pspmat.Local_.nnz = pspmat.Local_.colptr[numColLocal] - pspmat.Local_.colptr[0];
//
//
//  pspmat.Local_.rowind.resize( pspmat.Local_.nnz );
//  pspmat.nzvalLocal.resize ( pspmat.Local_.nnz );
//
//  //read rowIdx
//  MPI_Offset myRowIdxOffset = (3 + ((mpirank==0)?-1:0) )*sizeof(int) + (pspmat.size+1 + pspmat.Local_.colptr[0])*sizeof(Int);
//
//  lens[0] = (mpirank==0)?1:0;
//  lens[1] = pspmat.Local_.nnz;
//
//  MPI_Address(&np1, &disps[0]);
//  MPI_Address(pspmat.Local_.rowind.Data(), &disps[1]);
//
//  MPI_Type_hindexed(2, lens, disps, MPI_INT, &type);
//  MPI_Type_commit(&type);
//
//  err= MPI_File_read_at_all(fin, myRowIdxOffset, MPI_BOTTOM, 1, type,&status);
//
//  if (err != MPI_SUCCESS) {
//    throw std::logic_error( "error reading rowind" );
//  }
//  MPI_Type_free(&type);
//
//
//  //read nzval
//  MPI_Offset myNzValOffset = (3 + ((mpirank==0)?-1:0) )*sizeof(int) + (pspmat.size+1 + pspmat.nnz)*sizeof(Int) + pspmat.Local_.colptr[0]*sizeof(double);
//
//  lens[0] = (mpirank==0)?1:0;
//  lens[1] = pspmat.Local_.nnz;
//
//  MPI_Address(&np1, &disps[0]);
//  MPI_Address(pspmat.nzvalLocal.Data(), &disps[1]);
//
//  types[0] = MPI_INT;
//  types[1] = MPI_DOUBLE;
//
//  MPI_Type_create_struct(2, lens, disps, types, &type);
//  MPI_Type_commit(&type);
//
//  err = MPI_File_read_at_all(fin, myNzValOffset, MPI_BOTTOM, 1, type,&status);
//
//  if (err != MPI_SUCCESS) {
//    throw std::logic_error( "error reading nzval" );
//  }
//
//  MPI_Type_free(&type);
//
//
//  //convert to local references
//  for( Int i = 1; i < numColLocal + 1; i++ ){
//    pspmat.Local_.colptr[i] = pspmat.Local_.colptr[i] -  pspmat.Local_.colptr[0] + 1;
//  }
//  pspmat.Local_.colptr[0]=1;
//
//  MPI_Barrier( comm );
//
//  MPI_File_close(&fin);
//#ifndef _RELEASE_
//  PopCallStack();
//#endif
//
//  return ;
}		// -----  end of function ParaReadDistSparseMatrix  ----- 


//TODO we should do the same than for Real 
void ParaReadDistSparseMatrix ( const char* filename, DistSparseMatrix<Complex>& pspmat, MPI_Comm comm )
{





  // Get the processor information within the current communicator
  MPI_Barrier( comm );
  Int mpirank;  MPI_Comm_rank(comm, &mpirank);
  Int mpisize;  MPI_Comm_size(comm, &mpisize);
  MPI_Status mpistat;
  MPI_Datatype type;
  int lens[3];
  MPI_Aint disps[3];
  MPI_Datatype types[3];





  Int err = 0;

  int filemode = MPI_MODE_RDONLY | MPI_MODE_UNIQUE_OPEN;

  MPI_File fin;
  MPI_Status status;

  // FIXME Maybe change to MPI_Comm_Dup.
  pspmat.comm = comm;

  err = MPI_File_open(comm,(char*) filename, filemode, MPI_INFO_NULL,  &fin);

  if (err != MPI_SUCCESS) {
    #ifdef USE_ABORT
abort();
#endif
throw std::logic_error( "File cannot be opened!" );
  }

  // FIXME Note that nnz uses the Int data type for consistency of writing / reading
  // Read header
  if( mpirank == 0 ){
    err = MPI_File_read_at(fin, 0,(char*)&pspmat.size, 1, MPI_INT, &status);
    err = MPI_File_read_at(fin, sizeof(Int),(char*)&pspmat.nnz, 1, MPI_INT, &status);
  }


  /* define a struct that describes all our data */
  lens[0] = 1;
  lens[1] = 1;
  MPI_Address(&pspmat.size, &disps[0]);
  MPI_Address(&pspmat.nnz, &disps[1]);
  types[0] = MPI_INT;
  types[1] = MPI_INT;
  MPI_Type_struct(2, lens, disps, types, &type);
  MPI_Type_commit(&type);


  /* broadcast the header data to everyone */
  MPI_Bcast(MPI_BOTTOM, 1, type, 0, comm);

  MPI_Type_free(&type);

  pspmat.Local_.size = pspmat.size;
  pspmat.Global_.size = pspmat.size;
  pspmat.Global_.nnz = pspmat.nnz;
 
  // Compute the number of columns on each processor
  SYMPACK::vector<Int> numColLocalVec(mpisize);
  Int numColLocal, numColFirst;
  numColFirst = pspmat.size / mpisize;
  SetValue( numColLocalVec, numColFirst );
  numColLocalVec[mpisize-1] = pspmat.size - numColFirst * (mpisize-1);  // Modify the last entry  
  numColLocal = numColLocalVec[mpirank];
  pspmat.Local_.colptr.resize( numColLocal + 1 );



  MPI_Offset myColPtrOffset = (2 + ((mpirank==0)?0:1) )*sizeof(int) + (mpirank*numColFirst)*sizeof(Int);

  Int np1 = 0;
  lens[0] = (mpirank==0)?1:0;
  lens[1] = numColLocal + 1;

  MPI_Address(&np1, &disps[0]);
  MPI_Address(&pspmat.Local_.colptr[0], &disps[1]);

  MPI_Type_hindexed(2, lens, disps, MPI_INT, &type);
  MPI_Type_commit(&type);

  err= MPI_File_read_at_all(fin, myColPtrOffset, MPI_BOTTOM, 1, type, &status);

  if (err != MPI_SUCCESS) {
#ifdef USE_ABORT
    abort();
#endif
    throw std::logic_error( "error reading colptr" );
  }
  MPI_Type_free(&type);

  // Calculate nnz_loc on each processor
  pspmat.Local_.nnz = pspmat.Local_.colptr[numColLocal] - pspmat.Local_.colptr[0];


  pspmat.Local_.rowind.resize( pspmat.Local_.nnz );
  pspmat.nzvalLocal.resize ( pspmat.Local_.nnz );

  //read rowIdx
  MPI_Offset myRowIdxOffset = (3 + ((mpirank==0)?0:1) )*sizeof(int) + (pspmat.size+1 + (pspmat.Local_.colptr[0]-1))*sizeof(int);

  lens[0] = (mpirank==0)?1:0;
  lens[1] = pspmat.Local_.nnz;

  MPI_Address(&np1, &disps[0]);
  MPI_Address(&pspmat.Local_.rowind[0], &disps[1]);

  MPI_Type_hindexed(2, lens, disps, MPI_INT, &type);
  MPI_Type_commit(&type);

  err= MPI_File_read_at_all(fin, myRowIdxOffset, MPI_BOTTOM, 1, type,&status);

  if (err != MPI_SUCCESS) {
    #ifdef USE_ABORT
abort();
#endif
throw std::logic_error( "error reading rowind" );
  }
  MPI_Type_free(&type);


  //read nzval
  MPI_Offset myNzValOffset = (4 + ((mpirank==0)?0:1) )*sizeof(int) + (pspmat.size+1 + pspmat.nnz)*sizeof(Int) + (pspmat.Local_.colptr[0]-1)*sizeof(Complex);

  lens[0] = (mpirank==0)?1:0;
  lens[1] = pspmat.Local_.nnz;

  MPI_Address(&np1, &disps[0]);
  MPI_Address(&pspmat.nzvalLocal[0], &disps[1]);

  types[0] = MPI_INT;
  types[1] = MPI_COMPLEX;

  MPI_Type_create_struct(2, lens, disps, types, &type);
  MPI_Type_commit(&type);

  err = MPI_File_read_at_all(fin, myNzValOffset, MPI_BOTTOM, 1, type,&status);

  if (err != MPI_SUCCESS) {
    #ifdef USE_ABORT
abort();
#endif
throw std::logic_error( "error reading nzval" );
  }

  MPI_Type_free(&type);


  //convert to local references
  for( Int i = 1; i < numColLocal + 1; i++ ){
    pspmat.Local_.colptr[i] = pspmat.Local_.colptr[i] -  pspmat.Local_.colptr[0] + 1;
  }
  pspmat.Local_.colptr[0]=1;

  MPI_Barrier( comm );

  MPI_File_close(&fin);

  return ;





//#ifndef _RELEASE_
//  PushCallStack("ParaReadDistSparseMatrix");
//#endif
//  // Get the processor information within the current communicator
//  MPI_Barrier( comm );
//  Int mpirank;  MPI_Comm_rank(comm, &mpirank);
//  Int mpisize;  MPI_Comm_size(comm, &mpisize);
//  MPI_Status mpistat;
//  MPI_Datatype type;
//  int lens[3];
//  MPI_Aint disps[3];
//  MPI_Datatype types[3];
//  Int err = 0;
//
//
//
//  int filemode = MPI_MODE_RDONLY | MPI_MODE_UNIQUE_OPEN;
//
//  MPI_File fin;
//  MPI_Status status;
//
//
//  err = MPI_File_open(comm,(char*) filename, filemode, MPI_INFO_NULL,  &fin);
//
//  if (err != MPI_SUCCESS) {
//    throw std::logic_error( "File cannot be openeded!" );
//  }
//
//  // Read header
//  if( mpirank == 0 ){
//    err = MPI_File_read_at(fin, 0,(char*)&pspmat.size, 1, MPI_INT, &status);
//    err = MPI_File_read_at(fin, sizeof(Int),(char*)&pspmat.nnz, 1, MPI_INT, &status);
//  }
//
//  pspmat.comm = comm;
//
//  /* define a struct that describes all our data */
//  lens[0] = 1;
//  lens[1] = 1;
//  MPI_Address(&pspmat.size, &disps[0]);
//  MPI_Address(&pspmat.nnz, &disps[1]);
//  types[0] = MPI_INT;
//  types[1] = MPI_INT;
//  MPI_Type_struct(2, lens, disps, types, &type);
//  MPI_Type_commit(&type);
//
//  /* broadcast the header data to everyone */
//  MPI_Bcast(MPI_BOTTOM, 1, type, 0, comm);
//
//  MPI_Type_free(&type);
//
//  // Compute the number of columns on each processor
//  SYMPACK::vector<Int> numColLocalVec(mpisize);
//  Int numColLocal, numColFirst;
//  numColFirst = pspmat.size / mpisize;
//  SetValue( numColLocalVec, numColFirst );
//  numColLocalVec[mpisize-1] = pspmat.size - numColFirst * (mpisize-1);  // Modify the last entry	
//  numColLocal = numColLocalVec[mpirank];
//  pspmat.Local_.colptr.resize( numColLocal + 1 );
//
//
//
//  MPI_Offset myColPtrOffset = (2 + ((mpirank==0)?0:1) )*sizeof(int) + (mpirank*numColFirst)*sizeof(Int);
//
//  Int np1 = 0;
//  lens[0] = (mpirank==0)?1:0;
//  lens[1] = numColLocal + 1;
//
//  MPI_Address(&np1, &disps[0]);
//  MPI_Address(pspmat.Local_.colptr.Data(), &disps[1]);
//
//  MPI_Type_hindexed(2, lens, disps, MPI_INT, &type);
//  MPI_Type_commit(&type);
//
//  err= MPI_File_read_at_all(fin, myColPtrOffset, MPI_BOTTOM, 1, type, &status);
//
//  if (err != MPI_SUCCESS) {
//    throw std::logic_error( "error reading colptr" );
//  }
//  MPI_Type_free(&type);
//
//  // Calculate nnz_loc on each processor
//  pspmat.Local_.nnz = pspmat.Local_.colptr[numColLocal] - pspmat.Local_.colptr[0];
//
//
//  pspmat.Local_.rowind.resize( pspmat.Local_.nnz );
//  pspmat.nzvalLocal.resize ( pspmat.Local_.nnz );
//
//  //read rowIdx
//  MPI_Offset myRowIdxOffset = (3 + ((mpirank==0)?-1:0) )*sizeof(int) + (pspmat.size+1 + pspmat.Local_.colptr[0])*sizeof(Int);
//
//  lens[0] = (mpirank==0)?1:0;
//  lens[1] = pspmat.Local_.nnz;
//
//  MPI_Address(&np1, &disps[0]);
//  MPI_Address(pspmat.Local_.rowind.Data(), &disps[1]);
//
//  MPI_Type_hindexed(2, lens, disps, MPI_INT, &type);
//  MPI_Type_commit(&type);
//
//  err= MPI_File_read_at_all(fin, myRowIdxOffset, MPI_BOTTOM, 1, type,&status);
//
//  if (err != MPI_SUCCESS) {
//    throw std::logic_error( "error reading rowind" );
//  }
//  MPI_Type_free(&type);
//
//
//  //read nzval
//  MPI_Offset myNzValOffset = (3 + ((mpirank==0)?-1:0) )*sizeof(int) + (pspmat.size+1 + pspmat.nnz)*sizeof(Int) + pspmat.Local_.colptr[0]*sizeof(double);
//
//  lens[0] = (mpirank==0)?1:0;
//  lens[1] = pspmat.Local_.nnz;
//
//  MPI_Address(&np1, &disps[0]);
//  MPI_Address(pspmat.nzvalLocal.Data(), &disps[1]);
//
//  types[0] = MPI_INT;
//  types[1] = MPI_DOUBLE;
//
//  MPI_Type_create_struct(2, lens, disps, types, &type);
//  MPI_Type_commit(&type);
//
//  err = MPI_File_read_at_all(fin, myNzValOffset, MPI_BOTTOM, 1, type,&status);
//
//  if (err != MPI_SUCCESS) {
//    throw std::logic_error( "error reading nzval" );
//  }
//
//  MPI_Type_free(&type);
//
//
//  //convert to local references
//  for( Int i = 1; i < numColLocal + 1; i++ ){
//    pspmat.Local_.colptr[i] = pspmat.Local_.colptr[i] -  pspmat.Local_.colptr[0] + 1;
//  }
//  pspmat.Local_.colptr[0]=1;
//
//  MPI_Barrier( comm );
//
//  MPI_File_close(&fin);
//#ifndef _RELEASE_
//  PopCallStack();
//#endif
//
//  return ;
}		// -----  end of function ParaReadDistSparseMatrix  ----- 



void ReadDistSparseMatrixFormatted ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm )
{
	// Get the processor information within the current communicator
  MPI_Barrier( comm );
  Int mpirank;  MPI_Comm_rank(comm, &mpirank);
  Int mpisize;  MPI_Comm_size(comm, &mpisize);
	MPI_Status mpistat;
	std::ifstream fin;

  // Read basic information
	if( mpirank == 0 ){
		fin.open(filename);
		if( !fin.good() ){
			throw std::logic_error( "File cannot be openeded!" );
		}
		Int dummy;
		fin >> pspmat.size >> dummy;
		fin >> pspmat.nnz;
		// FIXME this is temporary and only applies to 4*4 matrix.
//	  fin	>> dummy;
	}
	
	MPI_Bcast(&pspmat.size, 1, MPI_INT, 0, comm);
	MPI_Bcast(&pspmat.nnz,  1, MPI_INT, 0, comm);

	// Read colptr

	SYMPACK::vector<Int>  colptr(pspmat.size+1);
	if( mpirank == 0 ){
		Int* ptr = &colptr[0];
		for( Int i = 0; i < pspmat.size+1; i++ )
			fin >> *(ptr++);
	}

	MPI_Bcast(&colptr[0], pspmat.size+1, MPI_INT, 0, comm);

	// Compute the number of columns on each processor
	SYMPACK::vector<Int> numColLocalVec(mpisize);
	Int numColLocal, numColFirst;
	numColFirst = pspmat.size / mpisize;
  SetValue( numColLocalVec, numColFirst );
  numColLocalVec[mpisize-1] = pspmat.size - numColFirst * (mpisize-1);  // Modify the last entry	
	numColLocal = numColLocalVec[mpirank];

	pspmat.Local_.colptr.resize( numColLocal + 1 );
	for( Int i = 0; i < numColLocal + 1; i++ ){
		pspmat.Local_.colptr[i] = colptr[mpirank * numColFirst+i] - colptr[mpirank * numColFirst] + 1;
	}

	// Calculate nnz_loc on each processor
	pspmat.Local_.nnz = pspmat.Local_.colptr[numColLocal] - pspmat.Local_.colptr[0];

  pspmat.Local_.rowind.resize( pspmat.Local_.nnz );
	pspmat.nzvalLocal.resize ( pspmat.Local_.nnz );

	// Read and distribute the row indices
	if( mpirank == 0 ){
		Int tmp;
		SYMPACK::vector<Idx> buf;
		Ptr numRead;
		for( Int ip = 0; ip < mpisize; ip++ ){
			numRead = colptr[ip*numColFirst + numColLocalVec[ip]] - 
				colptr[ip*numColFirst];
			buf.resize(numRead);
			Idx *ptr = &buf[0];
			for( Int i = 0; i < numRead; i++ ){
				fin >> *(ptr++);
			}
			if( ip > 0 ){
				MPI_Send(&numRead, sizeof(numRead), MPI_BYTE, ip, 0, comm);
				MPI_Send(&buf[0], numRead*sizeof(Idx), MPI_BYTE, ip, 1, comm);
			}
			else{
        pspmat.Local_.rowind = buf;
			}
		}
	}
	else{
		Ptr numRead;
		MPI_Recv(&numRead, sizeof(numRead), MPI_BYTE, 0, 0, comm, &mpistat);
		if( numRead != pspmat.Local_.nnz ){
			std::ostringstream msg;
			msg << "The number of columns in row indices do not match." << std::endl
				<< "numRead  = " << numRead << std::endl
				<< "nnzLocal = " << pspmat.Local_.nnz << std::endl;
			throw std::logic_error( msg.str().c_str() );
		}

    pspmat.Local_.rowind.resize( numRead );
		MPI_Recv( &pspmat.Local_.rowind[0], numRead*sizeof(Idx), MPI_BYTE, 0, 1, comm, &mpistat );
	}
		
//	std::cout << "Proc " << mpirank << " outputs Local_.rowind.size() = " 
//		<< pspmat.Local_.rowind.size() << endl;


	// Read and distribute the nonzero values
	if( mpirank == 0 ){
		Int tmp;
		SYMPACK::vector<Real> buf;
		Int numRead;
		for( Int ip = 0; ip < mpisize; ip++ ){
			numRead = colptr[ip*numColFirst + numColLocalVec[ip]] - 
				colptr[ip*numColFirst];
			buf.resize(numRead);
			Real *ptr = &buf[0];
			for( Int i = 0; i < numRead; i++ ){
				fin >> *(ptr++);
			}
			if( ip > 0 ){
				MPI_Send(&numRead, 1, MPI_INT, ip, 0, comm);
				MPI_Send(&buf[0], numRead, MPI_DOUBLE, ip, 1, comm);
			}
			else{
        pspmat.nzvalLocal = buf;
			}
		}
	}
	else{
		Int numRead;
		MPI_Recv(&numRead, 1, MPI_INT, 0, 0, comm, &mpistat);
		if( numRead != pspmat.Local_.nnz ){
			std::ostringstream msg;
			msg << "The number of columns in values do not match." << std::endl
				<< "numRead  = " << numRead << std::endl
				<< "nnzLocal = " << pspmat.Local_.nnz << std::endl;
			throw std::logic_error( msg.str().c_str() );
		}

    pspmat.nzvalLocal.resize( numRead );
		MPI_Recv( &pspmat.nzvalLocal[0], numRead, MPI_DOUBLE, 0, 1, comm, &mpistat );
	}

	// Close the file
	if( mpirank == 0 ){
    fin.close();
	}



  MPI_Barrier( comm );


	return ;
}		// -----  end of function ReadDistSparseMatrixFormatted  ----- 

//void
//GetDiagonal ( const DistSparseMatrix<Complex>& A, 
//		NumVec<Complex>& diag )
//{
//#ifndef _RELEASE_
//	PushCallStack("GetDiagonal");
//#endif
//	Int mpirank, mpisize;
//	MPI_Comm_rank( A.comm, &mpirank );
//	MPI_Comm_size( A.comm, &mpisize );
//
//  NumVec<Complex>	 diagLocal( A.size );
//	SetValue( diagLocal, Z_ZERO );
//	diag.resize( A.size );
//	SetValue( diag, Z_ZERO );
//
//	Int numColFirst = A.size / mpisize;
//	Int firstCol    = mpirank * numColFirst;
//	Int numColLocal = A.Local_.colptr.size() - 1;
//
//#if ( _DEBUGlevel_ >= 1 )
//	statusOFS << "numColFirst = " << numColFirst << std::endl;
//	statusOFS << "A.nzvalLocal.size = " << A.nzvalLocal.size() << std::endl;
//	statusOFS << "A.Local_.nnz = " << A.Local_.nnz << std::endl;
//#endif
//
//	// Note that the indices in DistSparseMatrix follows the FORTRAN convention
//  for( Int j = 0; j < numColLocal; j++ ){
//		Int jcol = j + firstCol + 1;
//		Int numRow = A.Local_.colptr(j+1) - A.Local_.colptr(j);
//		const Int* rowPtr = &A.Local_.rowind( A.Local_.colptr(j) - 1 );
//		// NOTE: The rows in DistSparseMatrix are not necessarily ordered.
//		// So lower_bound cannot be used here for fast searching. find has to be used. 
//		const Int* ptr = find( rowPtr, rowPtr + numRow, jcol ); 
//		if( ptr == rowPtr + numRow ){
//			std::ostringstream msg;
//			msg << "Serious problem. Did not find the row corresponding to the column." << std::endl
//				<< "This happens when j = " << j << ", jcol = " << jcol << ", and the row indices are " << std::endl
//				<< SYMPACK::vector<Int>( numRow, false, const_cast<Int*>(rowPtr) ) << std::endl;
//			throw std::logic_error( msg.str().c_str() );
//		}
//		Int diagIdx = ptr - A.Local_.rowind.Data();
//    diagLocal( jcol - 1 ) = A.nzvalLocal( diagIdx );
//	}
//
//	mpi::Allreduce( &diagLocal[0], &diag[0], A.size, MPI_SUM, A.comm );
//
//#ifndef _RELEASE_
//	PopCallStack();
//#endif
//
//	return ;
//}		// -----  end of function GetDiagonal  ----- 


void SetValue( SYMPACK::vector<char>& vec, bool val ){
  fill(vec.begin(),vec.end(),val);
}



#if 0
template <typename SCALAR, typename INSCALAR >
int ReadHB_PARA(std::string & filename, DistSparseMatrix<SCALAR> & HMat){
  MPI_Comm & workcomm = HMat.comm;

    int mpirank;
    MPI_Comm_rank(workcomm,&mpirank);
  
    int mpisize;
    MPI_Comm_size(workcomm,&mpisize);

  ifstream infile;
  infile.open(filename.c_str());

  string line;
  //read xadj on the first line of the input file
  stringstream iss;
  //skip 1st line
  if(getline(infile, line)){}
  Idx colptrCnt;
  Ptr rowindCnt;
  Ptr nzvalCnt;
  if(getline(infile, line)){
    iss.str("");
    iss.clear();
    iss<<line;
    Ptr dummy;
    iss>>dummy;
    iss>>colptrCnt>>rowindCnt>>nzvalCnt;
  }
  //read from third line

  auto & n = HMat.size;
  auto & nnz = HMat.nnz;


  auto m = n;
  if(getline(infile, line))
  {
    iss.str("");
    iss.clear();
    iss<<line;
    string type;
    iss>>type;
    iss>>m>>n>>nnz;
  }

  //compute local number of columns
  Idx nlocal = (mpirank<mpisize-1)?n/mpisize:n-mpirank*(Idx)(n/mpisize);
  Idx firstNode = mpirank*(int)(n/mpisize) + 1;
  //initialize local arrays
  HMat.Local_.colptr.resize(nlocal+1);

 
  //read from 4th line
  int colptrWidth = 0;
  int rowindWidth = 0;
  int nzvalWidth = 0;
  int colptrCntPerRow = 0;
  int rowindCntPerRow = 0;
  int nzvalCntPerRow = 0;
  if(getline(infile, line))
  {
    iss.str("");
    iss.clear();
    iss<<line;
    string format;
    iss>>format;
    int dummy;
    sscanf(format.c_str(),"(%dI%d)",&colptrCntPerRow,&colptrWidth);
    iss>>format;
    sscanf(format.c_str(),"(%dI%d)",&rowindCntPerRow,&rowindWidth);
    iss>>format;
    sscanf(format.c_str(),"(%dE%d.%d)",&nzvalCntPerRow,&nzvalWidth,&dummy);
  }

  //now compute the actual number of rows
  //colptr
  {
    size_t curPos = infile.tellg();
    size_t lineLastNode = std::ceil(( double(firstNode+nlocal) / double(colptrCntPerRow)));
    size_t lineFirstNode = std::ceil(( double(firstNode) / double(colptrCntPerRow)));
    size_t skip = (firstNode - 1)*colptrWidth + (lineFirstNode - 1);
    size_t readBytes = (nlocal+1)*colptrWidth + (lineLastNode - lineFirstNode);
    size_t skipAfter = (n+1 - (firstNode+nlocal))*colptrWidth + (colptrCnt - lineLastNode +1) ;

    infile.seekg(skip,ios_base::cur);

    {
      std::string rdStr;
      rdStr.resize(readBytes);

      infile.read(&rdStr[0], readBytes);

      istringstream iss(rdStr);
      Ptr j;
      Idx locPos = 0;
      while(iss>> j){
        HMat.Local_.colptr[locPos++]=j;
      }
    }

    infile.seekg(skipAfter,ios_base::cur);
    size_t curEnd = infile.tellg();
  }

  //convert to local colptr and compute local nnz
  Ptr first_idx = HMat.Local_.colptr.front();
  Ptr last_idx = HMat.Local_.colptr.back();
  for(Idx i=nlocal+1;i>=1;i--){
    HMat.Local_.colptr[i-1] = HMat.Local_.colptr[i-1] - HMat.Local_.colptr[0] + 1;
  }
  Ptr nnzLocal = HMat.Local_.colptr.back()-1;
  
  HMat.Local_.rowind.resize(nnzLocal);

  Ptr elem_idx;
  //rowind
  {
    size_t curPos = infile.tellg();
    size_t lineLastEdge = std::ceil(( double(last_idx-1) / double(rowindCntPerRow)));
    size_t lineFirstEdge = std::ceil(( double(first_idx) / double(rowindCntPerRow)));
    size_t skip = (first_idx - 1)*rowindWidth + (lineFirstEdge - 1);
    size_t readBytes = (last_idx - first_idx)*rowindWidth + (lineLastEdge - lineFirstEdge);
    size_t skipAfter = (nnz+1 - last_idx)*rowindWidth + (rowindCnt - lineLastEdge +1) ;

    infile.seekg(skip,ios_base::cur);

    {
      std::string rdStr;
      rdStr.resize(readBytes);

      infile.read(&rdStr[0], readBytes);
      istringstream iss(rdStr);
      Idx j;
      Ptr locPos = 0;
      while(iss>> j){
        HMat.Local_.rowind[locPos++]=j;
      }
    }

    infile.seekg(skipAfter,ios_base::cur);
    size_t curEnd = infile.tellg();
  }


  HMat.nzvalLocal.resize(nnzLocal);

  //nzval
  {
    size_t curPos = infile.tellg();
    size_t lineLastEdge = std::ceil(( double(last_idx-1) / double(nzvalCntPerRow)));
    size_t lineFirstEdge = std::ceil(( double(first_idx) / double(nzvalCntPerRow)));
    size_t skip = (first_idx - 1)*nzvalWidth + (lineFirstEdge - 1);
    size_t readBytes = (last_idx - first_idx)*nzvalWidth + (lineLastEdge - lineFirstEdge);
    size_t skipAfter = (nnz+1 - last_idx)*nzvalWidth + (nzvalCnt - lineLastEdge +1) ;

    infile.seekg(skip,ios_base::cur);

    {
      std::string rdStr;
      rdStr.resize(readBytes);

      infile.read(&rdStr[0], readBytes);
//      logfileptr->OFS()<<"nzval read string is"<<endl<<rdStr<<endl;

      istringstream iss(rdStr);
      INSCALAR j;
      Ptr locPos = 0;
      while(iss>> j){
        HMat.nzvalLocal[locPos++]=(SCALAR)j;
      }
    }

    infile.seekg(skipAfter,ios_base::cur);
    size_t curEnd = infile.tellg();
  }


  infile.close();


  HMat.globalAllocated = false;
  HMat.Local_.size = HMat.size;
  HMat.Local_.nnz = nnzLocal;

  return 0;
}


template <>
int ReadHB_PARA<double,double>(std::string & filename, DistSparseMatrix<double> & HMat){

//logfileptr->OFS()<<HMat.Local_.colptr<<endl;
//logfileptr->OFS()<<HMat.Local_.rowind<<endl;
//logfileptr->OFS()<<HMat.nzvalLocal<<endl;

  MPI_Comm & workcomm = HMat.comm;

    int mpirank;
    MPI_Comm_rank(workcomm,&mpirank);
  
    int mpisize;
    MPI_Comm_size(workcomm,&mpisize);

  ifstream infile;
  infile.open(filename.c_str());

  string line;
  //read xadj on the first line of the input file
  stringstream iss;
  //skip 1st line
  if(getline(infile, line)){}
  int64_t colptrCnt;
  int64_t rowindCnt;
  int64_t nzvalCnt;
  if(getline(infile, line)){
    iss.str("");
    iss.clear();
    iss<<line;
    int64_t dummy;
    iss>>dummy;
    iss>>colptrCnt>>rowindCnt>>nzvalCnt;
  }
  //read from third line

  auto & n = HMat.size;
  auto & nnz = HMat.nnz;


  int64_t m;
  if(getline(infile, line))
  {
    iss.str("");
    iss.clear();
    iss<<line;
    string type;
    iss>>type;
    iss>>m>>n>>nnz;
  }

  //compute local number of columns
  Idx nlocal = (mpirank<mpisize-1)?n/mpisize:n-mpirank*(int64_t)(n/mpisize);
  Idx firstNode = mpirank*(int64_t)(n/mpisize) + 1;
  //initialize local arrays
  HMat.Local_.colptr.resize(nlocal+1);

 
  //read from 4th line
  int colptrWidth = 0;
  int rowindWidth = 0;
  int nzvalWidth = 0;
  int colptrCntPerRow = 0;
  int rowindCntPerRow = 0;
  int nzvalCntPerRow = 0;
  if(getline(infile, line))
  {
    iss.str("");
    iss.clear();
    iss<<line;
    string format;
    iss>>format;
    int dummy;
    sscanf(format.c_str(),"(%dI%d)",&colptrCntPerRow,&colptrWidth);
    iss>>format;
    sscanf(format.c_str(),"(%dI%d)",&rowindCntPerRow,&rowindWidth);
    iss>>format;
    sscanf(format.c_str(),"(%dE%d.%d)",&nzvalCntPerRow,&nzvalWidth,&dummy);
  }

  //now compute the actual number of rows
  //colptr
  {
    size_t curPos = infile.tellg();
    size_t lineLastNode = std::ceil(( double(firstNode+nlocal) / double(colptrCntPerRow)));
    size_t lineFirstNode = std::ceil(( double(firstNode) / double(colptrCntPerRow)));
    size_t skip = (firstNode - 1)*colptrWidth + (lineFirstNode - 1);
    size_t readBytes = (nlocal+1)*colptrWidth + (lineLastNode - lineFirstNode);
    size_t skipAfter = (n+1 - (firstNode+nlocal))*colptrWidth + (colptrCnt - lineLastNode +1) ;

    infile.seekg(skip,ios_base::cur);

    {
      std::string rdStr;
      rdStr.resize(readBytes);

      infile.read(&rdStr[0], readBytes);
//      logfileptr->OFS()<<"read string is"<<endl<<rdStr<<endl;

      istringstream iss(rdStr);
      Ptr j;
      Idx locPos = 0;
      while(iss>> j){
        HMat.Local_.colptr[locPos++]=j;
      }
    }

    infile.seekg(skipAfter,ios_base::cur);
    size_t curEnd = infile.tellg();
  }


  //convert to local colptr and compute local nnz
  Ptr first_idx = HMat.Local_.colptr.front();
  Ptr last_idx = HMat.Local_.colptr.back();
  for(Idx i=nlocal+1;i>=1;i--){
    HMat.Local_.colptr[i-1] = HMat.Local_.colptr[i-1] - HMat.Local_.colptr[0] + 1;
  }
  Ptr nnzLocal = HMat.Local_.colptr.back()-1;
  
  HMat.Local_.rowind.resize(nnzLocal);

  Ptr elem_idx;
  //rowind
  {
    size_t curPos = infile.tellg();
    size_t lineLastEdge = std::ceil(( double(last_idx-1) / double(rowindCntPerRow)));
    size_t lineFirstEdge = std::ceil(( double(first_idx) / double(rowindCntPerRow)));
    size_t skip = (first_idx - 1)*rowindWidth + (lineFirstEdge - 1);
    size_t readBytes = (last_idx - first_idx)*rowindWidth + (lineLastEdge - lineFirstEdge);
    size_t skipAfter = (nnz+1 - last_idx)*rowindWidth + (rowindCnt - lineLastEdge +1) ;

    infile.seekg(skip,ios_base::cur);

    {
      std::string rdStr;
      rdStr.resize(readBytes);

      infile.read(&rdStr[0], readBytes);
//      logfileptr->OFS()<<"rowind read string is"<<endl<<rdStr<<endl;
      istringstream iss(rdStr);
      Idx j;
      Ptr locPos = 0;
      while(iss>> j){
        HMat.Local_.rowind[locPos++]=j;
      }
    }

    infile.seekg(skipAfter,ios_base::cur);
    size_t curEnd = infile.tellg();
  }



  HMat.nzvalLocal.resize(nnzLocal);

  //nzval
  {
    size_t curPos = infile.tellg();
    size_t lineLastEdge = std::ceil(( double(last_idx-1) / double(nzvalCntPerRow)));
    size_t lineFirstEdge = std::ceil(( double(first_idx) / double(nzvalCntPerRow)));
    size_t skip = (first_idx - 1)*nzvalWidth + (lineFirstEdge - 1);
    size_t readBytes = (last_idx - first_idx)*nzvalWidth + (lineLastEdge - lineFirstEdge);
    size_t skipAfter = (nnz+1 - last_idx)*nzvalWidth + (nzvalCnt - lineLastEdge +1) ;

    infile.seekg(skip,ios_base::cur);

    {
      std::string rdStr;
      Ptr locPos = 0;
      ////Read MB per MB
      //size_t already_read = 0;
      //while(already_read<readBytes){
      //  size_t toRead = std::min(readBytes-already_read,(size_t)1024*1024);
      //  rdStr.resize(toRead);
      //  infile.read(&rdStr[0], toRead);

        rdStr.resize(readBytes);
        infile.read(&rdStr[0], readBytes);

        istringstream iss(rdStr);
        double j;
        while(iss>> j){
          HMat.nzvalLocal[locPos++]=(double)j;
        }
      //  already_read+=toRead;
      //}
    }

    infile.seekg(skipAfter,ios_base::cur);
    size_t curEnd = infile.tellg();
  }

  infile.close();

  HMat.globalAllocated = false;
  HMat.Local_.size = HMat.size;
  HMat.Local_.nnz = nnzLocal;


  return 0;

}
#endif















}  // namespace SYMPACK
