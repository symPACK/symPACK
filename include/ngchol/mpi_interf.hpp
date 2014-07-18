#ifndef _NGCHOL_MPI_HPP_
#define _NGCHOL_MPI_HPP_

#include <iostream>

#include  "ngchol/Environment.hpp"
#include  "ngchol/CommTypes.hpp"

#include <mpi.h>

namespace LIBCHOLESKY{

  /// @namespace mpi
  ///
  /// @brief Interface with MPI to facilitate communication.
  namespace mpi{

    inline void
      Gatherv ( 
          Icomm& localIcomm, 
          Icomm& allIcomm,
          Int root,
          MPI_Comm          comm )
      {
        Int mpirank, mpisize;
        MPI_Comm_rank( comm, &mpirank );
        MPI_Comm_size( comm, &mpisize );

        Int localSize = localIcomm.size();
        std::vector<Int>  localSizeVec( mpisize );
        std::vector<Int>  localSizeDispls( mpisize );
        MPI_Gather( &localSize, sizeof(Int), MPI_BYTE, &localSizeVec[0], sizeof(Int), MPI_BYTE,root, comm );
        localSizeDispls[0] = 0;
        for( Int ip = 1; ip < mpisize; ip++ ){
          localSizeDispls[ip] = (localSizeDispls[ip-1] + localSizeVec[ip-1]);
        }
        Int totalSize = (localSizeDispls[mpisize-1] + localSizeVec[mpisize-1]);

        allIcomm.clear();
        allIcomm.resize( totalSize );

        MPI_Gatherv( localIcomm.front(), localSize, MPI_BYTE, allIcomm.front(), 
            &localSizeVec[0], &localSizeDispls[0], MPI_BYTE,root, comm	);

        //Mark the Icomm as "full"
        allIcomm.setHead(totalSize);
        return ;
      };		// -----  end of function Gatherv  ----- 


inline void 
Allreduce( Real* sendbuf, Real* recvbuf, Int count, MPI_Op op, MPI_Comm comm)
{
      MPI_Allreduce( sendbuf, recvbuf, count, MPI_DOUBLE, op, comm );
}


inline void 
Allreduce( Complex* sendbuf, Complex* recvbuf, Int count, MPI_Op op, MPI_Comm comm)
{
      MPI_Allreduce( sendbuf, recvbuf, 2*count, MPI_DOUBLE, op, comm );
}


//template<typename T>
//void 
//Allreduce( T* sendbuf, T* recvbuf, Int count, MPI_Op op, MPI_Comm comm)
//{
//      MPI_Allreduce( sendbuf, recvbuf, count*sizeof(T), MPI_BYTE, op, comm );
//}



  }
}



#endif

