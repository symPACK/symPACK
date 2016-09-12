#ifndef _SYMPACK_MPI_HPP_
#define _SYMPACK_MPI_HPP_

#include <iostream>

#include  "sympack/Environment.hpp"
#include  "sympack/CommTypes.hpp"

#include <mpi.h>
#include <functional>
#include <numeric>
#include <limits>

namespace SYMPACK{

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
        SYMPACK::vector<Int>  localSizeVec( mpisize );
        SYMPACK::vector<Int>  localSizeDispls( mpisize );
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



template<typename _Container, typename _Size>
int Alltoallv(_Container & sendbuf, const _Size *stotcounts, const _Size *stotdispls, MPI_Datatype sendtype,
                _Container & recvbuf, MPI_Comm & comm, std::function< void(_Container &, size_t)> & resize_func){
 
//gdb_lock();
  Int mpirank, mpisize;
  MPI_Comm_rank( comm, &mpirank );
  MPI_Comm_size( comm, &mpisize );

  MPI_Datatype size_type;
  MPI_Type_contiguous( sizeof(_Size), MPI_BYTE, &size_type );
  MPI_Type_commit(&size_type);
  //do the all to all
  std::vector<_Size> rtotcounts(mpisize,0);
  MPI_Alltoall(&stotcounts[0],1,size_type,&rtotcounts[0],1,size_type,comm);
  MPI_Type_free(&size_type);

  std::vector<_Size> rtotdispls(mpisize,0);
  rtotdispls[0] = 0;
  std::partial_sum(rtotcounts.begin(),rtotcounts.end()-1,&rtotdispls[1]);

  int maxchunk = 1300;//std::numeric_limits<int>::max()/mpisize;
  size_t total_send_size = 0;
  size_t total_recv_size = 0;
  size_t max_sr_size = 0;
  size_t max_s_gap = 0;
  size_t max_r_gap = 0;
  for(Int i = 0; i< mpisize; i++){
    max_sr_size = std::max(max_sr_size,std::max(stotcounts[i],rtotcounts[i]));
    max_s_gap = std::max(max_s_gap, stotcounts[i]>maxchunk?stotcounts[i] - maxchunk:0);
    max_r_gap = std::max(max_r_gap, rtotcounts[i]>maxchunk?rtotcounts[i] - maxchunk:0);

    total_send_size += size_t(stotcounts[i]);
    total_recv_size += size_t(rtotcounts[i]);
  }

  if(total_send_size>0 || total_recv_size>0){
    int chunk_atav= int(std::min(max_sr_size,size_t(maxchunk)));
    int split_atav = (int)std::ceil((double)max_sr_size/(double)chunk_atav);

    //resize the receive container
    resize_func(recvbuf,total_recv_size);

    //declare pointers because the offsets and displacements
    // can only be expressed in int, so we need to move the pointers instead
    //and may need to make copies if maximum gap is larger than maxint
    char* sptrs = (char*)&sendbuf[0];
    char* rptrs = (char*)&recvbuf[0];

    //get the byte size from the MPI_Datatype to know how much we need to advance the pointers
    int typesize = 0;
    MPI_Type_size(sendtype, &typesize); 


    std::vector<char> * tmpsbuf; 
    std::vector<char> * tmprbuf; 
    if(max_s_gap>0){
      tmpsbuf = new std::vector<char>();
    }

    if(max_r_gap>0){
      tmprbuf = new std::vector<char>();
    }

    std::vector<_Size> stotpositions(mpisize);
    std::copy(&stotdispls[0],&stotdispls[0]+mpisize,&stotpositions[0]);
    std::vector<_Size> rtotpositions(mpisize);
    std::copy(&rtotdispls[0],&rtotdispls[0]+mpisize,&rtotpositions[0]);

    std::vector<int> sdispls(mpisize,0);
    std::vector<int> scounts(mpisize,0);
    std::vector<int> rdispls(mpisize,0);
    std::vector<int> rcounts(mpisize,0);


    //now we know how many alltoallv are required
    for(int step_atav = 0; step_atav<split_atav;step_atav++){
      size_t offset_buffer = step_atav*chunk_atav;

      size_t tot_send = 0;
      for(Int i = 0; i< mpisize; i++){
        if(offset_buffer<=size_t(stotcounts[i])){
          scounts[i] = std::min(chunk_atav,int(stotcounts[i]-offset_buffer));
        }
        else{
          scounts[i]=0;
        }
        tot_send+=scounts[i];
      }

      size_t tot_recv = 0;
      for(Int i = 0; i< mpisize; i++){
        if(offset_buffer<=size_t(rtotcounts[i])){
          rcounts[i] = std::min(chunk_atav,int(rtotcounts[i]-offset_buffer));
        }
        else{
          rcounts[i]=0;
        }
        tot_recv+=rcounts[i];
      }


      if(max_s_gap>0){
        //compute send displacements
        sdispls[0] = 0;
        std::partial_sum(scounts.begin(),scounts.end()-1,&sdispls[1]);
        //make a copy of the appropriate chunk in the temp buffer
        tmpsbuf->resize(typesize*tot_send);
        for(Int i = 0; i< mpisize; i++){
          if(scounts[i]>0){
            std::copy(&sendbuf[stotpositions[i]],  &sendbuf[stotpositions[i]]+scounts[i],&(tmpsbuf->at(sdispls[i]*typesize)));
          }
        }
        //set the pointer
        sptrs = &tmpsbuf->front();
      }
      else{
        //compute send displacements
        sdispls[0] = 0;
        for(Int i = 1; i< mpisize; i++){
          if(scounts[i]>0){
            sdispls[i] = ((char*)&sendbuf[stotpositions[i]]-sptrs)/typesize;
          }
          else{
            sdispls[i] = sdispls[i-1];
          }
        }
      }

      if(max_r_gap>0){
        //compute recv displacements
        rdispls[0] = 0;
        std::partial_sum(rcounts.begin(),rcounts.end()-1,&rdispls[1]);
        //make a copy of the appropriate chunk in the temp buffer
        tmprbuf->resize(typesize*tot_recv);
        //set the pointer
        rptrs = &tmprbuf->front();
      }
      else{
        //compute recv displacements
        rdispls[0] = 0;
        for(Int i = 1; i< rcounts.size(); i++){
          if(rcounts[i]>0){
            rdispls[i] = ((char*)&recvbuf[rtotpositions[i]]-rptrs)/typesize;
          }
          else{
            rdispls[i] = rdispls[i-1];
          }
        }
      }


      MPI_Alltoallv((void*)sptrs, &scounts[0], &sdispls[0], sendtype,
          (void*)rptrs, &rcounts[0], &rdispls[0], sendtype,
          comm);

      if(max_s_gap==0){
        //advance the pointer
        sptrs += typesize*scounts[0]; 
      }

      if(max_r_gap==0){
        //advance the pointer
        rptrs += typesize*rcounts[0]; 
      }
      else{
        //need to copy back in the main recv buffer
        for(Int i = 0; i< mpisize; i++){
          if(rcounts[i]>0){
            std::copy(&(tmprbuf->at(rdispls[i]*typesize)),&(tmprbuf->at(rdispls[i]*typesize))+rcounts[i]*typesize, (char*)&recvbuf[rtotpositions[i]]);
          }
        }
      }

      //advance the send and recv positions
      for(Int i = 0; i< mpisize; i++){
        stotpositions[i]+=scounts[i];
        rtotpositions[i]+=rcounts[i];
      }
    }

    if(max_r_gap>0){
      delete tmprbuf;
    }
    if(max_s_gap>0){
      delete tmpsbuf;
    }
  }
  return 0;
}





template<typename T>
void Allreduce( T* sendbuf, T* recvbuf, Int count, MPI_Op op, MPI_Comm comm){};

template<>
void Allreduce<double>( double* sendbuf, double* recvbuf, Int count, MPI_Op op, MPI_Comm comm){
      MPI_Allreduce( sendbuf, recvbuf, count, MPI_DOUBLE, op, comm );
}

template<>
void Allreduce<float>( float* sendbuf, float* recvbuf, Int count, MPI_Op op, MPI_Comm comm){
      MPI_Allreduce( sendbuf, recvbuf, count, MPI_FLOAT, op, comm );
}

template<>
void Allreduce<std::complex<float> >( std::complex<float>* sendbuf, std::complex<float>* recvbuf, Int count, MPI_Op op, MPI_Comm comm){
      MPI_Allreduce( sendbuf, recvbuf, count, MPI_COMPLEX, op, comm );
}

template<>
void Allreduce<std::complex<double> >( std::complex<double>* sendbuf, std::complex<double>* recvbuf, Int count, MPI_Op op, MPI_Comm comm){
      MPI_Allreduce( sendbuf, recvbuf, count, MPI_DOUBLE_COMPLEX, op, comm );
}

  }
}



#endif

