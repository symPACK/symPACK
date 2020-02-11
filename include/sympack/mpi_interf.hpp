#ifndef _SYMPACK_MPI_HPP_
#define _SYMPACK_MPI_HPP_

#include <iostream>

#include  "sympack/Environment.hpp"
#include  "sympack/CommTypes.hpp"
#include  "sympack/utility.hpp"

#include <mpi.h>
#include <functional>
#include <numeric>
#include <limits>

namespace symPACK{

  template<typename T> void GenericMPIMax( void *in, void *inout, int *len, MPI_Datatype *dptr );

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
        int mpirank, mpisize;
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
      }		// -----  end of function Gatherv  ----- 








template<typename _Container, typename _Size>
inline int Alltoallv(_Container & sendbuf, const _Size *stotcounts, const _Size *stotdispls, MPI_Datatype sendtype,
                _Container & recvbuf, _Size *rtotcounts, _Size *rtotdispls,MPI_Comm & comm, std::function< void(_Container &, size_t)> & resize_func){

  int mpirank, mpisize, err;
  MPI_Comm_rank( comm, &mpirank );
  MPI_Comm_size( comm, &mpisize );
  err = MPI_SUCCESS;

  //if _Size is signed, we can save an alltoall, reducing the overhead if there is an error
  bool is_signed = std::numeric_limits<_Size>::is_signed;
  
  //first check that everything is positive
  if ( is_signed &&  // cast to intmax avoids compile warnings when _Size is unsigned
      (std::any_of(stotcounts,stotcounts+mpisize,[](const _Size & a){ return std::intmax_t(a) < 0;}) || 
                std::any_of(stotdispls,stotdispls+mpisize,[](const _Size & a){ return std::intmax_t(a) < 0;}) ) ) {
    return -1;
  }

  //convert to size_t
  std::vector< std::size_t > stotsizes(mpisize);
  for ( std::size_t i = 0; i < mpisize; i++ ) { 
    stotsizes[i] = std::size_t(stotcounts[i]);
  }

  //now check if everything is below threshold ( limit of integers )
  bool need_split = false;
  need_split = std::any_of(stotsizes.begin(),stotsizes.end(),
        [](const std::size_t & a){ return a > std::numeric_limits<int>::max();});

  
  if ( !need_split ) {
    std::vector< std::size_t > stotdispls_test(mpisize+1);
    stotdispls_test[0] = 0;
    std::partial_sum(stotsizes.begin(),stotsizes.end(),&stotdispls_test[1]);

    
    for ( std::size_t i = 0; i < mpisize+1; i++ ) {
      if ( stotdispls_test[i] / std::numeric_limits<int>::max() > 0 ) {
        need_split = true;
      }
    }
  }

  //now do an Allgather of the flag
  std::vector<char> need_splits(mpisize);
  err = MPI_Allgather(&need_split,sizeof(bool),MPI_BYTE,need_splits.data(),sizeof(bool),MPI_BYTE,comm);
  if ( err != MPI_SUCCESS ) return err;
  need_split = std::any_of(need_splits.begin(),need_splits.end(),
        [](const char & a){ return (bool)a;});


  size_t myGCD = 1;
  if ( need_split ) {

    //try to see if transferring larger chunks makes the error go away
    myGCD = findGCD(stotsizes.data(),stotsizes.size());
    //now we need to to an allgather of these
    std::vector<size_t> gcds(mpisize,0);
    err = MPI_Allgather(&myGCD,sizeof(size_t),MPI_BYTE,gcds.data(),sizeof(size_t),MPI_BYTE,comm);
    if ( err != MPI_SUCCESS ) return err;
    myGCD = findGCD(gcds.data(),gcds.size());

    need_split = false;
    for ( std::size_t i = 0; i < mpisize; i++ ) { 
      if (std::size_t(stotcounts[i])/std::size_t(myGCD) > std::numeric_limits<int>::max() ) {
        need_split = true;
        break;
      }
      else {
        stotsizes[i] = std::size_t(stotcounts[i])/std::size_t(myGCD);
      }
    }

    if ( !need_split ) {
      std::vector< std::size_t > stotdispls_test(mpisize+1);
      stotdispls_test[0] = 0;
      std::partial_sum(stotsizes.begin(),stotsizes.end(),&stotdispls_test[1]);

    
      for ( std::size_t i = 0; i < mpisize+1; i++ ) {
        if ( stotdispls_test[i] / std::numeric_limits<int>::max() > 0 ) {
          need_split = true;
          break;
        }
      }

      //now do an Allgather of the flag
      std::vector<char> need_splits(mpisize);
      err = MPI_Allgather(&need_split,sizeof(bool),MPI_BYTE,need_splits.data(),sizeof(bool),MPI_BYTE,comm);
      if ( err != MPI_SUCCESS ) return err;
      need_split = std::any_of(need_splits.begin(),need_splits.end(),
                                    [](const char & a){ return (bool)a;});
    }
  }

  if ( need_split ) {
    //TODO handle that later
    gdb_lock();
  }
  else {

    //get the byte size from the MPI_Datatype to know how much we need to advance the pointers
    int typesize = 0;
    MPI_Type_size(sendtype, &typesize); 

    MPI_Datatype fused_type;
    MPI_Type_contiguous( typesize*myGCD, MPI_BYTE, &fused_type );
    MPI_Type_commit(&fused_type);

    vector<int> ssizes2(mpisize,0);
    for (int i = 0; i< mpisize; i++) { ssizes2[i] = int( std::size_t(stotcounts[i]) / myGCD );}

    vector<int> sdispls2(mpisize+1,0);
    sdispls2[0] = 0;
    std::partial_sum(ssizes2.begin(),ssizes2.end(),&sdispls2[1]);

    vector<int> rsizes2(mpisize,0);
    err = MPI_Alltoall(&ssizes2[0],sizeof(int),MPI_BYTE,&rsizes2[0],sizeof(int),MPI_BYTE,comm);
    if ( err != MPI_SUCCESS ) return err;

    //compute receive displacements
    vector<int> rdispls2(mpisize+1,0);
    rdispls2[0] = 0;
    std::partial_sum(rsizes2.begin(),rsizes2.end(),&rdispls2[1]);

    //resize the receive container
    resize_func(recvbuf,rdispls2.back()*myGCD);

    //Now do the alltoallv
    err = MPI_Alltoallv(&sendbuf[0],&ssizes2[0],&sdispls2[0],fused_type,&recvbuf[0],&rsizes2[0],&rdispls2[0],fused_type,comm);
    if ( err != MPI_SUCCESS ) return err;

    MPI_Type_free(&fused_type);

    for (int i = 0; i< mpisize; i++) { rtotcounts[i] = _Size( std::size_t(rsizes2[i]) * myGCD );}
    for (int i = 0; i< mpisize+1; i++) { rtotdispls[i] = _Size( std::size_t(rdispls2[i]) * myGCD );}
  }

  return err;
}


















template<typename _Container, typename _Size>
inline int Alltoallv2(_Container & sendbuf, const _Size *stotcounts, const _Size *stotdispls, MPI_Datatype sendtype,
                _Container & recvbuf, _Size *rtotcounts, _Size *rtotdispls,MPI_Comm & comm, std::function< void(_Container &, size_t)> & resize_func){

  int mpirank, mpisize;
  MPI_Comm_rank( comm, &mpirank );
  MPI_Comm_size( comm, &mpisize );

  //if _Size is signed, we can save an alltoall, reducing the overhead if there is an error
  bool is_signed = std::numeric_limits<_Size>::is_signed;
  


  //first compare the localMax to the maximum chunk size. If something is bigger, set everything to -1 
  int maxchunk = std::numeric_limits<int>::max()/mpisize;
  const _Size * maxLocal = std::max_element(&stotcounts[0],&stotcounts[0]+mpisize);
  bool need_split = *maxLocal>maxchunk;


  MPI_Datatype size_type;
  MPI_Type_contiguous( sizeof(_Size), MPI_BYTE, &size_type );
  MPI_Type_commit(&size_type);
  //do the all to all


  if(need_split){
    std::vector<_Size> serrorsizes;
    if(is_signed){
      serrorsizes.reserve(mpisize);
      for(int i = 0; i<mpisize;i++){ serrorsizes[i] = -stotcounts[i]; }
    }
    else{
      serrorsizes.assign(mpisize,std::numeric_limits<_Size>::max());
    }
    MPI_Alltoall(&serrorsizes[0],1,size_type,&rtotcounts[0],1,size_type,comm);
  }
  else{
    MPI_Alltoall(&stotcounts[0],1,size_type,&rtotcounts[0],1,size_type,comm);
  }

  if(is_signed){
    const _Size * errorLocal = std::min_element(&rtotcounts[0],&rtotcounts[0]+mpisize);
    need_split = *errorLocal<0;
  }
  else{
    const _Size * errorLocal = std::max_element(&rtotcounts[0],&rtotcounts[0]+mpisize);
    need_split = *errorLocal==std::numeric_limits<_Size>::max();
  }










  size_t max_sr_size = 0;
  if(need_split){
    _Size globalMax = 0;
    MPI_Op MPI_SYMPACK_MAX; 
    MPI_Op_create( GenericMPIMax<_Size>, true, &MPI_SYMPACK_MAX ); 
    MPI_Allreduce(maxLocal,&globalMax,1,size_type,MPI_SYMPACK_MAX,comm);
    MPI_Op_free(&MPI_SYMPACK_MAX);
    logfileptr->OFS()<<"mpi::Alltoallv ERROR has been reported by one of the nodes, reduced global max is "<<globalMax<<std::endl;
    max_sr_size = std::max(max_sr_size,globalMax);

    if(is_signed){
        for(int i = 0; i<mpisize;i++){ rtotcounts[i] = rtotcounts[i]>0?rtotcounts[i]:-rtotcounts[i]; }
    }
    else{
      //redo the alltoall to get rid of the error flags
      MPI_Alltoall(&stotcounts[0],1,size_type,&rtotcounts[0],1,size_type,comm);
    }
  }

  MPI_Type_free(&size_type);










  rtotdispls[0] = 0;
  std::partial_sum(rtotcounts,rtotcounts+mpisize-1,&rtotdispls[1]);

  size_t total_send_size = 0;
  size_t total_recv_size = 0;
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
    logfileptr->OFS()<<"mpi::Alltoallv collective will be split in "<<split_atav<<" calls of size "<<chunk_atav<<std::endl;

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
          scounts[i] = int(std::min(_Size(chunk_atav),stotcounts[i]-offset_buffer));
        }
        else{
          scounts[i]=0;
        }
        tot_send+=scounts[i];
      }

      size_t tot_recv = 0;
      for(Int i = 0; i< mpisize; i++){
        if(offset_buffer<=size_t(rtotcounts[i])){
          rcounts[i] = int(std::min(_Size(chunk_atav),rtotcounts[i]-offset_buffer));
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


      int errmpi = MPI_Alltoallv((void*)sptrs, &scounts[0], &sdispls[0], sendtype,
          (void*)rptrs, &rcounts[0], &rdispls[0], sendtype,
          comm);


      bassert( scounts[mpirank] == rcounts[mpirank] );
      for ( size_t i = 0; i < scounts[mpirank] ; i++) {
        bassert( sptrs[ sdispls[mpirank] + i ] == rptrs[ rdispls[mpirank] + i ] );
      }



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

  rtotdispls[0] = 0;
  std::partial_sum(rtotcounts,rtotcounts+mpisize,&rtotdispls[1]);
  return 0;
}



template<typename _Container, typename _Size>
inline int Alltoallv(_Container & sendbuf, const _Size *stotcounts, const _Size *stotdispls, MPI_Datatype sendtype,
                _Container & recvbuf, MPI_Comm & comm, std::function< void(_Container &, size_t)> & resize_func){

  int mpisize;
  MPI_Comm_size( comm, &mpisize );
  std::vector<_Size> rtotdispls(mpisize+1,0);
  std::vector<_Size> rtotcounts(mpisize,0);
  return Alltoallv(sendbuf, stotcounts, stotdispls, sendtype,
                recvbuf, rtotcounts.data(), rtotdispls.data(), comm, resize_func);
}





template<typename T>
inline void Allreduce( T* sendbuf, T* recvbuf, Int count, MPI_Op op, MPI_Comm comm){}

template<>
inline void Allreduce<double>( double* sendbuf, double* recvbuf, Int count, MPI_Op op, MPI_Comm comm){
      MPI_Allreduce( sendbuf, recvbuf, count, MPI_DOUBLE, op, comm );
}

template<>
inline void Allreduce<float>( float* sendbuf, float* recvbuf, Int count, MPI_Op op, MPI_Comm comm){
      MPI_Allreduce( sendbuf, recvbuf, count, MPI_FLOAT, op, comm );
}

template<>
inline void Allreduce<std::complex<float> >( std::complex<float>* sendbuf, std::complex<float>* recvbuf, Int count, MPI_Op op, MPI_Comm comm){
      MPI_Allreduce( sendbuf, recvbuf, count, MPI_COMPLEX, op, comm );
}

template<>
inline void Allreduce<std::complex<double> >( std::complex<double>* sendbuf, std::complex<double>* recvbuf, Int count, MPI_Op op, MPI_Comm comm){
      MPI_Allreduce( sendbuf, recvbuf, count, MPI_DOUBLE_COMPLEX, op, comm );
}

  }
}



#endif //_SYMPACK_MPI_HPP_

