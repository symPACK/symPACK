#ifndef _SUPERNODAL_MATRIX_IMPL_HP_
#define _SUPERNODAL_MATRIX_IMPL_HPP_

#include "ngchol/SupernodalMatrix.hpp"

#include "ngchol/blas.hpp"
#include "ngchol/Ordering.hpp"
#include "ngchol/mpi_interf.hpp"

#include <queue>


#define BLOCKSIZE src_snode.Size()

#define TAG_INDEX 0
#define TAG_NZVAL 1
#define TAG_COUNT 2



namespace LIBCHOLESKY{





  template <typename T> void SupernodalMatrix<T>::Init(const DistSparseMatrix<T> & pMat, Int maxSnode,Mapping * pMapping, Int maxIsend, Int maxIrecv, MPI_Comm & pComm ){
  //Create the CommEnvironment object if necessary
    globalAllocated = false;
    if(CommEnv_!=NULL){
      delete CommEnv_;
    }    
    CommEnv_ = new CommEnvironment(pComm);

    Mapping_ = pMapping;

    //Options
    maxIsend_ = maxIsend;
    maxIrecv_ = maxIrecv;

    iSize_ = pMat.size;
    Local_ = pMat.GetLocalStructure();

    Local_.ToGlobal(Global_);
    Global_.ExpandSymmetric();
    //Create an Ordering object to hold the permutation
    Order_.SetStructure(Global_);


    //Reoder the matrix with MMD
    Order_.MMD();

//logfileptr->OFS()<<"Order.perm "<<Order_.perm<<endl;
//logfileptr->OFS()<<"Order.invp "<<Order_.invp<<endl;



    ETree_.ConstructETree(Global_,Order_);

//logfileptr->OFS()<<"ETree "<<ETree_<<endl;
    ETree_.PostOrderTree(Order_);
//logfileptr->OFS()<<"POETree "<<ETree_<<endl;

//logfileptr->OFS()<<"Order.perm "<<Order_.perm<<endl;
//logfileptr->OFS()<<"Order.invp "<<Order_.invp<<endl;




    IntNumVec cc,rc;
    Global_.GetLColRowCount(ETree_,Order_,cc,rc);
//    logfileptr->OFS()<<"colcnt "<<cc<<std::endl;
    ETree_.SortChildren(cc,Order_);

//logfileptr->OFS()<<"Sorted POETree "<<ETree_<<endl;
//    logfileptr->OFS()<<"sorted colcnt "<<cc<<std::endl;
//logfileptr->OFS()<<"Order.perm "<<Order_.perm<<endl;
//logfileptr->OFS()<<"Order.invp "<<Order_.invp<<endl;




#ifdef _DEBUG_
    logfileptr->OFS()<<"colcnt "<<cc<<std::endl;
    logfileptr->OFS()<<"rowcnt "<<rc<<std::endl;
#endif
double flops = 0.0;
for(Int i = 0; i<cc.m();++i){
  flops+= (double)pow((double)cc[i],2.0);
}

if(iam==0){
  cout<<"Flops: "<<flops<<endl;
}

    Global_.FindSupernodes(ETree_,Order_,cc,SupMembership_,Xsuper_,maxSnode);


#ifdef RELAXED_SNODE
    Global_.RelaxSupernodes(ETree_, cc,SupMembership_, Xsuper_, maxSnode );
    Global_.SymbolicFactorizationRelaxed(ETree_,Order_,cc,Xsuper_,SupMembership_,xlindx_,lindx_);

#else
    Global_.SymbolicFactorization(ETree_,Order_,cc,Xsuper_,SupMembership_,xlindx_,lindx_);
#endif



#ifdef REFINED_SNODES
    IntNumVec permRefined;
//    IntNumVec newPerm(Size());
    Global_.RefineSupernodes(ETree_,Order_, SupMembership_, Xsuper_, xlindx_, lindx_, permRefined);

//      Perm_.Resize(Size());
//      for(Int i =0; i<Perm_.m();++i){
//        Perm_[permRefined[i]-1] = permChild[i];
//      }

//    Perm_ = permRefined;
#else
    //Combine permMMD with Postorder
//    Perm_.Resize(Size());
//    for(Int i =0; i<Perm_.m();++i){
//      Perm_[i] = ETree_.FromPostOrder(i+1);
//    }
#endif

#ifdef _DEBUG_
    logfileptr->OFS()<<"Membership list is "<<SupMembership_<<std::endl;
    logfileptr->OFS()<<"xsuper "<<Xsuper_<<std::endl;
#endif

    GetUpdatingSupernodeCount(UpdateCount_,UpdateWidth_);




    Int iam = CommEnv_->MPI_Rank();
    Int np  = CommEnv_->MPI_Size();

    DistSparseMatrix<T> ExpA(CommEnv_->MPI_GetComm());
{
double timeSta = get_time();
    ExpA.size = pMat.size;
    ExpA.nnz = 2*pMat.nnz-pMat.size;
	  Int numColFirst = std::max(1,ExpA.size / np);
	  std::vector<Int> numColLocalVec(np,numColFirst);
    numColLocalVec[np-1] = ExpA.size - numColFirst * (np-1);  // Modify the last entry	
  	Int numColLocal = numColLocalVec[iam];
    //Expand A to symmetric storage
    Int localFirstCol = iam*numColFirst+1;
    Int localLastCol = localFirstCol+numColLocal-1;

    //Initialize the Local structure
    ExpA.Local_.size = ExpA.size;
    ExpA.Local_.colptr.Resize(numColLocal+1);

	  for( Int i = 0; i < numColLocal + 1; i++ ){
		  ExpA.Local_.colptr[i] = Global_.expColptr[iam * numColFirst+i] - Global_.expColptr[iam * numColFirst] + 1;
  	}

    ExpA.Local_.nnz = ExpA.Local_.colptr[numColLocal] - ExpA.Local_.colptr[0];

    ExpA.Local_.rowind.Resize(ExpA.Local_.nnz);

    Int globalColBeg = Global_.expColptr[iam*numColFirst];
    std::copy(&Global_.expRowind[globalColBeg-1],&Global_.expRowind[globalColBeg-1]+ExpA.Local_.nnz,ExpA.Local_.rowind.Data());
    ExpA.nzvalLocal.Resize(ExpA.Local_.nnz);

#ifdef _DEBUG_
logfileptr->OFS()<<"Global_.colptr: "<<Global_.colptr<<endl;
logfileptr->OFS()<<"Global_.rowind: "<<Global_.rowind<<endl;
logfileptr->OFS()<<"Global_.expColptr: "<<Global_.expColptr<<endl;
logfileptr->OFS()<<"Global_.expRowind: "<<Global_.expRowind<<endl;
logfileptr->OFS()<<"ExpA.colptr: "<<ExpA.Local_.colptr<<endl;
logfileptr->OFS()<<"ExpA.rowind: "<<ExpA.Local_.rowind<<endl;
logfileptr->OFS()<<"pMat.colptr: "<<pMat.Local_.colptr<<endl;
logfileptr->OFS()<<"pMat.rowind: "<<pMat.Local_.rowind<<endl;
#endif

    IntNumVec localColHead = pMat.Local_.colptr;
    IntNumVec localDestColHead = ExpA.Local_.colptr;
    
    NumVec<T> recvNzval;
    IntNumVec recvColptr;
    IntNumVec recvRowind;
    for(Int proc = 0; proc<np; ++proc){
      //communicate colptr
      //Broadcast size
      Int size = iam==proc?pMat.Local_.colptr.m():0;
      MPI_Bcast(&size,sizeof(size),MPI_BYTE,proc,ExpA.comm);
      //Broadcast colptr
      if(iam==proc){
        MPI_Bcast(pMat.Local_.colptr.Data(),size*sizeof(Int),MPI_BYTE,proc,ExpA.comm);
      }
      else{
        recvColptr.Resize(size);
        MPI_Bcast(recvColptr.Data(),size*sizeof(Int),MPI_BYTE,proc,ExpA.comm);
      }

      //communicate rowind
      size = iam==proc?pMat.Local_.rowind.m():0;
      MPI_Bcast(&size,sizeof(size),MPI_BYTE,proc,ExpA.comm);
      //Broadcast rowind
      if(iam==proc){
        MPI_Bcast(pMat.Local_.rowind.Data(),size*sizeof(Int),MPI_BYTE,proc,ExpA.comm);
      }
      else{
        recvRowind.Resize(size);
        MPI_Bcast(recvRowind.Data(),size*sizeof(Int),MPI_BYTE,proc,ExpA.comm);
      }


      //communicate nzvalLocal
      //Broadcast size
      size = iam==proc?pMat.nzvalLocal.m():0;
      MPI_Bcast(&size,sizeof(size),MPI_BYTE,proc,ExpA.comm);
      //Broadcast nzvalLocal
      if(iam==proc){
        MPI_Bcast(pMat.nzvalLocal.Data(),size*sizeof(T),MPI_BYTE,proc,ExpA.comm);
      }
      else{
        recvNzval.Resize(size);
        MPI_Bcast(recvNzval.Data(),size*sizeof(T),MPI_BYTE,proc,ExpA.comm);
      }

      int colptrSize = iam==proc?pMat.Local_.colptr.m()-1:recvColptr.m()-1;
      int * colptr = iam==proc?pMat.Local_.colptr.Data():recvColptr.Data();
      int * rowind = iam==proc?pMat.Local_.rowind.Data():recvRowind.Data();
      T * nzval = iam==proc?pMat.nzvalLocal.Data():recvNzval.Data();
      //Parse the received data and copy it in the ExpA.nzvalLocal
      Int recvFirstCol = proc*numColFirst+1;
      Int recvLastCol = recvFirstCol+numColLocalVec[proc]-1;
      for(Int colIdx = 1; colIdx<=colptrSize;++colIdx){
        Int col = recvFirstCol + colIdx-1;
        Int colBeg = colptr[colIdx-1];
        Int colEnd = colptr[colIdx]-1;
        for(Int rowIdx = colBeg; rowIdx<=colEnd;++rowIdx){
          Int row = rowind[rowIdx-1];
          if(row>=localFirstCol && row<=localLastCol){
            //copy only upper triangular part
            if(col<row){
            Int localCol = row - localFirstCol +1;
            Int & localDestColIdx = localDestColHead[localCol-1];
            ExpA.nzvalLocal[localDestColIdx-1] = nzval[rowIdx-1];
            localDestColIdx++;
            }
          }
        }
      }
    }
    //copy upper triangular part
      for(Int colIdx = 1; colIdx<pMat.Local_.colptr.m();++colIdx){
        Int & localDestColIdx = localDestColHead[colIdx-1];
        Int colBeg = pMat.Local_.colptr[colIdx-1];
        Int colEnd = pMat.Local_.colptr[colIdx]-1;
        for(Int rowIdx = colBeg; rowIdx<=colEnd;++rowIdx){
          Int row = pMat.Local_.rowind[rowIdx-1];
          ExpA.nzvalLocal[localDestColIdx-1] = pMat.nzvalLocal[rowIdx-1];
          localDestColIdx++;
        }
      }
double timeEnd = get_time();


#ifdef _DEBUG_
logfileptr->OFS()<<"ExpA.nzvalLocal: "<<ExpA.nzvalLocal<<endl;
logfileptr->OFS()<<"pMat.nzvalLocal: "<<pMat.nzvalLocal<<endl;
#endif

      //now we can do the copy by doing sends of whole columns to the dest processor
if(iam==0){
  cout<<"Time for expanding matrix to asymmetric: "<<timeEnd-timeSta<<endl;
}
}

//exit(-1);
    //copy

      Icomm buffer;
      std::vector<T> denseA;
//      std::vector<char> recv_buffer;
      Icomm recv_buffer;
#if 1
    for(Int I=1;I<Xsuper_.m();I++){
      Int fc = Xsuper_(I-1);
      Int lc = Xsuper_(I)-1;
      Int iWidth = lc-fc+1;
      Int fi = xlindx_(I-1);
      Int li = xlindx_(I)-1;
      Int iHeight = li-fi+1;

      Int iDest = Mapping_->Map(I-1,I-1);

#ifdef _DEBUG_
      logfileptr->OFS()<<"Supernode "<<I<<" is owned by P"<<iDest<<std::endl;
#endif
      
      //parse the first column to create the supernode structure
      if(iam==iDest){
        LocalSupernodes_.push_back( new SuperNode<T>(I,fc,lc,iHeight));
        SuperNode<T> & snode = *LocalSupernodes_.back();


        for(Int idx = fi; idx<=li;idx++){
          Int iStartRow = lindx_(idx-1);
          Int iPrevRow = iStartRow;
          Int iContiguousRows = 1;
          for(Int idx2 = idx+1; idx2<=li;idx2++){
            Int iCurRow = lindx_(idx2-1);
            if(iStartRow == lindx_(fi-1)){
              if(iCurRow>iStartRow+iWidth-1){
                //enforce the first block to be a square diagonal block
                break;
              }
            }

            if(iCurRow==iPrevRow+1){
              idx++;
              ++iContiguousRows;
              iPrevRow=iCurRow;
            }
            else{
              break;
            }
          }

          Int iCurNZcnt = iContiguousRows * iWidth;
#ifdef _DEBUG_
//          logfileptr->OFS()<<"Creating a new NZBlock for rows "<<iStartRow<<" to "<<iStartRow + iContiguousRows-1<<" with "<<iCurNZcnt<<" nz."<<std::endl;
#endif
          snode.AddNZBlock(iContiguousRows , iWidth,iStartRow);

        }

        snode.Shrink();
      }

      //Distribute the data

      //look at the owner of the first column of the supernode
      Int numColFirst = std::max(1,iSize_ / np);
      //copy the data from A into this Block structure
      std::vector<T> recvNzval;
      std::vector<Int> recvColptr;
      std::vector<Int> recvRowind;
      for(Int col = fc;col<=lc;col++){
        //corresponding column in the unsorted matrix A
        Int orig_col = Order_.perm[col-1];

//        //determine last contiguous col in the original matrix A
//        Int lccol = col;
//        while(++lccol<=lc){
//          Int orig_lccol = Order_.perm[lccol-1];
//          if(orig_lccol != orig_col+1){
//            break;
//          }
//        }

        Int iOwnerCol = std::min((orig_col-1)/numColFirst,np-1);
#if 1




        //column is local
        if(iam == iDest){
          SuperNode<T> & snode = *LocalSupernodes_.back();
          if(iam != iOwnerCol){
            Int size = 0;
            //receive sizeof rowind
            MPI_Recv(&size,sizeof(size),MPI_BYTE,iOwnerCol,col,CommEnv_->MPI_GetComm(),MPI_STATUS_IGNORE);
            //init dummy colptr
            recvColptr.resize(2);
            recvColptr[0] = 1;
            recvColptr[1] = size+1;

            recvRowind.resize(size);
            MPI_Recv(&recvRowind[0],size*sizeof(Int),MPI_BYTE,iOwnerCol,col,CommEnv_->MPI_GetComm(),MPI_STATUS_IGNORE);
            recvNzval.resize(size);
            MPI_Recv(&recvNzval[0],size*sizeof(T),MPI_BYTE,iOwnerCol,col,CommEnv_->MPI_GetComm(),MPI_STATUS_IGNORE);
          }

          //int colptrSize = iOwnerCol==iDest?pMat.Local_.colptr.m()-1:1;
          int * colptr = iOwnerCol==iDest?ExpA.Local_.colptr.Data():&recvColptr[0];
          int * rowind = iOwnerCol==iDest?ExpA.Local_.rowind.Data():&recvRowind[0];
          T * nzvalA = iOwnerCol==iDest?ExpA.nzvalLocal.Data():&recvNzval[0];

          Int local_col = iOwnerCol==iDest?(orig_col-(numColFirst)*iOwnerCol):1;
          //Int offset = iOwnerCol==iDest?0:colptr[0];

          Int colbeg = colptr[local_col-1];
          Int colend = colptr[local_col]-1;
          for(Int rowidx = colbeg; rowidx<=colend; ++rowidx){
            Int orig_row = rowind[rowidx-1];
            Int row = Order_.invp[orig_row-1];

            if(row>=col){
//gdb_lock();
              Int blkidx = snode.FindBlockIdx(row);
              NZBlockDesc & blk_desc = snode.GetNZBlockDesc(blkidx);
              Int local_row = row - blk_desc.GIndex + 1;
              Int local_col = col - fc + 1;
              T * nzval = snode.GetNZval(blk_desc.Offset);
              nzval[(local_row-1)*iWidth+local_col-1] = nzvalA[rowidx-1];
            }
          }
        }
        else{
          if(iam == iOwnerCol){
            Int local_col = (orig_col-(numColFirst)*iOwnerCol);
            Int colbeg = ExpA.Local_.colptr[local_col-1];
            Int colend = ExpA.Local_.colptr[local_col]-1;
            Int size = colend-colbeg+1;
            //MPI_Send
            //send sizeof rowind
            MPI_Send(&size,sizeof(size),MPI_BYTE,iDest,col,CommEnv_->MPI_GetComm());
            MPI_Send(&ExpA.Local_.rowind[colbeg-1],size*sizeof(Int),MPI_BYTE,iDest,col,CommEnv_->MPI_GetComm());
            MPI_Send(&ExpA.nzvalLocal[colbeg-1],size*sizeof(T),MPI_BYTE,iDest,col,CommEnv_->MPI_GetComm());
          }
        }

#else
        buffer.clear();
        buffer.resize(iHeight*(sizeof(Int) + sizeof(T)));
        //If the column comes later in the permuted matrix,
        //NZvals are already in the LT part of the original matrix
        //Otherwise we need to look at the structure of row orig_col


        //If I own the column in the original matrix
        //Copy the appropriate values in the string stream



        //NZvals are in row orig_col of original matrix
        Int prev_col = col;
        do{
          Int orig_prev_col = Order_.perm[prev_col-1];
          Int iOwnerCurCol = std::min((orig_prev_col-1)/numColFirst,np-1);

          if(iam == iOwnerCurCol){
            Int local_col = (orig_prev_col-(numColFirst)*iOwnerCurCol);
            Int colbeg = Local_.colptr[local_col-1];
            Int colend = Local_.colptr[local_col]-1;
            const T * pdData = &pMat.nzvalLocal[0];

            for(Int rowidx = colbeg; rowidx<=colend; ++rowidx){
              Int orig_row = Local_.rowind[rowidx-1];
              Int row = Order_.invp[orig_row-1];
              //serialize the value
              if(orig_prev_col == orig_col){
                if(row>=col){
                  buffer<<row;
                  assert(row>0);
                  T val = pdData[rowidx-1];
                  buffer<<val;
                }
              }
              else if(row == col){
                buffer<<prev_col;
                assert(prev_col>0);
                T val = pdData[rowidx-1];
                buffer<<val;
                break;
              }
            }
          }

          prev_col = ETree_.PostParent(prev_col-1);
        }while(prev_col!=0);



#if 1
          mpi::Gatherv(buffer,recv_buffer,iDest,CommEnv_->MPI_GetComm());
          //now I know if I need to send something and who I'm receiving from
          if(iam == iDest){
            denseA.resize(iSize_);
                char * pRecvData;
                Int size = recv_buffer.size();
                pRecvData = recv_buffer.front();

                //Unpack
                Int pos = 0;
                while(pos<size){
                  Int row = *reinterpret_cast<Int *>(&pRecvData[pos]);
assert(row>0);
                  pos += sizeof(Int);
                  T val = *reinterpret_cast<T *>(&pRecvData[pos]);
                  pos += sizeof(T);

                  //put this into denseA
                  denseA.at(row-1) = val;
                }

            //now parse A, pick the values in denseA and put it in L

            SuperNode<T> & snode = *LocalSupernodes_.back();
            {
              Int colbeg = Global_.expColptr[orig_col-1];
              Int colend = Global_.expColptr[orig_col]-1;
              for(Int rowidx = colbeg; rowidx<=colend; ++rowidx){
                Int orig_row = Global_.expRowind[rowidx-1];
                Int row = Order_.invp[orig_row-1];

                if(row>=col){
                  Int blkidx = snode.FindBlockIdx(row);
                  NZBlockDesc & blk_desc = snode.GetNZBlockDesc(blkidx);
                  Int local_row = row - blk_desc.GIndex + 1;
                  Int local_col = col - fc + 1;
                  T * nzval = snode.GetNZval(blk_desc.Offset);
                  nzval[(local_row-1)*iWidth+local_col-1] = denseA.at(row-1);
                }
              }
            }
          }
#else
        std::vector<Int> isSender(np,0);
        if(iam==iDest){
          //If I'm the destination, I should count the number of senders

          Int prev_col = col;
          do{
            Int orig_prev_col = Order_.perm[prev_col-1];
            Int iOwnerCurCol = std::min((orig_prev_col-1)/numColFirst,np-1);

            Int colbeg = Global_.expColptr[orig_prev_col-1];
            Int colend = Global_.expColptr[orig_prev_col]-1;

            for(Int rowidx = colbeg; rowidx<=colend; ++rowidx){
              Int orig_row = Global_.expRowind[rowidx-1];
              Int row = Order_.invp[orig_row-1];
              //Look only at the lower triangular part
              if(orig_row>= orig_prev_col){
                //serialize the value
                if(orig_prev_col == orig_col){
                  if(row>=col){
                    isSender[iOwnerCurCol] = 1;//true;      
                  }
                }
                else if(row == col){
                  isSender[iOwnerCurCol] = 1;//true;      
                  break;
                }
              }
            }

            prev_col = ETree_.PostParent(prev_col-1);
          }while(prev_col!=0);

          //              logfileptr->OFS()<<isSender<<endl;
        }


          //now I know if I need to send something and who I'm receiving from
          if(iam == iDest){

            denseA.resize(iSize_);
            std::vector<Icomm> recvBuffers(np);

            std::vector<MPI_Request> requestSizes(np,MPI_REQUEST_NULL);
            std::vector<MPI_Request> requestContents(np,MPI_REQUEST_NULL);
            std::vector<Int> sizes(np,0);
            Int senderCnt = 0;
            for(Int psrc = 0; psrc<np; ++psrc){
              if(isSender[psrc]){
                if(psrc!=iam){
                  senderCnt++;
                  MPI_Irecv(&sizes[psrc],sizeof(Int),MPI_BYTE,psrc,2*col+0,CommEnv_->MPI_GetComm(),&requestSizes[psrc]);
                }
              }
            }

            //MPI_Waitall(np,&requestSizes[0],MPI_STATUSES_IGNORE);

            Int indx = MPI_UNDEFINED;
            MPI_Status recv_status;
            if(senderCnt>0){
              MPI_Waitany(np,&requestSizes[0],&indx,&recv_status);
              while(indx!=MPI_UNDEFINED){
                //Post the corresponding Irecv for content
                //MPI_Request & req = requestSizes[indx];
                Int psrc = recv_status.MPI_SOURCE;
                Int size = sizes[psrc];
                recvBuffers[psrc].resize(size);
                //Do the Irecv for content
                assert(isSender[psrc]);
                MPI_Irecv(recvBuffers[psrc].front(),size,MPI_BYTE,psrc,2*col+1,CommEnv_->MPI_GetComm(),&requestContents[psrc]);

                //Wait for another size
                indx = MPI_UNDEFINED;
                MPI_Waitany(np,&requestSizes[0],&indx,&recv_status);
              }
            }

            //If there is some local data, unpackit while we receive the remote content
            if(isSender[iam]){
              Int size = buffer.size();
              char * pRecvData = buffer.front();
              //Unpack
              Int pos = 0;
              while(pos<size){
                Int row = *reinterpret_cast<Int *>(&pRecvData[pos]);
                assert(row>0);
                pos += sizeof(Int);
                T val = *reinterpret_cast<T *>(&pRecvData[pos]);
                pos += sizeof(T);

                //put this into denseA
                denseA.at(row-1) = val;
              }
            }



            if(senderCnt>0){
              //While there is a new content
              indx = MPI_UNDEFINED;
              MPI_Waitany(np,&requestContents[0],&indx,&recv_status);
              while(indx!=MPI_UNDEFINED){
                MPI_Request & req = requestContents[indx];
                Int psrc = recv_status.MPI_SOURCE;

                Icomm & buffer = recvBuffers[psrc];
                Int size = buffer.size();
                char * pRecvData = buffer.front();
                //Unpack
                Int pos = 0;
                while(pos<size){
                  Int row = *reinterpret_cast<Int *>(&pRecvData[pos]);
                  assert(row>0);
                  pos += sizeof(Int);
                  T val = *reinterpret_cast<T *>(&pRecvData[pos]);
                  pos += sizeof(T);

                  //put this into denseA
                  denseA.at(row-1) = val;
                }


                //Wait for another content
                indx = MPI_UNDEFINED;
                MPI_Waitany(np,&requestContents[0],&indx,&recv_status);
              }
            }

            //now parse A, pick the values in denseA and put it in L

            SuperNode<T> & snode = *LocalSupernodes_.back();
            {
              Int colbeg = Global_.expColptr[orig_col-1];
              Int colend = Global_.expColptr[orig_col]-1;
              for(Int rowidx = colbeg; rowidx<=colend; ++rowidx){
                Int orig_row = Global_.expRowind[rowidx-1];
                Int row = Order_.invp[orig_row-1];

                if(row>=col){
                  Int blkidx = snode.FindBlockIdx(row);
                  NZBlockDesc & blk_desc = snode.GetNZBlockDesc(blkidx);
                  Int local_row = row - blk_desc.GIndex + 1;
                  Int local_col = col - fc + 1;
                  T * nzval = snode.GetNZval(blk_desc.Offset);
                  nzval[(local_row-1)*iWidth+local_col-1] = denseA.at(row-1);
                }
              }
            }


          } 
          else{
            if(buffer.size()>0){
              //Send the size first
              Int size = buffer.size();
              MPI_Send(&size,sizeof(Int),MPI_BYTE,iDest,2*col+0,CommEnv_->MPI_GetComm());
              MPI_Send(buffer.front(),size,MPI_BYTE,iDest,2*col+1,CommEnv_->MPI_GetComm());
            }
          }

//        MPI_Barrier(CommEnv_->MPI_GetComm()); 
#endif
#endif





          //now I know if I need to send something and who I'm receiving from
////          if(iam == iDest){
////            denseA.resize(iSize_);
////            recv_buffer.resize(iHeight*(sizeof(Int)+sizeof(T)));
////            for(Int psrc = 0; psrc<np; ++psrc){
////              if(isSender[psrc]){
////                Int size = 0;
////                char * pRecvData;
////                if(psrc!=iam){
////                  MPI_Status recv_status;
////                  MPI_Recv(&recv_buffer[0],recv_buffer.size(),MPI_BYTE,psrc,col,CommEnv_->MPI_GetComm(),&recv_status);
////                  MPI_Get_count(&recv_status,MPI_BYTE,&size);
//////logfileptr->OFS()<<"MPI_Recv of "<<size<<" bytes of col "<<col<<" from P"<<psrc<<endl;
////                  pRecvData = &recv_buffer[0];
////                }
////                else{
////                  size = buffer.size();
////                  pRecvData = buffer.front();
////                }
////
////                //Unpack
////                Int pos = 0;
////                while(pos<size){
////                  Int row = *reinterpret_cast<Int *>(&pRecvData[pos]);
////assert(row>0);
////                  pos += sizeof(Int);
////                  T val = *reinterpret_cast<T *>(&pRecvData[pos]);
////                  pos += sizeof(T);
////
////                  //put this into denseA
////                  denseA.at(row-1) = val;
////                }
////
////              }
////            }
////
////            //now parse A, pick the values in denseA and put it in L
////
////            SuperNode<T> & snode = *LocalSupernodes_.back();
////            {
////              Int colbeg = Global_.expColptr[orig_col-1];
////              Int colend = Global_.expColptr[orig_col]-1;
////              for(Int rowidx = colbeg; rowidx<=colend; ++rowidx){
////                Int orig_row = Global_.expRowind[rowidx-1];
////                Int row = Order_.invp[orig_row-1];
////
////                if(row>=col){
////                  Int blkidx = snode.FindBlockIdx(row);
////                  NZBlockDesc & blk_desc = snode.GetNZBlockDesc(blkidx);
////                  Int local_row = row - blk_desc.GIndex + 1;
////                  Int local_col = col - fc + 1;
////                  T * nzval = snode.GetNZval(blk_desc.Offset);
////                  nzval[(local_row-1)*iWidth+local_col-1] = denseA.at(row-1);
////                }
////              }
////            }
////          }
////          else{
////            if(buffer.size()>0){
//////gdb_lock();
////              //Send the size first
////              Int size = buffer.size();
//////              MPI_Send(&size,sizeof(Int),MPI_BYTE,iDest,iam,CommEnv_->MPI_GetComm());
//////logfileptr->OFS()<<"MPI_Send of "<<size<<" bytes of col "<<col<<" to P"<<iDest<<endl;
////              MPI_Send(buffer.front(),size,MPI_BYTE,iDest,col,CommEnv_->MPI_GetComm());
////            }
////          }

//        MPI_Barrier(CommEnv_->MPI_GetComm()); 
 
        }
#ifdef _DEBUG_
      logfileptr->OFS()<<"--------------------------------------------------"<<std::endl;
#endif
    }
#endif

#ifdef _DEBUG_
    for(Int I=1;I<Xsuper_.m();I++){
      Int src_first_col = Xsuper_(I-1);
      Int src_last_col = Xsuper_(I)-1;
      Int iOwner = Mapping_->Map(I-1,I-1);
      //If I own the column, factor it
      if( iOwner == iam ){
        Int iLocalI = (I-1) / np +1 ;
        SuperNode<T> & src_snode = *LocalSupernodes_[iLocalI -1];
          

        logfileptr->OFS()<<"+++++++++++++"<<I<<"++++++++"<<std::endl;

          logfileptr->OFS()<<"cols: ";
          for(Int i = 0; i< src_snode.Size(); ++i){
              logfileptr->OFS()<<" "<<Order_.perm[src_first_col+i-1];
          }
            logfileptr->OFS()<<std::endl;
            logfileptr->OFS()<<std::endl;
        for(int blkidx=0;blkidx<src_snode.NZBlockCnt();++blkidx){

          NZBlockDesc & desc = src_snode.GetNZBlockDesc(blkidx);
          T * val = src_snode.GetNZval(desc.Offset);
          Int nRows = src_snode.NRows(blkidx);

          Int row = desc.GIndex;
          for(Int i = 0; i< nRows; ++i){
              logfileptr->OFS()<<row+i<<" | "<<Order_.perm[row+i-1]<<":   ";
            for(Int j = 0; j< src_snode.Size(); ++j){
              logfileptr->OFS()<<val[i*src_snode.Size()+j]<<" ";
            }
            logfileptr->OFS()<<std::endl;
          }

        logfileptr->OFS()<<"_______________________________"<<std::endl;
        }
      }
    }
#endif




  }

  template <typename T> SupernodalMatrix<T>::SupernodalMatrix(){
    CommEnv_=NULL;
  }

  template <typename T> SupernodalMatrix<T>::SupernodalMatrix(const DistSparseMatrix<T> & pMat, Int maxSnode,Mapping * pMapping, Int maxIsend, Int maxIrecv, MPI_Comm & pComm ){
    CommEnv_ = NULL;
    Init(pMat, maxSnode,pMapping, maxIsend, maxIrecv, pComm );
  }

  template <typename T> SupernodalMatrix<T>::~SupernodalMatrix(){
    for(Int i=0;i<LocalSupernodes_.size();++i){
      delete LocalSupernodes_[i];
    }

    for(Int i=0;i<Contributions_.size();++i){
      delete Contributions_[i];
    }
    if(CommEnv_!=NULL){
      delete CommEnv_;
    }
  }

  template <typename T> void SupernodalMatrix<T>::GetUpdatingSupernodeCount(IntNumVec & sc,IntNumVec & mw){
    sc.Resize(Xsuper_.m());
    SetValue(sc,I_ZERO);
    IntNumVec marker(Xsuper_.m());
    SetValue(marker,I_ZERO);
    mw.Resize(Xsuper_.m());
    SetValue(mw,I_ZERO);

    for(Int s = 1; s<Xsuper_.m(); ++s){
      Int first_col = Xsuper_(s-1);
      Int last_col = Xsuper_(s)-1;

      Int fi = xlindx_(s-1);
      Int li = xlindx_(s)-1;

//#ifndef _DEBUG_
//  #define nodebugtmp
//  #define _DEBUG_
//#endif



//#ifdef nodebugtmp
//  #undef _DEBUG_
//#endif



      Int iOwner = Mapping_->Map(s-1,s-1);
#ifdef _DEBUG_UPDATES_
      logfileptr->OFS()<<"Supernode "<<s<<" on P"<<iOwner<<" updates: ";
#endif

      for(Int row_idx = fi; row_idx<=li;++row_idx){
        Int row = lindx_(row_idx-1);
        Int supno = SupMembership_(row-1);

        if(marker(supno-1)!=s && supno!=s){

#ifdef _DEBUG_UPDATES_
          logfileptr->OFS()<<supno<<" ";
#endif
          ++sc(supno-1);
          marker(supno-1) = s;

          mw(supno-1) = max(mw(supno-1),last_col - first_col+1);

        }
      }

#ifdef _DEBUG_UPDATES_
      logfileptr->OFS()<<std::endl;
#endif

//#ifdef nodebugtmp
//  #undef _DEBUG_
//#endif
    }
  }


  template <typename T> SparseMatrixStructure SupernodalMatrix<T>::GetLocalStructure() const {
    return Local_;
  }

  template <typename T> SparseMatrixStructure SupernodalMatrix<T>::GetGlobalStructure(){
    if(!globalAllocated){
      Local_.ToGlobal(Global_);
      globalAllocated= true;
    }
    return Global_;
  }


//Routines related to the packing of data
  template<typename T> void SupernodalMatrix<T>::AddOutgoingComm(AsyncComms & outgoingSend, Int src_snode_id, Int src_snode_size, Int src_first_row, NZBlockDesc & pivot_desc, Int nzblk_cnt, T * nzval_ptr, Int nz_cnt){

    outgoingSend.push_back(new Icomm(3*sizeof(Int) + nzblk_cnt*sizeof(NZBlockDesc) + nz_cnt*sizeof(T),MPI_REQUEST_NULL));

    *outgoingSend.back()<<src_snode_id;
    *outgoingSend.back()<<nzblk_cnt;
    NZBlockDesc * dest_blocks_ptr = reinterpret_cast<NZBlockDesc *>(outgoingSend.back()->back());
    Serialize(*outgoingSend.back(),&pivot_desc,nzblk_cnt);
    dest_blocks_ptr->GIndex = src_first_row;
    dest_blocks_ptr->Offset = pivot_desc.Offset + (src_first_row - pivot_desc.GIndex)*src_snode_size;
    *outgoingSend.back()<<nz_cnt;
    Serialize(*outgoingSend.back(),nzval_ptr,nz_cnt);
  }

  template<typename T> void SupernodalMatrix<T>::AddOutgoingComm(AsyncComms & outgoingSend, Icomm * send_buffer){
    outgoingSend.push_back(send_buffer);
  }





template <typename T> void SupernodalMatrix<T>::AdvanceOutgoing(AsyncComms & outgoingSend){
  //Check for completion of outgoing communication
  if(!outgoingSend.empty()){
    AsyncComms::iterator it = outgoingSend.begin();
    while(it != outgoingSend.end()){
      int flag = 0;
      int error_code = MPI_Test(&(*it)->Request,&flag,MPI_STATUS_IGNORE);
      if(flag){
        it = outgoingSend.erase(it);
      }
      else{
        it++;
      }
    }
  }
}








/// Routines used to find the next update in the current supernode
/// They must be called in a while loop style
  template <typename T> inline bool SupernodalMatrix<T>::FindNextUpdate(SuperNode<T> & src_snode, Int & tgt_snode_id, Int & f_ur, Int & f_ub, Int & n_ur, Int & n_ub){

    Int iam = CommEnv_->MPI_Rank();
    Int np  = CommEnv_->MPI_Size();
//    TIMER_START(FIND_UPDATE);

    if(tgt_snode_id==0){   
      Int iOwner = Mapping_->Map(src_snode.Id()-1,src_snode.Id()-1);
      f_ub = iOwner==iam?1:0;
      n_ub = f_ub;
      n_ur = 0;
    }
    else{
      f_ub = n_ub;
    }

    if(src_snode.NZBlockCnt()>f_ub){
      NZBlockDesc * cur_desc = &src_snode.GetNZBlockDesc(f_ub); 
      f_ur = max(n_ur,cur_desc->GIndex); 
      //find the snode updated by that row
      tgt_snode_id = SupMembership_[f_ur-1];
      Int tgt_lc = Xsuper_[tgt_snode_id]-1;

      //or use FindBlockIdx
//      if(f_ub<src_snode.NZBlockCnt()-1){
        Int src_tgt_lb = src_snode.FindBlockIdx(tgt_lc);
        //if tgt_lc not found in the current column we need to resume at the next block 
        if(src_tgt_lb==-1){
          for(n_ub;n_ub<src_snode.NZBlockCnt();++n_ub){
            cur_desc = &src_snode.GetNZBlockDesc(n_ub);
            if(cur_desc->GIndex > tgt_lc){
              break;
            }
          }
          if(n_ub<src_snode.NZBlockCnt()){
            n_ur = cur_desc->GIndex;
          }
          else{
            n_ur = -1;
          }
        }
        else{
            n_ub = src_tgt_lb;
            cur_desc = &src_snode.GetNZBlockDesc(n_ub);
            if(cur_desc->GIndex + src_snode.NRows(n_ub)-1>tgt_lc){
              n_ur = tgt_lc+1;
            }
            else{
              ++n_ub;
              n_ur = (n_ub<src_snode.NZBlockCnt())?src_snode.GetNZBlockDesc(n_ub).GIndex:-1;
            }
        }
      //src_snode updates tgt_snode_id. Then we need to look from row n_ur and block l_ub
      return true;
    }
    else{
      return false;
    }
  }















#include "SupernodalMatrix_impl_FO.hpp"

#include "SupernodalMatrix_impl_FB.hpp"











template <typename T> void SupernodalMatrix<T>::Factorize(){
  TIMER_START(FACTORIZATION);
  FanOut();
  TIMER_STOP(FACTORIZATION);
}


//Solve related routines

  template <typename T> void SupernodalMatrix<T>::forward_update(SuperNode<T> * src_contrib,SuperNode<T> * tgt_contrib){

    Int iam = CommEnv_->MPI_Rank();
    Int np  = CommEnv_->MPI_Size();

    Int iOwner = Mapping_->Map(src_contrib->Id()-1,src_contrib->Id()-1);
    Int src_ncols = src_contrib->Size();
    Int tgt_ncols = tgt_contrib->Size();

    Int startBlock = (iam==iOwner)?1:0;
    Int tgt_blkidx = -1;

    Int src_blkidx = startBlock;
    while(src_blkidx<src_contrib->NZBlockCnt()){
      NZBlockDesc & src_desc = src_contrib->GetNZBlockDesc(src_blkidx);
      Int src_nrows = src_contrib->NRows(src_blkidx);
  
      if(tgt_blkidx<1){
          tgt_blkidx = tgt_contrib->FindBlockIdx(src_desc.GIndex);
      }

      NZBlockDesc & tgt_desc = tgt_contrib->GetNZBlockDesc(tgt_blkidx);
      Int tgt_nrows = tgt_contrib->NRows(tgt_blkidx);

      Int src_local_fr = max(tgt_desc.GIndex - src_desc.GIndex,0);
      Int src_lr = src_desc.GIndex+src_nrows-1;

      Int tgt_local_fr = max(src_desc.GIndex - tgt_desc.GIndex,0);
      Int tgt_lr = tgt_desc.GIndex+tgt_nrows-1;
      Int tgt_local_lr = min(src_lr,tgt_lr) - tgt_desc.GIndex;

      T * src = &src_contrib->GetNZval(src_desc.Offset)[src_local_fr*src_ncols];
      T * tgt = &tgt_contrib->GetNZval(tgt_desc.Offset)[tgt_local_fr*tgt_ncols];


//TODO understand why this axpy makes the code crash
//      for(Int i=0; i<(tgt_local_lr - tgt_local_fr +1)*src_ncols;++i,++tgt,++src){
//        *tgt+=*src;
//      }
      for(Int i=0; i<(tgt_local_lr - tgt_local_fr +1)*src_ncols;++i){
        tgt[i]+=src[i];
      }
//      blas::Axpy((tgt_local_lr - tgt_local_fr +1)*src_ncols,
//          ONE<T>(),src,1,tgt,1);

      if(src_lr>tgt_lr){
        //the src block hasn't been completely used and is
        // updating some lines in the nz block just below the diagonal block
        //this case happens only locally for the diagonal block
//        assert(tgt_blkidx==0);
        //skip to the next tgt nz block
        ++tgt_blkidx;
      }
      else{
        //skip to the next src nz block
        ++src_blkidx;
        if(src_blkidx<src_contrib->NZBlockCnt()){
          NZBlockDesc & next_src_desc = src_contrib->GetNZBlockDesc(src_blkidx);
          tgt_blkidx = tgt_contrib->FindBlockIdx(next_src_desc.GIndex);
        }
      }
    }
  }

  template <typename T> void SupernodalMatrix<T>::back_update(SuperNode<T> * src_contrib,SuperNode<T> * tgt_contrib){
    Int nrhs = tgt_contrib->Size();
    for(Int blkidx = 1; blkidx<tgt_contrib->NZBlockCnt();++blkidx){
      NZBlockDesc & tgt_desc = tgt_contrib->GetNZBlockDesc(blkidx);
      Int tgt_nrows = tgt_contrib->NRows(blkidx);

      Int src_nzblk_idx = src_contrib->FindBlockIdx(tgt_desc.GIndex);
      Int tgt_lr = tgt_desc.GIndex+tgt_nrows-1;
      Int src_lr; 
      do
      {
        NZBlockDesc & src_desc = src_contrib->GetNZBlockDesc(src_nzblk_idx);
        Int src_nrows = src_contrib->NRows(src_nzblk_idx);
        src_lr = src_desc.GIndex+src_nrows-1;

        Int src_local_fr = max(tgt_desc.GIndex - src_desc.GIndex,0);

        Int tgt_local_fr = max(src_desc.GIndex - tgt_desc.GIndex,0);
        Int tgt_local_lr = min(src_lr,tgt_lr) - tgt_desc.GIndex;

        T * src = &src_contrib->GetNZval(src_desc.Offset)[src_local_fr*nrhs];
        T * tgt = &tgt_contrib->GetNZval(tgt_desc.Offset)[tgt_local_fr*nrhs];

        std::copy(src,src+(tgt_local_lr - tgt_local_fr +1)*nrhs,tgt);
        //              lapack::Lacpy('N',nrhs,(tgt_local_lr - tgt_local_fr +1),
        //                  src,nrhs,  
        //                  tgt,nrhs);

        //do other block
        if(tgt_lr>src_lr){
//          assert(src_nzblk_idx==0);
          src_nzblk_idx++;
        }

      } while(tgt_lr>src_lr);
    }
  }

#ifdef _CHECK_RESULT_SEQ_
  template <typename T> void SupernodalMatrix<T>::Solve(NumMat<T> * RHS,  NumMat<T> & forwardSol, NumMat<T> * Xptr)
#else
  template <typename T> void SupernodalMatrix<T>::Solve(NumMat<T> * RHS,  NumMat<T> * Xptr)
#endif
{
    TIMER_START(SPARSE_SOLVE);

    NumMat<T> & B = *RHS;
    Int nrhs = B.n();

    Int iam = CommEnv_->MPI_Rank();
    Int np  = CommEnv_->MPI_Size();

    IntNumVec children(Xsuper_.m());
    SetValue(children,0);

    for(Int I=1;I<Xsuper_.m()-1;I++){
      Int fc = Xsuper_(I-1);
      Int lc = Xsuper_(I)-1;

      Int parent = ETree_.PostParent(lc-1);
      if(parent!=0){
        ++children(SupMembership_(parent-1)-1);
      }
    }

#ifdef _DEBUG_
    logfileptr->OFS()<<"Children vector is"<<children<<std::endl;
#endif


    IntNumVec UpdatesToDo = children;

    Contributions_.resize(LocalSupernodes_.size());
    std::vector<std::stack<Int> > LocalUpdates(LocalSupernodes_.size());

    AsyncComms outgoingSend;

    //This corresponds to the k loop in dtrsm
    for(Int I=1;I<Xsuper_.m();I++){
      Int iOwner = Mapping_->Map(I-1,I-1);
      //If I own the column, factor it
      if( iOwner == iam ){
        //Create the blocks of my contrib with the same nz structure as L
        //and the same width as the final solution
        Int iLocalI = (I-1) / np +1 ;
        SuperNode<T> * cur_snode = LocalSupernodes_[iLocalI-1];
        SuperNode<T> * contrib = new SuperNode<T>(I,1,nrhs, cur_snode->NRowsBelowBlock(0) );
        Contributions_[iLocalI-1] = contrib;


        for(Int blkidx = 0; blkidx<cur_snode->NZBlockCnt();++blkidx){
          NZBlockDesc & cur_desc = cur_snode->GetNZBlockDesc(blkidx);
          contrib->AddNZBlock(cur_snode->NRows(blkidx),nrhs,cur_desc.GIndex);
        }
        Int nRows = contrib->NRowsBelowBlock(0);
        std::fill(contrib->GetNZval(0),contrib->GetNZval(0)+nRows*nrhs,ZERO<T>());

        contrib->Shrink();
      }
    }





    CommList ContribsToSend; 
    DownCommList ContribsToSendDown; 


    std::vector<char> src_blocks;


    //forward-substitution phase
    //Sending contrib up the tree
    //Start from the leaves of the tree
    TIMER_START(SPARSE_FWD_SUBST);

    Int I =1;
    Int iLocalI =1;
    while(iLocalI<=LocalSupernodes_.size() || !ContribsToSend.empty() || !outgoingSend.empty()){
      //Check for completion of outgoing communication
      AdvanceOutgoing(outgoingSend);

      //process some of the delayed send

      if(iLocalI>0 && iLocalI<=LocalSupernodes_.size()){

      //If I own the column, factor it
        SuperNode<T> * cur_snode = LocalSupernodes_[iLocalI-1];
        I = cur_snode->Id();
        Int parent = ETree_.PostParent(cur_snode->LastCol()-1);
        SuperNode<T> * contrib = Contributions_[iLocalI-1];


        //Do all my updates (Local and remote)
        //Local updates
        while(!LocalUpdates[iLocalI-1].empty()){
          Int contrib_snode_id = LocalUpdates[iLocalI-1].top();
          LocalUpdates[iLocalI-1].pop();

          SuperNode<T> * dist_contrib = Contributions_[(contrib_snode_id-1) / np];

#ifdef _DEBUG_
          logfileptr->OFS()<<"LOCAL Supernode "<<I<<" is updated by contrib of Supernode "<<contrib_snode_id<<std::endl;
#endif

          forward_update(dist_contrib,contrib);
          //delete contributions[(contrib_snode_id-1) / np];
          --UpdatesToDo(I-1);
        }

        //do remote updates
        size_t max_bytes;
        Int nz_cnt;
        while(UpdatesToDo(I-1)>0){
          //receive children contrib
#ifdef _DEBUG_
          logfileptr->OFS()<<UpdatesToDo(I-1)<<" contribs left"<<endl;
#endif


          MPI_Status recv_status;
          int bytes_received = 0;
#if 1
          MPI_Probe(MPI_ANY_SOURCE,I,CommEnv_->MPI_GetComm(),&recv_status);
          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
          src_blocks.resize(bytes_received);
          MPI_Recv(&src_blocks[0],bytes_received,MPI_BYTE,recv_status.MPI_SOURCE,I,CommEnv_->MPI_GetComm(),&recv_status);
#else
          TIMER_START(RECV_MALLOC);
          max_bytes = 5*sizeof(Int); 
          Int nrows = cur_snode->NRowsBelowBlock(0);
          Int ncols = nrhs;
          nz_cnt = nrows * ncols;

          Int nblocks = nrows;//std::max((Int)ceil(nrows/2)+1,cur_snode->NZBlockCnt());
          max_bytes += (nblocks)*sizeof(NZBlockDesc);
          max_bytes += nz_cnt*sizeof(T); 

          src_blocks.resize(max_bytes);
          TIMER_STOP(RECV_MALLOC);

          TIMER_START(RECV_MPI);
#ifdef PROBE_FIRST
          MPI_Probe(MPI_ANY_SOURCE,I,CommEnv_->MPI_GetComm(),&recv_status);
          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);

          logfileptr->OFS()<<"Preparing to receive "<<bytes_received<<" bytes"<<endl;

          bool doabort = false;
          int prev_size = 0;
          if(src_blocks.size()<bytes_received){

            cout<<"We have a problem !!!! on P"<<iam<<"\n";
            gdb_lock();
            prev_size = src_blocks.size();
            doabort = true;
            //receive anyway
            src_blocks.resize(bytes_received);
          }
#endif

#ifdef PROBE_FIRST
          logfileptr->OFS()<<"Receiving from P"<<recv_status.MPI_SOURCE<<endl;
          MPI_Recv(&src_blocks[0],src_blocks.size(),MPI_BYTE,recv_status.MPI_SOURCE,I,CommEnv_->MPI_GetComm(),&recv_status);
#else
          MPI_Recv(&src_blocks[0],src_blocks.size(),MPI_BYTE,MPI_ANY_SOURCE,I,CommEnv_->MPI_GetComm(),&recv_status);
#endif
#endif

          SuperNode<T> dist_contrib;
          Deserialize(&src_blocks[0],dist_contrib);
//#ifdef PROBE_FIRST
//          if(doabort){
//
//            abort();
//          }
//#endif

#ifdef _DEBUG_
          logfileptr->OFS()<<"RECV contrib of Supernode "<<dist_contrib.Id()<<std::endl;
#endif

          forward_update(&dist_contrib,contrib);

          --UpdatesToDo(I-1);

        }

        assert(UpdatesToDo(I-1)==0);

        if(UpdatesToDo(I-1)==0){
//          logfileptr->OFS()<<"FORW Processing contrib "<<I<<std::endl;

          //This corresponds to the i loop in dtrsm
          for(Int blkidx = 0; blkidx<cur_snode->NZBlockCnt();++blkidx){

            NZBlockDesc & cur_desc = contrib->GetNZBlockDesc(blkidx);
            NZBlockDesc & chol_desc = cur_snode->GetNZBlockDesc(blkidx);
            NZBlockDesc & diag_desc = contrib->GetNZBlockDesc(0);

            Int cur_nrows = contrib->NRows(blkidx);
            Int chol_nrows = cur_snode->NRows(blkidx);
            Int diag_nrows = contrib->NRows(0);

            T * cur_nzval = contrib->GetNZval(cur_desc.Offset);
            T * chol_nzval = cur_snode->GetNZval(chol_desc.Offset);
            T * diag_nzval = contrib->GetNZval(diag_desc.Offset);

            //compute my contribution
            //Handle the diagonal block
            if(blkidx==0){
              //TODO That's where we can use the selective inversion
              //if we are processing the "pivot" block
              for(Int kk = 0; kk<cur_snode->Size(); ++kk){
                for(Int j = 0; j<nrhs;++j){
                  //logfileptr->OFS()<<"                      B = "<<B(Order_.perm[diag_desc.GIndex+kk-1]-1,j)<<std::endl;
                  Int srcRow = Order_.perm[diag_desc.GIndex+kk-1];
                  diag_nzval[kk*nrhs+j] = (B(srcRow-1,j) + diag_nzval[kk*nrhs+j]) / chol_nzval[kk*cur_snode->Size()+kk];
                  //diag_nzval[kk*nrhs+j] = (B(diag_desc.GIndex-1+kk,j) + diag_nzval[kk*nrhs+j]) / chol_nzval[kk*cur_snode->Size()+kk];
                  for(Int i = kk+1; i<cur_nrows;++i){
                    diag_nzval[i*nrhs+j] += -diag_nzval[kk*nrhs+j]*chol_nzval[i*cur_snode->Size()+kk];
                  }
                }
              }
            }
            else{
              for(Int kk = 0; kk<cur_snode->Size(); ++kk){
                for(Int j = 0; j<nrhs;++j){
                  for(Int i = 0; i<cur_nrows;++i){
                    cur_nzval[i*nrhs+j] += -diag_nzval[kk*nrhs+j]*chol_nzval[i*cur_snode->Size()+kk];
                  }
                }
              }
            }
          }


          //send to my parent
          if(parent!=0){
            Int parent_snode_id = SupMembership_[parent-1];

            Int iTarget = Mapping_->Map(parent_snode_id-1,parent_snode_id-1);

            if(iTarget!=iam){
#ifdef _DEBUG_
              logfileptr->OFS()<<"Remote Supernode "<<parent_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
#endif
              Int tgt_first_col = Xsuper_(parent_snode_id-1);
              Int tgt_last_col = Xsuper_(parent_snode_id)-1;
              Int src_nzblk_idx = 1;
              NZBlockDesc & pivot_desc = contrib->GetNZBlockDesc(src_nzblk_idx);

              Int src_first_row = pivot_desc.GIndex;

              bool isSkipped= false;

              Int next_local_contrib = (iLocalI < Contributions_.size())?Contributions_[iLocalI]->Id():Xsuper_.m();
                  if(next_local_contrib< parent_snode_id){
                  //need to push the prev src_last_row
                  ContribsToSend.push(DelayedComm(contrib->Id(),parent_snode_id,1,src_first_row));
#ifdef _DEBUG_DELAY_
                                      cout<<"P"<<iam<<" has delayed update from Contrib "<<I<<" to "<<parent_snode_id<<" from row "<<src_first_row<<endl;
#endif
                    isSkipped= true;
                }
              
              if(!isSkipped){
                    //Create a new Icomm buffer, serialize the contribution
                    // in it and add it to the outgoing comm list
                    Icomm * send_buffer = new Icomm();
                    Serialize(*send_buffer,*contrib,src_nzblk_idx,src_first_row);
                    AddOutgoingComm(outgoingSend,send_buffer);

                    if( outgoingSend.size() > maxIsend_){
                      TIMER_START(SEND_MPI);
                      MPI_Send(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,parent_snode_id,CommEnv_->MPI_GetComm());
                      TIMER_STOP(SEND_MPI);
                      outgoingSend.pop_back();
                    }
                    else{
                      TIMER_START(SEND_MPI);
                      MPI_Isend(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,parent_snode_id,CommEnv_->MPI_GetComm(),&outgoingSend.back()->Request);
                      TIMER_STOP(SEND_MPI);
                    }
#ifdef _DEBUG_            
                logfileptr->OFS()<<"     Send contribution "<<I<<" to Supernode "<<parent_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
//                logfileptr->OFS()<<"Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
//                logfileptr->OFS()<<"Sending "<<nzblk_cnt*sizeof(NZBlockDesc)+nz_cnt*sizeof(T) + 3*sizeof(Int)<<" bytes"<< std::endl;
#endif
              }
            }
            else{
#ifdef _DEBUG_
              logfileptr->OFS()<<"Local Supernode "<<parent_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
#endif
              Int iLocalJ = (parent_snode_id-1) / np +1 ;
              LocalUpdates[iLocalJ-1].push((Int)I);
            }
          }

      }



          }


      SendDelayedMessagesUp(iLocalI,ContribsToSend,outgoingSend,Contributions_);

#if 0
#ifdef _CHECK_RESULT_SEQ_
      //if(iLocalI>0 && iLocalI<=LocalSupernodes_.size())
      {
      MPI_Barrier(CommEnv_->MPI_GetComm());
      NumMat<T> tmp = B;
      GetSolution(tmp);


      Int nrows = 0;
      for(Int ii=1; ii<=I;++ii){ nrows+= Xsuper_[ii] - Xsuper_[ii-1];}

      NumMat<T> tmp3(tmp.m(),tmp.n());
      NumMat<T> tmpFwd(tmp.m(),tmp.n());
      for(Int i = 0; i<tmp.m();++i){
        for(Int j = 0; j<tmp.n();++j){
//          tmp3(Order_.perm[i]-1,j) = tmp(i,j);
//          tmpFwd(Order_.perm[i]-1,j) = forwardSol(i,j);

          tmp3(i,j) = tmp(Order_.perm[i]-1,j);
          tmpFwd(i,j) = forwardSol(Order_.perm[i]-1,j);
        }
      }

      NumMat<T> tmp2 = tmp3;
      
      blas::Axpy(tmp.m()*tmp.n(),-1.0,&tmpFwd(0,0),1,&tmp2(0,0),1);
      double norm = lapack::Lange('F',nrows,tmp.n(),&tmp2(0,0),tmp.m());
      logfileptr->OFS()<<"Norm after SuperNode "<<I<<" is "<<norm<<std::endl; 

        if(abs(norm)>=1e-1){
          for(Int i = 0;i<nrows/*tmp.m()*/;++i){
            logfileptr->OFS()<<tmpFwd(i,0)<<"       "<<tmp3(i,0)<<std::endl;
          }
        }
      }
#endif
#endif
          ++iLocalI;
    }

    while(!outgoingSend.empty()){
      AdvanceOutgoing(outgoingSend);
    } 
    MPI_Barrier(CommEnv_->MPI_GetComm());

    TIMER_STOP(SPARSE_FWD_SUBST);


    //Back-substitution phase
    TIMER_START(SPARSE_BACK_SUBST);

    //start from the root of the tree
    iLocalI = LocalSupernodes_.size() ;
    while(iLocalI>0|| !ContribsToSend.empty() || !outgoingSend.empty()){

      //Check for completion of outgoing communication
      AdvanceOutgoing(outgoingSend);

      if(iLocalI>0 && iLocalI<=LocalSupernodes_.size()){
        SuperNode<T> * cur_snode = LocalSupernodes_[iLocalI-1];
        I = cur_snode->Id();

        SuperNode<T> * contrib = Contributions_[iLocalI-1];

        Int parent = ETree_.PostParent(cur_snode->LastCol()-1);

        std::vector<Int> isBlockUpdated(contrib->NZBlockCnt(),0);

        if(parent!=0){
          Int parent_snode_id = SupMembership_[parent-1];
          Int iTarget = Mapping_->Map(parent_snode_id-1,parent_snode_id-1);
          //Do all my updates (Local and remote)
          //Local updates
          SuperNode<T> * dist_contrib;
          if(!LocalUpdates[iLocalI-1].empty()){
            Int contrib_snode_id = LocalUpdates[iLocalI-1].top();
            LocalUpdates[iLocalI-1].pop();

            dist_contrib = Contributions_[(contrib_snode_id-1) / np];
            back_update(dist_contrib,contrib);
          }
          else{
            //Receive parent contrib
            MPI_Status recv_status;
            Int bytes_received = 0;
            MPI_Probe(iTarget,I,CommEnv_->MPI_GetComm(),&recv_status);
            MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
            src_blocks.resize(bytes_received);

            MPI_Recv(&src_blocks[0],bytes_received,MPI_BYTE,iTarget,I,CommEnv_->MPI_GetComm(),&recv_status);


            dist_contrib = new SuperNode<T>();
            Deserialize(&src_blocks[0],*dist_contrib); 


            #ifdef _DEBUG_
            logfileptr->OFS()<<"RECV contrib of Supernode "<<dist_contrib->Id()<<std::endl;
            #endif

            back_update(dist_contrib,contrib);
            delete dist_contrib;
          }
        }

        //now compute MY contribution
#ifdef _DEBUG_
        logfileptr->OFS()<<"BACK Processing contrib "<<I<<std::endl;
#endif

        NZBlockDesc & diag_desc = cur_snode->GetNZBlockDesc(0);
        NZBlockDesc & tgt_desc = contrib->GetNZBlockDesc(0);

        T* diag_nzval = cur_snode->GetNZval(diag_desc.Offset);
        T* tgt_nzval = contrib->GetNZval(tgt_desc.Offset);

        for(Int j = 0; j<nrhs;++j){
          for(Int ii = cur_snode->Size()-1; ii>=0; --ii){
            T temp = tgt_nzval[ii*nrhs+j];

            //This corresponds to the k loop in dtrsm
            for(Int blkidx = 0; blkidx<cur_snode->NZBlockCnt();++blkidx){
              NZBlockDesc & chol_desc = cur_snode->GetNZBlockDesc(blkidx);
              Int chol_nrows = cur_snode->NRows(blkidx);

              Int src_blkidx = contrib->FindBlockIdx(chol_desc.GIndex);
              NZBlockDesc & cur_desc = contrib->GetNZBlockDesc(src_blkidx);
              Int cur_nrows = contrib->NRows(src_blkidx);

              T* chol_nzval = cur_snode->GetNZval(chol_desc.Offset);
              T* cur_nzval = contrib->GetNZval(cur_desc.Offset);

              for(Int kk = 0; kk< chol_nrows; ++kk){
                if(chol_desc.GIndex+kk>cur_snode->FirstCol()+ii){
                  Int src_row = chol_desc.GIndex - cur_desc.GIndex +kk;
                  if(src_row< cur_nrows){
                    temp += -chol_nzval[kk*cur_snode->Size()+ii]*cur_nzval[src_row*nrhs+j];
                  }
                }
              }
            }

            temp = temp / diag_nzval[ii*cur_snode->Size()+ii];
            tgt_nzval[ii*nrhs+j] = temp;
          }
        }


        //send to my children
        Int colIdx = cur_snode->FirstCol()-1;
        if(colIdx>0){
          Int children_found = 0;
          while(children_found<children(I-1)){
            Int child_snode_id = SupMembership_[colIdx-1];

            Int parent = ETree_.PostParent(colIdx-1);
            if(parent!=0){
              if(SupMembership_[parent-1]==cur_snode->Id()){
                Int iTarget = Mapping_->Map(child_snode_id-1,child_snode_id-1);

                if(iTarget!=iam){

                  bool isSkipped= false;

                  Int next_local_contrib = (iLocalI >1)?Contributions_[iLocalI-2]->Id():0;
                  if(next_local_contrib > child_snode_id){
                    //need to push the prev src_last_row
                    //Send
                    Int src_nzblk_idx = 0;
                    NZBlockDesc & pivot_desc = contrib->GetNZBlockDesc(src_nzblk_idx);
                    Int src_first_row = pivot_desc.GIndex;
                    ContribsToSendDown.push(DelayedComm(contrib->Id(),child_snode_id,0,src_first_row));
#ifdef _DEBUG_DELAY_
                    cout<<"P"<<iam<<" has delayed update from Contrib "<<I<<" to "<<child_snode_id<<" from row "<<src_first_row<<endl;
#endif
                    isSkipped= true;
                  }


                  if(!isSkipped){
#ifdef _DEBUG_
                    logfileptr->OFS()<<"Remote Supernode "<<child_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
#endif

                    Int src_nzblk_idx = 0;
                    NZBlockDesc & pivot_desc = contrib->GetNZBlockDesc(src_nzblk_idx);
                    Int src_first_row = pivot_desc.GIndex;
                    //Create a new Icomm buffer, serialize the contribution
                    // in it and add it to the outgoing comm list
                    Icomm * send_buffer = new Icomm();
                    Serialize(*send_buffer,*contrib,src_nzblk_idx,src_first_row);
                    AddOutgoingComm(outgoingSend,send_buffer);


                    if( outgoingSend.size() > maxIsend_){
                      TIMER_START(SEND_MPI);
                      MPI_Send(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,child_snode_id,CommEnv_->MPI_GetComm());
                      TIMER_STOP(SEND_MPI);
                      outgoingSend.pop_back();
                    }
                    else{
                      TIMER_START(SEND_MPI);
                      MPI_Isend(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,child_snode_id,CommEnv_->MPI_GetComm(),&outgoingSend.back()->Request);
                      TIMER_STOP(SEND_MPI);
                    }

#ifdef _DEBUG_            
                    logfileptr->OFS()<<"     Send contribution "<<I<<" to Supernode "<<child_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
//                    logfileptr->OFS()<<"Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
//                    logfileptr->OFS()<<"Sending "<<nzblk_cnt*sizeof(NZBlockDesc)<<" and "<<nz_cnt*sizeof(T)<<" bytes during BS"<<std::endl;
#endif
                  }
                }
                else{

#ifdef _DEBUG_
                  logfileptr->OFS()<<"Local Supernode "<<child_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
#endif
                  Int iLocalJ = (child_snode_id-1) / np +1 ;
                  LocalUpdates[iLocalJ-1].push((Int)I);
                }
                children_found++;
              }
            }
            //last column of the prev supernode
            colIdx = Xsuper_[child_snode_id-1]-1;
            if(colIdx==0){
              break;
            }
          }


        }
      }

      SendDelayedMessagesDown(iLocalI,ContribsToSendDown,outgoingSend,Contributions_);
      --iLocalI;
    }
    TIMER_STOP(SPARSE_BACK_SUBST);

    MPI_Barrier(CommEnv_->MPI_GetComm());
    TIMER_STOP(SPARSE_SOLVE);

}

template <typename T> void SupernodalMatrix<T>::GetFullFactors( NumMat<T> & fullMatrix){
    Int iam = CommEnv_->MPI_Rank();
    Int np  = CommEnv_->MPI_Size();
    if(iam==0){
      fullMatrix.Resize(this->Size(),this->Size());
      SetValue(fullMatrix,ZERO<T>());
    }



    //output L
    for(Int I=1;I<Xsuper_.m();I++){
      Int src_first_col = Xsuper_(I-1);
      Int src_last_col = Xsuper_(I)-1;
      Int iOwner = Mapping_->Map(I-1,I-1);
      //If I own the column, factor it

if( iOwner == iam){
        Int iLocalI = (I-1) / np +1 ;
        SuperNode<T> & src_snode = *LocalSupernodes_[iLocalI -1];
#ifdef _DEBUG_
        logfileptr->OFS()<<"Supernode "<<I<<"("<<src_snode.Id()<<") is on P"<<iOwner<<" local index is "<<iLocalI<<std::endl; 
        logfileptr->OFS()<<src_snode<<std::endl;
#endif
}
      if( iOwner == iam  && iam != 0){



        Int iLocalI = (I-1) / np +1 ;
        SuperNode<T> & src_snode = *LocalSupernodes_[iLocalI -1];


        NZBlockDesc * nzblk_desc = &src_snode.GetNZBlockDesc(0);
        Int size_blocks = src_snode.NZBlockCnt();
        MPI_Send((void*)&size_blocks,sizeof(Int),MPI_BYTE,0,I,CommEnv_->MPI_GetComm());
        MPI_Send((void*)nzblk_desc,size_blocks*sizeof(NZBlockDesc),MPI_BYTE,0,I,CommEnv_->MPI_GetComm());

        T * nzblk_nzval = src_snode.GetNZval(0);
        Int size_nzval = src_snode.NRowsBelowBlock(0)*src_snode.Size();
        MPI_Send((void*)&size_nzval,sizeof(Int),MPI_BYTE,0,I,CommEnv_->MPI_GetComm());
        MPI_Send((void*)nzblk_nzval,size_nzval*sizeof(T),MPI_BYTE,0,I,CommEnv_->MPI_GetComm());



      }

      if(iam==0){
        if(iOwner != iam){ 
          Int snode_size = Xsuper_[I] - Xsuper_[I-1];

          Int size_blocks, size_nzval;
          std::vector<NZBlockDesc> blocks;
          std::vector<T> nzval;


          MPI_Recv(&size_blocks,sizeof(Int),MPI_BYTE,iOwner,I,CommEnv_->MPI_GetComm(),MPI_STATUS_IGNORE);
          blocks.resize(size_blocks); 
          MPI_Recv(&blocks[0],size_blocks*sizeof(NZBlockDesc),MPI_BYTE,iOwner,I,CommEnv_->MPI_GetComm(),MPI_STATUS_IGNORE);


          MPI_Recv(&size_nzval,sizeof(Int),MPI_BYTE,iOwner,I,CommEnv_->MPI_GetComm(),MPI_STATUS_IGNORE);
          nzval.resize(size_nzval); 
          MPI_Recv(&nzval[0],size_nzval*sizeof(T),MPI_BYTE,iOwner,I,CommEnv_->MPI_GetComm(),MPI_STATUS_IGNORE);

          SuperNode<T> src_snode(I,Xsuper_[I-1],Xsuper_[I]-1,&blocks[0],size_blocks,&nzval[0],size_nzval);

          for(int blkidx=0;blkidx<src_snode.NZBlockCnt();++blkidx){
            NZBlockDesc & nzblk_desc = src_snode.GetNZBlockDesc(blkidx);
            Int nRows = src_snode.NRows(blkidx);
            T * nzblk_nzval = src_snode.GetNZval(nzblk_desc.Offset);

            T * dest = fullMatrix.Data();
            for(Int i = 0; i<nRows;++i){
              for(Int j = 0; j<src_snode.Size();++j){
                if(nzblk_desc.GIndex -1 + i >= src_snode.FirstCol()-1+j){
                  dest[nzblk_desc.GIndex -1 + i + (src_snode.FirstCol()-1+j)*fullMatrix.m()] = nzblk_nzval[i * src_snode.Size() + j];
                }
              }
            }
          } 
        }
        else{

          Int iLocalI = (I-1) / np +1 ;
          SuperNode<T> & src_snode = *LocalSupernodes_[iLocalI -1];
          for(int blkidx=0;blkidx<src_snode.NZBlockCnt();++blkidx){
            NZBlockDesc & nzblk_desc = src_snode.GetNZBlockDesc(blkidx);
            Int nRows = src_snode.NRows(blkidx);
            T * nzblk_nzval = src_snode.GetNZval(nzblk_desc.Offset);

            T * dest = fullMatrix.Data();
            for(Int i = 0; i<nRows;++i){
              for(Int j = 0; j<src_snode.Size();++j){
                if(nzblk_desc.GIndex -1 + i >= src_snode.FirstCol()-1+j){
                  dest[nzblk_desc.GIndex -1 + i + (src_snode.FirstCol()-1+j)*fullMatrix.m()] = nzblk_nzval[i * src_snode.Size() + j];
                }
              }
            }
          }
        }
      }
      MPI_Barrier(CommEnv_->MPI_GetComm());
    }

    
}

template<typename T> void SupernodalMatrix<T>::GetSolution(NumMat<T> & B){
    Int iam = CommEnv_->MPI_Rank();
    Int np  = CommEnv_->MPI_Size();
    Int nrhs = B.n();
    //Gather B from everybody and put it in the original matrix order
    std::vector<T> tmp_nzval;
    for(Int I=1; I<Xsuper_.m();++I){
      Int iOwner = Mapping_->Map(I-1,I-1);
      T * data;
      Int snode_size = Xsuper_[I] - Xsuper_[I-1];
      Int nzcnt = snode_size * nrhs;
      tmp_nzval.resize(nzcnt);

      if( iOwner == iam ){
        Int iLocalI = (I-1) / np +1 ;
        SuperNode<T> * contrib = Contributions_[iLocalI-1];
        data = contrib->GetNZval(0);
      }
      else{
        data = &tmp_nzval[0];
      }

      MPI_Bcast(data,nzcnt*sizeof(T),MPI_BYTE,iOwner,CommEnv_->MPI_GetComm());

      for(Int i = 0; i<snode_size; ++i){ 
        for(Int j = 0; j<nrhs; ++j){
          Int destRow = Xsuper_[I-1] + i;
          destRow = Order_.perm[destRow - 1];
          B(destRow-1, j) = data[i*nrhs + j];
        }
      }
    }
      MPI_Barrier(CommEnv_->MPI_GetComm());
}


template <typename T> void SupernodalMatrix<T>::SendDelayedMessagesDown(Int iLocalI, DownCommList & MsgToSend, AsyncComms & OutgoingSend, std::vector<SuperNode<T> *> & snodeColl){
  //typedef volatile LIBCHOLESKY::Int Int;

  if(snodeColl.empty() || MsgToSend.empty()) { return;}

  //Index of the first local supernode
  Int first_local_id = snodeColl.front()->Id();
  //Index of the last PROCESSED supernode
  Int prev_snode_id = iLocalI>=1?snodeColl[iLocalI-1]->Id():first_local_id;
  //Index of the next local supernode
  Int next_snode_id = iLocalI<=1?0:snodeColl[iLocalI-2]->Id();

  bool is_last = prev_snode_id<=1;

#ifdef _DEBUG_DELAY_
  {
    DownCommList tmp = MsgToSend;
    logfileptr->OFS()<<"Queue : ";
    while( tmp.size()>0){
      //Pull the highest priority message
      const DelayedComm & comm = tmp.top();
      Int src_snode_id = comm.src_snode_id;
      Int tgt_id = comm.tgt_snode_id;
      Int src_nzblk_idx = comm.src_nzblk_idx;
      Int src_first_row = comm.src_first_row;
      logfileptr->OFS()<<" { "<<src_snode_id<<" -> "<<tgt_id<<" }";
      tmp.pop();
    }
    logfileptr->OFS()<<endl;


  }
#endif

  while( MsgToSend.size()>0){
    //Pull the highest priority message
    const DelayedComm & comm = MsgToSend.top();
    Int src_snode_id = comm.src_snode_id;
    Int tgt_snode_id = comm.tgt_snode_id;
    Int src_nzblk_idx = comm.src_nzblk_idx;
    Int src_first_row = comm.src_first_row;

#ifdef _DEBUG_DELAY_
    logfileptr->OFS()<<"Picked { "<<src_snode_id<<" -> "<<tgt_snode_id<<" }"<<endl;
#endif

    if(tgt_snode_id>next_snode_id || is_last){

      Int iLocalSrc = (src_snode_id-1) / np +1 ;
      SuperNode<T> & prev_src_snode = *snodeColl[iLocalSrc -1];
      //this can be sent now

      Int iTarget = Mapping_->Map(tgt_snode_id-1,tgt_snode_id-1);
      if(iTarget != iam){
#ifdef _DEBUG_DELAY_
        logfileptr->OFS()<<"P"<<iam<<" has sent update from Supernode "<<prev_src_snode.Id()<<" to Supernode "<<tgt_snode_id<<endl;
        cout<<"P"<<iam<<" has sent update from Supernode "<<prev_src_snode.Id()<<" to Supernode "<<tgt_snode_id<<endl;
#endif

#ifdef _DEBUG_
        logfileptr->OFS()<<"Remote Supernode "<<tgt_snode_id<<" is updated by Supernode "<<prev_src_snode.Id()<<" rows "<<src_first_row/*<<" to "<<src_last_row*/<<" "<<src_nzblk_idx<<std::endl;
#endif

        NZBlockDesc & pivot_desc = prev_src_snode.GetNZBlockDesc(src_nzblk_idx);
        //Create a new Icomm buffer, serialize the contribution
        // in it and add it to the outgoing comm list
        Icomm * send_buffer = new Icomm();
        Serialize(*send_buffer,prev_src_snode,src_nzblk_idx,src_first_row);
        AddOutgoingComm(OutgoingSend,send_buffer);


        if( OutgoingSend.size() > maxIsend_){
          TIMER_START(SEND_MPI);
          MPI_Send(OutgoingSend.back()->front(),OutgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,CommEnv_->MPI_GetComm());
          TIMER_STOP(SEND_MPI);
          OutgoingSend.pop_back();
        }
        else{
          TIMER_START(SEND_MPI);
          MPI_Isend(OutgoingSend.back()->front(),OutgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,CommEnv_->MPI_GetComm(),&OutgoingSend.back()->Request);
          TIMER_STOP(SEND_MPI);
        }
      }

      //remove from the list
      MsgToSend.pop();
    }
    else{
      break;
    }
  }
}



template <typename T> void SupernodalMatrix<T>::SendDelayedMessagesUp(Int iLocalI, CommList & MsgToSend, AsyncComms & OutgoingSend, std::vector<SuperNode<T> *> & snodeColl){
//  typedef volatile LIBCHOLESKY::Int Int;

  if(snodeColl.empty() || MsgToSend.empty()) { return;}

  //Index of the last global snode to do
  Int last_snode_id = Xsuper_.m()-1;
  //Index of the last local supernode
  Int last_local_id = snodeColl.back()->Id();
  //Index of the last PROCESSED supernode
  Int prev_snode_id = iLocalI<=snodeColl.size()?snodeColl[iLocalI-1]->Id():last_local_id;
  //Index of the next local supernode
  Int next_snode_id = prev_snode_id>=last_local_id?last_snode_id+1:snodeColl[iLocalI]->Id();

  bool is_last = prev_snode_id>=last_local_id;


#ifdef _DEBUG_DELAY_
  {
    CommList tmp = MsgToSend;
    logfileptr->OFS()<<"Queue : ";
    while( tmp.size()>0){
      //Pull the highest priority message
      const DelayedComm & comm = tmp.top();
      Int src_snode_id = comm.src_snode_id;
      Int tgt_id = comm.tgt_snode_id;
      Int src_nzblk_idx = comm.src_nzblk_idx;
      Int src_first_row = comm.src_first_row;
      logfileptr->OFS()<<" { "<<src_snode_id<<" -> "<<tgt_id<<" }";
      tmp.pop();
    }
    logfileptr->OFS()<<endl;


  }
#endif
  while( MsgToSend.size()>0){
    //Pull the highest priority message
    const DelayedComm & comm = MsgToSend.top();
    Int src_snode_id = comm.src_snode_id;
    Int tgt_snode_id = comm.tgt_snode_id;
    Int src_nzblk_idx = comm.src_nzblk_idx;
    Int src_first_row = comm.src_first_row;

#ifdef _DEBUG_DELAY_
    logfileptr->OFS()<<"Picked { "<<src_snode_id<<" -> "<<tgt_snode_id<<" }"<<endl;
#endif

    if(tgt_snode_id < next_snode_id || is_last){
      Int iLocalSrc = (src_snode_id-1) / np +1 ;
      SuperNode<T> & prev_src_snode = *snodeColl[iLocalSrc -1];

      //this can be sent now
      Int iTarget = Mapping_->Map(tgt_snode_id-1,tgt_snode_id-1);
      if(iTarget != iam){
#ifdef _DEBUG_DELAY_
        logfileptr->OFS()<<"P"<<iam<<" has sent update from Supernode "<<prev_src_snode.Id()<<" to Supernode "<<tgt_snode_id<<endl;
        cout<<"P"<<iam<<" has sent update from Supernode "<<prev_src_snode.Id()<<" to Supernode "<<tgt_snode_id<<endl;
#endif

#ifdef _DEBUG_
        logfileptr->OFS()<<"Remote Supernode "<<tgt_snode_id<<" is updated by Supernode "<<prev_src_snode.Id()<<" rows "<<src_first_row/*<<" to "<<src_last_row*/<<" "<<src_nzblk_idx<<std::endl;
#endif
        NZBlockDesc & pivot_desc = prev_src_snode.GetNZBlockDesc(src_nzblk_idx);
        //Create a new Icomm buffer, serialize the contribution
        // in it and add it to the outgoing comm list
        Icomm * send_buffer = new Icomm();
        Serialize(*send_buffer,prev_src_snode,src_nzblk_idx,src_first_row);
        AddOutgoingComm(OutgoingSend,send_buffer);


        if( OutgoingSend.size() > maxIsend_){
          TIMER_START(SEND_MPI);
          MPI_Send(OutgoingSend.back()->front(),OutgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,CommEnv_->MPI_GetComm());
          TIMER_STOP(SEND_MPI);
          OutgoingSend.pop_back();
        }
        else{
          TIMER_START(SEND_MPI);
          MPI_Isend(OutgoingSend.back()->front(),OutgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,CommEnv_->MPI_GetComm(),&OutgoingSend.back()->Request);
          TIMER_STOP(SEND_MPI);
        }
      }
      //remove from the list
      MsgToSend.pop();
    }
    else{
      break;
    }
  }
}



#include "SupernodalMatrix_impl_deprecated.hpp"


} // namespace LIBCHOLESKY






#endif 
