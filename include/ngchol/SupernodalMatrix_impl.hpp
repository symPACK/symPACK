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





  template <typename T> void SupernodalMatrix<T>::Init(const DistSparseMatrix<T> & pMat, Int maxSnode,MAPCLASS & pMapping, Int maxIsend, Int maxIrecv, MPI_Comm & pComm ){
  //Create the CommEnvironment object if necessary
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

#ifdef _DEBUG_
    logfileptr->OFS()<<"Membership list is "<<SupMembership_<<std::endl;
    logfileptr->OFS()<<"xsuper "<<Xsuper_<<std::endl;
#endif

#ifdef RELAXED_SNODE
    Global_.RelaxSupernodes(ETree_, cc,SupMembership_, Xsuper_, maxSnode );
    Global_.SymbolicFactorizationRelaxed(ETree_,Order_,cc,Xsuper_,SupMembership_,xlindx_,lindx_);

#else
    Global_.SymbolicFactorization(ETree_,Order_,cc,Xsuper_,SupMembership_,xlindx_,lindx_);
#ifdef REFINED_SNODE
if(0){
    IntNumVec permRefined;
    IntNumVec newPerm(Size());
    Global_.RefineSupernodes(ETree_,Order_, SupMembership_, Xsuper_, xlindx_, lindx_, permRefined);


//newPerm = perm3;
//      for(Int i =0; i<perm3.m();++i){
//        newPerm[perm2[i]-1] = poperm[i];
//      }
//      logfileptr->OFS()<<"new perm == po perm 2 ?: "<<newPerm<<std::endl;
//      poperm2 = newPerm;
//
//      for(Int i =0; i<perm3.m();++i){
//        newPerm[perm3[i]-1] = perm2[i];
//      }
//      logfileptr->OFS()<<"perm2 and perm3: "<<newPerm<<std::endl;
//
//      perm3 = newPerm;
//    Perm_ = permRefined;
}
#endif

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


    GetUpdatingSupernodeCount(UpdateCount_,UpdateWidth_);

    Int iam = CommEnv_->MPI_Rank();
    Int np  = CommEnv_->MPI_Size();

    //copy

      Icomm buffer;
      std::vector<T> denseA;
//      std::vector<char> recv_buffer;
      Icomm recv_buffer;
    for(Int I=1;I<Xsuper_.m();I++){
      Int fc = Xsuper_(I-1);
      Int lc = Xsuper_(I)-1;
      Int iWidth = lc-fc+1;
      Int fi = xlindx_(I-1);
      Int li = xlindx_(I)-1;
      Int iHeight = li-fi+1;

      Int iDest = pMapping.Map(I-1,I-1);

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
      Int numColFirst = iSize_ / np;

      //copy the data from A into this Block structure
      for(Int col = fc;col<=lc;col++){
        //corresponding column in the unsorted matrix A
        Int orig_col = Order_.perm[col-1];

        //determine last contiguous col in the original matrix A
        Int lccol = col;
        while(++lccol<=lc){
          Int orig_lccol = Order_.perm[lccol-1];
          if(orig_lccol != orig_col+1){
            break;
          }
        }


        buffer.clear();
        buffer.resize(iHeight*(sizeof(Int) + sizeof(T)));
        //If the column comes later in the permuted matrix,
        //NZvals are already in the LT part of the original matrix
        //Otherwise we need to look at the structure of row orig_col

        Int iOwnerCol = std::min((orig_col-1)/numColFirst,np-1);

        //If I own the column in the original matrix
        //Copy the appropriate values in the string stream



        //NZvals are in row orig_col of original matrix
        std::vector<Int> isSender(np,0);

        Int prev_col = col;
        do{
          Int orig_prev_col = Order_.perm[prev_col-1];
          Int iOwnerCurCol = std::min((orig_prev_col-1)/numColFirst,np-1);

          Int local_col = (orig_prev_col-(numColFirst)*iOwnerCurCol);
          Int colbeg = Local_.colptr[local_col-1];
          Int colend = Local_.colptr[local_col]-1;
          const T * pdData = &pMat.nzvalLocal[0];

          for(Int rowidx = colbeg; rowidx<=colend; ++rowidx){
            Int orig_row = Local_.rowind[rowidx-1];
            Int row = Order_.invp[orig_row-1];
            //serialize the value
            if(orig_prev_col == orig_col){
              if(iam == iDest){
                isSender[iOwnerCurCol] = 1;
              }

              if(iam == iOwnerCurCol){
                if(row>=col){
                  buffer<<row;
                  assert(row>0);
                  T val = pdData[rowidx-1];
                  buffer<<val;
                }
              }
            }
            else if(row == col){
              if(iam == iDest){
                isSender[iOwnerCurCol] = 1;
              }

              if(iam == iOwnerCurCol){
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

////          //If I'm the destination, I should count the number of senders
////          if( iam == iDest ){
////
////          Int prev_col = col;
////          do{
////            Int orig_prev_col = Order_.perm[prev_col-1];
////            Int iOwnerCurCol = std::min((orig_prev_col-1)/numColFirst,np-1);
////
////              Int colbeg = Global_.expColptr[orig_prev_col-1];
////              Int colend = Global_.expColptr[orig_prev_col]-1;
////
////                for(Int rowidx = colbeg; rowidx<=colend; ++rowidx){
////                  Int orig_row = Global_.expRowind[rowidx-1];
////                  Int row = Order_.invp[orig_row-1];
////                  //Look only at the lower triangular part
////                  if(orig_row>= orig_prev_col){
////                    //serialize the value
////                    if(orig_prev_col == orig_col){
////                      if(row>=col){
////                        isSender[iOwnerCurCol] = 1;//true;      
////                      }
////                    }
////                    else if(row == col){
////                      isSender[iOwnerCurCol] = 1;//true;      
////                      break;
////                    }
////                  }
////                }
////
////            prev_col = ETree_.PostParent(prev_col-1);
////          }while(prev_col!=0);
////
//////              logfileptr->OFS()<<isSender<<endl;
////            }
////



#if 0
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
          //now I know if I need to send something and who I'm receiving from
          if(iam == iDest){

            denseA.resize(iSize_);
            std::vector<Icomm> recvBuffers(np);

            std::vector<MPI_Request> requestSizes(np,MPI_REQUEST_NULL);
            std::vector<MPI_Request> requestContents(np,MPI_REQUEST_NULL);
            std::vector<Int> sizes(np,0);
            for(Int psrc = 0; psrc<np; ++psrc){
              if(isSender[psrc]){
                if(psrc!=iam){
                  MPI_Irecv(&sizes[psrc],sizeof(Int),MPI_BYTE,psrc,2*col+0,CommEnv_->MPI_GetComm(),&requestSizes[psrc]);
                }
              }
            }

            //MPI_Waitall(np,&requestSizes[0],MPI_STATUSES_IGNORE);

            Int indx = MPI_UNDEFINED;
            MPI_Status recv_status;
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

#ifdef _DEBUG_
    for(Int I=1;I<Xsuper_.m();I++){
      Int src_first_col = Xsuper_(I-1);
      Int src_last_col = Xsuper_(I)-1;
      Int iOwner = Mapping_.Map(I-1,I-1);
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

  template <typename T> SupernodalMatrix<T>::SupernodalMatrix(const DistSparseMatrix<T> & pMat, Int maxSnode,MAPCLASS & pMapping, Int maxIsend, Int maxIrecv, MPI_Comm & pComm ){
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



#ifdef _DEBUG_UPDATES_
      logfileptr->OFS()<<"Supernode "<<s<<" updates: ";
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
      Int iOwner = Mapping_.Map(src_snode.Id()-1,src_snode.Id()-1);
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

    Int iOwner = Mapping_.Map(src_contrib->Id()-1,src_contrib->Id()-1);
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

      for(Int i=0; i<(tgt_local_lr - tgt_local_fr +1)*src_ncols;++i,++tgt,++src){
        *tgt+=*src;
      }
//      for(Int i=0; i<(tgt_local_lr - tgt_local_fr +1)*src_ncols;++i){
//        tgt[i]+=src[i];
//      }


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
#if 0
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
      Int iOwner = Mapping_.Map(I-1,I-1);
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

        contrib->Shrink();
      }
    }


    //forward-substitution phase
    //Sending contrib up the tree
    //Start from the leaves of the tree
    TIMER_START(SPARSE_FWD_SUBST);
    for(Int I=1;I<Xsuper_.m();I++){
      Int iOwner = Mapping_.Map(I-1,I-1);
      //If I own the column, factor it
      if( iOwner == iam ){
        Int iLocalI = (I-1) / np +1 ;
        SuperNode<T> * cur_snode = LocalSupernodes_[iLocalI-1];
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
        std::vector<char> src_blocks;
        std::vector<T> src_nzval;
        size_t max_bytes;
        Int nz_cnt;
        while(UpdatesToDo(I-1)>0){
          //receive children contrib
#ifdef _DEBUG_
          logfileptr->OFS()<<UpdatesToDo(I-1)<<" contribs left"<<endl;
#endif


    TIMER_START(RECV_MALLOC);
          if(src_blocks.size()==0){
            max_bytes = 3*sizeof(Int); 
            Int nrows = cur_snode->NRowsBelowBlock(0);
            Int ncols = nrhs;
            nz_cnt = nrows * ncols;

            Int nblocks = nrows;//std::max((Int)ceil(nrows/2)+1,cur_snode->NZBlockCnt());
            max_bytes += (nblocks)*sizeof(NZBlockDesc);
            max_bytes += nz_cnt*sizeof(T); 

            src_blocks.resize(max_bytes);
          }
    TIMER_STOP(RECV_MALLOC);

    TIMER_START(RECV_MPI);
          MPI_Status recv_status;
          int bytes_received = 0;

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


          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);

          Int src_snode_id = *(Int*)&src_blocks[0];
          Int src_nzblk_cnt = *(((Int*)&src_blocks[0])+1);
          NZBlockDesc * src_blocks_ptr = 
                    reinterpret_cast<NZBlockDesc*>(&src_blocks[2*sizeof(Int)]);
          Int src_nzval_cnt = *(Int*)(src_blocks_ptr + src_nzblk_cnt);
          T * src_nzval_ptr = (T*)((Int*)(src_blocks_ptr + src_nzblk_cnt)+1);
    TIMER_STOP(RECV_MPI);

          //Create the dummy supernode for that data
          SuperNode<T> dist_contrib(src_snode_id,1,nrhs, src_blocks_ptr, src_nzblk_cnt, src_nzval_ptr, src_nzval_cnt);

#ifdef PROBE_FIRST
          if(doabort){

            abort();
          }
#endif

#ifdef _DEBUG_
          logfileptr->OFS()<<"RECV contrib of Supernode "<<dist_contrib.Id()<<std::endl;
#endif

          forward_update(&dist_contrib,contrib);

          --UpdatesToDo(I-1);

        }


        if(UpdatesToDo(I-1)==0){

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
                    Int srcRow = Order_.perm[diag_desc.GIndex+kk-1];
//logfileptr->OFS()<<"                      B = "<<B(Order_.perm[diag_desc.GIndex+kk-1]-1,j)<<std::endl;
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

            Int iTarget = Mapping_.Map(parent_snode_id-1,parent_snode_id-1);

            if(iTarget!=iam){
#ifdef _DEBUG_
              logfileptr->OFS()<<"Remote Supernode "<<parent_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
#endif





        Int tgt_first_col = Xsuper_(parent_snode_id-1);
              Int tgt_last_col = Xsuper_(parent_snode_id)-1;
              Int src_nzblk_idx = 1;
              NZBlockDesc & pivot_desc = contrib->GetNZBlockDesc(src_nzblk_idx);

              Int src_first_row = pivot_desc.GIndex;
              Int local_first_row = 0;
              Int nzblk_cnt = contrib->NZBlockCnt()-src_nzblk_idx;

              T * nzval_ptr = contrib->GetNZval(pivot_desc.Offset
                  +local_first_row*nrhs);

              Int nz_cnt = (contrib->NRowsBelowBlock(src_nzblk_idx)
                  - local_first_row )*nrhs;





    TIMER_START(SEND_MALLOC);
              AddOutgoingComm(outgoingSend, contrib->Id(), nrhs, src_first_row, pivot_desc, nzblk_cnt, nzval_ptr, nz_cnt);
    TIMER_STOP(SEND_MALLOC);


    TIMER_START(SEND_MPI);
//                  assert(iTarget<np);
    MPI_Send(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,parent_snode_id,CommEnv_->MPI_GetComm());
    TIMER_STOP(SEND_MPI);
      outgoingSend.pop_back();


#ifdef _DEBUG_            
              logfileptr->OFS()<<"     Send contribution "<<I<<" to Supernode "<<parent_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
              logfileptr->OFS()<<"Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
              logfileptr->OFS()<<"Sending "<<nzblk_cnt*sizeof(NZBlockDesc)+nz_cnt*sizeof(T) + 3*sizeof(Int)<<" bytes"<< std::endl;
#endif
            }
            else{
#ifdef _DEBUG_
              logfileptr->OFS()<<"Local Supernode "<<parent_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
#endif
              Int iLocalJ = (parent_snode_id-1) / np +1 ;
              LocalUpdates[iLocalJ-1].push(I);
            }
          }
        }
      }

#ifdef _CHECK_RESULT_SEQ_
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
          tmp3(Order_.perm[i]-1,j) = tmp(i,j);
          tmpFwd(Order_.perm[i]-1,j) = forwardSol(i,j);
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

      //MPI_Barrier(CommEnv_->MPI_GetComm());
    }
    TIMER_STOP(SPARSE_FWD_SUBST);




    //Back-substitution phase

    TIMER_START(SPARSE_BACK_SUBST);

    //start from the root of the tree
    for(Int I=Xsuper_.m()-1;I>=1;--I){
      Int iOwner = Mapping_.Map(I-1,I-1);
      //If I own the column, factor it
      if( iOwner == iam ){
        Int iLocalI = (I-1) / np +1 ;
        SuperNode<T> * cur_snode = LocalSupernodes_[iLocalI-1];
        SuperNode<T> * contrib = Contributions_[iLocalI-1];

        Int parent = ETree_.PostParent(cur_snode->LastCol()-1);
        //Extend the contribution.

        std::vector<Int> isBlockUpdated(contrib->NZBlockCnt(),0);

        if(parent!=0){
          Int parent_snode_id = SupMembership_[parent-1];
          Int iTarget = Mapping_.Map(parent_snode_id-1,parent_snode_id-1);
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

//            assert(iTarget!=iam);

            //Receive parent contrib
            std::vector<char> src_blocks;
            std::vector<T> src_nzval;


            //MPI_Recv
            MPI_Status recv_status;

            //Receive the size of the blocks array
            Int max_bytes = 0;
            MPI_Recv(&max_bytes,sizeof(Int),MPI_BYTE,iTarget,I,CommEnv_->MPI_GetComm(),&recv_status);
            src_blocks.resize(max_bytes);


            //receive the index array
            MPI_Recv(&src_blocks[0],max_bytes,MPI_BYTE,iTarget,I,CommEnv_->MPI_GetComm(),&recv_status);
            //there an aditional integer storing the src id which needs to be removed
            Int src_nzblk_cnt = (max_bytes - sizeof(Int) ) / sizeof(NZBlockDesc);



            //Receive the size of the nzval array
            Int nz_cnt = 0;
            MPI_Recv(&nz_cnt,sizeof(Int),MPI_BYTE,iTarget,I,CommEnv_->MPI_GetComm(),&recv_status);
            src_nzval.resize(nz_cnt);

            //receive the nzval array
            MPI_Recv(&src_nzval[0],nz_cnt*sizeof(T),MPI_BYTE,iTarget,I,CommEnv_->MPI_GetComm(),&recv_status);

            Int src_snode_id = *(Int*)&src_blocks[0];
            NZBlockDesc * src_blocks_ptr = 
              reinterpret_cast<NZBlockDesc*>(&src_blocks[sizeof(Int)]);
            T * src_nzval_ptr = &src_nzval[0];

            //Create the dummy supernode for that data
            dist_contrib = new SuperNode<T>(src_snode_id,1,nrhs, src_blocks_ptr, src_nzblk_cnt, src_nzval_ptr, nz_cnt);

#ifdef _DEBUG_
            logfileptr->OFS()<<"RECV contrib of Supernode "<<dist_contrib->Id()<<std::endl;
#endif

            back_update(dist_contrib,contrib);
            delete dist_contrib;
          }
        }

        //now compute MY contribution



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
//            Int firstCol = Xsuper_[child_snode_id-1];
//            for(Int col = colIdx; col>=firstCol; --col){
//            }

            Int parent = ETree_.PostParent(colIdx-1);
            if(parent!=0){
            if(SupMembership_[parent-1]==cur_snode->Id()){
              Int iTarget = Mapping_.Map(child_snode_id-1,child_snode_id-1);

              if(iTarget!=iam){



#ifdef _DEBUG_
                logfileptr->OFS()<<"Remote Supernode "<<child_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
#endif



                std::vector<char> * pNewDesc = new std::vector<char>();

                Int tgt_first_col = Xsuper_(child_snode_id-1);
                Int tgt_last_col = Xsuper_(child_snode_id)-1;

                //Send
                Int src_nzblk_idx = 0;
                NZBlockDesc & pivot_desc = contrib->GetNZBlockDesc(src_nzblk_idx);

                Int src_first_row = pivot_desc.GIndex;
                Int local_first_row = 0;
                Int nzblk_cnt = contrib->NZBlockCnt()-src_nzblk_idx;
                pNewDesc->resize(sizeof(Int) + nzblk_cnt*sizeof(NZBlockDesc));

                char * send_ptr = &(*pNewDesc)[0];

                Int * id_ptr = reinterpret_cast<Int *>(&(*pNewDesc)[0]);
                NZBlockDesc * block_desc_ptr = 
                  reinterpret_cast<NZBlockDesc *>(&(*pNewDesc)[sizeof(Int)]);

                *id_ptr = contrib->Id();
                //copy the block descriptors
                std::copy(&pivot_desc,&pivot_desc+nzblk_cnt, block_desc_ptr);
                //change the first one
                block_desc_ptr->Offset += (src_first_row - block_desc_ptr->GIndex)*nrhs;
                block_desc_ptr->GIndex = src_first_row;

                Int bytes_size = nzblk_cnt*sizeof(NZBlockDesc)+sizeof(Int);


                T * nzval_ptr = contrib->GetNZval(pivot_desc.Offset
                    +local_first_row*nrhs);

                Int nz_cnt = (contrib->NRowsBelowBlock(src_nzblk_idx)
                    - local_first_row )*nrhs;

                //send the block descriptors
//                  assert(iTarget<np);
                MPI_Send(&bytes_size,sizeof(bytes_size),MPI_BYTE,iTarget,child_snode_id,CommEnv_->MPI_GetComm());
                MPI_Send(send_ptr,bytes_size, MPI_BYTE,iTarget,child_snode_id,CommEnv_->MPI_GetComm());
                //send the nzvals
                MPI_Send(&nz_cnt,sizeof(nz_cnt),MPI_BYTE,iTarget,child_snode_id,CommEnv_->MPI_GetComm());
                MPI_Send(nzval_ptr,nz_cnt*sizeof(T),MPI_BYTE,iTarget,child_snode_id,CommEnv_->MPI_GetComm());

                delete pNewDesc;

#ifdef _DEBUG_            
                logfileptr->OFS()<<"     Send contribution "<<I<<" to Supernode "<<child_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
                logfileptr->OFS()<<"Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
                logfileptr->OFS()<<"Sending "<<nzblk_cnt*sizeof(NZBlockDesc)<<" and "<<nz_cnt*sizeof(T)<<" bytes during BS"<<std::endl;
#endif
              }
              else{

#ifdef _DEBUG_
                logfileptr->OFS()<<"Local Supernode "<<child_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
#endif
                Int iLocalJ = (child_snode_id-1) / np +1 ;
                LocalUpdates[iLocalJ-1].push(I);
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

//      MPI_Barrier(CommEnv_->MPI_GetComm());

#ifdef _CHECK_RESULT_SEQ_
      {
//      MPI_Barrier(CommEnv_->MPI_GetComm());
//      NumMat<T> tmp = B;
//      GetSolution(tmp);
//            logfileptr->OFS()<<tmp<<std::endl;
//      Int nrows = 0;
//      for(Int ii=1; ii<=I;++ii){ nrows+= Xsuper_[ii] - Xsuper_[ii-1];}
//
//      NumMat<T> tmp2 = tmp;
//      blas::Axpy(tmp.m()*tmp.n(),-1.0,&forwardSol(0,0),1,&tmp2(0,0),1);
//      double norm = lapack::Lange('F',nrows,tmp.n(),&tmp2(0,0),tmp.m());
//      logfileptr->OFS()<<"Norm after SuperNode "<<I<<" is "<<norm<<std::endl; 
//
//        if(abs(norm)>=1e-1){
//          for(Int i = 0;i<tmp.m();++i){
//            logfileptr->OFS()<<forwardSol(i,0)<<"       "<<tmp(i,0)<<std::endl;
//          }
//        }
      }
#endif







    }
    TIMER_STOP(SPARSE_BACK_SUBST);

    MPI_Barrier(CommEnv_->MPI_GetComm());
    TIMER_STOP(SPARSE_SOLVE);

}
#else
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
      Int iOwner = Mapping_.Map(I-1,I-1);
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


    std::vector<char> src_blocks;


    //forward-substitution phase
    //Sending contrib up the tree
    //Start from the leaves of the tree
    TIMER_START(SPARSE_FWD_SUBST);

    volatile Int I =1;
    volatile Int iLocalI =1;
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


          TIMER_START(RECV_MALLOC);
          max_bytes = 3*sizeof(Int); 
          Int nrows = cur_snode->NRowsBelowBlock(0);
          Int ncols = nrhs;
          nz_cnt = nrows * ncols;

          Int nblocks = nrows;//std::max((Int)ceil(nrows/2)+1,cur_snode->NZBlockCnt());
          max_bytes += (nblocks)*sizeof(NZBlockDesc);
          max_bytes += nz_cnt*sizeof(T); 

          src_blocks.resize(max_bytes);
          TIMER_STOP(RECV_MALLOC);

          TIMER_START(RECV_MPI);
          MPI_Status recv_status;
          int bytes_received = 0;

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

#if 1
          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);

          Int src_snode_id = *(Int*)&src_blocks[0];
          Int src_nzblk_cnt = *(((Int*)&src_blocks[0])+1);
          NZBlockDesc * src_blocks_ptr = 
            reinterpret_cast<NZBlockDesc*>(&src_blocks[2*sizeof(Int)]);
          Int src_nzval_cnt = *(Int*)(src_blocks_ptr + src_nzblk_cnt);
          T * src_nzval_ptr = (T*)((Int*)(src_blocks_ptr + src_nzblk_cnt)+1);
          TIMER_STOP(RECV_MPI);

          //Create the dummy supernode for that data
          SuperNode<T> dist_contrib(src_snode_id,1,nrhs, src_blocks_ptr, src_nzblk_cnt, src_nzval_ptr, src_nzval_cnt);
#else

          SuperNode<T> dist_contrib;
          Deserialize(&src_blocks[0],dist_contrib);
#endif
#ifdef PROBE_FIRST
          if(doabort){

            abort();
          }
#endif

#ifdef _DEBUG_
          logfileptr->OFS()<<"RECV contrib of Supernode "<<dist_contrib.Id()<<std::endl;
#endif

          forward_update(&dist_contrib,contrib);

          --UpdatesToDo(I-1);

        }

        assert(UpdatesToDo(I-1)==0);

        if(UpdatesToDo(I-1)==0){
          logfileptr->OFS()<<"FORW Processing contrib "<<I<<std::endl;

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

            Int iTarget = Mapping_.Map(parent_snode_id-1,parent_snode_id-1);

            if(iTarget!=iam){







#ifdef _DEBUG_
              logfileptr->OFS()<<"Remote Supernode "<<parent_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
#endif





              Int tgt_first_col = Xsuper_(parent_snode_id-1);
              Int tgt_last_col = Xsuper_(parent_snode_id)-1;
              Int src_nzblk_idx = 1;
              NZBlockDesc & pivot_desc = contrib->GetNZBlockDesc(src_nzblk_idx);

              Int src_first_row = pivot_desc.GIndex;
              Int local_first_row = 0;
              Int nzblk_cnt = contrib->NZBlockCnt()-src_nzblk_idx;

              T * nzval_ptr = contrib->GetNZval(pivot_desc.Offset
                  +local_first_row*nrhs);

              Int nz_cnt = (contrib->NRowsBelowBlock(src_nzblk_idx)
                  - local_first_row )*nrhs;


              bool isSkipped= false;

              volatile Int next_local_contrib = (iLocalI < Contributions_.size())?Contributions_[iLocalI]->Id():Xsuper_.m();
                  if(next_local_contrib< parent_snode_id){
                  //need to push the prev src_last_row
                  ContribsToSend.push(DelayedComm(contrib->Id(),parent_snode_id,1,src_first_row));
#ifdef _DEBUG_DELAY_
                  //                    logfileptr->OFS()<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<tgt_snode_id<<" from row "<<src_first_row<<endl;
                                      cout<<"P"<<iam<<" has delayed update from Contrib "<<I<<" to "<<parent_snode_id<<" from row "<<src_first_row<<endl;
#endif
                    isSkipped= true;
                }
              
              if(!isSkipped){

//               if(!ContribsToSend.empty()){
//                  if(ContribsToSend.top().tgt_snode_id<parent_snode_id){
//                    SendDelayedMessages(I,ContribsToSend,outgoingSend,Contributions_);
//                  }
//                }




                TIMER_START(SEND_MALLOC);
                AddOutgoingComm(outgoingSend, contrib->Id(), nrhs, src_first_row, pivot_desc, nzblk_cnt, nzval_ptr, nz_cnt);
                TIMER_STOP(SEND_MALLOC);


                TIMER_START(SEND_MPI);
                //                  assert(iTarget<np);



                if( outgoingSend.size() > maxIsend_){
                  logfileptr->OFS()<<"Remote Supernode "<<parent_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
                  MPI_Send(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,parent_snode_id,CommEnv_->MPI_GetComm());
                  TIMER_STOP(SEND_MPI);
                  outgoingSend.pop_back();

                  logfileptr->OFS()<<"Remote Supernode "<<parent_snode_id<<" got the contribution of Supernode "<<I<<std::endl;
                }
                else{

                  TIMER_START(SEND_MPI);
                  MPI_Isend(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,parent_snode_id,CommEnv_->MPI_GetComm(),&outgoingSend.back()->Request);
                  TIMER_STOP(SEND_MPI);


                }
#ifdef _DEBUG_            
                logfileptr->OFS()<<"     Send contribution "<<I<<" to Supernode "<<parent_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
                logfileptr->OFS()<<"Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
                logfileptr->OFS()<<"Sending "<<nzblk_cnt*sizeof(NZBlockDesc)+nz_cnt*sizeof(T) + 3*sizeof(Int)<<" bytes"<< std::endl;
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

#ifdef _CHECK_RESULT_SEQ_
      {
//      MPI_Barrier(CommEnv_->MPI_GetComm());
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

          }


      SendDelayedMessages(iLocalI,ContribsToSend,outgoingSend,Contributions_);

          ++iLocalI;
      //MPI_Barrier(CommEnv_->MPI_GetComm());
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
          Int iTarget = Mapping_.Map(parent_snode_id-1,parent_snode_id-1);
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


            //MPI_Recv
            MPI_Status recv_status;



#if 1
            Int bytes_received = 0;
            MPI_Probe(iTarget,I,CommEnv_->MPI_GetComm(),&recv_status);
            MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
            src_blocks.resize(bytes_received);

            MPI_Recv(&src_blocks[0],bytes_received,MPI_BYTE,iTarget,I,CommEnv_->MPI_GetComm(),&recv_status);



            Int src_snode_id = *(Int*)&src_blocks[0];
            Int src_nzblk_cnt = *(((Int*)&src_blocks[0])+1);
            NZBlockDesc * src_blocks_ptr = 
              reinterpret_cast<NZBlockDesc*>(&src_blocks[2*sizeof(Int)]);
            Int src_nzval_cnt = *(Int*)(src_blocks_ptr + src_nzblk_cnt);
            T * src_nzval_ptr = (T*)((Int*)(src_blocks_ptr + src_nzblk_cnt)+1);




            dist_contrib = new SuperNode<T>(src_snode_id,1,nrhs, src_blocks_ptr, src_nzblk_cnt, src_nzval_ptr, src_nzval_cnt);


            //            dist_contrib = new SuperNode<T>();
            //            Deserialize(&src_blocks[0],*dist_contrib,nrhs); 


#else
            std::vector<T> src_nzval;
            //Receive the size of the blocks array
            Int max_bytes = 0;
            MPI_Recv(&max_bytes,sizeof(Int),MPI_BYTE,iTarget,I,CommEnv_->MPI_GetComm(),&recv_status);
            src_blocks.resize(max_bytes);


            //receive the index array
            MPI_Recv(&src_blocks[0],max_bytes,MPI_BYTE,iTarget,I,CommEnv_->MPI_GetComm(),&recv_status);
            //there an aditional integer storing the src id which needs to be removed
            Int src_nzblk_cnt = (max_bytes - sizeof(Int) ) / sizeof(NZBlockDesc);



            //Receive the size of the nzval array
            Int nz_cnt = 0;
            MPI_Recv(&nz_cnt,sizeof(Int),MPI_BYTE,iTarget,I,CommEnv_->MPI_GetComm(),&recv_status);
            src_nzval.resize(nz_cnt);

            logfileptr->OFS()<<"RECVing contrib to Supernode "<<contrib->Id()<<std::endl;
            //receive the nzval array
            MPI_Recv(&src_nzval[0],nz_cnt*sizeof(T),MPI_BYTE,iTarget,I,CommEnv_->MPI_GetComm(),&recv_status);

            Int src_snode_id = *(Int*)&src_blocks[0];
            NZBlockDesc * src_blocks_ptr = 
              reinterpret_cast<NZBlockDesc*>(&src_blocks[sizeof(Int)]);
            T * src_nzval_ptr = &src_nzval[0];

            //Create the dummy supernode for that data
            dist_contrib = new SuperNode<T>(src_snode_id,1,nrhs, src_blocks_ptr, src_nzblk_cnt, src_nzval_ptr, nz_cnt);
#endif


            //#ifdef _DEBUG_
            logfileptr->OFS()<<"RECV contrib of Supernode "<<dist_contrib->Id()<<std::endl;
            //#endif

            back_update(dist_contrib,contrib);
            delete dist_contrib;
          }
        }

        //now compute MY contribution


        logfileptr->OFS()<<"BACK Processing contrib "<<I<<std::endl;

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
            //            Int firstCol = Xsuper_[child_snode_id-1];
            //            for(Int col = colIdx; col>=firstCol; --col){
            //            }

            Int parent = ETree_.PostParent(colIdx-1);
            if(parent!=0){
              if(SupMembership_[parent-1]==cur_snode->Id()){
                Int iTarget = Mapping_.Map(child_snode_id-1,child_snode_id-1);

                if(iTarget!=iam){

                  bool isSkipped= false;

                  volatile Int next_local_contrib = (iLocalI >1)?Contributions_[iLocalI-2]->Id():0;
                  if(next_local_contrib > child_snode_id){
                    //need to push the prev src_last_row

                    //Send
                    Int src_nzblk_idx = 0;
                    NZBlockDesc & pivot_desc = contrib->GetNZBlockDesc(src_nzblk_idx);

                    Int src_first_row = pivot_desc.GIndex;
                    ContribsToSend.push(DelayedComm(contrib->Id(),child_snode_id,1,src_first_row));
#ifdef _DEBUG_DELAY_
                    //                    logfileptr->OFS()<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<tgt_snode_id<<" from row "<<src_first_row<<endl;
                    cout<<"P"<<iam<<" has delayed update from Contrib "<<I<<" to "<<child_snode_id<<" from row "<<src_first_row<<endl;
#endif
                    isSkipped= true;
                  }

                  if(!isSkipped){




#ifdef _DEBUG_
                    logfileptr->OFS()<<"Remote Supernode "<<child_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
#endif



                    std::vector<char> * pNewDesc = new std::vector<char>();

                    Int tgt_first_col = Xsuper_(child_snode_id-1);
                    Int tgt_last_col = Xsuper_(child_snode_id)-1;


                    Int src_nzblk_idx = 0;
                    NZBlockDesc & pivot_desc = contrib->GetNZBlockDesc(src_nzblk_idx);
                    Int src_first_row = pivot_desc.GIndex;
                    Int nzblk_cnt = contrib->NZBlockCnt()-src_nzblk_idx;
                    T * nzval_ptr = contrib->GetNZval(pivot_desc.Offset);
                    Int nz_cnt = (contrib->NRowsBelowBlock(src_nzblk_idx))*nrhs;


                    AddOutgoingComm(outgoingSend,contrib->Id(),nrhs,src_first_row,pivot_desc,nzblk_cnt,nzval_ptr,nz_cnt);

                    if( outgoingSend.size() > maxIsend_){
                      //logfileptr->OFS()<<"Remote Supernode "<<parent_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
                      MPI_Send(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,child_snode_id,CommEnv_->MPI_GetComm());
                      TIMER_STOP(SEND_MPI);
                      outgoingSend.pop_back();

                      //logfileptr->OFS()<<"Remote Supernode "<<parent_snode_id<<" got the contribution of Supernode "<<I<<std::endl;
                    }
                    else{

                      TIMER_START(SEND_MPI);
                      MPI_Isend(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,child_snode_id,CommEnv_->MPI_GetComm(),&outgoingSend.back()->Request);
                      TIMER_STOP(SEND_MPI);


                    }

#if 0
                    //Send
                    Int src_nzblk_idx = 0;
                    NZBlockDesc & pivot_desc = contrib->GetNZBlockDesc(src_nzblk_idx);

                    Int src_first_row = pivot_desc.GIndex;
                    Int local_first_row = 0;
                    Int nzblk_cnt = contrib->NZBlockCnt()-src_nzblk_idx;
                    pNewDesc->resize(sizeof(Int) + nzblk_cnt*sizeof(NZBlockDesc));

                    char * send_ptr = &(*pNewDesc)[0];

                    Int * id_ptr = reinterpret_cast<Int *>(&(*pNewDesc)[0]);
                    NZBlockDesc * block_desc_ptr = 
                      reinterpret_cast<NZBlockDesc *>(&(*pNewDesc)[sizeof(Int)]);

                    *id_ptr = contrib->Id();
                    //copy the block descriptors
                    std::copy(&pivot_desc,&pivot_desc+nzblk_cnt, block_desc_ptr);
                    //change the first one
                    block_desc_ptr->Offset += (src_first_row - block_desc_ptr->GIndex)*nrhs;
                    block_desc_ptr->GIndex = src_first_row;

                    Int bytes_size = nzblk_cnt*sizeof(NZBlockDesc)+sizeof(Int);


                    T * nzval_ptr = contrib->GetNZval(pivot_desc.Offset
                        +local_first_row*nrhs);

                    Int nz_cnt = (contrib->NRowsBelowBlock(src_nzblk_idx)
                        - local_first_row )*nrhs;

                    //send the block descriptors
                    //                  assert(iTarget<np);

                    logfileptr->OFS()<<"     Sending contribution "<<I<<" to Supernode "<<child_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
                    MPI_Send(&bytes_size,sizeof(bytes_size),MPI_BYTE,iTarget,child_snode_id,CommEnv_->MPI_GetComm());
                    MPI_Send(send_ptr,bytes_size, MPI_BYTE,iTarget,child_snode_id,CommEnv_->MPI_GetComm());
                    //send the nzvals
                    MPI_Send(&nz_cnt,sizeof(nz_cnt),MPI_BYTE,iTarget,child_snode_id,CommEnv_->MPI_GetComm());
                    MPI_Send(nzval_ptr,nz_cnt*sizeof(T),MPI_BYTE,iTarget,child_snode_id,CommEnv_->MPI_GetComm());
                    logfileptr->OFS()<<"     Sent contribution "<<I<<" to Supernode "<<child_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;


                    delete pNewDesc;
#endif

#ifdef _DEBUG_            
                    logfileptr->OFS()<<"     Send contribution "<<I<<" to Supernode "<<child_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
                    logfileptr->OFS()<<"Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
                    logfileptr->OFS()<<"Sending "<<nzblk_cnt*sizeof(NZBlockDesc)<<" and "<<nz_cnt*sizeof(T)<<" bytes during BS"<<std::endl;
#endif
                  }
                }
                else{

#ifdef _DEBUG_
                  logfileptr->OFS()<<"Local Supernode "<<child_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
#endif
                  Int iLocalJ = (child_snode_id-1) / np +1 ;
                  LocalUpdates[iLocalJ-1].push(I);
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

      SendDelayedMessages(iLocalI,ContribsToSend,outgoingSend,Contributions_,true);
      --iLocalI;
    }
    TIMER_STOP(SPARSE_BACK_SUBST);

    MPI_Barrier(CommEnv_->MPI_GetComm());
    TIMER_STOP(SPARSE_SOLVE);

}
#endif

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
      Int iOwner = Mapping_.Map(I-1,I-1);
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
      Int iOwner = Mapping_.Map(I-1,I-1);
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

////
////
////    //Print B
#ifdef _DEBUG_
//    logfileptr->OFS()<<"Solution is "<<B<<std::endl;
#endif
////

      MPI_Barrier(CommEnv_->MPI_GetComm());
}




#include "SupernodalMatrix_impl_deprecated.hpp"


} // namespace LIBCHOLESKY






#endif 
