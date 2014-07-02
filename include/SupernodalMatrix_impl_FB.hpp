#ifndef _SUPERNODAL_MATRIX_IMPL_FB_HPP_
#define _SUPERNODAL_MATRIX_IMPL_FB_HPP_

template <typename T> void SupernodalMatrix<T>::FanBoth( MPI_Comm & pComm ){

  TIMER_START(FACTORIZATION_FO);


  MPI_Comm_rank(pComm, &iam);
  MPI_Comm_size(pComm, &np);
  IntNumVec UpdatesToDo,LastUpdate;
  FBGetUpdateCount(UpdatesToDo,LastUpdate);
  IntNumVec AggregatesToRecv = UpdateCount_;


  CommList UpdatesToSend; 
  CommList UpdatesToRecv; 
  CommList FactorsToSend; 
  AsyncComms outgoingSend;



  std::vector<std::queue<LocalUpdate> > LocalUpdates(LocalSupernodes_.size());

  //std::unordered_map<Int,SuperNode*> aggVectors(Xsuper.m()-1); 
  std::vector< SuperNode<T> * > aggVectors(Xsuper_.m()-1,NULL);


#ifdef UPDATE_LIST
  std::list<SnodeUpdate> updates;
#endif

  Int maxwidth = 0;
  for(Int i = 1; i<Xsuper_.m(); ++i){
    Int width =Xsuper_(i) - Xsuper_(i-1);
    if(width>=maxwidth){
      maxwidth = width;
    }
  }

  NumMat<T> tmpBuf(iSize_,maxwidth);


  //dummy right looking cholesky factorization
  Int I =1;
  Int iLocalI = -1;
  while(I<Xsuper_.m() || !FactorsToSend.empty() || !outgoingSend.empty()){

//////    //Check for completion of outgoing communication
//////    if(!outgoingSend.empty()){
//////      std::list<Icomm *>::iterator it = outgoingSend.begin();
//////      while(it != outgoingSend.end()){
//////        int flag = 0;
//////        MPI_Test(&(*it)->Request,&flag,MPI_STATUS_IGNORE);
//////        if(flag){
//////          it = outgoingSend.erase(it);
//////        }
//////        else{
//////          it++;
//////        }
//////      }
//////    }
//////
//////    //process some of the delayed send
//////    SendDelayedMessages(I,FactorsToSend,outgoingSend);
////////    SendDelayedMessages2(I,AggregatesToSend,outgoingSend);
//////

    if(I<Xsuper_.m()){
      Int src_first_col = Xsuper_(I-1);
      Int src_last_col = Xsuper_(I)-1;
      Int iOwner = Mapping_.Map(I-1,I-1);
      //If I own the column, factor it
      if( iOwner == iam ){
        iLocalI = (I-1) / np +1 ;
        SuperNode<T> & src_snode = *LocalSupernodes_[iLocalI -1];
#ifdef _DEBUG_
        logfileptr->OFS()<<"Processing Supernode "<<I<<std::endl;
#endif


        //Aggregate all the updates

//////        //Do all my updates (Local and remote)
//////        //Local updates
//////        while(!LocalUpdates[iLocalI-1].empty()){
//////          Int src_snode_id = LocalUpdates[iLocalI-1].front();
//////          LocalUpdates[iLocalI-1].pop();
//////          Int src_nzblk_idx = LocalUpdates[iLocalI-1].front();
//////          LocalUpdates[iLocalI-1].pop();
//////          Int src_first_row = LocalUpdates[iLocalI-1].front();
//////          LocalUpdates[iLocalI-1].pop();
//////
//////          SuperNode<T> & local_src_snode = *LocalSupernodes_[(src_snode_id-1) / np];
//////
//////#ifdef _DEBUG_
//////          logfileptr->OFS()<<"LOCAL Supernode "<<src_snode.Id()<<" is updated by Supernode "<<src_snode_id<<" from row "<<src_first_row<<" "<<src_nzblk_idx<<std::endl;
//////#endif
//////
//////#ifdef SINGLE_BLAS
//////          UpdateSuperNode(local_src_snode,src_snode,src_nzblk_idx,tmpBuf, src_first_row);
//////#else
//////          UpdateSuperNode(local_src_snode,src_snode,src_nzblk_idx, src_first_row);
//////#endif
//////          //        logfileptr->OFS()<<"After "<<src_snode<<std::endl;
//////
//////          --UpdatesToDo(I-1);
//////#ifdef _DEBUG_
//////          logfileptr->OFS()<<UpdatesToDo(I-1)<<" updates left"<<endl;
//////#endif
//////        }
//////
        //Remote Aggregates
        std::vector<T> src_nzval;
        std::vector<char> src_blocks;

        Int nz_cnt;
        Int max_bytes;
        while(AggregatesToRecv(I-1)>0){
#ifdef _DEBUG_
          logfileptr->OFS()<<UpdatesToDo(I-1)<<" updates left"<<endl;
#endif


          TIMER_START(RECV_MALLOC);
          if(src_blocks.size()==0){
            max_bytes = 3*sizeof(Int); 
            //The upper bound must be of the width of the "largest" child
#ifdef _DEBUG_
            logfileptr->OFS()<<"Maximum width is "<<UpdateWidth_(I-1)<<std::endl;
#endif

            Int nrows = src_snode.NRowsBelowBlock(0);
            Int ncols = UpdateWidth_(I-1);
            nz_cnt = nrows * ncols;

            max_bytes += (std::max((Int)ceil(nrows/2)+1,src_snode.NZBlockCnt()))*sizeof(NZBlockDesc);
            max_bytes += nz_cnt*sizeof(T); 

            src_blocks.resize(max_bytes);
          }
          TIMER_STOP(RECV_MALLOC);

          TIMER_START(RECV_MPI);
          MPI_Status recv_status;
          int bytes_received = 0;

#ifdef PROBE_FIRST
          MPI_Probe(MPI_ANY_SOURCE,I,pComm,&recv_status);
          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);

          bool doabort = false;
          int prev_size = 0;
          if(src_blocks.size()<bytes_received){
            prev_size = src_blocks.size();
            doabort = true;
            //receive anyway
            src_blocks.resize(bytes_received);
          }
#endif


#ifdef PROBE_FIRST
          MPI_Recv(&src_blocks[0],max_bytes,MPI_BYTE,recv_status.MPI_SOURCE,I,pComm,&recv_status);
#else
          //receive the index array
          MPI_Recv(&src_blocks[0],max_bytes,MPI_BYTE,MPI_ANY_SOURCE,I,pComm,&recv_status);
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
          SuperNode<T> dist_src_snode(src_snode_id,Xsuper_[src_snode_id-1],Xsuper_[src_snode_id]-1, Size(), src_blocks_ptr, src_nzblk_cnt, src_nzval_ptr, src_nzval_cnt);

#ifdef PROBE_FIRST
          if(doabort){
            cout<<"We have a problem !!!! on P"<<iam<<"\n";
            gdb_lock();

            abort();
          }
#endif
#ifdef _DEBUG_
          logfileptr->OFS()<<"RECV Supernode "<<dist_src_snode.Id()<<std::endl;
#endif

#ifdef UPDATE_LIST
          TIMER_START(UPDATE_ANCESTORS);
          FindUpdates(dist_src_snode,updates);
          //now traverse the list
          for(std::list<SnodeUpdate>::iterator it = updates.begin(); it!=updates.end();it++){
            Int tgt_snode_id = it->tgt_snode_id;
            Int src_first_row = it->src_fr;
            Int src_nzblk_idx = dist_src_snode.FindBlockIdx(src_first_row);
            Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);

            if(iTarget == iam){
#ifdef _DEBUG_
              logfileptr->OFS()<<"RECV Supernode "<<tgt_snode_id<<" is updated by Supernode "<<dist_src_snode.Id()<<" from row "<<src_first_row<<" "<<src_nzblk_idx<<std::endl;
#endif

              Int iLocalJ = (tgt_snode_id-1) / np +1 ;
              SuperNode<T> & tgt_snode = *LocalSupernodes_[iLocalJ -1];

              AggregateSuperNode(dist_src_snode,tgt_snode,src_nzblk_idx, src_first_row);

              --AggregatesToRecv(tgt_snode_id-1);
#ifdef _DEBUG_
              logfileptr->OFS()<<AggregatesToRecv(tgt_snode_id-1)<<" aggregates left for Supernode "<<tgt_snode_id<<endl;
#endif
            }
          }
          TIMER_STOP(UPDATE_ANCESTORS);
#else
          //Update everything I own with that factor
          //Update the ancestors
          Int tgt_snode_id = 0;
          Int src_first_row = 0;
          Int src_last_row = 0;
          Int src_nzblk_idx = 0;


          TIMER_START(UPDATE_ANCESTORS);
          while(FindNextUpdate(dist_src_snode, src_nzblk_idx, src_first_row, src_last_row, tgt_snode_id)){ 
            Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);
            if(iTarget == iam){
#ifdef _DEBUG_
              logfileptr->OFS()<<"RECV Supernode "<<tgt_snode_id<<" is updated by Supernode "<<dist_src_snode.Id()<<" rows "<<src_first_row<<" to "<<src_last_row<<" "<<src_nzblk_idx<<std::endl;
#endif

              Int iLocalJ = (tgt_snode_id-1) / np +1 ;
              SuperNode<T> & tgt_snode = *LocalSupernodes_[iLocalJ -1];


            gdb_lock();
              AggregateSuperNode(dist_src_snode,tgt_snode,src_nzblk_idx, src_first_row);

              --AggregatesToRecv(tgt_snode_id-1);
#ifdef _DEBUG_
              logfileptr->OFS()<<AggregatesToRecv(tgt_snode_id-1)<<" aggregates left for Supernode "<<tgt_snode_id<<endl;
#endif

            }
          }
          TIMER_STOP(UPDATE_ANCESTORS);
#endif
        }
        //clear the buffer
        { vector<char>().swap(src_blocks);  }
        { vector<T>().swap(src_nzval);  }







#ifdef _DEBUG_
        assert(AggregatesToRecv(src_snode.Id()-1)==0);
        logfileptr->OFS()<<"  Factoring Supernode "<<I<<std::endl;
#endif
        logfileptr->OFS()<<"  Factoring Supernode "<<I<<std::endl;

        TIMER_START(FACTOR_PANEL);
        //Factorize Diagonal block
        NZBlockDesc & diag_desc = src_snode.GetNZBlockDesc(0);
        for(Int col = 0; col<src_snode.Size();col+=BLOCKSIZE){
          Int bw = min(BLOCKSIZE,src_snode.Size()-col);
          T * diag_nzval = &src_snode.GetNZval(diag_desc.Offset)[col+col*src_snode.Size()];
          lapack::Potrf( 'U', bw, diag_nzval, src_snode.Size());
          T * nzblk_nzval = &src_snode.GetNZval(diag_desc.Offset)[col+(col+bw)*src_snode.Size()];
          blas::Trsm('L','U','T','N',bw, src_snode.NRowsBelowBlock(0)-(col+bw), ONE<T>(),  diag_nzval, src_snode.Size(), nzblk_nzval, src_snode.Size());

          //update the rest !!! (next blocks columns)
          T * tgt_nzval = &src_snode.GetNZval(diag_desc.Offset)[col+bw+(col+bw)*src_snode.Size()];
          blas::Gemm('T','N',src_snode.Size()-(col+bw), src_snode.NRowsBelowBlock(0)-(col+bw),bw,MINUS_ONE<T>(),nzblk_nzval,src_snode.Size(),nzblk_nzval,src_snode.Size(),ONE<T>(),tgt_nzval,src_snode.Size());
        }

#ifdef _DEBUG_
        //        logfileptr->OFS()<<src_snode<<std::endl;
#endif

        TIMER_STOP(FACTOR_PANEL);

        //Send my factor to my ancestors. 
#ifdef _DEBUG_
        IntNumVec is_factor_sent(np);
        SetValue(is_factor_sent,0);
#else
        BolNumVec is_factor_sent(np);
        SetValue(is_factor_sent,false);
#endif

        BolNumVec is_skipped(np);
        SetValue(is_skipped,false);

        Int tgt_snode_id = 0;
        Int src_nzblk_idx = 0;
        Int src_first_row = 0;
        Int src_last_row = 0;


#ifdef UPDATE_LIST
        TIMER_START(FIND_UPDATED_ANCESTORS);
        FindUpdates(src_snode,updates);
        //now traverse the list
        for(std::list<SnodeUpdate>::iterator it = updates.begin(); it!=updates.end();it++){
          Int tgt_snode_id = it->tgt_snode_id;
          Int src_first_row = it->src_fr;
          Int src_nzblk_idx = src_snode.FindBlockIdx(src_first_row);
          Int iTarget = Mapping_.Map(tgt_snode_id-1,src_snode.Id()-1);

          if(iTarget != iam){

#ifdef _DEBUG_
            logfileptr->OFS()<<"Remote Supernode "<<tgt_snode_id<<" is updated by Supernode "<<I<<" rows "<<src_first_row<<std::endl;

            if(is_factor_sent[iTarget]){
              logfileptr->OFS()<<"Remote Supernode "<<tgt_snode_id<<" is updated Supernode "<<I<<" in factor "<<is_factor_sent[iTarget]<<std::endl;
            }
#endif

            if(!is_factor_sent[iTarget] && !is_skipped[iTarget] ){

#ifdef DELAY_SNODES
              //need a std::unordered_set to check whether 
              if(iLocalI < LocalSupernodes_.size()){
                if(LocalSupernodes_[iLocalI]->Id()< tgt_snode_id){
                  FactorsToSend.push_back(DelayedComm(src_snode.Id(),tgt_snode_id,src_nzblk_idx,src_first_row));
#ifdef _DEBUG_
                  logfileptr->OFS()<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<tgt_snode_id<<" from row "<<src_first_row<<endl;
                  cout<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<tgt_snode_id<<" from row "<<src_first_row<<endl;
#endif
                  is_skipped[iTarget] = true;
                  continue;
                }
              }
#endif

#ifdef _DEBUG_
              is_factor_sent[iTarget] = I;
#else
              is_factor_sent[iTarget] = true;
#endif

              Int tgt_first_col = Xsuper_(tgt_snode_id-1);
              Int tgt_last_col = Xsuper_(tgt_snode_id)-1;

              //Send
              NZBlockDesc & pivot_desc = src_snode.GetNZBlockDesc(src_nzblk_idx);

              Int local_first_row = src_first_row-pivot_desc.GIndex;
              Int nzblk_cnt = src_snode.NZBlockCnt()-src_nzblk_idx;
              Int nz_cnt = (src_snode.NRowsBelowBlock(src_nzblk_idx) - local_first_row )*src_snode.Size();
              assert(nz_cnt>0);

              T * nzval_ptr = src_snode.GetNZval(pivot_desc.Offset
                  +local_first_row*src_snode.Size());
              TIMER_START(SEND_MALLOC);
              AddOutgoingComm(outgoingSend, src_snode.Id(), src_snode.Size(), src_first_row, pivot_desc, nzblk_cnt, nzval_ptr, nz_cnt);
              TIMER_STOP(SEND_MALLOC);



              if(outgoingSend.size() > maxIsend_){
                TIMER_START(SEND_MPI);
                MPI_Send(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,pComm);
                TIMER_STOP(SEND_MPI);

                outgoingSend.pop_back();
              }
              else{
                TIMER_START(SEND_MPI);
                MPI_Isend(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,pComm,&outgoingSend.back()->Request);
                TIMER_STOP(SEND_MPI);
              }



#ifdef _DEBUG_            
              logfileptr->OFS()<<"Sending "<<nz_cnt*sizeof(T)<<" bytes to P"<<iTarget<<std::endl;
              logfileptr->OFS()<<"     Sent factor "<<src_snode.Id()<<" to Supernode "<<tgt_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
              logfileptr->OFS()<<"Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
              logfileptr->OFS()<<"Sending "<<src_blocks.size()<<" bytes to P"<<iTarget<<std::endl;
#endif

            }
          }
          else{
            Int iLocalJ = (tgt_snode_id-1) / np +1 ;

#ifdef _DEBUG_
            assert(LocalSupernodes_[iLocalJ-1]->Id()==tgt_snode_id);
#endif
            LocalUpdates[iLocalJ-1].push(LocalUpdate(src_snode.Id(),src_nzblk_idx,src_first_row));
          }
        }
        TIMER_STOP(FIND_UPDATED_ANCESTORS);
#else
        TIMER_START(FIND_UPDATED_ANCESTORS);
        while(FindNextUpdate(src_snode, src_nzblk_idx, src_first_row, src_last_row, tgt_snode_id)){ 
          Int iTarget = Mapping_.Map(tgt_snode_id-1,src_snode.Id()-1);

          if(iTarget != iam){


#ifdef _DEBUG_
            logfileptr->OFS()<<"Remote Supernode "<<tgt_snode_id<<" is updated by Supernode "<<I<<" rows "<<src_first_row<<" to "<<src_last_row<<" "<<src_nzblk_idx<<std::endl;

            if(is_factor_sent[iTarget]){
              logfileptr->OFS()<<"Remote Supernode "<<tgt_snode_id<<" is updated Supernode "<<I<<" in factor "<<is_factor_sent[iTarget]<<std::endl;
            }
#endif

            if(!is_factor_sent[iTarget] && !is_skipped[iTarget] ){

#ifdef DELAY_SNODES
              //need a std::unordered_set to check whether 
              if(iLocalI < LocalSupernodes_.size()){
                if(LocalSupernodes_[iLocalI]->Id()< tgt_snode_id){
                  //need to push the prev src_last_row
                  FactorsToSend.push_back(DelayedComm(src_snode.Id(),tgt_snode_id,src_nzblk_idx,src_first_row));
#ifdef _DEBUG_
                  logfileptr->OFS()<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<tgt_snode_id<<" from row "<<src_first_row<<endl;
                  cout<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<tgt_snode_id<<" from row "<<src_first_row<<endl;
#endif
                  is_skipped[iTarget] = true;
                  continue;
                }
              }
#endif



#ifdef _DEBUG_
              is_factor_sent[iTarget] = I;
#else
              is_factor_sent[iTarget] = true;
#endif

              Int tgt_first_col = Xsuper_(tgt_snode_id-1);
              Int tgt_last_col = Xsuper_(tgt_snode_id)-1;

              //Send
              NZBlockDesc & pivot_desc = src_snode.GetNZBlockDesc(src_nzblk_idx);

              Int local_first_row = src_first_row-pivot_desc.GIndex;



              //              logfileptr->OFS()<<src_first_row<<std::endl;
              //              logfileptr->OFS()<<local_first_row<<std::endl;
              //              logfileptr->OFS()<<pivot_desc.GIndex<<std::endl;
              //              assert(src_first_row < pivot_desc.GIndex + src_snode.NRows(src_nzblk_idx));

              Int nzblk_cnt = src_snode.NZBlockCnt()-src_nzblk_idx;
              Int nz_cnt = (src_snode.NRowsBelowBlock(src_nzblk_idx) - local_first_row )*src_snode.Size();
#ifdef _DEBUG_
              assert(nz_cnt>0);
#endif
              T * nzval_ptr = src_snode.GetNZval(pivot_desc.Offset
                  +local_first_row*src_snode.Size());


              TIMER_START(SEND_MALLOC);
              AddOutgoingComm(outgoingSend, src_snode.Id(), src_snode.Size(), src_first_row, pivot_desc, nzblk_cnt, nzval_ptr, nz_cnt);
              TIMER_STOP(SEND_MALLOC);

//gdb_lock();
              if(outgoingSend.size() > maxIsend_){
                TIMER_START(SEND_MPI);
                MPI_Send(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,pComm);
                TIMER_STOP(SEND_MPI);

                outgoingSend.pop_back();
              }
              else{
                TIMER_START(SEND_MPI);
                MPI_Isend(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,pComm,&outgoingSend.back()->Request);
                TIMER_STOP(SEND_MPI);
              }


#ifdef _DEBUG_            
              logfileptr->OFS()<<"Sending "<<nz_cnt*sizeof(T)<<" bytes to P"<<iTarget<<std::endl;
              logfileptr->OFS()<<"     Sent factor "<<src_snode.Id()<<" to Supernode "<<tgt_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
              logfileptr->OFS()<<"Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
              logfileptr->OFS()<<"Sending "<<outgoingSend.size()<<" bytes to P"<<iTarget<<std::endl;
#endif

            }
          }
//          else{
//            Int iLocalJ = (tgt_snode_id-1) / np +1 ;
//
//            assert(LocalSupernodes_[iLocalJ-1]->Id()==tgt_snode_id);
//
//            LocalUpdates[iLocalJ-1].push(I);
//            LocalUpdates[iLocalJ-1].push(src_nzblk_idx);
//            LocalUpdates[iLocalJ-1].push(src_first_row);
//          }

        }
        TIMER_STOP(FIND_UPDATED_ANCESTORS);
#endif
      }



      Int src_snode_id = I;
      //TODO replace this by a PLINDX ? 
      Int tgt_snode_id = FBUpdate(src_snode_id); 

      //If I am involved in updating tgt_snode_id
      if(tgt_snode_id!= -1){
        bool skipped = false;
#ifdef DELAY_SNODES
        if(LocalSupernodes_.size()>0){
          Int nextLocalI = max(1,iLocalI);
          Int next_snode_id = LocalSupernodes_[nextLocalI-1]->Id();
          if(next_snode_id<=src_snode_id){
            nextLocalI = min((Int)LocalSupernodes_.size()-1,nextLocalI+1);
            next_snode_id = LocalSupernodes_[nextLocalI-1]->Id();
          }

          if(next_snode_id < tgt_snode_id){
            skipped = true;

            //push this thing to be done later
            UpdatesToRecv.push_back(DelayedComm(src_snode_id,tgt_snode_id,-1,-1));
#ifdef _DEBUG_
                  logfileptr->OFS()<<"P"<<iam<<" has delayed reception of factor from Supernode "<<src_snode_id<<" to "<<tgt_snode_id<<endl;
#endif

          }
        }
#endif

        if(!skipped){
          std::vector<char> src_blocks;
          SuperNode<T> * cur_src_snode = FBRecvFactor(src_snode_id,tgt_snode_id,src_blocks);

          Int iSrcOwner = Mapping_.Map(src_snode_id-1,src_snode_id-1);

          //Compute update to tgt_snode_id and put it in my aggregate vector 

          //Update everything src_snode_id own with that factor
          //Update the ancestors
          tgt_snode_id = 0;
          Int src_first_row = 0;
          Int src_last_row = 0;
          Int src_nzblk_idx = 0;

          TIMER_START(UPDATE_ANCESTORS);
          while(FindNextUpdate(*cur_src_snode, src_nzblk_idx, src_first_row, src_last_row, tgt_snode_id)){

            Int iTarget = Mapping_.Map(tgt_snode_id-1,cur_src_snode->Id()-1);
            if(iTarget == iam){
              SuperNode<T> * tgt_aggreg;

              Int iTgtOwner = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);
              if(iTgtOwner == iam){
                //the aggregate vector is directly the target snode
                Int iLocalJ = (tgt_snode_id-1) / np +1 ;
                tgt_aggreg = LocalSupernodes_[iLocalJ -1];

                assert(tgt_snode_id == tgt_aggreg->Id());

              }
              else{
                //Check if src_snode_id already have an aggregate vector
                if(aggVectors[tgt_snode_id-1]==NULL){
                  aggVectors[tgt_snode_id-1] = new SuperNode<T>(tgt_snode_id, Xsuper_[tgt_snode_id-1], Xsuper_[tgt_snode_id]-1, Size(), xlindx_, lindx_);
                }
                tgt_aggreg = aggVectors[tgt_snode_id-1];
              }

#ifdef _DEBUG_
              logfileptr->OFS()<<"RECV Supernode "<<tgt_snode_id<<" is updated by Supernode "<<cur_src_snode->Id()<<" rows "<<src_first_row<<" to "<<src_last_row<<" "<<src_nzblk_idx<<std::endl;
#endif




#ifdef SINGLE_BLAS
              UpdateSuperNode(*cur_src_snode,*tgt_aggreg,src_nzblk_idx, tmpBuf,src_first_row);
#else
              UpdateSuperNode(*cur_src_snode,*tgt_aggreg,src_nzblk_idx, src_first_row);
#endif


              --UpdatesToDo(tgt_snode_id-1);
#ifdef _DEBUG_
              logfileptr->OFS()<<UpdatesToDo(tgt_snode_id-1)<<" updates left for Supernode "<<tgt_snode_id<<endl;
#endif


              if(iTgtOwner==iam){
                --AggregatesToRecv(tgt_snode_id-1);
#ifdef _DEBUG_
                logfileptr->OFS()<<AggregatesToRecv(tgt_snode_id-1)<<" aggregates left for Supernode "<<tgt_snode_id<<endl;
#endif
              }









              //If this is my last update sent it to tgt_snode_id
              if(src_snode_id == LastUpdate[tgt_snode_id-1]){
                if(iTgtOwner != iam){

                  gdb_lock();
#ifdef _DEBUG_
                  logfileptr->OFS()<<"Remote Supernode "<<tgt_snode_id<<" is updated by Supernode "<<src_snode_id<<std::endl;
#endif


                  //need a std::unordered_set to check whether 
                  //              if(iLocalI < LocalSupernodes_.size()){
                  //                if(LocalSupernodes_[iLocalI]->Id()< tgt_snode_id){
                  //                  //need to push the prev src_last_row
                  //                  AggregatesToSend.push_back(DelayedComm(src_snode_id,tgt_snode_id,src_nzblk_idx,src_first_row));
                  //#ifdef _DEBUG_
                  //                  logfileptr->OFS()<<"P"<<iam<<" has delayed update from Supernode "<<src_snode_id<<" to "<<tgt_snode_id<<" from row "<<src_first_row<<endl;
                  //                  cout<<"P"<<iam<<" has delayed update from Supernode "<<src_snode_id<<" to "<<tgt_snode_id<<" from row "<<src_first_row<<endl;
                  //#endif
                  //                  is_skipped[iTarget] = true;
                  //                  continue;
                  //                }
                  //              }


                  Int tgt_first_col = Xsuper_(tgt_snode_id-1);
                  Int tgt_last_col = Xsuper_(tgt_snode_id)-1;

                  //Send
                  NZBlockDesc & pivot_desc = tgt_aggreg->GetNZBlockDesc(0);
                  Int local_first_row = 0;//src_first_row-pivot_desc.GIndex;

                  Int nzblk_cnt = tgt_aggreg->NZBlockCnt();//-src_nzblk_idx;
                  Int nz_cnt = (tgt_aggreg->NRowsBelowBlock(0) - local_first_row )*tgt_aggreg->Size();
#ifdef _DEBUG_
                  assert(nz_cnt>0);
#endif
                  T * nzval_ptr = tgt_aggreg->GetNZval(pivot_desc.Offset
                      +local_first_row*tgt_aggreg->Size());


                  TIMER_START(SEND_MALLOC);
                  AddOutgoingComm(outgoingSend, src_snode_id, tgt_aggreg->Size(), tgt_aggreg->FirstCol(), pivot_desc, nzblk_cnt, nzval_ptr, nz_cnt);
                  TIMER_STOP(SEND_MALLOC);

                  //gdb_lock();
                  //              if(outgoingSend.size() > maxIsend_){
                  TIMER_START(SEND_MPI);
                  MPI_Send(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTgtOwner,tgt_snode_id,pComm);
                  TIMER_STOP(SEND_MPI);

                  outgoingSend.pop_back();
                  //              }
                  //              else{
                  //                TIMER_START(SEND_MPI);
                  //                MPI_Isend(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,pComm,&outgoingSend.back()->Request);
                  //                TIMER_STOP(SEND_MPI);
                  //              }


#ifdef _DEBUG_            
                  logfileptr->OFS()<<"Sending "<<nz_cnt*sizeof(T)<<" bytes to P"<<iTgtOwner<<std::endl;
                  logfileptr->OFS()<<"     Sent aggregate "<<tgt_snode_id<<" to Supernode "<<tgt_snode_id<<" on P"<<iTgtOwner<<std::endl;
                  logfileptr->OFS()<<"Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
                  logfileptr->OFS()<<"Sending "<<outgoingSend.size()<<" bytes to P"<<iTgtOwner<<std::endl;
#endif

                }








              }









            }
          }
          TIMER_STOP(UPDATE_ANCESTORS);

          if(iSrcOwner!=iam){
            delete cur_src_snode;
            //clear the buffer
            { vector<char>().swap(src_blocks);  }
          }

        }

      }


    }
    I++;


  MPI_Barrier(pComm);

  }



  MPI_Barrier(pComm);

  TIMER_STOP(FACTORIZATION_FO);
}

template<typename T> Int SupernodalMatrix<T>::FBUpdate(Int I){
      //Check if I have anything to update with that supernode
      //look at lindx_
      Int fi = xlindx_(I-1);
      Int li = xlindx_(I)-1;
       

 
      Int J = -1; 
      for(Int idx = fi; idx<=li;++idx){
        Int row = lindx_[idx-1];
        J = SupMembership_[row-1];
        Int iUpdater = Mapping_.Map(J-1,I-1);
        if(iUpdater == iam && J>I){
          break;
        }
      }
#ifdef _DEBUG_
assert(J>=1);
#endif

      return J;
      }


  template <typename T> void SupernodalMatrix<T>::FBGetUpdateCount(IntNumVec & sc,IntNumVec & lu){
    sc.Resize(Xsuper_.m());
    SetValue(sc,I_ZERO);

    lu.Resize(Xsuper_.m());
    SetValue(lu,I_ZERO);

    IntNumVec marker(Xsuper_.m());
    SetValue(marker,I_ZERO);

    for(Int s = 1; s<Xsuper_.m(); ++s){
      Int first_col = Xsuper_(s-1);
      Int last_col = Xsuper_(s)-1;

      Int fi = xlindx_(s-1);
      Int li = xlindx_(s)-1;

#ifndef _DEBUG_
  #define nodebugtmp
  #define _DEBUG_
#endif



#ifdef nodebugtmp
  #undef _DEBUG_
#endif



#ifdef _DEBUG_
      logfileptr->OFS()<<"Supernode "<<s<<" updates: ";
#endif

      for(Int row_idx = fi; row_idx<=li;++row_idx){
        Int row = lindx_(row_idx-1);
        Int supno = SupMembership_(row-1);

        if(marker(supno-1)!=s && supno!=s){

#ifdef _DEBUG_
          logfileptr->OFS()<<supno<<" ";
#endif


          Int iUpdater = Mapping_.Map(supno-1,s-1);

          if(iam == iUpdater){
            ++sc[supno-1];
            lu[supno-1] = s;
          }

          marker(supno-1) = s;
        }
      }

#ifdef _DEBUG_
      logfileptr->OFS()<<std::endl;
#endif

#ifdef nodebugtmp
  #undef _DEBUG_
#endif
    }
  }


template<typename T> SuperNode<T> * SupernodalMatrix<T>::FBRecvFactor(Int src_snode_id,Int tgt_snode_id, std::vector<char> & src_blocks){
  SuperNode<T> * cur_src_snode = NULL;
  Int iSrcOwner = Mapping_.Map(src_snode_id-1,src_snode_id-1);

  if(iSrcOwner==iam){
    Int iLocalI = (src_snode_id-1) / np +1 ;
    cur_src_snode = LocalSupernodes_[iLocalI -1];
  }
  else{
    Int nz_cnt;

    TIMER_START(RECV_MPI);
    MPI_Status recv_status;
    int bytes_received = 0;

    //gdb_lock();
    MPI_Probe(iSrcOwner,tgt_snode_id,pComm,&recv_status);
    MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);

    TIMER_START(RECV_MALLOC);
    src_blocks.resize(bytes_received);
    TIMER_STOP(RECV_MALLOC);


    MPI_Recv(&src_blocks[0],bytes_received,MPI_BYTE,recv_status.MPI_SOURCE,recv_status.MPI_TAG,pComm,&recv_status);
    MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);


    //unpack the data
    Int recv_snode_id = *(Int*)&src_blocks[0];
#ifdef _DEBUG_
    assert(src_snode_id == recv_snode_id);
#endif
    Int src_nzblk_cnt = *(((Int*)&src_blocks[0])+1);
    NZBlockDesc * src_blocks_ptr = 
      reinterpret_cast<NZBlockDesc*>(&src_blocks[2*sizeof(Int)]);
    Int src_nzval_cnt = *(Int*)(src_blocks_ptr + src_nzblk_cnt);
    T * src_nzval_ptr = (T*)((Int*)(src_blocks_ptr + src_nzblk_cnt)+1);
    TIMER_STOP(RECV_MPI);

    //Create the dummy supernode for that data
    SuperNode<T>* dist_src_snode = new SuperNode<T>(src_snode_id,Xsuper_[src_snode_id-1],Xsuper_[src_snode_id]-1, Size(), src_blocks_ptr, src_nzblk_cnt, src_nzval_ptr, src_nzval_cnt);



#ifdef _DEBUG_
    logfileptr->OFS()<<"RECV Supernode "<<dist_src_snode->Id()<<std::endl;
#endif
    cur_src_snode = dist_src_snode;

  }
  return cur_src_snode;
}



#endif
