#ifndef _SUPERNODAL_MATRIX_IMPL_HPP_
#define _SUPERNODAL_MATRIX_IMPL_HPP_

#include "SupernodalMatrix.hpp"



namespace LIBCHOLESKY{

template <typename T> SupernodalMatrix<T>::SupernodalMatrix(){
}

template <typename T> SupernodalMatrix<T>::SupernodalMatrix(const DistSparseMatrix<T> & pMat, FBMatrix* Afactptr ){
  Int iam = Afactptr->iam;
  Int np = Afactptr->np;


  pMat.ConstructETree(ETree_);

  IntNumVec cc,rc;
  pMat.GetLColRowCount(ETree_,cc,rc);

  pMat.FindSupernodes(ETree_,cc,Xsuper_);

  IntNumVec xlindx,lindx;
  IntNumVec xlnz;
  DblNumVec lnz;

  //TODO might need to be moved to SparseMatrixStructure ?
  pMat.SymbolicFactorization(ETree_,cc,Xsuper_,xlindx,xlnz,lindx,lnz);

  //copy
    for(Int I=1;I<Xsuper_.m();I++){
      Int fc = Xsuper_(I-1);
      Int lc = Xsuper_(I)-1;
      Int iWidth = lc-fc+1;
      Int fi = xlindx(I-1);
      Int li = xlindx(I)-1;

      Int iDest = Afactptr->MAP(I-1,I-1);

      //parse the first column to create the NZBlock
      if(iam==iDest){
        LocalSupernodes_.push_back(SuperNode(I,fc,lc));
        SuperNode & snode = LocalSupernodes_.back();


        for(Int idx = fi; idx<=li;idx++){
          Int iStartRow = lindx(idx-1);
          Int iPrevRow = iStartRow;
          Int iContiguousRows = 1;
          for(Int idx2 = idx+1; idx2<=li;idx2++){
            Int iCurRow = lindx(idx2-1);
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
          logfileptr->OFS()<<"Creating a new NZBlock for rows "<<iStartRow<<" to "<<iStartRow + iContiguousRows-1<<" with "<<iCurNZcnt<<" nz."<<std::endl;
#endif
          snode.AddNZBlock(iContiguousRows , iWidth,iStartRow);
        }

        
        //Add an extra past last row
//        Int iLastRow = lindx(li-1);
//        LocalIndices.push_back(iLastRow+1);

//DEBUG
//        logfileptr->OFS()<<"***********************************"<<std::endl;
//      logfileptr->OFS()<<LocalIndices<<std::endl;
//      for(int blkidx=0;blkidx<Lcol.size();++blkidx){
//        NZBlock<double> & nzblk = *Lcol[blkidx];
////        Int lastRow = nzblk.GIndex() + nzblk.NRows() -1;
//        logfileptr->OFS()<<nzblk<<std::endl;
//      }
//        logfileptr->OFS()<<"***********************************"<<std::endl;
//DEBUG


      }

      //Distribute the data

      //look at the owner of the first column
      Int numColFirst = pMat.size / np;
      Int prevOwner = -1;
      DblNumVec aRemoteCol;
      double * pdNzVal = NULL;
      Int iStartIdxCopy = 0;
      //copy the data from A into this Block structure
      for(Int i = fc;i<=lc;i++){

        Int iOwner = std::min((i-1)/numColFirst,np-1);

        if(iOwner != prevOwner){

          prevOwner = iOwner;
          //we need to transfer the bulk of columns
          //first compute the number of cols we need to transfer
          Int iNzTransfered = 0;
          Int iLCTransfered = lc;

          for(Int j =i;j<=lc;++j){
            iOwner = std::min((j-1)/numColFirst,np-1);
            if(iOwner == prevOwner){
              Int nrows = pMat.Global_.colptr(j) - pMat.Global_.colptr(j-1);
              iNzTransfered+=nrows;
              //              logfileptr->OFS()<<"Looking at col "<<j<<" which is on P"<<iOwner<<std::endl;
            }
            else{
              iLCTransfered = j-1;
              break;
            }
          } 


          //#ifdef _DEBUG_
          //logfileptr->OFS()<<"Column "<<i<<" to "<<iLCTransfered<<" are owned by P"<<prevOwner<<" and should go on P"<<iDest<<std::endl;
          //logfileptr->OFS()<<"They contain "<<iNzTransfered<<" nz"<<std::endl;
          //#endif

          //if data needs to be transfered
          if(iDest!=prevOwner){
            if(iam == iDest){
              aRemoteCol.Resize(iNzTransfered);
              //MPI_Recv
              MPI_Recv(aRemoteCol.Data(),iNzTransfered*sizeof(double),MPI_BYTE,prevOwner,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

              //logfileptr->OFS()<<"Received col "<<i<<" to "<<iLCTransfered<<" from P"<<iOwner<<std::endl;
              iStartIdxCopy = 0;
            } 
            else if (iam == prevOwner){

              Int local_i = (i-(numColFirst)*iam);
              Int iColptrLoc = pMat.Local_.colptr(local_i-1);
              Int iRowIndLoc = pMat.Local_.rowind(iColptrLoc-1);
              Int iLastRowIndLoc = pMat.Local_.rowind(pMat.Local_.colptr(local_i)-1-1);

              //logfileptr->OFS()<<"cols "<<i<<" to "<<iLCTransfered<<": rows "<<iRowIndLoc<<" to "<<iLastRowIndLoc<<std::endl;






              double * pdData = &pMat.nzvalLocal(iColptrLoc-1);

              //MPI_send
              MPI_Send(pdData,iNzTransfered*sizeof(double),MPI_BYTE,iDest,0,MPI_COMM_WORLD);

              //logfileptr->OFS()<<"Sent col "<<i<<" to "<<iLCTransfered<<" to P"<<iDest<<std::endl;
            }
          }
        }

        //copy the data if I own it
        if(iam == iDest){
          SuperNode & snode = LocalSupernodes_.back();

          //isData transfered or local
          if(iam!=prevOwner){
            pdNzVal = aRemoteCol.Data();
            //logfileptr->OFS()<<"pdNzVal is the remote buffer"<<std::endl;
            //logfileptr->OFS()<<aRemoteCol<<std::endl;
          }
          else{
            Int local_i = (i-(numColFirst)*iam);
            Int iColptrLoc = pMat.Local_.colptr(local_i-1);
            pdNzVal = &pMat.nzvalLocal(iColptrLoc-1);
            //logfileptr->OFS()<<"pdNzVal is the local pMat"<<std::endl;
            iStartIdxCopy = 0;
          }

          //Copy the data from pdNzVal in the appropriate NZBlock

          Int iGcolptr = pMat.Global_.colptr(i-1);
          Int iRowind = pMat.Global_.rowind(iGcolptr-1);
          Int iNrows = pMat.Global_.colptr(i) - pMat.Global_.colptr(i-1);
          Int idxA = 0;
          Int iLRow = 0;
          Int firstrow = fi + i-fc;
          for(Int idx = firstrow; idx<=li;idx++){
            iLRow = lindx(idx-1);
            //logfileptr->OFS()<<"Looking at L("<<iLRow<<","<<i<<")"<<std::endl;
            if( iLRow == iRowind){
              Int iNZBlockIdx = snode.FindBlockIdx(iLRow);

              double elem = pdNzVal[iStartIdxCopy + idxA];

              NZBlock<double> & pDestBlock = snode.GetNZBlock(iNZBlockIdx);
              //logfileptr->OFS()<<*pDestBlock<<std::endl;
              
              //find where we should put it
              Int localCol = i - fc;
              Int localRow = iLRow - pDestBlock.GIndex();
#ifdef _DEBUG_
              logfileptr->OFS()<<"Elem is A("<<iRowind<<","<<i<<") = "<<elem<<" at ("<<localRow<<","<<localCol<<") in "<<iNZBlockIdx<<"th NZBlock of L"<< std::endl;
#endif
              pDestBlock.Nzval(localRow,localCol) = elem;
              if(idxA<iNrows){
                idxA++;
                iRowind = pMat.Global_.rowind(iGcolptr+idxA-1);
              }
            }
          }

        }
        
      }
#ifdef _DEBUG_
      logfileptr->OFS()<<"--------------------------------------------------"<<std::endl;
#endif
    }

}

} // namespace LIBCHOLESKY


#endif 
