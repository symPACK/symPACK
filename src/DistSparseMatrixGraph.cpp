#include "sympack/DistSparseMatrixGraph.hpp"
#include "sympack/SparseMatrixStructure.hpp"
#include "sympack/ETree.hpp"
#include "sympack/utility.hpp"
#include <limits>       // std::numeric_limits

#include <iterator>
#include <set>
#include <list>
#include <vector>
#include <algorithm>
#include <numeric>

#ifdef _USE_COREDUMPER_
#include <google/coredumper.h>
#endif

#define USE_REDUCE

namespace SYMPACK{


  SparseMatrixGraph::SparseMatrixGraph(){
    baseval = 1;
    keepDiag = 1;
    sorted = 1;
    nnz = 0;
    size = 0;
  }


  void SparseMatrixGraph::SetBaseval(int aBaseval){
    if(aBaseval!= baseval){
      for(auto it = colptr.begin(); it!=colptr.end();it++){ *it = *it - baseval + aBaseval;}
      for(auto it = rowind.begin(); it!=rowind.end();it++){ *it = *it - baseval + aBaseval;}
    }
    baseval = aBaseval; 
  }

  void SparseMatrixGraph::SetKeepDiag(int aKeepDiag){
    if(aKeepDiag != keepDiag){
      if(aKeepDiag){
        //add the diagonal entry
        SYMPACK::vector<Idx> newRowind(rowind.size()+VertexCount());
        Ptr pos = 0;
        for(Idx locCol = 0;locCol < VertexCount(); locCol++){
          Idx col = baseval + locCol; // baseval based
          Ptr colbeg = colptr[locCol] - baseval;
          Ptr colend = colptr[locCol+1] - baseval;
          colptr[locCol] = pos + baseval; //baseval based
          bool diagAdded = false;
          for(Ptr eptr = colbeg; eptr < colend; eptr++){
            Idx row = rowind[eptr]; //baseval based
            if(row<col){
              newRowind[pos++] = row; 
            }
            if(row>col){
              if(!diagAdded){
                newRowind[pos++] = col;
                diagAdded = true;
              }
              newRowind[pos++] = row; 
            }
          }

          if(!diagAdded){
            newRowind[pos++] = col;
            diagAdded = true;
          }
        } 
        if(colptr.size()>0){
          colptr.back() = newRowind.size()+baseval;
        }
        rowind.swap(newRowind);
      }
      else{
        //remove the diagonal entry
        Ptr pos = 0;
        for(Idx locCol = 0;locCol < VertexCount(); locCol++){
          Idx col = baseval + locCol; // baseval based
          Ptr colbeg = colptr[locCol] - baseval;
          Ptr colend = colptr[locCol+1] - baseval;
          colptr[locCol] = pos + baseval; //baseval based
          for(Ptr eptr = colbeg; eptr < colend; eptr++){
            Idx row = rowind[eptr]; //baseval based
            if(row!=col){
              rowind[pos++] = row; 
            }
          }
        } 
        rowind.resize(pos);
        if(colptr.size()>0){
          colptr.back() = rowind.size()+baseval;
        }
      }
    }
    keepDiag = aKeepDiag; 
  }

  void SparseMatrixGraph::SetSorted(int aSorted){
    if(aSorted!=sorted){
      if(aSorted){
        SortEdges();
      }
      sorted = aSorted;
    }
  }


  void SparseMatrixGraph::SortEdges(){
    for(Idx locCol = 0 ; locCol< VertexCount(); locCol++){
      Ptr colbeg = colptr[locCol]-baseval; //now 0 based
      Ptr colend = colptr[locCol+1]-baseval; // now 0 based 
      sort(&rowind[0]+colbeg,&rowind[0]+colend,std::less<Ptr>());
    }
  }


  void SparseMatrixGraph::DistributeGraph(DistSparseMatrixGraph & dg){

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
      pinout[i] = max(pinout[i], pin[i]);
    } 
  } 






  DistSparseMatrixGraph::DistSparseMatrixGraph(){
    bIsExpanded = false;
    baseval = 1;
    keepDiag = 1;
    sorted = 1;
    nnz = 0;
    size = 0;
    comm = MPI_COMM_NULL;
    mpirank = 0;
    mpisize = 1;
    vertexDist.resize(mpisize+1,baseval);
  }

  DistSparseMatrixGraph::~DistSparseMatrixGraph(){
    if(comm != MPI_COMM_NULL){
      MPI_Comm_free(&comm);
    }
  }

  DistSparseMatrixGraph& DistSparseMatrixGraph::operator=( const DistSparseMatrixGraph& g ) {
    bIsExpanded = g.bIsExpanded;
    size = g.size;
    nnz = g.nnz;
    colptr = g.colptr;
    rowind = g.rowind;
    vertexDist = g.vertexDist;
    SetComm(g.comm);
    baseval = g.baseval;
    keepDiag = g.keepDiag;
    sorted = g.sorted;
    return *this;
  }

  DistSparseMatrixGraph::DistSparseMatrixGraph( const DistSparseMatrixGraph& g ):DistSparseMatrixGraph(){
    (*this) = g;
  }

  DistSparseMatrixGraph::DistSparseMatrixGraph(const SparseMatrixStructure & A):DistSparseMatrixGraph(){
    baseval = A.baseval;
    keepDiag = A.keepDiag;
    sorted = A.sorted;
    FromStructure(A);
  }



  void DistSparseMatrixGraph::SetComm(const MPI_Comm & aComm){
    if(aComm!=comm){
      if(comm!= MPI_COMM_NULL){
        MPI_Comm_free(&comm);
        //what do we do with vertex dist ?
      }
      MPI_Comm_dup(aComm,&comm);
      MPI_Comm_size(comm,&mpisize);
      MPI_Comm_rank(comm,&mpirank);

      //recompute a "balanced" distribution of vertices
      Int colPerProc = size / mpisize;
      vertexDist.assign(mpisize+1,colPerProc);
      vertexDist[0] = baseval;
      std::partial_sum(vertexDist.begin(),vertexDist.end(),vertexDist.begin());
      vertexDist.back() = size + baseval;
    }
  }

  void DistSparseMatrixGraph::SetBaseval(int aBaseval){
    if(aBaseval!= baseval){
      for(auto it = colptr.begin(); it!=colptr.end();it++){ *it = *it - baseval + aBaseval;}
      for(auto it = rowind.begin(); it!=rowind.end();it++){ *it = *it - baseval + aBaseval;}
      for(auto it = vertexDist.begin(); it!=vertexDist.end();it++){ *it = *it - baseval + aBaseval;}
    }
    baseval = aBaseval; 
  }

  void DistSparseMatrixGraph::SetKeepDiag(int aKeepDiag){
    if(aKeepDiag != keepDiag){
      if(aKeepDiag){
        //add the diagonal entry
        SYMPACK::vector<Idx> newRowind(rowind.size()+LocalVertexCount());
        Ptr pos = 0;
        for(Idx locCol = 0;locCol < LocalVertexCount(); locCol++){
          Idx col = LocalFirstVertex() + locCol; // baseval based
          Ptr colbeg = colptr[locCol] - baseval;
          Ptr colend = colptr[locCol+1] - baseval;
          colptr[locCol] = pos + baseval; //baseval based
          bool diagAdded = false;
          for(Ptr eptr = colbeg; eptr < colend; eptr++){
            Idx row = rowind[eptr]; //baseval based
            if(row<col){
              newRowind[pos++] = row; 
            }
            if(row>col){
              if(!diagAdded){
                newRowind[pos++] = col;
                diagAdded = true;
              }
              newRowind[pos++] = row; 
            }
          }

          if(!diagAdded){
            newRowind[pos++] = col;
            diagAdded = true;
          }
        } 
        if(colptr.size()>0){
          colptr.back() = newRowind.size()+baseval;
        }
        rowind.swap(newRowind);
      }
      else{
        //remove the diagonal entry
        Ptr pos = 0;
        for(Idx locCol = 0;locCol < LocalVertexCount(); locCol++){
          Idx col = LocalFirstVertex() + locCol; // baseval based
          Ptr colbeg = colptr[locCol] - baseval;
          Ptr colend = colptr[locCol+1] - baseval;
          colptr[locCol] = pos + baseval; //baseval based
          for(Ptr eptr = colbeg; eptr < colend; eptr++){
            Idx row = rowind[eptr]; //baseval based
            if(row!=col){
              rowind[pos++] = row; 
            }
          }
        } 
        rowind.resize(pos);
        if(colptr.size()>0){
          colptr.back() = rowind.size()+baseval;
        }
      }
    }
    keepDiag = aKeepDiag; 
  }

  void DistSparseMatrixGraph::SetSorted(int aSorted){
    if(aSorted!=sorted){
      if(aSorted){
        SortEdges();
      }
      sorted = aSorted;
    }
  }



  void DistSparseMatrixGraph::FromStructure(const SparseMatrixStructure & A){

    if(A.bIsGlobal){
      throw std::logic_error( "DistSparseMatrixGraph from Global (non distributed) SparseMatrixStructure not implemented yet\n" );
    }
    else{
      bIsExpanded = A.bIsExpanded;

      size = A.size;
      Idx locColCnt = A.IsExpanded()?A.expColptr.size()-1:A.colptr.size()-1;
      Idx firstCol = 0;



      Int colPerProc = size / mpisize;
      vertexDist.assign(mpisize+1,colPerProc);
      vertexDist[0] = baseval;
      std::partial_sum(vertexDist.begin(),vertexDist.end(),vertexDist.begin());
      vertexDist.back() = size + baseval;

      firstCol = LocalFirstVertex()-baseval; //0 based

      colptr.resize(locColCnt+1);

      const Ptr * acolptr = A.IsExpanded()?&A.expColptr[0]:&A.colptr[0];
      const Idx * arowind = A.IsExpanded()?&A.expRowind[0]:&A.rowind[0];

      nnz = keepDiag?(A.keepDiag?A.nnz:A.nnz+locColCnt):(A.keepDiag?A.nnz-locColCnt:A.nnz);
      rowind.resize(nnz);

      colptr[0] = baseval;
      for(Idx locCol = 0; locCol < locColCnt; locCol++){
        Idx col = firstCol + locCol; // 0 -based
        Ptr colbeg = acolptr[locCol] - A.baseval;
        Ptr colend = acolptr[locCol+1] - A.baseval;

        Ptr & pos = colptr[locCol+1];
        pos = colptr[locCol];

        //if we need to add the diagonal entry
        if(!A.keepDiag && keepDiag){
          rowind[pos++ - baseval] = col + baseval;
        }
        for(Ptr jptr = colbeg; jptr<colend; jptr++){
          Idx row = arowind[jptr] - A.baseval;
          if(keepDiag || row!=col){
            rowind[pos++ - baseval] = row + baseval;
          }
        }
      }

      if(sorted){
        SortEdges();
      }
    }
  }


  void DistSparseMatrixGraph::SortEdges(){
    for(Idx locCol = 0 ; locCol< LocalVertexCount(); locCol++){
      Ptr colbeg = colptr[locCol]-baseval; //now 0 based
      Ptr colend = colptr[locCol+1]-baseval; // now 0 based 
      sort(&rowind[0]+colbeg,&rowind[0]+colend,std::less<Ptr>());
    }
  }

  void DistSparseMatrixGraph::Permute(Int * invp){
    permute_(invp, NULL, 1);
  }
  void DistSparseMatrixGraph::Permute(Int * invp, Idx * newVertexDist){
    permute_(invp, newVertexDist, 1);
  }
  void DistSparseMatrixGraph::Permute(Int * invp, Int invpbaseval){
    permute_(invp, &vertexDist[0], invpbaseval);
  }
  void DistSparseMatrixGraph::Permute(Int * invp, Idx * newVertexDist, Int invpbaseval){
    permute_(invp, newVertexDist, invpbaseval);
  }

  void DistSparseMatrixGraph::permute_(Int * invp, Idx * newVertexDist, Int invpbaseval){
    TIMER_START(PERMUTE);

    //handle default parameter values
    if(newVertexDist==NULL){
      newVertexDist = &vertexDist[0];
    }

    //Idx colPerProc = size / mpisize;
    Int firstCol = LocalFirstVertex()-baseval; //0 based
    Int newFirstCol = newVertexDist[mpirank]-baseval; //0 based
    Ptr newVtxCount = newVertexDist[mpirank+1] - newVertexDist[mpirank];


    SYMPACK::vector<int> sizes(mpisize,0);

    for(Idx locCol = 0; locCol<LocalVertexCount(); locCol++){
      Idx col = firstCol + locCol; //0-based
      //Ptr colbeg = colptr[locCol] - baseval;
      //Ptr colend = colptr[locCol+1] - baseval;
      Ptr colbeg = colptr.at(locCol) - baseval;
      Ptr colend = colptr.at(locCol+1) - baseval;
      Idx permCol = invp[col]-invpbaseval; // 0 based;
      //find destination processors
      Idx pdest; for(pdest = 0; pdest<mpisize; pdest++){ if(permCol>=newVertexDist[pdest]-baseval && permCol < newVertexDist[pdest+1]-baseval){ break;} }
      //Idx pdest = min( (Idx)mpisize-1, permCol / colPerProc);
      sizes[pdest] += (colend - colbeg)*sizeof(Idx) + sizeof(Ptr) + sizeof(Idx); //extra 2 are for the count of elements and the permuted columns

      //now permute rows
      for(Ptr jptr = colbeg; jptr<colend; jptr++){
        //Idx row = rowind[jptr] - baseval; //0 based
        Idx row = rowind.at(jptr) - baseval; //0 based
        Idx permRow = invp[row] - invpbaseval; // 0 based
        //rowind[jptr] = permRow + baseval;
        rowind.at(jptr) = permRow + baseval;
      }
    }


    //First allgatherv to get the receive sizes
    SYMPACK::vector<int> displs(mpisize+1,0);
    SYMPACK::vector<int> rsizes(mpisize,0);
    SYMPACK::vector<int> rdispls(mpisize+1,0);

    MPI_Alltoall(&sizes[0],sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,comm);
    //logfileptr->OFS()<<rsizes<<endl;
    //    MPI_Alltoallv(&sizes[0],&displs[0],&rdispls[0],MPI_BYTE,&rsizes[0],&displs[0],&rdispls[0],MPI_BYTE,comm);
    //logfileptr->OFS()<<rsizes<<endl;



    displs[0] = 0;
    std::partial_sum(sizes.begin(),sizes.end(),displs.begin()+1);
    Int totSend = displs.back();//std::accumulate(sizes.begin(),sizes.end(),0);
    SYMPACK::vector<char> sbuf(totSend);

    //pack
    for(Idx locCol = 0; locCol<LocalVertexCount(); locCol++){
      Idx col = firstCol + locCol; 
      //Ptr colbeg = colptr[locCol]-baseval;
      //Ptr colend = colptr[locCol+1]-baseval;
      Ptr colbeg = colptr.at(locCol) - baseval;
      Ptr colend = colptr.at(locCol+1) - baseval;
      Idx permCol = invp[col]-invpbaseval; // perm is 1 based;
      //find destination processors
      Idx pdest; for(pdest = 0; pdest<mpisize; pdest++){ if(permCol>=newVertexDist[pdest]-baseval && permCol < newVertexDist[pdest+1]-baseval){ break;} }
      //Idx pdest = min((Idx)mpisize-1, permCol / colPerProc);

      int & pos = displs[pdest];
      Idx * pPermCol = (Idx*)&sbuf[pos];
      pos+=sizeof(Idx);
      Ptr * pRowsCnt = (Ptr*)&sbuf[pos];
      pos+=sizeof(Ptr);
      Idx * pPermRows = (Idx*)&sbuf[pos];
      pos+=(colend-colbeg)*sizeof(Idx);


      *pPermCol = permCol;
      *pRowsCnt = (colend - colbeg);
      //std::copy(&rowind[0]+colbeg ,&rowind[0]+colend, pPermRows );
      std::copy(&rowind.at(colbeg) ,&rowind.at(colend-1)+1, pPermRows );


      //      logfileptr->OFS()<<*pPermCol<<" ("<<locCol<<"): ";
      //      for(Ptr jptr = 0; jptr<*pRowsCnt; jptr++){
      //        logfileptr->OFS()<<pPermRows[jptr]<<" ";
      //      }
      //      logfileptr->OFS()<<endl;


    }


    //re compute send displs
    displs[0] = 0;
    std::partial_sum(sizes.begin(),sizes.end(),&displs[1]);

    //now recompute receiv displs with actual sizes 
    rdispls[0] = 0;
    std::partial_sum(rsizes.begin(),rsizes.end(),&rdispls[1]);
    Ptr totRecv = rdispls.back();//std::accumulate(rsizes.begin(),rsizes.end(),0);
    SYMPACK::vector<char> rbuf(totRecv);

    MPI_Alltoallv(&sbuf[0],&sizes[0],&displs[0],MPI_BYTE,&rbuf[0],&rsizes[0],&rdispls[0],MPI_BYTE,comm);

    sbuf.clear();
    sizes.clear(); 
    displs.clear();
    rdispls.clear();


    colptr.resize(newVtxCount+1);
    //unpack first by overwriting colptr and then rowind
    Ptr rpos = 0; 
    colptr[0]=baseval;
    while(rpos<totRecv){
      Idx * permCol = (Idx*)&rbuf[rpos];
      rpos+=sizeof(Idx);
      Ptr * rowsCnt = (Ptr*)&rbuf[rpos];
      rpos+=sizeof(Ptr);
      Idx * permRows = (Idx*)&rbuf[rpos];
      rpos += (*rowsCnt)*sizeof(Idx);

      Idx locCol = *permCol - newFirstCol;
      //colptr[locCol+1] = *rowsCnt; 
      colptr.at(locCol+1) = *rowsCnt; 
    }
    std::partial_sum(colptr.begin(),colptr.end(),colptr.begin());

    //now fill rowind
    Ptr nnzLoc = colptr.back()-baseval;
    rowind.resize(nnzLoc);
    SYMPACK::vector<Ptr> colpos = colptr;
    rpos = 0;
    while(rpos<totRecv){
      Idx * permCol = (Idx*)&rbuf[rpos];
      rpos+=sizeof(Idx);
      Ptr * rowsCnt = (Ptr*)&rbuf[rpos];
      rpos+=sizeof(Ptr);
      Idx * permRows = (Idx*)&rbuf[rpos];
      rpos += (*rowsCnt)*sizeof(Idx);

      Idx locCol = *permCol - newFirstCol;

      //logfileptr->OFS()<<*permCol<<" ("<<locCol<<"): ";
      //for(Ptr jptr = 0; jptr<*rowsCnt; jptr++){
      //  logfileptr->OFS()<<permRows[jptr]<<" ";
      //}
      //logfileptr->OFS()<<endl;

      std::copy(permRows,permRows + *rowsCnt, &rowind[colpos[locCol]-baseval]);
      //colpos[locCol] += *rowsCnt;
      colpos.at(locCol) += *rowsCnt;
    }

    //copy newVertexDist into vertexDist
    std::copy(newVertexDist,newVertexDist+mpisize+1,&vertexDist[0]);

    if(sorted){
      SortEdges();
    }

    TIMER_STOP(PERMUTE);
  }




  void DistSparseMatrixGraph::RedistributeSupernodal(Int nsuper, Int * xsuper, Int * xsuperdist, Int * supMembership){
    int ismpi=0;
    MPI_Initialized( &ismpi);
    int isnull= (comm == MPI_COMM_NULL);

    if(!ismpi || isnull){
      //throw an exception
      throw std::logic_error("MPI communicator needs to be set.");
    }

    SparseMatrixGraph sg;
    AllGatherStructure(sg);

    Int nsuperLocal = xsuperdist[mpirank+1]-xsuperdist[mpirank];//(iam!=np-1)?nsuper/np:nsuper-iam*(nsuper/np);
    Int firstSnode = xsuperdist[mpirank]-1;//0 based;//iam*(nsuper/np)+1;
    Int lastSnode = firstSnode + nsuperLocal-1;

    Int firstCol = LocalFirstVertex()-baseval; //0 based
    Idx nLocal = LocalVertexCount();
    Ptr nnzLocal = LocalEdgeCount();

    //These numbers are inclusive
    //What is the last column I should own
    //What is the last column I currently own
    Idx supLastColumn = xsuper[firstSnode+nsuperLocal]-1-1;//0-based
    Idx curLastColumn = firstCol + nLocal -1; //0-based
    Idx supFirstColumn = xsuper[firstSnode]-1;//0-based
    Idx curFirstColumn = firstCol; //0-based


    //    logfileptr->OFS()<<"["<<curFirstColumn<<" - "<<curLastColumn<<"] to ["<<supFirstColumn<<" - "<<supLastColumn<<"]"<<endl;


    Int numColSentBefore = supFirstColumn>curFirstColumn?min(curLastColumn,supFirstColumn) - curFirstColumn +(curLastColumn<supFirstColumn?1:0):0;
    Ptr numRowSentBefore = numColSentBefore>0?colptr[min(curLastColumn,supFirstColumn)-firstCol]-colptr[curFirstColumn-firstCol]:0;
    Int numColSentAfter = curLastColumn>supLastColumn?curLastColumn-max(curFirstColumn,supLastColumn)+(curFirstColumn>supLastColumn?1:0):0;

    Idx startColAfter = max(curFirstColumn,supLastColumn)+(curFirstColumn>supLastColumn?0:1);
    Ptr numRowSentAfter = numColSentAfter>0?colptr[curLastColumn+1-firstCol]-colptr[startColAfter-firstCol]:0;


    //logfileptr->OFS()<<"numColSentBefore: "<<numColSentBefore<<endl;
    //logfileptr->OFS()<<" numColSentAfter: "<<numColSentAfter <<endl;
    //logfileptr->OFS()<<"numRowSentBefore: "<<numRowSentBefore<<endl;
    //logfileptr->OFS()<<" numRowSentAfter: "<<numRowSentAfter <<endl;



    //    logfileptr->OFS()<<"colptr was: "<<colptr<<endl;
    //    logfileptr->OFS()<<"rowind was: "<<rowind<<endl;


    //first transfer column pointers
    vector<int> ssizes(mpisize,0);
    vector<int> ssizesR(mpisize,0);
    vector<int> ssizesAfter(mpisize,0);
    vector<int> sdispls(mpisize+1,-1);
    vector<int> sdisplsR(mpisize+1,-1);
    int totSent = 0;
    int totSentR = 0;

    if(numColSentBefore>0){
      //logfileptr->OFS()<<"Sending BEFORE columns: { ";
      for(Idx col = curFirstColumn; col<curFirstColumn+numColSentBefore; col++){
        bassert(col>=curFirstColumn && col<=curLastColumn);
        Int snode = supMembership[col];
        //            Idx pdest = min(mpisize-1, (snode-1) / supPerProc);
        Int pdest = 0;
        for(pdest = 0; pdest<mpisize;pdest++){
          if(xsuperdist[pdest]<=snode && snode<xsuperdist[pdest+1]){
            break;
          }
        }

        //logfileptr->OFS()<<col+1<<" ";
        bassert(pdest!=mpirank);

        if(sdispls[pdest]==-1){
          sdispls[pdest] = (col-firstCol)*sizeof(Ptr);
        }
        totSent+=sizeof(Ptr);
        ssizes[pdest]+=sizeof(Ptr);

        Ptr nrows = (colptr[col-firstCol +1] - colptr[col-firstCol]);
        if(sdisplsR[pdest]==-1){
          sdisplsR[pdest] = (colptr[col-firstCol] - baseval)*sizeof(Idx);
        }
        totSentR+=nrows*sizeof(Idx);
        ssizesR[pdest]+=nrows*sizeof(Idx);
      }
      //logfileptr->OFS()<<"}"<<endl;

      for(int p =0; p<mpisize;p++){
        if(ssizes[p]>0){
          //add an extra colptr entry to send
          ssizes[p]+=sizeof(Ptr); 
          totSent += sizeof(Ptr);
        }
      }
    }

    if(numColSentAfter>0){
      //If I have some columns to send
      //logfileptr->OFS()<<"Sending AFTER columns: { ";
      for(Idx col = startColAfter; col<startColAfter+numColSentAfter; col++){
        bassert(col>=curFirstColumn && col<=curLastColumn);
        Int snode = supMembership[col];
        //            Idx pdest = min(mpisize-1, (snode-1) / supPerProc);
        Int pdest = 0;
        for(pdest = 0; pdest<mpisize;pdest++){
          if(xsuperdist[pdest]<=snode && snode<xsuperdist[pdest+1]){
            break;
          }
        }


        //logfileptr->OFS()<<col+1<<" ";
        bassert(pdest!=mpirank);
        if(sdispls[pdest]==-1){
          sdispls[pdest] = (col-firstCol)*sizeof(Ptr);
        }
        totSent+=sizeof(Ptr);
        ssizesAfter[pdest]+=sizeof(Ptr);

        Ptr nrows = (colptr[col-firstCol +1] - colptr[col-firstCol]);
        if(sdisplsR[pdest]==-1){
          sdisplsR[pdest] = (colptr[col-firstCol] - baseval)*sizeof(Idx);
        }
        totSentR+=nrows*sizeof(Idx);
        ssizesR[pdest]+=nrows*sizeof(Idx);
      }
      //logfileptr->OFS()<<"}"<<endl;
      for(int p =0; p<mpisize;p++){
        if(ssizesAfter[p]>0){
          //add an extra colptr entry to send
          ssizesAfter[p]+=sizeof(Ptr); 
          totSent += sizeof(Ptr);
        }
      }

      for(int p =0; p<mpisize;p++){
        ssizes[p]+=ssizesAfter[p];
      }


    }

    //    logfileptr->OFS()<<"ssizes: "<<ssizes<<endl;
    //    logfileptr->OFS()<<"sdispls orig : "<<sdispls<<endl;


    vector<int> rsizes(mpisize,0);
    vector<int> rdispls(mpisize+1,0);
    MPI_Alltoall(&ssizes[0],sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,comm);
    rdispls[0] = 0;
    std::partial_sum(rsizes.begin(),rsizes.end(),&rdispls[1]);

    vector<int> rsizesR(mpisize,0);
    vector<int> rdisplsR(mpisize+1,0);
    MPI_Alltoall(&ssizesR[0],sizeof(int),MPI_BYTE,&rsizesR[0],sizeof(int),MPI_BYTE,comm);
    rdisplsR[0] = 0;
    std::partial_sum(rsizesR.begin(),rsizesR.end(),&rdisplsR[1]);

    //    logfileptr->OFS()<<"rsizes: "<<rsizes<<endl;
    //    logfileptr->OFS()<<"rdispls orig : "<<rdispls<<endl;


    SYMPACK::vector<char> recvBuf(rdispls.back()); 
    MPI_Alltoallv(&colptr[0], &ssizes[0], &sdispls[0], MPI_BYTE,
        &recvBuf[0], &rsizes[0], &rdispls[0], MPI_BYTE,
        comm);

    Int recvColCntBefore = 0;
    Int recvColCnt = 0;
    for(int p=0;p<mpisize;p++){
      if(rsizes[p]>0){
        if(p<mpirank){
          recvColCntBefore+=rsizes[p]/sizeof(Ptr)-1;
        }
        recvColCnt+=rsizes[p]/sizeof(Ptr)-1;
      }
    }

    //now shift the content if needed
    Ptr colptrSize = nLocal-numColSentBefore-numColSentAfter;
    if(numColSentBefore>0 && colptrSize>0){
      std::copy(colptr.begin()+numColSentBefore,colptr.begin()+nLocal+1,colptr.begin());
    }

    colptr.resize(colptrSize+recvColCnt+1);

    if(recvColCntBefore>0){
      std::copy_backward(colptr.begin(),colptr.begin()+colptrSize+1,colptr.begin()+recvColCntBefore+colptrSize+1);
    }


    //copy from recvBuf
    //if I have received something from ranks before me
    Idx colptrPos =0;
    Idx ptrBufPos = 0;
    Ptr * ptrBuf = (Ptr*)&recvBuf[0];


    //        logfileptr->OFS()<<"Recv colptr: "<<endl;
    //    for(int p=0;p<mpisize;p++){
    //      if(rsizes[p]>0){
    //        int numCols = (int)(rsizes[p]/sizeof(Ptr))-1;
    //        logfileptr->OFS()<<"P"<<p<<": ";
    //        for(int col = 0;col<=numCols;col++){
    //          logfileptr->OFS()<<ptrBuf[col]<<" ";
    //        }
    //        logfileptr->OFS()<<endl;
    //      }
    //    }

    Ptr offset = 1;
    for(int p=0;p<mpirank;p++){
      if(rsizes[p]>0){
        int numCols = (int)(rsizes[p]/sizeof(Ptr))-1;
        for(int col = numCols-1;col>=0;col--){
          colptr[colptrPos+col] = ptrBuf[ptrBufPos+col] - ptrBuf[ptrBufPos] + offset;
        }
        offset = ptrBuf[ptrBufPos+numCols] - ptrBuf[ptrBufPos] + offset;

        ptrBufPos+=numCols+1;
        colptrPos+=numCols;
      }
    }

    //apply the offset to my columns. 1-based because Idx may be unsigned
    for(Idx col = colptrSize+1;col>0;col--){
      colptr[colptrPos+col-1] = colptr[colptrPos+col-1] - colptr[colptrPos] + offset;
    }
    offset = colptr[colptrPos+colptrSize] - colptr[colptrPos] + offset;
    colptrPos+=colptrSize;

    //if I have something from ranks after me
    for(int p=mpirank+1;p<mpisize;p++){
      if(rsizes[p]>0){
        int numCols = (int)(rsizes[p]/sizeof(Ptr))-1;
        for(int col = numCols-1;col>=0;col--){
          colptr[colptrPos+col] = ptrBuf[ptrBufPos+col] - ptrBuf[ptrBufPos] + offset;
        }
        offset = ptrBuf[ptrBufPos+numCols] - ptrBuf[ptrBufPos] + offset;

        ptrBufPos+=numCols+1;
        colptrPos+=numCols;
      }
    }
    bassert(colptrPos == colptrSize+recvColCnt);
    //colptr[colptrPos] = offset;




    //#ifdef _USE_COREDUMPER_
    //    std::stringstream corename;
    //    corename << "core.sympack." << iam;
    //    WriteCoreDump(corename.str().c_str());
    //#endif




    recvBuf.resize(rdisplsR.back()); 
    MPI_Alltoallv(&rowind[0], &ssizesR[0], &sdisplsR[0], MPI_BYTE,
        &recvBuf[0], &rsizesR[0], &rdisplsR[0], MPI_BYTE,
        comm);
    Int recvCntBefore = 0;
    Int recvCnt = 0;
    for(int p=0;p<mpisize;p++){
      if(rsizesR[p]>0){
        if(p<mpirank){
          recvCntBefore+=rsizesR[p]/sizeof(Idx);
        }
        recvCnt+=rsizesR[p]/sizeof(Idx);
      }
    }



    //now shift the content if needed
    Ptr rowindSize = nnzLocal-numRowSentBefore-numRowSentAfter;
    if(numRowSentBefore>0 && rowindSize>0){
      std::copy(rowind.begin()+numRowSentBefore,rowind.begin()+nnzLocal,rowind.begin());
    }

    rowind.resize(rowindSize+recvCnt);

    //finish colptr
    if(colptr.size()>0){
      colptr.back()=rowind.size()+baseval;
    }
    //logfileptr->OFS()<<"colptr now is: "<<colptr<<endl;

    if(recvCntBefore>0){
      std::copy_backward(rowind.begin(),rowind.begin()+rowindSize,rowind.begin()+recvCntBefore+rowindSize);
    }


    //copy from recvBuf

    //if I have received something from ranks before me
    Idx rowindPos =0;
    Idx idxBufPos = 0;
    Idx * idxBuf = (Idx*)&recvBuf[0];



    //logfileptr->OFS()<<"Recv rowind: "<<endl;
    //for(int p=0;p<mpisize;p++){
    //  if(rsizesR[p]>0){
    //    int numRows = (int)(rsizesR[p]/sizeof(Idx));
    //    logfileptr->OFS()<<"P"<<p<<": ";
    //    for(int row = 0;row<=numRows;row++){
    //      logfileptr->OFS()<<idxBuf[row]<<" ";
    //    }
    //    logfileptr->OFS()<<endl;
    //  }
    //}




    for(int p=0;p<mpirank;p++){
      if(rsizesR[p]>0){
        int numRows = (int)(rsizesR[p]/sizeof(Idx));
        std::copy(&idxBuf[idxBufPos],&idxBuf[idxBufPos]+numRows,&rowind[rowindPos]);
        idxBufPos+=numRows;
        rowindPos+=numRows;
      }
    }

    rowindPos+=rowindSize;

    //if I have something from ranks after me
    for(int p=mpirank+1;p<mpisize;p++){
      if(rsizesR[p]>0){
        int numRows = (int)(rsizesR[p]/sizeof(Idx));
        std::copy(&idxBuf[idxBufPos],&idxBuf[idxBufPos]+numRows,&rowind[rowindPos]);
        idxBufPos+=numRows;
        rowindPos+=numRows;
      }
    }
    bassert(rowindPos == rowindSize+recvCnt);
    //logfileptr->OFS()<<"rowind now is: "<<rowind<<endl;

    //logfileptr->OFS()<<"Vertex distribution was: "<<vertexDist<<endl;
    supLastColumn = supLastColumn + baseval +1;
    MPI_Allgather(&supLastColumn,sizeof(supLastColumn),MPI_BYTE,&vertexDist[1],sizeof(supLastColumn),MPI_BYTE,comm);
    //logfileptr->OFS()<<"and is now: "<<vertexDist<<endl;


    //SparseMatrixGraph sgr;
    //AllGatherStructure(sgr);

    //logfileptr->OFS()<<sg.colptr<<endl;
    //logfileptr->OFS()<<sgr.colptr<<endl;
    //logfileptr->OFS()<<sg.rowind<<endl;
    //logfileptr->OFS()<<sgr.rowind<<endl;

  }







  void DistSparseMatrixGraph::ExpandSymmetric(){
    TIMER_START(EXPAND);
    if(!bIsExpanded){

      int ismpi=0;
      MPI_Initialized( &ismpi);
      int isnull= (comm == MPI_COMM_NULL);

      if(!ismpi || isnull){
        //throw an exception
        throw std::logic_error("MPI communicator needs to be set.");
      }


      if(!sorted){
        SetSorted(1);
      }

      MPI_Op MPI_SYMPACK_PTR_SUM; 
      MPI_Op MPI_SYMPACK_PTR_MAX; 
      MPI_Datatype MPI_SYMPACK_PTR; 
      MPI_Type_contiguous( sizeof(Ptr), MPI_BYTE, &MPI_SYMPACK_PTR ); 
      MPI_Type_commit( &MPI_SYMPACK_PTR ); 
      MPI_Op_create( PtrSum, true, &MPI_SYMPACK_PTR_SUM ); 
      MPI_Op_create( PtrMax, true, &MPI_SYMPACK_PTR_MAX ); 

      Idx N = size; 

      Idx locColCnt = LocalVertexCount();
      //make copies first
      SYMPACK::vector<Ptr> prevColptr;
      SYMPACK::vector<Idx> prevRowind;
      prevColptr.swap(colptr);
      prevRowind.swap(rowind);

      Ptr * pcolptr = &prevColptr[0];
      Idx * prowind = &prevRowind[0];
      Ptr locNNZ = prevRowind.size();

      //first generate the expanded distributed structure
      //            Int firstLocCol = (mpirank)*colPerProc; //0 based
      Idx firstLocCol = LocalFirstVertex()-baseval; //0 based
      Idx maxLocN = 0;
      for(int p = 0; p<mpisize;p++){maxLocN = max(maxLocN, vertexDist[p+1]-vertexDist[p]);}
      //               max(locColCnt,N-(mpisize-1)*colPerProc); // can be 0
      SYMPACK::vector<Ptr> remote_colptr(maxLocN+1);
      SYMPACK::vector<Idx> remote_rowind;
      SYMPACK::vector<Ptr> remote_rowindPos(maxLocN+1);
      SYMPACK::vector<Ptr> curPos(locColCnt);
      SYMPACK::vector<Ptr> prevPos(locColCnt);

      std::copy(pcolptr,pcolptr+locColCnt,curPos.begin());
      for(Int prow = 0; prow<mpisize; prow++){
        Int firstRemoteCol = vertexDist[prow]-baseval;//(prow)*colPerProc; // 0 based
        Int pastLastRemoteCol = vertexDist[prow+1]-baseval;//prow==mpisize-1?N:(prow+1)*colPerProc; // 0 based
        Ptr remColCnt = pastLastRemoteCol - firstRemoteCol;
        Ptr maxExtraNNZ = 0; 
        Ptr extraNNZ = 0;
        //receive from all the previous processor
        //receive extra nnzcnt...

        std::fill(remote_colptr.begin(),remote_colptr.end(),0);
        if(mpirank<=prow){
          //backup the position in each column
          std::copy(curPos.begin(),curPos.end(),prevPos.begin());
          //use remote_colptr to store the number of nnz per col first
          //loop through my local columns and figure out the extra nonzero per col on prow
          for(Idx locCol = 0 ; locCol< locColCnt; locCol++){
            Idx col = firstLocCol + locCol;  // 0 based
            Ptr colbeg = curPos[locCol]-baseval; //now 0 based
            Ptr colend = pcolptr[locCol+1]-baseval; // now 0 based
            for(Ptr rptr = colbeg ; rptr< colend ; rptr++ ){
              Idx row = prowind[rptr]-baseval; //0 based
              assert(row>=firstRemoteCol);
              if(row>col){
                //this goes onto prow
                if(row<pastLastRemoteCol){
                  //this is shifted by one to compute the prefix sum without copying afterwards
                  remote_colptr[row-firstRemoteCol+1]++;
                  extraNNZ++;
                }
                else{
                  break;
                }
              }
              curPos[locCol]++; // baseval based
            }
          }

          //put the local sum into the first element of the array
          remote_colptr[0] = 0;
          for(Idx p = 1; p<remColCnt+1;p++){ remote_colptr[0] += remote_colptr[p];}


          if(mpirank==prow){
            //we can now receive the number of NZ per col into the expanded pcolptr
            assert(locColCnt == remColCnt);
            colptr.assign(remColCnt+1,0);
            //this array will contain the max element in our custom reduce
            colptr[0] = 0;
            TIMER_START(REDUCE);
            MPI_Reduce(&remote_colptr[1],&colptr[1],remColCnt,MPI_SYMPACK_PTR,MPI_SYMPACK_PTR_SUM,prow,comm);
            MPI_Reduce(&remote_colptr[0],&colptr[0],1,MPI_SYMPACK_PTR,MPI_SYMPACK_PTR_MAX,prow,comm);
            TIMER_STOP(REDUCE);

            maxExtraNNZ = colptr[0];
            remote_rowind.resize(maxExtraNNZ);

            for(Idx locCol = 0 ; locCol< locColCnt; locCol++){
              Ptr colbeg = pcolptr[locCol]-baseval; //now 0 based
              Ptr colend = pcolptr[locCol+1]-baseval; // now 0 based
              colptr[locCol+1] += colend - colbeg; //- 1+ keepDiag;
              //At this point expColptr[locCol+1] contains the total number of NNZ of locCol 
              //we can now compute the expanded pcolptr
            }

            colptr[0] = baseval;
            for(Idx col = 1;col<remColCnt+1;col++){ colptr[col]+=colptr[col-1]; }
          }
          else{
            TIMER_START(REDUCE);
            MPI_Reduce(&remote_colptr[1],NULL,remColCnt,MPI_SYMPACK_PTR,MPI_SYMPACK_PTR_SUM,prow,comm);
            MPI_Reduce(&remote_colptr[0],NULL,1,MPI_SYMPACK_PTR,MPI_SYMPACK_PTR_MAX,prow,comm);
            TIMER_STOP(REDUCE);
            remote_rowind.resize(extraNNZ);
          }




///TODO


          /**************     Compute remote_colptr from the local nnz per col ***********/
          //compute a prefix sum of the nnz per column to get the new pcolptr
          remote_colptr[0] = baseval;
          for(Idx col = 1;col<=remColCnt;col++){ remote_colptr[col]+=remote_colptr[col-1]; }

          /**************     Fill remote_rowind now ****************/
          //make a copy of pcolptr in order to track the current position in each column
          std::copy(&remote_colptr[0],&remote_colptr[0]+remColCnt+1,remote_rowindPos.begin());
          //loop through my local columns and figure out the extra nonzero per col on remote processors
          for(Idx locCol = 0 ; locCol< locColCnt; locCol++){
            Idx col = firstLocCol + locCol;  // 0 based
            Ptr colbeg = prevPos[locCol]-baseval; //now 0 based
            Ptr colend = curPos[locCol]-baseval; // now 0 based

            for(Ptr rptr = colbeg ; rptr< colend ; rptr++ ){
              Idx row = prowind[rptr]-baseval; //0 based
              if(row>col){
                //this goes onto prow
                Idx locRow = row - firstRemoteCol;      
                remote_rowind[ remote_rowindPos[locRow]++ - baseval ] = col + baseval;
              }
            }
          }

#ifdef DEBUG
          logfileptr->OFS()<<"remote_colptr ";
          for(Idx col = 0;col<=remColCnt;col++){
            logfileptr->OFS()<<remote_colptr[col]<<" ";
          }
          logfileptr->OFS()<<endl;

          logfileptr->OFS()<<"remote_rowind ";
          for(Ptr col = 0;col<extraNNZ;col++){
            logfileptr->OFS()<<remote_rowind[col]<<" ";
          }
          logfileptr->OFS()<<endl;
#endif

          if(prow==mpirank){
#ifdef DEBUG
            logfileptr->OFS()<<"expColptr "<<colptr<<endl;
#endif
            if(colptr.size()>0){
              rowind.resize(colptr.back()-baseval);
            }
            std::copy(colptr.begin(),colptr.end(),remote_rowindPos.begin()); 
            for(Idx locCol = 0 ; locCol< locColCnt; locCol++){
              Idx col = firstLocCol + locCol;  // 0 based
              Ptr colbeg = pcolptr[locCol]-baseval; //now 0 based
              Ptr colend = pcolptr[locCol+1]-baseval; // now 0 based
              Ptr & pos = remote_rowindPos[locCol];

              //copy the local lower triangular part into the expanded structure
              for(Ptr rptr = colbeg ; rptr< colend ; rptr++ ){
                Idx row = prowind[rptr]-baseval; //0 based
                if(col!=row || keepDiag){
                  rowind[pos++ - baseval] = row + baseval;  
                }
              }

              //copy the local extra NNZ into the expanded structure
              colbeg = remote_colptr[locCol]-baseval; //now 0 based
              colend = remote_colptr[locCol+1]-baseval; // now 0 based 
              for(Ptr rptr = colbeg ; rptr< colend ; rptr++ ){
                Idx row = remote_rowind[rptr]-baseval; //0 based
                rowind[pos++ - baseval] = row + baseval;  
              }
            }

            TIMER_START(RECV);
            for(Int pcol = 0; pcol<prow; pcol++){
              //Use an MPI_Gatherv instead ? >> memory usage : p * n/p
              //Do mpi_recv from any, anytag for pcolptr and then do the matching rowind ?

              //logfileptr->OFS()<<"P"<<mpirank<<" receives pcolptr from P"<<pcol<<endl;
              //receive colptrs...
              MPI_Status status;
              MPI_Recv(&remote_colptr[0],(remColCnt+1)*sizeof(Ptr),MPI_BYTE,MPI_ANY_SOURCE,prow,comm,&status);

              //logfileptr->OFS()<<"P"<<mpirank<<" receives rowind from P"<<pcol<<endl;
              //receive rowinds...
              MPI_Recv(&remote_rowind[0],maxExtraNNZ*sizeof(Idx),MPI_BYTE,status.MPI_SOURCE,prow,comm,MPI_STATUS_IGNORE);

              TIMER_START(PROCESSING_RECV_DATA);
              //logfileptr->OFS()<<"P"<<mpirank<<" done receiving from P"<<pcol<<endl;
              for(Idx locCol = 0 ; locCol< locColCnt; locCol++){
                Idx col = firstLocCol + locCol;  // 0 based
                //copy the extra NNZ into the expanded structure
                Ptr colbeg = remote_colptr[locCol]-baseval; //now 0 based
                Ptr colend = remote_colptr[locCol+1]-baseval; // now 0 based 
                Ptr & pos = remote_rowindPos[locCol];
                for(Ptr rptr = colbeg ; rptr< colend ; rptr++ ){
                  Idx row = remote_rowind[rptr]-baseval; //0 based
                  rowind[pos++ - baseval] = row + baseval;  
                }
              }
              TIMER_STOP(PROCESSING_RECV_DATA);
            }
            TIMER_STOP(RECV);

            if(sorted){ 
              SortEdges();
            }

#ifdef DEBUG
            logfileptr->OFS()<<"expRowind "<<rowind<<endl;
            //logfileptr->OFS()<<"true expRowind "<<Global.expRowind<<endl;
#endif

          }
          else{
            TIMER_START(SEND);
            MPI_Send(&remote_colptr[0],(remColCnt+1)*sizeof(Ptr),MPI_BYTE,prow,prow,comm);
            MPI_Send(&remote_rowind[0],extraNNZ*sizeof(Idx),MPI_BYTE,prow,prow,comm);
            TIMER_STOP(SEND);
          }

        }
        else{
          TIMER_START(REDUCE);
            MPI_Reduce(&remote_colptr[1],NULL,remColCnt,MPI_SYMPACK_PTR,MPI_SYMPACK_PTR_SUM,prow,comm);
            MPI_Reduce(&remote_colptr[0],NULL,1,MPI_SYMPACK_PTR,MPI_SYMPACK_PTR_MAX,prow,comm);
          TIMER_STOP(REDUCE);
        }

      }

      MPI_Op_free(&MPI_SYMPACK_PTR_MAX);
      MPI_Op_free(&MPI_SYMPACK_PTR_SUM);
      MPI_Type_free(&MPI_SYMPACK_PTR);

      //nnz = LocalEdgeCount();
      bIsExpanded =true;
    }
    TIMER_STOP(EXPAND);
  }

  void DistSparseMatrixGraph::AllGatherStructure(SparseMatrixGraph & g){
    int ismpi=0;
    MPI_Initialized( &ismpi);
    int isnull= (comm == MPI_COMM_NULL);

    if(!ismpi || isnull){
      //throw an exception
      throw std::logic_error("MPI communicator needs to be set.");
    }


    g.size = size;
    g.SetBaseval(baseval);  
    g.SetSorted(sorted);
    g.SetKeepDiag(keepDiag);
    g.bIsExpanded = bIsExpanded;

    //get other proc vertex counts
    Idx localVertexCnt = LocalVertexCount();
    SYMPACK::vector<Idx> remoteVertexCnt(mpisize,0);
    MPI_Allgather(&localVertexCnt,sizeof(localVertexCnt),MPI_BYTE,&remoteVertexCnt[0],sizeof(localVertexCnt),MPI_BYTE,comm);
    Idx totalVertexCnt = std::accumulate(remoteVertexCnt.begin(),remoteVertexCnt.end(),0,std::plus<Idx>());
    g.colptr.resize(totalVertexCnt+1);
    //compute receive displacements
    SYMPACK::vector<int> rsizes(mpisize,0);
    for(int p = 0; p<mpisize;p++){rsizes[p] = (int)remoteVertexCnt[p]*sizeof(Ptr);}
    SYMPACK::vector<int> rdispls(mpisize+1,0);
    rdispls[0]=0;
    std::partial_sum(rsizes.begin(),rsizes.end(),rdispls.begin()+1);
    MPI_Allgatherv(&colptr[0],localVertexCnt*sizeof(Ptr),MPI_BYTE,&g.colptr[0],&rsizes[0],&rdispls[0],MPI_BYTE,comm);


    Ptr localEdgeCnt = LocalEdgeCount();
    SYMPACK::vector<Ptr> remoteEdgeCnt(mpisize,0);
    MPI_Allgather(&localEdgeCnt,sizeof(localEdgeCnt),MPI_BYTE,&remoteEdgeCnt[0],sizeof(localEdgeCnt),MPI_BYTE,comm);
    Ptr totalEdgeCnt = std::accumulate(remoteEdgeCnt.begin(),remoteEdgeCnt.end(),0,std::plus<Ptr>());

    g.rowind.resize(totalEdgeCnt);

    //compute receive displacements
    rsizes.assign(mpisize,0);
    for(int p = 0; p<mpisize;p++){rsizes[p] = (int)remoteEdgeCnt[p]*sizeof(Idx);}
    rdispls[0]=0;
    std::partial_sum(rsizes.begin(),rsizes.end(),rdispls.begin()+1);
    MPI_Allgatherv(&rowind[0],localEdgeCnt*sizeof(Idx),MPI_BYTE,&g.rowind[0],&rsizes[0],&rdispls[0],MPI_BYTE,comm);




    //fix colptr
    Idx pos = remoteVertexCnt[0];
    Ptr offset = 0;
    for(int p=1;p<mpisize;p++){
      offset+=remoteEdgeCnt[p-1]; 
      for(Idx lv = 0; lv < remoteVertexCnt[p]; lv++){
        g.colptr[pos++] += offset;//remoteEdgeCnt[p-1];//(vertexDist[p] - baseval);
      }
    }
    if(g.colptr.size()>0){
      g.colptr.back()=totalEdgeCnt + baseval;
    }
  }



  void DistSparseMatrixGraph::GatherStructure(SparseMatrixGraph & g, int proot){
    int ismpi=0;
    MPI_Initialized( &ismpi);
    int isnull= (comm == MPI_COMM_NULL);

    if(!ismpi || isnull){
      //throw an exception
      throw std::logic_error("MPI communicator needs to be set.");
    }

    g.size = size;
    g.SetBaseval(baseval);  
    if(iam==proot){
    g.SetSorted(sorted);
    g.SetKeepDiag(keepDiag);
    g.bIsExpanded = bIsExpanded;
    }

    //get other proc vertex counts
    Idx localVertexCnt = LocalVertexCount();
    SYMPACK::vector<Idx> remoteVertexCnt(mpisize,0);
    MPI_Allgather(&localVertexCnt,sizeof(localVertexCnt),MPI_BYTE,&remoteVertexCnt[0],sizeof(localVertexCnt),MPI_BYTE,comm);
    if(iam==proot){
      Idx totalVertexCnt = std::accumulate(remoteVertexCnt.begin(),remoteVertexCnt.end(),0,std::plus<Idx>());
      g.colptr.resize(totalVertexCnt+1);


      //compute receive displacements
      SYMPACK::vector<int> rsizes(mpisize,0);
      for(int p = 0; p<mpisize;p++){rsizes[p] = (int)remoteVertexCnt[p]*sizeof(Ptr);}
      SYMPACK::vector<int> rdispls(mpisize+1,0);
      rdispls[0]=0;
      std::partial_sum(rsizes.begin(),rsizes.end(),rdispls.begin()+1);
      MPI_Gatherv(&colptr[0],localVertexCnt*sizeof(Ptr),MPI_BYTE,&g.colptr[0],&rsizes[0],&rdispls[0],MPI_BYTE,proot,comm);
    }
    else{
    MPI_Gatherv(&colptr[0],localVertexCnt*sizeof(Ptr),MPI_BYTE,NULL,NULL,NULL,MPI_BYTE,proot,comm);
    }


    Ptr localEdgeCnt = LocalEdgeCount();
    SYMPACK::vector<Ptr> remoteEdgeCnt(mpisize,0);
    MPI_Gather(&localEdgeCnt,sizeof(localEdgeCnt),MPI_BYTE,&remoteEdgeCnt[0],sizeof(localEdgeCnt),MPI_BYTE,proot,comm);
    Ptr totalEdgeCnt = std::accumulate(remoteEdgeCnt.begin(),remoteEdgeCnt.end(),0,std::plus<Ptr>());


    if(iam==proot){
    g.rowind.resize(totalEdgeCnt);
    //compute receive displacements
      SYMPACK::vector<int> rsizes(mpisize,0);
      SYMPACK::vector<int> rdispls(mpisize+1,0);
    for(int p = 0; p<mpisize;p++){rsizes[p] = (int)remoteEdgeCnt[p]*sizeof(Idx);}
    rdispls[0]=0;
    std::partial_sum(rsizes.begin(),rsizes.end(),rdispls.begin()+1);
    MPI_Gatherv(&rowind[0],localEdgeCnt*sizeof(Idx),MPI_BYTE,&g.rowind[0],&rsizes[0],&rdispls[0],MPI_BYTE,proot,comm);
    }
    else{
    MPI_Gatherv(&rowind[0],localEdgeCnt*sizeof(Idx),MPI_BYTE,NULL,NULL,NULL,MPI_BYTE,proot,comm);
    }

    if(iam==proot){
      //fix colptr
      Idx pos = remoteVertexCnt[0];
      Ptr offset = 0;
      for(int p=1;p<mpisize;p++){
        offset+=remoteEdgeCnt[p-1]; 
        for(Idx lv = 0; lv < remoteVertexCnt[p]; lv++){
          g.colptr[pos++] += offset;//remoteEdgeCnt[p-1];//(vertexDist[p] - baseval);
        }
      }
      if(g.colptr.size()>0){
        g.colptr.back()=totalEdgeCnt + baseval;
      }
    }
  }






















  //  void DistSparseMatrixGraph::GetLColRowCount(ETree & tree, Ordering & aOrder, SYMPACK::vector<Int> & cc, SYMPACK::vector<Int> & rc){
  //    //The tree need to be postordered
  //    if(!tree.IsPostOrdered()){
  //      tree.PostOrderTree(aOrder);
  //    }
  //
  //    ExpandSymmetric();
  //
  //    TIMER_START(GetColRowCount_Classic);
  //    cc.resize(size);
  //    rc.resize(size);
  //
  //    SYMPACK::vector<Idx> level(size+1);
  //    SYMPACK::vector<Int> weight(size+1);
  //    SYMPACK::vector<Idx> fdesc(size+1);
  //    SYMPACK::vector<Idx> nchild(size+1);
  //    SYMPACK::vector<Idx> set(size);
  //    SYMPACK::vector<Idx> prvlf(size);
  //    SYMPACK::vector<Idx> prvnbr(size);
  //
  //
  //    Idx xsup = 1;
  //    level[0] = 0;
  //    for(Idx k = size; k>=1; --k){
  //      rc[k-1] = 1;
  //      cc[k-1] = 0;
  //      set[k-1] = k;
  //      prvlf[k-1] = 0;
  //      prvnbr[k-1] = 0;
  //      level[k] = level[tree.PostParent(k-1)] + 1;
  //      weight[k] = 1;
  //      fdesc[k] = k;
  //      nchild[k] = 0;
  //    }
  //
  //    nchild[0] = 0;
  //    fdesc[0] = 0;
  //    for(Idx k =1; k<size; ++k){
  //      Idx parent = tree.PostParent(k-1);
  //      weight[parent] = 0;
  //      ++nchild[parent];
  //      Idx ifdesc = fdesc[k];
  //      if  ( ifdesc < fdesc[parent] ) {
  //        fdesc[parent] = ifdesc;
  //      }
  //    }
  //
  //    for(Idx lownbr = 1; lownbr<=size; ++lownbr){
  //      Int lflag = 0;
  //      Idx ifdesc = fdesc[lownbr];
  //      Idx oldnbr = aOrder.perm[lownbr-1];
  //      Ptr jstrt = expColptr[oldnbr-1];
  //      Ptr jstop = expColptr[oldnbr] - 1;
  //
  //
  //      //           -----------------------------------------------
  //      //           for each ``high neighbor'', hinbr of lownbr ...
  //      //           -----------------------------------------------
  //      for(Ptr j = jstrt; j<=jstop;++j){
  //        Idx hinbr = expRowind[j-1];
  //        hinbr = aOrder.invp[hinbr-1];
  //        if  ( hinbr > lownbr )  {
  //          if  ( ifdesc > prvnbr[hinbr-1] ) {
  //            //                       -------------------------
  //            //                       increment weight[lownbr].
  //            //                       -------------------------
  //            ++weight[lownbr];
  //            Idx pleaf = prvlf[hinbr-1];
  //            //                       -----------------------------------------
  //            //                       if hinbr has no previous ``low neighbor'' 
  //            //                       then ...
  //            //                       -----------------------------------------
  //            if  ( pleaf == 0 ) {
  //              //                           -----------------------------------------
  //              //                           ... accumulate lownbr-->hinbr path length 
  //              //                               in rowcnt[hinbr].
  //              //                           -----------------------------------------
  //              rc[hinbr-1] += level[lownbr] - level[hinbr];
  //            }
  //            else{
  //              //                           -----------------------------------------
  //              //                           ... otherwise, lca <-- find[pleaf], which 
  //              //                               is the least common ancestor of pleaf 
  //              //                               and lownbr.
  //              //                               (path halving.)
  //              //                           -----------------------------------------
  //              Idx last1 = pleaf;
  //              Idx last2 = set[last1-1];
  //              Idx lca = set[last2-1];
  //              while(lca != last2){
  //                set[last1-1] = lca;
  //                last1 = lca;
  //                last2 = set[last1-1];
  //                lca = set[last2-1];
  //              }
  //              //                           -------------------------------------
  //              //                           accumulate pleaf-->lca path length in 
  //              //                           rowcnt[hinbr].
  //              //                           decrement weight(lca).
  //              //                           -------------------------------------
  //              rc[hinbr-1] += level[lownbr] - level[lca];
  //              --weight[lca];
  //            }
  //            //                       ----------------------------------------------
  //            //                       lownbr now becomes ``previous leaf'' of hinbr.
  //            //                       ----------------------------------------------
  //            prvlf[hinbr-1] = lownbr;
  //            lflag = 1;
  //          }
  //          //                   --------------------------------------------------
  //          //                   lownbr now becomes ``previous neighbor'' of hinbr.
  //          //                   --------------------------------------------------
  //          prvnbr[hinbr-1] = lownbr;
  //        }
  //      }
  //      //           ----------------------------------------------------
  //      //           decrement weight ( parent[lownbr] ).
  //      //           set ( p[lownbr] ) <-- set ( p[lownbr] ) + set[xsup].
  //      //           ----------------------------------------------------
  //      Idx parent = tree.PostParent(lownbr-1);
  //      --weight[parent];
  //
  //
  //      //merge the sets
  //      if  ( lflag == 1  || nchild[lownbr] >= 2 ) {
  //        xsup = lownbr;
  //      }
  //      set[xsup-1] = parent;
  //    }
  //
  //
  //
  //#ifdef _DEBUG_
  //    logfileptr->OFS()<<"deltas "<<weight<<std::endl;
  //#endif
  //
  //    for(Int k = 1; k<=size; ++k){
  //      Int temp = cc[k-1] + weight[k];
  //      cc[k-1] = temp;
  //      Int parent = tree.PostParent(k-1);
  //      if  ( parent != 0 ) {
  //        cc[parent-1] += temp;
  //      }
  //    }
  //
  //    //      logfileptr->OFS()<<"column counts "<<cc<<std::endl;
  //
  //    TIMER_STOP(GetColRowCount_Classic);
  //  }
  //
  //
  //

  //
  //#ifdef REFINE_SNODE
  //  //EXPERIMENTAL STUFF
  //  typedef std::set<Int> nodeset;
  //  typedef std::list<nodeset*> partitions;
  //  typedef SYMPACK::vector<nodeset*> vecset;
  //  void SparseMatrixStructure::RefineSupernodes(ETree& tree, Ordering & aOrder, SYMPACK::vector<Int> & supMembership, SYMPACK::vector<Int> & xsuper, PtrVec & xlindx, IdxVec & lindx, SYMPACK::vector<Int> & perm){
  //
  //    perm.resize(size);
  //
  //    Int nsuper = xsuper.size()-1;
  //    //Using std datatypes first
  //
  //    vecset adj(nsuper);
  //    vecset snodes(nsuper);
  //
  //    partitions L;
  //
  //    SYMPACK::vector<Int> origPerm(size);
  //
  //    //init L with curent supernodal partition
  //    Int pos = 1;
  //    for(Int i = 1; i<=nsuper; ++i){
  //      adj[i-1] = new nodeset();
  //      Ptr fi = xlindx[i-1];
  //      Ptr li = xlindx[i]-1;
  //      for(Ptr idx = fi; idx<=li;idx++){
  //        adj[i-1]->insert(lindx[idx-1]);
  //      }
  //
  //      L.push_back(new nodeset());
  //      snodes[i-1] = L.back();
  //      nodeset * curL = L.back();
  //      Int fc = xsuper[i-1];
  //      Int lc = xsuper[i]-1;
  //      for(Int node = fc; node<=lc;++node){
  //        curL->insert(node);
  //        origPerm[node-1] = pos++;
  //      }
  //    }
  //
  //    Int K = nsuper;
  //    partitions::reverse_iterator Kit;
  //    for(Kit = L.rbegin();Kit != L.rend(); ++Kit){
  //      logfileptr->OFS()<<"Looking at snode "<<K<<std::endl;
  //
  //      assert( snodes[K-1] == *Kit);
  //
  //      //        logfileptr->OFS()<<"Adj is "<<*adj[K-1]<<std::endl;
  //
  //      partitions::reverse_iterator Jit;
  //      //    partitions tmp;
  //      Jit = L.rbegin();
  //      Int count = L.size() - K;
  //      while(count>0){
  //        //        logfileptr->OFS()<<"L is "<<**Jit<<std::endl;
  //        //for(Jit = L.rbegin();Jit!=Kit;++Jit){
  //        nodeset * inter = new nodeset();
  //
  //        std::set_intersection((*Jit)->begin(),(*Jit)->end(),
  //            adj[K-1]->begin(),adj[K-1]->end(),
  //            std::inserter(*inter, inter->begin()));
  //        //        logfileptr->OFS()<<"Intersect is "<<*inter<<std::endl;
  //
  //        if(inter->size()>0){
  //          nodeset * diff = new nodeset();
  //          std::set_difference((*Jit)->begin(),(*Jit)->end(),
  //              adj[K-1]->begin(),adj[K-1]->end(),
  //              std::inserter(*diff, diff->begin()));
  //
  //          //        logfileptr->OFS()<<"Diff is "<<*diff<<std::endl;
  //
  //          if(diff->size()>0){
  //            //          tmp.push_back(diff);
  //            //          tmp.push_back(inter);
  //            //replace Jit by inter and diff
  //            (*Jit)->swap(*diff);
  //            L.insert(Jit.base(),inter);
  //            delete diff;
  //            //          (*Jit)->swap(*inter);
  //            //          L.insert(Jit.base(),diff);
  //            //          delete inter;
  //            ++Jit; 
  //            ++Jit; 
  //          }
  //          else{
  //            //          tmp.push_back(*Jit);
  //            delete diff;
  //          }
  //        }
  //        else{
  //          //        tmp.push_back(*Jit);
  //          delete inter;
  //        }
  //        //}
  //        ++Jit;
  //        --count;
  //      }
  //
  //
  //      partitions::iterator it;
  //      for(it = L.begin();it != L.end(); ++it){
  //        logfileptr->OFS()<<"[ ";
  //        for(nodeset::iterator nit = (*it)->begin();nit != (*it)->end(); ++nit){
  //          logfileptr->OFS()<<*nit<<" ";
  //        }
  //        logfileptr->OFS()<<"] ";
  //      }
  //      logfileptr->OFS()<<std::endl;
  //      --K;
  //    }
  //
  //    partitions::iterator it;
  //    for(it = L.begin();it != L.end(); ++it){
  //      logfileptr->OFS()<<"[ ";
  //      for(nodeset::iterator nit = (*it)->begin();nit != (*it)->end(); ++nit){
  //        logfileptr->OFS()<<*nit<<" ";
  //      }
  //      logfileptr->OFS()<<"] ";
  //    }
  //    logfileptr->OFS()<<std::endl;
  //
  //
  //    //construct perm
  //    pos = 1;
  //    for(it = L.begin();it != L.end(); ++it){
  //      for(nodeset::iterator nit = (*it)->begin();nit != (*it)->end(); ++nit){
  //        perm[*nit-1] = pos++;
  //      }
  //    }
  //
  //    logfileptr->OFS()<<"Orig col order "<<origPerm<<std::endl;
  //    logfileptr->OFS()<<"Refined col order "<<perm<<std::endl;
  //
  //    //  SYMPACK::vector<Int> lindxTemp = lindx;
  //    //  //change lindx to reflect the new ordering
  //    //  for(Int i = 1; i<=nsuper; ++i){
  //    //    Int fi = xlindx(i-1);
  //    //    Int li = xlindx(i)-1;
  //    //    for(Int idx = fi; idx<=li;idx++){
  //    //       lindx[idx-1] = perm[lindxTemp[idx-1]-1];
  //    //    }
  //    //  }
  //    //
  //    //    logfileptr->OFS()<<"Previous lindx "<<lindxTemp<<std::endl;
  //    //    logfileptr->OFS()<<"Refined lindx "<<lindx<<std::endl;
  //
  //
  //    for(it = L.begin();it != L.end(); ++it){
  //      delete (*it);
  //    }
  //
  //    for(Int i = 1; i<=nsuper; ++i){
  //      delete adj[i-1];
  //    }
  //
  //    aOrder.Compose(perm);
  //
  //  }
  //#endif
  //
  //
  //  void SparseMatrixStructure::RelaxSupernodes(ETree& tree, SYMPACK::vector<Int> & cc,SYMPACK::vector<Int> & supMembership, SYMPACK::vector<Int> & xsuper, RelaxationParameters & params ){
  //
  //    Int nsuper = xsuper.size()-1;
  //
  //    DisjointSet sets;
  //    sets.Initialize(nsuper);
  //    SYMPACK::vector<Int> ncols(nsuper);
  //    SYMPACK::vector<Int> zeros(nsuper);
  //    SYMPACK::vector<Int> newCC(nsuper);
  //    for(Int ksup=nsuper;ksup>=1;--ksup){
  //      Int cset = sets.makeSet(ksup);
  //      sets.Root(cset-1)=ksup;
  //
  //      Int fstcol = xsuper[ksup-1];
  //      Int lstcol = xsuper[ksup]-1;
  //      Int width = lstcol - fstcol +1;
  //      Int length = cc[fstcol-1];
  //      ncols[ksup-1] = width;
  //      zeros[ksup-1] = 0;
  //      newCC[ksup-1] = length;
  //    }
  //
  //
  //
  //    for(Int ksup=nsuper;ksup>=1;--ksup){
  //      Int fstcol = xsuper[ksup-1];
  //      Int lstcol = xsuper[ksup]-1;
  //      Int width = ncols[ksup-1];
  //      Int length = cc[fstcol-1];
  //
  //      Int parent_fstcol = tree.PostParent(lstcol-1);
  //      if(parent_fstcol!=0){
  //        Int parent_snode = supMembership[parent_fstcol-1];
  //        Int pset = sets.find(parent_snode);
  //        parent_snode = sets.Root(pset-1);
  //
  //        bool merge = (parent_snode == ksup+1);
  //
  //
  //        if(merge){
  //          Int parent_width = ncols[parent_snode-1];
  //
  //          Int parent_fstcol = xsuper[parent_snode-1];
  //          Int parent_lstcol = xsuper[parent_snode]-1;
  //          Int totzeros = zeros[parent_snode-1];
  //          Int fused_cols = width + parent_width;
  //
  //          merge = false;
  //          if(fused_cols <= params.nrelax0){
  //            merge = true;
  //          }
  //          else if(fused_cols <=params.maxSize){
  //            double child_lnz = cc[fstcol-1];
  //            double parent_lnz = cc[parent_fstcol-1];
  //            double xnewzeros = width * (parent_lnz + width  - child_lnz);
  //
  //            if(xnewzeros == 0){
  //              merge = true;
  //            }
  //            else{
  //              //all these values are the values corresponding to the merged snode
  //              double xtotzeros = (double)totzeros + xnewzeros;
  //              double xfused_cols = (double) fused_cols;
  //              //new number of nz
  //              double xtotsize = (xfused_cols * (xfused_cols+1)/2) + xfused_cols * (parent_lnz - parent_width);
  //              //percentage of explicit zeros
  //              double z = xtotzeros / xtotsize;
  //
  //              Int totsize = (fused_cols * (fused_cols+1)/2) + fused_cols * ((Int)parent_lnz - parent_width);
  //              totzeros += (Int)xnewzeros;
  //
  //              merge = ((fused_cols <= params.nrelax1 && z < params.zrelax0) 
  //                  || (fused_cols <= params.nrelax2 && z < params.zrelax1)
  //                  || (z<params.zrelax2)) &&
  //                (xtotsize < std::numeric_limits<Int>::max() / sizeof(double));
  //            }
  //
  //          }
  //
  //          if(merge){
  //            //              std::cout<<"merge "<<ksup<<" and "<<parent_snode<<std::endl;
  //            ncols[ksup-1] += ncols[parent_snode-1]; 
  //            zeros[ksup-1] = totzeros;
  //            newCC[ksup-1] = width + newCC[parent_snode-1];
  //            sets.Union(ksup,parent_snode,ksup);
  //          }
  //        } 
  //
  //      }
  //    }
  //
  //    SYMPACK::vector<Int> relXSuper(nsuper+1);
  //    Int nrSuper = 0;
  //    for(Int ksup=1;ksup<=nsuper;++ksup){
  //      Int kset = sets.find(ksup);
  //      if(ksup == sets.Root(kset-1)){
  //        Int fstcol = xsuper[ksup-1];
  //        relXSuper[nrSuper] = fstcol;
  //        newCC[nrSuper] = newCC[ksup-1];
  //        ++nrSuper;
  //      }
  //    }
  //    relXSuper[nrSuper] = xsuper[nsuper];
  //    relXSuper.resize(nrSuper+1);
  //
  //    for(Int ksup=1;ksup<=nrSuper;++ksup){
  //      Int fstcol = relXSuper[ksup-1];
  //      Int lstcol = relXSuper[ksup]-1;
  //      for(Int col = fstcol; col<=lstcol;++col){
  //        supMembership[col-1] = ksup;
  //        cc[col-1] = newCC[ksup-1] + col-fstcol;
  //      }
  //    }
  //
  //    xsuper = relXSuper;
  //    ///      //adjust the column counts
  //    ///      for(Int col=i-2;col>=i-supsize;--col){
  //    ///        cc[col-1] = cc[col]+1;
  //    ///      }
  //
  //
  //  }
  //
  //  void SparseMatrixStructure::SymbolicFactorizationRelaxed(ETree& tree,Ordering & aOrder, const SYMPACK::vector<Int> & cc,const SYMPACK::vector<Int> & xsuper,const SYMPACK::vector<Int> & SupMembership, PtrVec & xlindx, IdxVec & lindx){
  //    TIMER_START(SymbolicFactorization);
  //
  //
  //    if(!bIsGlobal){
  //      throw std::logic_error( "SparseMatrixStructure must be global in order to call SymbolicFactorization\n" );
  //    }
  //
  //    Int nsuper = xsuper.size()-1;
  //
  //
  //
  //
  //    Ptr nzbeg = 0;
  //    //nzend points to the last used slot in lindx
  //    Ptr nzend = 0;
  //
  //    //tail is the end of list indicator (in rchlnk, not mrglnk)
  //    Idx tail = size +1;
  //
  //    Idx head = 0;
  //
  //    //Array of length nsuper containing the children of 
  //    //each supernode as a linked list
  //    SYMPACK::vector<Idx> mrglnk(nsuper,0);
  //
  //    //Array of length n+1 containing the current linked list 
  //    //of merged indices (the "reach" set)
  //    SYMPACK::vector<Idx> rchlnk(size+1);
  //
  //    //Array of length n used to mark indices as they are introduced
  //    // into each supernode's index set
  //    SYMPACK::vector<Int> marker(size,0);
  //
  //
  //
  //    xlindx.resize(nsuper+1);
  //
  //    //Compute the sum of the column count and resize lindx accordingly
  //    Ptr nofsub = 1;
  //    for(Int i =0; i<cc.size();++i){
  //      nofsub+=cc[i];
  //    }
  //
  //    lindx.resize(nofsub);
  //
  //
  //    Ptr point = 1;
  //    for(Int ksup = 1; ksup<=nsuper; ++ksup){
  //      Int fstcol = xsuper[ksup-1];
  //      xlindx[ksup-1] = point;
  //      point += cc[fstcol-1]; 
  //    } 
  //    xlindx[nsuper] = point;
  //
  //
  //
  //    for(Int ksup = 1; ksup<=nsuper; ++ksup){
  //      Int fstcol = xsuper[ksup-1];
  //      Int lstcol = xsuper[ksup]-1;
  //      Int width = lstcol - fstcol +1;
  //      Int length = cc[fstcol-1];
  //      Ptr knz = 0;
  //      rchlnk[head] = tail;
  //      Int jsup = mrglnk[ksup-1];
  //
  //      //If ksup has children in the supernodal e-tree
  //      if(jsup>0){
  //        //copy the indices of the first child jsup into 
  //        //the linked list, and mark each with the value 
  //        //ksup.
  //        Int jwidth = xsuper[jsup]-xsuper[jsup-1];
  //        Ptr jnzbeg = xlindx[jsup-1] + jwidth;
  //        Ptr jnzend = xlindx[jsup] -1;
  //        for(Ptr jptr = jnzend; jptr>=jnzbeg; --jptr){
  //          Idx newi = lindx[jptr-1];
  //          ++knz;
  //          marker[newi-1] = ksup;
  //          rchlnk[newi] = rchlnk[head];
  //          rchlnk[head] = newi;
  //        }
  //
  //        //for each subsequent child jsup of ksup ...
  //        jsup = mrglnk[jsup-1];
  //        while(jsup!=0 && knz < length){
  //          //merge the indices of jsup into the list,
  //          //and mark new indices with value ksup.
  //
  //          jwidth = xsuper[jsup]-xsuper[jsup-1];
  //          jnzbeg = xlindx[jsup-1] + jwidth;
  //          jnzend = xlindx[jsup] -1;
  //          Int nexti = head;
  //          for(Ptr jptr = jnzbeg; jptr<=jnzend; ++jptr){
  //            Idx newi = lindx[jptr-1];
  //            Idx i;
  //            do{
  //              i = nexti;
  //              nexti = rchlnk[i];
  //            }while(newi > nexti);
  //
  //            if(newi < nexti){
  //#ifdef _DEBUG_
  //              logfileptr->OFS()<<jsup<<" is a child of "<<ksup<<" and "<<newi<<" is inserted in the structure of "<<ksup<<std::endl;
  //#endif
  //              ++knz;
  //              rchlnk[i] = newi;
  //              rchlnk[newi] = nexti;
  //              marker[newi-1] = ksup;
  //              nexti = newi;
  //            }
  //          }
  //          jsup = mrglnk[jsup-1];
  //        }
  //      }
  //
  //      //structure of a(*,fstcol) has not been examined yet.  
  //      //"sort" its structure into the linked list,
  //      //inserting only those indices not already in the
  //      //list.
  //      if(knz < length){
  //        for(Int row = fstcol; row<=lstcol; ++row){
  //          Idx newi = row;
  //          if(newi > fstcol && marker[newi-1] != ksup){
  //            //position and insert newi in list and
  //            // mark it with kcol
  //            Idx nexti = head;
  //            Idx i;
  //            do{
  //              i = nexti;
  //              nexti = rchlnk[i];
  //            }while(newi > nexti);
  //            ++knz;
  //            rchlnk[i] = newi;
  //            rchlnk[newi] = nexti;
  //            marker[newi-1] = ksup;
  //          }
  //        }
  //
  //
  //        for(Int col = fstcol; col<=lstcol; ++col){
  //          Int node = aOrder.perm[col-1];
  //
  //          Ptr knzbeg = expColptr[node-1];
  //          Ptr knzend = expColptr[node]-1;
  //          for(Ptr kptr = knzbeg; kptr<=knzend;++kptr){
  //            Idx newi = expRowind[kptr-1];
  //            newi = aOrder.invp[newi-1];
  //
  //            if(newi > fstcol && marker[newi-1] != ksup){
  //              //position and insert newi in list and
  //              // mark it with kcol
  //              Idx nexti = head;
  //              Idx i;
  //              do{
  //                i = nexti;
  //                nexti = rchlnk[i];
  //              }while(newi > nexti);
  //              ++knz;
  //              rchlnk[i] = newi;
  //              rchlnk[newi] = nexti;
  //              marker[newi-1] = ksup;
  //            }
  //          }
  //        }
  //
  //      } 
  //
  //      //if ksup has no children, insert fstcol into the linked list.
  //      if(rchlnk[head] != fstcol){
  //        rchlnk[fstcol] = rchlnk[head];
  //        rchlnk[head] = fstcol;
  //        ++knz;
  //      }
  //
  //      assert(knz == cc[fstcol-1]);
  //
  //
  //      //copy indices from linked list into lindx(*).
  //      nzbeg = nzend+1;
  //      nzend += knz;
  //      assert(nzend+1 == xlindx[ksup]);
  //      Idx i = head;
  //      for(Ptr kptr = nzbeg; kptr<=nzend;++kptr){
  //        i = rchlnk[i];
  //        lindx[kptr-1] = i;
  //      } 
  //
  //      //if ksup has a parent, insert ksup into its parent's 
  //      //"merge" list.
  //      if(length > width){
  //        Idx pcol = lindx[xlindx[ksup-1] + width -1];
  //        Int psup = SupMembership[pcol-1];
  //        mrglnk[ksup-1] = mrglnk[psup-1];
  //        mrglnk[psup-1] = ksup;
  //      }
  //    }
  //
  //    lindx.resize(nzend+1);
  //
  //    TIMER_STOP(SymbolicFactorization);
  //  }
  //
  //  void SparseMatrixStructure::SymbolicFactorization(ETree& tree, Ordering & aOrder, const SYMPACK::vector<Int> & cc,const SYMPACK::vector<Int> & xsuper,const SYMPACK::vector<Int> & SupMembership, PtrVec & xlindx, IdxVec & lindx){
  //    TIMER_START(SymbolicFactorization);
  //
  //
  //    if(!bIsGlobal){
  //      throw std::logic_error( "SparseMatrixStructure must be global in order to call SymbolicFactorization\n" );
  //    }
  //
  //    Int nsuper = xsuper.size()-1;
  //
  //    Ptr nzbeg = 0;
  //    //nzend points to the last used slot in lindx
  //    Ptr nzend = 0;
  //
  //    //tail is the end of list indicator (in rchlnk, not mrglnk)
  //    Idx tail = size +1;
  //
  //    Idx head = 0;
  //
  //    //Array of length nsuper containing the children of 
  //    //each supernode as a linked list
  //    SYMPACK::vector<Idx> mrglnk(nsuper,0);
  //
  //    //Array of length n+1 containing the current linked list 
  //    //of merged indices (the "reach" set)
  //    SYMPACK::vector<Idx> rchlnk(size+1);
  //
  //    //Array of length n used to mark indices as they are introduced
  //    // into each supernode's index set
  //    SYMPACK::vector<Int> marker(size,0);
  //
  //
  //
  //    xlindx.resize(nsuper+1);
  //
  //    //Compute the sum of the column count and resize lindx accordingly
  //    Ptr nofsub = 1;
  //    for(Int i =0; i<cc.size();++i){
  //      nofsub+=cc[i];
  //    }
  //
  //    lindx.resize(nofsub);
  //
  //
  //    Ptr point = 1;
  //    for(Int ksup = 1; ksup<=nsuper; ++ksup){
  //      Int fstcol = xsuper[ksup-1];
  //      xlindx[ksup-1] = point;
  //      point += cc[fstcol-1]; 
  //    } 
  //    xlindx[nsuper] = point;
  //
  //    for(Int ksup = 1; ksup<=nsuper; ++ksup){
  //      Int fstcol = xsuper[ksup-1];
  //      Int lstcol = xsuper[ksup]-1;
  //      Int width = lstcol - fstcol +1;
  //      Int length = cc[fstcol-1];
  //      Idx knz = 0;
  //      rchlnk[head] = tail;
  //      Int jsup = mrglnk[ksup-1];
  //
  //      //If ksup has children in the supernodal e-tree
  //      if(jsup>0){
  //        //copy the indices of the first child jsup into 
  //        //the linked list, and mark each with the value 
  //        //ksup.
  //        Int jwidth = xsuper[jsup]-xsuper[jsup-1];
  //        Ptr jnzbeg = xlindx[jsup-1] + jwidth;
  //        Ptr jnzend = xlindx[jsup] -1;
  //        for(Ptr jptr = jnzend; jptr>=jnzbeg; --jptr){
  //          Idx newi = lindx[jptr-1];
  //          ++knz;
  //          marker[newi-1] = ksup;
  //          rchlnk[newi] = rchlnk[head];
  //          rchlnk[head] = newi;
  //        }
  //
  //        //for each subsequent child jsup of ksup ...
  //        jsup = mrglnk[jsup-1];
  //        while(jsup!=0 && knz < length){
  //          //merge the indices of jsup into the list,
  //          //and mark new indices with value ksup.
  //
  //          jwidth = xsuper[jsup]-xsuper[jsup-1];
  //          jnzbeg = xlindx[jsup-1] + jwidth;
  //          jnzend = xlindx[jsup] -1;
  //          Idx nexti = head;
  //          for(Ptr jptr = jnzbeg; jptr<=jnzend; ++jptr){
  //            Idx newi = lindx[jptr-1];
  //            Idx i;
  //            do{
  //              i = nexti;
  //              nexti = rchlnk[i];
  //            }while(newi > nexti);
  //
  //            if(newi < nexti){
  //#ifdef _DEBUG_
  //              logfileptr->OFS()<<jsup<<" is a child of "<<ksup<<" and "<<newi<<" is inserted in the structure of "<<ksup<<std::endl;
  //#endif
  //              ++knz;
  //              rchlnk[i] = newi;
  //              rchlnk[newi] = nexti;
  //              marker[newi-1] = ksup;
  //              nexti = newi;
  //            }
  //          }
  //          jsup = mrglnk[jsup-1];
  //        }
  //      }
  //
  //      //structure of a(*,fstcol) has not been examined yet.  
  //      //"sort" its structure into the linked list,
  //      //inserting only those indices not already in the
  //      //list.
  //      if(knz < length){
  //        Idx node = aOrder.perm[fstcol-1];
  //        Ptr knzbeg = expColptr[node-1];
  //        Ptr knzend = expColptr[node]-1;
  //        for(Ptr kptr = knzbeg; kptr<=knzend;++kptr){
  //          Idx newi = expRowind[kptr-1];
  //          newi = aOrder.invp[newi-1];
  //          if(newi > fstcol && marker[newi-1] != ksup){
  //            //position and insert newi in list and
  //            // mark it with kcol
  //            Idx nexti = head;
  //            Idx i;
  //            do{
  //              i = nexti;
  //              nexti = rchlnk[i];
  //            }while(newi > nexti);
  //            ++knz;
  //            rchlnk[i] = newi;
  //            rchlnk[newi] = nexti;
  //            marker[newi-1] = ksup;
  //          }
  //        }
  //      }
  //
  //      //if ksup has no children, insert fstcol into the linked list.
  //      if(rchlnk[head] != fstcol){
  //        rchlnk[fstcol] = rchlnk[head];
  //        rchlnk[head] = fstcol;
  //        ++knz;
  //      }
  //
  //      assert(knz == cc[fstcol-1]);
  //
  //
  //      //copy indices from linked list into lindx(*).
  //      nzbeg = nzend+1;
  //      nzend += knz;
  //      assert(nzend+1 == xlindx[ksup]);
  //      Int i = head;
  //      for(Int kptr = nzbeg; kptr<=nzend;++kptr){
  //        i = rchlnk[i];
  //        lindx[kptr-1] = i;
  //      } 
  //
  //      //if ksup has a parent, insert ksup into its parent's 
  //      //"merge" list.
  //      if(length > width){
  //        Int pcol = lindx[xlindx[ksup-1] + width -1];
  //        Int psup = SupMembership[pcol-1];
  //        mrglnk[ksup-1] = mrglnk[psup-1];
  //        mrglnk[psup-1] = ksup;
  //      }
  //    }
  //
  //    lindx.resize(nzend+1);
  //
  //    TIMER_STOP(SymbolicFactorization);
  //  }
}
