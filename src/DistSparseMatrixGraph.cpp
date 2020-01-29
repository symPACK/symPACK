#include "sympack/DistSparseMatrixGraph.hpp"
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

namespace symPACK{


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
        std::vector<Idx> newRowind(rowind.size()+VertexCount());
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
      bassert(colbeg<=colend);
      sort(&rowind[0]+colbeg,&rowind[0]+colend,std::less<Idx>());
    }
  }


  void SparseMatrixGraph::DistributeGraph(DistSparseMatrixGraph & dg){
    throw std::logic_error( "SparseMatrixGraph::DistributeGraph(DistSparseMatrixGraph & ) not implemented\n" );
  }

  void SparseMatrixGraph::BroadcastGraph(MPI_Comm comm, Int root){

    MPI_Datatype graphType;
    int blocklen[8];
    MPI_Aint disps[8];

    MPI_Datatype Idxtype;
    MPI_Type_contiguous( sizeof(Idx), MPI_BYTE, &Idxtype );
    MPI_Type_commit(&Idxtype);
    MPI_Datatype Ptrtype;
    MPI_Type_contiguous( sizeof(Ptr), MPI_BYTE, &Ptrtype );
    MPI_Type_commit(&Ptrtype);

    MPI_Datatype types[8] = {MPI_BYTE,MPI_BYTE,MPI_BYTE,MPI_BYTE,MPI_BYTE,
      MPI_BYTE,Ptrtype,Idxtype};
    blocklen[0] = sizeof(size);
    blocklen[1] = sizeof(nnz);
    blocklen[2] = sizeof(baseval);
    blocklen[3] = sizeof(keepDiag);
    blocklen[4] = sizeof(sorted);
    blocklen[5] = sizeof(expanded);
    blocklen[6] = colptr.size();
    blocklen[7] = rowind.size();

    MPI_Get_address( (void *)&size,  &disps[0]);
    MPI_Get_address( (void *)&nnz,  &disps[1]);
    MPI_Get_address( (void *)&baseval,  &disps[2]);
    MPI_Get_address( (void *)&keepDiag,  &disps[3]);
    MPI_Get_address( (void *)&sorted,  &disps[4]);
    MPI_Get_address( (void *)&expanded,  &disps[5]);
    MPI_Get_address( (void *)&colptr[0],  &disps[6]);
    MPI_Get_address( (void *)&rowind[0],  &disps[7]);

    MPI_Type_create_struct(8, blocklen, disps, types, &graphType);
    MPI_Type_commit(&graphType);


    MPI_Bcast(MPI_BOTTOM,1,graphType,root,comm); 

    MPI_Type_free(&graphType);
    MPI_Type_free(&Ptrtype);
    MPI_Type_free(&Idxtype);
  }

  void SparseMatrixGraph::ExpandSymmetric(){



    SYMPACK_TIMER_START(EXPAND);
    if(!expanded){
      int mpisize = 1;
      int mpirank = 0;

      Idx N = size; 


      Idx firstLocCol = LocalFirstVertex()-baseval;
      Idx locColCnt = LocalVertexCount();

      std::vector<duet > Isend,Irecv;
      std::vector<int> ssizes(mpisize,0);
      std::vector<int> rsizes(mpisize,0);
      std::vector<int> sdispls(mpisize+1,0);
      std::vector<int> rdispls(mpisize+1,0);

      std::vector<Ptr> curPos;
      std::vector<Ptr> prevPos;

      //loop through my local columns and figure out the extra nonzero per col on prow
      SYMPACK_TIMER_START(DistMat_Expand_count);
      for(Idx locCol = 0 ; locCol< locColCnt; locCol++){
        Idx col = firstLocCol + locCol;  // 0 based
        Ptr colbeg = colptr[locCol]-baseval; //now 0 based
        Ptr colend = colptr[locCol+1]-baseval; // now 0 based
        for(Ptr rptr = colbeg ; rptr< colend ; rptr++ ){
          Idx row = rowind[rptr]-baseval; //0 based
          if(row>col){
            //where should it go ?
            Idx pdest = 0;
            ssizes[pdest]++;
          }
        }
      }
      SYMPACK_TIMER_STOP(DistMat_Expand_count);


      sdispls[0] = 0;
      std::partial_sum(ssizes.begin(),ssizes.end(),sdispls.begin()+1);
      int totSend = sdispls.back();
      Isend.resize(totSend);
      //serialize
      SYMPACK_TIMER_START(DistMat_Expand_pack);
      for(Idx locCol = 0 ; locCol< locColCnt; locCol++){
        Idx col = firstLocCol + locCol;  // 0 based
        Ptr colbeg = colptr[locCol]-baseval; //now 0 based
        Ptr colend = colptr[locCol+1]-baseval; // now 0 based
        for(Ptr rptr = colbeg ; rptr< colend ; rptr++ ){
          Idx row = rowind[rptr]-baseval; //0 based
          if(row>col){
            Idx pdest = 0;
            duet & trip = Isend[sdispls[pdest]++];
            trip.row = col;
            trip.col = row;
          }
        }
      }  
      SYMPACK_TIMER_STOP(DistMat_Expand_pack);

      sdispls[0] = 0;
      std::partial_sum(ssizes.begin(),ssizes.end(),sdispls.begin()+1);
      totSend = sdispls.back();

      SYMPACK_TIMER_START(DistMat_Expand_communication);
      rsizes = ssizes;
      rdispls = sdispls;


      rdispls[0] = 0;
      std::partial_sum(rsizes.begin(),rsizes.end(),rdispls.begin()+1);
      int totRecv = rdispls.back();
      Irecv.resize(totRecv);

      std::copy(&Isend[0],&Isend[0] + totRecv, &Irecv[0]);
      SYMPACK_TIMER_STOP(DistMat_Expand_communication);

      //now parse

      std::vector<Ptr> newColptr(colptr.size(),0);
      for(int col=colptr.size()-1;col>0;col--){
        //convert to count instead
        newColptr[col] = colptr[col] - colptr[col-1];//baseval-based
      }

      //update the column counts
      for(auto it = Irecv.begin(); it!=Irecv.end(); it++){
        //triplets are 0-based but everything is currently shifted by one
        newColptr[it->col-firstLocCol+1]++;
      }

      //turn it back to colptr
      newColptr[0]=baseval;
      std::partial_sum(newColptr.begin(),newColptr.end(),newColptr.begin());


      Ptr newNNZ = newColptr.back()-baseval;
      rowind.resize(newNNZ);
      this->nnz = nnz*2-size;// newNNZ;


      SYMPACK_TIMER_START(DistMat_Expand_shift);
      //shift the content
      for(int col=colptr.size()-2;col>=0;col--){
        std::copy_backward(&rowind[colptr[col]-baseval],&rowind[colptr[col+1]-baseval],&rowind[newColptr[col+1]-baseval]);
      }
      SYMPACK_TIMER_STOP(DistMat_Expand_shift);

      //add the new content, using colptr as a position backup
      //turn colptr to count
      for(int col=colptr.size()-1;col>0;col--){
        colptr[col] = colptr[col] - colptr[col-1];//baseval-based
      }

      SYMPACK_TIMER_START(DistMat_Expand_unpack);
      for(auto it = Irecv.begin(); it!=Irecv.end(); it++){
        duet & trip = *it;
        Idx locCol = trip.col-firstLocCol;//0-based
        //colptr contains column count shifted by one
        Ptr pos = newColptr[locCol + 1] - colptr[locCol+1]-1;//baseval-based
        colptr[locCol+1]++;

        rowind[pos-baseval] = trip.row+baseval;
      }
      //exchange content for the colptr array
      colptr.swap(newColptr);
      SYMPACK_TIMER_STOP(DistMat_Expand_unpack);

      expanded = 1;
      if(!sorted){
        SetSorted(1);
      }

    }
    SYMPACK_TIMER_STOP(EXPAND);
  }








  DistSparseMatrixGraph::DistSparseMatrixGraph(){
    expanded = 0;
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
    expanded = g.expanded;
    size = g.size;
    nnz = g.nnz;
    colptr = g.colptr;
    rowind = g.rowind;
    SetComm(g.comm);
    vertexDist = g.vertexDist;
    baseval = g.baseval;
    keepDiag = g.keepDiag;
    sorted = g.sorted;
    return *this;
  }

  DistSparseMatrixGraph::DistSparseMatrixGraph( const DistSparseMatrixGraph& g ):DistSparseMatrixGraph(){
    (*this) = g;
  }

  void DistSparseMatrixGraph::SetComm(const MPI_Comm & aComm){
    if(aComm!=comm){
      bool doDup = true;
      int prevSize = -1;
      if(comm!= MPI_COMM_NULL){
        int isSame = 0;
        MPI_Comm_compare(comm, aComm, &isSame);
        if(isSame==MPI_IDENT){
          doDup=false;
        }
        else{
          MPI_Comm_size(comm,&prevSize);
          MPI_Comm_free(&comm);
          //what do we do with vertex dist ?
        }
      }
      if(doDup){
        MPI_Comm_dup(aComm,&comm);
        MPI_Comm_size(comm,&mpisize);
        MPI_Comm_rank(comm,&mpirank);

        if(prevSize!=mpisize){
          //recompute a "balanced" distribution of vertices
          Int colPerProc = size / mpisize;
          vertexDist.assign(mpisize+1,colPerProc);
          vertexDist[0] = baseval;
          std::partial_sum(vertexDist.begin(),vertexDist.end(),vertexDist.begin());
          vertexDist.back() = size + baseval;
          //redistribute the graph ?
        }
      }
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

  void DistSparseMatrixGraph::SetExpanded(int aExpanded) {
    if(expanded != aExpanded){
      if(aExpanded){
        //expand to unsymmetric storage
        ExpandSymmetric();
      }
      else{
        //restrict to symmetric storage
        ToSymmetric();
      }
    }
  }



  void DistSparseMatrixGraph::SetKeepDiag(int aKeepDiag){
    if(aKeepDiag != keepDiag){
      if(aKeepDiag){
        //add the diagonal entry
        std::vector<Idx> newRowind(rowind.size()+LocalVertexCount());
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




  void DistSparseMatrixGraph::SortEdges(){
    for(Idx locCol = 0 ; locCol< LocalVertexCount(); locCol++){
      Ptr colbeg = colptr[locCol]-baseval; //now 0 based
      Ptr colend = colptr[locCol+1]-baseval; // now 0 based 
      bassert(colbeg<=colend);
      sort(&rowind[0]+colbeg,&rowind[0]+colend,std::less<Idx>());
    }
  }

  void DistSparseMatrixGraph::Permute(Int * invp){
    permute_(invp, nullptr, 1);
  }
  void DistSparseMatrixGraph::Permute(Int * invp, Idx * newVertexDist){
    permute_(invp, newVertexDist, 1);
  }

  void DistSparseMatrixGraph::Redistribute(Idx * newVertexDist){
    permute_(nullptr, newVertexDist, 1);
  }


  void DistSparseMatrixGraph::Permute(Int * invp, Int invpbaseval){
    permute_(invp, &vertexDist[0], invpbaseval);
  }
  void DistSparseMatrixGraph::Permute(Int * invp, Idx * newVertexDist, Int invpbaseval){
    permute_(invp, newVertexDist, invpbaseval);
  }

  void DistSparseMatrixGraph::permute_(Int * invp, Idx * newVertexDist, Int invpbaseval){
    SYMPACK_TIMER_START(PERMUTE);

    //handle default parameter values
    if(newVertexDist==nullptr){
      newVertexDist = &vertexDist[0];
    }

    Int firstCol = LocalFirstVertex()-baseval; //0 based
    Int newFirstCol = newVertexDist[mpirank]-baseval; //0 based
    Ptr newVtxCount = newVertexDist[mpirank+1] - newVertexDist[mpirank];

    std::vector<int> sizes(mpisize,0);

    for(Idx locCol = 0; locCol<LocalVertexCount(); locCol++){
      Idx col = firstCol + locCol; //0-based
      Ptr colbeg = colptr[locCol] - baseval;
      Ptr colend = colptr[locCol+1] - baseval;
      Idx permCol = invp!=nullptr?invp[col]-invpbaseval:col; // 0 based;
      //find destination processors
      Idx pdest; for(pdest = 0; pdest<mpisize; pdest++){ if(permCol>=newVertexDist[pdest]-baseval && permCol < newVertexDist[pdest+1]-baseval){ break;} }

      //now permute rows
      for(Ptr jptr = colbeg; jptr<colend; jptr++){
        Idx row = rowind[jptr] - baseval; //0 based
        Idx permRow = invp!=nullptr?invp[row] - invpbaseval:row; // 0 based

        if(permRow<permCol && !expanded){
          Idx pdestR; for(pdestR = 0; pdestR<mpisize; pdestR++){ if(permRow>=newVertexDist[pdestR]-baseval && permRow < newVertexDist[pdestR+1]-baseval){ break;} }
          sizes[pdestR]++;
        }
        else if(permRow>permCol || (permRow==permCol && keepDiag) || expanded){
          sizes[pdest]++;
        }

      }
    }



    //First allgatherv to get the receive sizes
    std::vector<int> displs(mpisize+1,0);
    displs[0] = 0;
    std::partial_sum(sizes.begin(),sizes.end(),displs.begin()+1);
    Int totSend = displs.back();
    std::vector< std::pair<Idx,Idx> > sbuf(totSend);

    //pack
    for(Idx locCol = 0; locCol<LocalVertexCount(); locCol++){
      Idx col = firstCol + locCol; 
      Ptr colbeg = colptr[locCol]-baseval;
      Ptr colend = colptr[locCol+1]-baseval;
      Idx permCol = invp!=nullptr?invp[col]-invpbaseval:col; // 0 based;
      //find destination processors
      Idx pdest; for(pdest = 0; pdest<mpisize; pdest++){ if(permCol>=newVertexDist[pdest]-baseval && permCol < newVertexDist[pdest+1]-baseval){ break;} }

      for(Ptr jptr = colbeg; jptr<colend; jptr++){
        Idx row = rowind[jptr] - baseval; //0 based
        Idx permRow = invp!=nullptr?invp[row] - invpbaseval:row; // 0 based
        if(permRow<permCol && !expanded){
          Idx pdestR; for(pdestR = 0; pdestR<mpisize; pdestR++){ if(permRow>=newVertexDist[pdestR]-baseval && permRow < newVertexDist[pdestR+1]-baseval){ break;} }
          assert(displs[pdestR]<sbuf.size());
          sbuf[displs[pdestR]++] = std::make_pair(permRow,permCol);
        }
        else if(permRow>permCol || (permRow==permCol && keepDiag) || expanded){
          assert(displs[pdest]<sbuf.size());
          sbuf[displs[pdest]++] = std::make_pair(permCol,permRow);
        }
      }
    }

    MPI_Datatype type;
    MPI_Type_contiguous( sizeof( std::pair<Idx,Idx> ), MPI_BYTE, &type );
    MPI_Type_commit(&type);




    //re compute send displs in bytes
    displs[0] = 0;
    std::partial_sum(sizes.begin(),sizes.end(),&displs[1]);


    std::vector<int> rsizes(mpisize,0);
    std::vector<int> rdispls(mpisize+1,0);
    MPI_Alltoall(sizes.data(),sizeof(int),MPI_BYTE,rsizes.data(),sizeof(int),MPI_BYTE,comm);

    //now compute receiv displs with actual sizes 
    rdispls[0] = 0;
    std::partial_sum(rsizes.begin(),rsizes.end(),&rdispls[1]);
    Ptr totRecv = rdispls.back();
    std::vector< std::pair<Idx,Idx> > rbuf(totRecv);
    MPI_Alltoallv(&sbuf[0],&sizes[0],&displs[0],type,&rbuf[0],&rsizes[0],&rdispls[0],type,comm);

    MPI_Type_free(&type);

    std::vector<Ptr> newColptr(newVtxCount+1,0);
    //compute column counts
    for(auto it = rbuf.begin(); it!=rbuf.end(); it++){
      //pairs are 0-based but everything is currently shifted by one
      newColptr[it->first-newFirstCol+1]++;
    }
    //turn it back to colptr
    newColptr[0]=baseval;
    std::partial_sum(newColptr.begin(),newColptr.end(),newColptr.begin());

    Ptr newNNZ = newColptr.back()-baseval;
    rowind.resize(newNNZ);

    colptr = newColptr;

    //now use newColptr as a position backup
    for(auto it = rbuf.begin(); it!=rbuf.end(); it++){
      Idx locCol = it->first-newFirstCol;//0-based
      //colptr contains column count shifted by one
      Ptr pos = newColptr[locCol]++;//baseval-based
      rowind[pos-baseval] = it->second+baseval;
    }

    //copy newVertexDist into vertexDist
    std::copy(newVertexDist,newVertexDist+mpisize+1,&vertexDist[0]);


    if(sorted){
      SortEdges();
    }

    SYMPACK_TIMER_STOP(PERMUTE);
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

    Int nsuperLocal = xsuperdist[mpirank+1]-xsuperdist[mpirank];
    Int firstSnode = xsuperdist[mpirank]-1;//0 based
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


    Int numColSentBefore = supFirstColumn>curFirstColumn?std::min(curLastColumn,supFirstColumn) - curFirstColumn +(curLastColumn<supFirstColumn?1:0):0;
    Ptr numRowSentBefore = numColSentBefore>0?colptr[std::min(curLastColumn,supFirstColumn)-firstCol]-colptr[curFirstColumn-firstCol]:0;
    Int numColSentAfter = curLastColumn>supLastColumn?curLastColumn-std::max(curFirstColumn,supLastColumn)+(curFirstColumn>supLastColumn?1:0):0;

    Idx startColAfter = std::max(curFirstColumn,supLastColumn)+(curFirstColumn>supLastColumn?0:1);
    Ptr numRowSentAfter = numColSentAfter>0?colptr[curLastColumn+1-firstCol]-colptr[startColAfter-firstCol]:0;

    //first transfer column pointers
    vector<int> ssizes(mpisize,0);
    vector<int> ssizesR(mpisize,0);
    vector<int> ssizesAfter(mpisize,0);
    vector<int> sdispls(mpisize+1,-1);
    vector<int> sdisplsR(mpisize+1,-1);
    int totSent = 0;
    int totSentR = 0;

    if(numColSentBefore>0){
      for(Idx col = curFirstColumn; col<curFirstColumn+numColSentBefore; col++){
        bassert(col>=curFirstColumn && col<=curLastColumn);
        Int snode = supMembership[col];
        Int pdest = 0;
        for(pdest = 0; pdest<mpisize;pdest++){
          if(xsuperdist[pdest]<=snode && snode<xsuperdist[pdest+1]){
            break;
          }
        }

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
      for(Idx col = startColAfter; col<startColAfter+numColSentAfter; col++){
        bassert(col>=curFirstColumn && col<=curLastColumn);
        Int snode = supMembership[col];
        Int pdest = 0;
        for(pdest = 0; pdest<mpisize;pdest++){
          if(xsuperdist[pdest]<=snode && snode<xsuperdist[pdest+1]){
            break;
          }
        }

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



    std::vector<char> recvBuf(rdispls.back()); 
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

    if(recvCntBefore>0){
      std::copy_backward(rowind.begin(),rowind.begin()+rowindSize,rowind.begin()+recvCntBefore+rowindSize);
    }


    //copy from recvBuf

    //if I have received something from ranks before me
    Idx rowindPos =0;
    Idx idxBufPos = 0;
    Idx * idxBuf = (Idx*)&recvBuf[0];


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

    supLastColumn = supLastColumn + baseval +1;
    MPI_Allgather(&supLastColumn,sizeof(supLastColumn),MPI_BYTE,&vertexDist[1],sizeof(supLastColumn),MPI_BYTE,comm);
  }







  void DistSparseMatrixGraph::ExpandSymmetric(){
    SYMPACK_TIMER_START(EXPAND);
    if(!expanded){

      int ismpi=0;
      MPI_Initialized( &ismpi);
      int isnull= (comm == MPI_COMM_NULL);

      if(!ismpi || isnull){
        //throw an exception
        throw std::logic_error("MPI communicator needs to be set.");
      }


      int mpisize,mpirank;
      MPI_Comm_size(comm,&mpisize);
      MPI_Comm_rank(comm,&mpirank);

      Idx N = size; 


      Idx firstLocCol = LocalFirstVertex()-baseval;
      Idx locColCnt = LocalVertexCount();

      std::vector<duet > Isend,Irecv;
      std::vector<int> ssizes(mpisize,0);
      std::vector<int> rsizes(mpisize,0);
      std::vector<int> sdispls(mpisize+1,0);
      std::vector<int> rdispls(mpisize+1,0);

      std::vector<Ptr> curPos;
      std::vector<Ptr> prevPos;

      //loop through my local columns and figure out the extra nonzero per col on prow
      SYMPACK_TIMER_START(DistMat_Expand_count);
      for(Idx locCol = 0 ; locCol< locColCnt; locCol++){
        Idx col = firstLocCol + locCol;  // 0 based
        Ptr colbeg = colptr[locCol]-baseval; //now 0 based
        Ptr colend = colptr[locCol+1]-baseval; // now 0 based
        for(Ptr rptr = colbeg ; rptr< colend ; rptr++ ){
          Idx row = rowind[rptr]-baseval; //0 based
          if(row>col){
            //where should it go ?
            Idx pdest; for(pdest = 0; pdest<mpisize; pdest++){ if(row>=vertexDist[pdest]-baseval && row < vertexDist[pdest+1]-baseval){ break;} }
            ssizes[pdest]++;
          }
        }
      }
      SYMPACK_TIMER_STOP(DistMat_Expand_count);


      sdispls[0] = 0;
      std::partial_sum(ssizes.begin(),ssizes.end(),sdispls.begin()+1);
      int totSend = sdispls.back();
      Isend.resize(totSend);
      //serialize
      SYMPACK_TIMER_START(DistMat_Expand_pack);
      for(Idx locCol = 0 ; locCol< locColCnt; locCol++){
        Idx col = firstLocCol + locCol;  // 0 based
        Ptr colbeg = colptr[locCol]-baseval; //now 0 based
        Ptr colend = colptr[locCol+1]-baseval; // now 0 based
        for(Ptr rptr = colbeg ; rptr< colend ; rptr++ ){
          Idx row = rowind[rptr]-baseval; //0 based
          if(row>col){
            Idx pdest; for(pdest = 0; pdest<mpisize; pdest++){ if(row>=vertexDist[pdest]-baseval && row < vertexDist[pdest+1]-baseval){ break;} }
            duet & trip = Isend[sdispls[pdest]++];
            trip.row = col;
            trip.col = row;
          }
        }
      }  
      SYMPACK_TIMER_STOP(DistMat_Expand_pack);

      sdispls[0] = 0;
      std::partial_sum(ssizes.begin(),ssizes.end(),sdispls.begin()+1);
      totSend = sdispls.back();

      MPI_Datatype type;
      MPI_Type_contiguous( sizeof(duet), MPI_BYTE, &type );
      MPI_Type_commit(&type);

      SYMPACK_TIMER_START(DistMat_Expand_communication);
      MPI_Alltoall(&ssizes[0],sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,comm);


      rdispls[0] = 0;
      std::partial_sum(rsizes.begin(),rsizes.end(),rdispls.begin()+1);
      int totRecv = rdispls.back();
      Irecv.resize(totRecv);

      MPI_Alltoallv(&Isend[0],&ssizes[0],&sdispls[0],type,&Irecv[0],&rsizes[0],&rdispls[0],type,comm);
      SYMPACK_TIMER_STOP(DistMat_Expand_communication);

      MPI_Type_free(&type);
      //now parse


      std::vector<Ptr> newColptr(colptr.size(),0);
      for(int col=colptr.size()-1;col>0;col--){
        //convert to count instead
        newColptr[col] = colptr[col] - colptr[col-1];//baseval-based
      }

      //update the column counts
      for(auto it = Irecv.begin(); it!=Irecv.end(); it++){
        //triplets are 0-based but everything is currently shifted by one
        newColptr[it->col-firstLocCol+1]++;
      }

      //turn it back to colptr
      newColptr[0]=baseval;
      std::partial_sum(newColptr.begin(),newColptr.end(),newColptr.begin());


      Ptr newNNZ = newColptr.back()-baseval;
      rowind.resize(newNNZ);
      this->nnz = newNNZ;


      SYMPACK_TIMER_START(DistMat_Expand_shift);
      //shift the content
      for(int col=colptr.size()-2;col>=0;col--){
        std::copy_backward(&rowind[colptr[col]-baseval],&rowind[colptr[col+1]-baseval],&rowind[newColptr[col+1]-baseval]);
      }
      SYMPACK_TIMER_STOP(DistMat_Expand_shift);

      //add the new content, using colptr as a position backup
      //turn colptr to count
      for(int col=colptr.size()-1;col>0;col--){
        colptr[col] = colptr[col] - colptr[col-1];//baseval-based
      }

      SYMPACK_TIMER_START(DistMat_Expand_unpack);
      for(auto it = Irecv.begin(); it!=Irecv.end(); it++){
        duet & trip = *it;
        Idx locCol = trip.col-firstLocCol;//0-based
        //colptr contains column count shifted by one
        Ptr pos = newColptr[locCol + 1] - colptr[locCol+1]-1;//baseval-based
        colptr[locCol+1]++;

        rowind[pos-baseval] = trip.row+baseval;
      }
      //exchange content for the colptr array
      colptr.swap(newColptr);
      SYMPACK_TIMER_STOP(DistMat_Expand_unpack);

      expanded = 1;

      if(!sorted){
        SetSorted(1);
      }
    }
    SYMPACK_TIMER_STOP(EXPAND);
  }



  void DistSparseMatrixGraph::ToSymmetric(){
    scope_timer(a,DistSparseMatrixGraph::ToSymmetric);
    if(expanded){
      Idx N = size; 
      Idx firstLocCol = LocalFirstVertex()-baseval;
      Idx locColCnt = LocalVertexCount();

      std::vector<Ptr> nnzcnt(locColCnt+1,0);
      //Loop through local columns and count new nnz per col
      //loop through my local columns and figure out the extra nonzero per col on prow
      SYMPACK_TIMER_START(DistMat_ToSymm_count);
      for(Idx locCol = 0 ; locCol< locColCnt; locCol++){
        Idx col = firstLocCol + locCol;  // 0 based
        Ptr colbeg = colptr[locCol]-baseval; //now 0 based
        Ptr colend = colptr[locCol+1]-baseval; // now 0 based
        for(Ptr rptr = colbeg ; rptr< colend ; rptr++ ){
          Idx row = rowind[rptr]-baseval; //0 based
          if(row>col || (keepDiag && row == col)){
            nnzcnt[locCol+1]++;
          }
        }
      }
      SYMPACK_TIMER_STOP(DistMat_ToSymm_count);

      //turn nnzcnt into the future colptr
      nnzcnt.front() = baseval;
      std::partial_sum(nnzcnt.begin(),nnzcnt.end(),nnzcnt.begin());

      std::vector<Idx> newRowind;
      newRowind.reserve(nnzcnt.back());

      for(Idx locCol = 0 ; locCol< locColCnt; locCol++){
        Idx col = firstLocCol + locCol;  // 0 based
        Ptr colbeg = colptr[locCol]-baseval; //now 0 based
        Ptr colend = colptr[locCol+1]-baseval; // now 0 based
        for(Ptr rptr = colbeg ; rptr< colend ; rptr++ ){
          Idx row = rowind[rptr]-baseval; //0 based
          if(row>=col || (keepDiag && row == col)){
            newRowind.push_back(row+baseval);
          }
        }
      }

      colptr.swap(nnzcnt);
      rowind.swap(newRowind);
      nnz = (nnz - size)/2 + size;
      expanded = 0;
    }
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
    g.expanded = expanded;

    //get other proc vertex counts
    Idx localVertexCnt = LocalVertexCount();
    std::vector<Idx> remoteVertexCnt(mpisize,0);
    MPI_Allgather(&localVertexCnt,sizeof(localVertexCnt),MPI_BYTE,&remoteVertexCnt[0],sizeof(localVertexCnt),MPI_BYTE,comm);
    Idx totalVertexCnt = std::accumulate(remoteVertexCnt.begin(),remoteVertexCnt.end(),0,std::plus<Idx>());
    g.colptr.resize(totalVertexCnt+1);
    //compute receive displacements
    std::vector<int> rsizes(mpisize,0);
    for(int p = 0; p<mpisize;p++){rsizes[p] = (int)remoteVertexCnt[p];}
    std::vector<int> rdispls(mpisize+1,0);
    rdispls[0]=0;
    std::partial_sum(rsizes.begin(),rsizes.end(),rdispls.begin()+1);

    MPI_Datatype Ptrtype;
    MPI_Type_contiguous( sizeof(Ptr), MPI_BYTE, &Ptrtype );
    MPI_Type_commit(&Ptrtype);
    MPI_Allgatherv(&colptr[0],localVertexCnt,Ptrtype,&g.colptr[0],&rsizes[0],&rdispls[0],Ptrtype,comm);
    MPI_Type_free(&Ptrtype);


    Ptr localEdgeCnt = LocalEdgeCount();
    std::vector<Ptr> remoteEdgeCnt(mpisize,0);
    MPI_Allgather(&localEdgeCnt,sizeof(localEdgeCnt),MPI_BYTE,&remoteEdgeCnt[0],sizeof(localEdgeCnt),MPI_BYTE,comm);
    Ptr totalEdgeCnt = std::accumulate(remoteEdgeCnt.begin(),remoteEdgeCnt.end(),0,std::plus<Ptr>());

    g.rowind.resize(totalEdgeCnt);

    //compute receive displacements
    rsizes.assign(mpisize,0);
    for(int p = 0; p<mpisize;p++){rsizes[p] = (int)remoteEdgeCnt[p];}
    rdispls[0]=0;
    std::partial_sum(rsizes.begin(),rsizes.end(),rdispls.begin()+1);
    MPI_Datatype Idxtype;
    MPI_Type_contiguous( sizeof(Idx), MPI_BYTE, &Idxtype );
    MPI_Type_commit(&Idxtype);
    MPI_Allgatherv(&rowind[0],localEdgeCnt,Idxtype,&g.rowind[0],&rsizes[0],&rdispls[0],Idxtype,comm);
    MPI_Type_free(&Idxtype);


    //fix colptr
    Idx pos = remoteVertexCnt[0];
    Ptr offset = 0;
    for(int p=1;p<mpisize;p++){
      offset+=remoteEdgeCnt[p-1]; 
      for(Idx lv = 0; lv < remoteVertexCnt[p]; lv++){
        g.colptr[pos++] += offset;
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


    int iam =0;
    int np =1;
    MPI_Comm_rank(comm,&iam);
    MPI_Comm_size(comm,&np);

    g.size = size;
    g.SetBaseval(baseval);  
    if(iam==proot){
      g.SetSorted(sorted);
      g.SetKeepDiag(keepDiag);
      g.expanded = expanded;
    }

    //get other proc vertex counts
    Idx localVertexCnt = LocalVertexCount();
    std::vector<Idx> remoteVertexCnt(mpisize,0);
    MPI_Allgather(&localVertexCnt,sizeof(localVertexCnt),MPI_BYTE,&remoteVertexCnt[0],sizeof(localVertexCnt),MPI_BYTE,comm);


    MPI_Datatype type;
    MPI_Type_contiguous( sizeof(Ptr), MPI_BYTE, &type );
    MPI_Type_commit(&type);


    if(iam==proot){
      Idx totalVertexCnt = std::accumulate(remoteVertexCnt.begin(),remoteVertexCnt.end(),0,std::plus<Idx>());
      g.colptr.resize(totalVertexCnt+1);


      //compute receive displacements
      std::vector<int> rsizes(mpisize,0);
      for(int p = 0; p<mpisize;p++){rsizes[p] = (int)remoteVertexCnt[p];}
      std::vector<int> rdispls(mpisize+1,0);
      rdispls[0]=0;
      std::partial_sum(rsizes.begin(),rsizes.end(),rdispls.begin()+1);
      MPI_Gatherv(&colptr[0],localVertexCnt,type,&g.colptr[0],&rsizes[0],&rdispls[0],type,proot,comm);
    }
    else{
      MPI_Gatherv(&colptr[0],localVertexCnt,type,nullptr,nullptr,nullptr,type,proot,comm);
    }

    MPI_Type_free(&type);

    Ptr localEdgeCnt = LocalEdgeCount();
    std::vector<Ptr> remoteEdgeCnt(mpisize,0);
    MPI_Gather(&localEdgeCnt,sizeof(localEdgeCnt),MPI_BYTE,&remoteEdgeCnt[0],sizeof(localEdgeCnt),MPI_BYTE,proot,comm);
    Ptr totalEdgeCnt = std::accumulate(remoteEdgeCnt.begin(),remoteEdgeCnt.end(),0,std::plus<Ptr>());

    MPI_Type_contiguous( sizeof(Idx), MPI_BYTE, &type );
    MPI_Type_commit(&type);

    if(iam==proot){
      g.rowind.resize(totalEdgeCnt);
      //compute receive displacements
      std::vector<int> rsizes(mpisize,0);
      std::vector<int> rdispls(mpisize+1,0);
      for(int p = 0; p<mpisize;p++){rsizes[p] = (int)remoteEdgeCnt[p];}
      rdispls[0]=0;
      std::partial_sum(rsizes.begin(),rsizes.end(),rdispls.begin()+1);
      MPI_Gatherv(&rowind[0],localEdgeCnt,type,&g.rowind[0],&rsizes[0],&rdispls[0],type,proot,comm);
    }
    else{
      MPI_Gatherv(&rowind[0],localEdgeCnt,type,nullptr,nullptr,nullptr,type,proot,comm);
    }
    MPI_Type_free(&type);

    if(iam==proot){
      //fix colptr
      Idx pos = remoteVertexCnt[0];
      Ptr offset = 0;
      for(int p=1;p<mpisize;p++){
        offset+=remoteEdgeCnt[p-1]; 
        for(Idx lv = 0; lv < remoteVertexCnt[p]; lv++){
          g.colptr[pos++] += offset;
        }
      }
      if(g.colptr.size()>0){
        g.colptr.back()=totalEdgeCnt + baseval;
      }
    }
  }


}
