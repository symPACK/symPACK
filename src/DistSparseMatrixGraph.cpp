/*
	 Copyright (c) 2016 The Regents of the University of California,
	 through Lawrence Berkeley National Laboratory.  

   Author: Mathias Jacquelin
	 
   This file is part of symPACK. All rights reserved.

	 Redistribution and use in source and binary forms, with or without
	 modification, are permitted provided that the following conditions are met:

	 (1) Redistributions of source code must retain the above copyright notice, this
	 list of conditions and the following disclaimer.
	 (2) Redistributions in binary form must reproduce the above copyright notice,
	 this list of conditions and the following disclaimer in the documentation
	 and/or other materials provided with the distribution.
	 (3) Neither the name of the University of California, Lawrence Berkeley
	 National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
	 be used to endorse or promote products derived from this software without
	 specific prior written permission.

	 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
	 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
	 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
	 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
	 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
	 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
	 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
	 ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
	 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
	 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

	 You are under no obligation whatsoever to provide any bug fixes, patches, or
	 upgrades to the features, functionality or performance of the source code
	 ("Enhancements") to anyone; however, if you choose to make your Enhancements
	 available either publicly, or directly to Lawrence Berkeley National
	 Laboratory, without imposing a separate written license agreement for such
	 Enhancements, then you hereby grant the following license: a non-exclusive,
	 royalty-free perpetual license to install, use, modify, prepare derivative
	 works, incorporate into other computer software, distribute, and sublicense
	 such enhancements or derivative works thereof, in binary and source code form.
*/
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
      sort(&rowind[0]+colbeg,&rowind[0]+colend,std::less<Ptr>());
    }
  }


  void SparseMatrixGraph::DistributeGraph(DistSparseMatrixGraph & dg){
    throw std::logic_error( "SparseMatrixGraph::DistributeGraph(DistSparseMatrixGraph & ) not implemented\n" );
  }

  void SparseMatrixGraph::BroadcastGraph(MPI_Comm comm, Int root){

    MPI_Datatype graphType;
    int blocklen[8];
    MPI_Aint disps[8];

    MPI_Datatype types[8] = {MPI_BYTE,MPI_BYTE,MPI_BYTE,MPI_BYTE,MPI_BYTE,
      MPI_BYTE,MPI_BYTE,MPI_BYTE};
    blocklen[0] = sizeof(size);
    blocklen[1] = sizeof(nnz);
    blocklen[2] = sizeof(baseval);
    blocklen[3] = sizeof(keepDiag);
    blocklen[4] = sizeof(sorted);
    blocklen[5] = sizeof(expanded);
    blocklen[6] = colptr.size()*sizeof(Ptr);
    blocklen[7] = rowind.size()*sizeof(Idx);

    MPI_Address( (void *)&size,  &disps[0]);
    MPI_Address( (void *)&nnz,  &disps[1]);
    MPI_Address( (void *)&baseval,  &disps[2]);
    MPI_Address( (void *)&keepDiag,  &disps[3]);
    MPI_Address( (void *)&sorted,  &disps[4]);
    MPI_Address( (void *)&expanded,  &disps[5]);
    MPI_Address( (void *)&colptr[0],  &disps[6]);
    MPI_Address( (void *)&rowind[0],  &disps[7]);

    MPI_Type_create_struct(8, blocklen, disps, types, &graphType);
    MPI_Type_commit(&graphType);


    MPI_Bcast(MPI_BOTTOM,1,graphType,root,comm); 
      
    MPI_Type_free(&graphType);
  }

  void SparseMatrixGraph::ExpandSymmetric(){
    throw std::logic_error( "SparseMatrixGraph::ExpandSymmetric() not implemented\n" );
  }








  DistSparseMatrixGraph::DistSparseMatrixGraph(){
    expanded = false;
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
if(colbeg>colend){logfileptr->OFS()<<colptr<<std::endl; gdb_lock();}
      sort(&rowind[0]+colbeg,&rowind[0]+colend,std::less<Ptr>());
    }
  }

  void DistSparseMatrixGraph::Permute(Int * invp){
    permute_(invp, NULL, 1);
  }
  void DistSparseMatrixGraph::Permute(Int * invp, Idx * newVertexDist){
    permute_(invp, newVertexDist, 1);
  }

  void DistSparseMatrixGraph::Redistribute(Idx * newVertexDist){
    permute_(NULL, newVertexDist, 1);
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
    if(newVertexDist==NULL){
      newVertexDist = &vertexDist[0];
    }

    //Idx colPerProc = size / mpisize;
    Int firstCol = LocalFirstVertex()-baseval; //0 based
    Int newFirstCol = newVertexDist[mpirank]-baseval; //0 based
    Ptr newVtxCount = newVertexDist[mpirank+1] - newVertexDist[mpirank];


    std::vector<int> sizes(mpisize,0);

    for(Idx locCol = 0; locCol<LocalVertexCount(); locCol++){
      Idx col = firstCol + locCol; //0-based
      Ptr colbeg = colptr[locCol] - baseval;
      Ptr colend = colptr[locCol+1] - baseval;
      //Ptr colbeg = colptr.at(locCol) - baseval;
      //Ptr colend = colptr.at(locCol+1) - baseval;
      Idx permCol = invp!=NULL?invp[col]-invpbaseval:col; // 0 based;
      //find destination processors
      Idx pdest; for(pdest = 0; pdest<mpisize; pdest++){ if(permCol>=newVertexDist[pdest]-baseval && permCol < newVertexDist[pdest+1]-baseval){ break;} }
      //Idx pdest = min( (Idx)mpisize-1, permCol / colPerProc);
      sizes[pdest] += (colend - colbeg)*sizeof(Idx) + sizeof(Ptr) + sizeof(Idx); //extra 2 are for the count of elements and the permuted columns

      //now permute rows
      for(Ptr jptr = colbeg; jptr<colend; jptr++){
        Idx row = rowind[jptr] - baseval; //0 based
        //Idx row = rowind.at(jptr) - baseval; //0 based
        Idx permRow = invp!=NULL?invp[row] - invpbaseval:row; // 0 based
        rowind[jptr] = permRow + baseval;
        //rowind.at(jptr) = permRow + baseval;
      }
    }


    //First allgatherv to get the receive sizes
    std::vector<int> displs(mpisize+1,0);
    std::vector<int> rsizes(mpisize,0);
    std::vector<int> rdispls(mpisize+1,0);

    MPI_Alltoall(&sizes[0],sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,comm);
    //logfileptr->OFS()<<rsizes<<std::endl;
    //    MPI_Alltoallv(&sizes[0],&displs[0],&rdispls[0],MPI_BYTE,&rsizes[0],&displs[0],&rdispls[0],MPI_BYTE,comm);
    //logfileptr->OFS()<<rsizes<<std::endl;



    displs[0] = 0;
    std::partial_sum(sizes.begin(),sizes.end(),displs.begin()+1);
    Int totSend = displs.back();//std::accumulate(sizes.begin(),sizes.end(),0);
    std::vector<char> sbuf(totSend);

    //pack
    for(Idx locCol = 0; locCol<LocalVertexCount(); locCol++){
      Idx col = firstCol + locCol; 
      Ptr colbeg = colptr[locCol]-baseval;
      Ptr colend = colptr[locCol+1]-baseval;
      //Ptr colbeg = colptr.at(locCol) - baseval;
      //Ptr colend = colptr.at(locCol+1) - baseval;
      Idx permCol = invp!=NULL?invp[col]-invpbaseval:col; // perm is 1 based;
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
      std::copy(&rowind[0]+colbeg ,&rowind[0]+colend, pPermRows );
      //std::copy(&rowind.at(colbeg) ,&rowind.at(colend-1)+1, pPermRows );


      //      logfileptr->OFS()<<*pPermCol<<" ("<<locCol<<"): ";
      //      for(Ptr jptr = 0; jptr<*pRowsCnt; jptr++){
      //        logfileptr->OFS()<<pPermRows[jptr]<<" ";
      //      }
      //      logfileptr->OFS()<<std::endl;


    }


    //re compute send displs
    displs[0] = 0;
    std::partial_sum(sizes.begin(),sizes.end(),&displs[1]);

    //now recompute receiv displs with actual sizes 
    rdispls[0] = 0;
    std::partial_sum(rsizes.begin(),rsizes.end(),&rdispls[1]);
    Ptr totRecv = rdispls.back();//std::accumulate(rsizes.begin(),rsizes.end(),0);
    std::vector<char> rbuf(totRecv);

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
      colptr[locCol+1] = *rowsCnt; 
      //colptr.at(locCol+1) = *rowsCnt; 
    }
    std::partial_sum(colptr.begin(),colptr.end(),colptr.begin());

    //now fill rowind
    Ptr nnzLoc = colptr.back()-baseval;
    rowind.resize(nnzLoc);
    std::vector<Ptr> colpos = colptr;
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
      //logfileptr->OFS()<<std::endl;

      std::copy(permRows,permRows + *rowsCnt, &rowind[colpos[locCol]-baseval]);
      colpos[locCol] += *rowsCnt;
      //colpos.at(locCol) += *rowsCnt;
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


    //    logfileptr->OFS()<<"["<<curFirstColumn<<" - "<<curLastColumn<<"] to ["<<supFirstColumn<<" - "<<supLastColumn<<"]"<<std::endl;


    Int numColSentBefore = supFirstColumn>curFirstColumn?std::min(curLastColumn,supFirstColumn) - curFirstColumn +(curLastColumn<supFirstColumn?1:0):0;
    Ptr numRowSentBefore = numColSentBefore>0?colptr[std::min(curLastColumn,supFirstColumn)-firstCol]-colptr[curFirstColumn-firstCol]:0;
    Int numColSentAfter = curLastColumn>supLastColumn?curLastColumn-std::max(curFirstColumn,supLastColumn)+(curFirstColumn>supLastColumn?1:0):0;

    Idx startColAfter = std::max(curFirstColumn,supLastColumn)+(curFirstColumn>supLastColumn?0:1);
    Ptr numRowSentAfter = numColSentAfter>0?colptr[curLastColumn+1-firstCol]-colptr[startColAfter-firstCol]:0;


    //logfileptr->OFS()<<"numColSentBefore: "<<numColSentBefore<<std::endl;
    //logfileptr->OFS()<<" numColSentAfter: "<<numColSentAfter <<std::endl;
    //logfileptr->OFS()<<"numRowSentBefore: "<<numRowSentBefore<<std::endl;
    //logfileptr->OFS()<<" numRowSentAfter: "<<numRowSentAfter <<std::endl;



    //    logfileptr->OFS()<<"colptr was: "<<colptr<<std::endl;
    //    logfileptr->OFS()<<"rowind was: "<<rowind<<std::endl;


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
      //logfileptr->OFS()<<"}"<<std::endl;

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
      //logfileptr->OFS()<<"}"<<std::endl;
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

    //    logfileptr->OFS()<<"ssizes: "<<ssizes<<std::endl;
    //    logfileptr->OFS()<<"sdispls orig : "<<sdispls<<std::endl;


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

    //    logfileptr->OFS()<<"rsizes: "<<rsizes<<std::endl;
    //    logfileptr->OFS()<<"rdispls orig : "<<rdispls<<std::endl;


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


    //        logfileptr->OFS()<<"Recv colptr: "<<std::endl;
    //    for(int p=0;p<mpisize;p++){
    //      if(rsizes[p]>0){
    //        int numCols = (int)(rsizes[p]/sizeof(Ptr))-1;
    //        logfileptr->OFS()<<"P"<<p<<": ";
    //        for(int col = 0;col<=numCols;col++){
    //          logfileptr->OFS()<<ptrBuf[col]<<" ";
    //        }
    //        logfileptr->OFS()<<std::endl;
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
    //logfileptr->OFS()<<"colptr now is: "<<colptr<<std::endl;

    if(recvCntBefore>0){
      std::copy_backward(rowind.begin(),rowind.begin()+rowindSize,rowind.begin()+recvCntBefore+rowindSize);
    }


    //copy from recvBuf

    //if I have received something from ranks before me
    Idx rowindPos =0;
    Idx idxBufPos = 0;
    Idx * idxBuf = (Idx*)&recvBuf[0];



    //logfileptr->OFS()<<"Recv rowind: "<<std::endl;
    //for(int p=0;p<mpisize;p++){
    //  if(rsizesR[p]>0){
    //    int numRows = (int)(rsizesR[p]/sizeof(Idx));
    //    logfileptr->OFS()<<"P"<<p<<": ";
    //    for(int row = 0;row<=numRows;row++){
    //      logfileptr->OFS()<<idxBuf[row]<<" ";
    //    }
    //    logfileptr->OFS()<<std::endl;
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
    //logfileptr->OFS()<<"rowind now is: "<<rowind<<std::endl;

    //logfileptr->OFS()<<"Vertex distribution was: "<<vertexDist<<std::endl;
    supLastColumn = supLastColumn + baseval +1;
    MPI_Allgather(&supLastColumn,sizeof(supLastColumn),MPI_BYTE,&vertexDist[1],sizeof(supLastColumn),MPI_BYTE,comm);
    //logfileptr->OFS()<<"and is now: "<<vertexDist<<std::endl;


    //SparseMatrixGraph sgr;
    //AllGatherStructure(sgr);

    //logfileptr->OFS()<<sg.colptr<<std::endl;
    //logfileptr->OFS()<<sgr.colptr<<std::endl;
    //logfileptr->OFS()<<sg.rowind<<std::endl;
    //logfileptr->OFS()<<sgr.rowind<<std::endl;

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

        //for(auto it = ssizes.begin();it!=ssizes.end();it++){  (*it)*=sizeof(triplet<F>);}
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
        int totRecv = rdispls.back();///sizeof(triplet<F>);
        Irecv.resize(totRecv);

        MPI_Alltoallv(&Isend[0],&ssizes[0],&sdispls[0],type,&Irecv[0],&rsizes[0],&rdispls[0],type,comm);
        SYMPACK_TIMER_STOP(DistMat_Expand_communication);

        MPI_Type_free(&type);
        //now parse
        //for(auto it = rsizes.begin();it!=rsizes.end();it++){  (*it)/=sizeof(triplet<F>);}


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

        expanded =true;
        //keepDiag = 1;

//        if(GetSorted()){
//          scope_timer(a,DistMat_Expand_sort);
//          SetSorted(false);
//          SortGraph(); 
//        }
      if(!sorted){
        SetSorted(1);
      }



/////      if(!sorted){
/////        SetSorted(1);
/////      }
/////
/////      MPI_Op MPI_SYMPACK_PTR_SUM; 
/////      MPI_Op MPI_SYMPACK_PTR_MAX; 
/////      MPI_Datatype MPI_SYMPACK_PTR; 
/////      MPI_Type_contiguous( sizeof(Ptr), MPI_BYTE, &MPI_SYMPACK_PTR ); 
/////      MPI_Type_commit( &MPI_SYMPACK_PTR ); 
/////      MPI_Op_create( PtrSum, true, &MPI_SYMPACK_PTR_SUM ); 
/////      MPI_Op_create( PtrMax, true, &MPI_SYMPACK_PTR_MAX ); 
/////
/////      Idx N = size; 
/////
/////      Idx locColCnt = LocalVertexCount();
/////      //make copies first
/////      std::vector<Ptr> prevColptr;
/////      std::vector<Idx> prevRowind;
/////      prevColptr.swap(colptr);
/////      prevRowind.swap(rowind);
/////
/////      Ptr * pcolptr = &prevColptr[0];
/////      Idx * prowind = &prevRowind[0];
/////      Ptr locNNZ = prevRowind.size();
/////
/////      //first generate the expanded distributed structure
/////      //            Int firstLocCol = (mpirank)*colPerProc; //0 based
/////      Idx firstLocCol = LocalFirstVertex()-baseval; //0 based
/////      Idx maxLocN = 0;
/////      for(int p = 0; p<mpisize;p++){maxLocN = max(maxLocN, vertexDist[p+1]-vertexDist[p]);}
/////      //               max(locColCnt,N-(mpisize-1)*colPerProc); // can be 0
/////      std::vector<Ptr> remote_colptr(maxLocN+1);
/////      std::vector<Idx> remote_rowind;
/////      std::vector<Ptr> remote_rowindPos(maxLocN+1);
/////      std::vector<Ptr> curPos(locColCnt);
/////      std::vector<Ptr> prevPos(locColCnt);
/////
/////      std::copy(pcolptr,pcolptr+locColCnt,curPos.begin());
/////      for(Int prow = 0; prow<mpisize; prow++){
/////        Int firstRemoteCol = vertexDist[prow]-baseval;//(prow)*colPerProc; // 0 based
/////        Int pastLastRemoteCol = vertexDist[prow+1]-baseval;//prow==mpisize-1?N:(prow+1)*colPerProc; // 0 based
/////        Ptr remColCnt = pastLastRemoteCol - firstRemoteCol;
/////        Ptr maxExtraNNZ = 0; 
/////        Ptr extraNNZ = 0;
/////        //receive from all the previous processor
/////        //receive extra nnzcnt...
/////
/////        std::fill(remote_colptr.begin(),remote_colptr.end(),0);
/////        if(mpirank<=prow){
/////          //backup the position in each column
/////          std::copy(curPos.begin(),curPos.end(),prevPos.begin());
/////          //use remote_colptr to store the number of nnz per col first
/////          //loop through my local columns and figure out the extra nonzero per col on prow
/////          for(Idx locCol = 0 ; locCol< locColCnt; locCol++){
/////            Idx col = firstLocCol + locCol;  // 0 based
/////            Ptr colbeg = curPos[locCol]-baseval; //now 0 based
/////            Ptr colend = pcolptr[locCol+1]-baseval; // now 0 based
/////            for(Ptr rptr = colbeg ; rptr< colend ; rptr++ ){
/////              Idx row = prowind[rptr]-baseval; //0 based
/////              assert(row>=firstRemoteCol);
/////              if(row>col){
/////                //this goes onto prow
/////                if(row<pastLastRemoteCol){
/////                  //this is shifted by one to compute the prefix sum without copying afterwards
/////                  remote_colptr[row-firstRemoteCol+1]++;
/////                  extraNNZ++;
/////                }
/////                else{
/////                  break;
/////                }
/////              }
/////              curPos[locCol]++; // baseval based
/////            }
/////          }
/////
/////          //put the local sum into the first element of the array
/////          remote_colptr[0] = 0;
/////          for(Idx p = 1; p<remColCnt+1;p++){ remote_colptr[0] += remote_colptr[p];}
/////
/////
/////          if(mpirank==prow){
/////            //we can now receive the number of NZ per col into the expanded pcolptr
/////            assert(locColCnt == remColCnt);
/////            colptr.assign(remColCnt+1,0);
/////            //this array will contain the max element in our custom reduce
/////            colptr[0] = 0;
/////            SYMPACK_TIMER_START(REDUCE);
/////            MPI_Reduce(&remote_colptr[1],&colptr[1],remColCnt,MPI_SYMPACK_PTR,MPI_SYMPACK_PTR_SUM,prow,comm);
/////            MPI_Reduce(&remote_colptr[0],&colptr[0],1,MPI_SYMPACK_PTR,MPI_SYMPACK_PTR_MAX,prow,comm);
/////            SYMPACK_TIMER_STOP(REDUCE);
/////
/////            maxExtraNNZ = colptr[0];
/////            remote_rowind.resize(maxExtraNNZ);
/////
/////            for(Idx locCol = 0 ; locCol< locColCnt; locCol++){
/////              Ptr colbeg = pcolptr[locCol]-baseval; //now 0 based
/////              Ptr colend = pcolptr[locCol+1]-baseval; // now 0 based
/////              colptr[locCol+1] += colend - colbeg; //- 1+ keepDiag;
/////              //At this point expColptr[locCol+1] contains the total number of NNZ of locCol 
/////              //we can now compute the expanded pcolptr
/////            }
/////
/////            colptr[0] = baseval;
/////            for(Idx col = 1;col<remColCnt+1;col++){ colptr[col]+=colptr[col-1]; }
/////          }
/////          else{
/////            SYMPACK_TIMER_START(REDUCE);
/////            MPI_Reduce(&remote_colptr[1],NULL,remColCnt,MPI_SYMPACK_PTR,MPI_SYMPACK_PTR_SUM,prow,comm);
/////            MPI_Reduce(&remote_colptr[0],NULL,1,MPI_SYMPACK_PTR,MPI_SYMPACK_PTR_MAX,prow,comm);
/////            SYMPACK_TIMER_STOP(REDUCE);
/////            remote_rowind.resize(extraNNZ);
/////          }
/////
/////
/////
/////
////////TODO
/////
/////
/////          /**************     Compute remote_colptr from the local nnz per col ***********/
/////          //compute a prefix sum of the nnz per column to get the new pcolptr
/////          remote_colptr[0] = baseval;
/////          for(Idx col = 1;col<=remColCnt;col++){ remote_colptr[col]+=remote_colptr[col-1]; }
/////
/////          /**************     Fill remote_rowind now ****************/
/////          //make a copy of pcolptr in order to track the current position in each column
/////          std::copy(&remote_colptr[0],&remote_colptr[0]+remColCnt+1,remote_rowindPos.begin());
/////          //loop through my local columns and figure out the extra nonzero per col on remote processors
/////          for(Idx locCol = 0 ; locCol< locColCnt; locCol++){
/////            Idx col = firstLocCol + locCol;  // 0 based
/////            Ptr colbeg = prevPos[locCol]-baseval; //now 0 based
/////            Ptr colend = curPos[locCol]-baseval; // now 0 based
/////
/////            for(Ptr rptr = colbeg ; rptr< colend ; rptr++ ){
/////              Idx row = prowind[rptr]-baseval; //0 based
/////              if(row>col){
/////                //this goes onto prow
/////                Idx locRow = row - firstRemoteCol;      
/////                remote_rowind[ remote_rowindPos[locRow]++ - baseval ] = col + baseval;
/////              }
/////            }
/////          }
/////
/////#ifdef DEBUG
/////          logfileptr->OFS()<<"remote_colptr ";
/////          for(Idx col = 0;col<=remColCnt;col++){
/////            logfileptr->OFS()<<remote_colptr[col]<<" ";
/////          }
/////          logfileptr->OFS()<<std::endl;
/////
/////          logfileptr->OFS()<<"remote_rowind ";
/////          for(Ptr col = 0;col<extraNNZ;col++){
/////            logfileptr->OFS()<<remote_rowind[col]<<" ";
/////          }
/////          logfileptr->OFS()<<std::endl;
/////#endif
/////
/////          if(prow==mpirank){
/////#ifdef DEBUG
/////            logfileptr->OFS()<<"expColptr "<<colptr<<std::endl;
/////#endif
/////            if(colptr.size()>0){
/////              rowind.resize(colptr.back()-baseval);
/////            }
/////            std::copy(colptr.begin(),colptr.end(),remote_rowindPos.begin()); 
/////            for(Idx locCol = 0 ; locCol< locColCnt; locCol++){
/////              Idx col = firstLocCol + locCol;  // 0 based
/////              Ptr colbeg = pcolptr[locCol]-baseval; //now 0 based
/////              Ptr colend = pcolptr[locCol+1]-baseval; // now 0 based
/////              Ptr & pos = remote_rowindPos[locCol];
/////
/////              //copy the local lower triangular part into the expanded structure
/////              for(Ptr rptr = colbeg ; rptr< colend ; rptr++ ){
/////                Idx row = prowind[rptr]-baseval; //0 based
/////                if(col!=row || keepDiag){
/////                  rowind[pos++ - baseval] = row + baseval;  
/////                }
/////              }
/////
/////              //copy the local extra NNZ into the expanded structure
/////              colbeg = remote_colptr[locCol]-baseval; //now 0 based
/////              colend = remote_colptr[locCol+1]-baseval; // now 0 based 
/////              for(Ptr rptr = colbeg ; rptr< colend ; rptr++ ){
/////                Idx row = remote_rowind[rptr]-baseval; //0 based
/////                rowind[pos++ - baseval] = row + baseval;  
/////              }
/////            }
/////
/////            SYMPACK_TIMER_START(RECV);
/////            for(Int pcol = 0; pcol<prow; pcol++){
/////              //Use an MPI_Gatherv instead ? >> memory usage : p * n/p
/////              //Do mpi_recv from any, anytag for pcolptr and then do the matching rowind ?
/////
/////              //logfileptr->OFS()<<"P"<<mpirank<<" receives pcolptr from P"<<pcol<<std::endl;
/////              //receive colptrs...
/////              MPI_Status status;
/////              MPI_Recv(&remote_colptr[0],(remColCnt+1)*sizeof(Ptr),MPI_BYTE,MPI_ANY_SOURCE,prow,comm,&status);
/////
/////              //logfileptr->OFS()<<"P"<<mpirank<<" receives rowind from P"<<pcol<<std::endl;
/////              //receive rowinds...
/////              MPI_Recv(&remote_rowind[0],maxExtraNNZ*sizeof(Idx),MPI_BYTE,status.MPI_SOURCE,prow,comm,MPI_STATUS_IGNORE);
/////
/////              SYMPACK_TIMER_START(PROCESSING_RECV_DATA);
/////              //logfileptr->OFS()<<"P"<<mpirank<<" done receiving from P"<<pcol<<std::endl;
/////              for(Idx locCol = 0 ; locCol< locColCnt; locCol++){
/////                Idx col = firstLocCol + locCol;  // 0 based
/////                //copy the extra NNZ into the expanded structure
/////                Ptr colbeg = remote_colptr[locCol]-baseval; //now 0 based
/////                Ptr colend = remote_colptr[locCol+1]-baseval; // now 0 based 
/////                Ptr & pos = remote_rowindPos[locCol];
/////                for(Ptr rptr = colbeg ; rptr< colend ; rptr++ ){
/////                  Idx row = remote_rowind[rptr]-baseval; //0 based
/////                  rowind[pos++ - baseval] = row + baseval;  
/////                }
/////              }
/////              SYMPACK_TIMER_STOP(PROCESSING_RECV_DATA);
/////            }
/////            SYMPACK_TIMER_STOP(RECV);
/////
/////            if(sorted){ 
/////              SortEdges();
/////            }
/////
/////#ifdef DEBUG
/////            logfileptr->OFS()<<"expRowind "<<rowind<<std::endl;
/////            //logfileptr->OFS()<<"true expRowind "<<Global.expRowind<<std::endl;
/////#endif
/////
/////          }
/////          else{
/////            SYMPACK_TIMER_START(SEND);
/////            MPI_Send(&remote_colptr[0],(remColCnt+1)*sizeof(Ptr),MPI_BYTE,prow,prow,comm);
/////            MPI_Send(&remote_rowind[0],extraNNZ*sizeof(Idx),MPI_BYTE,prow,prow,comm);
/////            SYMPACK_TIMER_STOP(SEND);
/////          }
/////
/////        }
/////        else{
/////          SYMPACK_TIMER_START(REDUCE);
/////            MPI_Reduce(&remote_colptr[1],NULL,remColCnt,MPI_SYMPACK_PTR,MPI_SYMPACK_PTR_SUM,prow,comm);
/////            MPI_Reduce(&remote_colptr[0],NULL,1,MPI_SYMPACK_PTR,MPI_SYMPACK_PTR_MAX,prow,comm);
/////          SYMPACK_TIMER_STOP(REDUCE);
/////        }
/////
/////      }
/////
/////      MPI_Op_free(&MPI_SYMPACK_PTR_MAX);
/////      MPI_Op_free(&MPI_SYMPACK_PTR_SUM);
/////      MPI_Type_free(&MPI_SYMPACK_PTR);
/////
/////      MPI_Barrier(comm);
/////
/////      //nnz = LocalEdgeCount();
/////      expanded =true;
    }
    SYMPACK_TIMER_STOP(EXPAND);
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
    for(int p = 0; p<mpisize;p++){rsizes[p] = (int)remoteVertexCnt[p]*sizeof(Ptr);}
    std::vector<int> rdispls(mpisize+1,0);
    rdispls[0]=0;
    std::partial_sum(rsizes.begin(),rsizes.end(),rdispls.begin()+1);
    MPI_Allgatherv(&colptr[0],localVertexCnt*sizeof(Ptr),MPI_BYTE,&g.colptr[0],&rsizes[0],&rdispls[0],MPI_BYTE,comm);


    Ptr localEdgeCnt = LocalEdgeCount();
    std::vector<Ptr> remoteEdgeCnt(mpisize,0);
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
    if(iam==proot){
      Idx totalVertexCnt = std::accumulate(remoteVertexCnt.begin(),remoteVertexCnt.end(),0,std::plus<Idx>());
      g.colptr.resize(totalVertexCnt+1);


      //compute receive displacements
      std::vector<int> rsizes(mpisize,0);
      for(int p = 0; p<mpisize;p++){rsizes[p] = (int)remoteVertexCnt[p]*sizeof(Ptr);}
      std::vector<int> rdispls(mpisize+1,0);
      rdispls[0]=0;
      std::partial_sum(rsizes.begin(),rsizes.end(),rdispls.begin()+1);
      MPI_Gatherv(&colptr[0],localVertexCnt*sizeof(Ptr),MPI_BYTE,&g.colptr[0],&rsizes[0],&rdispls[0],MPI_BYTE,proot,comm);
    }
    else{
    MPI_Gatherv(&colptr[0],localVertexCnt*sizeof(Ptr),MPI_BYTE,NULL,NULL,NULL,MPI_BYTE,proot,comm);
    }


    Ptr localEdgeCnt = LocalEdgeCount();
    std::vector<Ptr> remoteEdgeCnt(mpisize,0);
    MPI_Gather(&localEdgeCnt,sizeof(localEdgeCnt),MPI_BYTE,&remoteEdgeCnt[0],sizeof(localEdgeCnt),MPI_BYTE,proot,comm);
    Ptr totalEdgeCnt = std::accumulate(remoteEdgeCnt.begin(),remoteEdgeCnt.end(),0,std::plus<Ptr>());


    if(iam==proot){
    g.rowind.resize(totalEdgeCnt);
    //compute receive displacements
      std::vector<int> rsizes(mpisize,0);
      std::vector<int> rdispls(mpisize+1,0);
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


}
