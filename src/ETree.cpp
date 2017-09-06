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
/// @file etree.cpp
/// @brief Implementation of the elimination-tree related algorithms.
/// @author Mathias Jacquelin
/// @date 2012-08-31
#include "sympack/ETree.hpp"
#include "sympack/utility.hpp"

#include <upcxx.h>

namespace symPACK{


  void DisjointSet::Initialize(Int n){
    pp_.resize(n,0);
    root_.resize(n);
  }

  void DisjointSet::Finalize(){
    pp_.resize(0);
  }
}

namespace symPACK{
  ETree::ETree(){
    bIsPostOrdered_=false;
  }


  void ETree::BTreeToPO(std::vector<Int> & fson, std::vector<Int> & brother, std::vector<Int> & invpos){
    //Do a depth first search to construct the postordered tree
    std::vector<Int> stack(n_);
    invpos.resize(n_);

    Int stacktop=0, vertex=n_,m=0;
    bool exit = false;
    while( m<n_){
      do{
        stacktop++;
        stack[stacktop-1] = vertex;
        vertex = fson[vertex-1];
      }while(vertex>0);

      while(vertex==0){

        if(stacktop<=0){
          exit = true;
          break;
        }
        vertex = stack[stacktop-1];
        stacktop--;
        m++;

        invpos[vertex-1] = m;
        vertex = brother[vertex-1];
      }

      if(exit){
        break;
      }
    }
  }


  void ETree::PostOrderTree(Ordering & aOrder,Int * relinvp){
    if(n_>0 && !bIsPostOrdered_){

      SYMPACK_TIMER_START(PostOrder);

      std::vector<Int> fson(n_,0);
      std::vector<Int> & brother = poparent_;
      brother.resize(n_,0);

      Int lroot = n_;
      for(Int vertex=n_-1; vertex>0; vertex--){
        Int curParent = parent_[vertex-1];
        if(curParent==0 || curParent == vertex){
          brother[lroot-1] = vertex;
          lroot = vertex;
        }
        else{
          brother[vertex-1] = fson[curParent-1];
          fson[curParent-1] = vertex;
        }
      }


#ifdef _DEBUG_
      logfileptr->OFS()<<"parent "<<parent_<<std::endl;
      logfileptr->OFS()<<"fson "<<fson<<std::endl;
      logfileptr->OFS()<<"brother "<<brother<<std::endl;
#endif

      std::vector<Int> invpos;
      BTreeToPO(fson,brother,invpos);

      //modify the parent list ?
      // node i is now node invpos(i-1)
      for(Int i=1; i<=n_;i++){
        Int nunode = invpos[i-1];
        Int ndpar = parent_[i-1];
        if(ndpar>0){
          ndpar = invpos[ndpar-1];
        }
        poparent_[nunode-1] = ndpar;
      }

      //we need to compose aOrder.invp and invpos
      aOrder.Compose(invpos);

      if(relinvp!=nullptr){
        std::copy(invpos.begin(),invpos.end(),relinvp);
      }

      bIsPostOrdered_ = true;


#ifdef _DEBUG_
      logfileptr->OFS()<<"new parent: "<<brother<<std::endl;
#endif

      SYMPACK_TIMER_STOP(PostOrder);
    }

  }


  void ETree::SortChildren(std::vector<Int> & cc, Ordering & aOrder){
    if(!bIsPostOrdered_){
      this->PostOrderTree(aOrder);
    }

    std::vector<Int> fson(n_,0);
    std::vector<Int> brother(n_,0);
    std::vector<Int> lson(n_,0);

    //Get Binary tree representation
    Int lroot = n_;
    for(Int vertex=n_-1; vertex>0; vertex--){
      Int curParent = PostParent(vertex-1);
      //Int curParent = parent_(vertex-1);
      if(curParent==0 || curParent == vertex){
        brother[lroot-1] = vertex;
        lroot = vertex;
      }
      else{
        Int ndlson = lson[curParent-1];
        if(ndlson > 0){
          if  ( cc[vertex-1] >= cc[ndlson-1] ) {
            brother[vertex-1] = fson[curParent-1];
            fson[curParent-1] = vertex;
          }
          else{                                                                                                                                            
            brother[ndlson-1] = vertex;
            lson[curParent-1] = vertex;
          }                                                                                                                                                               
        }
        else{
          fson[curParent-1] = vertex;
          lson[curParent-1] = vertex;
        }
      }
    }
    brother[lroot-1]=0;


    std::vector<Int> invpos;
    //      std::vector<Int> invperm;

    //Compute the parent permutation and update postNumber_
    //Do a depth first search to construct the postordered tree
    std::vector<Int> stack(n_);
    invpos.resize(n_);
    //      invperm.Resize(n_);

    Int stacktop=0, vertex=n_,m=0;
    bool exit = false;
    while( m<n_){
      do{
        stacktop++;
        stack[stacktop-1] = vertex;
        vertex = fson[vertex-1];
      }while(vertex>0);

      while(vertex==0){

        if(stacktop<=0){
          exit = true;
          break;
        }
        vertex = stack[stacktop-1];
        stacktop--;
        m++;

        invpos[vertex-1] = m;

        vertex = brother[vertex-1];
      }

      if(exit){
        break;
      }
    }


    //modify the parent list ?
    // node i is now node invpos(i-1)
    for(Int i=1; i<=n_;i++){
      Int nunode = invpos[i-1];
      Int ndpar = poparent_[i-1];
      if(ndpar>0){
        ndpar = invpos[ndpar-1];
      }
      brother[nunode-1] = ndpar;
    }
    poparent_ = brother;




    //Permute CC     
    for(Int node = 1; node <= n_; ++node){
      Int nunode = invpos[node-1];
      stack[nunode-1] = cc[node-1];
    }

    for(Int node = 1; node <= n_; ++node){
      cc[node-1] = stack[node-1];
    }

    //Compose the two permutations
    aOrder.Compose(invpos);

  }

  void ETree::ConstructETree(DistSparseMatrixGraph & aDistExp, Ordering & aOrder){
    bIsPostOrdered_=false;
    n_ = aDistExp.size;

    int mpisize;
    MPI_Comm_size(aDistExp.GetComm(),&mpisize);

    int mpirank;
    MPI_Comm_rank(aDistExp.GetComm(),&mpirank);

#if 1
    int expandedGraph = aDistExp.expanded;
    
    if ( !expandedGraph ){
      throw std::logic_error( "DistSparseMatrixGraph must be expanded and permuted\n" );
    }

    parent_.assign(n_,0);
    std::vector<Int> ancstr(n_,0);

    Idx fc = aDistExp.LocalFirstVertex(); //1 - based

    MPI_Datatype Inttype;
    MPI_Type_contiguous( sizeof(Int), MPI_BYTE, &Inttype );
    MPI_Type_commit(&Inttype);

    if(mpirank>0){
      MPI_Recv(&parent_[0],n_,Inttype,mpirank-1,mpirank-1,aDistExp.GetComm(),MPI_STATUS_IGNORE);
      MPI_Recv(&ancstr[0],n_,Inttype,mpirank-1,mpirank-1,aDistExp.GetComm(),MPI_STATUS_IGNORE);
    }
    for(Idx locCol = 0; locCol< aDistExp.LocalVertexCount(); locCol++){
      Idx i = fc + locCol; // 1 - based;
      parent_[i-1] = 0;
      ancstr[i-1] = 0;

      //logfileptr->OFS()<<i<<" parent_ "<<parent_<<std::endl;
      //logfileptr->OFS()<<i<<" ancstr "<<ancstr<<std::endl;
      //logfileptr->OFS()<<"i = "<<i<<std::endl;

      Ptr jstrt = aDistExp.colptr[locCol]; //1-based
      Ptr jstop = aDistExp.colptr[locCol+1] -1;//1-based
      if (jstrt<jstop){
        for(Ptr j = jstrt; j<=jstop; ++j){
          Idx nbr = aDistExp.rowind[j-1]; //1-based
//            logfileptr->OFS()<<"   nbr = "<<nbr<<std::endl;
          if  ( nbr < i ){
            //                       -------------------------------------------
            //                       for each nbr, find the root of its current
            //                       elimination tree.  perform path compression
            //                       as the subtree is traversed.
            //                       -------------------------------------------
            Int break_loop = 0;
            if  ( ancstr[nbr-1] == i ){
              break_loop = 1;
            }
            else{

              while(ancstr[nbr-1] >0){
                if  ( ancstr[nbr-1] == i ){
                  break_loop = 1;
                  break;
                }
                Int next = ancstr[nbr-1];
                ancstr[nbr-1] = i;
                nbr = next;
              }
              //              logfileptr->OFS()<<std::endl;
              //                       --------------------------------------------
              //                       now, nbr is the root of the subtree.  make i
              //                       the parent node of this root.
              //                       --------------------------------------------
              if(!break_loop){
                parent_[nbr-1] = i;
                ancstr[nbr-1] = i;
              }
            }
          }
        }
      }
    }

    if(mpirank<mpisize-1){
      //        logfileptr->OFS()<<"my parent now is: "<<myParent<<std::endl;
      //        logfileptr->OFS()<<"my ancstr now is: "<<myAncstr<<std::endl;
      //logfileptr->OFS()<<"parent now is: "<<parent_<<std::endl;
      //logfileptr->OFS()<<"ancstr now is: "<<ancstr<<std::endl;
      MPI_Send(&parent_[0],n_,Inttype,mpirank+1,mpirank,aDistExp.GetComm());
      MPI_Send(&ancstr[0],n_ ,Inttype,mpirank+1,mpirank,aDistExp.GetComm());
    }

    //Now proc mpisize-1 bcast the parent_ array
    MPI_Bcast(&parent_[0],n_,Inttype,mpisize-1,aDistExp.GetComm());


    MPI_Type_free( &Inttype );
#else
    throw std::logic_error( "ETree::ConstructETree(DistSparseMatrixGraph & , Ordering & ) not implemented\n" );
    DistSparseMatrixGraph tmpGraph = aDistExp;

    //Expand to unsymmetric storage
    tmpGraph.ExpandSymmetric();

    SYMPACK_TIMER_START(Construct_Etree);
    //std::fill(parent_.begin(),parent_.end(),0);



#if 1 
    //first permute locally then redistribute with alltoallv then do the etree
    tmpGraph.Permute(&aOrder.invp[0]);

    parent_.assign(n_,0);
    std::vector<Int> ancstr(n_,-1);

    Idx fc = tmpGraph.LocalFirstVertex()-tmpGraph.GetBaseval(); //0 - based

    if(mpirank>0){
      MPI_Recv(&parent_[0],fc*sizeof(Int),MPI_BYTE,mpirank-1,mpirank-1,tmpGraph.GetComm(),MPI_STATUS_IGNORE);
      MPI_Recv(&ancstr[0],fc*sizeof(Int),MPI_BYTE,mpirank-1,mpirank-1,tmpGraph.GetComm(),MPI_STATUS_IGNORE);
    }
    for(Idx locCol = 0; locCol< tmpGraph.LocalVertexCount(); locCol++){

      Idx i = fc + locCol; // 0 - based;
      parent_[i] = 0;
      ancstr[i] = 0;

      //logfileptr->OFS()<<"i = "<<i+1<<std::endl;

      Ptr jstrt = tmpGraph.colptr[locCol] - tmpGraph.GetBaseval(); //0-based
      Ptr jstop = tmpGraph.colptr[locCol+1] - tmpGraph.GetBaseval();//0-based
      if(jstrt<jstop-1){
        for(Ptr j = jstrt; j<jstop; ++j){
          Idx nbr = tmpGraph.rowind[j] - tmpGraph.GetBaseval(); //0-based
          //logfileptr->OFS()<<"   nbr = "<<nbr+1<<std::endl;
          if  ( nbr < i ){
            // -------------------------------------------
            // for each nbr, find the root of its current
            // elimination tree.  perform path compression
            // as the subtree is traversed.
            // -------------------------------------------
            //column i (unpermuted) is not the parent of column nbr
            if  ( ancstr[nbr] != i+1 ){
              //logfileptr->OFS()<<"path: "<<nbr<<" ";
              Int break_loop = 0;
              while(ancstr[nbr] >0){
                if  ( ancstr[nbr] == i+1 ){
                  break_loop = 1;
                  break;
                }
                Int next = ancstr[nbr];
                ancstr[nbr] = i+1;
                nbr = next - 1;
                //logfileptr->OFS()<<nbr+1<<" ";
              }
              //logfileptr->OFS()<<std::endl;

              // --------------------------------------------
              // now, nbr is the root of the subtree.  make i
              // the parent node of this root.
              // --------------------------------------------
              if(!break_loop){
                parent_[nbr] = i + 1; // 1-based
                ancstr[nbr] = i + 1; // 1-based
              }
            }
          }
        }
      }

    }

    if(mpirank<mpisize-1){
      //        logfileptr->OFS()<<"my parent now is: "<<myParent<<std::endl;
      //        logfileptr->OFS()<<"my ancstr now is: "<<myAncstr<<std::endl;
      //logfileptr->OFS()<<"parent now is: "<<parent_<<std::endl;
      //logfileptr->OFS()<<"ancstr now is: "<<ancstr<<std::endl;
      MPI_Send(&parent_[0],(fc+tmpGraph.LocalVertexCount())*sizeof(Int),MPI_BYTE,mpirank+1,mpirank,tmpGraph.GetComm());
      MPI_Send(&ancstr[0],(fc+tmpGraph.LocalVertexCount())*sizeof(Int),MPI_BYTE,mpirank+1,mpirank,tmpGraph.GetComm());
    }

    //Now proc mpisize-1 bcast the parent_ array
    MPI_Bcast(&parent_[0],n_*sizeof(Int),MPI_BYTE,mpisize-1,tmpGraph.GetComm());


#else
#endif
#endif

    SYMPACK_TIMER_STOP(Construct_Etree);

  }


  void ETree::ConstructETree(SparseMatrixGraph & sgraph, Ordering & aOrder, MPI_Comm & aComm){
    int iam =0;
    int np =1;
    MPI_Comm_rank(aComm,&iam);
    MPI_Comm_size(aComm,&np);


    bIsPostOrdered_=false;
    n_ = sgraph.size;

    if(iam == 0 && (!sgraph.IsExpanded() ) ){
      throw std::logic_error( "SparseMatrixGraph must be expanded\n" );
    }

    SYMPACK_TIMER_START(Construct_Etree_Classic);
    parent_.assign(n_,0);

    if(iam==0){
      std::vector<Int> ancstr(n_);

      for(Int i = 1; i<=n_; ++i){
        parent_[i-1] = 0;
        ancstr[i-1] = 0;
      //logfileptr->OFS()<<i<<" parent_ "<<parent_<<std::endl;
      //logfileptr->OFS()<<i<<" ancstr "<<ancstr<<std::endl;

        Int node = aOrder.perm[i-1];
        //          logfileptr->OFS()<<"i = "<<node<<std::endl;
        Ptr jstrt = sgraph.colptr[node-1];
        Ptr jstop = sgraph.colptr[node] - 1;
        if  ( jstrt < jstop ){
          for(Ptr j = jstrt; j<=jstop; ++j){
            Idx nbr = sgraph.rowind[j-1];
//            logfileptr->OFS()<<"   nbr = "<<nbr<<std::endl;
            //logfileptr->OFS()<<"   nbr = "<<nbr<<"  ";
            nbr = aOrder.invp[nbr-1];
            //logfileptr->OFS()<<"|  nbr = "<<nbr<<std::endl;
            if  ( nbr < i ){
              //                       -------------------------------------------
              //                       for each nbr, find the root of its current
              //                       elimination tree.  perform path compression
              //                       as the subtree is traversed.
              //                       -------------------------------------------
              Int break_loop = 0;
              if  ( ancstr[nbr-1] == i ){
                break_loop = 1;
              }
              else{

                while(ancstr[nbr-1] >0){
                  if  ( ancstr[nbr-1] == i ){
                    break_loop = 1;
                    break;
                  }
                  Int next = ancstr[nbr-1];
                  ancstr[nbr-1] = i;
                  nbr = next;
                }
                //              logfileptr->OFS()<<std::endl;
                //                       --------------------------------------------
                //                       now, nbr is the root of the subtree.  make i
                //                       the parent node of this root.
                //                       --------------------------------------------
                if(!break_loop){
                  parent_[nbr-1] = i;
                  ancstr[nbr-1] = i;
                }
              }
            }
          }
        }
      }
    }
    
        MPI_Datatype type;
        MPI_Type_contiguous( sizeof(Int), MPI_BYTE, &type );
        MPI_Type_commit(&type);
    //Broadcast  
    MPI_Bcast(&parent_[0],n_,type,0,aComm);
        MPI_Type_free(&type);

    SYMPACK_TIMER_STOP(Construct_Etree_Classic);

  }






  ETree ETree::ToSupernodalETree(std::vector<Int> & aXsuper,std::vector<Int> & aSupMembership,Ordering & aOrder) const{
    ETree newTree;
    newTree.n_ = aXsuper.size()-1;
    newTree.parent_.resize(aXsuper.size()-1);


    assert(bIsPostOrdered_);

    for(Int snode=1; snode<=newTree.n_; ++snode){
      Int fc = aXsuper[snode-1];
      Int lc = aXsuper[snode]-1;
      Int parent_col = this->PostParent(lc-1);
      Int parentSnode = ( parent_col == 0) ? 0:aSupMembership[parent_col-1];

      newTree.parent_[snode-1] = parentSnode;
#ifdef _DEBUG_
      logfileptr->OFS()<<"parent of curSnode "<<snode<<" is "<<parentSnode<<std::endl;
#endif
    } 

    newTree.poparent_ = newTree.parent_;
    newTree.bIsPostOrdered_ = true;


    return newTree;

  }





  void ETree::DeepestFirst(Ordering & aOrder)
  {
    assert(bIsPostOrdered_);









    std::vector<Int> treesize(n_,0);
    std::vector<Int> depths(n_);
    //first, compute the depth of each node
    for(Int col=n_; col>=1; --col){
      Int parent = PostParent(col-1);
      if(parent==0){
        depths[col-1]=0;
      }
      else{
        treesize[parent-1]++;
        depths[col-1]=depths[parent-1]+1;
      }
    }


    for(Int col=n_; col>=1; --col){
      Int parent = PostParent(col-1);
      if(parent!=0){
        depths[parent-1]=std::max(depths[col-1],depths[parent-1]);
      }
    }




    //relabel the nodes within a subtree based on their depths

  }


}

