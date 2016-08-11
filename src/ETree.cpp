/// @file etree.cpp
/// @brief Implementation of the elimination-tree related algorithms.
/// @author Mathias Jacquelin
/// @date 2012-08-31
#include "sympack/ETree.hpp"
#include "sympack/utility.hpp"

#include <upcxx.h>

namespace SYMPACK{


  void DisjointSet::Initialize(Int n){
    pp_.resize(n,0);
    root_.resize(n);
  }

  void DisjointSet::Finalize(){
    pp_.resize(0);
  }
}

namespace SYMPACK{
  ETree::ETree(){
    bIsPostOrdered_=false;
  }

  ETree::ETree(SparseMatrixStructure & aGlobal, Ordering & aOrder){
    bIsPostOrdered_ = false;
    ConstructETree(aGlobal,aOrder);
  }



  void ETree::BTreeToPO(SYMPACK::vector<Int> & fson, SYMPACK::vector<Int> & brother, SYMPACK::vector<Int> & invpos){
    //Do a depth first search to construct the postordered tree
    SYMPACK::vector<Int> stack(n_);
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


  void ETree::PostOrderTree(Ordering & aOrder){
    if(n_>0 && !bIsPostOrdered_){

      TIMER_START(PostOrder);

      SYMPACK::vector<Int> fson(n_,0);
      SYMPACK::vector<Int> & brother = poparent_;
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

      SYMPACK::vector<Int> invpos;
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



      bIsPostOrdered_ = true;


#ifdef _DEBUG_
      logfileptr->OFS()<<"new parent: "<<brother<<std::endl;
#endif

      TIMER_STOP(PostOrder);
    }

  }


  void ETree::SortChildren(SYMPACK::vector<Int> & cc, Ordering & aOrder){
    if(!bIsPostOrdered_){
      this->PostOrderTree(aOrder);
    }

    SYMPACK::vector<Int> fson(n_,0);
    SYMPACK::vector<Int> brother(n_,0);
    SYMPACK::vector<Int> lson(n_,0);

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


    SYMPACK::vector<Int> invpos;
    //      SYMPACK::vector<Int> invperm;

    //Compute the parent permutation and update postNumber_
    //Do a depth first search to construct the postordered tree
    SYMPACK::vector<Int> stack(n_);
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

    //Sort the parent
  }

  void ETree::ConstructETree(DistSparseMatrixGraph & aDistExp, Ordering & aOrder){
    throw std::logic_error( "ETree::ConstructETree(DistSparseMatrixGraph & , Ordering & ) not implemented\n" );
    bIsPostOrdered_=false;
    n_ = aDistExp.size;

    DistSparseMatrixGraph tmpGraph = aDistExp;

    //Expand to unsymmetric storage
    tmpGraph.ExpandSymmetric();

    TIMER_START(Construct_Etree);
    //std::fill(parent_.begin(),parent_.end(),0);


    int mpisize;
    MPI_Comm_size(tmpGraph.GetComm(),&mpisize);


#if 1 
    int mpirank;
    MPI_Comm_rank(tmpGraph.GetComm(),&mpirank);
    //first permute locally then redistribute with alltoallv then do the etree
    tmpGraph.Permute(&aOrder.invp[0]);

    parent_.assign(n_,0);
    SYMPACK::vector<Int> ancstr(n_,-1);

    Idx fc = tmpGraph.LocalFirstVertex()-tmpGraph.GetBaseval(); //0 - based

    if(mpirank>0){
      MPI_Recv(&parent_[0],fc*sizeof(Int),MPI_BYTE,mpirank-1,mpirank-1,tmpGraph.GetComm(),MPI_STATUS_IGNORE);
      MPI_Recv(&ancstr[0],fc*sizeof(Int),MPI_BYTE,mpirank-1,mpirank-1,tmpGraph.GetComm(),MPI_STATUS_IGNORE);
    }
    for(Idx locCol = 0; locCol< tmpGraph.LocalVertexCount(); locCol++){

      Idx i = fc + locCol; // 0 - based;
      parent_[i] = 0;
      ancstr[i] = 0;

      //logfileptr->OFS()<<"i = "<<i+1<<endl;

      Ptr jstrt = tmpGraph.colptr[locCol] - tmpGraph.GetBaseval(); //0-based
      Ptr jstop = tmpGraph.colptr[locCol+1] - tmpGraph.GetBaseval();//0-based
      if(jstrt<jstop-1){
        for(Ptr j = jstrt; j<jstop; ++j){
          Idx nbr = tmpGraph.rowind[j] - tmpGraph.GetBaseval(); //0-based
          //logfileptr->OFS()<<"   nbr = "<<nbr+1<<endl;
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
              //logfileptr->OFS()<<endl;

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
      //        logfileptr->OFS()<<"my parent now is: "<<myParent<<endl;
      //        logfileptr->OFS()<<"my ancstr now is: "<<myAncstr<<endl;
      logfileptr->OFS()<<"parent now is: "<<parent_<<endl;
      logfileptr->OFS()<<"ancstr now is: "<<ancstr<<endl;
      MPI_Send(&parent_[0],(fc+tmpGraph.LocalVertexCount())*sizeof(Int),MPI_BYTE,mpirank+1,mpirank,tmpGraph.GetComm());
      MPI_Send(&ancstr[0],(fc+tmpGraph.LocalVertexCount())*sizeof(Int),MPI_BYTE,mpirank+1,mpirank,tmpGraph.GetComm());
    }

    //Now proc mpisize-1 bcast the parent_ array
    MPI_Bcast(&parent_[0],n_*sizeof(Int),MPI_BYTE,mpisize-1,tmpGraph.GetComm());


#else
#endif

    TIMER_STOP(Construct_Etree);

  }


  void ETree::ConstructETree(SparseMatrixGraph & sgraph, Ordering & aOrder, upcxx::team * team){
    bIsPostOrdered_=false;
    n_ = sgraph.size;

    if(iam == 0 && (!sgraph.IsExpanded() ) ){
      throw std::logic_error( "SparseMatrixGraph must be expanded\n" );
    }

    TIMER_START(Construct_Etree_Classic);
    parent_.resize(n_,0);

    if(iam==0){
    SYMPACK::vector<Int> ancstr(n_);



    for(Int i = 1; i<=n_; ++i){
      parent_[i-1] = 0;
      ancstr[i-1] = 0;
      Int node = aOrder.perm[i-1];
      //          logfileptr->OFS()<<"i = "<<node<<endl;
      Ptr jstrt = sgraph.colptr[node-1];
      Ptr jstop = sgraph.colptr[node] - 1;
      if  ( jstrt < jstop ){
        for(Ptr j = jstrt; j<=jstop; ++j){
          Idx nbr = sgraph.rowind[j-1];
          //logfileptr->OFS()<<"   nbr = "<<nbr<<"  ";
          nbr = aOrder.invp[nbr-1];
          //          logfileptr->OFS()<<"   nbr = "<<nbr<<endl;
          //logfileptr->OFS()<<"|  nbr = "<<nbr<<endl;
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

              //              logfileptr->OFS()<<"path: "<<nbr<<" ";
              while(ancstr[nbr-1] >0){
                if  ( ancstr[nbr-1] == i ){
                  break_loop = 1;
                  break;
                }
                Int next = ancstr[nbr-1];
                ancstr[nbr-1] = i;
                nbr = next;
                //              logfileptr->OFS()<<nbr<<" ";
              }
              //              logfileptr->OFS()<<endl;
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
    
    //Broadcast  
    team->bcast(&parent_[0], &parent_[0], n_*sizeof(Int), 0);

    TIMER_STOP(Construct_Etree_Classic);

  }




  void ETree::ConstructETree(SparseMatrixStructure & aGlobal, Ordering & aOrder){
    bIsPostOrdered_=false;
    n_ = aGlobal.size;

    //Expand to symmetric storage
    aGlobal.ExpandSymmetric();

    TIMER_START(Construct_Etree_Classic);
    parent_.resize(n_,0);


    SYMPACK::vector<Int> ancstr(n_);



    for(Int i = 1; i<=n_; ++i){
      parent_[i-1] = 0;
      ancstr[i-1] = 0;
      Int node = aOrder.perm[i-1];
      //          logfileptr->OFS()<<"i = "<<node<<endl;
      Ptr jstrt = aGlobal.expColptr[node-1];
      Ptr jstop = aGlobal.expColptr[node] - 1;
      if  ( jstrt < jstop ){
        for(Ptr j = jstrt; j<=jstop; ++j){
          Idx nbr = aGlobal.expRowind[j-1];
          //logfileptr->OFS()<<"   nbr = "<<nbr<<"  ";
          nbr = aOrder.invp[nbr-1];
          //          logfileptr->OFS()<<"   nbr = "<<nbr<<endl;
          //logfileptr->OFS()<<"|  nbr = "<<nbr<<endl;
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

              //              logfileptr->OFS()<<"path: "<<nbr<<" ";
              while(ancstr[nbr-1] >0){
                if  ( ancstr[nbr-1] == i ){
                  break_loop = 1;
                  break;
                }
                Int next = ancstr[nbr-1];
                ancstr[nbr-1] = i;
                nbr = next;
                //              logfileptr->OFS()<<nbr<<" ";
              }
              //              logfileptr->OFS()<<endl;
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


#if 0
    for(Int node = 1; node<=n_; ++node){
      Int i = aOrder.invp[node-1];
      parent_[i-1] = 0;
      ancstr[i-1] = 0;
      logfileptr->OFS()<<"node = "<<node<<endl; 
      Ptr jstrt = aGlobal.expColptr[node-1];
      Ptr jstop = aGlobal.expColptr[node] - 1;
      if  ( jstrt < jstop ){
        for(Ptr j = jstrt; j<=jstop; ++j){
          Idx nbr = aGlobal.expRowind[j-1];
          logfileptr->OFS()<<"   nbr = "<<nbr<<"  ";
          nbr = aOrder.invp[nbr-1];
          logfileptr->OFS()<<"|  nbr = "<<nbr<<endl;
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

              logfileptr->OFS()<<"path: "<<nbr<<" ";
              while(ancstr[nbr-1] >0){
                if  ( ancstr[nbr-1] == i ){
                  break_loop = 1;
                  break;
                }
                Int next = ancstr[nbr-1];
                ancstr[nbr-1] = i;
                nbr = next;
                logfileptr->OFS()<<nbr<<" ";
              }
              logfileptr->OFS()<<endl;
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
    logfileptr->OFS()<<parent_<<endl;
#endif






    TIMER_STOP(Construct_Etree_Classic);

  }





  ETree ETree::ToSupernodalETree(SYMPACK::vector<Int> & aXsuper,SYMPACK::vector<Int> & aSupMembership,Ordering & aOrder) const{
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









    SYMPACK::vector<Int> treesize(n_,0);
    SYMPACK::vector<Int> depths(n_);
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
        depths[parent-1]=max(depths[col-1],depths[parent-1]);
      }
    }




    //relabel the nodes within a subtree based on their depths

  }


}

