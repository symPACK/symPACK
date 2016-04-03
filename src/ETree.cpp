/// @file etree.cpp
/// @brief Implementation of the elimination-tree related algorithms.
/// @author Mathias Jacquelin
/// @date 2012-08-31
#include "sympack/ETree.hpp"
#include "sympack/utility.hpp"

#include <upcxx.h>

namespace SYMPACK{


  void DisjointSet::Initialize(Int n){
    pp_.resize(n);
    SetValue(pp_,I_ZERO);
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
      //postNumber_.Resize(n_);
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
//          perm[m-1] = vertex;

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

      SYMPACK::vector<Int> fson(n_);
      SetValue(fson, I_ZERO);
      SYMPACK::vector<Int> & brother = poparent_;
      brother.resize(n_);
      SetValue(brother, I_ZERO);

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

      SYMPACK::vector<Int> fson(n_);
      SetValue(fson, I_ZERO);
      SYMPACK::vector<Int> brother(n_);
      SetValue(brother, I_ZERO);
      SYMPACK::vector<Int> lson(n_);
      SetValue(lson, I_ZERO);

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
             //if  ( cc(ToPostOrder(vertex)-1) >= cc(ToPostOrder(ndlson)-1) ) {
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

////      SYMPACK::vector<Int> poParent(n_+1);
////            for(Int i=1; i<=n_;i++){
////              Int nunode = perm_(i-1);
////              Int ndpar = parent_(i-1);
////              if(ndpar>0){
////                ndpar = postNumber_(ndpar-1);
////              }
////              poParent(nunode-1) = ndpar;
////            }
////#ifdef _DEBUG_
////      logfileptr->OFS()<<"ORDERED fson "<<fson<<std::endl;
////      logfileptr->OFS()<<"ORDERED brother "<<brother<<std::endl;
////
////      SYMPACK::vector<Int> poParent(n_+1);
////            for(Int i=1; i<=n_;i++){
////              Int nunode = postNumber_(i-1);
////              Int ndpar = parent_(i-1);
////              if(ndpar>0){
////                ndpar = postNumber_(ndpar-1);
////              }
////              poParent(nunode-1) = ndpar;
////            }
////
////
////      logfileptr->OFS()<<"ORDERED new parent: "<<poParent<<std::endl;
////      logfileptr->OFS()<<"ORDERED postNumber: "<<postNumber_<<std::endl;
////      logfileptr->OFS()<<"ORDERED invPostNumber: "<<invPostNumber_<<std::endl;
////#endif
////
//////      logfileptr->OFS()<<"perm: "<<perm<<std::endl;
//////      logfileptr->OFS()<<"invperm: "<<invperm<<std::endl;



//    return perm;

  }

  void ETree::ConstructETree(DistSparseMatrixGraph & aDistExp, Ordering & aOrder){
    bIsPostOrdered_=false;
    n_ = aDistExp.size;

    //Expand to unsymmetric storage
    aDistExp.ExpandSymmetric();

    TIMER_START(Construct_Etree);
    //std::fill(parent_.begin(),parent_.end(),0);


    int mpisize;
    MPI_Comm_size(aDistExp.comm,&mpisize);


#if 1 
    int mpirank;
    MPI_Comm_rank(aDistExp.comm,&mpirank);
    //first permute locally then redistribute with alltoallv then do the etree
    aDistExp.Permute(&aOrder.invp[0]);

    //TODO This is not ok for memory scaling....
    parent_.resize(n_);

//    SYMPACK::vector<Int> myParent(aDistExp.LocalVertexCount(),-1);
//    SYMPACK::vector<Int> myAncstr(aDistExp.LocalVertexCount(),-1);
    //TODO This is not ok for memory scaling....
    SYMPACK::vector<Int> ancstr(n_,-1);

    Idx fc = (iam)*(n_/mpisize); //0 - based

      if(mpirank>0){
        MPI_Recv(&parent_[0],fc*sizeof(Int),MPI_BYTE,mpirank-1,mpirank-1,aDistExp.comm,MPI_STATUS_IGNORE);
        MPI_Recv(&ancstr[0],fc*sizeof(Int),MPI_BYTE,mpirank-1,mpirank-1,aDistExp.comm,MPI_STATUS_IGNORE);
      }
      for(Idx locCol = 0; locCol< aDistExp.LocalVertexCount(); locCol++){

        Idx i = fc + locCol; // 0 - based;
        parent_[i] = 0;
        ancstr[i] = 0;
//        myParent[locCol] = 0;
//        myAncstr[locCol] = 0;

        Ptr jstrt = aDistExp.colptr[locCol] - aDistExp.baseval; //0-based
        Ptr jstop = aDistExp.colptr[locCol+1] - aDistExp.baseval;//0-based
        if(jstrt<jstop-1){
          for(Ptr j = jstrt; j<jstop; ++j){
            Idx nbr = aDistExp.rowind[j] - aDistExp.baseval; //0-based
            if  ( nbr < i ){
              // -------------------------------------------
              // for each nbr, find the root of its current
              // elimination tree.  perform path compression
              // as the subtree is traversed.
              // -------------------------------------------
              //column i (unpermuted) is not the parent of column nbr
              if  ( ancstr[nbr] != i+1 ){
                Int break_loop = 0;
                while(ancstr[nbr] >0){
                  if  ( ancstr[nbr] == i+1 ){
                    break_loop = 1;
                    break;
                  }
                  Int next = ancstr[nbr];
                  ancstr[nbr] = i+1;
                  nbr = next - 1;
                }

                // --------------------------------------------
                // now, nbr is the root of the subtree.  make i
                // the parent node of this root.
                // --------------------------------------------
                if(!break_loop){
                  parent_[nbr] = i + 1; // 1-based
                  ancstr[nbr] = i + 1; // 1-based
                }
              }

              
//              if  ( myAncstr[nbr-fc] != i+1 ){
//                Int break_loop = 0;
//                while(myAncstr[nbr-fc] >0){
//                  if  ( myAncstr[nbr-fc] == i+1 ){
//                    break_loop = 1;
//                    break;
//                  }
//                  Int next = myAncstr[nbr-fc];
//                  myAncstr[nbr-fc] = i+1;
//                  nbr = next - 1;
//                }
//
//                // --------------------------------------------
//                // now, nbr is the root of the subtree.  make i
//                // the parent node of this root.
//                // --------------------------------------------
//                if(!break_loop){
//                  myParent[nbr-fc] = i + 1; // 1-based
//                  myAncstr[nbr-fc] = i + 1; // 1-based
//                }
//              }

            }
          }
        }

      }

      if(mpirank<mpisize-1){
//        logfileptr->OFS()<<"my parent now is: "<<myParent<<endl;
//        logfileptr->OFS()<<"my ancstr now is: "<<myAncstr<<endl;
        logfileptr->OFS()<<"parent now is: "<<parent_<<endl;
        logfileptr->OFS()<<"ancstr now is: "<<ancstr<<endl;
        MPI_Send(&parent_[0],(fc+aDistExp.LocalVertexCount())*sizeof(Int),MPI_BYTE,mpirank+1,mpirank,aDistExp.comm);
        MPI_Send(&ancstr[0],(fc+aDistExp.LocalVertexCount())*sizeof(Int),MPI_BYTE,mpirank+1,mpirank,aDistExp.comm);
      }

    //Now proc mpisize-1 bcast the parent_ array
    MPI_Bcast(&parent_[0],n_*sizeof(Int),MPI_BYTE,mpisize-1,aDistExp.comm);


#else

    //allocate a upcxx shared array ?
    upcxx::shared_array<Int> sh_parent;
    upcxx::shared_array<Int> sh_ancstr;
    sh_ancstr.init(n_,n_/mpisize); 
    sh_parent.init(n_,n_/mpisize); 

    //for(int b = 0; b < n_; b+=n_/mpisize){
    //  upcxx::global_ptr<Int> ptr = &sh_parent[b];
    //  if(ptr.where()==iam){
    //    Int * lptr = (Int*)ptr;
    //    Int cnt = min(n_/mpisize,n_ - b);
    //    std::fill(lptr,lptr+cnt,0);
    //  }
    //}

    //MPI_Barrier(aDistExp.comm);
    workteam->barrier();


    //Int * pLocalAncstr = (int64_t*)workloads[iam].raw_ptr();


    Idx fc = (iam)*(n_/mpisize); //0 - based

    //for(Idx locCol = 0; locCol< aDistExp.LocalVertexCount(); locCol++){
    for(Idx locCol = 0; locCol< n_; locCol++){
      //Idx i = fc + locCol; // 0 - based;
      Idx i = locCol;
      Idx node = aOrder.perm[i] - 1; // 0-based (perm is 1 based)
   
//      logfileptr->OFS()<<"node = "<<node+1<<endl; 
      if(node >= fc && node<fc+aDistExp.LocalVertexCount()){     
        sh_parent[i] = 0;
        sh_ancstr[i] = 0;
        //parent_[i] = 0;
        //ancstr[i] = 0;
        //Idx node = fc + locCol; // 0 - based;
        //Idx i = aOrder.invp[node] - 1; // 0-based (perm is 1 based)
        //parent_[i] = 0;
        //ancstr[i] = 0;

        Idx locNode = node - fc;
        Ptr jstrt = aDistExp.colptr[locNode] - aDistExp.baseval; //0-based
        Ptr jstop = aDistExp.colptr[locNode+1] - aDistExp.baseval;//0-based
        for(Ptr j = jstrt; j<jstop; ++j){
          Idx nbr = aDistExp.rowind[j] - aDistExp.baseval; //0-based
          //logfileptr->OFS()<<"   nbr = "<<nbr+1<<"  ";
          nbr = aOrder.invp[nbr] - 1; // 0-based (invp is 1 based);
          //assert(nbr<aOrder.invp.size());
          //logfileptr->OFS()<<"|  nbr = "<<nbr+1<<endl;
          //upper triangular part ?
          if  ( nbr < i ){
            // -------------------------------------------
            // for each nbr, find the root of its current
            // elimination tree.  perform path compression
            // as the subtree is traversed.
            // -------------------------------------------
            //column i (unpermuted) is not the parent of column nbr
            Int anc = sh_ancstr[nbr];
            //Int anc = ancstr[nbr];
            if  ( anc != i+1 ){
              Int break_loop = 0;
              //logfileptr->OFS()<<"path: "<<nbr+1<<" ";
              while(anc >0){
                if  ( anc == i+1 ){
                  break_loop = 1;
                  break;
                }
                Int next = sh_ancstr[nbr];
                sh_ancstr[nbr] = i+1;
                //Int next = ancstr[nbr];
                //ancstr[nbr] = i+1;
                nbr = next - 1;
                anc = sh_ancstr[nbr];
                //anc = ancstr[nbr];
                //logfileptr->OFS()<<nbr+1<<" ";
              }
              //logfileptr->OFS()<<endl;

              // --------------------------------------------
              // now, nbr is the root of the subtree.  make i
              // the parent node of this root.
              // --------------------------------------------
              if(!break_loop){
                sh_parent[nbr] = i + 1; // 1-based
                sh_ancstr[nbr] = i + 1; // 1-based
                //parent[nbr] = i + 1; // 1-based
                //ancstr[nbr] = i + 1; // 1-based
              }
            }
          }
        }
      }
      workteam->barrier();
      //workteam->reduce();
    }

  //if(np>1){
  //  //Do an allgatherv 
  //  //SYMPACK::vector<int> sizes(np,(n_/np)*sizeof(Int));
  //  //sizes.back() = (n_ - (np-1)*(n_/np))*sizeof(Int);
  //  //SYMPACK::vector<int>displs(np,0);
  //  //std::copy(&sizes.front(),&sizes.back(),&displs[1]);
  //  //std::partial_sum(displs.begin(),displs.end(),displs.begin());
  //  //MPI_Allgatherv(&parent_[


  //  //do an allreduce
  //  //logfileptr->OFS()<<"parent before: "<<parent_<<endl;
  //  MPI_Allreduce(MPI_IN_PLACE,&parent_[0],n_,MPI_INT,MPI_MAX,aDistExp.comm);
  //}



//    int b = 0;
//    while(1){
//      Int p = (b)%mpisize;
//      //for(int p = 0; p< mpisize; p++){
//      Idx fc = (b)*(n_/mpisize); //0 - based
//      Idx cnt = n_/mpisize;
//      if(fc+cnt>n_){ cnt = n_ - fc; }
//      upcxx::copy(sh_parent[fc].raw_ptr(),&parent_[fc],cnt);
//      if(fc+cnt==n_){break;}
//      b++;
//    }

    //MPI_Barrier(aDistExp.comm);
    workteam->barrier();
    parent_.resize(n_);
    for(Idx i = 0; i<n_; i++){ parent_[i] = sh_parent[i]; }

    //for(int b = 0; b < n_; b+=n_/mpisize){
    //  upcxx::global_ptr<Int> ptr = &sh_parent[b];
    ////  if(ptr.where()==iam){
    ////    Int * lptr = (Int*)ptr;
    //  Int cnt = min(n_/mpisize,n_ - b);
    //  upcxx::copy(ptr,upcxx::global_ptr<Int>(&parent_[b]),cnt);
    ////    std::copy(lptr,lptr+cnt,&parent_[b]);
    ////  }
    //}

    //MPI_Barrier(aDistExp.comm);
    workteam->barrier();
#endif

TIMER_STOP(Construct_Etree);

  }




  void ETree::ConstructETree(SparseMatrixStructure & aGlobal, Ordering & aOrder){
    bIsPostOrdered_=false;
    n_ = aGlobal.size;

    //Expand to symmetric storage
    aGlobal.ExpandSymmetric();

TIMER_START(Construct_Etree_Classic);
    parent_.resize(n_);
    SetValue(parent_,I_ZERO );


    SYMPACK::vector<Int> ancstr(n_);



    for(Int i = 1; i<=n_; ++i){
            parent_[i-1] = 0;
            ancstr[i-1] = 0;
            Int node = aOrder.perm[i-1];
            Ptr jstrt = aGlobal.expColptr[node-1];
            Ptr jstop = aGlobal.expColptr[node] - 1;
            if  ( jstrt < jstop ){
              for(Ptr j = jstrt; j<=jstop; ++j){
                    Idx nbr = aGlobal.expRowind[j-1];
          //logfileptr->OFS()<<"   nbr = "<<nbr<<"  ";
                    nbr = aOrder.invp[nbr-1];
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

              //logfileptr->OFS()<<"path: "<<nbr<<" ";
                        while(ancstr[nbr-1] >0){
                          if  ( ancstr[nbr-1] == i ){
                            break_loop = 1;
                            break;
                          }
                          Int next = ancstr[nbr-1];
                          ancstr[nbr-1] = i;
                          nbr = next;
                //logfileptr->OFS()<<nbr<<" ";
                        }
              //logfileptr->OFS()<<endl;
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

