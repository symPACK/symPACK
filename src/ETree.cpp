#include "sympack/ETree.hpp"
#include "sympack/utility.hpp"


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

    //Compute the parent permutation and update postNumber_
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

      Ptr jstrt = aDistExp.colptr[locCol]; //1-based
      Ptr jstop = aDistExp.colptr[locCol+1] -1;//1-based
      if (jstrt<jstop){
        for(Ptr j = jstrt; j<=jstop; ++j){
          Idx nbr = aDistExp.rowind[j-1]; //1-based
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
      MPI_Send(&parent_[0],n_,Inttype,mpirank+1,mpirank,aDistExp.GetComm());
      MPI_Send(&ancstr[0],n_ ,Inttype,mpirank+1,mpirank,aDistExp.GetComm());
    }

    //Now proc mpisize-1 bcast the parent_ array
    MPI_Bcast(&parent_[0],n_,Inttype,mpisize-1,aDistExp.GetComm());


    MPI_Type_free( &Inttype );

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

        Int node = aOrder.perm[i-1];
        Ptr jstrt = sgraph.colptr[node-1];
        Ptr jstop = sgraph.colptr[node] - 1;
        if  ( jstrt < jstop ){
          for(Ptr j = jstrt; j<=jstop; ++j){
            Idx nbr = sgraph.rowind[j-1];
            nbr = aOrder.invp[nbr-1];
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

