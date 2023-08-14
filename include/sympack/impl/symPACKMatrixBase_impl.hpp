#ifndef _SYMPACK_MATRIX_BASE_IMPL_HPP_
#define _SYMPACK_MATRIX_BASE_IMPL_HPP_

#include <sympack/symPACKMatrixBase.hpp>

#include "sympack/Ordering.hpp"

namespace symPACK{
  extern "C" {
    void FORTRAN(ordsup) (int * ordflag, int *  altflag, int *  NEQNS, int *  nofsub, int *  nsuper, 
        int * xsuper, int *  xlindx, int *  lindx , int *  snode , int *  perm  , 
        int * invp  , int *  freeforw, int *  freeback, int *  sforw, int *  sback, 
        int * setseg_forw, int *  setseg_back, int *  nodehead, 
        int * nodeforw, int *  nodeback, 
        int *  setsnode, int *  supperm, int *  mark, int *  set  , int *  compset,
        int *  invp2 , int *  heap                             );

    void FORTRAN(ordsup_ind_tsp_paths2)
      (  int * nadj  , int * neqns , int * nofsub, int * nsuper, int * supsiz,
         int * xsuper, int * xlindx, int * lindx , int * snode , int * xadj  , 
         int * adjncy, int * etpar , int * perm  , int * invp  , int * iflag , 
         int * xskadj, int * sklenf, int * sklenb, int * skadj , int * invp2 , 
         int * link  , int * fstloc, int * sperm , int * fstloc2, int * dist1, 
         int * suppar, int * iwsiz , int * iwork , int * rep             );

    void FORTRAN(ordsup_ind_tsp_paths)
      (  int * nadj  , int * neqns , int * nofsub, int * nsuper, int * supsiz,
         int * xsuper, int * xlindx, int * lindx , int * snode , int * xadj  , 
         int * adjncy, int * etpar , int * perm  , int * invp  , int * iflag , 
         int * xskadj, int * sklenf, int * sklenb, int * skadj , int * invp2 , 
         int * link  , int * fstloc, int * sperm , int * fstloc2, int * dist1, 
         int * suppar, int * iwsiz , int * iwork              );


  }
}


namespace symPACK {

  template <typename T> 
    inline void symPACKMatrixMeta<T>::findSupernodes(ETree& tree, Ordering & aOrder, std::vector<Int> & cc,std::vector<Int> & supMembership, std::vector<Int> & xsuper, Int maxSize ){
      SYMPACK_TIMER_START(FindSupernodes);
      Int size = this->iSize_;
      //TODO: tree order cc supmembership xsuper are all members of the class. no need for argument
      supMembership.resize(size);

      Int nsuper = 1;
      Int supsize = 1;
      supMembership[0] = 1;

      for(Int i =2; i<=size;i++){
        Int prev_parent = tree.PostParent(i-2);
        if(prev_parent == i){
          if(cc[i-2] == cc[i-1]+1 ) {
            if(supsize<maxSize || maxSize==0){
              ++supsize;
              supMembership[i-1] = nsuper;
              continue;
            }
          }
        }

        nsuper++;
        supsize = 1;
        supMembership[i-1] = nsuper;
      }

      xsuper.resize(nsuper+1);
      Int lstsup = nsuper+1;
      for(Int i = size; i>=1;--i){
        Int ksup = supMembership[i-1];
        if(ksup!=lstsup){
          xsuper[lstsup-1] = i + 1; 
        }
        lstsup = ksup;
      }
      xsuper[0]=1;
      SYMPACK_TIMER_STOP(FindSupernodes);
    }

  template <typename T> 
    inline void symPACKMatrixMeta<T>::getLColRowCount(DistSparseMatrixGraph & dgraph, std::vector<Int> & cc, std::vector<Int> & rc){
      scope_timer(q,GetColRowCount_Classic);

      //The tree need to be postordered
      if(!this->ETree_.IsPostOrdered()){
        this->ETree_.PostOrderTree(this->Order_);
      }

      if(this->iam == 0 && (!dgraph.IsExpanded() ) ){
        throw std::logic_error( "DistSparseMatrixGraph must be expanded and permuted\n" );
      }

      dgraph.SetBaseval(1);
      dgraph.SetKeepDiag(1);


      int mpisize;
      MPI_Comm_size(dgraph.GetComm(),&mpisize);

      int mpirank;
      MPI_Comm_rank(dgraph.GetComm(),&mpirank);

      Int size = dgraph.size;

      MPI_Datatype Idxtype;
      MPI_Type_contiguous( sizeof(Idx), MPI_BYTE, &Idxtype );
      MPI_Type_commit(&Idxtype);

      MPI_Datatype Inttype;
      MPI_Type_contiguous( sizeof(Int), MPI_BYTE, &Inttype );
      MPI_Type_commit(&Inttype);

      {
        std::vector<Idx> level(size+1);
        std::vector<Idx> nchild(size+1,0);
        std::vector<Idx> fdesc(size+1);

        std::vector<Idx> storage;
        Int * weight = nullptr;
        Idx * set = nullptr;
        Idx * prvlf = nullptr;
        Idx * prvnbr = nullptr;
        Int * drc = nullptr;
        Idx * pxsup = nullptr;

        if(mpirank==0){
          storage.resize(5*size+1+1);
          weight = (Int*)&storage[0];
          set = &storage[size+1];
          prvlf = &storage[2*size+1];
          prvnbr = &storage[3*size+1];
          drc = (Int*)&storage[4*size+1];
          pxsup = &storage[5*size+1];
          *pxsup = 1;
        }

        level[0] = 0;
        for(Idx k = size; k>=1; --k){
          if(mpirank==0){
            drc[k-1] = 1;
            set[k-1] = k;
            prvlf[k-1] = 0;
            prvnbr[k-1] = 0;
            weight[k] = 1;
          }
          fdesc[k] = k;
          level[k] = level[this->ETree_.PostParent(k-1)] + 1;
          nchild[k] = 0;
        }

        nchild[0] = 0;
        fdesc[0] = 0;
        for(Idx k =1; k<size; ++k){
          Idx parent = this->ETree_.PostParent(k-1);

          ++nchild[parent];
          if(mpirank==0){
            weight[parent] = 0;
          }

            Idx ifdesc = fdesc[k];
            if  ( ifdesc < fdesc[parent] ) {
              fdesc[parent] = ifdesc;
            }
        }




        Idx firstLocCol = dgraph.LocalFirstVertex()-1;

        if(mpirank>0){
          //If something is coming, we can resize          
          MPI_Probe(mpirank-1,mpirank-1,this->graph_.GetComm(),MPI_STATUS_IGNORE);

          storage.resize(5*size+1+1);
          weight = (Int*)&storage[0];
          set = &storage[size+1];
          prvlf = &storage[2*size+1];
          prvnbr = &storage[3*size+1];
          drc = (Int*)&storage[4*size+1];
          pxsup = &storage[5*size+1];
          MPI_Recv(storage.data(),storage.size(),Idxtype,mpirank-1,mpirank-1,this->graph_.GetComm(),MPI_STATUS_IGNORE);
        }

        for(Idx loclownbr = 1; loclownbr<=dgraph.LocalVertexCount(); ++loclownbr){
          Idx lownbr = firstLocCol + loclownbr;

          Int lflag = 0;
          Idx ifdesc = fdesc[lownbr];
          Ptr jstrt = dgraph.colptr[loclownbr-1];
          Ptr jstop = dgraph.colptr[loclownbr] - 1;
          //           -----------------------------------------------
          //           for each ``high neighbor'', hinbr of lownbr ...
          //           -----------------------------------------------

          for(Ptr j = jstrt; j<=jstop;++j){
            Idx hinbr = dgraph.rowind[j-1];
            if  ( hinbr > lownbr )  {
              if  ( ifdesc > prvnbr[hinbr-1] ) {
                //                       -------------------------
                //                       increment weight[lownbr].
                //                       -------------------------
                ++weight[lownbr];
                Idx pleaf = prvlf[hinbr-1];
                //                       -----------------------------------------
                //                       if hinbr has no previous ``low neighbor'' 
                //                       then ...
                //                       -----------------------------------------
                if  ( pleaf == 0 ) {
                  //                           -----------------------------------------
                  //                           ... accumulate lownbr-->hinbr path length 
                  //                               in rowcnt[hinbr].
                  //                           -----------------------------------------
                  drc[hinbr-1] += level[lownbr] - level[hinbr];
                }
                else{
                  //                           -----------------------------------------
                  //                           ... otherwise, lca <-- find[pleaf], which 
                  //                               is the least common ancestor of pleaf 
                  //                               and lownbr.
                  //                               (path halving.)
                  //                           -----------------------------------------
                  Idx last1 = pleaf;
                  Idx last2 = set[last1-1];
                  Idx lca = set[last2-1];
                  while(lca != last2){
                    set[last1-1] = lca;
                    last1 = lca;
                    last2 = set[last1-1];
                    lca = set[last2-1];
                  }
                  //                           -------------------------------------
                  //                           accumulate pleaf-->lca path length in 
                  //                           rowcnt[hinbr].
                  //                           decrement weight(lca).
                  //                           -------------------------------------
                  drc[hinbr-1] += level[lownbr] - level[lca];
                  --weight[lca];
                }
                //                       ----------------------------------------------
                //                       lownbr now becomes ``previous leaf'' of hinbr.
                //                       ----------------------------------------------
                prvlf[hinbr-1] = lownbr;
                lflag = 1;
              }
              //                   --------------------------------------------------
              //                   lownbr now becomes ``previous neighbor'' of hinbr.
              //                   --------------------------------------------------
              prvnbr[hinbr-1] = lownbr;
            }
          }
          //           ----------------------------------------------------
          //           decrement weight ( parent[lownbr] ).
          //           set ( p[lownbr] ) <-- set ( p[lownbr] ) + set[xsup].
          //           ----------------------------------------------------
          Idx parent = this->ETree_.PostParent(lownbr-1);
          --weight[parent];

          if  ( lflag == 1  || nchild[lownbr] >= 2 ) {
            *pxsup = lownbr;
          }
          set[*pxsup-1] = parent;
        }

        if(mpirank<mpisize-1){
          MPI_Send(storage.data(),storage.size(),Idxtype,mpirank+1,mpirank,this->graph_.GetComm());
        }
        else{

          cc.assign(size,0);
          for(Int k = 1; k<=size; ++k){
            Int temp = cc[k-1] + weight[k];
            cc[k-1] = temp;
            Int parent = this->ETree_.PostParent(k-1);
            if  ( parent != 0 ) {
              cc[parent-1] += temp;
            }
          }

          rc.resize(size);
          for(size_t i=0;i<size;i++){ rc[i] = drc[i];}
        }
      }

      //Broadcast to everyone

      if (mpirank<mpisize-1){
        cc.resize(size);
      }

      MPI_Bcast(&cc[0],size,Inttype,mpisize-1,this->fullcomm_);

      rc.resize(size);
      MPI_Bcast(&rc[0],size,Inttype,mpisize-1,this->fullcomm_);



      MPI_Type_free(&Inttype);
      MPI_Type_free(&Idxtype);

    }




  template <typename T> 
    inline void symPACKMatrixMeta<T>::getLColRowCount(SparseMatrixGraph & sgraph, std::vector<Int> & cc, std::vector<Int> & rc){
      scope_timer(q,GetColRowCount_Classic);
      //The tree need to be postordered
      if(!this->ETree_.IsPostOrdered()){
        this->ETree_.PostOrderTree(this->Order_);
      }

      if(this->iam == 0 && (!sgraph.IsExpanded() ) ){
        throw std::logic_error( "SparseMatrixGraph must be expanded\n" );
      }

      sgraph.SetBaseval(1);
      sgraph.SetKeepDiag(1);

      MPI_Datatype Inttype;
      MPI_Type_contiguous( sizeof(Int), MPI_BYTE, &Inttype );
      MPI_Type_commit(&Inttype);

      Int size = sgraph.size;

      if(this->iam==0){
        cc.resize(size);
        rc.resize(size);
        std::vector<Idx> level(size+1);
        std::vector<Int> weight(size+1,0);
        std::vector<Idx> fdesc(size+1);
        std::vector<Idx> nchild(size+1);
        std::vector<Idx> set(size);
        std::vector<Idx> prvlf(size);
        std::vector<Idx> prvnbr(size);

        Idx xsup = 1;
        level[0] = 0;
        for(Idx k = size; k>=1; --k){
          rc[k-1] = 1;
          cc[k-1] = 0;
          set[k-1] = k;
          prvlf[k-1] = 0;
          prvnbr[k-1] = 0;
          level[k] = level[this->ETree_.PostParent(k-1)] + 1;
          weight[k] = 1;
          fdesc[k] = k;
          nchild[k] = 0;
        }

        nchild[0] = 0;
        fdesc[0] = 0;
        for(Idx k =1; k<size; ++k){
          Idx parent = this->ETree_.PostParent(k-1);
          weight[parent] = 0;
          ++nchild[parent];
          Idx ifdesc = fdesc[k];
          if  ( ifdesc < fdesc[parent] ) {
            fdesc[parent] = ifdesc;
          }
        }

        for(Idx lownbr = 1; lownbr<=size; ++lownbr){
          Int lflag = 0;
          Idx ifdesc = fdesc[lownbr];
          Idx oldnbr = this->Order_.perm[lownbr-1];
          Ptr jstrt = sgraph.colptr[oldnbr-1];
          Ptr jstop = sgraph.colptr[oldnbr] - 1;

          //           -----------------------------------------------
          //           for each ``high neighbor'', hinbr of lownbr ...
          //           -----------------------------------------------
          for(Ptr j = jstrt; j<=jstop;++j){
            Idx hinbr = sgraph.rowind[j-1];
            hinbr = this->Order_.invp[hinbr-1];
            if  ( hinbr > lownbr )  {
              if  ( ifdesc > prvnbr[hinbr-1] ) {
                //                       -------------------------
                //                       increment weight[lownbr].
                //                       -------------------------
                ++weight[lownbr];
                Idx pleaf = prvlf[hinbr-1];
                //                       -----------------------------------------
                //                       if hinbr has no previous ``low neighbor'' 
                //                       then ...
                //                       -----------------------------------------
                if  ( pleaf == 0 ) {
                  //                           -----------------------------------------
                  //                           ... accumulate lownbr-->hinbr path length 
                  //                               in rowcnt[hinbr].
                  //                           -----------------------------------------
                  rc[hinbr-1] += level[lownbr] - level[hinbr];
                }
                else{
                  //                           -----------------------------------------
                  //                           ... otherwise, lca <-- find[pleaf], which 
                  //                               is the least common ancestor of pleaf 
                  //                               and lownbr.
                  //                               (path halving.)
                  //                           -----------------------------------------
                  Idx last1 = pleaf;
                  Idx last2 = set[last1-1];
                  Idx lca = set[last2-1];
                  while(lca != last2){
                    set[last1-1] = lca;
                    last1 = lca;
                    last2 = set[last1-1];
                    lca = set[last2-1];
                  }
                  //                           -------------------------------------
                  //                           accumulate pleaf-->lca path length in 
                  //                           rowcnt[hinbr].
                  //                           decrement weight(lca).
                  //                           -------------------------------------
                  rc[hinbr-1] += level[lownbr] - level[lca];
                  --weight[lca];
                }
                //                       ----------------------------------------------
                //                       lownbr now becomes ``previous leaf'' of hinbr.
                //                       ----------------------------------------------
                prvlf[hinbr-1] = lownbr;
                lflag = 1;
              }
              //                   --------------------------------------------------
              //                   lownbr now becomes ``previous neighbor'' of hinbr.
              //                   --------------------------------------------------
              prvnbr[hinbr-1] = lownbr;
            }
          }
          //           ----------------------------------------------------
          //           decrement weight ( parent[lownbr] ).
          //           set ( p[lownbr] ) <-- set ( p[lownbr] ) + set[xsup].
          //           ----------------------------------------------------
          Idx parent = this->ETree_.PostParent(lownbr-1);
          --weight[parent];

          if  ( lflag == 1  || nchild[lownbr] >= 2 ) {
            xsup = lownbr;
          }
          set[xsup-1] = parent;
        }

        for(Int k = 1; k<=size; ++k){
          Int temp = cc[k-1] + weight[k];
          cc[k-1] = temp;
          Int parent = this->ETree_.PostParent(k-1);
          if  ( parent != 0 ) {
            cc[parent-1] += temp;
          }
        }
      }

      if(this->iam!=0){
        cc.resize(size);
        rc.resize(size);
      }
      //Broadcast to everyone 
      MPI_Bcast(&cc[0],size,Inttype,0,this->fullcomm_);
      MPI_Bcast(&rc[0],size,Inttype,0,this->fullcomm_);

      MPI_Type_free(&Inttype);
    }





  template <typename T> 
    inline void symPACKMatrixMeta<T>::relaxSupernodes(ETree& tree, std::vector<Int> & cc,std::vector<Int> & supMembership, std::vector<Int> & xsuper, RelaxationParameters & params ){
      //todo tree cc supmembership xsuper and relax params are members, no need for arguments
      Int nsuper = xsuper.size()-1;

      DisjointSet sets;
      sets.Initialize(nsuper);
      std::vector<Int> ncols(nsuper);
      std::vector<Int> zeros(nsuper);
      std::vector<Int> newCC(nsuper);
      for(Int ksup=nsuper;ksup>=1;--ksup){
        Int cset = sets.makeSet(ksup);
        sets.Root(cset-1)=ksup;

        Int fstcol = xsuper[ksup-1];
        Int lstcol = xsuper[ksup]-1;
        Int width = lstcol - fstcol +1;
        Int length = cc[fstcol-1];
        ncols[ksup-1] = width;
        zeros[ksup-1] = 0;
        newCC[ksup-1] = length;
      }



      for(Int ksup=nsuper;ksup>=1;--ksup){
        Int fstcol = xsuper[ksup-1];
        Int lstcol = xsuper[ksup]-1;
        Int width = ncols[ksup-1];
        Int length = cc[fstcol-1];

        Int parent_fc = tree.PostParent(lstcol-1);
        if(parent_fc!=0){
          Int parent_snode = supMembership[parent_fc-1];
          Int pset = sets.find(parent_snode);
          parent_snode = sets.Root(pset-1);

          bool merge = (parent_snode == ksup+1);


          if(merge){
            Int parent_width = ncols[parent_snode-1];

            Int parent_fc = xsuper[parent_snode-1];
            Int totzeros = zeros[parent_snode-1];
            Int merged_snode_size = width + parent_width;

            //TODO rename nrelax0 as maxSnodeSize
            //TODO rename nrelax1 as ???
            //TODO rename nrelax2 as ???
            //TODO rename zrelax0 as ???
            //TODO rename zrelax1 as ???
            //TODO rename zrelax2 as ???
            merge = false;

            //Supernode is extremely small > merge it
            //TODO TRY THIS
            if( (merged_snode_size <=params.maxSize || params.maxSize==0 ) && params.nrelax0>0){
              if(merged_snode_size <= params.nrelax0){
                merge = true;
              }
              else if(merged_snode_size <=params.maxSize || params.maxSize==0){
                double nnzchild = cc[fstcol-1];
                double nnzparent = cc[parent_fc-1];
                double xnewzeros = width * (nnzparent + width  - nnzchild);

                //The merge is not creating extra fill=in, proceed safely
                if(xnewzeros == 0){
                  merge = true;
                }
                else{
                  //candidate merged supernode characteristics
                  double xtotzeros = (double)totzeros + xnewzeros;
                  double xmerged_snode_size = (double) merged_snode_size;
                  //new number of nz
                  double xtotsize = (xmerged_snode_size * (xmerged_snode_size+1)/2) + xmerged_snode_size * (nnzparent - parent_width);
                  //percentage of explicit zeros
                  double z = xtotzeros / xtotsize;

                  Int totsize = (merged_snode_size * (merged_snode_size+1)/2) + merged_snode_size * ((Int)nnzparent - parent_width);
                  totzeros += (Int)xnewzeros;

                  //check that we will not have Integer overflow issues with the Ptr type
                  Ptr ptr_max = std::numeric_limits<Ptr>::max();
                  if (xtotsize * sizeof(double) < ptr_max){
                    if (merged_snode_size <= params.nrelax1 && z < params.zrelax0){
                      merge = true;
                    }
                    else if (merged_snode_size <= params.nrelax2 && z < params.zrelax1){
                      merge = true;
                    }
                    else if (z<params.zrelax2){
                      merge = true;
                    }
                  }
                }

              }
            }

            //Merge the two supernodes
            if(merge){
              ncols[ksup-1] += ncols[parent_snode-1]; 
              zeros[ksup-1] = totzeros;
              newCC[ksup-1] = width + newCC[parent_snode-1];
              sets.Union(ksup,parent_snode,ksup);
            }
          } 

        }
      }

      std::vector<Int> relXSuper(nsuper+1);
      Int nrSuper = 0;
      for(Int ksup=1;ksup<=nsuper;++ksup){
        Int kset = sets.find(ksup);
        if(ksup == sets.Root(kset-1)){
          Int fstcol = xsuper[ksup-1];
          relXSuper[nrSuper] = fstcol;
          newCC[nrSuper] = newCC[ksup-1];
          ++nrSuper;
        }
      }
      relXSuper[nrSuper] = xsuper[nsuper];
      relXSuper.resize(nrSuper+1);

      for(Int ksup=1;ksup<=nrSuper;++ksup){
        Int fstcol = relXSuper[ksup-1];
        Int lstcol = relXSuper[ksup]-1;
        for(Int col = fstcol; col<=lstcol;++col){
          supMembership[col-1] = ksup;
          cc[col-1] = newCC[ksup-1] - (col-fstcol);
        }
      }


      xsuper = relXSuper;
    }


  template <typename T> 
    inline void symPACKMatrixMeta<T>::symbolicFactorizationRelaxedDist(std::vector<Int> & cc){
      scope_timer(a,SymbolicFactorization);
      Int size = this->iSize_;
      ETree& tree = this->ETree_;
      Ordering & aOrder = this->Order_;
      DistSparseMatrixGraph & graph = this->graph_;
      std::vector<Int> & xsuper = this->Xsuper_;
      std::vector<Int> & SupMembership = this->SupMembership_;
      PtrVec & xlindx = this->locXlindx_;
      IdxVec & lindx = this->locLindx_;
      MPI_Comm & comm = this->graph_.comm;

      //permute the graph
      DistSparseMatrixGraph pGraph = graph;

      {
        double tstart = get_time();

        //recompute vertexDist based on XsuperDist
        std::vector<Idx> newVertexDist;
        newVertexDist.resize(this->all_np+1,0);
        newVertexDist[this->all_np] = pGraph.size+1;
        for(int p = 0; p < this->all_np; p++){
          Int S = this->XsuperDist_[p];
          newVertexDist[p] = this->Xsuper_[S-1];
        }

        pGraph.Permute(&this->Order_.invp[0],&newVertexDist[0]);
        double tstop = get_time();
        logfileptr->OFS()<<"Permute time: "<<tstop-tstart<<std::endl;
      }
      Int nsuper = xsuper.size()-1;
      std::list<MPI_Request> mpirequests;

      Ptr nzbeg = 0;
      //nzend points to the last used slot in lindx
      Ptr nzend = 0;

      //tail is the end of list indicator (in rchlnk, not mrglnk)
      Idx tail = size +1;
      Idx head = 0;

      //Array of length nsuper containing the children of 
      //each supernode as a linked list
      std::vector<Idx> mrglnk(nsuper,0);

      //Array of length n+1 containing the current linked list 
      //of merged indices (the "reach" set)
      std::vector<Idx> rchlnk(size+1);

      //Array of length n used to mark indices as they are introduced
      // into each supernode's index set
      std::vector<Int> marker(size,0);

      Int nsuperLocal = this->XsuperDist_[this->iam+1]-this->XsuperDist_[this->iam];
      Int firstSnode = this->XsuperDist_[this->iam];
      Int lastSnode = firstSnode + nsuperLocal-1;
      xlindx.resize(nsuperLocal+1);


      Int firstColumn = this->Xsuper_[firstSnode-1];
      Int numColumnsLocal = this->Xsuper_[lastSnode] - firstColumn;


      //Compute the sum of the column count and resize lindx accordingly
      //nofsub will be the local nnz now
      Ptr nofsub = 0;
      for(Int ksup = 1; ksup<=nsuper; ++ksup){
        if(ksup>=firstSnode && ksup<=lastSnode){
          Int fstcol = xsuper[ksup-1];
          nofsub+=cc[fstcol-1];
        }
      }
      lindx.resize(nofsub);

      Ptr point = 1;
      for(Int ksup = 1; ksup<=nsuper; ++ksup){
        if(ksup>=firstSnode && ksup<=lastSnode){
          Int locSnode = ksup - firstSnode +1;
          Int fstcol = xsuper[ksup-1];
          xlindx[locSnode-1] = point;
          point += cc[fstcol-1]; 
        }
      } 
      xlindx[nsuperLocal] = point;

      std::vector<Ptr> recvXlindx;
      std::vector<Idx> recvLindx;

      if(this->iam>0){
        //build the mrglnk array


        for(Int ksup = 1; ksup<firstSnode; ++ksup){
          Int fstcol = xsuper[ksup-1];
          Int lstcol = xsuper[ksup]-1;
          Int width = lstcol - fstcol +1;
          Int length = cc[fstcol-1];

          //if ksup has a parent, insert ksup into its parent's 
          //"merge" list.
          if(length > width){
            Idx pcol = tree.PostParent(fstcol+width-1-1);  
            Int psup = SupMembership[pcol-1];
            mrglnk[ksup-1] = mrglnk[psup-1];
            mrglnk[psup-1] = ksup;
          }
        }

      }


      point = 1;
      for(Int locksup = 1; locksup<=nsuperLocal; ++locksup){
        Int ksup = locksup + firstSnode - 1;
        Int fstcol = xsuper[ksup-1];
        Int lstcol = xsuper[ksup]-1;
        Int width = lstcol - fstcol +1;
        Int length = cc[fstcol-1];
        Ptr knz = 0;
        rchlnk[head] = tail;

        //If ksup has children in the supernodal e-tree
        Int jsup = mrglnk[ksup-1];
        if(jsup>0){
          //copy the indices of the first child jsup into 
          //the linked list, and mark each with the value 
          //ksup.
          Int parentJ = ksup;
          do{
            Int jwidth = xsuper[jsup]-xsuper[jsup-1];

            Ptr * jxlindx = nullptr;
            Idx * jlindx = nullptr;
            Int locjsup = -1;
            if(jsup>=firstSnode && jsup<=lastSnode){
              locjsup = jsup - firstSnode +1;
              jxlindx = &xlindx[0];
              jlindx = &lindx[0];
            }
            else {
              MPI_Status status;
              recvLindx.resize(size);
              //receive jsup lindx
              Int psrc = 0; for(psrc = 0; psrc<this->iam;psrc++){ if(this->XsuperDist_[psrc]<=jsup && jsup<this->XsuperDist_[psrc+1]){ break; } }
              MPI_Recv(&recvLindx[0],recvLindx.size()*sizeof(Idx),MPI_BYTE,psrc,jsup+this->all_np,comm,&status);
              //get actual number of received elements
              int count = 0;
              MPI_Get_count(&status,MPI_BYTE,&count);
              count/=sizeof(Idx);

              //compute jsup xlindx
              recvXlindx.resize(2);
              recvXlindx[0] = 1;
              recvXlindx[1] = count +1; 

              locjsup = 1;
              jxlindx = &recvXlindx[0];
              jlindx = &recvLindx[0];
            }

            Ptr jnzbeg = jxlindx[locjsup-1] + jwidth;
            Ptr jnzend = jxlindx[locjsup] -1;
            if(parentJ == ksup){
              for(Ptr jptr = jnzend; jptr>=jnzbeg; --jptr){
                Idx newi = jlindx[jptr-1];
                ++knz;
                marker[newi-1] = ksup;
                rchlnk[newi] = rchlnk[head];
                rchlnk[head] = newi;
              }
            }
            else{
              Int nexti = head;
              for(Ptr jptr = jnzbeg; jptr<=jnzend; ++jptr){
                Idx newi = jlindx[jptr-1];
                Idx i;
                do{
                  i = nexti;
                  nexti = rchlnk[i];
                }while(newi > nexti);

                if(newi < nexti){
                  ++knz;
                  rchlnk[i] = newi;
                  rchlnk[newi] = nexti;
                  marker[newi-1] = ksup;
                  nexti = newi;
                }
              }
            }

            parentJ = jsup;
            jsup = mrglnk[jsup-1];
          } while(jsup!=0 && knz < length);


          //TODO do better than this:need to ainline void sending unnecessary data
          //receive the speculative sends
          jsup = mrglnk[ksup-1];
          //get the next element of the list
          jsup = mrglnk[jsup-1];
          while(jsup>0){
            Int psrc = 0; for(psrc = 0; psrc<this->iam;psrc++){ if(this->XsuperDist_[psrc]<=jsup && jsup<this->XsuperDist_[psrc+1]){ break; } }
            if(psrc!=this->iam){
              MPI_Status status;
              MPI_Request request;
              MPI_Irecv(&recvLindx[0],recvLindx.size()*sizeof(Idx),MPI_BYTE,psrc,jsup+this->all_np,comm,&request);
              MPI_Cancel(&request);
            }
            jsup = mrglnk[jsup-1];
          }
        }

        //structure of a(*,fstcol) has not been examined yet.  
        //"sort" its structure into the linked list,
        //inserting only those indices not already in the
        //list.
        if(knz < length){
          //loop on local columns instead for LOCAL EXPANDED structure
          for(Int row = fstcol; row<=lstcol; ++row){
            Idx newi = row;
            if(newi > fstcol && marker[newi-1] != ksup){
              //position and insert newi in list and
              // mark it with kcol
              Idx nexti = head;
              Idx i;
              do{
                i = nexti;
                nexti = rchlnk[i];
              }while(newi > nexti);
              ++knz;
              rchlnk[i] = newi;
              rchlnk[newi] = nexti;
              marker[newi-1] = ksup;
            }
          }

            for(Int col = fstcol; col<=lstcol; ++col){
              Int local_col = col - firstColumn + 1;
              Ptr knzbeg = pGraph.colptr[local_col-1];
              Ptr knzend = pGraph.colptr[local_col]-1;
              for(Ptr kptr = knzbeg; kptr<=knzend;++kptr){
                Idx newi = pGraph.rowind[kptr-1];

                if(newi > fstcol && marker[newi-1] != ksup){
                  //position and insert newi in list and
                  // mark it with kcol
                  Idx nexti = head;
                  Idx i;
                  do{
                    i = nexti;
                    nexti = rchlnk[i];
                  }while(newi > nexti);
                  ++knz;
                  rchlnk[i] = newi;
                  rchlnk[newi] = nexti;
                  marker[newi-1] = ksup;
                }
              }

              if(this->options_.relax.nrelax0==0 && this->options_.order_refinement_str == "NO") {
                break;
              }
            }
        } 

        //if ksup has no children, insert fstcol into the linked list.
        if(rchlnk[head] != fstcol){
          rchlnk[fstcol] = rchlnk[head];
          rchlnk[head] = fstcol;
          ++knz;
        }

        {
            Idx i = head;
            for(Int col = fstcol; col<=lstcol; ++col){
              Int local_col = col - firstColumn + 1;
              Ptr knzbeg = pGraph.colptr[local_col-1];
              Ptr knzend = pGraph.colptr[local_col]-1;
              for(Ptr kptr = knzbeg; kptr<=knzend;++kptr){
                Idx newi = pGraph.rowind[kptr-1];
              }
              if(this->options_.relax.nrelax0==0 && this->options_.order_refinement_str == "NO") {
                break;
              }
            }
        }

        bassert(knz == cc[fstcol-1]);


        //copy indices from linked list into lindx(*).
        nzbeg = nzend+1;
        nzend += knz;
        xlindx[locksup] = nzend+1;
        bassert(nzend+1 == xlindx[locksup]);
        Idx i = head;
        for(Ptr kptr = nzbeg; kptr<=nzend;++kptr){
          i = rchlnk[i];
          lindx[kptr-1] = i;
        } 

        //if ksup has a parent, insert ksup into its parent's 
        //"merge" list.
        if(length > width){
          Idx pcol = tree.PostParent(fstcol+width-1-1);  
          Int psup = SupMembership[pcol-1];
          mrglnk[ksup-1] = mrglnk[psup-1];
          mrglnk[psup-1] = ksup;

          //send L asap
          Int pdest = 0; for(pdest = 0; pdest<this->all_np;pdest++){ if(this->XsuperDist_[pdest]<=psup && psup<this->XsuperDist_[pdest+1]){ break; } }
          //if remote
          if(pdest!=this->iam){
            mpirequests.push_back(MPI_REQUEST_NULL);
            MPI_Request & request = mpirequests.back();
            MPI_Isend(&lindx[xlindx[locksup-1]-1],length*sizeof(Idx),MPI_BYTE,pdest,ksup+this->all_np,comm,&request);
          }
        }
      }

      bassert(nzend==0 || this->iam<this->np);
      lindx.resize(nzend);
      MPI_Barrier(comm);
    }




  template <typename T> 
    inline void symPACKMatrixMeta<T>::gatherLStructure(std::vector<Ptr>& xlindx, std::vector<Idx> & lindx){
      //Gather this->locXlindx_ and this->locLindx_
      //get other proc vertex counts
      Idx localVertexCnt = this->locXlindx_.size()-1;
      std::vector<Idx> remoteVertexCnt(this->np,0);
      MPI_Gather(&localVertexCnt,sizeof(localVertexCnt),MPI_BYTE,&remoteVertexCnt[0],sizeof(localVertexCnt),MPI_BYTE,0,this->workcomm_);
      Idx totalVertexCnt = std::accumulate(remoteVertexCnt.begin(),remoteVertexCnt.end(),0,std::plus<Idx>());
      if(this->iam==0){
        xlindx.resize(totalVertexCnt+1);
      }
      //compute receive displacements

    MPI_Datatype type;
    MPI_Type_contiguous( sizeof(Ptr), MPI_BYTE, &type );
    MPI_Type_commit(&type);
    MPI_Datatype typeIdx;
    MPI_Type_contiguous( sizeof(Idx), MPI_BYTE, &typeIdx );
    MPI_Type_commit(&typeIdx);

      std::vector<int> rsizes(this->np,0);
      for(int p = 0; p<this->np;p++){rsizes[p] = (int)remoteVertexCnt[p];}
      std::vector<int> rdispls(this->np+1,0);
      rdispls[0]=0;
      std::partial_sum(rsizes.begin(),rsizes.end(),rdispls.begin()+1);
      MPI_Gatherv(&this->locXlindx_[0],localVertexCnt,type,&xlindx[0],&rsizes[0],&rdispls[0],type,0,this->workcomm_);


      Ptr localEdgeCnt = this->locLindx_.size();
      std::vector<Ptr> remoteEdgeCnt(this->np,0);
      MPI_Gather(&localEdgeCnt,sizeof(localEdgeCnt),MPI_BYTE,&remoteEdgeCnt[0],sizeof(localEdgeCnt),MPI_BYTE,0,this->workcomm_);
      Ptr totalEdgeCnt = std::accumulate(remoteEdgeCnt.begin(),remoteEdgeCnt.end(),0,std::plus<Ptr>());

      //fix xlindx
      if(this->iam==0){

        Idx pos = remoteVertexCnt[0];
        Ptr offset = 0;
        for(int p=1;p<this->np;p++){
          offset+=remoteEdgeCnt[p-1]; 
          for(Idx lv = 0; lv < remoteVertexCnt[p]; lv++){
            xlindx[pos++] += offset;
          }
        }
        xlindx.back()=totalEdgeCnt + 1;

        lindx.resize(totalEdgeCnt);
      }

      //compute receive displacements
      rsizes.assign(this->np,0);
      for(int p = 0; p<this->np;p++){rsizes[p] = (int)remoteEdgeCnt[p];}
      rdispls[0]=0;
      std::partial_sum(rsizes.begin(),rsizes.end(),rdispls.begin()+1);
      MPI_Gatherv(&this->locLindx_[0],localEdgeCnt,typeIdx,&lindx[0],&rsizes[0],&rdispls[0],typeIdx,0,this->workcomm_);
    MPI_Type_free(&typeIdx);
    MPI_Type_free(&type);
    }




  template <typename T> 
    inline void symPACKMatrixMeta<T>::refineSupernodes(int ordflag,int altflag,DistSparseMatrix<T> * pMat){
      ETree& tree = this->ETree_;
      Ordering & aOrder = this->Order_;
      std::vector<Int> & supMembership = this->SupMembership_; 
      std::vector<Int> & xsuper = this->Xsuper_; 

      std::vector<int> ixlindx;
      std::vector<int> ilindx;

      std::vector<int>  new_invp;

      //Gather this->locXlindx_ and this->locLindx_
      {
        std::vector<Ptr> xlindx;
        std::vector<Idx> lindx;
        this->gatherLStructure(xlindx, lindx);

        if(this->iam==0){
          ixlindx.resize(xlindx.size());
          for(int i = 0;i<xlindx.size();i++){
            ixlindx[i] = xlindx[i];
          }
          ilindx.resize(lindx.size());
          for(int i = 0;i<lindx.size();i++){
            ilindx[i] = lindx[i];
          }
        }
      }

      if(this->iam==0){
        int neqns = this->iSize_;
        int nofsub =ilindx.size();
        int nsuper = xsuper.size()-1;

        std::vector<int>  freeforw(neqns,0);
        std::vector<int>  freeback(neqns,0);
        std::vector<int>  sforw(neqns,0);
        std::vector<int>  sback(neqns,0);
        std::vector<int> setseg_forw(neqns,0);
        std::vector<int>  setseg_back(neqns,0);
        std::vector<int>  nodehead(neqns,0);
        std::vector<int> nodeforw(neqns,0);
        std::vector<int>  nodeback(neqns,0);
        std::vector<int>  setsnode(neqns,0);
        std::vector<int>  supperm(nsuper,0);
        std::vector<int>  mark(neqns+1,0);
        std::vector<int>  set (neqns,0);
        std::vector<int>  compset(neqns,0);
        std::vector<int>  invp2(neqns,0);
        std::vector<int>  heap (2*nsuper,0);

        new_invp.assign(neqns,0);


        std::vector<int>  new_perm(neqns,0);
        std::iota(new_invp.begin(),new_invp.end(),1);
        std::iota(new_perm.begin(),new_perm.end(),1);

        FORTRAN(ordsup) (
            &ordflag, &altflag, &neqns , &nofsub, &nsuper, 
            &xsuper[0], &ixlindx[0], &ilindx[0], &supMembership[0], 
            &new_perm[0], &new_invp[0], 
            &freeforw[0], &freeback[0], &sforw[0], &sback[0], 
            &setseg_forw[0], &setseg_back[0], &nodehead[0], 
            &nodeforw[0], &nodeback[0], 
            &setsnode[0], &supperm[0], &mark[0], &set[0], &compset[0],
            &invp2[0], &heap[0]);

        this->Order_.Compose(new_invp);

      }

      // broadcast invp
      Int N = aOrder.invp.size();
      MPI_Bcast(&aOrder.invp[0],N*sizeof(int),MPI_BYTE,0,this->fullcomm_);
      MPI_Bcast(&aOrder.perm[0],N*sizeof(int),MPI_BYTE,0,this->fullcomm_);
    }


}

#endif // _SYMPACK_MATRIX_BASE_IMPL_HPP_
