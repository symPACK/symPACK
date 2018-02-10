/*
   "symPACK" Copyright (c) 2016, The Regents of the University of California,
   through Lawrence Berkeley National Laboratory (subject to receipt of any
   required approvals from the U.S. Dept. of Energy).  All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

   (1) Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

   (2) Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

   (3) Neither the name of the University of California, Lawrence Berkeley
   National Laboratory, U.S. Dept. of Energy nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
   FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

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
#ifndef _SYMPACK_MATRIX2D_DECL_HPP_
#define _SYMPACK_MATRIX2D_DECL_HPP_

#include "sympack/Environment.hpp"
#include "sympack/symPACKMatrix.hpp"




namespace symPACK{


  template <typename colptr_t, typename rowind_t, typename T> 
    class symPACKMatrix2D: public symPACKMatrixMeta<T>{
    public:

    class block_t {
      public:
        std::vector<T> nzval;
    };

    class blockCell_t :public std::vector<block_t> {
      using intrank_t = int;
      public:
        Int i;
        Int j;
        intrank_t owner;
        blockCell_t():std::vector<block_t>(),i(0),j(0),owner(0) {}
    };

    //using blockCell_t = std::vector<block_t>;

    protected:
//      std::vector<blockCell_t> cells_;
      std::unordered_map<Int, blockCell_t> cells_;
      //TODO ideally, this should be DistSparseMatrixGraph<colptr_t,rowind_t>
#ifndef NO_MPI
#endif
    public:

      symPACKMatrix2D();
      ~symPACKMatrix2D();



      void Init(symPACKOptions & options );
      void Init(DistSparseMatrix<T> & pMat, symPACKOptions & options );

      void SymbolicFactorization(DistSparseMatrix<T> & pMat);
      void DistributeMatrix(DistSparseMatrix<T> & pMat);




  };

  template <typename colptr_t, typename rowind_t, typename T>
  symPACKMatrix2D<colptr_t,rowind_t,T>::symPACKMatrix2D():
    symPACKMatrixMeta<T>()
  {
#ifndef NO_MPI
#endif
  }

  template <typename colptr_t, typename rowind_t, typename T>
  symPACKMatrix2D<colptr_t,rowind_t,T>::~symPACKMatrix2D(){
  }

  template <typename colptr_t, typename rowind_t, typename T>
   void symPACKMatrix2D<colptr_t,rowind_t,T>::Init(symPACKOptions & options ){
    scope_timer(a,symPACKMatrix2D::Init);
    this->options_ = options;
    logfileptr->verbose = this->options_.verbose>0;

#ifndef NO_MPI
    this->all_np = 0;
    MPI_Comm_size(this->options_.MPIcomm,&this->all_np);
    MPI_Comm_rank(this->options_.MPIcomm,&this->iam);
#endif

#ifdef NEW_UPCXX
    this->iam = upcxx::rank_me();
    this->all_np = upcxx::rank_n();
#endif

    this->np = this->options_.used_procs(this->all_np);


#ifndef NO_MPI
    MPI_Comm_split(this->options_.MPIcomm,this->iam<this->np,this->iam,&this->workcomm_);

    //do another split to contain P0 and all the non working processors
    if(this->all_np!=this->np){
      this->non_workcomm_ = MPI_COMM_NULL;
      MPI_Comm_split(this->options_.MPIcomm,this->iam==0||this->iam>=this->np,this->iam==0?0:this->iam-this->np+1,&this->non_workcomm_);
    }
#else
//    upcxx::team_all.split(this->iam<this->np,new_rank, this->team_);
#endif

   } 

  template <typename colptr_t, typename rowind_t, typename T>
   void symPACKMatrix2D<colptr_t,rowind_t,T>::Init(DistSparseMatrix<T> & pMat,symPACKOptions & options ){
   } 

  template <typename colptr_t, typename rowind_t, typename T>
    void symPACKMatrix2D<colptr_t,rowind_t,T>::SymbolicFactorization(DistSparseMatrix<T> & pMat ){
      scope_timer(a,symPACKMatrix2D::SymbolicFactorization);

      //This has to be declared here to be able to debug ...
      std::vector<int, Mallocator<int> > xadj;
      std::vector<int, Mallocator<int> > adj;
      Idx row;
      int fc,lc,colbeg,colend,col;
      Ptr supbeg,supend,rowidx;
      Int I;

#ifndef NO_MPI
      if(this->fullcomm_!=MPI_COMM_NULL){
        MPI_Comm_free(&this->fullcomm_);
      }
      MPI_Comm_dup(pMat.comm,&this->fullcomm_);
#endif

      this->iSize_ = pMat.size;

      this->graph_ = pMat.GetLocalGraph();
      this->graph_.SetBaseval(1);
      this->graph_.SetSorted(1);
      this->graph_.ExpandSymmetric();
      logfileptr->OFS()<<"Matrix structure expanded"<<std::endl;

      SparseMatrixGraph * sgraph = NULL;
      {
        {
          double timeSta = get_time();
          scope_timer(c,symPACKMatrix2D::ordering);

          if(this->options_.orderingStr=="MMD"){
            this->options_.ordering = symPACK::MMD;
          }
          else if(this->options_.orderingStr=="RCM"){
            this->options_.ordering = symPACK::RCM;
          }
          else if(this->options_.orderingStr=="AMD"){
            this->options_.ordering = symPACK::AMD;
          }
#ifdef USE_METIS
          else if(this->options_.orderingStr=="METIS"){
            this->options_.ordering = symPACK::METIS;
          }
#endif
#ifdef USE_SCOTCH
          else if(this->options_.orderingStr=="SCOTCH"){
            this->options_.ordering = symPACK::SCOTCH;
          }
#endif
#ifdef USE_PARMETIS
          else if(this->options_.orderingStr=="PARMETIS"){
            this->options_.ordering = symPACK::PARMETIS;
          }
#endif
#ifdef USE_PTSCOTCH
          else if(this->options_.orderingStr=="PTSCOTCH"){
            this->options_.ordering = symPACK::PTSCOTCH;
          }
#endif
          else if(this->options_.orderingStr=="NATURAL"){
            this->options_.ordering = symPACK::NATURAL;
          }
          else if(this->options_.orderingStr=="USER"){
            if(this->options_.perm ==NULL){
              throw std::logic_error( "When using USER, symPACKOptions.perm must be provided.\n" );
            }
            else{
              this->options_.ordering = symPACK::USER;
            }
          }
          else{
            std::stringstream sstr;
            sstr<<"This ordering method is not supported by symPACK. Valid options are:";
            sstr<<"NATURAL MMD AMD RCM NDBOX NDGRID USER ";
#ifdef USE_SCOTCH
            sstr<<"SCOTCH ";
#endif
#ifdef USE_PTSCOTCH
            sstr<<"PTSCOTCH ";
#endif
#ifdef USE_METIS
            sstr<<"METIS ";
#endif
#ifdef USE_PARMETIS
            sstr<<"PARMETIS ";
#endif
            sstr<<std::endl;
            throw std::logic_error( sstr.str()  );

          }

          this->Order_.NpOrdering = this->options_.NpOrdering;
          switch(this->options_.ordering){
            case MMD:
              {
                sgraph = new SparseMatrixGraph();
                this->graph_.GatherStructure(*sgraph,0);
                sgraph->SetBaseval(1);
                sgraph->SetKeepDiag(0);
                this->Order_.MMD(*sgraph, this->graph_.comm);
              }
              break;

            case RCM:
              {
                sgraph = new SparseMatrixGraph();
                this->graph_.GatherStructure(*sgraph,0);
                sgraph->SetBaseval(1);
                sgraph->SetKeepDiag(0);
                this->Order_.RCM(*sgraph, this->graph_.comm);
              }
              break;

            case AMD:
              {
                sgraph = new SparseMatrixGraph();
                this->graph_.GatherStructure(*sgraph,0);
                sgraph->SetBaseval(1);
                sgraph->SetKeepDiag(0);
                this->Order_.AMD(*sgraph, this->graph_.comm);
              }
              break;

            case USER:
              {
                this->Order_.perm.resize(this->iSize_);
                //find baseval
                auto baseval = std::min_element(&this->options_.perm[0],&this->options_.perm[0]+this->iSize_);
                //now compute everything in 1 based 
                for(int i=0;i<this->Order_.perm.size();i++){this->Order_.perm[i] = this->options_.perm[i] - *baseval + 1; /*1 based*/}
                this->Order_.invp.resize(this->iSize_);
                for(int i=0;i<this->Order_.perm.size();i++){this->Order_.invp[this->Order_.perm[i]-1] = i;} 
              }
              break;

            case NDBOX:
              this->Order_.NDBOX(this->Size(), this->graph_.comm);
              break;

            case NDGRID:
              this->Order_.NDGRID(this->Size(), this->graph_.comm);
              break;

#ifdef USE_SCOTCH
            case SCOTCH:
              {

                sgraph = new SparseMatrixGraph();
                this->graph_.GatherStructure(*sgraph,0);
                sgraph->SetKeepDiag(0);
                this->Order_.SCOTCH(*sgraph, this->graph_.comm);
              }
              break;
#endif
#ifdef USE_METIS
            case METIS:
              {
                sgraph = new SparseMatrixGraph();
                this->graph_.GatherStructure(*sgraph,0);
                sgraph->SetKeepDiag(0);
                this->Order_.METIS(*sgraph, this->graph_.comm);
              }
              break;
#endif
#ifdef USE_PARMETIS
            case PARMETIS:
              {
                this->graph_.SetKeepDiag(0);
                this->Order_.PARMETIS(this->graph_);
              }
              break;
#endif
#ifdef USE_PTSCOTCH
            case PTSCOTCH:
              {
                this->graph_.SetKeepDiag(0);
                this->Order_.PTSCOTCH(this->graph_);
              }
              break;
#endif
            case NATURAL:
              this->Order_.perm.resize(this->iSize_);
              for(int i=0;i<this->Order_.perm.size();i++){this->Order_.perm[i]=i+1;} 
              this->Order_.invp = this->Order_.perm;
              break;
            default:
              {
                std::stringstream sstr;
                sstr<<"This ordering method is not supported by symPACK. Valid options are:";
                sstr<<"MMD AMD RCM NDBOX NDGRID USER ";
#ifdef USE_SCOTCH
                sstr<<"SCOTCH ";
#endif
#ifdef USE_PTSCOTCH
                sstr<<"PTSCOTCH ";
#endif
#ifdef USE_METIS
                sstr<<"METIS ";
#endif
#ifdef USE_PARMETIS
                sstr<<"PARMETIS ";
#endif
                sstr<<std::endl;
                throw std::logic_error( sstr.str()  );
                //do nothing: either natural or user provided ordering
              }
              break;
          }

          //The ordering is available on every processor of the full communicator
          double timeStop = get_time();
          if(this->iam==0 && this->options_.verbose){
            std::cout<<"Ordering time: "<<timeStop - timeSta<<std::endl;
          }
          logfileptr->OFS()<<"Ordering done"<<std::endl;
        }
      }

      if(this->options_.dumpPerm>0){
        logfileptr->OFS()<<"perm = [";
        for (auto i : this->Order_.perm){ 
          logfileptr->OFS()<<i<<" ";
        }
        logfileptr->OFS()<<"]"<<std::endl;
      }

      std::vector<Int> cc;
      std::vector<Int> rc;


      if(sgraph==NULL){ 
        logfileptr->OFS()<<"copying graph"<<std::endl;
        auto graph = this->graph_;
        double timeSta = get_time();
        graph.SetKeepDiag(1);
        graph.SetBaseval(1);
        logfileptr->OFS()<<"graph copied"<<std::endl;
        //Expand to unsymmetric storage
        graph.ExpandSymmetric();
        logfileptr->OFS()<<"graph expanded"<<std::endl;

        auto backinvp = this->Order_.invp;
        graph.Permute(this->Order_.invp.data());
        logfileptr->OFS()<<"graph permuted 1/2"<<std::endl;
        this->ETree_.ConstructETree(graph,this->Order_);
        this->ETree_.PostOrderTree(this->Order_);
        logfileptr->OFS()<<"ETree computed and postordered"<<std::endl;

        std::vector<Int> relinvp;
        this->Order_.GetRelativeInvp(backinvp,relinvp);
        graph.Permute(relinvp.data());
        logfileptr->OFS()<<"graph permuted 2/2"<<std::endl;

        double timeSta_cc = get_time();
        this->getLColRowCount(graph,cc,rc);
        double timeStop_cc = get_time();
        if(this->iam==0 && this->options_.verbose){
          std::cout<<"Column count (distributed) construction time: "<<timeStop_cc - timeSta_cc<<std::endl;
        }

        if(this->options_.ordering != NATURAL){
          this->ETree_.SortChildren(cc,this->Order_);
          if(this->options_.dumpPerm>0){
            logfileptr->OFS()<<"perm = "<<this->Order_.perm<<std::endl;
          }
        }

        double timeStop = get_time();
        if(this->iam==0 && this->options_.verbose){
          std::cout<<"Elimination tree construction time: "<<timeStop - timeSta<<std::endl;
        }
      }
      else
      {
        //gather the graph if necessary to build the elimination tree
        if(sgraph==NULL){ 
          sgraph = new SparseMatrixGraph();
          this->graph_.GatherStructure(*sgraph,0);
        }
        this->graph_.SetKeepDiag(1);
        sgraph->SetBaseval(1);
        sgraph->SetKeepDiag(1);

        double timeSta = get_time();
        this->ETree_.ConstructETree(*sgraph,this->Order_,this->fullcomm_);
        this->ETree_.PostOrderTree(this->Order_);

        double timeSta_cc = get_time();
        this->getLColRowCount(*sgraph,cc,rc);
        double timeStop_cc = get_time();
        if(this->iam==0 && this->options_.verbose){
          std::cout<<"Column count (gather + serial + bcast) construction time: "<<timeStop_cc - timeSta_cc<<std::endl;
        }

        if(this->options_.ordering != NATURAL){
          this->ETree_.SortChildren(cc,this->Order_);
          if(this->options_.dumpPerm>0){
            logfileptr->OFS()<<"perm = "<<this->Order_.perm<<std::endl;
          }
        }

        double timeStop = get_time();
        if(this->iam==0 && this->options_.verbose){
          std::cout<<"Elimination tree construction time: "<<timeStop - timeSta<<std::endl;
        }
      }

      //compute some statistics
      if(this->iam==0 && this->options_.verbose){
        double flops = 0.0;
        int64_t NNZ = 0;
        for(Int i = 0; i<cc.size();++i){
          flops+= (double)pow((double)cc[i],2.0);
          NNZ+=cc[i];
        }
        std::cout<<"Flops: "<<flops<<std::endl;
        std::cout<<"NNZ in L factor: "<<NNZ<<std::endl;
      }

    //get rid of the sequential graph
    if(sgraph!=NULL && this->options_.order_refinement_str.substr(0,4) != "TSPB"){delete sgraph;}

    { 
      double timeSta = get_time();
      this->findSupernodes(this->ETree_,this->Order_,cc,this->SupMembership_,this->Xsuper_,this->options_.relax.maxSize);
      logfileptr->OFS()<<"Supernodes found"<<std::endl;

      if(this->options_.relax.nrelax0>0)
      {
        this->relaxSupernodes(this->ETree_, cc,this->SupMembership_, this->Xsuper_, this->options_.relax );
        logfileptr->OFS()<<"Relaxation done"<<std::endl;
      }

      //modify this->np since it cannot be greater than the number of supernodes
      this->np = std::min(this->np,this->options_.used_procs(this->Xsuper_.size()-1)); 

      if(this->workcomm_!=MPI_COMM_NULL){
        MPI_Comm_free(&this->workcomm_);
      }
      MPI_Comm_split(this->options_.MPIcomm,this->iam<this->np,this->iam,&this->workcomm_);
      this->group_.reset( new RankGroup( this->workcomm_ ) ); 

      //do another split to contain P0 and all the non working processors
      if(this->all_np!=this->np){
        if(this->non_workcomm_!=MPI_COMM_NULL){
          MPI_Comm_free(&this->non_workcomm_);
        }
        this->non_workcomm_ = MPI_COMM_NULL;
        MPI_Comm_split(this->options_.MPIcomm,this->iam==0||this->iam>=this->np,this->iam==0?0:this->iam-this->np+1,&this->non_workcomm_);
      }

      //Compute this->XsuperDist_
      std::vector<Idx> newVertexDist;
      {
        Idx supPerProc = std::max((size_t)1,(this->Xsuper_.size()-1) / this->np);
        this->XsuperDist_.resize(this->all_np+1,0);
        this->XsuperDist_[this->np] = this->Xsuper_.size();
        for(int p = 0; p<this->np; p++){
          this->XsuperDist_[p]= std::min(this->XsuperDist_[this->np],(Int)(p*supPerProc+1));
        }
        for(int p = this->np+1; p<this->all_np; p++){
          this->XsuperDist_[p]= this->XsuperDist_[p-1];
        }
        this->XsuperDist_[this->all_np] = this->Xsuper_.size();


        newVertexDist.resize(this->all_np+1,0);
        newVertexDist[this->all_np] = this->iSize_+1;
        for(int p = 0; p < this->all_np; p++){
          Int S = this->XsuperDist_[p];
          newVertexDist[p] = this->Xsuper_[S-1];
        }
      }

//      //now create the mapping
//      if(this->Mapping_ != NULL){
//        delete this->Mapping_;
//      }
//
//      Int pmapping = this->np;
//      Int pr = (Int)sqrt((double)pmapping);
//      if(this->options_.mappingTypeStr ==  "MODWRAP2DTREE"){
//        this->Mapping_ = new ModWrap2DTreeMapping(pmapping, pr, pr, 1);
//      }
//      else if(this->options_.mappingTypeStr ==  "WRAP2D"){
//        this->Mapping_ = new Wrap2D(pmapping, pr, pr, 1);
//      }
//      else if(this->options_.mappingTypeStr ==  "WRAP2DFORCED"){
//        this->Mapping_ = new Wrap2DForced(pmapping, pr, pr, 1);
//      }
//      else if(this->options_.mappingTypeStr ==  "MODWRAP2DNS"){
//        this->Mapping_ = new Modwrap2DNS(pmapping, pr, pr, 1);
//      }
//      else if(this->options_.mappingTypeStr ==  "ROW2D"){
//        this->Mapping_ = new Row2D(pmapping, pmapping, pmapping, 1);
//      }
//      else if(this->options_.mappingTypeStr ==  "COL2D"){
//        this->Mapping_ = new Col2D(pmapping, pmapping, pmapping, 1);
//      }
//      else{
//        this->Mapping_ = new Modwrap2D(pmapping, pr, pr, 1);
//      }


      double timeStaSymb = get_time();
      this->symbolicFactorizationRelaxedDist(cc);

      double timeStopSymb = get_time();
      if(this->iam==0 && this->options_.verbose){
        std::cout<<"Symbolic factorization time: "<<timeStopSymb - timeStaSymb<<std::endl;
      }
      logfileptr->OFS()<<"Symbfact done"<<std::endl;


      if(this->options_.order_refinement_str != "NO") {
        double timeSta = get_time();
        if(this->options_.order_refinement_str == "SET"){ 
          this->refineSupernodes(3,1,&pMat);
        }
        else if(this->options_.order_refinement_str == "SET10"){ 
          this->refineSupernodes(1,0,&pMat);
        }
        else if(this->options_.order_refinement_str == "SET11"){ 
          this->refineSupernodes(1,1,&pMat);
        }
        else if(this->options_.order_refinement_str == "SET20"){ 
          this->refineSupernodes(2,0,&pMat);
        }
        else if(this->options_.order_refinement_str == "SET21"){ 
          this->refineSupernodes(2,1,&pMat);
        }
        else if(this->options_.order_refinement_str == "SET30"){ 
          this->refineSupernodes(3,0,&pMat);
        }
        else if(this->options_.order_refinement_str == "SET31"){ 
          this->refineSupernodes(3,1,&pMat);
        }
        else if(this->options_.order_refinement_str == "TSP"){
          auto SupETree = this->ETree_.ToSupernodalETree(this->Xsuper_,this->SupMembership_,this->Order_);

          std::vector<Ptr> xlindx;
          std::vector<Idx> lindx;
          this->gatherLStructure(xlindx, lindx);

          if(this->iam==0){
            TSP::SymbolMatrix * symbmtx = TSP::GetPastixSymbolMatrix(this->Xsuper_,this->SupMembership_, xlindx, lindx);
            TSP::Order * psorder = TSP::GetPastixOrder(symbmtx,this->Xsuper_, SupETree, &this->Order_.perm[0], &this->Order_.invp[0]);

            double timeSta = get_time();
            TSP::symbolReordering( symbmtx, psorder, 0, std::numeric_limits<int>::max(), 0 );
            double timeStop = get_time();
            if(this->iam==0 && this->options_.verbose){
              std::cout<<"TSP reordering done in "<<timeStop-timeSta<<std::endl;
            }

            //overwrite order
            for(int i = 0; i < this->Order_.perm.size(); ++i){
              this->Order_.perm[i] = psorder->peritab[i]+1;
              this->Order_.invp[i] = psorder->permtab[i]+1;
            }
          }

          MPI_Bcast(this->Order_.perm.data(),this->Order_.perm.size()*sizeof(Int),MPI_BYTE,0,this->fullcomm_);
          MPI_Bcast(this->Order_.invp.data(),this->Order_.invp.size()*sizeof(Int),MPI_BYTE,0,this->fullcomm_);
        } 
        else if(this->options_.order_refinement_str.substr(0,4) == "TSPB"){

          if(sgraph==NULL){ 
            sgraph = new SparseMatrixGraph();
            this->graph_.GatherStructure(*sgraph,0);
          }

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
            bassert(sgraph!=nullptr);
            int neqns = this->iSize_;

            int nofsub =ilindx.size();
            int nsuper = xsuper.size()-1;


            new_invp.assign(neqns,0);


            std::vector<int>  new_perm(neqns,0);

            for(size_t i =0; i<new_perm.size(); i++){ new_invp[i] = this->Order_.invp[i];}
            for(size_t i =0; i<new_perm.size(); i++){ new_perm[i] = this->Order_.perm[i];}

            int supsiz = 0;
            for(I = 1; I <= nsuper; I++){
              Int fc = this->Xsuper_[I-1];
              Int lc = this->Xsuper_[I]-1;
              supsiz = std::max(supsiz,lc-fc+1);
            }

            sgraph->SetKeepDiag(0);
            xadj.resize(sgraph->colptr.size());
            adj.resize(sgraph->rowind.size());
            int nadj = adj.size();
            for(size_t i = 0; i< sgraph->colptr.size(); i++){ xadj[i] = int(sgraph->colptr[i]); }
            for(size_t i = 0; i< sgraph->rowind.size(); i++){ adj[i] = int(sgraph->rowind[i]); }

            if(sgraph!=NULL){delete sgraph;}

            std::vector<int, Mallocator<int> > etpar(neqns);
            for(int i = 0; i<neqns; i++){ etpar[i] = tree.PostParent(i); }


            std::vector<int, Mallocator<int> > xskadj(neqns+1);
            std::vector<int, Mallocator<int> > sklenf(neqns);
            std::vector<int, Mallocator<int> > sklenb(neqns);
            std::vector<int, Mallocator<int> > skadj(nadj);
            std::vector<int, Mallocator<int> > invp2(neqns);
            std::vector<int, Mallocator<int> > link(neqns);
            std::vector<int, Mallocator<int> > fstloc(nsuper);
            std::vector<int, Mallocator<int> > sperm(nsuper);
            std::vector<int, Mallocator<int> > fstloc2(neqns);
            std::vector<int, Mallocator<int> > dist1((supsiz+1)*(supsiz+1));
            std::vector<int, Mallocator<int> > suppar(nsuper);
            int iwsiz = 8*neqns+3;
            std::vector<int, Mallocator<int> > iwork(iwsiz);

            int iflag = 0;
            double timeSta = get_time();
            if(this->options_.order_refinement_str == "TSPB"){
              std::vector<int, Mallocator<int> > rep(neqns);
              FORTRAN(ordsup_ind_tsp_paths2)
                ( 
                 &nadj, &neqns , &nofsub, &nsuper, &supsiz, 
                 xsuper.data(), ixlindx.data(), ilindx.data(), supMembership.data(), 
                 xadj.data(), adj.data(), 
                 etpar.data(),
                 new_perm.data(), new_invp.data(), 
                 &iflag , 
                 xskadj.data(), sklenf.data(), sklenb.data(), 
                 skadj.data() , invp2.data() , link.data()  , 
                 fstloc.data(), sperm.data() , fstloc2.data(), dist1.data(), 
                 suppar.data(), &iwsiz , iwork.data() , rep.data()             );
            }
            else{
              FORTRAN(ordsup_ind_tsp_paths)
                ( 
                 &nadj, &neqns , &nofsub, &nsuper, &supsiz, 
                 xsuper.data(), ixlindx.data(), ilindx.data(), supMembership.data(), 
                 xadj.data(), adj.data(), 
                 etpar.data(),
                 new_perm.data(), new_invp.data(), 
                 &iflag , 
                 xskadj.data(), sklenf.data(), sklenb.data(), 
                 skadj.data() , invp2.data() , link.data()  , 
                 fstloc.data(), sperm.data() , fstloc2.data(), dist1.data(), 
                 suppar.data(), &iwsiz , iwork.data() );
            }

            double timeStop = get_time();
            if(this->iam==0 && this->options_.verbose){
              std::cout<<"TSPB reordering done in "<<timeStop-timeSta<<std::endl;
            }


            //this->Order_.Compose(new_invp);
            for(size_t i =0; i<new_perm.size(); i++){ this->Order_.invp[i] = new_invp[i];}
            for(size_t i =0; i<new_perm.size(); i++){ this->Order_.perm[i] = new_perm[i];}
          }

          // broadcast invp
          Int N = aOrder.invp.size();
          MPI_Bcast(&aOrder.invp[0],N*sizeof(Int),MPI_BYTE,0,this->fullcomm_);
          MPI_Bcast(&aOrder.perm[0],N*sizeof(Int),MPI_BYTE,0,this->fullcomm_);

        }

        double timeStop = get_time();

        if(this->iam==0 && this->options_.verbose){
          std::cout<<"Supernode reordering done in "<<timeStop-timeSta<<std::endl;
        }

        {
          double timeSta = get_time();
          this->symbolicFactorizationRelaxedDist(cc);
          double timeStop = get_time();
          if(this->iam==0 && this->options_.verbose){
            std::cout<<"Symbolic factorization time: "<<timeStop - timeSta<<std::endl;
          }
          logfileptr->OFS()<<"Symbfact done"<<std::endl;
        }
      }

      double timeStop = get_time();
      if(this->iam==0 && this->options_.verbose){
        std::cout<<"Total symbolic factorization time: "<<timeStop - timeSta<<std::endl;
      }

      //Print statistics
      if(this->options_.print_stats){
        OrderStats stats;
        stats.get(this->Xsuper_, this->XsuperDist_, this->locXlindx_, this->locLindx_, this->fullcomm_);
        if (this->iam==0){
          stats.print();
        }
      }
    }

#ifdef _DEBUG_
    logfileptr->OFS()<<"Membership list is "<<this->SupMembership_<<std::endl;
    logfileptr->OFS()<<"xsuper "<<this->Xsuper_<<std::endl;
#endif


#ifdef _OUTPUT_ETREE_
    logfileptr->OFS()<<"ETree is "<<this->ETree_<<std::endl;
    {
      auto supETree = this->ETree_.ToSupernodalETree(this->Xsuper_,this->SupMembership_,this->Order_);
      logfileptr->OFS()<<"Supernodal ETree is "<<supETree<<std::endl;
      logfileptr->OFS()<<"Xsuper is "<<this->Xsuper_<<std::endl;
    } 
#endif

    //compute a mapping
    // Do a 2D cell-cyclic distribution for now.
    auto nsuper = this->Xsuper_.size()-1;
    auto iam = this->iam;
    auto np = this->np;
    Int npcol = std::floor(std::sqrt(np));
    Int nprow = std::floor(np / npcol);


    auto coord2supidx = [nsuper](Int i, Int j){ return (i-j) + j*(nsuper+1) - (j+1)*(j)/2;  };
    //TODO implement this
    auto idx2coord = [nsuper](Int idx){ return 0; };

#if 1    
    {
      using cell_tuple_t = std::tuple<Idx,Idx>;
      std::list< cell_tuple_t> cellsToSend;

      Int numLocSnode = this->XsuperDist_[this->iam+1]-this->XsuperDist_[this->iam];
      Int firstSnode = this->XsuperDist_[this->iam];
      for(Int locsupno = 1; locsupno<this->locXlindx_.size(); ++locsupno){
        Idx I = locsupno + firstSnode-1;
        Int first_col = this->Xsuper_[I-1];
        Int last_col = this->Xsuper_[I]-1;

        Ptr lfi = this->locXlindx_[locsupno-1];
        Ptr lli = this->locXlindx_[locsupno]-1;

        Int K = -1; 
        Idx K_prevSnode = -1;
        for(Ptr K_sidx = lfi; K_sidx<=lli;K_sidx++){
          Idx K_row = this->locLindx_[K_sidx-1];
          K = this->SupMembership_[K_row-1];

          //Split at boundary or after diagonal block
          if(K!=K_prevSnode){
            if(K>=I){
              cellsToSend.push_back(std::make_tuple(K,I));
            }
          }
          K_prevSnode = K;
        }
      }

      MPI_Datatype type;
      MPI_Type_contiguous( sizeof(cell_tuple_t), MPI_BYTE, &type );
      MPI_Type_commit(&type);

      //then do an allgatherv
      //compute send sizes
      int ssize = cellsToSend.size();

      //Build the contiguous array
      vector<cell_tuple_t> sendbuf;
      sendbuf.reserve(ssize);
      sendbuf.insert(sendbuf.end(),cellsToSend.begin(),cellsToSend.end());

      //allgather receive sizes
      vector<int> rsizes(this->np,0);
      MPI_Allgather(&ssize,sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,this->workcomm_);

      //compute receive displacements
      vector<int> rdispls(this->np+1,0);
      rdispls[0] = 0;
      std::partial_sum(rsizes.begin(),rsizes.end(),&rdispls[1]);


      //Now do the allgatherv
      vector<cell_tuple_t> recvbuf(rdispls.back());
      MPI_Allgatherv(&sendbuf[0],ssize,type,&recvbuf[0],&rsizes[0],&rdispls[0],type,this->workcomm_);
      MPI_Type_free(&type);

      for(auto it = recvbuf.begin();it!=recvbuf.end();it++){
        auto & cur_cell = (*it);
        Int i = std::get<0>(cur_cell);
        Int j = std::get<1>(cur_cell);
        auto pcol = j % npcol;
        auto idx = coord2supidx(i,j);
        auto prow = i % nprow;
        auto p = pcol*nprow+prow;
        cells_[idx].owner = p;
        cells_[idx].i = i;
        cells_[idx].j = j;
      }
    }

#else

    cells_.resize(nsuper*(nsuper-1)/2.0+nsuper);

    //should create cell only if non empty
    for(Int j = 0; j < nsuper; j++){
      auto pcol = j % npcol;
      for(Int i = j; i < nsuper; i++){
        auto idx = coord2supidx(i,j);
        auto prow = i % nprow;
        auto p = pcol*nprow+prow;
        cells_[idx].owner = p;
        cells_[idx].i = i;
        cells_[idx].j = j;
      }
    }
#endif

        if(iam==0){
          std::cout<<"#supernodes = "<<nsuper<<" #cells = "<<cells_.size()<< " which is "<<cells_.size()*sizeof(blockCell_t)<<" bytes"<<std::endl;
        }

    //generate task graph
    {
      auto supETree = this->ETree_.ToSupernodalETree(this->Xsuper_,this->SupMembership_,this->Order_);

      //we will need to communicate if only partial xlindx_, lindx_
      //idea: build tasklist per processor and then exchange
      //tuple is: src_snode,tgt_snode,op_type,lt_first_row = updated fc, facing_first_row = updated fr,
      using upd_tuple_t = std::tuple<Idx,Idx,Factorization::op_type,Idx,Idx>;
      std::map<Idx, std::list< upd_tuple_t> > Updates;
      std::vector<int> marker(this->np,0);

      std::vector< std::map< Idx, std::pair<Idx,std::set<Idx> > > > updCnt(this->Xsuper_.size() -1 );

      Int numLocSnode = this->XsuperDist_[this->iam+1]-this->XsuperDist_[this->iam];
      Int firstSnode = this->XsuperDist_[this->iam];
      for(Int locsupno = 1; locsupno<this->locXlindx_.size(); ++locsupno){
        Idx I = locsupno + firstSnode-1;
        Int first_col = this->Xsuper_[I-1];
        Int last_col = this->Xsuper_[I]-1;
        auto & fcell = cells_[coord2supidx(I-1,I-1)];
        Int iOwner = fcell.owner;

        //Create the factor task on the owner
        Updates[iOwner].push_back(std::make_tuple(I,I,Factorization::op_type::FACTOR,first_col,first_col));


        Ptr lfi = this->locXlindx_[locsupno-1];
        Ptr lli = this->locXlindx_[locsupno]-1;
        std::list<Int> ancestor_rows;

        Int K = -1; 
        Idx K_prevSnode = -1;
        for(Ptr K_sidx = lfi; K_sidx<=lli;K_sidx++){
          Idx K_row = this->locLindx_[K_sidx-1];
          K = this->SupMembership_[K_row-1];

          //Split at boundary or after diagonal block
          if(K!=K_prevSnode){
            if(K>I){
              ancestor_rows.push_back(K_row);
            }
          }
          K_prevSnode = K;
        }

        for(auto J_row : ancestor_rows){
          Int J = this->SupMembership_[J_row-1];

          auto & fodcell = cells_[coord2supidx(J-1,I-1)];
          Int iFODOwner = fodcell.owner;

//          if(iOwner!=iFODOwner){
            Updates[iOwner].push_back(std::make_tuple(I,I,Factorization::op_type::TRSM_SEND,J_row,J_row));
            Updates[iFODOwner].push_back(std::make_tuple(I,I,Factorization::op_type::TRSM_RECV,J_row,J_row));
//          }

          Updates[iFODOwner].push_back(std::make_tuple(I,I,Factorization::op_type::TRSM,J_row,J_row));

          //add the UPDATE_SEND from FODOwner tasks
          for(auto K_row : ancestor_rows){
            K = this->SupMembership_[K_row-1];

            if(K>=J){
              auto & tgtupdcell = cells_[coord2supidx(K-1,J-1)];
              Int iTgtOwner = tgtupdcell.owner;
              //TODO update this for non fan-out mapping
              Int iUpdOwner = tgtupdcell.owner;

//              if(iFODOwner!=iUpdOwner){
                //cell(J,I) to cell(K,J)
                Updates[iFODOwner].push_back(std::make_tuple(I,J,Factorization::op_type::UPDATE2D_SEND_OD,K_row,J_row));
                Updates[iUpdOwner].push_back(std::make_tuple(I,J,Factorization::op_type::UPDATE2D_RECV_OD,K_row,J_row));
//              }
            }

            if(K<=J){
              auto & tgtupdcell = cells_[coord2supidx(J-1,K-1)];
              Int iTgtOwner = tgtupdcell.owner;
              //TODO update this for non fan-out mapping
              Int iUpdOwner = tgtupdcell.owner;

              auto & facingcell = cells_[coord2supidx(J-1,I-1)];
              Int iFacingOwner = facingcell.owner;

              //sender point of view
              //cell(J,I) to cell(J,K)
//              if(iFacingOwner!=iUpdOwner){
                Updates[iFacingOwner].push_back(std::make_tuple(I,K,Factorization::op_type::UPDATE2D_SEND,K_row,J_row));
                Updates[iUpdOwner].push_back(std::make_tuple(I,K,Factorization::op_type::UPDATE2D_RECV,K_row,J_row));
//              }
            }
          }

          for(auto K_row : ancestor_rows){
            K = this->SupMembership_[K_row-1];
            if(K>=J){
              auto & tgtupdcell = cells_[coord2supidx(K-1,J-1)];
              Int iTgtOwner = tgtupdcell.owner;
              //TODO update this for non fan-out mapping
              Int iUpdOwner = tgtupdcell.owner;

              //update on cell(K,J)
              Updates[iUpdOwner].push_back(std::make_tuple(I,J,Factorization::op_type::UPDATE2D_COMP,J_row,K_row));

//              if(iUpdOwner!=iTgtOwner){
                //cell(K,J) to cell (K,J)
                Updates[iUpdOwner].push_back(std::make_tuple(I,J,Factorization::op_type::AGGREGATE2D_SEND,J_row,K_row));
                Updates[iTgtOwner].push_back(std::make_tuple(I,J,Factorization::op_type::AGGREGATE2D_RECV,J_row,K_row));
//              }
            }
          }
        }
      }

      MPI_Datatype type;
      MPI_Type_contiguous( sizeof(upd_tuple_t), MPI_BYTE, &type );
      MPI_Type_commit(&type);

      //then do an alltoallv
      //compute send sizes
      vector<int> ssizes(this->np,0);
      for(auto itp = Updates.begin();itp!=Updates.end();itp++){
        ssizes[itp->first] = itp->second.size();//*sizeof(std::pair<Idx,Idx>);
      }

      //compute send displacements
      vector<int> sdispls(this->np+1,0);
      sdispls[0] = 0;
      std::partial_sum(ssizes.begin(),ssizes.end(),&sdispls[1]);

      //Build the contiguous array of pairs
      vector<upd_tuple_t> sendbuf;
      sendbuf.reserve(sdispls.back());

      for(auto itp = Updates.begin();itp!=Updates.end();itp++){
        sendbuf.insert(sendbuf.end(),itp->second.begin(),itp->second.end());
      }

      //gather receive sizes
      vector<int> rsizes(this->np,0);
      MPI_Alltoall(&ssizes[0],sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,this->workcomm_);

      //compute receive displacements
      vector<int> rdispls(this->np+1,0);
      rdispls[0] = 0;
      std::partial_sum(rsizes.begin(),rsizes.end(),&rdispls[1]);

      //Now do the alltoallv
      vector<upd_tuple_t> recvbuf(rdispls.back());///sizeof(std::pair<Idx,Idx>));
      MPI_Alltoallv(&sendbuf[0],&ssizes[0],&sdispls[0],type,&recvbuf[0],&rsizes[0],&rdispls[0],type,this->workcomm_);

      MPI_Type_free(&type);

      using taskGraph2D_t = std::unordered_map< Int,  std::list< SparseTask2D * > > ;
      taskGraph2D_t task_graph;

      //do a top down traversal of my local task list
      for(auto it = recvbuf.begin();it!=recvbuf.end();it++){
        auto & cur_op = (*it);
        auto & src_snode = std::get<0>(cur_op);
        auto & tgt_snode = std::get<1>(cur_op);
        auto & type = std::get<2>(cur_op);
        auto & lt_first_row = std::get<3>(cur_op);
        auto & facing_first_row = std::get<4>(cur_op);

SparseTask2D * ptask = nullptr;
        switch(type){
          case Factorization::op_type::TRSM:
          case Factorization::op_type::FACTOR:
          case Factorization::op_type::UPDATE2D_COMP:
            {
              ptask = new SparseTask2D;
              size_t tuple_size = sizeof(upd_tuple_t);
              ptask->meta.resize(tuple_size);
              auto meta = reinterpret_cast<upd_tuple_t*>(ptask->meta.data());
              meta[0] = *it;
            }
            break;
          default:
            break;
        }


        auto I = src_snode;
        auto J = tgt_snode;
        switch(type){
          case Factorization::op_type::TRSM:
            {
              auto J = this->SupMembership_[facing_first_row-1];
              logfileptr->OFS()<<"TRSM"<<" from "<<I<<" to ("<<I<<") cell ("<<J<<","<<I<<")"<<std::endl;
              task_graph[coord2supidx(J-1,I-1)].push_back(ptask);
            }
            break;
          case Factorization::op_type::FACTOR:
            {
              logfileptr->OFS()<<"FACTOR"<<" cell ("<<J<<","<<I<<")"<<std::endl;
              task_graph[coord2supidx(I-1,I-1)].push_back(ptask);
            }
            break;
          case Factorization::op_type::UPDATE2D_COMP:
            {
              auto K = this->SupMembership_[facing_first_row-1];
              logfileptr->OFS()<<"UPDATE"<<" from "<<I<<" to "<<J<<" facing cell ("<<K<<","<<I<<") and cell ("<<J<<","<<I<<") to cell ("<<K<<","<<J<<")"<<std::endl;
              task_graph[coord2supidx(K-1,J-1)].push_back(ptask);
            }
            break;
          default:
            break;
        }
      }

      auto task_lookup = [&task_graph,this,&coord2supidx](Int cellI_, Int cellJ_, Int I_, Int J_, Int K_, Factorization::op_type type_){
        static std::unordered_map< Factorization::op_type , SparseTask2D * > last_ptr;


        //first check last_ptr
        auto last_it = last_ptr.find(type_);
        if ( last_it != last_ptr.end() ) {
          auto meta = reinterpret_cast<upd_tuple_t*>( last_it->second->meta.data() );
          auto & src_snode = std::get<0>(meta[0]);
          auto & tgt_snode = std::get<1>(meta[0]);
          auto & facing_row = std::get<4>(meta[0]);
          auto K = this->SupMembership_[facing_row-1];
          auto & type = std::get<2>(meta[0]);

          bassert(type==type_);

          if( (src_snode == I_) && (tgt_snode == J_) && ( K == K_) ){
            return last_it->second;
          }
        }

        auto idx = coord2supidx(cellI_-1,cellJ_-1);
        auto it = std::find_if(task_graph[idx].begin(),task_graph[idx].end(),
            [I_,J_,K_,type_,this]( SparseTask2D * ptask)->bool{
            auto meta = reinterpret_cast<upd_tuple_t*>(ptask->meta.data());
            auto & src_snode = std::get<0>(meta[0]);
            auto & tgt_snode = std::get<1>(meta[0]);
            auto & facing_row = std::get<4>(meta[0]);
            auto K = this->SupMembership_[facing_row-1];
            auto & type = std::get<2>(meta[0]);
            return (src_snode == I_) && (tgt_snode == J_) && (K == K_) && (type == type_) ;
            });

        if (it!=task_graph[idx].end()){
          //backup
          last_ptr[type_] = *it;
          return *it;
        }
        else{
          return (SparseTask2D*)nullptr;
        }
      };



      //process dependency tasks
      for(auto it = recvbuf.begin();it!=recvbuf.end();it++){
        auto & cur_op = (*it);
        auto & src_snode = std::get<0>(cur_op);
        auto & tgt_snode = std::get<1>(cur_op);
        auto & type = std::get<2>(cur_op);
        auto & lt_first_row = std::get<3>(cur_op);
        auto & facing_first_row = std::get<4>(cur_op);


        auto I = src_snode;
        auto J = tgt_snode;

        switch(type){
          case Factorization::op_type::TRSM_SEND:
            {
              auto J = this->SupMembership_[facing_first_row-1];

              auto taskptr = task_lookup(I,I,I,I,I,Factorization::op_type::FACTOR);
              bassert(taskptr!=nullptr); 
              logfileptr->OFS()<<"TRSM_SEND"<<" from "<<I<<" to "<<I<<" cell ("<<J<<","<<I<<")"<<std::endl;
              taskptr->out_dependencies.push_back( std::make_tuple(J,I) );
              logfileptr->OFS()<<"        out dep added to FACTOR"<<" from "<<I<<" to "<<I<<std::endl;
            }
            break;
          case Factorization::op_type::TRSM_RECV:
            {
              auto J = this->SupMembership_[facing_first_row-1];

              //find the TRSM and add cell I,I as incoming dependency
              auto taskptr = task_lookup(J,I,I,I,J,Factorization::op_type::TRSM);
              bassert(taskptr!=nullptr); 
              logfileptr->OFS()<<"TRSM_RECV"<<" from "<<I<<" to "<<I<<" cell ("<<J<<","<<I<<")"<<std::endl;
              taskptr->in_dependencies.push_back( std::make_tuple(I,I) );
              logfileptr->OFS()<<"        in dep added to TRSM"<<" from "<<I<<" to "<<I<<std::endl;
            }
            break;
          case Factorization::op_type::AGGREGATE2D_RECV:
            {
              auto K = this->SupMembership_[facing_first_row-1];

              logfileptr->OFS()<<"AGGREGATE2D_RECV"<<" from "<<I<<" to "<<J<<" cell ("<<K<<","<<J<<") to cell("<<K<<","<<J<<")"<<std::endl;
              //find the FACTOR or TRSM and add cell K,J as incoming dependency
              if(K==J){
                auto taskptr = task_lookup(K,J,J,J,J,Factorization::op_type::FACTOR);
                bassert(taskptr!=nullptr); 

                taskptr->in_dependencies.push_back( std::make_tuple(K,J) );
                logfileptr->OFS()<<"        in dep added to FACTOR"<<" from "<<J<<" to "<<J<<std::endl;
              }
              else{
                auto taskptr = task_lookup(K,J,J,J,K,Factorization::op_type::TRSM);
                bassert(taskptr!=nullptr); 
                taskptr->in_dependencies.push_back( std::make_tuple(K,J) );
                logfileptr->OFS()<<"        in dep added to TRSM"<<" from "<<J<<" to "<<J<<std::endl;
              }
            }
            break;
          case Factorization::op_type::AGGREGATE2D_SEND:
            {
              //TODO this might be useless
              auto K = this->SupMembership_[facing_first_row-1];

              logfileptr->OFS()<<"AGGREGATE2D_SEND"<<" from "<<I<<" to "<<J<<" cell ("<<K<<","<<J<<") to cell("<<K<<","<<J<<")"<<std::endl;
              //find the UPDATE2D_COMP and add cell K,J as outgoing dependency
              auto taskptr = task_lookup(K,J,I,J,K,Factorization::op_type::UPDATE2D_COMP);
              bassert(taskptr!=nullptr); 

              taskptr->out_dependencies.push_back( std::make_tuple(K,J) );
              logfileptr->OFS()<<"        out dep added to UPDATE_COMP"<<" from "<<I<<" to "<<J<<std::endl;
            }
            break;

          case Factorization::op_type::UPDATE2D_RECV:
            {
              //TODO this might be useless
              Idx K = this->SupMembership_[facing_first_row-1];
              std::swap(K,J);
              logfileptr->OFS()<<"UPDATE_RECV"<<" from "<<I<<" to "<<K<<" cell ("<<J<<","<<I<<") to cell ("<<J<<","<<K<<")"<<std::endl;
              //find the UPDATE2D_COMP and add cell J,I as incoming dependency
              auto taskptr = task_lookup(J,K,I,K,J,Factorization::op_type::UPDATE2D_COMP);
              bassert(taskptr!=nullptr);
              taskptr->in_dependencies.push_back( std::make_tuple(J,I) );
              logfileptr->OFS()<<"        in dep added to UPDATE_COMP"<<" from "<<I<<" to "<<K<<std::endl;
            }
            break;
          case Factorization::op_type::UPDATE2D_SEND:
            {
              Idx K = this->SupMembership_[facing_first_row-1];
              std::swap(K,J);
              logfileptr->OFS()<<"UPDATE_SEND"<<" from "<<I<<" to "<<K<<" cell ("<<J<<","<<I<<") to cell ("<<J<<","<<K<<")"<<std::endl;
              //find the TRSM and add cell J,K as outgoing dependency
              auto taskptr = task_lookup(J,I,I,I,J,Factorization::op_type::TRSM);
              bassert(taskptr!=nullptr); 
              taskptr->out_dependencies.push_back( std::make_tuple(J,K) );
              logfileptr->OFS()<<"        out dep added to TRSM"<<" from "<<I<<" to "<<I<<std::endl;
            }
            break;
          case Factorization::op_type::UPDATE2D_SEND_OD:
            {
              auto K = this->SupMembership_[lt_first_row-1];
              logfileptr->OFS()<<"UPDATE_SEND OD"<<" from "<<I<<" to "<<J<<" cell ("<<J<<","<<I<<") to cell ("<<K<<","<<J<<")"<<std::endl;

              //find the TRSM and add cell K,J as outgoing dependency
              auto taskptr = task_lookup(J,I,I,I,J,Factorization::op_type::TRSM);
              bassert(taskptr!=nullptr); 
              taskptr->out_dependencies.push_back( std::make_tuple(K,J) );
              logfileptr->OFS()<<"        out dep DOWN added to TRSM"<<" from "<<I<<" to "<<I<<std::endl;
            }
            break;
          case Factorization::op_type::UPDATE2D_RECV_OD:
            {
              auto K = this->SupMembership_[lt_first_row-1];
              logfileptr->OFS()<<"UPDATE_RECV OD"<<" from "<<I<<" to "<<J<<" cell ("<<J<<","<<I<<") to cell ("<<K<<","<<J<<")"<<std::endl;

              //find the UPDATE2D_COMP and add cell J,I as incoming dependency
              auto taskptr = task_lookup(K,J,I,J,K,Factorization::op_type::UPDATE2D_COMP);
              bassert(taskptr!=nullptr); 
              taskptr->in_dependencies.push_back( std::make_tuple(J,I) );
              logfileptr->OFS()<<"        in dep DOWN added to UPDATE2D_COMP"<<" from "<<I<<" to "<<J<<std::endl;
            }
            break;
          default:
            break;
        }

      }


//Now we have our local part of the task graph
for(auto it = task_graph.begin(); it != task_graph.end(); it++){
  auto idx = it->first;
  auto tasklist = it->second;

  for( auto && ptask: tasklist){
     auto meta = reinterpret_cast<upd_tuple_t*>(ptask->meta.data());
     auto & src_snode = std::get<0>(meta[0]);
     auto & tgt_snode = std::get<1>(meta[0]);
     auto & lt_first_row = std::get<3>(meta[0]);
     auto & facing_row = std::get<4>(meta[0]);
     auto & type = std::get<2>(meta[0]);

     auto I = src_snode;
     auto J = tgt_snode;
     auto K = this->SupMembership_[facing_row-1];

        switch(type){
          case Factorization::op_type::TRSM:
            {
              logfileptr->OFS()<<"TRSM"<<" from "<<I<<" to ("<<I<<") cell ("<<K<<","<<I<<")"<<std::endl;
            }
            break;
          case Factorization::op_type::FACTOR:
            {
              logfileptr->OFS()<<"FACTOR"<<" cell ("<<I<<","<<I<<")"<<std::endl;
            }
            break;
          case Factorization::op_type::UPDATE2D_COMP:
            {
              logfileptr->OFS()<<"UPDATE"<<" from "<<I<<" to "<<J<<" facing cell ("<<K<<","<<I<<") and cell ("<<J<<","<<I<<") to cell ("<<K<<","<<J<<")"<<std::endl;
            }
            break;
          default:
            delete ptask;
            break;
        }

      logfileptr->OFS()<<"      input dependencies: ";
      for(auto &&tpl : ptask->in_dependencies){
        logfileptr->OFS()<<"("<<std::get<0>(tpl)<<","<<std::get<1>(tpl)<<") ";
      }
      logfileptr->OFS()<<std::endl;
      logfileptr->OFS()<<"      output dependencies: ";
      for(auto &&tpl : ptask->out_dependencies){
        logfileptr->OFS()<<"("<<std::get<0>(tpl)<<","<<std::get<1>(tpl)<<") ";
      }
      logfileptr->OFS()<<std::endl;


  }
}




////      //do a top down traversal of my local task list
////      for(auto it = recvbuf.begin();it!=recvbuf.end();it++){
////        auto & cur_op = (*it);
////        auto & src_snode = std::get<0>(cur_op);
////        auto & tgt_snode = std::get<1>(cur_op);
////        auto & type = std::get<2>(cur_op);
////        auto & lt_first_row = std::get<3>(cur_op);
////        auto & facing_first_row = std::get<4>(cur_op);
////
////        std::shared_ptr<GenericTask> pTask(new SparseTask2D);
////        SparseTask2D & Task = *(SparseTask2D*)pTask.get();
////
////
////        auto I = src_snode;
////        auto J = tgt_snode;
////
////        size_t tuple_size = sizeof(upd_tuple_t);
////        Task.meta.resize(tuple_size);
////        auto meta = reinterpret_cast<upd_tuple_t*>(Task.meta.data());
////        meta[0] = *it;
////
////
////        switch(type){
////          case Factorization::op_type::TRSM:
////            {
////              auto J = this->SupMembership_[facing_first_row-1];
////              last_trsm = pTask;
////              logfileptr->OFS()<<"TRSM"<<" from "<<I<<" to ("<<I<<") cell ("<<J<<","<<I<<")"<<std::endl;
////
////              task_graph[coord2supidx(J-1,I-1)].push_back(&Task);
////
////              auto & fcell = cells_[coord2supidx(I-1,I-1)];
////              Int iOwner = fcell.owner;
////              if ( iOwner == this->iam ) {
////                //enqueue Factor I,I as an incoming dependency
////                Task.in_dependencies.push_back( std::make_tuple(I,I) );
////                
////                SparseTask2D & fTask = *(SparseTask2D*)last_factor.get();
////                bassert(last_factor!=nullptr);
////                auto meta = reinterpret_cast<upd_tuple_t*>(fTask.meta.data());
////                bassert ( std::get<0>(meta[0]) == I && std::get<1>(meta[0]) == I) ;
////                fTask.out_dependencies.push_back( std::make_tuple(J,I) );
////                logfileptr->OFS()<<"        out dep added to FACTOR"<<" from "<<I<<" to "<<I<<std::endl;
////              }
////              else{
////                //enqueue TRSM_RECV J,I as an incoming dependency
////                Task.in_dependencies.push_back( std::make_tuple(J,I) );
////
////                bassert(last_trsm_recv!=nullptr);
////                SparseTask2D & trTask = *(SparseTask2D*)last_trsm_recv.get();
////                auto meta = reinterpret_cast<upd_tuple_t*>(trTask.meta.data());
////                bassert ( std::get<0>(meta[0]) == I && std::get<1>(meta[0]) == I 
////                    && this->SupMembership_[std::get<4>(meta[0])-1] == J);
////                trTask.out_dependencies.push_back( std::make_tuple(J,I) );
////                logfileptr->OFS()<<"        out dep added to TSRM_RECV"<<" from "<<I<<" to "<<I<<std::endl;
////              } 
////            }
////            break;
////          case Factorization::op_type::TRSM_SEND:
////            {
////              last_trsm_send = pTask;
////              auto J = this->SupMembership_[facing_first_row-1];
////
////              task_graph[coord2supidx(J-1,I-1)].push_back(&Task);
////
////              logfileptr->OFS()<<"TRSM_SEND"<<" from "<<I<<" to "<<I<<" cell ("<<J<<","<<I<<")"<<std::endl;
////              //enqueue Factor I,I as an incoming dependency
////              Task.in_dependencies.push_back( std::make_tuple(I,I) );
////              Task.out_dependencies.push_back( std::make_tuple(J,I) );
////
////              //add this as an outgoing dependency to the last factor
////              SparseTask2D & fTask = *(SparseTask2D*)last_factor.get();
////              auto meta = reinterpret_cast<upd_tuple_t*>(fTask.meta.data());
////              bassert(last_factor!=nullptr);
////              bassert( std::get<0>(meta[0]) == I && std::get<1>(meta[0]) == I);
////              fTask.out_dependencies.push_back( std::make_tuple(J,I) );
////              logfileptr->OFS()<<"        out dep added to FACTOR"<<" from "<<I<<" to "<<I<<std::endl;
////            }
////            break;
////          case Factorization::op_type::TRSM_RECV:
////            {
////              last_trsm_recv = pTask;
////              auto J = this->SupMembership_[facing_first_row-1];
////              task_graph[coord2supidx(J-1,I-1)].push_back(&Task);
////
////              logfileptr->OFS()<<"TRSM_RECV"<<" from "<<I<<" to "<<I<<" cell ("<<J<<","<<I<<")"<<std::endl;
////
////              //enqueue TRSM_SEND J,I as an incoming dependency
////              Task.in_dependencies.push_back( std::make_tuple(J,I) );
////              Task.out_dependencies.push_back( std::make_tuple(J,I) );
////            }
////            break;
////
////          case Factorization::op_type::AGGREGATE2D_SEND:
////            {
////              last_aggregate_send = pTask;
////              auto K = this->SupMembership_[facing_first_row-1];
////              task_graph[coord2supidx(K-1,J-1)].push_back(&Task);
////              logfileptr->OFS()<<"AGGREGATE2D_SEND"<<" from "<<I<<" to "<<J<<" cell ("<<K<<","<<J<<") to cell("<<K<<","<<K<<")"<<std::endl;
////
////              Task.in_dependencies.push_back( std::make_tuple(I,J) );
////
////              SparseTask2D & uTask = *(SparseTask2D*)last_update.get();
////              bassert(last_update!=nullptr);
////
////              auto meta = reinterpret_cast<upd_tuple_t*>(uTask.meta.data());
////              bassert( std::get<0>(meta[0]) == I && std::get<1>(meta[0]) == J);
////              uTask.out_dependencies.push_back( std::make_tuple(K,J) );
////              logfileptr->OFS()<<"        out dep added to UPDATE_COMP ("<<K<<","<<J<<") from "<<I<<" to "<<J<<std::endl;
////            }
////            break;
////          case Factorization::op_type::FACTOR:
////            {
////              last_factor = pTask;
////              task_graph[coord2supidx(I-1,I-1)].push_back(&Task);
////              logfileptr->OFS()<<"FACTOR"<<" cell ("<<J<<","<<I<<")"<<std::endl;
////
////
////
////            }
////            break;
////          case Factorization::op_type::UPDATE2D_COMP:
////            {
////              last_update = pTask;
////              auto K = this->SupMembership_[facing_first_row-1];
////              task_graph[coord2supidx(K-1,J-1)].push_back(&Task);
////              logfileptr->OFS()<<"UPDATE"<<" from "<<I<<" to "<<J<<" facing cell ("<<K<<","<<I<<") and cell ("<<J<<","<<I<<") to cell ("<<K<<","<<J<<")"<<std::endl;
////
////              Task.in_dependencies.push_back( std::make_tuple(J,I) );
////              Task.in_dependencies.push_back( std::make_tuple(K,I) );
////              Task.out_dependencies.push_back( std::make_tuple(K,J) );
////
////              auto & fcell = cells_[coord2supidx(J-1,I-1)];
////              Int iOwner = fcell.owner;
////              if ( iOwner == this->iam ) {
////                //add this as an outgoing dependency to the last TRSM OD
////                bassert(last_trsm!=nullptr);
////                SparseTask2D & tTask = *(SparseTask2D*)last_trsm.get();
////                auto meta = reinterpret_cast<upd_tuple_t*>(tTask.meta.data());
////                bassert( std::get<0>(meta[0]) == I && std::get<1>(meta[0]) == I 
////                          && this->SupMembership_[std::get<4>(meta[0])-1] == J);
////
////                tTask.out_dependencies.push_back( std::make_tuple(K,J) );
////                logfileptr->OFS()<<"        out dep added to TRSM ("<<J<<","<<I<<")"<<std::endl;
////
////              }
////              else{
////                //add this as an outgoing dependency to the last UPDATE2D_RECV_OD
////                bassert(last_update_od_recv!=nullptr);
////                SparseTask2D & tTask = *(SparseTask2D*)last_update_od_recv.get();
////                auto meta = reinterpret_cast<upd_tuple_t*>(tTask.meta.data());
////                bassert( std::get<0>(meta[0]) == I && std::get<1>(meta[0]) == J 
////                          && this->SupMembership_[std::get<4>(meta[0])-1] == J);
////                tTask.out_dependencies.push_back( std::make_tuple(K,J) );
////                logfileptr->OFS()<<"        out dep added to UPDATE2D_RECV OD ("<<J<<","<<I<<")"<<std::endl;
////              }
////
////              auto & odcell = cells_[coord2supidx(K-1,I-1)];
////              Int iFODOwner = odcell.owner;
////              if ( iOwner == this->iam ) {
////                /////add this as an outgoing dependency to the last top TRSM
////                ///bassert(last_top_trsm!=nullptr);
////                ///SparseTask2D & tTask = *(SparseTask2D*)last_top_trsm.get();
////                ///auto meta = reinterpret_cast<upd_tuple_t*>(tTask.meta.data());
////                ///bassert( std::get<0>(meta[0]) == I && std::get<1>(meta[0]) == I 
////                ///          && this->SupMembership_[std::get<4>(meta[0])-1] == J);
////                ///tTask.out_dependencies.push_back( std::make_tuple(J,K) );
////                ///logfileptr->OFS()<<"        out dep added to TRSM ("<<J<<","<<I<<")"<<std::endl;
////              }
////              else{
////                //add this as an outgoing dependency to the last UPDATE2D_RECV_OD
////                //bassert(last_update_recv!=nullptr);
////                //SparseTask2D & tTask = *(SparseTask2D*)last_update_recv.get();
////                //auto meta = reinterpret_cast<upd_tuple_t*>(tTask.meta.data());
////                //bassert( std::get<0>(meta[0]) == I && std::get<1>(meta[0]) == J 
////                //          && this->SupMembership_[std::get<4>(meta[0])-1] == K);
////                //tTask.out_dependencies.push_back( std::make_tuple(J,K) );
////                //logfileptr->OFS()<<"        out dep added to UPDATE2D_RECV ("<<K<<","<<I<<")"<<std::endl;
////              }
////
////            }
////            break;
////          case Factorization::op_type::UPDATE2D_SEND_OD:
////            {
////              last_update_od_send = pTask;
////              auto K = this->SupMembership_[facing_first_row-1];
////              logfileptr->OFS()<<"UPDATE_SEND OD"<<" from "<<I<<" to "<<J<<" cell ("<<J<<","<<I<<") to cell ("<<K<<","<<J<<")"<<std::endl;
////              task_graph[coord2supidx(J-1,I-1)].push_back(&Task);
////
////              Task.in_dependencies.push_back( std::make_tuple(J,I) );
////              Task.out_dependencies.push_back( std::make_tuple(K,J) );
////
////              //add this as an outgoing dependency to the last TRSM if it corresponds
////              SparseTask2D & tTask = *(SparseTask2D*)last_trsm.get();
////              bassert(last_trsm!=nullptr);
////              auto meta = reinterpret_cast<upd_tuple_t*>(tTask.meta.data());
////              bassert( std::get<0>(meta[0]) == I && std::get<1>(meta[0]) == I 
////                  && this->SupMembership_[std::get<4>(meta[0])-1] == J);
////              tTask.out_dependencies.push_back( std::make_tuple(K,J) );
////              logfileptr->OFS()<<"        out dep DOWN added to TRSM ("<<J<<","<<I<<")"<<std::endl;
////            }
////            break;
////          case Factorization::op_type::UPDATE2D_RECV_OD:
////            {
////              last_update_od_recv = pTask;
////              auto K = this->SupMembership_[facing_first_row-1];
////              logfileptr->OFS()<<"UPDATE_RECV OD"<<" from "<<I<<" to "<<J<<" cell ("<<J<<","<<I<<") to cell ("<<K<<","<<J<<")"<<std::endl;
////              task_graph[coord2supidx(K-1,J-1)].push_back(&Task);
////
////              Task.in_dependencies.push_back( std::make_tuple(J,I) );
////              Task.out_dependencies.push_back( std::make_tuple(K,J) );
////            }
////            break;
////          case Factorization::op_type::UPDATE2D_SEND:
////            {
////              last_update_send = pTask;
////              auto K = this->SupMembership_[facing_first_row-1];
////                logfileptr->OFS()<<"UPDATE_SEND"<<" from "<<I<<" to "<<J<<" cell ("<<K<<","<<I<<") to cell ("<<K<<","<<J<<")"<<std::endl;
////              task_graph[coord2supidx(K-1,I-1)].push_back(&Task);
////
////              Task.in_dependencies.push_back( std::make_tuple(K,I) );
////              Task.out_dependencies.push_back( std::make_tuple(K,J) );
////
////              //add this as an outgoing dependency to the last TRSM if it corresponds
////              SparseTask2D & tTask = *(SparseTask2D*)last_trsm.get();
////              bassert(last_trsm!=nullptr);
////              auto meta = reinterpret_cast<upd_tuple_t*>(tTask.meta.data());
////              bassert( std::get<0>(meta[0]) == I && std::get<1>(meta[0]) == I 
////                          && this->SupMembership_[std::get<4>(meta[0])-1] == K);
////
////                tTask.out_dependencies.push_back( std::make_tuple(J,K) );
////                logfileptr->OFS()<<"        out dep added to TRSM ("<<K<<","<<I<<")"<<std::endl;
////            }
////            break;
////          case Factorization::op_type::UPDATE2D_RECV:
////            {
////              last_update_recv = pTask;
////              auto K = this->SupMembership_[facing_first_row-1];
////                logfileptr->OFS()<<"UPDATE_RECV"<<" from "<<I<<" to "<<J<<" cell ("<<K<<","<<I<<") to cell ("<<K<<","<<J<<")"<<std::endl;
////              task_graph[coord2supidx(K-1,J-1)].push_back(&Task);
////
////              Task.in_dependencies.push_back( std::make_tuple(K,I) );
////              Task.out_dependencies.push_back( std::make_tuple(K,J) );
////
////            }
////            break;
////        }
////
////        //push the task somewhere, in a per supernode bucket list for instance ?
////      }







    }



    //redistribute the matrix

    } 

  template <typename colptr_t, typename rowind_t, typename T>
   void symPACKMatrix2D<colptr_t,rowind_t,T>::DistributeMatrix(DistSparseMatrix<T> & pMat ){
   } 
}

//#include <sympack/impl/symPACKMatrix2D_impl.hpp>
#endif //_SYMPACK_MATRIX2D_DECL_HPP_

