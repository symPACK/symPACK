#ifndef _SUPERNODAL_MATRIX_DECL_HPP_
#define _SUPERNODAL_MATRIX_DECL_HPP_

#include "Environment.hpp"
#include "SuperNode.hpp"

#include "NumVec.hpp"
#include "DistSparseMatrix.hpp"
#include "ETree.hpp"
#include "FBMatrix.hpp"



#include <list>
#include <vector>
#include "Mapping.hpp"



namespace LIBCHOLESKY{

struct SnodeUpdate;

struct LocalUpdate{
  Int src_snode_id;
  Int src_nzblk_idx;
  Int src_first_row;
  LocalUpdate(Int snode_id,Int nzblk_idx,Int first_row):src_snode_id(snode_id),src_nzblk_idx(nzblk_idx),src_first_row(first_row){};
};

struct DelayedComm{
  Int tgt_snode_id;
  Int src_snode_id;

  Int src_nzblk_idx;
  Int src_first_row;
  Int src_last_row;

  DelayedComm(Int a_src_snode_id, Int a_tgt_snode_id, Int a_src_nzblk_idx, Int a_src_first_row):src_snode_id(a_src_snode_id),tgt_snode_id(a_tgt_snode_id),src_first_row(a_src_first_row),src_nzblk_idx(a_src_nzblk_idx){};
};

struct OutgoingComm{
  std::vector<char> * pSrcBlocks;
  Int head;
  MPI_Request Request;
  OutgoingComm(Int aSize, MPI_Request aRequest):Request(aRequest){
   pSrcBlocks = new std::vector<char>(aSize);
    head = 0;
  };
  ~OutgoingComm(){  
    delete pSrcBlocks; 
    if(Request !=MPI_REQUEST_NULL){
      MPI_Request_free(&Request);
    }
  };
  inline char * back(){ return &pSrcBlocks->at(head);}
  inline char * front(){ return &pSrcBlocks->front();}
  inline Int size(){ return pSrcBlocks->size();}

};

template <typename T> inline void Serialize( OutgoingComm& os,  const T * val, const Int count){
  T* dest = reinterpret_cast<T*>(&os.pSrcBlocks->at(os.head));
  std::copy(val,val+count,dest);
  os.head += count*sizeof(T);
}

template <typename T> inline OutgoingComm& operator<<( OutgoingComm& os,  const T & val){
  T* dest = reinterpret_cast<T*>(&os.pSrcBlocks->at(os.head));
  std::copy(&val,&val+1,dest);
  os.head += sizeof(T);

  return os;
}





typedef std::list<DelayedComm> CommList;
typedef std::list<OutgoingComm *> Isends;

template <typename T> class SupernodalMatrix{
  protected:




  bool globalAllocated = false;

  IntNumVec Xsuper_;
  IntNumVec SupMembership_;
  MAPCLASS Mapping_;
 
  SparseMatrixStructure Local_;
  SparseMatrixStructure Global_;

  IntNumVec UpdateCount_;
  IntNumVec UpdateWidth_;

  ETree ETree_;
  ETree SupETree_;
  Int iSize_;
  Int maxIsend_;
  std::vector<SuperNode<T> * > LocalSupernodes_;
  std::vector<SuperNode<T> *> Contributions_;
  IntNumVec xlindx_;
  IntNumVec lindx_;




  void AddOutgoingComm(Isends & outgoingSend, Int src_snode_id, Int src_snode_size ,Int src_first_row, NZBlockDesc & pivot_desc, Int nzblk_cnt, T * nzval_ptr, Int nz_cnt);


  void GetUpdatingSupernodeCount( IntNumVec & sc,IntNumVec & mw);


  inline void FindUpdates(SuperNode<T> & src_snode, std::list<SnodeUpdate> & updates  );
  inline bool FindNextUpdate(SuperNode<T> & src_snode, Int & src_nzblk_idx, Int & src_first_row,  Int & src_last_row, Int & tgt_snode_id);


#ifdef SINGLE_BLAS
  inline void UpdateSuperNode(SuperNode<T> & src_snode, SuperNode<T> & tgt_snode,Int & pivot_idx, NumMat<T> & tmpBuf, Int  pivot_fr = I_ZERO);
#else
  inline void UpdateSuperNode(SuperNode<T> & src_snode, SuperNode<T> & tgt_snode,Int & pivot_idx, Int  pivot_fr = I_ZERO);
#endif

  inline void AggregateSuperNode(SuperNode<T> & src_snode, SuperNode<T> & tgt_snode, Int &pivot_idx, Int  pivot_fr = I_ZERO);

    void SendDelayedMessages(Int cur_snode_id, CommList & MsgToSend, Isends & OutgoingSend);




  Int FBUpdate(Int I);
  void FBGetUpdateCount(IntNumVec & sc, IntNumVec & lu);

  SuperNode<T> * FBRecvFactor(Int src_snode_id,Int tgt_snode_id, std::vector<char> & src_blocks);

  public:

	MPI_Comm     pComm;        
  
  IntNumVec perm_;


  void forward_update(SuperNode<T> * src_contrib,SuperNode<T> * tgt_contrib);
  void back_update(SuperNode<T> * src_contrib,SuperNode<T> * tgt_contrib);


	/// @brief MPI communicator
	//MPI_Comm     comm = MPI_COMM_NULL;        


  SupernodalMatrix();
  SupernodalMatrix(const DistSparseMatrix<T> & pMat, Int maxSnode,MAPCLASS & pMapping, Int maxIsend, MPI_Comm & pComm );

  //Accessors
  Int Size(){return iSize_;}
  IntNumVec & GetSupernodes(){ return Xsuper_;}
  std::vector<SuperNode<T> *  > & GetLocalSupernodes(){ return LocalSupernodes_; } 
  SuperNode<T> & GetLocalSupernode(Int i){ return *LocalSupernodes_[i]; } 
  ETree & GetETree(){return ETree_;}

  SparseMatrixStructure GetGlobalStructure();
  SparseMatrixStructure GetLocalStructure() const;

  void Factorize(MPI_Comm & pComm);
  void FanOut( MPI_Comm & pComm );
  void FanBoth( MPI_Comm & pComm );



#ifdef _CHECK_RESULT_SEQ_
  void Solve(NumMat<T> * RHS, MPI_Comm & pComm,NumMat<T> & forwardSol, NumMat<T> * Xptr=NULL);
#else
  void Solve(NumMat<T> * RHS, MPI_Comm & pComm,/*NumMat<T> & forwardSol,*/ NumMat<T> * Xptr=NULL);
#endif

  void GetFullFactors( NumMat<T> & fullMatrix, MPI_Comm &pComm);

  void GetSolution(NumMat<T> & B, MPI_Comm pComm);
};



} // namespace LIBCHOLESKY

#include "SupernodalMatrix_impl.hpp"

#endif 
