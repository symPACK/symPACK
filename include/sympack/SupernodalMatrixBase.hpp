#ifndef _SUPERNODAL_MATRIX_BASE_DECL_HPP_
#define _SUPERNODAL_MATRIX_BASE_DECL_HPP_

#include <list>
#include <vector>

namespace SYMPACK{

  struct SnodeUpdateFB;

  class SupernodalMatrixBase{
    protected:
//      std::vector<std::list<SnodeUpdateFB> > taskLists_;
//      std::list<SnodeUpdateFB*> readyTasks_;
      
      //MAPCLASS describing the Mapping of the computations
//      Mapping * Mapping_;

      
//      virtual inline bool FindNextSupernodeUpdate(Int snodeId, SnodeUpdate & nextUpdate, bool isLocal=true) = 0;

    public:
      SupernodalMatrixBase(){}


//      inline void FactorReceived(Int src_snode, Int tgt_snode){
//        //go through the tasks in this supernode
//
//            //go through supernode column structure to see whether
//            //there are other update tasks using this factor
////this doesn't work as the column have more non zeroes than the factor
//    SnodeUpdate curUpdate;
//    while(src_snode.FindNextUpdate(curUpdate,Xsuper_,SupMembership_)){ 
//      Int iTarget = this->Mapping_->Map(curUpdate.tgt_snode_id-1,src_snode.Id()-1);
//    while(FindNextSupernodeUpdate(tgt,curUpdate)){  
//      Int iTarget = Mapping_->Map(curUpdate.tgt_snode_id-1,src_snode.Id()-1);
//      if(iam==iTarget){
//        for(auto it = taskLists_[curUpdate.tgt_snode-1].begin(); it!=taskLists_[curUpdate.tgt_snode-1].end(); it++){
//          if(it->src==src){
//          }
//        }
//      }
//    }
//
//
//
//        assert(taskLists_.size()>=tgt);
//        for(auto it = taskLists_[tgt-1].begin(); it!=taskLists_[tgt-1].end(); it++){
//          if(it->src==src){
//            //check if it's local or remote
//            assert(it->remote_deps>0);
//            it->remote_deps--;
//
//
//            //if task is ready, move it to the ready list
//            if(it->remote_deps==0 && it->local_deps==0){
//              readyTasks_.push_back(&*it);
//            }
//          }
//        }
//      }; 
  };


}
#endif
