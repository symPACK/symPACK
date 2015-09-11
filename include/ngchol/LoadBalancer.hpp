#ifndef _LOAD_BALANCER_HPP_
#define _LOAD_BALANCER_HPP_

#include "ngchol/Environment.hpp"
#include "ngchol/ETree.hpp"



namespace LIBCHOLESKY{


class LoadBalancer{
  protected:
    std::vector<Int> procMap_;
    Int np_;
    Int n_;

  public:
      void Dump(){
        logfileptr->OFS()<<"Resulting proc mapping: "<<endl;
        for(auto it = procMap_.begin(); it!= procMap_.end();++it){
            logfileptr->OFS()<<*it<<" ";
        }
        logfileptr->OFS()<<endl;
      }
  
      LoadBalancer(Int np, Int n){
        n_=n;
        np_=np;
        procMap_.resize(n);
      };

      virtual inline std::vector<Int> & GetMap() =0;
};



class NNZBalancer: public LoadBalancer{
  protected:
    IntNumVec & Xsuper_;
    IntNumVec & cc_;
  public:
      NNZBalancer(Int np, IntNumVec & Xsuper, IntNumVec & pCc):Xsuper_(Xsuper),cc_(pCc),LoadBalancer(np,Xsuper.m()-1){
      };

 
      virtual inline std::vector<Int> & GetMap(){
          std::vector<double> load(np_,0.0);

          for(Int i = 1; i< Xsuper_.m();  ++i){
            //find least loaded processor
            vector<double>::iterator it = std::min_element(load.begin(),load.end());
            Int proc = (Int)(it - load.begin());
            Int width = Xsuper_[i] - Xsuper_[i-1];
            Int height = cc_[i-1];
            *it += (double)(width*height);
            procMap_[i-1] = proc;
          } 


      }
};



class TreeLoadBalancer: public LoadBalancer{
  public:
      class ProcGroup{
        protected:
        vector<Int> ranks_;
        double total_load_;
//        Int worker_;
        public:
        vector<Int> & Ranks(){return ranks_;}
//        Int & Worker(){return worker_;}
        double & Load(){return total_load_;}


  friend std::ostream& operator<<( std::ostream& os,  ProcGroup& group){
//    os<<"Worker is P"<<group.Worker()<<std::endl;
    os<<"Members are: "<<std::endl;
    for(Int ranks =0; ranks<group.Ranks().size();++ranks){
      os<<"P"<<group.Ranks()[ranks]<<" ";
    }
    os<<std::endl;
    return os;
  }



      };

  protected:
    vector< ProcGroup > procGroups_;
    ETree supETree_;

    vector< Int > groupIdx_;
    vector< Int > groupWorker_;
    vector< ProcGroup > levelGroups_;

    virtual double factor_cost(Int m, Int n) = 0;
    virtual double update_cost(Int m, Int n, Int k)=0;


  public:
      TreeLoadBalancer(Int np, ETree& supETree):supETree_(supETree),LoadBalancer(np,supETree.Size()){
          groupIdx_.resize(n_+1,0);
          groupWorker_.resize(n_+1,0);
      };

      vector< ProcGroup > & ProcGroups(){return procGroups_;}

      vector< ProcGroup > & LevelGroups(){return levelGroups_;}
      vector< Int > & GroupIdx(){return groupIdx_;}
      vector< Int > & GroupWorker(){return groupWorker_;}


};

class SubtreeToSubcube: public TreeLoadBalancer{

  protected:
    bool fan_in_;
    IntNumVec & Xsuper_;
    IntNumVec & SupMembership_;
    PtrVec & Xlindx_;
    IdxVec & Lindx_;
    IntNumVec & cc_;


    double factor_cost(Int m, Int n){
      return CHOLESKY_COST(m,n);
    }

    double update_cost(Int m, Int n, Int k){
      return 2.0*m*n*k;
    }



  public:
      SubtreeToSubcube(Int np, ETree & supETree,IntNumVec & Xsuper, IntNumVec & SupMembership,PtrVec & Xlindx, IdxVec & Lindx, IntNumVec & pCc, bool fan_in = true):Xsuper_(Xsuper),SupMembership_(SupMembership),cc_(pCc),Xlindx_(Xlindx),Lindx_(Lindx),TreeLoadBalancer(np,supETree){
          fan_in_=fan_in;
      };


  virtual inline std::vector<Int> & GetMap(){
      if(procGroups_.size()==0){
          //compute number of children and load
          Int numLevels = 1;
          vector<double> SubTreeLoad(supETree_.Size()+1,0.0);
          vector<double> NodeLoad(supETree_.Size()+1,0.0);
          vector<Int> children(supETree_.Size()+1,0);
          for(Int I=1;I<=supETree_.Size();I++){
            Int parent = supETree_.Parent(I-1);
            ++children[parent];
            Idx fc = Xsuper_[I-1];
            Int width = Xsuper_[I] - Xsuper_[I-1];
            Int height = cc_[fc-1];
            //cost of factoring curent panel
            double local_load = factor_cost(height,width);
            SubTreeLoad[I]+=local_load;
            NodeLoad[I]+=local_load;
            //cost of updating ancestors
            {
              Ptr fi = Xlindx_[I-1];
              Ptr li = Xlindx_[I]-1;

              //parse rows
              Int tgt_snode_id = I;
              Idx tgt_fr = fc;
              Idx tgt_lr = fc;
              Ptr tgt_fr_ptr = 0;
              Ptr tgt_lr_ptr = 0;
              for(Ptr rowptr = fi; rowptr<=li;rowptr++){
                Idx row = Lindx_[rowptr-1]; 
                if(SupMembership_[row-1]==tgt_snode_id){ continue;}

                //we have a new tgt_snode_id
                tgt_lr_ptr = rowptr-1;
                tgt_lr = Lindx_[rowptr-1 -1];
                if(tgt_snode_id!=I){
                  Int m = li-tgt_fr_ptr+1;
                  Int n = width;
                  Int k = tgt_lr_ptr - tgt_fr_ptr+1; 
                  double cost = update_cost(m,n,k);
                  if(fan_in_){
                    SubTreeLoad[I]+=cost;
                    NodeLoad[I]+=cost;
                  }
                  else{
                    SubTreeLoad[tgt_snode_id]+=cost;
                    NodeLoad[tgt_snode_id]+=cost;
                  }
                }
                tgt_fr = row;
                tgt_fr_ptr = rowptr;
                tgt_snode_id = SupMembership_[row-1];
              }
              //do the last update supernode
              //we have a new tgt_snode_id
              tgt_lr_ptr = li;
              tgt_lr = Lindx_[li -1];
              if(tgt_snode_id!=I){
                Int m = li-tgt_fr_ptr+1;
                Int n = width;
                Int k = tgt_lr_ptr - tgt_fr_ptr+1; 
                double cost = update_cost(m,n,k);
                  if(fan_in_){
                    SubTreeLoad[I]+=cost;
                    NodeLoad[I]+=cost;
                  }
                  else{
                    SubTreeLoad[tgt_snode_id]+=cost;
                    NodeLoad[tgt_snode_id]+=cost;
                  }
                }
            }
            SubTreeLoad[parent]+=SubTreeLoad[I];
          }

          vector<Int> levels(supETree_.Size()+1,0);
          for(Int I=n_; I>= 1;I--){ 
            Int parent = supETree_.Parent(I-1);
            if(parent==0){levels[I]=0;}
            else{ levels[I] = levels[parent]+1; numLevels = max(numLevels,levels[I]);}
          }
          numLevels++;
            logfileptr->OFS()<<"levels : "; 
            for(Int i = 0; i<levels.size();++i){logfileptr->OFS()<<levels.at(i)<<" ";}
            logfileptr->OFS()<<endl;


          logfileptr->OFS()<<"SubTreeLoad is "<<SubTreeLoad<<endl;


          //procmaps[0]/pstart[0] represents the complete list
          procGroups_.resize(n_+1);
          vector<Int> pstart(n_+1,0);
          procGroups_[0].Ranks().reserve(np);
          for(Int p = 0;p<np;++p){ procGroups_[0].Ranks().push_back(p);}


          levelGroups_.push_back(ProcGroup());//reserve(numLevels);
          levelGroups_[0].Ranks().reserve(np);
          for(Int p = 0;p<np;++p){levelGroups_[0].Ranks().push_back(p);}


          for(Int I=n_; I>= 1;I--){ 
            Int parent = supETree_.Parent(I-1);

            //split parent's proc list
            double parent_load = 0.0;

            if(parent!=0){
              Int fc = Xsuper_[parent-1];
              Int width = Xsuper_[parent] - Xsuper_[parent-1];
              Int height = cc_[fc-1];
              parent_load = NodeLoad[parent];// factor_cost(height,width);
              procGroups_[parent].Load() = parent_load;
            }


            double proportion = min(1.0,SubTreeLoad[I]/(SubTreeLoad[parent]-parent_load));
            Int npParent = procGroups_[parent].Ranks().size();
            Int pFirstIdx = min(pstart[parent],npParent-1);
            Int npIdeal =(Int)std::round(npParent*proportion);
            Int numProcs = max(1,min(npParent-pFirstIdx,npIdeal));
            Int pFirst = procGroups_[parent].Ranks().at(pFirstIdx);

            logfileptr->OFS()<<I<<" npParent = "<<npParent<<" pstartParent = "<<pstart[parent]<<" childrenParent = "<<children[parent]<<" pFirst = "<<pFirst<<" numProcs = "<<numProcs<<" proportion = "<<proportion<<endl; 
            pstart[parent]+= numProcs;

            if(npParent!=numProcs){
            if(iam>=pFirst && iam<pFirst+numProcs){
              levelGroups_.push_back(ProcGroup());//reserve(numLevels);
              levelGroups_.back().Ranks().reserve(numProcs);
              for(Int p = pFirst; p<pFirst+numProcs;++p){ levelGroups_.back().Ranks().push_back(p);}
              groupIdx_[I] = levelGroups_.size()-1;
              logfileptr->OFS()<<"DEBUG "<<I<<" = "<<groupIdx_[I]<<" | ";for(Int p = pFirst; p<pFirst+numProcs;++p){ logfileptr->OFS()<<p<<" ";} logfileptr->OFS()<<endl;
            }
            }
            else{
              groupIdx_[I] = groupIdx_[parent];
              logfileptr->OFS()<<"DEBUG "<<I<<" = "<<groupIdx_[I]<<" | ";for(Int p = pFirst; p<pFirst+numProcs;++p){ logfileptr->OFS()<<p<<" ";} logfileptr->OFS()<<endl;
            }

            procGroups_[I].Ranks().reserve(numProcs);
            for(Int p = pFirst; p<pFirst+numProcs;++p){ procGroups_[I].Ranks().push_back(p);}

            pstart[parent]++;

            logfileptr->OFS()<<I<<": "; 
            for(Int i = 0; i<procGroups_[I].Ranks().size();++i){logfileptr->OFS()<<procGroups_[I].Ranks().at(i)<<" ";}
            logfileptr->OFS()<<endl;
          }


          //now choose which processor to get
          std::vector<double> load(np,0.0);
          for(Int I=1;I<=supETree_.Size();I++){
            Int minLoadP= -1;
            double minLoad = -1;
            for(Int i = 0; i<procGroups_[I].Ranks().size();++i){
              Int proc = procGroups_[I].Ranks().at(i);
              if(load[proc]<minLoad || minLoad==-1){
                minLoad = load[proc];
                minLoadP = proc;
              }
            }

            groupWorker_[I] = minLoadP;
            //procGroups_[I].Worker() = minLoadP;


            Int fc = Xsuper_[I-1];
            Int width = Xsuper_[I] - Xsuper_[I-1];
            Int height = cc_[fc-1];
            //cost of factoring current node + updating ancestors (how about fanboth ?)
            //this is exactly the cost of FAN-In while 
            //CHOLESKY_COST(height,width) + SUM_child SubtreeLoads - SUM_descendant CHOLESKY_COSTS
            // would be the cost of Fan-Out
            double local_load = NodeLoad[I];//width*height*height;

            load[minLoadP]+=local_load;
          }


          logfileptr->OFS()<<"Proc load: "<<load<<endl;


          for(Int I = 1; I<procGroups_.size();++I){ 
            //procMap_[I-1] = procGroups_[I].Worker();
            procMap_[I-1] = groupWorker_[I];
          }

          for(Int I = 0; I<levelGroups_.size();++I){ 
            logfileptr->OFS()<<"Group "<<I<<": "<<std::endl<<levelGroups_[I]<<endl;
          }


        }

        return procMap_;
  }

};

//TODO
class SubtreeToSubcubeVolume: public TreeLoadBalancer{
  protected:
    bool fan_in_;
    IntNumVec & Xsuper_;
    IntNumVec & SupMembership_;
    PtrVec & Xlindx_;
    IdxVec & Lindx_;
    IntNumVec & cc_;


    double factor_cost(Int m, Int n){
      return 0;//CHOLESKY_COST(m,n);
    }

    double update_cost(Int m, Int n, Int k){
      return m*n;
    }


  public:
      SubtreeToSubcubeVolume(Int np, ETree & supETree,IntNumVec & Xsuper, IntNumVec & SupMembership,PtrVec & Xlindx, IdxVec & Lindx, IntNumVec & pCc, bool fan_in = true):Xsuper_(Xsuper),SupMembership_(SupMembership),cc_(pCc),Xlindx_(Xlindx),Lindx_(Lindx),TreeLoadBalancer(np,supETree){
          fan_in_=fan_in;
      };


  virtual inline std::vector<Int> & GetMap(){
      if(procGroups_.size()==0){
          //compute number of children and load
          Int numLevels = 1;
          vector<double> SubTreeLoad(supETree_.Size()+1,0.0);
          vector<double> NodeLoad(supETree_.Size()+1,0.0);
          vector<Int> children(supETree_.Size()+1,0);
          for(Int I=1;I<=supETree_.Size();I++){
            Int parent = supETree_.Parent(I-1);
            ++children[parent];
            Idx fc = Xsuper_[I-1];
            Int width = Xsuper_[I] - Xsuper_[I-1];
            Int height = cc_[fc-1];
            //cost of factoring curent panel
            double local_load = factor_cost(height,width);
            SubTreeLoad[I]+=local_load;
            NodeLoad[I]+=local_load;
            //cost of updating ancestors
            {
              Ptr fi = Xlindx_[I-1];
              Ptr li = Xlindx_[I]-1;

              //parse rows
              Int tgt_snode_id = I;
              Idx tgt_fr = fc;
              Idx tgt_lr = fc;
              Ptr tgt_fr_ptr = 0;
              Ptr tgt_lr_ptr = 0;
              for(Ptr rowptr = fi; rowptr<=li;rowptr++){
                Idx row = Lindx_[rowptr-1]; 
                if(SupMembership_[row-1]==tgt_snode_id){ continue;}

                //we have a new tgt_snode_id
                tgt_lr_ptr = rowptr-1;
                tgt_lr = Lindx_[rowptr-1 -1];
                if(tgt_snode_id!=I){
                  Int m = li-tgt_fr_ptr+1;
                  Int n = width;
                  Int k = tgt_lr_ptr - tgt_fr_ptr+1; 
                  double cost = update_cost(m,n,k);
                  if(fan_in_){
                    SubTreeLoad[I]+=cost;
                    NodeLoad[I]+=cost;
                  }
                  else{
                    SubTreeLoad[tgt_snode_id]+=cost;
                    NodeLoad[tgt_snode_id]+=cost;
                  }
                }
                tgt_fr = row;
                tgt_fr_ptr = rowptr;
                tgt_snode_id = SupMembership_[row-1];
              }
              //do the last update supernode
              //we have a new tgt_snode_id
              tgt_lr_ptr = li;
              tgt_lr = Lindx_[li -1];
              if(tgt_snode_id!=I){
                Int m = li-tgt_fr_ptr+1;
                Int n = width;
                Int k = tgt_lr_ptr - tgt_fr_ptr+1; 
                double cost = update_cost(m,n,k);
                  if(fan_in_){
                    SubTreeLoad[I]+=cost;
                    NodeLoad[I]+=cost;
                  }
                  else{
                    SubTreeLoad[tgt_snode_id]+=cost;
                    NodeLoad[tgt_snode_id]+=cost;
                  }
                }
            }
            SubTreeLoad[parent]+=SubTreeLoad[I];
          }

          vector<Int> levels(supETree_.Size()+1,0);
          for(Int I=n_; I>= 1;I--){ 
            Int parent = supETree_.Parent(I-1);
            if(parent==0){levels[I]=0;}
            else{ levels[I] = levels[parent]+1; numLevels = max(numLevels,levels[I]);}
          }
          numLevels++;
            logfileptr->OFS()<<"levels : "; 
            for(Int i = 0; i<levels.size();++i){logfileptr->OFS()<<levels.at(i)<<" ";}
            logfileptr->OFS()<<endl;


          logfileptr->OFS()<<"SubTreeLoad is "<<SubTreeLoad<<endl;


          //procmaps[0]/pstart[0] represents the complete list
          procGroups_.resize(n_+1);
          vector<Int> pstart(n_+1,0);
          procGroups_[0].Ranks().reserve(np);
          for(Int p = 0;p<np;++p){ procGroups_[0].Ranks().push_back(p);}


          levelGroups_.push_back(ProcGroup());//reserve(numLevels);
          levelGroups_[0].Ranks().reserve(np);
          for(Int p = 0;p<np;++p){levelGroups_[0].Ranks().push_back(p);}


          for(Int I=n_; I>= 1;I--){ 
            Int parent = supETree_.Parent(I-1);

            //split parent's proc list
            double parent_load = 0.0;

            if(parent!=0){
              Int fc = Xsuper_[parent-1];
              Int width = Xsuper_[parent] - Xsuper_[parent-1];
              Int height = cc_[fc-1];
              parent_load = NodeLoad[parent];// factor_cost(height,width);
              procGroups_[parent].Load() = parent_load;
            }


            double proportion = min(1.0,SubTreeLoad[I]/(SubTreeLoad[parent]-parent_load));
            Int npParent = procGroups_[parent].Ranks().size();
            Int pFirstIdx = min(pstart[parent],npParent-1);
            Int npIdeal =(Int)std::round(npParent*proportion);
            Int numProcs = max(1,min(npParent-pFirstIdx,npIdeal));
            Int pFirst = procGroups_[parent].Ranks().at(pFirstIdx);

            logfileptr->OFS()<<I<<" npParent = "<<npParent<<" pstartParent = "<<pstart[parent]<<" childrenParent = "<<children[parent]<<" pFirst = "<<pFirst<<" numProcs = "<<numProcs<<" proportion = "<<proportion<<endl; 
            pstart[parent]+= numProcs;

            if(npParent!=numProcs){
            if(iam>=pFirst && iam<pFirst+numProcs){
              levelGroups_.push_back(ProcGroup());//reserve(numLevels);
              levelGroups_.back().Ranks().reserve(numProcs);
              for(Int p = pFirst; p<pFirst+numProcs;++p){ levelGroups_.back().Ranks().push_back(p);}
              groupIdx_[I] = levelGroups_.size()-1;
              logfileptr->OFS()<<"DEBUG "<<I<<" = "<<groupIdx_[I]<<" | ";for(Int p = pFirst; p<pFirst+numProcs;++p){ logfileptr->OFS()<<p<<" ";} logfileptr->OFS()<<endl;
            }
            }
            else{
              groupIdx_[I] = groupIdx_[parent];
              logfileptr->OFS()<<"DEBUG "<<I<<" = "<<groupIdx_[I]<<" | ";for(Int p = pFirst; p<pFirst+numProcs;++p){ logfileptr->OFS()<<p<<" ";} logfileptr->OFS()<<endl;
            }

            procGroups_[I].Ranks().reserve(numProcs);
            for(Int p = pFirst; p<pFirst+numProcs;++p){ procGroups_[I].Ranks().push_back(p);}

            pstart[parent]++;

            logfileptr->OFS()<<I<<": "; 
            for(Int i = 0; i<procGroups_[I].Ranks().size();++i){logfileptr->OFS()<<procGroups_[I].Ranks().at(i)<<" ";}
            logfileptr->OFS()<<endl;
          }


          //now choose which processor to get
          std::vector<double> load(np,0.0);
          for(Int I=1;I<=supETree_.Size();I++){
            Int minLoadP= -1;
            double minLoad = -1;
            for(Int i = 0; i<procGroups_[I].Ranks().size();++i){
              Int proc = procGroups_[I].Ranks().at(i);
              if(load[proc]<minLoad || minLoad==-1){
                minLoad = load[proc];
                minLoadP = proc;
              }
            }

            groupWorker_[I] = minLoadP;
            //procGroups_[I].Worker() = minLoadP;


            Int fc = Xsuper_[I-1];
            Int width = Xsuper_[I] - Xsuper_[I-1];
            Int height = cc_[fc-1];
            //cost of factoring current node + updating ancestors (how about fanboth ?)
            //this is exactly the cost of FAN-In while 
            //CHOLESKY_COST(height,width) + SUM_child SubtreeLoads - SUM_descendant CHOLESKY_COSTS
            // would be the cost of Fan-Out
            double local_load = NodeLoad[I];//width*height*height;

            load[minLoadP]+=local_load;
          }


          logfileptr->OFS()<<"Proc load: "<<load<<endl;


          for(Int I = 1; I<procGroups_.size();++I){ 
            //procMap_[I-1] = procGroups_[I].Worker();
            procMap_[I-1] = groupWorker_[I];
          }

          for(Int I = 0; I<levelGroups_.size();++I){ 
            logfileptr->OFS()<<"Group "<<I<<": "<<std::endl<<levelGroups_[I]<<endl;
          }


        }

        return procMap_;
  }


//  virtual inline std::vector<Int> & GetMap(){
//      if(procGroups_.size()==0){
//
////////gdb_lock();
//////
//////          //compute number of children and load
//////          Int numLevels = 1;
//////          vector<double> SubTreeLoad(supETree_.Size()+1,0.0);
//////          vector<double> NodeLoad(supETree_.Size()+1,0.0);
//////          vector<Int> children(supETree_.Size()+1,0);
//////          for(Int I=1;I<=supETree_.Size();I++){
//////            Int parent = supETree_.Parent(I-1);
//////            ++children[parent];
//////            Idx fc = Xsuper_[I-1];
//////            Int width = Xsuper_[I] - Xsuper_[I-1];
//////            Int height = cc_[fc-1];
//////            //cost of factoring curent panel
//////            double local_load = 0;//CHOLESKY_COST(height,width);
//////            SubTreeLoad[I]+=local_load;
//////            NodeLoad[I]+=local_load;
//////            //cost of updating ancestors
//////#if 1
//////            {
//////              Ptr fi = Xlindx_[I-1];
//////              Ptr li = Xlindx_[I]-1;
//////
//////              //parse rows
//////              Int tgt_snode_id = I;
//////              Idx tgt_fr = fc;
//////              Idx tgt_lr = fc;
//////              Ptr tgt_fr_ptr = 0;
//////              Ptr tgt_lr_ptr = 0;
//////              for(Ptr rowptr = fi; rowptr<=li;rowptr++){
//////                Idx row = Lindx_[rowptr-1]; 
//////                if(SupMembership_[row-1]==tgt_snode_id){ continue;}
//////
//////                //we have a new tgt_snode_id
//////                tgt_lr_ptr = rowptr-1;
//////                tgt_lr = Lindx_[rowptr-1 -1];
//////                if(tgt_snode_id!=I){
//////                  Int m = li-tgt_fr_ptr+1;
//////                  Int n = width;
//////                  Int k = tgt_lr_ptr - tgt_fr_ptr+1; 
//////                  //update cost is the volume
//////                  double update_cost = m*k;
//////                  SubTreeLoad[I]+=update_cost;
//////                  NodeLoad[I]+=update_cost;
//////                }
//////                tgt_fr = row;
//////                tgt_fr_ptr = rowptr;
//////                tgt_snode_id = SupMembership_[row-1];
//////              }
//////              //do the last update supernode
//////              //we have a new tgt_snode_id
//////              tgt_lr_ptr = li;
//////              tgt_lr = Lindx_[li -1];
//////              if(tgt_snode_id!=I){
//////                Int m = li-tgt_fr_ptr+1;
//////                Int n = width;
//////                Int k = tgt_lr_ptr - tgt_fr_ptr+1; 
//////                //update cost is the volume
//////                double update_cost = m*k;
//////                SubTreeLoad[I]+=update_cost;
//////                NodeLoad[I]+=update_cost;
//////              }
//////            }
//////            SubTreeLoad[parent]+=SubTreeLoad[I];
//////#else
//////            SubTreeLoad[parent]+=SubTreeLoad[I];
//////#endif
//////          }
//////
//////          vector<Int> levels(supETree_.Size()+1,0);
//////          for(Int I=n_; I>= 1;I--){ 
//////            Int parent = supETree_.Parent(I-1);
//////            if(parent==0){levels[I]=0;}
//////            else{ levels[I] = levels[parent]+1; numLevels = max(numLevels,levels[I]);}
//////          }
//////          numLevels++;
//////            logfileptr->OFS()<<"levels : "; 
//////            for(Int i = 0; i<levels.size();++i){logfileptr->OFS()<<levels.at(i)<<" ";}
//////            logfileptr->OFS()<<endl;
//////
//////
//////          logfileptr->OFS()<<"SubTreeLoad is "<<SubTreeLoad<<endl;
//////
//////
//////          //procmaps[0]/pstart[0] represents the complete list
//////          procGroups_.resize(n_+1);
//////          vector<Int> pstart(n_+1,0);
//////          procGroups_[0].Ranks().reserve(np);
//////          for(Int p = 0;p<np;++p){ procGroups_[0].Ranks().push_back(p);}
//////
//////
//////          levelGroups_.push_back(ProcGroup());//reserve(numLevels);
//////          levelGroups_[0].Ranks().reserve(np);
//////          for(Int p = 0;p<np;++p){levelGroups_[0].Ranks().push_back(p);}
//////
//////
//////          for(Int I=n_; I>= 1;I--){ 
//////            Int parent = supETree_.Parent(I-1);
//////
//////            //split parent's proc list
//////            double parent_load = 0.0;
//////
//////            if(parent!=0){
//////              Int fc = Xsuper_[parent-1];
//////              Int width = Xsuper_[parent] - Xsuper_[parent-1];
//////              Int height = cc_[fc-1];
//////              parent_load = 0;//width*height;//CHOLESKY_COST(height,width);
//////              procGroups_[parent].Load() = parent_load;
//////            }
//////
//////
//////            double proportion = min(1.0,SubTreeLoad[I]/(SubTreeLoad[parent]-parent_load));
//////            Int npParent = procGroups_[parent].Ranks().size();
//////            Int pFirstIdx = min(pstart[parent],npParent-1);
//////            Int npIdeal =(Int)std::round(npParent*proportion);
//////            Int numProcs = max(1,min(npParent-pFirstIdx,npIdeal));
//////            Int pFirst = procGroups_[parent].Ranks().at(pFirstIdx);
//////
//////            logfileptr->OFS()<<I<<" npParent = "<<npParent<<" pstartParent = "<<pstart[parent]<<" childrenParent = "<<children[parent]<<" pFirst = "<<pFirst<<" numProcs = "<<numProcs<<" proportion = "<<proportion<<endl; 
//////            pstart[parent]+= numProcs;
//////
//////            if(npParent!=numProcs){
//////            if(iam<pFirst || iam>=pFirst+numProcs){
//////
//////            }
//////            else{
//////              levelGroups_.push_back(ProcGroup());//reserve(numLevels);
//////              levelGroups_.back().Ranks().reserve(numProcs);
//////              for(Int p = pFirst; p<pFirst+numProcs;++p){ levelGroups_.back().Ranks().push_back(p);}
//////              groupIdx_[I] = levelGroups_.size()-1;
//////              logfileptr->OFS()<<"DEBUG "<<I<<" = "<<groupIdx_[I]<<" | ";for(Int p = pFirst; p<pFirst+numProcs;++p){ logfileptr->OFS()<<p<<" ";} logfileptr->OFS()<<endl;
//////            }
//////            }
//////            else{
//////              groupIdx_[I] = groupIdx_[parent];
//////              logfileptr->OFS()<<"DEBUG "<<I<<" = "<<groupIdx_[I]<<" | ";for(Int p = pFirst; p<pFirst+numProcs;++p){ logfileptr->OFS()<<p<<" ";} logfileptr->OFS()<<endl;
//////            }
//////
//////            procGroups_[I].Ranks().reserve(numProcs);
//////            for(Int p = pFirst; p<pFirst+numProcs;++p){ procGroups_[I].Ranks().push_back(p);}
//////
//////            pstart[parent]++;
//////
//////            logfileptr->OFS()<<I<<": "; 
//////            for(Int i = 0; i<procGroups_[I].Ranks().size();++i){logfileptr->OFS()<<procGroups_[I].Ranks().at(i)<<" ";}
//////            logfileptr->OFS()<<endl;
//////          }
//////
//////
//////          //now choose which processor to get
//////          std::vector<double> load(np,0.0);
//////          for(Int I=1;I<=supETree_.Size();I++){
//////            Int minLoadP= -1;
//////            double minLoad = -1;
//////            for(Int i = 0; i<procGroups_[I].Ranks().size();++i){
//////              Int proc = procGroups_[I].Ranks().at(i);
//////              if(load[proc]<minLoad || minLoad==-1){
//////                minLoad = load[proc];
//////                minLoadP = proc;
//////              }
//////            }
//////
//////            procGroups_[I].Worker() = minLoadP;
//////
//////
//////            Int fc = Xsuper_[I-1];
//////            Int width = Xsuper_[I] - Xsuper_[I-1];
//////            Int height = cc_[fc-1];
//////            //cost of factoring current node + updating ancestors (how about fanboth ?)
//////            //this is exactly the cost of FAN-In while 
//////            //CHOLESKY_COST(height,width) + SUM_child SubtreeLoads - SUM_descendant CHOLESKY_COSTS
//////            // would be the cost of Fan-Out
//////            double local_load = NodeLoad[I];//width*height*height;
//////
//////            load[minLoadP]+=local_load;
//////          }
//////
//////
//////          logfileptr->OFS()<<"Proc load: "<<load<<endl;
//////
//////
//////          for(Int I = 1; I<procGroups_.size();++I){ 
//////            procMap_[I-1] = procGroups_[I].Worker();
//////          }
//////
//////          for(Int I = 0; I<levelGroups_.size();++I){ 
//////            logfileptr->OFS()<<"Group "<<I<<": "<<std::endl<<levelGroups_[I]<<endl;
//////          }
//////
//////
//        }
//
//        return procMap_;
//  }

};


}

#endif //_LOAD_BALANCER_HPP_
