#include <iostream>
#include <cmath>


#include "ngchol/DistSparseMatrix.hpp"
#include "ngchol/SparseMatrixStructure.hpp"

extern "C" {
#include "bebop/util/config.h"
#include "bebop/smc/sparse_matrix.h"
#include "bebop/smc/csr_matrix.h"
#include "bebop/smc/csc_matrix.h"
#include "bebop/smc/sparse_matrix_ops.h"

#include "bebop/util/get_options.h"
#include "bebop/util/init.h"
#include "bebop/util/malloc.h"
#include "bebop/util/timer.h"
#include "bebop/util/util.h"
}


using namespace LIBCHOLESKY;

int main(int argc, char **argv) 
{
  MPI_Init(&argc,&argv);

  int iam=0;
  int np=1;
  MPI_Comm worldcomm;
  MPI_Comm_dup(MPI_COMM_WORLD,&worldcomm);
  //  int np, iam;
  MPI_Comm_size(worldcomm,&np);
  MPI_Comm_rank(worldcomm,&iam);





    sparse_matrix_file_format_t informat;
    DistSparseMatrix<double> HMat(worldcomm);
if(argc>2){
    informat = sparse_matrix_file_format_string_to_enum (argv[2]);
}
else{
    informat = sparse_matrix_file_format_string_to_enum ("HARWELL_BOEING");
}

    sparse_matrix_t* Atmp = load_sparse_matrix (informat, argv[1]);
    sparse_matrix_convert (Atmp, CSC);
    const csc_matrix_t * cscptr = (const csc_matrix_t *) Atmp->repr;
    HMat.CopyData(cscptr->n,cscptr->nnz,cscptr->colptr,cscptr->rowidx,(double *)cscptr->values);
      destroy_sparse_matrix (Atmp);

    if(iam==0){ cout<<"Matrix order is "<<HMat.size<<endl; }

    int maxSnode = -1;
    vector<int> SupMembership,Xsuper;
    vector<int64_t> Xlindx;
    vector<int32_t> Lindx;





    SparseMatrixStructure Local = HMat.GetLocalStructure();
    SparseMatrixStructure Global;
    Local.ToGlobal(Global,worldcomm);
    Global.ExpandSymmetric();
    if(iam==0){ cout<<"Matrix expanded"<<endl; }





    //Create an Ordering object to hold the permutation
    Ordering Order;
    Order.SetStructure(Global);
    if(iam==0){ cout<<"Structure set"<<endl; }

        //Reoder the matrix with MMD
        Order.MMD();
//        Order.AMD();
    if(iam==0){ cout<<"Ordering done"<<endl; }

    ETree Etree;
    Etree.ConstructETree(Global,Order);
    Etree.PostOrderTree(Order);
    if(iam==0){ cout<<"ETree done"<<endl; }
    vector<int> cc,rc;
    Global.GetLColRowCount(Etree,Order,cc,rc);
    Etree.SortChildren(cc,Order);

    if(iam==0){
    //  cout<<"colcnt "<<cc<<std::endl;
    //  cout<<"rowcnt "<<rc<<std::endl;
    }

    double flops = 0.0;
    for(int i = 0; i<cc.size();++i){
      flops+= (double)pow((double)cc[i],2.0);
    }

    if(iam==0){
      cout<<"Flops: "<<flops<<endl;
    }

    Global.FindSupernodes(Etree,Order,cc,SupMembership,Xsuper,maxSnode);
    if(iam==0){ cout<<"Supernodes found"<<endl; }

    Global.RelaxSupernodes(Etree, cc,SupMembership, Xsuper, maxSnode );

    if(iam==0){ 
      cout<<"1 Xsuper is:"<<endl;
      for(int i =0; i<Xsuper.size();++i){
        cout<<Xsuper[i]<<" "; 
      }
      cout<<endl;

      cout<<"1 perm is:"<<endl;
      for(int i =0; i<Order.perm.size();++i){
        cout<<Order.perm[i]<<" "; 
      }
      cout<<endl;
    }

    if(iam==0){ cout<<"Relaxation done"<<endl; }
//    Global.SymbolicFactorizationRelaxed(Etree,Order,cc,Xsuper,SupMembership,Xlindx,Lindx);
    Global.SymbolicFactorization(Etree,Order,cc,Xsuper,SupMembership,Xlindx,Lindx);
    if(iam==0){ cout<<"Symbfact done"<<endl; }
    vector<int> permRefined,origPerm,newXsuper;
    Global.RefineSupernodes(Etree,Order, SupMembership, Xsuper, Xlindx, Lindx, permRefined,origPerm,newXsuper);
    if(iam==0){ cout<<"Refinement done"<<endl; }


    if(iam==0){ 
      cout<<"2 Xsuper is:"<<endl;
      for(int i =0; i<Xsuper.size();++i){
        cout<<Xsuper[i]<<" "; 
      }
      cout<<endl;

      cout<<"2 origPerm is:"<<endl;
      for(int i =0; i<origPerm.size();++i){
        cout<<origPerm[i]<<" "; 
      }
      cout<<endl;


      cout<<"2 permRefined is:"<<endl;
      for(int i =0; i<permRefined.size();++i){
        cout<<permRefined[i]<<" "; 
      }
      cout<<endl;


      cout<<"2 Final perm is:"<<endl;
      for(int i =0; i<Order.perm.size();++i){
        cout<<Order.perm[i]<<" "; 
      }
      cout<<endl;

      cout<<"2 orig CC is:"<<endl;
      for(int i =0; i<origPerm.size();++i){
        cout<<cc[origPerm[i]-1]<<" "; 
      }
      cout<<endl;


      cout<<"2 refined CC is:"<<endl;
      for(int i =0; i<permRefined.size();++i){
        cout<<cc[permRefined[i]-1]<<" "; 
      }
      cout<<endl;




    }

//
//#ifdef _DEBUG_
//    logfileptr->OFS()<<"Membership list is "<<SupMembership_<<std::endl;
//    logfileptr->OFS()<<"xsuper "<<Xsuper_<<std::endl;
//#endif
//









MPI_Finalize();

}
