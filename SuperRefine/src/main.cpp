#include <iostream>
#include <cmath>
#include <complex>

#include "iohb.h"

#include "ngchol/DistSparseMatrix.hpp"
#include "ngchol/SparseMatrixStructure.hpp"


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





#if 0
    DistSparseMatrix<double> HMat(worldcomm);
    sparse_matrix_file_format_t informat;
if(argc>2){
    informat = sparse_matrix_file_format_string_to_enum (argv[2]);
}
else{
    informat = sparse_matrix_file_format_string_to_enum ("HARWELL_BOEING");
}

    sparse_matrix_t* Atmp = load_sparse_matrix (informat, argv[1]);
    sparse_matrix_convert (Atmp, CSC);
    const csc_matrix_t * cscptr = (const csc_matrix_t *) Atmp->repr;
    HMat.CopyData(cscptr->n,cscptr->nnz,cscptr->colptr,cscptr->rowidx,(double *)cscptr->values,true);


  for(int i=0;i<=cscptr->n;++i){
    cout<<cscptr->colptr[i]<<" ";
  }
  cout<<endl;

  for(int i=0;i<cscptr->nnz;++i){
    cout<<cscptr->rowidx[i]<<" ";
  }
  cout<<endl;




    destroy_sparse_matrix (Atmp);

    if(iam==0){ cout<<"Matrix order is "<<HMat.size<<endl; }

    int maxSnode = -1;
    vector<int> SupMembership,Xsuper;
    vector<int64_t> Xlindx;
    vector<int32_t> Lindx;





    SparseMatrixStructure Local = HMat.GetLocalStructure();


  for(int i=0;i<Local.colptr.size();++i){
    cout<<Local.colptr[i]<<" ";
  }
  cout<<endl;

  for(int i=0;i<Local.rowind.size();++i){
    cout<<Local.rowind[i]<<" ";
  }
  cout<<endl;




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
#else

{
  int M, N, nonzeros, Nrhs;
  int *colptr, *rowind;
  double *val;
  char *Type;

  readHB_info(argv[1], &M, &N, &nonzeros, &Type, &Nrhs);
  readHB_newmat_double(argv[1], &M, &N, &nonzeros, &colptr, &rowind, &val);


  SparseMatrixStructure Local;
  if ( Type[0] == 'R' ){

    DistSparseMatrix<double> HMat(worldcomm);
    HMat.CopyData(M,nonzeros,colptr,rowind,val,true);
    free(colptr);
    free(rowind);
    free(val);


    if(iam==0){ cout<<"Matrix order is "<<HMat.size<<endl; }

    Local = HMat.GetLocalStructure();
  }
  else{
    DistSparseMatrix<complex<double> > HMat(worldcomm);
    HMat.CopyData(N,nonzeros,colptr,rowind,(complex<double> *)val);
    free(colptr);
    free(rowind);
    free(val);

    if(iam==0){ cout<<"Matrix order is "<<HMat.size<<endl; }

    Local = HMat.GetLocalStructure();

  }


  int maxSnode = -1;
  vector<int> SupMembership,Xsuper;
  vector<int64_t> Xlindx;
  vector<int32_t> Lindx;


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

#ifdef RELAX_SNODE
  Global.RelaxSupernodes(Etree, cc,SupMembership, Xsuper, maxSnode );
  if(iam==0){ cout<<"Relaxation done"<<endl; }
  Global.SymbolicFactorizationRelaxed(Etree,Order,cc,Xsuper,SupMembership,Xlindx,Lindx);
#else
  Global.SymbolicFactorization(Etree,Order,cc,Xsuper,SupMembership,Xlindx,Lindx);
  if(iam==0){ cout<<"Symbfact done"<<endl; }
#endif

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




  vector<int> permRefined,origPerm,newXsuper;
  Global.RefineSupernodes(Etree,Order, SupMembership, Xsuper, Xlindx, Lindx, permRefined,origPerm);
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

  //safety check
  {
    ETree Etree2;
    Etree2.ConstructETree(Global,Order);
    Etree2.PostOrderTree(Order);
    vector<int> cc2,rc2;
    Global.GetLColRowCount(Etree2,Order,cc2,rc2);

    if(iam==0){

      int nnz = 0;
      for(int i =0; i<cc.size();++i){
        cout<<cc[i]<<" "; 
        nnz+=cc[i];
      }
      cout<<endl;
      cout<<nnz<<endl;

      int nnz2 = 0;
      for(int i =0; i<cc2.size();++i){
        cout<<cc2[i]<<" "; 
        nnz2+=cc2[i];
      }
      cout<<endl;
      cout<<nnz2<<endl;
    }
  }





}
#endif 











//
//#ifdef _DEBUG_
//    logfileptr->OFS()<<"Membership list is "<<SupMembership_<<std::endl;
//    logfileptr->OFS()<<"xsuper "<<Xsuper_<<std::endl;
//#endif
//









MPI_Finalize();

}
