#include <iostream>
#include <cmath>
#include <complex>
#include <cassert>

#include "iohb.h"

#include "ngchol/Utility.hpp"
#include "ngchol/DistSparseMatrix.hpp"
#include "ngchol/SparseMatrixStructure.hpp"

//int countBlock(vector<int> & Xsuper, vector<int64_t> & Xlindx, vector<int32_t> & Lindx){
//        int num_blocks = 0;
//        for(int I = 1; I<Xlindx.size();++I){
//          cout<<I<<"("<<Xsuper[I-1]<<".."<<Xsuper[I]-1<<"): ";
//
//                  //count number of contiguous blocks (at least one diagonal block)
//                  int32_t fc = Xsuper[I-1];
//                  int32_t lc = Xsuper[I]-1;
//                  int32_t fi = Xlindx[I-1];
//                  int32_t li = Xlindx[I]-1;
//                  int32_t nzBlockCnt = 1;
//                  int32_t iPrevRow = Lindx[fi-1]-1;
//                
//                  for(int64_t idx = fi; idx<=li;idx++){
//                    int32_t iRow = Lindx[idx-1];
//    
//                    //enforce the first block to be a square diagonal block
//                    if(nzBlockCnt==1 && iRow>lc){
//                      nzBlockCnt++;
//                      cout<<"| ";
//                    }
//                    else if(iRow!=iPrevRow+1){
//                      nzBlockCnt++;
//                      cout<<"| ";
//                    }
//                    cout<<iRow<<" ";
//
//                    iPrevRow=iRow;
//                  }
//                cout<<" || "<<nzBlockCnt<<endl;
//            num_blocks+=nzBlockCnt;
//          }
//
//
//    return num_blocks;
//}

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
      HMat.CopyData(N,nonzeros,colptr,rowind,(complex<double> *)val,true);
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
    cout<<"Adj: "<<Local.rowind<<endl;
    cout<<"Adj: "<<Global.rowind<<endl;

    Global.ExpandSymmetric();
    if(iam==0){ cout<<"Matrix expanded"<<endl; }
    //Create an Ordering object to hold the permutation
    Ordering Order;
    Order.SetStructure(Global);
    if(iam==0){ cout<<"Structure set"<<endl; }

    //Reoder the matrix with MMD
    Order.MMD();
    //Order.AMD();
    if(iam==0){ cout<<"Ordering done"<<endl; }
    cout<<"MMD perm: "<<Order.perm<<endl;


    ETree Etree;
    Etree.ConstructETree(Global,Order);
    Etree.PostOrderTree(Order);
    if(iam==0){ cout<<"ETree done"<<endl; }
    //cout<<"etree "<<Etree<<std::endl;
    vector<int> cc,rc;
    Global.GetLColRowCount(Etree,Order,cc,rc);
    //cout<<"colcnt "<<cc<<std::endl;
    Etree.SortChildren(cc,Order);
    cout<<"colcnt sorted "<<cc<<std::endl;


    double flops = 0.0;
    for(int i = 0; i<cc.size();++i){
      flops+= (double)pow((double)cc[i],2.0);
    }

    if(iam==0){
      cout<<"Flops: "<<flops<<endl;
    }

//    Global.FindSupernodes(Etree,Order,cc,SupMembership,Xsuper,maxSnode);
    Global.FindFundamentalSupernodes(Etree,Order,cc,SupMembership,Xsuper,maxSnode);
    if(iam==0){ cout<<"Supernodes found"<<endl; }

#ifdef RELAX_SNODE
    Global.RelaxSupernodes(Etree, cc,SupMembership, Xsuper, maxSnode );
    if(iam==0){ cout<<"Relaxation done"<<endl; }
    Global.SymbolicFactorizationRelaxed(Etree,Order,cc,Xsuper,SupMembership,Xlindx,Lindx);
#else
    Global.SymbolicFactorization(Etree,Order,cc,Xsuper,SupMembership,Xlindx,Lindx);
    if(iam==0){ cout<<"Symbfact done"<<endl; }
#endif

#if 1
//VERBOSE
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


      cout<<"perm is:"<<Order.perm<<endl;
      cout<<"Xsuper is "<<Xsuper<<endl;
      cout<<"Xlindx is "<<Xlindx<<endl;
      cout<<"Lindx is "<<Lindx<<endl;








    }
#endif

    Ordering OrderSave = Order;

    vector<int> permRefined,origPerm,newXsuper;
    Global.RefineSupernodes(Etree,Order, SupMembership, Xsuper, Xlindx, Lindx, permRefined,origPerm);
    if(iam==0){ cout<<"Refinement done"<<endl; }


    if(iam==0){ 
//#ifdef VERBOSE
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
//#endif

      if(iam==0){ cout<<"Original Supernodes: "<<Xsuper<<endl; }


    }

    {
#if 0
      //safety check
      {
        ETree Etree2;
        Etree2.ConstructETree(Global,OrderSave);
        Etree2.PostOrderTree(OrderSave);
        vector<int> cc2,rc2;
        Global.GetLColRowCount(Etree2,OrderSave,cc2,rc2);
      Etree2.SortChildren(cc2,Order);

        for(int i =0; i<cc.size();++i){
          assert(cc[i]==cc2[i]);
        }


      }
#endif

      vector<int> SupMembership2,Xsuper2;
      vector<int64_t> Xlindx2;
      vector<int32_t> Lindx2;

      ETree Etree2;
      Etree2.ConstructETree(Global,Order);
      Etree2.PostOrderTree(Order);
      vector<int> cc2,rc2;
      Global.GetLColRowCount(Etree2,Order,cc2,rc2);
      Etree2.SortChildren(cc2,Order);

      double flops2 = 0.0;
      for(int i = 0; i<cc2.size();++i){
        flops2+= (double)pow((double)cc2[i],2.0);
      }

      if(iam==0){
        cout<<"Flops 2: "<<flops2<<endl;
      }

//      Global.FindSupernodes(Etree2,Order,cc2,SupMembership2,Xsuper2,maxSnode);
      Global.FindFundamentalSupernodes(Etree2,Order,cc2,SupMembership2,Xsuper2,maxSnode);
      if(iam==0){ cout<<"Supernodes found"<<endl; }

#ifdef RELAX_SNODE
      Global.RelaxSupernodes(Etree2, cc2,SupMembership2, Xsuper2, maxSnode );
      if(iam==0){ cout<<"Relaxation done"<<endl; }
      if(iam==0){ cout<<"Supernodes: "<<Xsuper2<<endl; }
      Global.SymbolicFactorizationRelaxed(Etree2,Order,cc2,Xsuper2,SupMembership2,Xlindx2,Lindx2);
#else
      if(iam==0){ cout<<"Supernodes: "<<Xsuper2<<endl; }
      Global.SymbolicFactorization(Etree2,Order,cc2,Xsuper2,SupMembership2,Xlindx2,Lindx2);
      if(iam==0){ cout<<"Symbfact done"<<endl; }
#endif




      if(iam==0){


        int num_blocks = 0;
        for(int I = 1; I<Xlindx.size();++I){
          cout<<I<<"("<<Xsuper[I-1]<<".."<<Xsuper[I]-1<<"): ";

                  //count number of contiguous blocks (at least one diagonal block)
                  int32_t fc = Xsuper[I-1];
                  int32_t lc = Xsuper[I]-1;
                  int32_t fi = Xlindx[I-1];
                  int32_t li = Xlindx[I]-1;
                  int32_t nzBlockCnt = 1;
                  int32_t iPrevRow = Lindx[fi-1]-1;
                
                  for(int64_t idx = fi; idx<=li;idx++){
                    int32_t iRow = Lindx[idx-1];
    
                    //enforce the first block to be a square diagonal block
                    if(nzBlockCnt==1 && iRow>lc){
                      nzBlockCnt++;
                      cout<<"| ";
                    }
                    else if(iRow!=iPrevRow+1){
                      nzBlockCnt++;
                      cout<<"| ";
                    }
                    cout<<iRow<<" ";

                    iPrevRow=iRow;
                  }
                cout<<" || "<<nzBlockCnt<<endl;
            num_blocks+=nzBlockCnt;
          }


        int num_blocks2 = 0;
        for(int I = 1; I<Xlindx2.size();++I){
          cout<<I<<"("<<Xsuper2[I-1]<<".."<<Xsuper2[I]-1<<"): ";

                  //count number of contiguous blocks (at least one diagonal block)
                  int32_t fc = Xsuper2[I-1];
                  int32_t lc = Xsuper2[I]-1;
                  int32_t fi = Xlindx2[I-1];
                  int32_t li = Xlindx2[I]-1;
                  int32_t nzBlockCnt = 1;
                  int32_t iPrevRow = Lindx2[fi-1]-1;
                
                  for(int64_t idx = fi; idx<=li;idx++){
                    int32_t iRow = Lindx2[idx-1];
    
                    //enforce the first block to be a square diagonal block
                    if(nzBlockCnt==1 && iRow>lc){
                      nzBlockCnt++;
                      cout<<"| ";
                    }
                    else if(iRow!=iPrevRow+1){
                      nzBlockCnt++;
                      cout<<"| ";
                    }
                    cout<<iRow<<" ";

                    iPrevRow=iRow;
                  }
                cout<<" || "<<nzBlockCnt<<endl;
            num_blocks2+=nzBlockCnt;
          }


//    Global.SymbolicFactorizationRelaxed(Etree,Order,cc,Xsuper,SupMembership,Xlindx,Lindx);




//        for(int I = 1; I<Xlindx.size();++I){
//          cout<<I<<"("<<Xsuper[I-1]<<".."<<Xsuper[I]-1<<"): ";
//          int32_t prevRow = -1;
//          for(int64_t rowi = Xlindx[I-1]; rowi<=Xlindx[I]-1; ++rowi){
//            int32_t row = Lindx[rowi-1];
//            cout<<row<<" ";
//            if(prevRow == -1 || prevRow!=row-1){
//              num_blocks2++;
//            }
//            prevRow=row;
//          }
//          cout<<" | "<<num_blocks2<<endl;
//        }
//
//        int num_blocks = 0;
//        for(int I = 1; I<Xlindx2.size();++I){
//          cout<<I<<"("<<Xsuper2[I-1]<<".."<<Xsuper2[I]-1<<"): ";
//          int32_t prevRow = -1;
//          for(int64_t rowi = Xlindx2[I-1]; rowi<=Xlindx2[I]-1; ++rowi){
//            int32_t row = Lindx2[rowi-1];
//            cout<<row<<" ";
//            if(prevRow == -1 || prevRow!=row-1){
//              num_blocks++;
//            }
//            prevRow=row;
//          }
//          cout<<" | "<<num_blocks<<endl;
//        }

        cout<<"Contiguous blocks: "<<num_blocks<<" vs "<<num_blocks2<<endl;


        int nnz = 0;
        for(int i =0; i<cc.size();++i){
          nnz+=cc[i];
        }
        for(int i =0; i<Xsuper.size()-1;++i){
          int x = Xsuper[i+1]-Xsuper[i];
          nnz+= x*(x-1)/2;
        }
        //      cout<<endl;
        cout<<nnz<<endl;

        int nnz2 = 0;
        for(int i =0; i<cc2.size();++i){
          //        cout<<cc2[i]<<" "; 
          nnz2+=cc2[i];
        }
        for(int i =0; i<Xsuper2.size()-1;++i){
          int x = Xsuper2[i+1]-Xsuper2[i];
          nnz2+= x*(x-1)/2;
        }
        //      cout<<endl;
        cout<<nnz2<<endl;


      }
    }

  }
#endif 

  MPI_Finalize();

}

