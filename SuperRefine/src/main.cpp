#include <iostream>
#include <cmath>
#include <complex>
#include <cassert>
#include <limits>

#include "iohb.h"

#include "ngchol/Utility.hpp"
#include "ngchol/DistSparseMatrix.hpp"
#include "ngchol/SparseMatrixStructure.hpp"
#include "ngchol/PastixReorder.hpp"
#include "ngchol/timer.hpp"

#define FUNDAMENTAL
//#define RELAX_SNODE

struct Stats{
  int64_t totalSnodeBlocks = 0;
  int64_t totalBlocks = 0;
  double blocksPerSnode = 0;
  double blocksPerCol = 0;
  double avgSnodeBlockSize = 0;
  double avgBlockSize = 0;
  void reset(){
   totalSnodeBlocks = 0;
   totalBlocks = 0;
   blocksPerSnode = 0;
   blocksPerCol = 0;
   avgSnodeBlockSize = 0;
   avgBlockSize = 0;
  }
  void print(){
   cout<<"totalSnodeBlocks: "<<totalSnodeBlocks<<endl;
   cout<<"totalBlocks: "<<totalBlocks<<endl;
   cout<<"blocksPerSnode: "<<blocksPerSnode<<endl;
   cout<<"blocksPerCol: "<<blocksPerCol<<endl;
   cout<<"avgSnodeBlockSize: "<<avgSnodeBlockSize<<endl;
   cout<<"avgBlockSize: "<<avgBlockSize<<endl;
  }
};

void countBlock(vector<int> & Xsuper, vector<int64_t> & Xlindx, vector<int32_t> & Lindx, Stats & stats){
        stats.reset();

        int64_t supNrows = 0;
        for(int I = 1; I<Xlindx.size();++I){
#ifdef VERBOSE1
          cout<<I<<"("<<Xsuper[I-1]<<".."<<Xsuper[I]-1<<"): ";
#endif
                  //count number of contiguous blocks (at least one diagonal block)
                  int32_t fc = Xsuper[I-1];
                  int32_t lc = Xsuper[I]-1;
                  int32_t fi = Xlindx[I-1];
                  int32_t li = Xlindx[I]-1;
                  int32_t iPrevRow = Lindx[fi-1]-1;
                  int32_t iFirstRow = Lindx[fi-1];
              
                  int32_t width = lc - fc + 1; 
 
                  for(int col = fc; col<=lc;col++){ 
                    //1 to count the diagonal block, 0 to skip it
                    int32_t nzBlockCnt = 0;//1;
                    int32_t height = li - fi + 1;
                    for(int64_t idx = fi; idx<=li;idx++){
                      int32_t iRow = Lindx[idx-1];
                      if(iRow<col){
                        --height; 
                      }
                      //enforce the first block to be a square diagonal block
                      if(nzBlockCnt==1 && iRow>col){
                        nzBlockCnt++;
                        stats.avgBlockSize+=col-iFirstRow+1;
                        if(col==fc){
                          stats.avgSnodeBlockSize+=width;
                        }
                        iFirstRow=iRow;
#ifdef VERBOSE1
                        cout<<"| ";
#endif
                      }
                      else if(iRow!=iPrevRow+1){
                        nzBlockCnt++;
                        stats.avgBlockSize+=iPrevRow-iFirstRow+1;
                        if(col==fc){
                          stats.avgSnodeBlockSize+=iPrevRow-iFirstRow+1;
                        }
                        iFirstRow=iRow;
#ifdef VERBOSE1
                        cout<<"| ";
#endif
                      }
#ifdef VERBOSE1
                      cout<<iRow<<" ";
#endif

                      iPrevRow=iRow;
                    }
#ifdef VERBOSE1
                    cout<<" || "<<nzBlockCnt<<endl;
#endif


                    stats.totalBlocks+=nzBlockCnt;
                    if(col==lc){
                      stats.totalSnodeBlocks+=nzBlockCnt;
                      supNrows+=height;
                    }
                  }
          }

        stats.avgBlockSize/=stats.totalBlocks;
        stats.avgSnodeBlockSize/=stats.totalSnodeBlocks;

        stats.blocksPerCol=stats.totalBlocks/(Xsuper.back()-1);
        stats.blocksPerSnode=stats.totalSnodeBlocks/(Xsuper.size()-1);
        cout<<"N is "<<Xsuper.back()-1<<endl;
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

    Global.ExpandSymmetric();
    if(iam==0){ cout<<"Matrix expanded"<<endl; }
    //Create an Ordering object to hold the permutation
    Ordering order;
    order.SetStructure(Global);
    if(iam==0){ cout<<"Structure set"<<endl; }
{

TIMER_START(ORDERING);
double tstart = MPI_Wtime();
    //Reoder the matrix with MMD
    if(argc>2){
      string ordering(argv[2]);
      if(ordering == "AMD"){
        order.AMD();
      }
      else if (ordering == "METIS"){
        order.METIS();
      }
      else if (ordering == "SCOTCH"){
        order.SCOTCH();
      }
      else{
        order.MMD();
      }
    }
    else{
      order.MMD();
    }

double tstop = MPI_Wtime();
TIMER_STOP(ORDERING);
    if(iam==0){ cout<<"ordering done in "<<tstop-tstart<<endl; }
}
#ifdef VERBOSE
    cout<<"MMD perm: "<<order.perm<<endl;
#endif

    ETree Etree;
    Etree.ConstructETree(Global,order);
    Etree.PostOrderTree(order);
    if(iam==0){ cout<<"ETree done"<<endl; }
    //cout<<"etree "<<Etree<<std::endl;
    vector<int> cc,rc;
    Global.GetLColRowCount(Etree,order,cc,rc);
    //cout<<"colcnt "<<cc<<std::endl;
    Etree.SortChildren(cc,order);
#ifdef VERBOSE
    cout<<"colcnt sorted "<<cc<<std::endl;
#endif

    double flops = 0.0;
    for(int i = 0; i<cc.size();++i){
      flops+= (double)pow((double)cc[i],2.0);
    }

    if(iam==0){
      cout<<"Flops: "<<flops<<endl;
    }

#ifdef FUNDAMENTAL
    Global.FindFundamentalSupernodes(Etree,order,cc,SupMembership,Xsuper,maxSnode);
#else
    Global.FindSupernodes(Etree,order,cc,SupMembership,Xsuper,maxSnode);
#endif
    if(iam==0){ cout<<"Supernodes found"<<endl; }

#ifdef RELAX_SNODE
    Global.RelaxSupernodes(Etree, cc,SupMembership, Xsuper, maxSnode );
    if(iam==0){ cout<<"Relaxation done"<<endl; }
    Global.SymbolicFactorizationRelaxed(Etree,order,cc,Xsuper,SupMembership,Xlindx,Lindx);
#else
    Global.SymbolicFactorization(Etree,order,cc,Xsuper,SupMembership,Xlindx,Lindx);
    if(iam==0){ cout<<"Symbfact done"<<endl; }
#endif

      cout<<"Non refined perm is:"<<order.perm<<endl;
#ifdef VERBOSE
      cout<<"Xsuper is "<<Xsuper<<endl;
      //cout<<"Xlindx is "<<Xlindx<<endl;
      //cout<<"Lindx is "<<Lindx<<endl;
#endif




    Ordering orderSave = order;

    vector<int> permRefined,origPerm,newXsuper;
{
double tstart = MPI_Wtime();
#if 0
    Global.RefineSupernodes(Etree,order, SupMembership, Xsuper, Xlindx, Lindx, permRefined,origPerm);
#else
    ETree SupETree = Etree.ToSupernodalETree(Xsuper,SupMembership,order);
    SymbolMatrix * symbmtx = GetPastixSymbolMatrix(Xsuper,SupMembership, Xlindx, Lindx);
    Order * psorder = GetPastixOrder(symbmtx,Xsuper, SupETree, &order.perm[0], &order.invp[0]);
    symbolReordering( symbmtx, psorder, 0, std::numeric_limits<int>::max(), 0 );
#endif
double tstop = MPI_Wtime();
    if(iam==0){ cout<<"Refinement done in "<<tstop-tstart<<endl; }
}


    if(iam==0){ 
      cout<<"Refined perm is:"<<order.perm<<endl;
#ifdef VERBOSE
      cout<<"Original Supernodes: "<<Xsuper<<endl; 
#endif
    }

int nsuper = Xsuper.size()-1;
cout<<"=========================================="<<endl;
        Stats stats;
        countBlock(Xsuper, Xlindx, Lindx,stats);

{
vector<int64_t> dummy;
Xlindx.swap(dummy);
vector<int32_t> dummy2;
Lindx.swap(dummy2);
}

    {
#if 0
      //safety check
      {
        ETree Etree2;
        Etree2.ConstructETree(Global,orderSave);
        Etree2.PostorderTree(orderSave);
        vector<int> cc2,rc2;
        Global.GetLColRowCount(Etree2,orderSave,cc2,rc2);
      Etree2.SortChildren(cc2,order);

        for(int i =0; i<cc.size();++i){
          assert(cc[i]==cc2[i]);
        }


      }
#endif

      vector<int> SupMembership2,Xsuper2;
      vector<int64_t> Xlindx2;
      vector<int32_t> Lindx2;

      ETree Etree2;
      Etree2.ConstructETree(Global,order);
      Etree2.PostOrderTree(order);
      vector<int> cc2,rc2;
      Global.GetLColRowCount(Etree2,order,cc2,rc2);
      Etree2.SortChildren(cc2,order);

      double flops2 = 0.0;
      for(int i = 0; i<cc2.size();++i){
        flops2+= (double)pow((double)cc2[i],2.0);
      }

      if(iam==0){
        cout<<"Flops 2: "<<flops2<<endl;
      }

#ifdef FUNDAMENTAL
      Global.FindFundamentalSupernodes(Etree2,order,cc2,SupMembership2,Xsuper2,maxSnode);
#else
      Global.FindSupernodes(Etree2,order,cc2,SupMembership2,Xsuper2,maxSnode);
#endif

      if(iam==0){ cout<<"Supernodes found"<<endl; }

#ifdef RELAX_SNODE
      Global.RelaxSupernodes(Etree2, cc2,SupMembership2, Xsuper2, maxSnode );
      if(iam==0){ cout<<"Relaxation done"<<endl; }
#ifdef VERBOSE
      if(iam==0){ cout<<"New Supernodes: "<<Xsuper2<<endl; }
#endif
      Global.SymbolicFactorizationRelaxed(Etree2,order,cc2,Xsuper2,SupMembership2,Xlindx2,Lindx2);
#else

#ifdef VERBOSE
      if(iam==0){ cout<<"New Supernodes: "<<Xsuper2<<endl; }
#endif
      Global.SymbolicFactorization(Etree2,order,cc2,Xsuper2,SupMembership2,Xlindx2,Lindx2);
      if(iam==0){ cout<<"Symbfact done"<<endl; }
#endif




      if(iam==0){

cout<<"=========================================="<<endl;
        Stats statsRefined;
        countBlock(Xsuper2, Xlindx2, Lindx2,statsRefined);

{
vector<int64_t> dummy;
Xlindx2.swap(dummy);
vector<int32_t> dummy2;
Lindx2.swap(dummy2);
}

cout<<"=========================================="<<endl;
        stats.print();
cout<<"=========================================="<<endl;
        statsRefined.print();
cout<<"=========================================="<<endl;

        int64_t nnz = 0;
        for(int i =0; i<cc.size();++i){
          nnz+=cc[i];
        }
        for(int i =0; i<Xsuper.size()-1;++i){
          int x = Xsuper[i+1]-Xsuper[i];
          nnz+= x*(x-1)/2;
        }
        //      cout<<endl;
        cout<<nnz<<endl;

        int64_t nnz2 = 0;
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

