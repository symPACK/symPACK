#include "ngchol/PastixReorder.hpp"

#include "ngchol/Utility.hpp"
#include "ngchol/timer.hpp"
#include <cstdio>
#include <cstring>
#include <cassert>
#include <iostream>

#define VERBOSE1




void countBlock(const vector<int> & Xsuper, const vector<int64_t> & Xlindx, const vector<int32_t> & Lindx, const vector<int> & supMembership, vector<int> & blkCount){
  blkCount.assign(Xlindx.size(),0);

  int totalBlocks = 0;
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

    //only go through last column
    for(int col = lc; col<=lc;col++){ 
      //1 to count the diagonal block, 0 to skip it
      int32_t nzBlockCnt = 1;
      int prevFacing = I;
      for(int64_t idx = fi; idx<=li;idx++){
        int32_t iRow = Lindx[idx-1];
        int facing = supMembership[iRow-1];
        //enforce the first block to be a square diagonal block
        if(nzBlockCnt==1 && iRow>col){
          nzBlockCnt++;
          iFirstRow=iRow;
#ifdef VERBOSE1
                        cout<<"| ";
#endif
        }
        else if(iRow!=iPrevRow+1 || prevFacing != facing){
//          if(prevFacing != facing){
//            std::cout<<"BLOCK SPLITTED "<<prevFacing<<" "<<facing<<std::endl;
//          }
          nzBlockCnt++;
          iFirstRow=iRow;
#ifdef VERBOSE1
                        cout<<"| ";
#endif
        }

#ifdef VERBOSE1
                      cout<<iRow<<" ";
#endif
        iPrevRow=iRow;
        prevFacing = facing;
      }

#ifdef VERBOSE1
                    cout<<" || "<<nzBlockCnt<<endl;
#endif

      //totalBlocks+=nzBlockCnt;
      if(col==lc){
        blkCount[I-1] = nzBlockCnt;
        totalBlocks+=nzBlockCnt;
      }
    }
  }
  blkCount.back() = totalBlocks;
}


void symbolBuildRowtab(SymbolMatrix *symbptr) {
  SymbolCblk *cblk;
  SymbolBlok *blok;
  int *innbr, *intmp, *browtab;
  int  itercblk;
  int  cblknbr;
  int  edgenbr = symbptr->bloknbr - symbptr->cblknbr;

  cblknbr = symbptr->cblknbr;

  innbr = new int[cblknbr];
  memset( innbr, 0, cblknbr * sizeof(int) );

  /* Count the number of input edge per cblk */
  cblk = symbptr->cblktab;
  blok = symbptr->bloktab;
  for(itercblk=0; itercblk<cblknbr; itercblk++, cblk++)
  {
    int iterblok = cblk[0].bloknum + 1;
    int lbloknum = cblk[1].bloknum;

    /* Skip diagonal block */
    blok++;

    /* Off-diagonal blocks */
    for( ; iterblok < lbloknum; iterblok++, blok++)
    {
      innbr[ blok->fcblknm ]++;
    }
  }

  /* Initialize the brownum fields */
  cblk = symbptr->cblktab;
  cblk->brownum = 0;
  for(itercblk=0; itercblk<cblknbr; itercblk++, cblk++)
  {
    cblk[1].brownum = cblk[0].brownum + innbr[ itercblk ];
    innbr[itercblk] = cblk[0].brownum;
  }
  assert( cblk[0].brownum == edgenbr );

  /* Initialize the browtab */
  browtab = new int[edgenbr];

  cblk = symbptr->cblktab;
  blok = symbptr->bloktab;
  for(itercblk=0; itercblk<symbptr->cblknbr; itercblk++, cblk++)
  {
    int iterblok = cblk[0].bloknum + 1;
    int lbloknum = cblk[1].bloknum;

    /* Skip diagonal block */
    blok++;

    /* Off-diagonal blocks */
    for( ; iterblok < lbloknum; iterblok++, blok++)
    {
      intmp = innbr + blok->fcblknm;
      browtab[ *intmp ] = iterblok;
      (*intmp)++;
    }
  }
  //if (symbptr->browtab == NULL) {
  //   delete [] symbptr->browtab;
  //}
  symbptr->browtab = browtab;

  delete [] innbr;
  return;
}



Order * GetPastixOrder(SymbolMatrix * symbptr,const vector<int> & xsuper, const LIBCHOLESKY::ETree & etree, const int * perm, const int * iperm){

  Order * order = new Order();
  order->baseval = 1;
  order->vertnbr = symbptr->nodenbr;
  order->cblknbr = symbptr->cblknbr;

  order->permtab = new int[order->vertnbr];
  for(int i =0; i<order->vertnbr;++i){order->permtab[i] = iperm[i];}

  order->peritab = new int[order->vertnbr];
  for(int i =0; i<order->vertnbr;++i){order->peritab[i] = perm[i];}

  order->rangtab = new int[xsuper.size()];
  for(int i = 0;i<xsuper.size();++i){order->rangtab[i] = xsuper[i];} 

  order->treetab = new int[etree.n()];
  for(int i = 0;i<xsuper.size();++i){order->rangtab[i] = etree.PostParent(i);} 




  return order;
};


SymbolMatrix *  GetPastixSymbolMatrix(const vector<int> & xsuper,const vector<int> & supMembership, vector<int64_t> & xlindx, vector<int32_t> & lindx){
  //count the number of blocks and blocks per supernode
  vector<int> blkCount;
  countBlock(xsuper, xlindx, lindx, supMembership, blkCount);


  SymbolMatrix * symbmtx = new SymbolMatrix;

  symbmtx->baseval = 1;
  //ignore dof
  symbmtx->cblknbr = xsuper.size()-1;
  symbmtx->bloknbr = blkCount.back();
  //nodenbr should be n
  symbmtx->nodenbr = xsuper.back()-1;





  symbmtx->cblktab = new SymbolCblk[symbmtx->cblknbr+1];
  symbmtx->bloktab = new SymbolBlok[symbmtx->bloknbr];
  int blockIdx = 1;
  for(int I = 1; I<xlindx.size();++I){
    //count number of contiguous blocks (at least one diagonal block)
    int32_t fc = xsuper[I-1];
    int32_t lc = xsuper[I]-1;
    int32_t fi = xlindx[I-1];
    int32_t li = xlindx[I]-1;
    int32_t iPrevRow = lindx[fi-1]-1;
    int32_t iFirstRow = lindx[fi-1];

    int32_t width = lc - fc + 1; 

    //only go through last column
    for(int col = lc; col<=lc;col++){ 
      //1 to count the diagonal block, 0 to skip it
      int32_t nzBlockCnt = 1;
      int prevFacing = I;
      for(int64_t idx = fi; idx<=li;idx++){
        int32_t iRow = lindx[idx-1];
        //enforce the first block to be a square diagonal block
        int facing = supMembership[iRow-1];
        if(nzBlockCnt==1 && iRow>col){

          symbmtx->bloktab[blockIdx-1].frownum = fc;  /*< First row index            */
          symbmtx->bloktab[blockIdx-1].lrownum = lc;  /*< Last row index (inclusive) */
          symbmtx->bloktab[blockIdx-1].lcblknm = I;  /*< Local column block         */
          symbmtx->bloktab[blockIdx-1].fcblknm = I;  /*< Facing column block        */ //>> CBLK que tu update (ou cblk couran pr le bloc diagonal)


          symbmtx->cblktab[I-1].fcolnum = fc;
          symbmtx->cblktab[I-1].lcolnum = lc;
          symbmtx->cblktab[I-1].bloknum = blockIdx;  /*< First block in column (diagonal) */


          blockIdx++;
                        nzBlockCnt++;

          iFirstRow=iRow;
        }
        else if(iRow!=iPrevRow+1 || prevFacing != facing){

if(prevFacing != facing){
  std::cout<<"BLOCK SPLITTED "<<prevFacing<<" "<<facing<<std::endl;
}

          int facingSnode = supMembership[iFirstRow-1];
  std::cout<<"facing "<<facingSnode<<" "<<facing<<std::endl;
          symbmtx->bloktab[blockIdx-1].frownum = iFirstRow;  /*< First row index            */
          symbmtx->bloktab[blockIdx-1].lrownum = iPrevRow;  /*< Last row index (inclusive) */
          symbmtx->bloktab[blockIdx-1].lcblknm = I;  /*< Local column block         */
          symbmtx->bloktab[blockIdx-1].fcblknm = facingSnode;  /*< Facing column block        */ //>> CBLK que tu update (ou cblk couran pr le bloc diagonal)


          blockIdx++;
                        nzBlockCnt++;

          iFirstRow=iRow;
        }

        iPrevRow=iRow;
        prevFacing = facing;
      }

      if(col==lc){
          int facingSnode = supMembership[iFirstRow-1];
          symbmtx->bloktab[blockIdx-1].frownum = iFirstRow;  /*< First row index            */
          symbmtx->bloktab[blockIdx-1].lrownum = iPrevRow;  /*< Last row index (inclusive) */
          symbmtx->bloktab[blockIdx-1].lcblknm = I;  /*< Local column block         */
          symbmtx->bloktab[blockIdx-1].fcblknm = facingSnode;  /*< Facing column block        */ //>> CBLK que tu update (ou cblk couran pr le bloc diagonal)


          blockIdx++;
      }

    }
  }


  symbolBuildRowtab(symbmtx);




  return symbmtx;
}




///**
// *
// * @file symbol_reordering.c
// *
// *  PaStiX symbol structure routines
// *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
// *  LaBRI, University of Bordeaux 1 and IPB.
// *
// * @version 5.1.0
// * @author Gregoire Pichon
// * @author Mathieu Faverge
// * @author Pierre Ramet
// * @date 2015-04
// *
// **/
//#include "common.h"
//#include "symbol.h"
//#include "order.h"

  static inline int
compute_cblklevel( const int *treetab,
    int *levels,
    int  cblknum )
{
  /* cblknum level has already been computed */
  if ( levels[cblknum] != 0 ) {
    return levels[cblknum];
  }
  else {
    int father = treetab[cblknum];

    if ( father == -1 ) {
      return 1;
    }
    else {
      return compute_cblklevel( treetab, levels, father ) + 1;
    }
  }
}

  static inline int
hamming_distance(int **vectors,
    int  *vectors_size,
    int   xi,
    int   xj,
    int   stop)
{
  /* For the fictive vertex */
  if (xi == -1){
    return vectors_size[xj];
  }
  if (xj == -1){
    return vectors_size[xi];
  }

  int sum = 0;
  int *set1 = vectors[xi];
  int *set2 = vectors[xj];
  int *end1 = vectors[xi] + vectors_size[xi];
  int *end2 = vectors[xj] + vectors_size[xj];

  if (vectors_size[xi] - vectors_size[xj] >= stop){
    return stop;
  }
  if (vectors_size[xj] - vectors_size[xi] >= stop){
    return stop;
  }

  while((set1 < end1) && (set2 < end2))
  {
    if( *set1 == *set2)
    {
      set1++;
      set2++;
    }
    else if( *set1 < *set2 )
    {
      while (( set1 < end1 ) && ( *set1 < *set2 ))
      {
        sum ++;
        set1++;
      }
    }
    else if( *set1 > *set2 )
    {
      while (( set2 < end2 ) && ( *set1 > *set2 ))
      {
        sum ++;
        set2++;
      }
    }
    else
    {
      assert(0);
    }

    /* The computation is stopped if sum overlapped a given limit */
    if (sum >= stop){
      return stop;
    }
  }

  sum += end1-set1;
  sum += end2-set2;

  if (sum >= stop){
    return stop;
  }

  return sum;
}


  static inline void
symbol_reorder_tsp(int size, Order *order, int sn_id,
    int **lw_vectors, int *lw_vectors_size,
    int **up_vectors, int *up_vectors_size,
    int stop_criteria, int stop_when_fitting)
{

  if ( size < 3 ) {
    return;
  }

  int  i, j, k, l, elected;
  int *tmpinvp;
  int *tmplen;
  int distance;

  tmpinvp = (int*)malloc((size+1)*sizeof(int));
    tmplen = (int*)malloc((size+1)*sizeof(int));
    memset(tmplen, 0, (size+1)*sizeof(int));

  tmpinvp[0] = -1;
  tmpinvp[1] = 0;

  distance = hamming_distance(lw_vectors, lw_vectors_size, 0, -1, stop_criteria);

  tmplen[0] = distance;
  tmplen[1] = distance;

  int min_cut = -1;
  for(i=1; i<size; i++) {
    int first_pos;
    int last_pos;

    int lw_before_pos;
    int lw_after_pos;

    int up_before_pos;
    int up_after_pos;

    /* Start by adding the row in first position */
    lw_before_pos = hamming_distance(lw_vectors, lw_vectors_size, i, tmpinvp[0], stop_criteria);
    lw_after_pos  = hamming_distance(lw_vectors, lw_vectors_size, i, tmpinvp[1], stop_criteria);
    up_after_pos  = hamming_distance(up_vectors, up_vectors_size, i, tmpinvp[1], 1);

    int minl = lw_before_pos + lw_after_pos - tmplen[0];
    int mpos = 1;
    int min_cut = -1;

    for(j=1; j<i; j++ ){
      up_before_pos = up_after_pos;
      up_after_pos  = hamming_distance(up_vectors, up_vectors_size, i, tmpinvp[j+1], 1);

      if ( up_before_pos < 1 ||
          up_after_pos  < 1 )
      {

        /* If split was used previously, this first distance may not be already computed */
        if (lw_after_pos == -1)
          lw_before_pos = hamming_distance(lw_vectors, lw_vectors_size, i, tmpinvp[j], stop_criteria);
        else
          lw_before_pos = lw_after_pos;


        lw_after_pos = hamming_distance(lw_vectors, lw_vectors_size, i, tmpinvp[j+1], stop_criteria);

        l = lw_before_pos + lw_after_pos - tmplen[j];


        /* Minimize the cut between two lines, for the same TSP result */
        if ( l == minl ) {
          if (lw_before_pos < min_cut){
            min_cut = lw_before_pos;
            minl = l; mpos = j+1;
          }
          if (lw_after_pos < min_cut){
            min_cut = lw_after_pos;
            minl = l; mpos = j+1;
          }
        }

        /* Position that minimizes TSP */
        if ( l < minl ) {
          minl = l; mpos = j+1;

          min_cut = lw_before_pos;
          if (lw_after_pos < min_cut){
            min_cut = lw_after_pos;
          }
        }

        if ( l < minl ) {
          minl = l; mpos = j+1;
          min_cut = lw_before_pos;
          if (lw_after_pos < min_cut){
            min_cut = lw_after_pos;
          }
        }


        /* Stop if two lines are equal (already done tmpinvp[j]) */
        if (lw_after_pos == 0){
          min_cut = 0;
          minl = l; mpos = j+1;
          j = i;
        }
      }
      else{
        lw_after_pos = -1;
      }

    }

    /* Test between last and first */
    first_pos = hamming_distance(lw_vectors, lw_vectors_size, i, tmpinvp[0], stop_criteria);
    last_pos  = hamming_distance(lw_vectors, lw_vectors_size, i, tmpinvp[i], stop_criteria);

    lw_before_pos = hamming_distance(lw_vectors, lw_vectors_size, i, tmpinvp[mpos-1], stop_criteria);
    lw_after_pos  = hamming_distance(lw_vectors, lw_vectors_size, i, tmpinvp[mpos  ], stop_criteria);

    l = first_pos + last_pos - tmplen[i];
    if ( l < minl ) {
      minl = l; mpos = i+1;
    }

    if (mpos > 0){
      tmplen[mpos-1] = lw_before_pos;
    }

    if (mpos < (i+1))
    {
      int tmpi, tmpl;
      k = i;
      l = lw_after_pos;

      /* Insert the line in the tmpinvp/tmplen arrays */
      for(j=mpos; j<i+2; j++ )
      {
        tmpi = tmpinvp[j];
        tmpl = tmplen[j];

        tmpinvp[j] = k;
        tmplen[j]  = l;

        k = tmpi;
        l = tmpl;
      }
    }
    else {
      tmpinvp[i+1] = i;
      tmplen[i+1]  = first_pos;
    }
  }

  elected = 0;
  for (i=0; i<size; i++)
  {
    if (tmpinvp[i] == -1){
      elected = i;
    }
  }

  int *sn_connected;
  sn_connected = (int*)malloc(size*sizeof(int));
  {
    //TODO
    int *peritab = order->peritab + order->rangtab[sn_id];
    for (i=0; i<size; i++)
    {
      sn_connected[i] = peritab[ tmpinvp[(i + 1 + elected)%(size+1)] ];
    }
    memcpy( peritab, sn_connected, size * sizeof(int) );
  }

  free(sn_connected);
  free(tmpinvp);
  free(tmplen);
}

  static inline void
symbol_reorder_cblk( const SymbolMatrix *symbptr,
    const SymbolCblk   *cblk,
    Order              *order,
    const int *levels,
    int        cblklvl,
    int       *depthweight,
    int        depthmax,
    int        split_level,
    int                 stop_criteria,
    int                 stop_when_fitting,
    double             *time_compute_vectors,
    double             *time_update_perm)
{
  SymbolBlok *blok;
  int **up_vectors, *up_vectors_size;
  int **lw_vectors, *lw_vectors_size;
  int size = cblk->lcolnum - cblk->fcolnum + 1;
  int local_split_level = split_level;
  int i, iterblok;
  int *brow = symbptr->browtab;
  double timer;

  /**
   * Compute hamming vectors in two subsets:
   *   - The upper subset contains the cblk with level higher than the split_level
   *     in the elimination tree, (or depth lower than levels[cblk])
   *   - The lower subset contains the cblk with level lower than the split_level
   *     in the elimination tree, (or depth higher than levels[cblk])
   *
   * The delimitation between the lower and upper levels is made such that
   * the upper level represents 17% to 25% of the total number of cblk.
   */
//  clockStart(timer);
  {
    int weight = 0;

    /* Compute the weigth of each level */
    //MATHIAS: this is a loop through blocks in current column 
    for(iterblok=cblk[0].brownum; iterblok<cblk[1].brownum; iterblok++)
    {
      int blokweight;
      blok = symbptr->bloktab + brow[iterblok];
      blokweight = blok->lrownum - blok->frownum + 1;
      depthweight[ levels[ blok->lcblknm ] - 1 ] += blokweight;
      weight += blokweight;
    }

    /**
     * Compute the split_level:
     *    We start with the given split_level parameter
     *    and we try to correct it to minimize the following iterative process
     */
    {
      /* Current for each line within the current cblk the number of contributions */
      int up_total = 0;
      int lw_total = 0;
      int sign = 0;

split:
      up_total = 0;
      lw_total = 0;

      for(i=0; i<local_split_level; i++)
      {
        up_total += depthweight[i];
      }
      for(; i<depthmax; i++)
      {
        lw_total += depthweight[i];
      }

      /* If there are too many upper bloks */
      if ( (lw_total < (5 * up_total)) &&
          (lw_total > 10) && (up_total > 10) && (sign <= 0))
      {
        local_split_level--;
        sign--;
        goto split;
      }

      /* If there are too many lower bloks */
      if ( (lw_total > (3 * up_total)) &&
          (lw_total > 10) && (up_total > 10) && (sign >= 0) )
      {
        local_split_level++;
        sign++;
        goto split;
      }

      /* Adjust to depth of the level array */
      /* symbol_reorder_cblk( symbptr, cblk, order, */
      /*                      levels, levels[itercblk], */
      /*                      depthweight + levels[itercblk], maxdepth-levels[itercblk], */
      /*                      split_level, stop_criteria, stop_when_fitting, */
      /*                      &time_compute_vectors, &time_update_perm); */
      /* local_split_level += cblklvl; */
      /* for(i=0; (i<local_split_level) && (depthweight[i] != 0); i++) */
      /* for(; (i<depthmax) && (depthweight[i] != 0); i++) */
    }

    /* Compute the Hamming vector size for each row of the cblk */
    up_vectors_size = (int*)malloc(size*sizeof(int));
    memset(up_vectors_size, 0, size * sizeof(int));
    lw_vectors_size = (int*)malloc(size*sizeof(int));
    memset(lw_vectors_size, 0, size * sizeof(int));

    for(iterblok=cblk[0].brownum; iterblok<cblk[1].brownum; iterblok++)
    {
      blok = symbptr->bloktab + brow[iterblok];

      /* For upper levels in nested dissection */
      if (levels[blok->lcblknm] <= local_split_level){
        for (i=blok->frownum; i<=blok->lrownum; i++){
          int index = i - cblk->fcolnum;
          up_vectors_size[index]++;
        }
      }
      else{
        for (i=blok->frownum; i<=blok->lrownum; i++){
          int index = i - cblk->fcolnum;
          lw_vectors_size[index]++;
        }
      }
    }

    /* Initiate Hamming vectors structure */
    lw_vectors = (int**)malloc(size*sizeof(int*));
    up_vectors = (int**)malloc(size*sizeof(int*));
    for (i=0; i<size; i++) {
      lw_vectors[i] = (int*)malloc(lw_vectors_size[i]*sizeof(int));
      up_vectors[i]=(int*)malloc(up_vectors_size[i]*sizeof(int));
      memset(lw_vectors[i], 0, lw_vectors_size[i] * sizeof(int));
      memset(up_vectors[i], 0, up_vectors_size[i] * sizeof(int));
    }
    memset(lw_vectors_size, 0, size * sizeof(int));
    memset(up_vectors_size, 0, size * sizeof(int));

    /* Fill-in vectors structure with contributing cblks */
    for(iterblok=cblk[0].brownum; iterblok<cblk[1].brownum; iterblok++)
    {
      blok = symbptr->bloktab + brow[iterblok];

      /* For upper levels in nested dissection */
      if (levels[blok->lcblknm] <= local_split_level) {
        for (i=blok->frownum; i<=blok->lrownum; i++){
          int index = i - cblk->fcolnum;
          up_vectors[index][up_vectors_size[index]] = blok->lcblknm;
          up_vectors_size[index]++;
        }
      }
      else{
        for (i=blok->frownum; i<=blok->lrownum; i++){
          int index = i - cblk->fcolnum;
          lw_vectors[index][lw_vectors_size[index]] = blok->lcblknm;
          lw_vectors_size[index]++;
        }
      }
    }
  }

//  clockStop(timer);
//  *time_compute_vectors += clockVal(timer);

//  clockStart(timer);
  {
    /* Apply the pseudo-TSP algorithm to the rows in the current supernode */
    symbol_reorder_tsp(size, order, cblk - symbptr->cblktab,
        lw_vectors, lw_vectors_size,
        up_vectors, up_vectors_size,
        stop_criteria, stop_when_fitting);
  }
//  clockStop(timer);
//  *time_update_perm += clockVal(timer);

  for (i=0; i<size; i++){
    free(lw_vectors[i]);
    free(up_vectors[i]);
  }

  free(lw_vectors);
  free(up_vectors);
  free(lw_vectors_size);
  free(up_vectors_size);
}

/* For split_level parameter */
/* The chosen level to reduce computational cost: no effects if set to 0 */
/* A first comparison is computed according to upper levels */
/* If hamming distances are equal, the computation goes through lower levels */

/* For stop_criteria parameter */
/* Criteria to limit the number of comparisons when computing hamming distances */

/* For stop_when_fitting parameter */
/* Criteria to insert a line when no extra-blok is created */
/* If set to 0, the algorithm will minimize the cut between two lines */

  void
symbolReordering( const SymbolMatrix *symbptr,
    Order *order,
    int split_level,
    int stop_criteria,
    int stop_when_fitting )
{
  SymbolCblk  *cblk;
  int itercblk;
  int cblknbr = symbptr->cblknbr;

  double time_compute_vectors = 0.;
  double time_update_perm     = 0.;

  int i, maxdepth;
  int *levels, *depthweight;

  /* Create the level array to compute the depth of each cblk and the maximum depth */
  {
    maxdepth = 0;
    levels = new int[cblknbr];
    for(int i = 0; i<cblknbr; ++i){levels[i] = 0;}

    for (i=0; i<cblknbr; i++) {
      levels[i] = compute_cblklevel( order->treetab, levels, i );
      maxdepth = std::max( maxdepth, levels[i] );
    }
  }

  /**
   * Solves the Traveler Salesman Problem on each cblk to minimize the number
   * of off-diagonal blocks per row
   */
  cblk = symbptr->cblktab;
  depthweight = new int[maxdepth];
  for (itercblk=0; itercblk<cblknbr; itercblk++, cblk++) {

    memset( depthweight, 0, maxdepth * sizeof(int) );

    symbol_reorder_cblk( symbptr, cblk, order,
        levels, levels[itercblk],
        depthweight, maxdepth,
        split_level, stop_criteria, stop_when_fitting,
        &time_compute_vectors, &time_update_perm);
  }

  printf("Time to compute vectors  %lf s\n", time_compute_vectors);
  printf("Time to update  perm     %lf s\n", time_update_perm);

  /* Update the permutation */
  for (i=0; i<symbptr->nodenbr; i++) {
    order->permtab[ order->peritab[i] ] = i;
  }
  delete [] levels; levels = NULL;
  delete [] depthweight; depthweight = NULL;
}

  void
symbolReorderingPrintComplexity( const SymbolMatrix *symbptr )
{
  SymbolCblk  *cblk;
  int itercblk, iterblok;
  int cblknbr;
  int nbflops;

  cblk    = symbptr->cblktab;
  cblknbr = symbptr->cblknbr;
  nbflops = 0;

  /**
   * nbcblk is the number of non zeroes intersection between indivudal rows
   * and block columns.
   */
  for(itercblk=0; itercblk<cblknbr; itercblk++, cblk++)
  {
    int width;
    int nbcblk = 0;

    for(iterblok=cblk[0].brownum; iterblok<cblk[1].brownum; iterblok++)
    {
      SymbolBlok *blok = symbptr->bloktab + symbptr->browtab[iterblok];
      assert( blok->fcblknm == itercblk );

      nbcblk += blok->lrownum - blok->frownum + 1;
    }
    width = cblk->lcolnum - cblk->fcolnum + 1;
    nbflops += nbcblk * (width-1);

    if ( itercblk == (cblknbr-1) ) {
      int localflops = nbcblk * (width-1);
      fprintf(stdout, " Number of operations in reordering for last supernode: %ld (%lf)\n",
          localflops, (double)localflops / (double)(nbflops) * 100. );
    }
  }
  fprintf(stdout, " Number of operations in reordering: %ld\n", nbflops );
}
