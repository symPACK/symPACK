#ifndef __PASTIX_HEADER__
#define __PASTIX_HEADER__

#include "ngchol/ETree.hpp"
#include <vector>

typedef struct SymbolBlok_ {
  int frownum;  /*< First row index            */
  int lrownum;  /*< Last row index (inclusive) */
  int lcblknm;  /*< Local column block         */
  int fcblknm;  /*< Facing column block        */ //>> CBLK que tu update (ou cblk couran pr le bloc diagonal)
} SymbolBlok;

typedef struct SymbolCblk_ {
  int fcolnum;  /*< First column index               */
  int lcolnum;  /*< Last column index (inclusive)    */
  int bloknum;  /*< First block in column (diagonal) */
  int brownum;  /*< First block in row facing the diagonal block in browtab, 0-based */
} SymbolCblk;



typedef struct SymbolMatrix_ {
  int            baseval;  /*< Base value for numberings               */
  int            dof;      /*< Degrees of freedom per node
                             (constant if > 0, unconstant if 0 (not implemented)) */
  int            cblknbr;  /*< Number of column blocks                 */
  int            bloknbr;  /*< Number of blocks                        */
  int            nodenbr;  /*< Number of node in the compressed symbol */
  SymbolCblk   * cblktab;  /*< Array of column blocks [+1,based]       */
  SymbolBlok   * bloktab;  /*< Array of blocks [based]                 */
  int * browtab;  /*< Global index of a given block number in bloktab 0-based    */
} SymbolMatrix;

typedef struct Order_ {
  int  baseval;   /*< base value used for numbering       */
  int  vertnbr;   /*< Number of vertices                  */
  int  cblknbr;   /*< Number of column blocks             */
  int *permtab;   /*< Permutation array [based]           */ //scotch
  int *peritab;   /*< Inverse permutation array [based]   */
  int *rangtab;   /*< Column block range array [based,+1] */
  int *treetab;   /*< Partitioning tree [based]           */
} Order;






void countBlock(const vector<int> & Xsuper, const vector<int64_t> & Xlindx, const vector<int32_t> & Lindx, const vector<int> & supMembership, vector<int> & blkCount);
void symbolBuildRowtab(SymbolMatrix *symbptr);
Order * GetPastixOrder(SymbolMatrix * symbptr,const vector<int> & xsuper, const LIBCHOLESKY::ETree & etree, const int * perm, const int * iperm);
SymbolMatrix *  GetPastixSymbolMatrix(const vector<int> & xsuper,const vector<int> & supMembership, vector<int64_t> & xlindx, vector<int32_t> & lindx);


void symbolReordering( const SymbolMatrix *symbptr, Order *order, int split_level, int stop_criteria, int stop_when_fitting );
void symbolReorderingPrintComplexity( const SymbolMatrix *symbptr );
#endif
