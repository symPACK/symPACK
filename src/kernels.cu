#include "sympack/kernels.cuh"
#include <stdio.h>

/* These are cuda kernels that could theoretically replace some of the 
 * array operations that happen in blockCell_t::update.
 * Right now they are not used for anything, but I've decided to leave them here
 * in case they one day might be useful.
 */

namespace symPACK {
namespace cudaKernels {
__global__ void update_tgt(int src_nrows,
                            int first_col, int tgt_width,
                            int * col_array, int col_array_sz,
                            int * offset_array, 
                            double * tgt, double * buf) {
    //TODO: Parallelize this
    //int n_thread_rows = src_nrows / blockDim.x;
    //int start_row = (n_thread_rows * threadIdx.x);
    for (int rowidx = 0; rowidx < src_nrows; rowidx++) {
        for (int colidx = 0; colidx < col_array_sz; colidx++) {
            int col = col_array[colidx];
            int tgt_colidx = col - first_col;
            tgt[offset_array[rowidx] + tgt_colidx] += buf[rowidx*tgt_width+colidx];
        }
    }
}


__global__ void set_offset(int lr, int tgtOffset, int offset, int tgt_width,
	       		   int row, int rowidx, int tgt_snode_size,
		   	   int * offset_arr) {
	for (int cr = row; cr<=lr; cr++) {
	    offset+=tgt_width;
	    offset_arr[rowidx] = tgtOffset + (cr-row)*tgt_snode_size;
	    rowidx++;
	}
}	


__global__ void set_colindx(int colindx, int cur_src_fc, int cur_src_lc, int * colindx_arr) {
	
	for (int col = cur_src_fc; col <= cur_src_lc; col++) {
	     colindx_arr[colindx++] = col;
	}

}


void set_colindx_wrapper(int colindx, int cur_src_fc, int cur_src_lc, int * colindx_arr) {
	set_colindx<<<1, 1>>>(colindx, cur_src_fc, cur_src_lc, colindx_arr);
}


void set_offset_wrapper(int lr, int tgtOffset, int offset, int tgt_width,
	       		   int row, int rowidx, int tgt_snode_size,
		   	   int * offset_arr) {
     set_offset<<<1, 1>>>(lr, tgtOffset, offset, tgt_width,
		     	  row, rowidx, tgt_snode_size,
			  offset_arr);
}

void update_tgt_wrapper(int src_nrows,
        		int first_col, int tgt_width,
        		int * col_array, int col_array_sz,
        		int * offset_array, 
        		double * tgt, double * buf) {
	update_tgt<<<1, 1>>>(src_nrows, first_col, tgt_width,
               		     col_array, col_array_sz,
                             offset_array, tgt, buf);
}



}//namespace cudaKernels
}//namespace symPACK
