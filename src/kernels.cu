#include "sympack/kernels.cuh"
#include <stdio.h>

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
		   	   double * offset_arr) {
	for (int cr = row; cr<lr; cr++) {
	    offset+=tgt_width;
	    offset_arr[rowidx] = tgtOffset + (cr-row)*tgt_snode_size;
	    rowidx++;
	}
}	

void set_offset_wrapper(int lr, int tgtOffset, int offset, int tgt_width,
	       		   int row, int rowidx, int tgt_snode_size,
		   	   double * offset_arr) {
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
