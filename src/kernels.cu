#include "sympack/kernels.cuh"

namespace symPACK {
namespace cudaKernels {
__global__ void update_tgt(int src_nrows,
                            int first_col, int tgt_width,
                            int * col_array, int col_array_sz,
                            int * offset_array, 
                            double * tgt, double * buf) {
    
    int n_thread_rows = src_nrows / blockDim.x;
    int start_row = (n_thread_rows * threadIdx.x);
    for (int rowidx = start_row; rowidx < start_row + n_thread_rows; rowidx++) {
        for (int colidx = 0; colidx < col_array_sz; colidx++) {
            int col = col_array[colidx];
            int tgt_colidx = col - first_col;
            tgt[offset_array[rowidx] + tgt_colidx] += buf[rowidx * tgt_width + colidx];
        }
    }

}

    
    void update_tgt_wrapper(int src_nrows,
        int first_col, int tgt_width,
        int * col_array, int col_array_sz,
        int * offset_array, 
        double * tgt, double * buf) {
            update_tgt<<<1, src_nrows/4>>>(src_nrows, first_col, tgt_width,
                            col_array, col_array_sz,
                            offset_array, tgt, buf);
        }
}
}