# TODO
----------
This document contains a list of known issues with symPACK, as well as a list of potential directions for future enhancements to symPACK.

## Known Issues
----------
#### Reading Matrix Market (.mm) files
While symPACK does have a method that is meant to read Matrix Market files, this method does not work, leaving symPACK only capable of reading Rutherford-Boeing (.rb) files.

## Future Work
----------
#### Supernode/block amalgamation 
In certain cases, it may be beneficial to combine contiguous blocks or supernodes into a single entity. This could reduce communication costs and enable more effective use of BLAS routines on each block. 
For more information on supernode amalgamation, see **C. C. Ashcraft and R. G. Grimes.** "The inuence of relaxed supernode partitions on the multifrontal method." ACM Trans. Math. Software, 15(4):291{309, Dec. 1989.

#### Improved GPU offloading heuristic 
The default size thresholds for offloading computations to the GPU are not extensively optimized, and are likely only suitable for certain devices. It is not ideal for users to have to manually tune these parameters, so some sort of more intelligent heuristic to determine
when a computation should be offloaded to a GPU could be beneficial.

#### Add GPU operations to LDL^t decomposition
In addition to standard Cholesky factorization (LL^t), symPACK also supports LDL^t decomposition. This decomposition option does not currently support GPU operations.
It would be relatively easy to implement GPU operations for the LDL^t decomposition, since symPACK's code for this decomposition is largely analogous to symPACK's standard Cholesky factorization code, with some minor differences.

#### Improved intra-node scheduling
Right now, if multiple tasks are available for execution on the ready task queue (RTQ), symPACK simply selects whichever happens to be first on the queue. 
A more sophisticated scheduling policy could be of potential benefit to symPACK.
