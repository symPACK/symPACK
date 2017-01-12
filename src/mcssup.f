C***********************************************************************
C***********************************************************************
C
C   Version:        0.4
C   Last modified:  February 13, 1995
C   Authors:        Esmond G. Ng and Barry W. Peyton
C
C   Mathematical Sciences Section, Oak Ridge National Laboratory
C
C***********************************************************************
C***********************************************************************
C*************     mcssup ..... MCS Supernode ordering    **************
C***********************************************************************
C***********************************************************************
C
C   PURPOSE: 
C       this routine reorders the supernodes backwards by maximum
c       cardinality.  It is analogous to a maximum cardinality ordering.
c       The ordering can then be used to enhance the performance of 
c       ordsup.f, which attempts to reduce the number of segments
c       (within a supernode) for index lists in LINDX.
C
C   INPUT PARAMETERS:
C       (I) NEQNS       -   NUMBER OF EQUATIONS
C       (I) NOFSUB      -   NUMBER OF SUBSCRIPTS STORED IN LINDX(*).
c                           it is the size of lindx(*)
C       (I) NSUPER      -   NUMBER OF SUPERNODES.
C       (I) XSUPER(*)   -   ARRAY OF LENGTH NSUPER+1, CONTAINING THE
C                           FIRST COLUMN OF EACH SUPERNODE.
C       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING
C                           SUPERNODE MEMBERSHIP.
C       (I) XLINDX      -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS 
C                           INTO THE SUBSCRIPT VECTOR.
C       (I) LINDX       -   ARRAY OF LENGTH nofsub, CONTAINING THE
C                           COMPRESSED SUBSCRIPTS.
c
C   OUTPUT PARAMETERS:
c       (i) supperm(*)  -   array of length neqns, contains the new 
c                           ordering of the supernodes, maps new numbers 
c                           to old.
C       
C   WORKING PARAMETERS:
c       (i) suppar      _   array of length nsuper, contains the parent
c                           vector of the supernodal elimination tree.
c       (i) fchild      -   array of length nsuper, contains the first
c                           child vector for the supernodal elimination
c                           tree.
c       (i) siblng      -   array of length nsuper, contains the next
c                           sibling vector for the supernodal elimination
c                           tree.
c       (i) heap        -   array of length 2*nsuper, contains the heap
c                           for ordering the "active" supernodes by 
c                           maximum cardinality.
c       (i) heapinv     -   array of length nsuper, contains the heap
c                           inverse pointers from the supernode to its
c                           cell in the heap.
c       (i) ndesc       -   array of length nsuper, contains the number
c                           of descendants in the supernodal elimination
c                           tree (counting self as a descendant).
C
C***********************************************************************
C
      SUBROUTINE  mcssup (  ordflag, neqns , nofsub, nsuper, xsuper, 
     &                      snode , xlindx, lindx , supperm, suppar, 
     &                      fchild, siblng, heap  , heapinv, ndesc   )
C
C***********************************************************************
C
C       -----------
C       PARAMETERS.
C       -----------
        integer             neqns , nofsub, nsuper, ordflag
c
        INTEGER             fchild(nsuper), heap(2*nsuper),
     &                      heapinv(nsuper), lindx(nofsub),
     &                      ndesc(nsuper) , siblng(nsuper), 
     &                      snode(neqns)  , suppar(nsuper), 
     &                      supperm(nsuper), xlindx(nsuper+1), 
     &                      xsuper(nsuper+1)
c
c       ----------------
c       local variables.
c       ----------------
        INTEGER             card  , heapsize, istart, istop, jsup  ,
     &                      ksuper, nxtchild, offset, score
C
C***********************************************************************
C
c       -----------------------------------------------------------------
c       initialize null parent vector, and
c       account for jsup's contribution to its own number of descendants.
c       -----------------------------------------------------------------
        do  jsup = 1, nsuper
            suppar(jsup) = 0
            ndesc(jsup) = 1
        end do
c
c       ----------------------------
c       for every supernode jsup ...
c       ----------------------------
        do  jsup = 1, nsuper
            offset = xsuper(jsup+1) - xsuper(jsup)
            istart = xlindx(jsup) + offset
            istop = xlindx(jsup+1) - 1
c           -----------------------------------------------------------
c           if jsup is not a root in the supernode elimination tree ...
c           -----------------------------------------------------------
            if  ( istart .le. istop )  then
c               --------------------------------------------------
c               store parent supernode of jsup.
c               compute jsup's contribution to its parent's number 
c               of descendants.
c               --------------------------------------------------
                suppar(jsup) = snode(lindx(istart))
                ndesc(suppar(jsup)) = ndesc(suppar(jsup)) + ndesc(jsup)
            end if
        end do
c
c       ----------------------------------------------------------
c       construct the binary tree representation of the supernodal
c       elimination tree.
c       ----------------------------------------------------------
        call  betree ( nsuper, suppar, fchild, siblng )
c
c       ----------------------------
c       initially the heap is empty.
c       ----------------------------
        heapsize  = 0
c       ---------------------------
c       for each supernode jsup ...
c       ---------------------------
        do  jsup = 1, nsuper
c           --------------------------------------
c           if jsup is the root of its subtree ...
c           --------------------------------------
            if  ( suppar(jsup) .eq. 0 )  then
c               -----------------------------------------------
c               insert jsup onto the heap with cardinality zero
c               as its score.
c               -----------------------------------------------
                call  ins_heap ( heap, heapsize, heapinv, jsup, 
     &                           0                              )
            end if
        end do
c
c       ---------------------------------------------------
c       do while there are unnumbered supernodes remaining.
c       ---------------------------------------------------
        ksuper = nsuper
        do while  ( ksuper .gt. 0 ) 
c           ----------------------------------------------------------
c           get a supernode jsup of maximum cardinality from the heap,
c           and delete it from the heap.
c           ----------------------------------------------------------
            jsup = heap(2)
            call  del_heap ( heap, heapsize, heapinv, jsup )
c           -----------------------------------
c           for each child nxtchild of jsup ...
c           -----------------------------------
            nxtchild = fchild(jsup)
            do while  ( nxtchild .gt. 0 )  
c               ------------------------------------------------------
c               use cardinality or number of descendants as the score.
c               ------------------------------------------------------
                if  ( ordflag .eq. 2 )  then
                    offset = xsuper(nxtchild+1) - xsuper(nxtchild)
                    card = xlindx(nxtchild+1) - xlindx(nxtchild) 
     &                     - offset
                    score  = - card
                else if  ( ordflag .eq. 3 )  then
                    score = - ndesc(nxtchild)
                end if
c               ----------------------------------------------
c               insert nxtchild onto the heap using its score.
c               ----------------------------------------------
                call  ins_heap  ( heap, heapsize, heapinv, nxtchild,
     &                            score                               )
                nxtchild = siblng(nxtchild)
            end do
c           ------------------------------------------------------
c           number supernode jsup next (backward) in the ordering.
c           ------------------------------------------------------
            supperm(ksuper) = jsup
            ksuper = ksuper - 1
        end do
        RETURN
      END
