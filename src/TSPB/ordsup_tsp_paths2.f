C***********************************************************************
C***********************************************************************
C
C   Version:        0.4
C   Last modified:  January 12, 1995
C   Authors:        Esmond G. Ng and Barry W. Peyton
C
C   Mathematical Sciences Section, Oak Ridge National Laboratory
C
C***********************************************************************
C***********************************************************************
C**************     FCNTHN  ..... FIND NONZERO COUNTS    ***************
C***********************************************************************
C***********************************************************************
C
C   PURPOSE:
C       THIS SUBROUTINE DETERMINES THE ROW COUNTS AND COLUMN COUNTS IN
C       THE CHOLESKY FACTOR.  IT USES A DISJOINT SET UNION ALGORITHM.
C
C       TECHNIQUES:
C       1) SUPERNODE DETECTION.
C       2) PATH HALVING.
C       3) NO UNION BY RANK.
C
C   NOTES:
C       1) ASSUMES A POSTORDERING OF THE ELIMINATION TREE.
C
C   INPUT PARAMETERS:
C       (I) NEQNS       -   NUMBER OF EQUATIONS.
C       (I) ADJLEN      -   LENGTH OF ADJACENCY STRUCTURE.
C       (I) XADJ(*)     -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS
C                           TO THE ADJACENCY STRUCTURE.
C       (I) ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1, CONTAINING
C                           THE ADJACENCY STRUCTURE.
C       (I) PERM(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE
C                           POSTORDERING.
C       (I) INVP(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE
C                           INVERSE OF THE POSTORDERING.
C       (I) ETPAR(*)    -   ARRAY OF LENGTH NEQNS, CONTAINING THE
C                           ELIMINATION TREE OF THE POSTORDERED MATRIX.
C
C   OUTPUT PARAMETERS:
C       (I) ROWCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER
C                           OF NONZEROS IN EACH ROW OF THE FACTOR,
C                           INCLUDING THE DIAGONAL ENTRY.
C       (I) COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER
C                           OF NONZEROS IN EACH COLUMN OF THE FACTOR,
C                           INCLUDING THE DIAGONAL ENTRY.
C       (I) NLNZ        -   NUMBER OF NONZEROS IN THE FACTOR, INCLUDING
C                           THE DIAGONAL ENTRIES.
C
C   WORK PARAMETERS:
C       (I) SET(*)      -   ARRAY OF LENGTH NEQNS USED TO MAINTAIN THE
C                           DISJOINT SETS (I.E., SUBTREES).
C       (I) PRVLF(*)    -   ARRAY OF LENGTH NEQNS USED TO RECORD THE
C                           PREVIOUS LEAF OF EACH ROW SUBTREE.
C       (I) LEVEL(*)    -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE LEVEL
C                           (DISTANCE FROM THE ROOT).
C       (I) WEIGHT(*)   -   ARRAY OF LENGTH NEQNS+1 CONTAINING WEIGHTS
C                           USED TO COMPUTE COLUMN COUNTS.
C       (I) FDESC(*)    -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE
C                           FIRST (I.E., LOWEST-NUMBERED) DESCENDANT.
C       (I) NCHILD(*)   -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE
C                           NUMBER OF CHILDREN.
C       (I) PRVNBR(*)   -   ARRAY OF LENGTH NEQNS USED TO RECORD THE
C                           PREVIOUS ``LOWER NEIGHBOR'' OF EACH NODE.
C
C   FIRST CREATED ON    APRIL 12, 1990.
C   LAST UPDATED ON     JANUARY 12, 1995.
C
C***********************************************************************
C
        SUBROUTINE ordsup_tsp_paths2
     &                   (  nnodes, nadj  , neqns , xskadj, sklenf, 
     &                      skadj , invp  , perm  , nsuper, xsuper, 
     &                      snode , nsuper2, sperm, suppar, fstloc,
     &                      row_c , curlf , prvlf , prv_l , lback , 
     &                      lforw , rowcnt, fdesc , level , set     )
C
C       -----------
C       PARAMETERS.
C       -----------
c
        INTEGER             nadj  , neqns , nnodes, nsuper, nsuper2
c
        INTEGER             curlf(nnodes)   , fdesc(0:nsuper2),
     &                      fstloc(neqns)   , invp(neqns)     ,
     &                      lback(nnodes)   , level(0:nsuper2),
     &                      lforw(nnodes)   , perm(neqns)     , 
     &                      prv_l(0:nnodes,0:nnodes), prvlf(nnodes), 
     &                      row_c(0:nnodes,0:nnodes), rowcnt(nnodes), 
     &                      set(nsuper2)    , skadj(nadj)     , 
     &                      sklenf(neqns)   , snode(neqns)    , 
     &                      sperm(nsuper2)  , suppar(nsuper2) , 
     &                      xskadj(neqns+1) , xsuper(nsuper+1)
C
C       ----------------
C       LOCAL VARIABLES.
C       ----------------
        INTEGER             floc  , fstnod, i     , i0    , ifdesc, 
     &                      iloc  , inew  , iold  , j     , jleft ,
     &                      jright, ksuper, last1 , last2 , lca   ,
     &                      lca2  , lhead , lloc  , ltail , nleafji,
     &                      oldfst, oldsup, parent, pleafi,
     &                      pleafj, pleafji, rootsup 
c
        logical             skipflag
C
C***********************************************************************
C
c       ***************
c       intializations.
c       ***************

C       ----------------------------
C       COMPUTE LEVEL(*).
C       INITIALIZE SET(*), fdesc(*).
C       ----------------------------
        level(0) = 0
        do  ksuper = nsuper2, 1, -1
            fdesc(ksuper) = ksuper
            set(ksuper) = ksuper
            level(ksuper) = level(suppar(ksuper)) + 1
        end do
c       --------------------------------------
c       rootsup: root supernode --- global ID.
c       i0: last node in previous supernode
c       --------------------------------------
        rootsup = sperm(nsuper2)
        i0 = xsuper(rootsup) - 1
c       ---------------------------------------------------
c       initialize rowcnt(*), prvlf(*), curlf(*).
c       note: indexed by local node ID, not global node ID.
c       ---------------------------------------------------
        do  i = 1, nnodes
            rowcnt(i) = 1
            prvlf(i) = 0
            curlf(i) = 0
        end do
c       ---------------------------------------------------
c       initialize row_c(*,*), prv_l(*,*).
c       note: indexed by local node ID, not global node ID.
c       ---------------------------------------------------
        do  i = 1, nnodes
c           ---------------
c           lower triangle.
c           ---------------
            do  j = 1, i-1
                row_c(i,j) = 1
            end do
c           ---------------
c           upper triangle.
c           ---------------
            do  j = i+1, nnodes
                prv_l(i,j) = 0
            end do
        end do
c       ------------------------------------
c       compute fdesc(*): first descendants.
c       ------------------------------------
        FDESC(0) = 0
        DO  ksuper = 1, nsuper2
            PARENT = suppar(ksuper)
            IFDESC = FDESC(ksuper)
            IF  ( IFDESC .LT. FDESC(PARENT) )  THEN
                FDESC(PARENT) = IFDESC
            ENDIF
        end do   
c       --------------------------------------------------
c       initialize empty "visited" list.
c       initially all nodes are marked as not in the list.
c       --------------------------------------------------
        lhead = 0
        ltail = 0
        do  i = 1, nnodes
            lforw(i) = -1
            lback(i) = -1
        end do
c
c       ************************************************************************
c       main loop: for each supernode ksuper, compute new paths beginning
c       at ksuper in row subtrees and new paths in ij-intersection row subtrees.
c       ************************************************************************
c
c       --------------------------------------------
c       for each supernode ksuper (in postorder) ...
c       --------------------------------------------
        do  ksuper = 1, nsuper2
            oldsup = sperm(ksuper)
            fstnod = xsuper(oldsup)
            oldfst = perm(fstnod)
c           -------------------------------------------
c           for each node i (local ID) joined to ksuper
c           by a skeleton graph edge, update its 
c           "current leaf" to be ksuper.
c           also, determine where these skeleton graph
c           edges are (floc, lloc).
c           -------------------------------------------
            floc = fstloc(oldfst)
            lloc = xskadj(oldfst) + sklenf(oldfst) - 1
            if  ( floc .gt. lloc)  go to 500
            skipflag = .true.
            do  iloc = floc, lloc
                iold = skadj(iloc)
                inew = invp(iold)
                if  ( snode(inew) .ne. rootsup )  go to 50
                skipflag = .false.
                i = invp(iold) - i0
                curlf(i) = ksuper
            end do
            go to 90
   50       continue
            lloc = iloc - 1
c
c           -------------------------------------------
c           for each node i (local ID) joined to ksuper
c           by a skeleton graph edge ...
c           -------------------------------------------
   90       continue
            if  ( skipflag )  go to 500
            fstloc(oldfst) = lloc + 1
            do  iloc = floc, lloc
                iold = skadj(iloc)
                i = invp(iold) - i0
c
c               ***************************************
c               update row count for node i (local ID).
c               ***************************************
c
                pleafi = prvlf(i)
c               --------------------------------------------
c               if i has no previous leaf supernode then ...
c               --------------------------------------------
                if  ( pleafi .eq. 0 )  then
C                   -------------------------------------------
C                   ... ACCUMULATE ksuper-->nsuper2 PATH LENGTH 
C                       IN ROWCNT(i).
C                   -------------------------------------------
                    rowcnt(i) = rowcnt(i) + 
     &                          level(ksuper) - level(nsuper2)
                else
C                   ------------------------------------------
C                   ... OTHERWISE, LCA <-- FIND(PLEAFI), WHICH 
C                       IS THE LEAST COMMON ANCESTOR OF PLEAFI
C                       AND ksuper.
C                       (PATH HALVING.)
C                   ------------------------------------------
                    LAST1 = PLEAFI
                    LAST2 = SET(LAST1)
                    LCA = SET(LAST2)
  100               CONTINUE
                        IF  ( LCA .NE. LAST2 )  THEN
                            SET(LAST1) = LCA
                            LAST1 = LCA
                            LAST2 = SET(LAST1)
                            LCA = SET(LAST2)
                            GO TO 100
                        ENDIF
C                   --------------------------------------
C                   ACCUMULATE PLEAFI-->LCA PATH LENGTH IN 
C                   ROWCNT(i).
C                   --------------------------------------
                    ROWCNT(i) = ROWCNT(i)
     &                          + LEVEL(ksuper) - LEVEL(LCA)
                end if
c
c               ************************************************
c               update row counts for the ij-intersection trees,
c               when necessary.
c               ************************************************
c
                jright = lforw(i)
                jleft = lback(i)
c               -----------------------------------------
c               if i is not yet in the "visited" list ...
c               -----------------------------------------
                if  ( jleft .eq. -1 )  then
c                   ---------------------------------------
c                   ... start at the beginning of the list.
c                   ---------------------------------------
                    j = lhead
                else
c                   ---------------------------------------------
c                   ... otherwise, start immediately to the right
c                       of i.
c                   ---------------------------------------------
                    j = jright
                end if
c               ----------------------------------------------------
c               for each node j that we may need to find LCA for ...
c               ----------------------------------------------------
                do while  ( j .gt. 0 ) 
c                   -----------------------------------------
c                   if the "current leaf" of j is in ksuper's
c                   elimination subtree ...
c                   -----------------------------------------
                    if  ( fdesc(ksuper) .le. curlf(j) )  then
c                       ------------------------------------------------------
c                       ... then ksuper is a new ij-intersection subtree leaf.
c                       ------------------------------------------------------
                        nleafji = ksuper
                    else
C                       ------------------------------------------
C                       ... nleafji <-- FIND(PLEAFj), WHICH 
C                           IS THE LEAST COMMON ANCESTOR OF PLEAFj
C                           AND ksuper.
C                           (PATH HALVING.)
C                       ------------------------------------------
                        pleafj = prvlf(j)
                        LAST1 = PLEAFj
                        LAST2 = SET(LAST1)
                        nleafji = SET(LAST2)
  200                   CONTINUE
                            IF  ( nleafji .NE. LAST2 )  THEN
                                SET(LAST1) = nleafji
                                LAST1 = nleafji
                                LAST2 = SET(LAST1)
                                nleafji = SET(LAST2)
                                GO TO 200
                            ENDIF
                    end if
c                   --------------------------------------------------------
c                   get the "previous leaf" for the ij-intersection subtree.
c                   --------------------------------------------------------
                    if  ( j .lt. i )  then
                        pleafji = prv_l(j,i)
                    else
                        pleafji = prv_l(i,j)
                    end if
c                   --------------------------------------------
c                   if the previous ji-tree leaf and the new one
c                   are the same, there is nothing to do.
c                   --------------------------------------------
                    if  ( nleafji .eq. pleafji )  go to 400
c                   ----------------------------------------------------
c                   if the ij-intersection tree has no previous leaf ...
c                   ----------------------------------------------------
                    if  ( pleafji .eq. 0 )  then
C                       -------------------------------------------------
C                       ... ACCUMULATE nleafji-->nsuper2 PATH LENGTH 
C                           IN row_c(i,j) or row_c(j,i).
c                           make new ij-leaf nleafji the "previous leaf".
C                       -------------------------------------------------
                        if  ( j .lt. i )  then
                            row_c(i,j) = row_c(i,j) + 
     &                                   level(nleafji) - level(nsuper2)
                            prv_l(j,i) = nleafji
                        else
                            row_c(j,i) = row_c(j,i) + 
     &                                   level(nleafji) - level(nsuper2)
                            prv_l(i,j) = nleafji
                        end if
                    else
                        if  ( nleafji .lt. pleafji )  then
c                           ----------------------------------
c                           pleafji is an ancestor of nleafji.
c                           ----------------------------------
                            lca2 = pleafji
                        else
C                           --------------------------------------------
C                           ... otherwise, LCA2 <-- FIND(PLEAFji), WHICH 
C                               IS THE LEAST COMMON ANCESTOR OF PLEAFji
C                               AND nleafji.
C                               (PATH HALVING.)
C                           --------------------------------------------
                            LAST1 = PLEAFji
                            LAST2 = SET(LAST1)
                            LCA2 = SET(LAST2)
  300                       CONTINUE
                                IF  ( LCA2 .NE. LAST2 )  THEN
                                    SET(LAST1) = LCA2
                                    LAST1 = LCA2
                                    LAST2 = SET(LAST1)
                                    LCA2 = SET(LAST2)
                                    GO TO 300
                                ENDIF
                        end if
C                       ---------------------------------------------
C                       ACCUMULATE nleafji-->LCA2 PATH LENGTH IN 
C                       ROW_c(i,j) or row_c(j,i).
c                       make new ij-leaf nleafji the "previous leaf".
C                       ---------------------------------------------
                        if  ( j .lt. i )  then
                            row_c(i,j) = row_c(i,j) + 
     &                                   level(nleafji) - level(lca2)
                            prv_l(j,i) = nleafji
                        else
                            row_c(j,i) = row_c(j,i) + 
     &                                   level(nleafji) - level(lca2)
                            prv_l(i,j) = nleafji
                        end if
                    end if
c
  400               continue
c                   ------
c                   next j
c                   ------
                    j = lforw(j)
                end do
c               -----------------------------------------
c               if i is already in the "visited" list ...
c               -----------------------------------------
                if  ( jleft .ne. -1 )  then
c                   ---------------------------
c                   ... remove i from the list.
c                   ---------------------------
                    if  ( jleft .gt. 0 )  then
                        lforw(jleft) = jright
                    else
                        lhead = jright
                    end if
                    if  ( jright .gt. 0 )  then
                        lback(jright) = jleft
                    else
                        ltail = jleft
                    end if
                end if
c               ------------------------
c               add i to the tail of ...
c               ------------------------
                if  ( ltail .gt. 0 )  then
c                   --------------------
c                   ... a nonempty list.
c                   --------------------
                    lforw(ltail) = i
                    lback(i) = ltail
                    lforw(i) = 0
                    ltail = i
                else
c                   ------------------
c                   ... an empty list.
c                   ------------------
                    lhead = i
                    ltail = i
                    lforw(i) = 0
                    lback(i) = 0
                end if
c               ------
c               next i
c               ------
            end do
c           -------------------------------------------
c           for each node i (local ID) joined to ksuper
c           by a skeleton graph edge, update its 
c           "previous leaf".
c           -------------------------------------------
            do  iloc = floc, lloc
                iold = skadj(iloc)
                i = invp(iold) - i0
                prvlf(i) = ksuper
            end do
c           -----------------------------------
c           union operation for DSU operations.
c           -----------------------------------
  500       continue
            parent = suppar(ksuper)
            set(ksuper) = parent
c       -----------
c       next ksuper
c       -----------
        end do
c
c       ----------------------
c       compute the distances.
c       (make them symmetric.)
c       ----------------------
        do  i = 1, nnodes
            row_c(i,0) = rowcnt(i)
            row_c(0,i) = rowcnt(i)
            do  j = 1, i-1
                row_c(i,j) = rowcnt(i) + rowcnt(j) - 2 * row_c(i,j)
                row_c(j,i) = row_c(i,j)
            end do
        end do
C
        RETURN
      END
