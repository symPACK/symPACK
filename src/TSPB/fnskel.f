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
C**************     fnskel  ..... FIND skeleton graph    ***************
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
C
C   WORK PARAMETERS:
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
      SUBROUTINE fnskel  (  NEQNS , ADJLEN, XADJ  , ADJNCY, PERM  ,
     &                      INVP  , ETPAR , xskel , skel  , sklenf,
     &                      sklenb, FDESC , PRVNBR                  )
C
C       -----------
C       PARAMETERS.
C       -----------
        INTEGER             ADJLEN, NEQNS
        INTEGER             
     &                      ADJNCY(ADJLEN)  , ETPAR(NEQNS)    , 
     &                      FDESC(0:NEQNS)  , INVP(NEQNS)     ,
     &                      PERM(NEQNS)     , PRVNBR(NEQNS)   , 
     &                      skel(adjlen)    , sklenb(neqns)   , 
     &                      sklenf(neqns)   , XADJ(neqns+1)   , 
     &                      xskel(neqns+1) 
C
C       ----------------
C       LOCAL VARIABLES.
C       ----------------
        INTEGER             HINBR , IFDESC, J     , JSTOP , JSTRT , 
     &                      K     , LOWNBR, oldhi , oldlow, OLDNBR, 
     &                      PARENT
C
C***********************************************************************
C
C       ----------------------------
C       COMPUTE FDESC(*), NCHILD(*).
C       INITIALIZE PRVNBR(*).
C       ----------------------------
        DO  100  K = NEQNS, 1, -1
            xskel(k) = xadj(k)
            sklenf(k) = 0
            sklenb(k) = 0
            FDESC(K) = K
            PRVNBR(K) = 0
  100   CONTINUE
        xskel(neqns+1) = xadj(neqns+1)
        FDESC(0) = 0
        DO  200  K = 1, NEQNS
            PARENT = ETPAR(K)
            IFDESC = FDESC(K)
            IF  ( IFDESC .LT. FDESC(PARENT) )  THEN
                FDESC(PARENT) = IFDESC
            ENDIF
  200   CONTINUE
C       ------------------------------------
C       FOR EACH ``LOW NEIGHBOR'' LOWNBR ...
C       ------------------------------------
        DO  600  LOWNBR = 1, NEQNS
            IFDESC = FDESC(LOWNBR)
            OLDNBR = PERM(LOWNBR)
            JSTRT = XADJ(OLDNBR)
            JSTOP = XADJ(OLDNBR+1) - 1
C           -----------------------------------------------
C           FOR EACH ``HIGH NEIGHBOR'', HINBR OF LOWNBR ...
C           -----------------------------------------------
            DO  500  J = JSTRT, JSTOP
                oldhi = adjncy(j)
                HINBR = INVP(oldhi)
                IF  ( HINBR .GT. LOWNBR )  THEN
                    IF  ( IFDESC .GT. PRVNBR(HINBR) )  THEN
C                       ----------------------------------------
C                       lownbr is a leaf in hinbr's row subtree.
C                       ----------------------------------------
                        sklenb(oldhi) = sklenb(oldhi) + 1
                        skel(xskel(oldhi+1)-sklenb(oldhi)) = oldnbr
                    ENDIF
C                   --------------------------------------------------
C                   LOWNBR NOW BECOMES ``PREVIOUS NEIGHBOR'' OF HINBR.
C                   --------------------------------------------------
                    PRVNBR(HINBR) = LOWNBR
                ENDIF
  500       CONTINUE
  600   CONTINUE
c
c       *********************************************************
c       create the forward skeleton lists so that they are sorted
c       in ascending order.
c       *********************************************************
c
c       ---------------------------------------------
c       for each "high neighbor" hinbr (in order) ...
c       ---------------------------------------------
        do  700  hinbr = 1, neqns
            oldhi = perm(hinbr)
            jstrt = xskel(oldhi+1) - 1
            jstop = xskel(oldhi+1) - sklenb(oldhi)
            do  j = jstrt, jstop, -1
                oldlow = skel(j)
c               -------------------------------------------------------
c               hinbr is joined to its next "low neighbor" in its list.
c               -------------------------------------------------------
                sklenf(oldlow) = sklenf(oldlow) + 1
                skel(xskel(oldlow)+sklenf(oldlow)-1) = oldhi
            end do
  700   continue
C
        RETURN
      END
