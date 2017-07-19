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
C************* ordsup_ind_tsp_paths2                      **************
C***********************************************************************
C***********************************************************************
C
C   PURPOSE: 
c       This routine performs a reordering within supernodes using
c       the Bordeaux approach based on the Traveling Salesman Problem.
c       It has compression of the skeleton graph.
C
C   INPUT PARAMETERS:
c       (i) nadj        -   length of adjacency structure.
C       (I) NEQNS       -   NUMBER OF EQUATIONS.
C       (I) NOFSUB      -   NUMBER OF SUBSCRIPTS TO STORED IN
C                           LINDX(*).
C       (I) NSUPER      -   NUMBER OF SUPERNODES.
c       (i) supsiz      -   maximum supernode size.
C       (I) XSUPER(*)   -   ARRAY OF LENGTH NSUPER+1, CONTAINING THE
C                           FIRST COLUMN OF EACH SUPERNODE.
C       (I) XLINDX(*)   -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS 
C                           INTO THE SUBSCRIPT VECTOR.
C       (I) LINDX(*)    -   ARRAY OF LENGTH nofsub, CONTAINING THE
C                           COMPRESSED SUBSCRIPTS.
C       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING
C                           SUPERNODE MEMBERSHIP.
c       (i) xadj(*)     -   array of length neqns+1, containing pointers
c                           to the adjacency structure.
c       (i) adjncy(*)   -   array of length xadj(neqns+1)-1, containing
c                           the adjacency structure.
c       (i) etpar(*)    -   array of length neqns, containing the parent 
c                           vector of the elimination tree.
c
C   updated PARAMETERS:
C       (I) PERM(*)     -   ARRAY OF LENGTH NEQNS CONTAINING THE
C                           POSTORDERING initially.  On output,
c                           it contains a reordering of the nodes
c                           within each supernode.
C       (I) INVP(*)     -   ARRAY OF LENGTH NEQNS CONTAINING THE
C                           INVERSE OF THE POSTORDERING on input.
c                           on output, it contains the inverse of the
c                           new ordering in perm(*).
C
C   Output PARAMETER:
c       (i) iflag       -   error flag
c                              0 = sufficient working storage (>= 8*neqns+3 )
c                              -1 = insufficient working storage
c
C   WORKING PARAMETERS:
c       (i) xskadj(*)   -   array of length neqns+1, containing the 
c                           pointers for the skeleton graph.
c       (i) sklenf(*)   -   array of length neqns, containing the 
c                           lengths of the forward skeleton lists.
c       (i) sklenb(*)   -   array of length neqns, containing the 
c                           lengths of the back skeleton lists.
c       (i) skadj(*)    -   array of length nadj = xadj(neqns+1)-1, 
c                           containing the skeleton adjacency structure.
c       (i) invp2(*)    -   array of length neqns that records the inverse
c                           of the new permutation mapping current numbers to
c                           new numbers.
c       (i) link(*)     -   array of length neqns, contains the supernodal
c                           row subtrees as a linked list.
c       (i) fstloc(*)   -   array of length nsuper, contains current first
c                           location within index lists.
c       (i) sperm(*)    -   array of length nsuper, maps a supernode's 
c                           postordering number to it original number.
c       (i) fstloc2(*)  -   array of length neqns containing pointers to
c                           the next segments of skeleton neighbors to use.
c       (i) dist1(*,*)  -   two dimensional array dist1(0:supsiz,0:supsiz),
c                           used to store the TSP distances.
c       (i) suppar(*)   -   array of length nsuper containing the parent
c                           vector for the posordered supernodal 
c                           elimination subtree.
c       (i) iwsiz       -   This is an input parameter.  It is the size
c                           of working storage.  It must be no smaller 
c                           than 8*neqns+3.
c       (i) iwork       -   array of length iwsize, used as working storage
c                           by various subroutines in various ways.
c       (i) rep         -   array of length neqns, used to point to a
c                           node's representative node.
C
C***********************************************************************
C
      SUBROUTINE  ordsup_ind_tsp_paths2
     &                   (  nadj  , NEQNS , nofsub, nsuper, supsiz,
     &                      xsuper, xlindx, lindx , snode , xadj  , 
     &                      adjncy, etpar , perm  , invp  , iflag , 
     &                      xskadj, sklenf, sklenb, skadj , invp2 , 
     &                      link  , fstloc, sperm , fstloc2, dist1, 
     &                      suppar, iwsiz , iwork , rep             )
C
C***********************************************************************
C
C       -----------
C       PARAMETERS.
C       -----------
c
        INTEGER             iflag , iwsiz, nadj  , neqns, nofsub, 
     &                      nsuper, supsiz
c
        INTEGER             
     &                      adjncy(nadj)  ,
     &                      dist1(0:supsiz,0:supsiz),
     &                      etpar(neqns)  , fstloc(nsuper), 
     &                      fstloc2(neqns), invp(neqns)   , 
     &                      invp2(neqns)  , iwork(iwsiz)  , 
     &                      lindx(nofsub) , link(nsuper)  , 
     &                      perm(neqns)   , skadj(nadj)   , 
     &                      sklenb(neqns) , sklenf(neqns) , 
     &                      snode(neqns)  , sperm(nsuper) , 
     &                      suppar(nsuper), xadj(neqns+1) ,
     &                      xlindx(nsuper+1), xskadj(neqns+1),
     &                      xsuper(nsuper+1),
     &                      rep(neqns)
c
c       ----------------
c       local variables.
c       ----------------
c
        integer             colnum, i     , istart, istop , j     , 
     &                      jcol  , newsup, nnodes, nsuper2, nxtroot, 
     &                      nxtsup, offset, oldsup, rootsup
C
C***********************************************************************
C
        iflag = 0
        if  ( iwsiz .lt. 8*neqns + 3 )  then
            iflag = -1
            return
        end if
c
c       ---------------------------
c       compute the skeleton graph.
c       ---------------------------
        call  fnskel2 ( neqns , nadj  , xadj  , adjncy, perm  ,
     &                  invp  , etpar , xskadj, skadj , sklenf,
     &                  sklenb,
     &                  iwork(1),
     &                  iwork(neqns+2)
     &                )
c
c       --------------------------------------------
c       compress the skeleton graph (forward edges).
c       --------------------------------------------
        call  compress_skel2
     &                ( neqns , nadj  , snode , invp  ,
     &                  xskadj, sklenf, sklenb, skadj ,
     &                  rep   ,
     &                  iwork(1),
     &                  iwork(neqns+1),
     &                  iwork(2*neqns+1)
     &                )
c
c
c       ----------------------
c       initialize fstloc2(*).
c       ----------------------
        do  jcol = 1, neqns
            fstloc2(jcol) = xskadj(jcol)
        end do

c       -------------------
c       initialize link(*).
c       -------------------
        do  oldsup = 1, nsuper
            link(oldsup) = 0
        end do
c
c       ---------------------------------------------------------------
c       For each supernode rootsup from the bottom up in the supernodal
c       elimination forest ...
c       ---------------------------------------------------------------
c      
        colnum = 0
        do  rootsup = 1, nsuper
c           print *,'rootsup:',rootsup
c
c           -------------------------------------------------------
c           postorder the row supernodal elimination subtree rooted
c           at rootsup.
c           -------------------------------------------------------
            call  post_suptree2
     &               (  rootsup, NEQNS, nofsub, nsuper, xsuper, 
     &                  xlindx, lindx , snode , link  , fstloc, 
     &                  nsuper2, sperm, suppar, 
     &                  iwork(1),
     &                  iwork(nsuper+1),
     &                  iwork(2*nsuper+1),
     &                  iwork(3*nsuper+1),
     &                  iwork(4*nsuper+1)                        )
c
c           ----------------------------------------------
c           compute the distances for the TSP.
c           the subroutine partitions subtrees into paths.
c           ----------------------------------------------
            nnodes = xsuper(rootsup+1) - xsuper(rootsup)
            call  ordsup_tsp_paths2
     &               ( nnodes, nadj  , neqns , xskadj, sklenf, 
     &                 skadj , invp  , perm  , nsuper, xsuper, 
     &                 snode , nsuper2, sperm, suppar, fstloc2, 
     &                 dist1 ,
     &                 iwork(1),
     &                 iwork(nnodes+1),
     &                 dist1 ,
     &                 iwork(2*nnodes+1),
     &                 iwork(3*nnodes+1),
     &                 iwork(4*nnodes+1),
     &                 iwork(5*nnodes+1),
     &                 iwork(5*nnodes+nsuper2+2),
     &                 iwork(5*nnodes+2*nsuper2+3)               )
c
c           -----------------------------------------------------
c           approximate the solution to the TSP using the nearest
c           insertion heuristic.
c           -----------------------------------------------------
            call  near_ins2
     &               (  nnodes,
     &                  rootsup, nsuper2, NEQNS, dist1, nsuper, 
     &                  xsuper, colnum, invp2 , 
     &                  iwork(1),
     &                  iwork(nnodes+2),
     &                  iwork(2*nnodes+3),
     &                  iwork(2*nnodes+nsuper2+4),
     &                  iwork(3*nnodes+nsuper2+4),
     &                  iwork(4*nnodes+nsuper2+4),
     &                  iwork(5*nnodes+nsuper2+4),
     &                  iwork(6*nnodes+nsuper2+4), rep            )
c
c           -----------------------------------------------------
c           link the supernodes 1 through nsuper2 - 1 to the next
c           supernodal elimination subtree in which they appear.
c           also, update fstloc(*).
c           -----------------------------------------------------
            do  newsup = 1, nsuper2 - 1
                oldsup = sperm(newsup)
                if  ( fstloc(oldsup) .lt. xlindx(oldsup+1) )  then
                    nxtroot = snode(lindx(fstloc(oldsup)))
                    nxtsup = link(nxtroot)
                    link(oldsup) = nxtsup
                    link(nxtroot) = oldsup
                end if
            end do
c
c           ----------------------------------------------------------
c           link supernode rootsup to the first supernodal elimination
c           subtree in which it appears.  
c           Also initialize fstloc(rootsup).
c           ----------------------------------------------------------
            offset = xsuper(rootsup+1) - xsuper(rootsup)
            istart = xlindx(rootsup) + offset
            istop = xlindx(rootsup+1) - 1
            if  ( istart .le. istop )  then
                nxtroot = snode(lindx(istart))
                nxtsup = link(nxtroot)
                link(rootsup) = nxtsup
                link(nxtroot) = rootsup
                fstloc(rootsup) = istart
            end if
c
c           -----------------------
c           next supernode rootsup.
c           -----------------------
        end do
c
c       ----------------------------------------------------
c       compose the original ordering with the new ordering.
c       ----------------------------------------------------
        call  invinv ( neqns , invp  , invp2 , perm    )
c
        RETURN
      END
