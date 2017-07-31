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
C************* ordsup_ind ..... order within supernodes   **************
C***********************************************************************
C***********************************************************************
C
C   PURPOSE: 
c       This routine performs partition refinement within the current
c       supernode partition in an attempt to reduce the number of 
c       "segments" in the compressed indices.
C
C   INPUT PARAMETERS:
c       (i) altflag     -      0 - There is no alternating between moving the 
c                                  adjacent subset up and moving it down.
c                              1 - There is alternating between moving the 
c                                  adjacent subset up and moving it down.
C       (I) NEQNS       -   NUMBER OF EQUATIONS
C       (I) NOFSUB      -   NUMBER OF SUBSCRIPTS TO BE STORED IN
C                           LINDX(*).
C       (I) NSUPER      -   NUMBER OF SUPERNODES.
C       (I) XSUPER(*)   -   ARRAY OF LENGTH NSUPER+1, CONTAINING THE
C                           FIRST COLUMN OF EACH SUPERNODE.
C       (I) XLINDX      -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS 
C                           INTO THE SUBSCRIPT VECTOR.
C       (I) LINDX       -   ARRAY OF LENGTH MAXSUB, CONTAINING THE
C                           COMPRESSED SUBSCRIPTS.
C       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING
C                           SUPERNODE MEMBERSHIP.
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
C   WORKING PARAMETERS:
c       (i) freeforw(*) -   array of length neqns containing the forward
c                           links for the list of free empty sets.
c       (i) freeback(*) -   array of length neqns containing the backward
c                           links for the list of free empty sets.
c       (i) sforw(*)    -   array of length neqns containing the forward
c                           links for the list of sets in the current
c                           refined partition.
c       (i) sback(*)    -   array of length neqns containing the backward
c                           links for the list of sets in the current
c                           refined partition.
c       (i) setseg_forw(*) - array of length neqns containing the forward
c                            links in the list of "first sets" of segments
c                            of contiguous adjacent sets.
c       (i) setseg_back(*) - array of length neqns containing the backward
c                            links in the list of "first sets" of segments
c                            of contiguous adjacent sets.
c       (i) nodehead(*) -   array of length neqns containing the head of 
c                           the node lists for each set.
c       (i) nodeforw(*) -   array of length neqns containing the forward
c                           links for the node lists for each set.
c       (i) nodeback(*) -   array of length neqns containing the backward
c                           links for the node lists for each set.
c       (i) setsnode(*) -   array of length neqns containing the supernode
c                           number that a set is a subset of.
c       (i) mark(*)     -   array of length NEQNS+1 used to mark the sets
c                           adjacent to the current supernode being processed.
c       (i) set(*)      -   array of length neqns that records the set
c                           a node belongs to.
c       (i) compset(*)  -   array of length neqns used to record the companion
c                           set used to split a set adjacent to the current
c                           supernode being processed.
c       (i) invp2(*)    -   array of length neqns that record the inverse
c                           of the new permutation mapping current numbers to
c                           new numbers.
C
C***********************************************************************
C
      SUBROUTINE  near_ins
     &                   (  nnodes, rootsup, nsuper2, NEQNS, dist1, 
     &                      nsuper, xsuper, colnum, invp2 , lforw , 
     &                      lback , mdhead, mdforw, mdback, outforw,
     &                      outback, md                              )
C
C***********************************************************************
C
C       -----------
C       PARAMETERS.
C       -----------
c
        INTEGER             colnum, neqns , nnodes, nsuper, nsuper2, 
     &                      rootsup
c
        INTEGER             dist1(0:nnodes,0:nnodes),
     &                      invp2(neqns)     , lback(0:nnodes) ,
     &                      lforw(0:nnodes)  , md(nnodes)      ,
     &                      mdback(nnodes)   , mdforw(nnodes)  ,
     &                      mdhead(0:nsuper2), outback(nnodes) ,
     &                      outforw(nnodes)  , xsuper(nsuper+1)
c
c       ----------------
c       local variables.
c       ----------------
c
        integer             cost  , d     , i     , j     , j0    , 
     &                      lhead , mincost, mindist, next, nextj , 
     &                      nodecnt, nxtnod, oldd , outhead,
     &                      prev
C
C***********************************************************************
C
c       ******************************************
c       approximate the solution to the TSP by the 
c       nearest insersion method.
c       ******************************************
c
c       ---------------------------------------
c       initialize empty minimum distance bins.
c       ---------------------------------------
        do  d = 0, nsuper2
            mdhead(d) = 0
        end do
c
c       --------------------------------------------------------
c       place nodes 1,...,nnodes in the list of nodes not yet in
c       the insertion list.
c       --------------------------------------------------------
        do  i = nnodes-1, 1, -1
            outforw(i) = i+1
            outback(i+1) = i
        end do
        outforw(nnodes) = 0
        outback(1) = 0
        outhead = 1
c
c       ---------------------------------------------------
c       place vertex 0 alone in the nearest insersion list.
c       ---------------------------------------------------
c       print *,' '
c       print *,'insert 0 - begins insertion list'
        lhead = 0
        lback(0) = -1
        lforw(0) = -1
c
c       ----------------------------
c       initialize minimum distance.
c       ----------------------------
        mindist = nsuper2
c       -----------------------------------------------------
c       for each node i not in the nearest insertion list ...
c       -----------------------------------------------------
        do  i = 1, nnodes
c           ------------------------------------------
c           get distance from i to 0, and place node i
c           in its distance bin.
c           ------------------------------------------
            d = dist1(i,0)
            md(i) = d
            mindist = min ( mindist, d )
            nxtnod = mdhead(d)
            mdhead(d) = i
            mdforw(i) = nxtnod
            if  ( nxtnod .gt. 0 )  then
                mdback(nxtnod) = i
            end if
            mdback(i) = 0
        end do
c
c       i = 0
c
        nodecnt = nnodes
c       -------------------------------------------------------------------------
c       while there are still nodes to insert into the nearest insersion list ...
c       -------------------------------------------------------------------------
        do  while  ( nodecnt .ge. 1 )
c           i = i + 1
c           -----------------------------------------------------------
c           find a node i yet to be inserted and of minimum distance to 
c           the current list.
c           -----------------------------------------------------------
            do  while  ( mdhead(mindist) .eq. 0 )
                mindist = mindist + 1
            end do
            i = mdhead(mindist)
c           print *,'    next i to be inserted:',i
c           ---------------------------------------------------
c           remove node i from it current minimum distance bin.
c           ---------------------------------------------------
            nxtnod = mdforw(i)
            mdhead(mindist) = nxtnod
            if  ( nxtnod .gt. 0 )  then
                mdback(nxtnod) = 0
            end if
c           ---------------------------------------------------------
c           remove node i from the list of nodes currently not in the 
c           nearest insertion list.
c           ---------------------------------------------------------
            next = outforw(i)
            prev = outback(i)
            if  ( next .gt. 0 )  then
                outback(next) = prev
            end if
            if  ( prev .gt. 0 )  then
                outforw(prev) = next
            else
                outhead = next
            end if
c
c           -------------------------------------------------------
c           if i is the second node to be inserted into the nearest
c           insertion list ...
c           -------------------------------------------------------
            if  ( nodecnt .eq. nnodes )  then
c               -------------------------------------------------------
c               just append i to the end of the nearest insertion list.
c               -------------------------------------------------------
c               print *,'    i is second (on end)',i
                next = lhead
                lforw(next) = i
                lback(i) = next
                lforw(i) = -1
            else
c               ------------------------------------------------------
c               otherwise, we will insert i where it will increase the
c               the tour distance least.
c               ------------------------------------------------------
c               print *,' '
c               print *,'    two already in list - insert by cost'
                mincost = 2 * nsuper2
                j = lhead
c               -----------------------------------------------
c               for each successive pair j, nextj in the cycle.
c               -----------------------------------------------
                do  while  ( j .ge. 0 )
                    nextj = lforw(j)
                    if  ( nextj .eq. -1 )  then
c                       --------------------------
c                       j is last, nextj is first.
c                       --------------------------
                        nextj = lhead
                    end if
c                   print *,'        j,nextj:',j,nextj
c                   --------------------------------------------------
c                   compute the cost of inserting between j and jnext.
c                   --------------------------------------------------
                    cost = dist1(i,j) + dist1(i,nextj) 
     &                     - dist1(j,nextj)
c                   -------------------------------------
c                   if the new cost is an improvement ...
c                   -------------------------------------
                    if  ( cost .lt. mincost )  then
c                       -----------------------------------
c                       remember the cost and the location.
c                       -----------------------------------
                        mincost = cost
                        prev = j
                        next = nextj
c                       print *,'        cost:',cost
                    end if
c                   ---------------------------
c                   get next j to be processed.
c                   ---------------------------
                    if  ( nextj .eq. lhead )  then
                        j = -1
                    else
                        j = nextj
                    end if
                end do
c               ----------------------------------------------------
c               insert i into the least costly position found above.
c               ----------------------------------------------------
                lforw(prev) = i
                lback(i) = prev
                if  ( next .eq. lhead )  then
                    lforw(i) = -1
                else
                    lforw(i) = next
                    lback(next) = i
                end if
            end if
c           print *,' '
c
c           --------------------------------------------------------------------
c           for each node j not yet inserted into the nearest insertion iist ...
c           --------------------------------------------------------------------
c           print *,' '
c           print *,'    update distances'
            j = outhead
            do  while  ( j .gt. 0 )  
c               print *,'    j:',j
c               ---------------------------
c               get distance d from j to i.
c               ---------------------------
                d = dist1(i,j)
c               ------------------------------------------------
c               if the new distance improves the distance from j
c               to the current tour ...
c               ------------------------------------------------
                if  ( d .lt. md(j) )  then
c                   print *,'    d,oldd',d,md(j)
                    oldd = md(j)
c                   ---------------------------------------
c                   remove j from its minimum distance bin.
c                   ---------------------------------------
                    prev = mdback(j)
                    next = mdforw(j)
                    if  ( prev .eq. 0 )  then
                        mdhead(oldd) = next
                    else
                        mdforw(prev) = next
                    end if
                    if  ( next .gt. 0 )  then
                        mdback(next) = prev
                    end if
c                   -------------------------------------------
c                   insert j into its new minimum distance bin.
c                   -------------------------------------------
                    next = mdhead(d) 
                    mdforw(j) = next
                    if  ( next .gt. 0 )  then
                        mdback(next) = j
                    end if
                    mdback(j) = 0
                    mdhead(d) = j
c                   -------------------------
c                   update minimum distances.
c                   -------------------------
                    md(j) = d
                    mindist = min ( mindist, d )
                end if
c               ------
c               next j
c               ------
                j = outforw(j)
            end do
c           ------------------------------------------------
c           decrement counter for number of nodes processed.
c           ------------------------------------------------
            nodecnt = nodecnt - 1
        end do
c
c       -------------------------------------
c       order the nodes of supernode rootsup.
c       -------------------------------------
        j0 = xsuper(rootsup) - 1
        j = lforw(0)
        do  while  ( j .ge. 0 )  
c           --------------------------------------
c           order j next within supernode rootsup.
c           --------------------------------------
c           print *,'    j:',j
            colnum = colnum + 1
            invp2(j0+j) = colnum
            j = lforw(j)
        end do
c
        RETURN
      END
