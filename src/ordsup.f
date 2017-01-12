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
C*************     ordsup ..... order within supernodes   **************
C***********************************************************************
C***********************************************************************
C
C   PURPOSE: 
c       This routine performs partition refinement within the current
c       supernode partition in an attempt to reduce the number of 
c       "segments" in the compressed indices.
C
C   INPUT PARAMETERS:
c       (i) ordflag     -      1 - The supernodes will be processed in the reverse
c                                  of the current labelling (the natural ordering).
c                              2 - The supernodes will be processed in the reverse
c                                  of a maximum cardinality search (MCS).
c                              3 - The supernodes will be processed in the reverse
c                                  order by the number of descendants they have
c                                  in the supernodal elimination tree (larger
c                                  number of descendants before lesser number
c                                  of descendants).
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
c       (i) supperm(*)  -   array of length nsuper containing the permutation
c                           of the supernodes that controls the order in
c                           which they are processed.
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
c       (i) heap(*)     -   array of length 2*nsuper used to store a heap
c                           of supernodes if they are reordered by subroutine
c                           mcssup.
C
C***********************************************************************
C
      SUBROUTINE  ordsup (  ordflag, altflag, NEQNS, nofsub, nsuper, 
     &                      xsuper, xlindx, lindx , snode , perm  , 
     &                      invp  , freeforw, freeback, sforw, sback, 
     &                      setseg_forw, setseg_back, nodehead, 
     &                      nodeforw, nodeback, 
     &                      setsnode, supperm, mark, set  , compset,
     &                      invp2 , heap                             )
C
C***********************************************************************
C
C       -----------
C       PARAMETERS.
C       -----------
        INTEGER             altflag, neqns, nofsub, nsuper, ordflag
        INTEGER             compset(neqns), freeback(neqns),
     &                      freeforw(neqns), heap(2*nsuper),
     &                      invp(neqns)   , invp2(neqns)  ,
     &                      lindx(nofsub) , mark(0:neqns) , 
     &                      nodeback(neqns), nodeforw(neqns), 
     &                      nodehead(neqns), perm(neqns)  ,
     &                      set(neqns)    , setseg_back(neqns),
     &                      setseg_forw(neqns), sback(neqns),
     &                      setsnode(neqns), sforw(neqns) , 
     &                      snode(neqns)  , supperm(nsuper), 
     &                      xlindx(nsuper+1), xsuper(nsuper+1)
c
c       ----------------
c       local variables.
c       ----------------
        integer             a     , after , b     , before, cursnode,
     &                      freehead, fstfree, fstnod, fstset, i    , 
     &                      inode , iset  , istart, istop , nextfree, 
     &                      nextnod, nextset, num , newsup, offset, 
     &                      oldsup, prevsnode, s  , shead , sprime, 
     &                      setseg_head, t, tprime, up
C
C***********************************************************************
C
c       ----------------------------------
c       initialize various vectors.
c       mark(*): marker vector for sets.
c       set(*): set membership for nodes.
c       compset(*): companion sets of sets.
c       -----------------------------------
        mark(0) = 0
        do  i = 1, neqns
            mark(i) = 0
            set(i) = snode(i)
            compset(i) = 0
        end do

c       -----------------------------------------------
c       order the supernodes in a topological ordering.
c       -----------------------------------------------
        if  ( ordflag .eq. 1 )  then
c           -----------------
c           natural ordering.
c           -----------------
            do  oldsup = 1, nsuper
                supperm(oldsup) = oldsup
            end do
        else 
c           --------------------------------------------------------------------------
c           ordflag = 2:
c           maximum cardinality search ordering (much like a minimum degree ordering).
c           ordflag = 3:
c           order supernodes by the number of descendants the have in the supernodal
c           elimination tree.
c           --------------------------------------------------------------------------
            call  mcssup ( ordflag, neqns , nofsub, nsuper, xsuper, 
     &                     snode , xlindx, lindx , supperm, sforw, 
     &                     sback , freeforw, heap, freeback, nodehead )
        end if
c
c       ----------------------------------------------------
c       construct initial list of available free empty sets.
c       ----------------------------------------------------
        freehead  = 0
        do  iset = nsuper+1, neqns
            nodehead(iset) = 0
            setsnode(iset) = 0
            fstfree = freehead
            freeforw(iset) = fstfree
            if  ( fstfree .gt. 0 )  then
                freeback(fstfree) = iset
            end if
            freehead = iset
            freeback(iset) = 0
        end do
c
c       -------------------------------------------------------------
c       compute the list of supernodes (sets), where each supernode's
c       parent is after the supernode in the list.
c       also record the supernode the set is a subset of.
c       -------------------------------------------------------------
        shead = 0
        do  newsup = nsuper, 1, -1
            oldsup = supperm(newsup)
            fstset = shead
            sforw(oldsup) = fstset
            if  ( fstset .gt. 0 )  then
                sback(fstset) = oldsup
            end if
            shead = oldsup
            sback(oldsup) = 0
            setsnode(oldsup) = oldsup
        end do
c
c       ------------------------------------------------------------
c       insert the nodes of each supernode (set) into its node list.
c       ------------------------------------------------------------
        do  newsup = 1, nsuper
            oldsup = supperm(newsup)
            nodehead(oldsup) = 0
            do  inode = xsuper(oldsup+1)-1, xsuper(oldsup), -1
                fstnod = nodehead(oldsup)
                nodeforw(inode) = fstnod
                if  ( fstnod .gt. 0 )  then
                    nodeback(fstnod) = inode
                end if
                nodehead(oldsup) = inode
                nodeback(inode) = 0
            end do
c           ----------------------------------------------------------------
c           mark the supernode (set) as a singleton if it has only one node.
c           (It will never need a companion set for partitioning.)
c           ----------------------------------------------------------------
            if  ( nodeforw(nodehead(oldsup)) .eq. 0 )  then
                compset(oldsup) = -1
            end if
        end do
c
c       -------------------------------------------------------------
c       For each supernode oldsup from the top down in the supernodal
c       elimination forest ...
c       -------------------------------------------------------------
c      
        do  newsup = nsuper, 1, -1
            oldsup = supperm(newsup)
c           --------------------------------------------------------
c           The list of "set segments" of sets adjacent to oldsup is 
c           initially empty.
c           --------------------------------------------------------
            setseg_head = 0
c
c           ********************************************************
c           While processing the nodes adjacent to supernode oldsup,
c           we detect the sets in the partition that may possibly
c           be split in two.
c           ********************************************************
c           ------------------------------------------------------------
c           For each node inode in the higher adjacency set of supernode
c           oldsup ...
c           ------------------------------------------------------------
            offset = xsuper(oldsup+1) - xsuper(oldsup)
            istart = xlindx(oldsup) + offset
            istop = xlindx(oldsup+1) - 1
            do  i = istart, istop
                inode = lindx(i)
c               --------------------------------------------------------
c               Get the current set that contains inode.
c               Also record the supernode that the set s is a subset of.
c               --------------------------------------------------------
                s = set(inode)
c               --------------------------------------------------------
c               If inode is the first node from set s adjacent to oldsup
c               and set s contains more than one node ...
c               --------------------------------------------------------
                if  ( compset(s) .eq. 0 )  then
c                   --------------------------------------------------------
c                   Get a free empty set sprime to pair with s and remove it 
c                   from the list of free sets.
c                   --------------------------------------------------------
                    sprime = freehead
                    nextfree = freeforw(sprime)
                    freehead = nextfree
                    if  ( nextfree .gt. 0 )  then
                        freeback(nextfree) = 0
                    end if
c                   --------------------------------------
c                   associate the free set sprime with s.
c                   sprime will be the companion set of s.
c                   --------------------------------------
                    compset(s) = sprime
                end if
c               ------------------------------------------------------------
c               If inode is the first node from set s adjacent to oldsup ...
c               ------------------------------------------------------------
                if  ( mark(s) .eq. 0 )  then
c                   ------------------------------------------------------
c                   modify the set segment "first sets" to account for the
c                   inclusion of s.
c                   ------------------------------------------------------
                    after = sforw(s)
                    before = sback(s)
                    if  ( after .ne. 0 )  then
                        if  ( mark(after) .eq. 1  .and.
     &                        setsnode(after) .eq. setsnode(s) )  then
c                           ----------------------------------------------
c                           if the next set after s is already included
c                           and a subset of the same supernode, then the
c                           inclusion of s means that "after" will no
c                           longer be first in a segment of sets adjacent
c                           to supernode oldsup.
c                           Hence, remove "after" from the list of segment
c                           "first sets".
c                           ----------------------------------------------
                            a = setseg_forw(after)
                            b = setseg_back(after)
                            if  ( b .eq. 0 )  then
                                setseg_head = a
                            else
                                setseg_forw(b) = a
                            end if
                            if  ( a .gt. 0 )  then
                                setseg_back(a) = b
                            end if
                        end if
                    end if
                    if  ( mark(before) .eq. 0  .or.
     &                    setsnode(before) .ne. setsnode(s) )  then
c                       ----------------------------------------------
c                       if the next set before s is not already
c                       included or it is a subset of a different
c                       supernode, then the inclusion of s means that
c                       s currently begins a segment of sets adjacent
c                       to the supernode oldsup.
c                       Hence add set s to the current list of segment
c                       "first sets".
c                       ----------------------------------------------
                        setseg_forw(s) = setseg_head
                        setseg_back(s) = 0
                        if  ( setseg_head .gt. 0 )  then
                            setseg_back(setseg_head) = s
                        end if
                        setseg_head = s
                    end if
c                   -------------------------
c                   mark set s as "included".
c                   -------------------------
                    mark(s) = 1
                end if
c               ------------------------------------------
c               if set s has more than one node,
c               then move node inode to the companion set.
c               ------------------------------------------
                if  ( compset(s) .gt. 0 )  then
c                   ----------------------------------
c                   get the companion set sprime of s.
c                   ----------------------------------
                    sprime = compset(s)
c                   ------------------------
c                   remove inode from set s.
c                   ------------------------
                    after = nodeforw(inode)
                    before = nodeback(inode)
                    if  ( before .eq. 0 )  then
                        nodehead(s) = after
                    else
                        nodeforw(before) = after
                    end if
                    if  ( after .gt. 0 )  then
                        nodeback(after) = before
                    end if
c                   ------------------------
c                   add inode to set sprime.
c                   ------------------------
                    fstnod = nodehead(sprime)
                    nodeforw(inode) = fstnod
                    if  ( fstnod .gt. 0 )  then
                        nodeback(fstnod) = inode
                    end if
                    nodeback(inode) = 0
                    nodehead(sprime) = inode
                    set(inode) = sprime
                end if
c               -----------
c               next inode.
c               -----------
            end do
c
c           *********************************************************
c           while processing the sets adjacent to supernode oldsup,
c           contiguous segment of sets by contiguous segment of sets, 
c           we refine the current partition.
c           *********************************************************
c           --------------------------------------------------
c           for each "first set" s of a contiguous set of sets 
c           adjacent to supernode oldsup ...
c           --------------------------------------------------
            s = setseg_head
            do while  ( s .gt. 0 )  
                t = s
                if  ( altflag .eq. 0 )  then
                    up = 0
                end if
c               ---------------------------------------------
c               for each set t in this contiguous set of sets
c               adjacent to supernode oldsup ...
c               ---------------------------------------------
                prevsnode = 0
                do while  ( t .gt. 0  .and.  mark(t) .eq. 1 )
                    cursnode = setsnode(t)
                    if  ( altflag .eq. 1  .and.
     &                    cursnode .ne. prevsnode )  then
                        up = 1
                    end if
                    prevsnode = cursnode
c                   ------------------------------------------------------------
c                   "before" precedes and "after" follows t in the list of sets.
c                   ------------------------------------------------------------
                    before = sback(t)
                    after = sforw(t)
c                   --------------------------------------------------------------
c                   if t is not a singleton set ...
c                   (Singleton sets are skipped here because they were not split.)
c                   --------------------------------------------------------------
                    if  ( compset(t) .ne. -1 )  then
c                       -------------------------------------------------
c                       get t's companion set and unmark t in compset(*).
c                       -------------------------------------------------
                        tprime = compset(t)
                        compset(t)= 0
c                       ----------------------------------------------------------
c                       if every node of set t is adjacent to supernode oldsup ...
c                       ----------------------------------------------------------
                        if  ( nodehead(t) .eq. 0 )  then
c                           -----------------------------------------------------------
c                           set t is empty; put it back in the pool of free empty sets.
c                           -----------------------------------------------------------
                            fstset = freehead
                            freeforw(t) = fstset
                            if  ( fstset .gt. 0 )  then
                                freeback(fstset) = t
                            end if
                            freehead = t
                            freeback(t) = 0
c                           -------------------------------------------------------
c                           the sequence is to be changed to before, tprime, after.
c                           insert tprime between before and after.
c                           -------------------------------------------------------
                            if  ( before .eq. 0 )  then
                                shead = tprime
                                sback(tprime) = 0
                            else
                                sforw(before) = tprime
                                sback(tprime) = before
                            end if
                            sforw(tprime) = after
                            if  ( after .gt. 0 )  then
                                sback(after) = tprime
                            end if
                            setsnode(tprime) = setsnode(t)
                            setsnode(t) = 0
                        else
c                           -------------------------------
c                           both t and tprime are nonempty.
c                           -------------------------------
                            if  ( up .eq. 1 )  then
c                               ----------------------------------------------------------
c                               the sequence is to be changed to before, t, tprime, after.
c                               insert tprime between t and after.
c                               ----------------------------------------------------------
                                sforw(t) = tprime
                                sback(tprime) = t
                                sforw(tprime) = after
                                if  ( after .gt. 0 )  then
                                    sback(after) = tprime
                                end if
                                if  ( altflag .eq. 1 )  then
                                    up = 0
                                end if
                            else
c                               ----------------------------------------------------------
c                               the sequence is to be changed to before, tprime, t, after.
c                               insert tprime between before and t.
c                               ----------------------------------------------------------
                                if  ( before .eq. 0 )  then
                                    shead = tprime
                                    sback(tprime) = 0
                                else
                                    sforw(before) = tprime
                                    sback(tprime) = before
                                end if
                                sforw(tprime) = t
                                sback(t) = tprime
                                if  ( altflag .eq. 1 )  then
                                    up = 1
                                end if
                            end if
c                           ----------------------------------------------------------
c                           mark sets t and tprime as singleton sets if they have only
c                           one node.
c                           ----------------------------------------------------------
                            if  ( nodeforw(nodehead(t)) 
     &                            .eq. 0                )  then
                                compset(t) = -1
                            end if
                            if  ( nodeforw(nodehead(tprime)) 
     &                            .eq. 0                     )  then
                                compset(tprime) = -1
                            end if
                            setsnode(tprime) = setsnode(t)
                        end if
                    end if
c                   --------------------------------------------
c                   next (possible) set t in segment is "after".
c                   also, unmark set t.
c                   --------------------------------------------
                    mark(t) = 0
                    t = after
                end do
c               --------------------------------
c               next "first set" s of a segment.
c               --------------------------------
                s = setseg_forw(s)
            end do
c
c           ----------------------
c           next supernode oldsup.
c           ----------------------
        end do
c
c       *************************************************************
c       compute the ordering.
c       invp2(*) maps the current node number to the new node number.
c       *************************************************************
c
c       ------------------------------------------------------
c       record the first numbers for a node in each supernode.
c       ------------------------------------------------------
        do  oldsup = 1, nsuper
            mark(oldsup) = xsuper(oldsup)
        end do
c
c       ---------------------------------------------------------
c       for each set nextset of the refined ordered partition ...
c       ---------------------------------------------------------
        nextset = shead
        do while  ( nextset .gt. 0 )
c           ----------------------------------------
c           for each node nextnod of set nextset ...
c           ----------------------------------------
            nextnod = nodehead(nextset)
            do while  ( nextnod .gt. 0 )  
                oldsup = snode(nextnod)
c               -------------------------------------------------
c               order nextnod next within its original supernode.
c               -------------------------------------------------
                invp2(nextnod) = mark(oldsup)
                mark(oldsup) = mark(oldsup) +  1
                nextnod = nodeforw(nextnod)
            end do
            nextset = sforw(nextset)
        end do
c       ----------------------------------------------------
c       compose the original ordering with the new ordering.
c       ----------------------------------------------------
        call  invinv ( neqns , invp  , invp2 , perm    )
c
        RETURN
      END
