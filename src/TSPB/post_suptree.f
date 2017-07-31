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
      SUBROUTINE  post_suptree
     &                   (  rootsup, NEQNS, nofsub, nsuper, xsuper, 
     &                      xlindx, lindx , snode , link  , fstloc, 
     &                      nsuper2, sperm, suppar, parent, fchild, 
     &                      nxtsib, stack , sinvp                   )
C
C***********************************************************************
C
C       -----------
C       PARAMETERS.
C       -----------
c
        INTEGER             neqns, nofsub, nsuper, nsuper2, rootsup
c
        INTEGER             fchild(nsuper), fstloc(nsuper), 
     &                      lindx(nofsub) , link(nsuper)  ,
     &                      nxtsib(nsuper), parent(nsuper), 
     &                      sinvp(nsuper) , snode(neqns)  ,
     &                      sperm(nsuper) , stack(nsuper) ,
     &                      suppar(nsuper), xlindx(nsuper+1), 
     &                      xsuper(nsuper+1)
c
c       ----------------
c       local variables.
c       ----------------
c
        integer             fstsup, i     , inode , istart, istop ,
     &                      newpar, newsup, nsuper3, offset, oldsup,
     &                      oldpar, parsup 

C
C***********************************************************************
C
c       --------------------------------------------------
c       compute the parent vector for rootsup's supernodal
c       elimination subtree.
c       initialize fchild(*), nxtsib(*).
c       --------------------------------------------------
        parent(rootsup) = 0
        fchild(rootsup) = 0
        nxtsib(rootsup) = 0
        nsuper2 = 1
        oldsup = link(rootsup)
        do  while  ( oldsup .gt. 0 ) 
            nsuper2 = nsuper2 + 1
            fchild(oldsup) = 0
            nxtsib(oldsup) = 0
            offset = xsuper(oldsup+1) - xsuper(oldsup)
            istart = xlindx(oldsup) + offset
            istop = xlindx(oldsup+1) - 1
            if  ( istart .le. istop )  then
                parent(oldsup) = snode(lindx(istart))
            end if
            oldsup = link(oldsup)
        end do
c
c       ---------------------------------------
c       for every supernode oldsup in rootsup's
c       supernodal elimination subtree ...
c       ---------------------------------------
        oldsup = link(rootsup)
        do  while  ( oldsup .gt. 0 ) 
c           -------------------------------------------------
c           calculate "one beyond" the location of the last
c           node in oldsup's segment within rootsup (fstloc).
c           -------------------------------------------------
            istart = fstloc(oldsup)
            istop = xlindx(oldsup+1) - 1
            do  i = istart, istop
                inode = lindx(i)
                if  ( snode(inode) .ne. rootsup )  go to 100
            end do
            fstloc(oldsup) = istop + 1
            go to 200
  100       fstloc(oldsup) = i
  200       continue
c           -----------------------------------------------
c           insert oldsup into binary representation of the 
c           subtree.
c           -----------------------------------------------
            parsup = parent(oldsup)
            fstsup = fchild(parsup)
            fchild(parsup) = oldsup
            nxtsib(oldsup) = fstsup
            oldsup = link(oldsup)
        end do
c
c       ---------------------------------------------
c       postorder the supernodal elimination subtree.
c       ---------------------------------------------
        call sup_etpost
     &           (  rootsup, parent, fchild, nxtsib, stack ,
     &              nsuper3, sperm , sinvp                   )
c
c       ------------------
c       consistency check.
c       ------------------
c       if  ( nsuper3 .ne. nsuper2 )  then
c           print *, 'nsuper2,nsuper3:', nsuper2,nsuper3
c           stop
c       end if
c
c       -------------------------------------------------------------
c       translate supernodal elimination subtree to use the labels of
c       the postordering (in suppar).
c       -------------------------------------------------------------
        do  newsup = 1, nsuper2
            oldsup = sperm(newsup)
            oldpar = parent(oldsup)
            if  ( oldpar .eq. 0 )  then
                newpar = 0
            else
                newpar = sinvp(oldpar)
            end if
            suppar(newsup) = newpar
        end do
c
        RETURN
      END
