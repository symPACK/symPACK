C***********************************************************************
C***********************************************************************
C
C   Version:        0.4
C   Last modified:  December 27, 1994
C   Authors:        Esmond G. Ng and Barry W. Peyton
C
C   Mathematical Sciences Section, Oak Ridge National Laboratory
C
C***********************************************************************
C***********************************************************************
C****     ORDNAT ..... NATURAL ORDERING                     ************
C***********************************************************************
C***********************************************************************
C
C     PURPOSE - THIS ROUTINE RECORDS THE INITIAL ORDERING IN THE
C               ORDERING VECTORS PERM AND INVP.
C
C     INPUT PARAMETERS -
C        NEQNS  - NUMBER OF EQUATIONS.
C
C     OUTPUT PARAMETERS -
C        PERM   - THE "NATURAL" ORDERING; I.E., THE INITIAL
C                 ORDERING.
C        INVP   - THE INVERSE OF PERM.
C
C***********************************************************************
C
      subroutine ordnat  (  neqns , perm  , invp    )
c
c***********************************************************************
c
        integer    invp(neqns) ,    perm(neqns)
        integer    neqns
c
        integer     i
c
c***********************************************************************
c
        do  700  i = 1, neqns
            perm(i) = i
            invp(i) = i
  700   continue
        return
c
      end
