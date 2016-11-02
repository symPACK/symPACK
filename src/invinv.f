C***********************************************************************
C***********************************************************************
C
C   Version:        0.4
C   Last modified:  December 27, 1994
C   Authors:        Joseph W.H. Liu
C
C   Mathematical Sciences Section, Oak Ridge National Laboratory
C
C***********************************************************************
C***********************************************************************
C***********     INVINV ..... CONCATENATION OF TWO INVP     ************
C***********************************************************************
C***********************************************************************
C
C   WRITTEN BY JOSEPH LIU (JUL 17, 1985)
C
C   PURPOSE:
C       TO PERFORM THE MAPPING OF
C           ORIGINAL-INVP --> INTERMEDIATE-INVP --> NEW INVP
C       AND THE RESULTING ORDERING REPLACES INVP.  THE NEW PERMUTATION
C       VECTOR PERM IS ALSO COMPUTED.
C
C   INPUT PARAMETERS:
C       NEQNS           -   NUMBER OF EQUATIONS.
C       INVP2           -   THE SECOND INVERSE PERMUTATION VECTOR.
C
C   UPDATED PARAMETERS:
C       INVP            -   THE FIRST INVERSE PERMUTATION VECTOR.  ON
C                           OUTPUT, IT CONTAINS THE NEW INVERSE
C                           PERMUTATION.
C
C   OUTPUT PARAMETER:
C       PERM            -   NEW PERMUTATION VECTOR (CAN BE THE SAME AS
C                           INVP2).
C
C***********************************************************************
C
      SUBROUTINE  INVINV (  NEQNS , INVP  , INVP2 , PERM            )
C
C***********************************************************************
C
        INTEGER*4           INVP(*)       , INVP2(*)      ,
     &                      PERM(*)
C
        INTEGER*4           NEQNS
C
C***********************************************************************
C
        INTEGER*4           I     , INTERM, NODE
C
C***********************************************************************
C
        DO  100  I = 1, NEQNS
            INTERM = INVP(I)
            INVP(I) = INVP2(INTERM)
  100   CONTINUE
C
        DO  200  I = 1, NEQNS
            NODE = INVP(I)
            PERM(NODE) = I
  200   CONTINUE
C
        RETURN
      END
