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
C******     BETREE ..... BINARY TREE REPRESENTATION OF ETREE     *******
C***********************************************************************
C***********************************************************************
C
C   WRITTEN BY JOSEPH LIU (JUL 17, 1985)
C
C   PURPOSE:
C       TO DETERMINE THE BINARY TREE REPRESENTATION OF THE ELIMINATION
C       TREE GIVEN BY THE PARENT VECTOR.  THE RETURNED REPRESENTATION
C       WILL BE GIVEN BY THE FIRST-SON AND BROTHER VECTORS.  THE ROOT
C       OF THE BINARY TREE IS ALWAYS NEQNS.
C
C   INPUT PARAMETERS:
C       NEQNS           -   NUMBER OF EQUATIONS.
C       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE.
C                           IT IS ASSUMED THAT PARENT(I) > I EXCEPT OF
C                           THE ROOTS.
C
C   OUTPUT PARAMETERS:
C       FSON            -   THE FIRST SON VECTOR.
C       BROTHR          -   THE BROTHER VECTOR.
C
C***********************************************************************
C
      SUBROUTINE  BETREE (  NEQNS , PARENT, FSON  , BROTHR          )
C
C***********************************************************************
C
        INTEGER*4           BROTHR(*)     , FSON(*)       ,
     &                      PARENT(*)
C
        INTEGER*4           NEQNS
C
C***********************************************************************
C
        INTEGER*4           LROOT , NODE  , NDPAR
C
C***********************************************************************
C
        IF  ( NEQNS .LE. 0 )  RETURN
C
        DO  100  NODE = 1, NEQNS
            FSON(NODE) = 0
            BROTHR(NODE) = 0
  100   CONTINUE
        LROOT = NEQNS
C       ------------------------------------------------------------
C       FOR EACH NODE := NEQNS-1 STEP -1 DOWNTO 1, DO THE FOLLOWING.
C       ------------------------------------------------------------
        IF  ( NEQNS .LE. 1 )  RETURN
        DO  300  NODE = NEQNS-1, 1, -1
            NDPAR = PARENT(NODE)
            IF  ( NDPAR .LE. 0  .OR.  NDPAR .EQ. NODE )  THEN
C               -------------------------------------------------
C               NODE HAS NO PARENT.  GIVEN STRUCTURE IS A FOREST.
C               SET NODE TO BE ONE OF THE ROOTS OF THE TREES.
C               -------------------------------------------------
                BROTHR(LROOT) = NODE
                LROOT = NODE
            ELSE
C               -------------------------------------------
C               OTHERWISE, BECOMES FIRST SON OF ITS PARENT.
C               -------------------------------------------
                BROTHR(NODE) = FSON(NDPAR)
                FSON(NDPAR) = NODE
            ENDIF
  300   CONTINUE
        BROTHR(LROOT) = 0
C
        RETURN
      END
