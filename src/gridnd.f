* Esmond Ng, Oak Ridge National Laboratory
* Esmond Ng, Lawrence Berkeley National Laboratory
* Last modified:  September 24, 2006
* $Header$
* $Log$
************************************************************************
************************************************************************
************************************************************************
*****     GRIDND ..... NESTED DISSECTION ON A RECTANGULAR GRID     *****
************************************************************************
************************************************************************
*
*   PURPOSE -
*       This subroutine computes a theoretical nested dissection
*       ordering on a P-by-Q rectangular grid using depth-first
*       search.
*
*   INPUT PARAMETERS -
*       P       -   Number of rows in the grid.
*       Q       -   Number of columns in the grid.
*
*   OUTPUT PARAMETERS -
*       GRID    -   The inverse ordering.  It should be passed into
*                   this routine as a one-dimensional array of length
*                   P*Q.          
*
*   ERROR CODES -
*       IERROR  -   0 if no errors and 1 if STACK is not large enough.
*
*   WORKING STORAGE -
*       STACK   -   During execution, it is used as a stack.  Its size
*                   is P*Q.
*
*   NOTES -
*       There should be enough space in STACK if P*Q is larger than
*       25.
*
************************************************************************
*
       SUBROUTINE  GRIDND ( P, Q, GRID, STACK, STACKSZ, IERROR )
*
************************************************************************
*
        INTEGER*4       IERROR, P     , Q
        INTEGER*4       GRID(P,*)
        INTEGER*4       STACK(*)
*
        INTEGER*4       HALF  , I     , II    , ILAST , J     ,
     +                  JJ    , JLAST , LAST  , M     , MID   ,
     +                  N     , NNODES, TOP, STACKSZ
*
************************************************************************
*
        IERROR = 0
*       -----------------------------------------------------
*       Initialize the stack.
*
*       The first two items are the x and y coordinates
*       of the bottom left node of the subgrid.
*       the last two items are the numbers of nodes in
*       the x and y directions.
*
*       TOP points to the top of the stack.
*       LAST is the number of nodes remaining to be labelled.
*       -----------------------------------------------------
        LAST = P*Q
        NNODES = LAST
*
        TOP = 4
        IF  ( TOP .GT. STACKSZ )  THEN
            IERROR = 1
            RETURN
        END IF
        STACK(1) = 1
        STACK(2) = 1
        STACK(3) = P
        STACK(4) = Q
*
  100   CONTINUE
*           ------------------------------------
*           While the stack is not empty, do ...
*           ------------------------------------
            IF  ( TOP .NE. 0 )  THEN
*
*               -------------------------------
*               Pop the subgrid from the stack.
*               -------------------------------
                I = STACK(TOP-3)
                J = STACK(TOP-2)
                M = STACK(TOP-1)
                N = STACK(TOP)
                TOP = TOP - 4
*
                IF  ( M .GE. N )  THEN
*
                    IF  ( M .LE. 2 )  THEN
*
*                       -------------------------------
*                       The subgrid is too small.  No
*                       further dissection is possible.
*                       -------------------------------
                        ILAST = I + M - 1
                        JLAST = J + N - 1
                        DO  II = ILAST, I, -1
                            DO  JJ = JLAST, J, -1
                                GRID(II,JJ) = LAST
                                LAST = LAST - 1
                            END DO
                        END DO
*
                    ELSE
*
*                       ----------------------------
*                       Find a horizontal separator
*                       and dissect the subgrid into
*                       into two smaller subgrids.
*                       ----------------------------
                        HALF = (M - 1)/2
                        MID = I + HALF
*                       -------------------------
*                       Push the two new subgrids
*                       onto the stack.
*                       -------------------------
                        IF  ( TOP+8 .GT. STACKSZ )  THEN
                            IERROR = 1
                            RETURN
                        END IF
                        STACK(TOP+1) = I
                        STACK(TOP+2) = J
                        STACK(TOP+3) = HALF
                        STACK(TOP+4) = N
                        TOP = TOP + 4
                        STACK(TOP+1) = MID + 1
                        STACK(TOP+2) = J
                        STACK(TOP+3) = M - HALF - 1
                        STACK(TOP+4) = N
                        TOP = TOP + 4
*                       ---------------------------------
*                       Label the nodes on the separator.
*                       ---------------------------------
                        JLAST = J + N - 1
                        DO  JJ = JLAST, J, -1
                            GRID(MID,JJ) = LAST
                            LAST = LAST - 1
                        END DO
*
                    ENDIF
*
                ELSE
*
                    IF  ( N .LE. 2 )  THEN
*
*                       -------------------------------
*                       The subgrid is too small.  No
*                       further dissection is possible.
*                       -------------------------------
                        ILAST = I + M - 1
                        JLAST = J + N - 1
                        DO  II = ILAST, I, -1
                            DO  JJ = JLAST, J, -1
                                GRID(II,JJ) = LAST
                                LAST = LAST - 1
                            END DO
                        END DO
*
                    ELSE
*
*                       ----------------------------
*                       Find a vertical separator
*                       and dissect the subgrid into
*                       into two smaller subgrids.
*                       ----------------------------
                        HALF = (N - 1)/2
                        MID = J + HALF
*                       -------------------------
*                       Push the two new subgrids
*                       onto the stack.
*                       -------------------------
                        IF  ( TOP+8 .GT. STACKSZ )  THEN
                            IERROR = 1
                            RETURN
                        END IF
                        STACK(TOP+1) = I
                        STACK(TOP+2) = J
                        STACK(TOP+3) = M
                        STACK(TOP+4) = HALF
                        TOP = TOP + 4
                        STACK(TOP+1) = I
                        STACK(TOP+2) = MID + 1
                        STACK(TOP+3) = M
                        STACK(TOP+4) = N - HALF - 1
                        TOP = TOP + 4
*                       ---------------------------------
*                       Label the nodes on the separator.
*                       ---------------------------------
                        ILAST = I + M - 1
                        DO  II = ILAST, I, -1
                            GRID(II,MID) = LAST
                            LAST = LAST - 1
                        END DO
*
                    ENDIF
*
                ENDIF
*
                GO TO 100
*
            ENDIF
*
            RETURN
*
        END
