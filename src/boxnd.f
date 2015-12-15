* Esmond Ng, Oak Ridge National Laboratory
* Esmond Ng, Lawrence Berkeley National Laboratory
* Last modified:  September 24, 2006
* $Header$
* $Log$
************************************************************************
************************************************************************
************************************************************************
*****      BOXND ..... NESTED DISSECTION ON A 3D RECTANGULAR BOX   *****
************************************************************************
************************************************************************
*
*   PURPOSE -
*       This subroutine computes a theoretical nested dissection
*       ordering on a P-by-Q rectangular grid using depth-first
*       search.
*
*   INPUT PARAMETERS -
*       P       -   Number of rows in a plan of the box.
*       Q       -   Number of columns in a plan of the box.
*       R       -   Number of plans in the box.
*
*   OUTPUT PARAMETERS -
*       BOX     -   The inverse ordering.  It should be passed into
*                   this routine as a one-dimensional array of length
*                   P*Q*R.
*
*   ERROR CODES -
*       IERROR  -   0 if no errors and 1 if STACK is not large enough.
*
*   WORKING STORAGE -
*       STACK   -   During execution, it is used as a stack.  Its size
*                   is P*Q*R.
*
*   NOTES -
*       There should be enough space in STACK if P*Q is larger than
*       25.
*
************************************************************************
*
       SUBROUTINE  BOXND ( P, Q, R, BOX, STACK, STACKSZ, IERROR )
*
************************************************************************
*
        INTEGER*4       IERROR, P     , Q     , R
        INTEGER*4       BOX(P,Q,*)    , STACK(*)
*
        INTEGER*4       HALF  , I     , II    , ILAST , J     ,
     +                  JJ    , JLAST , LAST  , M     , MID   ,
     +                  N     , NNODES, TOP   , K     , KK    ,
     +                  KLAST , L     , STACKSZ
*
************************************************************************
*
        IERROR = 0
*       -----------------------------------------------------
*       Initialize the stack.
*
*       The first three items are the x, y and z coordinates
*       of the front bottom left node of the sub-box.
*       the last three items are the numbers of nodes in
*       the x, y and z directions.
*
*       TOP points to the top of the stack.
*       LAST is the number of nodes remaining to be labelled.
*       -----------------------------------------------------
        LAST = P*Q*R
        NNODES = LAST
*
        TOP = 6
        IF  ( TOP .GT. STACKSZ )  THEN
            IERROR = 1
            RETURN
        END IF
        STACK(1) = 1
        STACK(2) = 1
        STACK(3) = 1
        STACK(4) = P
        STACK(5) = Q
        STACK(6) = R
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
                I = STACK(TOP-5)
                J = STACK(TOP-4)
                K = STACK(TOP-3)
                M = STACK(TOP-2)
                N = STACK(TOP-1)
                L = STACK(TOP)
                TOP = TOP - 6
*
                IF  ( M .GE. N .AND. M .GE. L )  THEN
*
                    IF  ( M .LE. 2 )  THEN
*
*                       -------------------------------
*                       The subbox is too small.  No
*                       further dissection is possible.
*                       -------------------------------
                        ILAST = I + M - 1
                        JLAST = J + N - 1
                        KLAST = K + L - 1
                        DO  II = ILAST, I, -1
                            DO  JJ = JLAST, J, -1
                                DO  KK = KLAST, K, -1
                                  BOX(II,JJ,KK) = LAST
                                  LAST = LAST - 1
                                END DO
                            END DO
                        END DO
*
                    ELSE
*
*                       ----------------------------
*                       Find a separator 
*                       and dissect the sub-box into
*                       into two smaller sub-boxes.
*                       ----------------------------
                        HALF = (M - 1)/2
                        MID = I + HALF
*                       -------------------------
*                       Push the two new sub-boxess
*                       onto the stack.
*                       -------------------------
                        IF  ( TOP+12 .GT. STACKSZ )  THEN
                            IERROR = 1
                            RETURN
                        END IF
                        STACK(TOP+1) = I
                        STACK(TOP+2) = J
                        STACK(TOP+3) = K
                        STACK(TOP+4) = HALF
                        STACK(TOP+5) = N
                        STACK(TOP+6) = L
                        TOP = TOP + 6
                        STACK(TOP+1) = MID + 1
                        STACK(TOP+2) = J
                        STACK(TOP+3) = K
                        STACK(TOP+4) = M - HALF - 1
                        STACK(TOP+5) = N
                        STACK(TOP+6) = L
                        TOP = TOP + 6
*                       ---------------------------------
*                       Label the nodes on the separator.
*                       which is in the JK plane
*                       ---------------------------------
                        JLAST = J + N - 1
                        KLAST = K + L - 1
                        DO  JJ = JLAST, J, -1
                            DO  KK = KLAST, K, -1
                                BOX(MID,JJ,KK) = LAST
                                LAST = LAST - 1
                            END DO
                        END DO
*
                    ENDIF
*
                ELSE IF  ( N .GE. M .AND. N .GE. L )  THEN
*
                    IF  ( N .LE. 2 )  THEN
*
*                       -------------------------------
*                       The subbox is too small.  No
*                       further dissection is possible.
*                       -------------------------------
                        ILAST = I + M - 1
                        JLAST = J + N - 1
                        KLAST = K + L - 1
                        DO  II = ILAST, I, -1
                            DO  JJ = JLAST, J, -1
                                DO  KK = KLAST, K, -1
                                  BOX(II,JJ,KK) = LAST
                                  LAST = LAST - 1
                                END DO
                            END DO
                        END DO
*
                    ELSE
*
*                       ----------------------------
*                       Find a separator 
*                       and dissect the sub-box into
*                       into two smaller sub-boxes.
*                       ----------------------------
                        HALF = (N - 1)/2
                        MID = J + HALF
*                       -------------------------
*                       Push the two new sub-boxess
*                       onto the stack.
*                       -------------------------
                        IF  ( TOP+12 .GT. STACKSZ )  THEN
                            IERROR = 1
                            RETURN
                        END IF
                        STACK(TOP+1) = I
                        STACK(TOP+2) = J
                        STACK(TOP+3) = K
                        STACK(TOP+4) = M
                        STACK(TOP+5) = HALF
                        STACK(TOP+6) = L
                        TOP = TOP + 6
                        STACK(TOP+1) = I
                        STACK(TOP+2) = MID + 1
                        STACK(TOP+3) = K
                        STACK(TOP+4) = M
                        STACK(TOP+5) = N - HALF - 1
                        STACK(TOP+6) = L
                        TOP = TOP + 6
*                       ---------------------------------
*                       Label the nodes on the separator.
*                       which is in the IK plane
*                       ---------------------------------
                        ILAST = I + M - 1
                        KLAST = K + L - 1
                        DO  II = ILAST, I, -1
                            DO  KK = KLAST, K, -1
                                BOX(II,MID,KK) = LAST
                                LAST = LAST - 1
                            END DO
                        END DO
*
                    ENDIF
*
                ELSE
*
                    IF  ( L .LE. 2 )  THEN
*
*                       -------------------------------
*                       The subbox is too small.  No
*                       further dissection is possible.
*                       -------------------------------
                        ILAST = I + M - 1
                        JLAST = J + N - 1
                        KLAST = K + L - 1
                        DO  II = ILAST, I, -1
                            DO  JJ = JLAST, J, -1
                                DO  KK = KLAST, K, -1
                                  BOX(II,JJ,KK) = LAST
                                  LAST = LAST - 1
                                END DO
                            END DO
                        END DO
*
                    ELSE
*
*                       ----------------------------
*                       Find a separator 
*                       and dissect the sub-box into
*                       into two smaller sub-boxes.
*                       ----------------------------
                        HALF = (L - 1)/2
                        MID = K + HALF
*                       -------------------------
*                       Push the two new sub-boxess
*                       onto the stack.
*                       -------------------------
                        IF  ( TOP+12 .GT. STACKSZ )  THEN
                            IERROR = 1
                            RETURN
                        END IF
                        STACK(TOP+1) = I
                        STACK(TOP+2) = J
                        STACK(TOP+3) = K
                        STACK(TOP+4) = M
                        STACK(TOP+5) = N
                        STACK(TOP+6) = HALF
                        TOP = TOP + 6
                        STACK(TOP+1) = I
                        STACK(TOP+2) = J
                        STACK(TOP+3) = MID + 1
                        STACK(TOP+4) = M
                        STACK(TOP+5) = N
                        STACK(TOP+6) = L - HALF - 1
                        TOP = TOP + 6
*                       ---------------------------------
*                       Label the nodes on the separator.
*                       which is in the IJ plane
*                       ---------------------------------
                        ILAST = I + M - 1
                        JLAST = J + N - 1
                        DO  II = ILAST, I, -1
                            DO  JJ = JLAST, J, -1
                                BOX(II,JJ,MID) = LAST
                                LAST = LAST - 1
                            END DO
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
