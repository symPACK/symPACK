      SUBROUTINE GRIDNDWRAPPER(P,Q)
        INTEGER*4       P,Q,IERROR
        INTEGER*4       GRID(P,Q)    , STACK(P,Q)
        CALL GRIDND ( P, Q, GRID, STACK, IERROR )
        WRITE (*,*) GRID
      END
 
      PROGRAM MAINBOX
        INTEGER*4       P,Q
        CHARACTER(len=32):: arg
        P = 0
        Q = 0


        IF ( IARGC() .eq. 2) THEN
          CALL getarg(1, arg)
          READ(arg,'(I4)'), P
          CALL getarg(2, arg)
          READ(arg,'(I4)'), Q
        END IF

        CALL GRIDNDWRAPPER ( P, Q )
      END
