      SUBROUTINE BOXNDWRAPPER(P,Q,R)
        INTEGER*4       P,Q,R,IERROR
        INTEGER*4       BOX(P,Q,R)    , STACK(P,Q,R)
        CALL BOXND ( P, Q, R, BOX, STACK, IERROR )
        WRITE (*,*) BOX
      END
 
      PROGRAM MAINBOX
        INTEGER*4       P,Q,R
        CHARACTER(len=32):: arg
        P = 0
        Q = 0
        R = 0


        IF ( IARGC() .eq. 3) THEN
          CALL getarg(1, arg)
          READ(arg,'(I4)'), P
          CALL getarg(2, arg)
          READ(arg,'(I4)'), Q
          CALL getarg(3, arg)
          READ(arg,'(I4)'), R
        END IF

        CALL BOXNDWRAPPER ( P, Q, R )
      END
