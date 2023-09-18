SUBROUTINE db3con
  !-----------------------------------------------
  !   read drawbead-sheet pair control information
  !-----------------------------------------------

  !   m o d u l e s
  USE c_input

  IMPLICIT NONE

  !   l o c a l   v a r i a b l e s

  INTEGER :: j
  TYPE (db3pair_db), POINTER :: pair
  !-----------------------------------------------------

  WRITE (lures, 40,ERR=9999)

  DO
     CALL listen ('db3con')
     IF (exists('ENDINT')) EXIT
     CALL new_db3(pair) !ALLOCATE (pair)
     IF (exists('PAIR  ',j)) THEN
       pair%pname = TRIM(words(j+1))
       WRITE(lures,"(/'Drawbead Pair name:',A)") TRIM(pair%pname)
     ELSE
       CALL runend ('Db3con: Drawb. Pair name not given!')
     END IF

     pair%rparm (1) = getrea ('TMAX  ', 0d0, '!Maximum drawbead force ...........')
     pair%rparm (2) = getrea ('ELMOD ', 0d0, '!Elastic modulus of the drawbd ....')
     pair%iparm (5) = getint ('NNODE ', 4, '?Number of nodes per sheet elements')
     pair%iparm (6) = getint ('IDBFOU', 0, '!Control code of history output ...')

     CALL add_db3pair (pair, headp, tailp)
  END DO

  RETURN

!30 FORMAT(' Error read drawbead control card for drawbead-line #',i3)
40 FORMAT (/, /, 4x, 'D r a w b e a d s   i n f o r m a t i o n' /)
!70 FORMAT(/ &
!        '       Drawbead-line data',  &
!         /, /,  &
!     /, '       Maximum drawbead force/unit length......   ', e15.5,  &
!     /, '       "Elastic modulus" of the drawbead.......   ', e15.5,  &
!     /, '          [force/length/displacement] ', /,  &
!        '       Hist. output code of restraining force..      ', i4,  &
!     /, '          0 : no, 1 : yes  ' /,  &
!        '       Number of nodes/element in sheet def.(3 or 4)', i8,  &
!        /)
 9999 CALL runen2('')
END SUBROUTINE db3con
