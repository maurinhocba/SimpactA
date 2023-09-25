   SUBROUTINE rdmpcr (chead,ctail,nc,namex,namey)
   !
   ! read a serie of multi-point curves
   !
   IMPLICIT NONE
   TYPE (curve), POINTER :: chead,ctail      !head and tail pointers
   INTEGER (kind=4), INTENT(OUT) :: nc       !number of curves
   CHARACTER (len=*), INTENT(IN) :: namex, & !variable name
                                    namey    !variable name

   !local variables
   TYPE (curve), POINTER :: c
   INTEGER (kind=4) :: j,k
   REAL (kind=8) :: dx,dy

   nc = 0         !initializes number of curves
   DO
     CALL listen ('RDMPCR')                !read a line
     IF( exists('ENDCUR') )EXIT            !check end of curves
     IF( exists('DEFINI',k) )THEN          !check key-word
       ALLOCATE (c)                        !get memory for curve
       WRITE(lures,"(/)")
       c%name = get_name('DEFINI',posin=k,texts='!  CURVE NAME:') !curve name
       c%np=getint('NUMPTS',0,'!Number of Sample Points ........')
       ALLOCATE ( c%val(3,c%np) )          !get memory for curve values
       WRITE (lures,"(/,5x,'Point    Abscissa   Function Value ')",ERR=9999)
       DO j=1,c%np                         !read values
         CALL listen('RDCURV')
         IF (exists(TRIM(namex))) THEN
           c%val(1,j)=getrea(TRIM(namex),0d0,'!-Abscissa .........................')
           c%val(2,j)=getrea(TRIM(namey),0d0,'!-Function Value ...................')
           WRITE (lures,"(5x,i5,2e15.5)",ERR=9999)j,c%val(1:2,j)
         ELSE
           CALL runend('RDMPCR: Error in the curve input.  ')
         END IF
       END DO
       ! calculate slopes
       DO j = 1, c%np-1
         dx = c%val(1,j+1) - c%val(1,j)
         dy = c%val(2,j+1) - c%val(2,j)
         c%val(3,j) = dy/dx
       END DO
       c%val(3,c%np) = c%val(3,c%np-1)

       !  add to list of curves
       CALL add_cur (c, chead, ctail)
       nc = nc + 1          !increment counter
     ELSE
       CALL runend('RDMPCR: Curve definition expected ')
     END IF
   END DO

   RETURN
   9999 CALL runen2('')
   END SUBROUTINE rdmpcr
