 SUBROUTINE mtype_c (matdef,mtype)
 !read material characteristics
 !to compute MTYPE
 USE ctrl_db, ONLY : ndime
 IMPLICIT NONE
 !dummy arguments
 INTEGER (kind=4), INTENT(OUT) :: matdef(:),mtype
 !local variables
 CHARACTER(len=midn), PARAMETER :: pars(8,7) = (/ &
 'RIGID ','ISOTRO','ORTHOT','RUBBER','POINTS','NONCOM','------','------', &
 'ELASTI','MISES ','HILL48','HILL79','HILL90','BBC03 ','NORM  ','------', &
 'IHNONE','IHLINE','IHLUDW','IHSATU','IHPOIN','IHHOLO','IHVOCE','------', &
 'KHNONE','KHLINE','KHSATU','KHPOIN','KINDAM','------','------','------', &
 'VISCOU','ELAPLA','ELVIPL','------','------','------','------','------', &
 'THERMA','------','------','------','------','------','------','------', &
 'ASSOCI','NONASO','------','------','------','------','------','------' /)

 ! Word ELASTI appeared twice in above array (eliminated from row 6)

 INTEGER(kind=4) :: i,j,k,l,m,n,ls,kk
 CHARACTER (len=63) :: cc,dd

 ! material space dimension
 matdef(1) = getint('NDIME ',ndime,' Material Space Dimension ........')
 ! check data input

 dd = ''
 n = 1
 DO i=1,7    ! for each characteristic
   k = 1     ! initializes string
   kk = 0    ! initializes counter
   cc = ''
   l = matdef(i+1)   ! initializes to default value
   DO j=1,8  ! for each possible value
     IF( exists(TRIM(pars(j,i)),m) )THEN     !check if key-word exists
       ls =  LEN_TRIM(pars(j,i))
       cc(k:k+ls) = TRIM(pars(j,i))//' '
       l = j
       k = k+ls+1
       kk = kk + 1
     END IF
   END DO
   IF( kk > 1 )THEN
     WRITE(lures,"('Invalid combination for one material ',A)",ERR=9999) TRIM(cc)
     CALL runend('MTYPE_C:Invalid combination for one material')
   END IF
   matdef(i+1) = l
   IF( l > 0 )THEN
     dd(n:n+ls+2) = TRIM(pars(l,i))//' & '
     n = n + ls+ 2
   END IF
 END DO

 WRITE(lures,"('         Material Definition ',A)",ERR=9999) TRIM(dd)

 mtype = getint('MATYPE',0,'!MATERIAL TYPE ..................')

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE mtype_c
