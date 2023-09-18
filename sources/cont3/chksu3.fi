 SUBROUTINE chksu3(nsurm,codes,surname,tsurf)

!     Generates array CODES(2,*)
!         (1) = number of times the surface is used as MASTER
!         (2) = number of times the surface is used as SLAVE
!     surname(*) = label of the surface
!     tsurf   = number of total surfaces encountered
!
!     Input data for this generation is taken from pair data base

 IMPLICIT NONE

 INTEGER (kind=4), INTENT (IN) :: nsurm  !dimension of the array
 INTEGER (kind=4), INTENT (OUT) :: codes(2,nsurm), tsurf
 CHARACTER(len=*), INTENT (OUT) :: surname(nsurm)

 !Local Variables
 TYPE (pair3_db), POINTER :: pair
 CHARACTER (len=mnam) :: sname
 INTEGER (kind=4) :: i,j,ipair

 !initializes
 tsurf = 0            !Initializes number of surfaces
 codes = 0            !Initializes surface usage
 !surname = '                              '   !Clean surface names (unnecessary)
 pair => headp                                      !point to first pair

 DO ipair=1,npair                                   !for each pair
   DO j=1,2
     IF( j == 1 )THEN
       sname = pair%master                          !master surface label
     ELSE
       sname = pair%slave                           !slave surface label
     END IF
     i = 1                                          !search index
     DO
       IF(i > tsurf)THEN                            !unexisting surface
         tsurf = tsurf + 1                          !new surface
         surname(tsurf) = sname                     !store name
         codes(j,i) = 1                             !set as used
         EXIT
       ELSE IF( TRIM(surname(i)) == TRIM(sname))THEN            !existing
         codes(j,i) = codes(j,i) + 1                !increase counter
         EXIT
       END IF
       i = i+1                                    !increase index
     END DO
   END DO
   pair => pair%next                              !next pair
 END DO
 RETURN
 END SUBROUTINE chksu3
