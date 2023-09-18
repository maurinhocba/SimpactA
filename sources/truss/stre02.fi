 SUBROUTINE stre02 (sec, dstra, gausv, hyper, xl)

 !   stress computation for truss materials

 IMPLICIT NONE
 TYPE( section ), POINTER :: sec
 REAL (kind=8), INTENT(IN OUT) :: gausv(:),dstra
 REAL (kind=8), INTENT(IN) :: xl
 LOGICAL, INTENT(IN) :: hyper   !Hyper-elastic material (.TRUE.)
 REAL (kind=8) :: j

 INTERFACE
   INCLUDE 'stre02_1.h'
   INCLUDE 'stre02_6.h'
 END INTERFACE

 dstra = dstra + sec%rprop(2)    !adds initial strain

 SELECT CASE ( sec%mtbas%mtype )

 CASE (1,2,5) !standard elasto-plastic model
   CALL stre02_1 (sec%mtbas,dstra,gausv(2:),hyper,j)

 !CASE (3)
 !  ! 1D
 !  ! non-linear plastic
 !  CALL stre02_3 (sec%mtbas,dstra,gausv(2:),hype)

 CASE (6)  ! model for rubbers
   CALL stre02_6 (sec%mtbas,xl,gausv(2:),j)
 CASE DEFAULT
   WRITE(*,*)' matyp ', sec%mtbas%mtype
 END SELECT

 gausv(4) = gausv(3)*sec%rprop(1)  ! area

 !  compute force from stress * A0
 SELECT CASE (sec%iprop(1))
 CASE (0)
   gausv(4) = gausv(4)*xl    !2nd Piola-Kirchhoff
 CASE (1)
   gausv(4) = gausv(4)*j/xl  !Cauchy
 CASE (2)
   gausv(4) = gausv(4)/xl    !Kirchhoff
 CASE (3)
   !gausv(4) = gausv(4)       !1st P-K
 END SELECT

 RETURN
 END SUBROUTINE stre02
