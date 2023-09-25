 SUBROUTINE locla9(nelem,nnode,axesc,heade,iwrit,eule0,coord)
 !**********************************************************************
 !
 !      READ local systems for SHREV elements
 !
 !**********************************************************************
 IMPLICIT NONE
 INTEGER (kind=4) :: nelem,nnode,axesc,iwrit
 REAL (kind=8) eule0(:,:),coord(:,:)
 TYPE(ele09), POINTER :: heade

 ! local variables
 TYPE(ele09), POINTER :: e
 INTEGER (kind=4) i,j,k,m,n
 REAL    (kind=8) x(2),aux,f360
 REAL    (kind=8),PARAMETER :: facto = 1.74532925199433d-02
 LOGICAL :: first

 IF(axesc > 0) THEN
   CALL listen('LOCLA9')
   IF (.NOT.exists('SHRLOC')) CALL runend('ELMDA9: SHR_LOCAL_SYSTEM EXPECTED .')
   IF(iwrit == 1)WRITE(lures,"(/,'   local cartesian systems for SHREV elements')")

   e => heade       !point to first
   DO
     CALL listen('LOCLA9')
     IF (exists('ENDLOC')) EXIT

     i=getint('IELEM ',0,'!ELEMENT NUMBER ...................')
     first = .TRUE.
     DO
       IF(.NOT.ASSOCIATED( e ) .AND. .NOT.first)THEN
         WRITE(lures,"(' element label not found nelem=',i5)") i
         STOP
       ELSE IF(.NOT.ASSOCIATED( e ) .AND. first)THEN
         e => heade
         first = .FALSE.
       END IF
       IF( i == e%numel )EXIT
       e => e%next
     END DO

     DO j = 1,axesc
       ! Read angle between local and global system
       aux = param(j+1)
       IF(iwrit == 1) WRITE(lures,"(f10.5)")aux
       e%phil(j) = aux*facto
     END DO
     e => e%next
   END DO
 ELSE IF(axesc < 0) THEN
   axesc = -axesc
   f360 = 360d0*facto
   e => heade
   DO k=1,nelem
     m = e%lnods(1)     !first node
     n = e%lnods(nnode) !last node
     x(1:2) = coord(1:2,n) - coord(1:2,m)
     CALL vecuni(2,x,aux)
     IF(x(1) >= 0d0 .AND. x(2) >= 0d0)THEN !first
       aux = ASIN(x(2)) - (eule0(1,m)+eule0(1,n))/2d0
     ELSE IF(x(1) < 0d0 .AND. x(2) >= 0d0)THEN !second
       aux = ACOS(x(1)) - (eule0(1,m)+eule0(1,n))/2d0
     ELSE IF(x(1) <= 0d0 .AND. x(2) < 0d0)THEN !third
       aux = f360 - ACOS(x(1)) -(eule0(1,m)+eule0(1,n))/2d0
     ELSE IF(x(1) > 0d0 .AND. x(2) < 0d0)THEN !fourth
       aux = ASIN(x(2)) - (eule0(1,m)+eule0(1,m))/2d0
     END IF
     DO j = 1,axesc
       e%phil(j) = aux
     END DO
     e => e%next
   END DO
 END IF

 RETURN
 END SUBROUTINE locla9
