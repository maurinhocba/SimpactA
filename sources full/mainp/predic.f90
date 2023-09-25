 SUBROUTINE predic(istep,neq,disax,ddisp,dlamb,ipred)
 !***********************************************************************
 !
 !*** this routine predicts displacement and load step increments
 !    according to selected path and previous increments
 !
 !***********************************************************************
 IMPLICIT NONE
 !       routine arguments
 INTEGER (kind=4),INTENT(IN) :: istep, & !step to know how many increments exist
                                neq,   & !size of array ddisp and disax
                                ipred    !prediction order
 REAL (kind=8),INTENT(IN) :: dlamb       !increment in load/displacement factor
 REAL (kind=8),INTENT(IN) :: disax(:,:)  !previous incremental displacements
 REAL (kind=8),INTENT(OUT) :: ddisp(:)   !prediction of incremental displacements
 !       local variables
 INTEGER (kind=4) :: neq1
 REAL    (kind=8) :: r,s1,s2,f(3),l1,l2,l3,a(3,3),det,b(3)

 IF( istep == 1 )THEN ! no prediction
   ddisp = 0d0
 ELSE
   neq1 = neq+1       !size of array disax (position of load factor increment)
   !ip = 1                         !prediction candidate
   l1 = disax(neq1,1)/dlamb       !load factor increment from previous steps
   l2 = disax(neq1,2)/dlamb
   l3 = disax(neq1,3)/dlamb
   r = 1d0/l1                     !auxiliar
   ddisp = r  * disax(1:neq,1) ! linear prediction
   !!RETURN !higher predictions do not improve convergence
   !CALL residual(rfnorm(1))   !compute residual forces and norm
   !ibest = 1
   IF( istep > 2 .AND. ipred > 1)THEN                           ! quadratic prediction
     !ip = 2                      !prediction candidate
     s1 = (1d0-l2)/l1/(l1-l2)
     s2 = 1d0*(l1-1d0)/l2/(l1-l2)
     IF( ABS(s1+3d0) > 0.0001 .OR. ABS(s2-1d0) > 0.0001)PRINT *, 's1 ',s1,' s2 ',s2
     ddisp = disax(1:neq,1) * s1 + disax(1:neq,2) * s2
     !CALL residual(rfnorm(2))   !compute residual forces and norm
     !IF( rfnorm(2) < rfnorm(1) )   ibest = 2
     IF( istep > 3 .AND. ipred > 2)THEN                         ! cubic prediction
        !ip = 3                   !prediction candidate
        det = l1*l2*l3*(l1**2*l2-l1*l2**2+l2**2*l3-l2*l3**2+l3**2*l1-l3*l1**2)
        a(1,1) = (l2**2*l3-l2*l3**2)/det
        a(2,1) = (l3**2*l1-l3*l1**2)/det
        a(3,1) = (l1**2*l2-l1*l2**2)/det
        a(1,2) = (l2*l3**3-l2**3*l3)/det
        a(2,2) = (l3*l1**3-l3**3*l1)/det
        a(3,2) = (l1*l2**3-l1**3*l2)/det
        a(1,3) = (l2**3*l3**2-l2**2*l3**3)/det
        a(2,3) = (l3**3*l1**2-l3**2*l1**3)/det
        a(3,3) = (l1**3*l2**2-l1**2*l2**3)/det
        b = (/ 1d0,1d0,1d0 /)
        f = MATMUL(a,b)
        IF( ABS(f(1)+6d0) > 0.0001 .OR. ABS(f(2)-4d0) > 0.0001 .OR. ABS(f(3)+1d0) > 0.0001)THEN
           WRITE(55,"(4e12.4)")l3,l2,l1,dlamb
           WRITE(55,"(e12.4)")det
           WRITE(55,"(3f12.6)")a
           WRITE(55,"(3f10.4)")f
        END IF
        ddisp = disax(1:neq,1) * f(1) + disax(1:neq,2) * f(2) + disax(1:neq,3) * f(3)
        !CALL residual(rfnorm(3))   !compute residual forces and norm
        !IF( rfnorm(ibest) < rfnorm(3) )  ibest = 3
     END IF
     !IF( ip /= ibest )THEN
     !   IF( ibest == 1 ) ddisp = r  * disax(1:neq,1)
     !   IF( ibest == 2 ) ddisp = disax(1:neq,1) * s1 + disax(1:neq,2) * s2
     !END IF
   END IF

 END IF
 RETURN
 END SUBROUTINE predic
! !!!!!!!!!!!!!
! SUBROUTINE residual( neq,ddisp,dlamb,value )
! ! compute norm of residual forces to select prediction
! USE
! IMPLICIT NONE
! ! dummy arguments
! INTEGER(kind=4), INTENT(IN) :: neq
! REAL(kind=8), INTENT(IN) :: dlamb,ddisp(:)
! REAL(kind=8), INTENT(OUT) :: value
! ! local variables
! ! T A S K S
! !0- keep coordinates and euler systems to be restored
!   = coora
!   = euler
! !1- update coordinates COORA, COORC & EULER are used
! CALL coract(dlamb,velor(1:nvfix,nvelr+1),ddisp,coora,euler)
! !compute internal nodal forces
! CALL elemnt('RESVPL', deltc=1d0, istop=istop, ttime= ttime) !=> resid
! !2- compute contact nodal forces
! IF(numct > 0)THEN
!   fcont = 0d0
!   incdis = coora - coorc       ! check stack overflow
!   CALL contac('FORCES',iwrit,dtcal=1d0, ttime=ttime, neq=neq, &
!                velnp=incdis)    !toutd=toutd,
! END IF
! !3- Assemble internal, contact and external loads
! CALL ensve1(ndofn*npoin,ifpre(1,1),resid(1,1),eqres(1),npsdf(1),nesdf(1),ftsdf(1))
! !4- Add contact forces (FCONT)
! IF (numct > 0) CALL ensve1(ndimc*npoin,id(1,1),fcont(1,1),eqres(1),npsdf(1),nesdf(1),ftsdf(1))
! !5- add external equivalent forces (Force)
! IF (nload > 0) THEN
!   n1 = nload + 1
!   DO ieq=1,neq
!     eqres(ieq) = eqres(ieq) - force(ieq,n1)   !!!
!   END DO
! END IF
! !6- compute norm
! value = DOT_PRODUCT(eqres,eqres)
! !7- Restore coordinates and local systems
! coora =
! coorc =  !this is not changed
! euler =
! RETURN
! END SUBROUTINE residual
