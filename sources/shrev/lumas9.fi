 SUBROUTINE lumas9(iwrit,ntype,nelem,nnode,ngaus,shape,weigh, &
                   emass,sumat,heade)
 !********************************************************************
 !
 !***  evaluates mass matrix for beam/shell elements
 !
 !********************************************************************
 IMPLICIT NONE

 INTEGER (kind=4), INTENT(IN) :: ntype,nelem,nnode,ngaus,iwrit
 REAL (kind=8), INTENT(IN) :: shape(:,:),weigh(:)
 REAL (kind=8), INTENT(IN OUT) :: emass(:,:),sumat
 TYPE( ele09 ), POINTER :: heade
!          local variables
 INTEGER (kind=4) ielem,i,g,n,matty,isec,osec
 REAL    (kind=8) sume1,tmass,auxvl,area,iner,trmas(3),thick

 INTEGER (kind=4), PARAMETER :: ndime = 2
 REAL (kind=8),PARAMETER :: twopi=6.283185307179586d0
 TYPE( ele09 ), POINTER :: e
 TYPE(section), POINTER :: sec

 sume1=0.0
 e => heade
 osec = -1
 DO ielem=1,nelem
   !          lumped mass for each element
   isec = e%matno
   !         evaluates jacobians
   !         section properties
   IF( isec /= osec )THEN             !present and previous material diff
     sec => psecs(isec)%p
     matty = sec%secdef(3)
     thick = sec%rprop(1)
     IF(matty == 1 .OR. matty == 2) THEN
       IF(ntype == 1)THEN
         area = sec%mtbas%prope(5)*sec%rprop(4)
       ELSE
         area = sec%mtbas%prope(5)*thick
       END IF
     ELSE
       WRITE(lures,*,ERR=9999) ' material type not defined for beam/shell ',matty
     END IF
     osec  = isec
   END IF
   IF(e%jac(1)/thick < SQRT(1d0/6d0)) THEN
      iner = area*thick**2/12d0   !standard
   ELSE
      iner = area*e%jac(1)**2/2d0           !Hughes correction
   END IF
   !         numerical integration of nodal shape functions
   tmass=0.0                               !initializes total mass
   trmas(1:nnode) = 0d0                    !initializes translational mass
   DO g=1,ngaus                            !for each integration point
     auxvl = e%jac(g)*weigh(g)               !
     IF(ntype == 3) auxvl = auxvl*e%r0(g)*twopi  !for shell of revolution
     tmass = tmass+auxvl                       !sum on total mass
     DO n=1,nnode                          !for each node
       trmas(n) = trmas(n) + auxvl*shape(n,g)  !translational mass
     END DO
     IF(ntype == 3 .AND. nnode == 3) THEN  !shell of revol. and 3-node
       IF(trmas(3) > 20.d0 * trmas(1)) THEN   !very low mass in first node
          trmas(1) = trmas(1)+ 0.05*trmas(3)  !modify mass distribution
          trmas(3) = 0.95*trmas(3)
       ELSE IF(trmas(1) > 20.d0*trmas(3)) THEN !very low mass in third node
          trmas(3) = trmas(3)+ 0.05*trmas(1)   !modify mass distribution
          trmas(1) = 0.95*trmas(1)
       END IF
     END IF
   END DO
   !       sums over the global vector
   sume1 = sume1+tmass*area        !sum over the global mass
   DO n=1,nnode                    !for each node
     i = e%lnods(n)                !global node number
     emass(1:ndime,i) = emass(1:ndime,i)+trmas(n)*area !translational mass
      emass(ndime+1,i) = emass(ndime+1,i)+trmas(n)*iner !rotational mass
    END DO
    e => e%next
  END DO

  IF(iwrit == 1) &
  WRITE(lures,"(//'  Total mass for beam/shell element :',e15.7)",ERR=9999) sume1
  sumat = sumat+sume1

  RETURN
  9999 CALL runen2('')
  END SUBROUTINE lumas9
