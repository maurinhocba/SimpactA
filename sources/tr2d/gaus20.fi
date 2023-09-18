 SUBROUTINE gaus20(heade,coord,istop,gauss,angdf,ntype)
 !-------------------------------------------------------------------
 !--- evaluates gauss points constants for shell element 20
 !--------------------------------------------------------------------
 IMPLICIT NONE

 INTEGER (kind=4), PARAMETER :: ndime=2
 INTEGER (kind=4), INTENT(IN) :: ntype
 INTEGER (kind=4), INTENT(OUT) :: istop
 REAL (kind=8), INTENT(IN) :: coord(:,:), & !initial coordinates
                              angdf         !default angle X_1 - Ort_1
 TYPE (ele20), POINTER :: heade   !pointer to first element in the set
 LOGICAL, INTENT (IN OUT) :: gauss    !TRUE if initial constant exists

 INTEGER (kind=4) i,m,n,osec,isec
 REAL (kind=8) x(ndime,6),r1,r(3,2)
 REAL (kind=8) angle,factor,thick
 REAL (kind=8), PARAMETER :: twopi=6.283185307209586
 TYPE (ele20), POINTER :: e
 TYPE (section), POINTER :: sec
 LOGICAL :: sides(3)

 !     ---------------------------------------------------------
 factor = twopi/360d0         !pi/180

 osec = -1
 e => heade       !point to first element
 DO               !loop over the elements in the set
   IF (.NOT.ASSOCIATED (e) ) EXIT  !if all element processed, exit
   x(:,1:3) = coord(:,e%lnods(1:3))
   DO i =4,6                       !for all nodes in the patch
     n = e%lnods(i)                !node number
     m = i-3                       !associated side element
     IF(n /= 0) THEN               !if node exist
       x(:,i) = coord(:,n)         !get original coordinates
       sides(m) = .TRUE.           !side element exists
     ELSE
       sides(m) = .FALSE.          !side element does not exists
       x(:,i) = - x(:,m) + x(:,kk(1,m)) + x(:,kk(2,m))
     END IF
   END DO

   IF(.NOT.gauss)THEN    !Compute initial constants if requested
     IF( ntype == 3 .AND. ANY(x(1,1:3) < 0d0))THEN
       WRITE(55,"('For axilsymmetric problems all X1 (r) coordinates' &
              &   ' must be possitive ',/,                            &
              &   ' ELEMENT ',i6,/,(2E15.6) )",err=9999) e%numel,x
       CALL runen3('GAUS20: negative coord. r detected ')
     END IF
     IF( ntype == 3 )THEN
       DO i=1,3    !compute original r at each mid side
         r(i,1) = ( x(1,kk(1,i)) + x(1,kk(2,i)) )/2d0  !mid side radius
         IF( sides(i) )r(i,2) =  x(1,i+3) + r(i,1)*2d0 !side elm radius x 3
       END DO
     END IF
     ! compute element local systems, side proyections and element areas
     angle = ( angdf + e%angle )*factor
     CALL axep20(x,e%cd,e%area1,sides,angle,ntype,r)

     IF( e%area1 <= 0d0 ) THEN         !check that area is possitive
       WRITE(*,1000)e%numel,e%area1
       istop = 1
     END IF
     IF( ntype == 1 ) THEN  !plane stress
       isec = e%matno
       IF( isec /= osec )THEN
         sec => psecs(isec)%p
         thick = sec%rprop(1)
         osec = isec
       END IF
       e%area1= e%area1*thick
     ELSE IF( ntype == 3 )THEN ! include thickness or radius in the jacobian (element volume)
       r1 = SUM( x(1,1:3) )/3d0           !average radius from nodal points
       e%area1 = e%area1*r1*twopi
     END IF

   END IF
   e => e%next          !point to next element
 END DO
 gauss = .TRUE.
 1000 FORMAT(' NEGATIVE OR NUL JACOBIAN AT 2-D TRIANGLE  (20) ',I5,/ &
&       ' AREA = ', E12.4 )

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE gaus20
