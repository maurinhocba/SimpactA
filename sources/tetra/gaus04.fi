 SUBROUTINE gaus04(heade,coord,istop,gauss,angdf,nodvo)
 !-------------------------------------------------------------------
 !--- evaluates gauss points constants for solid element 04
 !--------------------------------------------------------------------
 IMPLICIT NONE

 INTEGER (kind=4), INTENT(OUT) :: istop
 REAL (kind=8), INTENT(IN) :: coord(:,:), & !initial coordinates
                              angdf(:)      !default angles Gl_X - Ort_Axes
 TYPE (ele04), POINTER :: heade    !pointer to first element in the set
 LOGICAL, INTENT (IN) :: nodvo    !TRUE if nodal volume is used
 LOGICAL, INTENT (IN OUT) :: gauss    !TRUE if initial constant exists

 INTEGER (kind=4) i,j,k,l,m,n
 REAL (kind=8) x(3,4)
 REAL (kind=8) :: ang(3),factor,rm(3,3),lb(3,3),vol
 TYPE (ele04), POINTER :: e

 INTERFACE
   INCLUDE 'lgst18.h'
 END INTERFACE
 !     ---------------------------------------------------------
 factor = ATAN(1d0)/45d0         !pi/180
 ang = angdf*factor
 CALL inrotm(ang,rm(1,1))
 e => heade       !point to first element
 DO               !loop over the elements in the set
   IF (.NOT.ASSOCIATED (e) ) EXIT  !if all element processed, exit
   x(:,1:4) = coord(:,e%lnods(1:4))
   IF(.NOT.gauss)THEN    !Compute initial constants if requested
     ! compute element local systems, cartesyan derivatives and element volume
     IF( ANY(e%angle /= 0d0))THEN
       ang = e%angle*factor
       CALL inrotm(ang,lb(1,1))
       lb = MATMUL(rm,lb)
     ELSE
       lb = rm
     END IF
     ! compute cartesian derivatives and volume factors
     CALL axep04(x,e%cd,e%dvol,lb)

     IF( e%dvol <= 0d0 ) THEN         !check that volume is possitive
       WRITE(*,1000)e%numel,e%dvol
       istop = 1
     END IF
     ! for nodal integrated elements compute averaged cartesian derivatives
     IF( nodvo )THEN
       vol = e%dvol     !element volume (weigthing factor
       DO i=1,nnode     !for each node
         n = e%lnods(i)        !node
         varn(0,n) = varn(0,n) + vol    !add element volume to nodal volume
         ncard(n)%p(1,:) = ncard(n)%p(1,:) + e%cd(i,:)*vol  !nodal contribution
         DO j=i+1,nnode                 !for every other node
           m = e%lnods(j)      !node
           k = 1               !search position of "m" in list at node i
           DO
             IF( nlnod(n)%p(k) == m )EXIT
             k = k + 1
           END DO
           l = 1
           DO                  !search position of "n" in list at node j
             IF( nlnod(m)%p(l) == n )EXIT
             l = l + 1
           END DO
           ncard(n)%p(k,:) = ncard(n)%p(k,:) + e%cd(j,:)*vol
           ncard(m)%p(l,:) = ncard(m)%p(l,:) + e%cd(i,:)*vol
         END DO
       END DO
     END IF

   END IF
   e => e%next          !point to next element
 END DO
 gauss = .TRUE.
 1000 FORMAT(' NEGATIVE OR NUL JACOBIAN AT 3-D TETRAHEDRA(04) ',I5,/ &
        ' JACOBIAN = ', E12.4 )

 RETURN
 END SUBROUTINE gaus04
