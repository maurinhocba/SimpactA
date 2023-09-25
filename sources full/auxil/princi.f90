 SUBROUTINE princi(ndime,i,ll)
 !
 ! Compute principal values and directions of a 3-D 2nd order tensor
 ! used to compute principal values and directions of INERTIA tensor
 ! NOT adequate for Cauchy-Green tensor
 ! not very efficient but called once in a run
 !
 IMPLICIT NONE
 INTEGER (kind=4),INTENT (IN) :: ndime !problem dimension, MUST be 3
 REAL (kind=8), INTENT(IN OUT) :: i(ndime,ndime), & ! inertia tensor in global system (IN)
                                  ll(ndime,ndime)   ! local system at the node (IN)

 REAL (kind=8),PARAMETER :: pi = 3.14159265358979
 INTEGER (kind=4) :: k,lg(1)
 REAL (kind=8) :: i1,i2,i3,i13,phi,lp(3),la,a(2,2),det,l(3,3),pr(3)

 !     compute mean value and deviator tensor
 i1   = (i(1,1) + i(2,2) + i(3,3))     ! I1
 i13  = i1/3d0                         ! I1/3
 i(1,1) = i(1,1) -i13                  ! deviator tensor
 i(2,2) = i(2,2) -i13
 i(3,3) = i(3,3) -i13
 !     compute principal values
 i2  =SQRT(4d0*( - i(1,1)*i(2,2) - i(1,1)*i(3,3) - i(2,2)*i(3,3)  &  ! (4/3 I2)^(1/2)
                 + i(1,2)**2     + i(1,3)**2     + i(2,3)**2 )/3d0)
 i3  = i(1,1)*i(2,2)*i(3,3) + 2*i(1,2)*i(2,3)*i(3,1)  &              !third inv (det)
      -i(1,1)*i(2,3)**2 - i(2,2)*i(1,3)**2 - i(3,3)*i(1,2)**2
 IF(i2/i1 < 1.e-6)THEN     !all eigenvalues are almost equal
   i(1,1) = i(1,1) + i13
   i(2,2) = i(2,2) + i13
   i(3,3) = i(3,3) + i13
   ll(:,1) = (/ 1d0, 0d0, 0d0 /)
   ll(:,2) = (/ 0d0, 1d0, 0d0 /)
   ll(:,3) = (/ 0d0, 0d0, 1d0 /)
   RETURN
 END IF
 phi = ASIN(-4d0*i3/(i2**3))/3d0          !angle
 lp(1) = i2*SIN(phi+2d0/3d0*pi)         !largest root
 lp(2) = i2*SIN(phi)                    !interm. root
 lp(3) = i2*SIN(phi+4d0/3d0*pi)         !smallest root
 !     compute firt and third principal directions
 DO k=1,3,2  !principal directions 1 & 3
   !  first option, use the first two equations
   a(1,1) = i(1,1) - lp(k)
   a(1,2) = i(1,2)
   a(2,1) = i(2,1)
   a(2,2) = i(2,2) - lp(k)
   det = a(1,1)*a(2,2) - a(1,2)*a(2,1)
   IF( ABS(det) > 1d-8)THEN
     i1  = (-i(1,3)*a(2,2)+i(2,3)*a(1,2) )/det
     i2  = ( i(1,3)*a(2,1)-i(2,3)*a(1,1) )/det
     i3  = 1d0
     la  = SQRT(i1**2+i2**2+i3**2)
     l(1,k)  = i1/la
     l(2,k)  = i2/la
     l(3,k)  = i3/la
   ELSE
     ! second option, use  the last two equations
     a(1,1) = i(2,2) - lp(k)
     a(1,2) = i(2,3)
     a(2,1) = i(3,2)
     a(2,2) = i(3,3) - lp(k)
     det = a(1,1)*a(2,2) - a(1,2)*a(2,1)
     IF( ABS(det) > 1d-8)THEN
       i1  = 1d0
       i2  = (-i(2,1)*a(2,2)+i(3,1)*a(1,2) )/det
       i3  = ( i(2,1)*a(2,1)-i(3,1)*a(1,1) )/det
       la  = SQRT(i1**2+i2**2+i3**2)
       l(1,k)  = i1/la
       l(2,k)  = i2/la
       l(3,k)  = i3/la
     ELSE
       ! Third option, use  the first and third equations
       a(1,1) = i(1,1) - lp(k)
       a(1,2) = i(1,3)
       a(2,1) = i(3,1)
       a(2,2) = i(3,3) - lp(k)
       det = a(1,1)*a(2,2) - a(1,2)*a(2,1)
       IF( ABS(det) > 1d-8)THEN
         i1  = (-i(1,2)*a(2,2)+i(3,2)*a(1,2) )/det
         i2  = 1d0
         i3  = ( i(1,2)*a(2,1)-i(3,2)*a(1,1) )/det
         la  = SQRT(i1**2+i2**2+i3**2)
         l(1,k)  = i1/la
         l(2,k)  = i2/la
         l(3,k)  = i3/la
       ELSE
         RETURN
       END IF
     END IF
   END IF
 END DO
 ! compute second direction using the vectorial product (det = 1)
 CALL vecpro(l(1,3),l(1,1),l(1,2))
 !Reorder directions to minimize rotations with respecto to input matrix LL
 ! select 1st direction
 pr(1) = DOT_PRODUCT(l(:,1),ll(:,1))     !cos 1-1
 pr(2) = DOT_PRODUCT(l(:,2),ll(:,1))     !cos 2-1
 pr(3) = DOT_PRODUCT(l(:,3),ll(:,1))     !cos 3-1
 lg = MAXLOC( ABS(pr) )           !direction with larger component along original X1 direction
 IF( lg(1) == 2 )THEN               !if it is the second
   lp =  (/  lp(2), lp(3), lp(1) /)                     !shift
   IF( pr(lg(1)) > 0d0 )THEN
     l  = RESHAPE((/ l(:,2),l(:,3),l(:,1) /),(/3,3/))
   ELSE
     l  =-RESHAPE((/ l(:,2),l(:,3),l(:,1) /),(/3,3/))
   END IF
 ELSE IF ( lg(1) == 3 )THEN         !if it is the third
   lp =         (/  lp(3), lp(1), lp(2) /)
   IF( pr(lg(1)) > 0d0 )THEN
     l  = RESHAPE((/ l(:,3),l(:,1),l(:,2) /),(/3,3/))
   ELSE
     l  =-RESHAPE((/ l(:,3),l(:,1),l(:,2) /),(/3,3/))
   END IF
 END IF
 ! select 2nd direction
 pr(2) = DOT_PRODUCT(l(:,2),ll(:,2))  !cos 2-2
 pr(3) = DOT_PRODUCT(l(:,3),ll(:,2))  !cos 3-2
 IF( ABS(pr(3)) > ABS(pr(2)) )THEN    !compare 2nd component of the 2-3 directions
   lp(2:3) =         (/   lp(3),  lp(2) /)          !swap if necessary
   IF( pr(3) > 0d0 )THEN
     l(:,2:3) = RESHAPE((/ l(:,3),-l(:,2) /),(/3,2/))
   ELSE
     l(:,2:3) = RESHAPE((/-l(:,3), l(:,2) /),(/3,2/))
   END IF
 END IF
 ! pass principal values
 i(1,1) = lp(1) + i13
 i(2,2) = lp(2) + i13
 i(3,3) = lp(3) + i13
 ll = l
 RETURN
 END SUBROUTINE princi
