 SUBROUTINE bfle11(ntype,fixit,l0,r0,lb,t,n,x,cab,sh,b,g)
   ! computes b-matrix for standard element
   ! 2-d beam/shell rotation-free elements
   IMPLICIT NONE
   INTEGER (kind=4), INTENT(IN) :: ntype,    & !problem type
                                   fixit(2), & !side type
                                   g           !gauss point
   REAL (kind=8), INTENT(IN) :: l0(3),  & !original element sizes
                                r0(2),     & !original X1 and curvature
                                lb(3),  & !stretchings
                                t(2),   & !tangent vector
                                n(2),   & !normal vector
                                x(2,4), & !element coordinates
                                cab(4,2),&! a0,1/li,gamma
                                sh(2)     ! shape function

   REAL (kind=8), INTENT(OUT) :: b(2,4,2) !B matrix
   REAL (kind=8) :: bn(2,4,2),  &  !nodal B matrix
                    l              !length at present element

   ! local variables
   REAL (kind=8) :: fac,fac1,fac2,fac3,fac4,ni(2),ti(2)

   !IF( g == 1)THEN
     bn   = 0d0       !initializes
     l    = l0(2)*lb(1)  !length of present element
     !              CURVATURE IN THE PLANE
     fac = 1d0/l0(2)     !inverse length
     ! curvature at first node
     SELECT CASE (fixit(1))
     CASE (0)  !sides exists
       ! compute ti, ni
       ti = (x(:,2) - x(:,1))*cab(2,1) !tangent at side element
       ni = (/ -ti(2), ti(1) /)        !normal at side element
       ! compute factors
       fac2 = fac*fac*cab(3,1)         !gamma/l0**2
       fac1 = fac*cab(4,1)*lb(1)       !Ci*ls/l0
       fac3 = fac1*cab(2,1)            !        /li
       fac4 = fac1/l                   !        /l
       bn(1:2,1,1) =        -fac3*ni
       bn(1:2,2,1) = -fac2*t+fac3*ni+fac4*n
       bn(1:2,3,1) = +fac2*t        -fac4*n
     CASE (1)  !first side simple supported
       !Nothing
     CASE (2)  !first side clamped
       fac1 = fac*cab(2,1)               !2/l0**2       2/l0 = cab(2,1)
       fac2 = fac*cab(3,1)*fac           !gam/l0**2
       bn(1:2,2,1) =  fac1*n-fac2*t
       bn(1:2,3,1) = -bn(1:2,2,1)
     CASE (3)  !first side branching (first part)
       fac1 = -cab(3,1)*fac*fac          !-Gam/l0**2
       fac2 = cab(2,1)*fac               !2/l0**2
       bn(1:2,2,1) =  fac1*t+fac2*n
       bn(1:2,3,1) = -bn(1:2,2,1)
     END SELECT

     ! curvature at second node
     SELECT CASE (fixit(2))
     CASE (0)  !side exists
       ! compute td, nd
       ti = (x(:,4) - x(:,3))*cab(2,2) !tangent at side element
       ni = (/ -ti(2), ti(1) /)        !normal at side element
       ! compute factors
       fac1 = fac*fac*cab(3,2)         !gamma/l0**2
       fac2 = fac*cab(4,2)*lb(1)       !Cd*ls/l0
       fac3 = fac2/l                   !        /l
       fac4 = fac2*cab(2,2)            !        /ld
       bn(1:2,2,2) = -fac1*t-fac3*n
       bn(1:2,3,2) =  fac1*t+fac3*n+fac4*ni
       bn(1:2,4,2) =               -fac4*ni
     CASE ( 1)  !second side simple supported
       !nothing
     CASE (2)  !second side clamped
       fac1 =  fac*cab(3,2)*fac        !gam/l0**2
       fac2 = -fac*cab(2,2)            !-2/(l0**2)    2/l0 = cab(2,2)
       bn(1:2,2,2) = fac2*n-fac1*t
       bn(1:2,3,2) = -bn(1:2,2,2)
     CASE (3)  !second side branching (first part)
       fac1 = -cab(3,2)*fac*fac          !-gam/l0**2
       fac2 = -fac*cab(2,2)              !-2/l0**2
       bn(1:2,2,2) =  fac1*t+fac2*n
       bn(1:2,3,2) = -bn(1:2,2,2)
     END SELECT

   !END IF
   IF( g > 0 )THEN
     ! Bb matrix at integration point
     b(:,:,1) = sh(1)*bn(:,:,1) + sh(2)*bn(:,:,2)

     !               2ND CURVATURE
     IF( ntype == 3 )THEN        !for axil-symmetryc problems
       fac = 1d0/r0(1)/r0(1)                   !-1/r^2
       b(1,2:3,2) = -fac*(t(2)-r0(2))*sh       !t(2)-r0(2) = sin(f) - sin(f^0)
       b(2,2:3,2) = 0d0
       fac = fac*(sh(1)*x(1,2)+sh(2)*x(1,3))*t(1)/l  !t(1) = cos(f)
       b(:,2,2) = b(:,2,2) + fac*n
       b(:,3,2) = b(:,3,2) - fac*n
     END IF
   ELSE
     ! Stabilization Bb matrix
     b(:,:,1) = -bn(:,:,1) + bn(:,:,2)
   END IF

   RETURN
 END SUBROUTINE bfle11
