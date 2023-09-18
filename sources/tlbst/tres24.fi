 SUBROUTINE tres24( heade, coora, tresi)
 !-------------------------------------------------------------------
 !
 !          EVALUATION OF HEAT RESIDUAL for BST++
 !           !
 !-------------------------------------------------------------------
 IMPLICIT NONE

 ! dummy arguments
 REAL(kind=8),INTENT (IN) :: coora(:,:)       !present coordinates

 REAL(kind=8),INTENT (IN OUT) :: tresi(:)   !internal nodal forces
 TYPE (ele24), POINTER  :: heade              !pointer to first element

 ! local variables
 INTEGER (kind=4), PARAMETER :: nn = 1000000

 LOGICAL :: newmt    ! TRUE if material constant computation necessary

 INTEGER (kind=4) ielem, & !element number
                  isec,  & !associated material
                  i,j,n,l, & !different indexes
                  mtype, & !associated material type
                  na,nk, & !number of points in curve
                  osec     !associated material of previous element

 REAL (kind=8) x(3,3),    & !nodal coordinates of the patch
               thick,     & !thickness (original)
               thnew,     & !present thickness
               alpha,     & !Thermal dilatation parameter
               bulk,      & !Bulk modulus
               condu,     & !conductivity
               a(3),b(3), & !side projections
               t1(3),t2(3), & !sides
               area2,     & !present area (twice)
               t(3,3),    & !nodal temperatures
               rl(3,3),   & !nodal internal equivalent forces
               k11,k12,k13,k22,k23,k33,k,f,ta,tb,tc,dt,td

 TYPE (ele24), POINTER :: e    !pointer to an element data
 TYPE (section), POINTER :: sec    !pointer to a section data
 TYPE (mater), POINTER :: mat    !pointer to a section data
 TYPE (curve), POINTER :: ca,ck

 !***  Computes nodal internal forces, integrating trought the thickness

 osec = -1                           !no previous material
 b(1) = 0d0                          !always
 td = tdtime
 IF( td == 0 ) td = 1d0              !to avoid an error due to initial plastic work
 e => heade                          !point to first element in the list
 DO                                  !loop over all the elements
   IF (.NOT.ASSOCIATED (e) ) EXIT    !last element processed => EXIT
   ielem = e%numel                   !element number
   isec = e%matno                    !associated material
   newmt = osec /= isec              !same material than previous ?
   IF(newmt)THEN                     !if new material, get properties
     sec => psecs(isec)%p
     mat => sec%mtbas
     osec = isec                    !keep present material
     mtype = mat%mtype              !type of base material
     thick = sec%rprop(1)           !original thickness
     bulk  = mat%prope(4)           !Bulk modulus
     alpha = mat%prope(6)           !Thermical Dilatation Coeff
     na    = INT(sec%mtbas%props(6))        !associated curve for thermal dilatation
     CALL cur_point( na,sec%mtbas%chead,ca) !point to the curve
     condu = mat%prope(8)           !conductivity
     nk    = INT(sec%mtbas%props(8))        !associated curve for Heat conductivity
     CALL cur_point( nk,sec%mtbas%chead,ck) !point to the curve
   END IF

   !*** computes actual configuration

   x(:,1:3)  = coora(:,e%lnods(1:3)) ! get coordinates of the main triangle
   ! get and transform temperatures in to hierachic DOFs
   dt = 0
   DO i=1,3
     n = e%lnods(i)
     SELECT CASE (ndoft)
     CASE (1)
       t(i,1) = tempe(1,n)        ! T0
       dt = dt+ dtemp(1,n)        ! dT0
     CASE (2)
       t(i,1) = (+tempe(1,n)+ tempe(2,n))/2d0    ! T0
       t(i,2) = (-tempe(1,n)+ tempe(2,n))/2d0    ! T1
       dt = dt+ (+dtemp(1,n)+ dtemp(2,n))/2d0    ! dT0
     CASE (3)
       t(i,1) =   tempe(1,n)                                 ! T0
       t(i,2) = (-tempe(2,n)+ tempe(3,n))/2d0                ! T1
       t(i,3) =  -tempe(1,n)+(tempe(2,n)+ tempe(3,n))/2d0    ! T2
       dt = dt+   dtemp(1,n)                                 ! dT0
     END SELECT
   END DO
   !  for temperature dependency of thermal parameters
   IF( na > 0 .OR. nk > 0 )THEN
     ! temperature at element center (middle surface)
     f = SUM(t(:,1))/3d0
     i = 1
     !compute thermal  parameters as a function of temperature
     IF( na > 0 ) alpha = inte_cr(ca%val,na,f,i)
     IF( nk > 0 ) condu = inte_cr(ck%val,nk,f,i)
   END IF
   !*** evaluate the side projections (derivatives of the in-plane mapping)

   t1 = x(:,3) - x(:,2)  !first side 2-3
   t2 = x(:,1) - x(:,3)  !second side 3-1
   a(1) = SQRT(DOT_PRODUCT(t1,t1))         !first side length (base length)
   a(2) = DOT_PRODUCT(t2,t1)/a(1)          !proyection of second side over first
   !b(1) = 0d0                              !no proyection
   b(2) = SQRT(DOT_PRODUCT(t2,t2)-a(2)**2) !triangle height
   a(3) = -a(1)-a(2)
   b(3) = -b(2) !-b(1)
   area2 = a(1)*b(2)                       !twice the area
   a = a/area2                             !cartesian derivatives of mapping (2)
   b = b/area2                             !-cartesian derivatives of mapping (1)
   thnew = thick*e%lb                  !thickness at new configuration

   ! residual heat due to conductivity
   ! (a(i)*a(j) + b(i)*b(j) )/(2A)^2
   k11 = a(1)*a(1) !+ b(1)*b(1)   !
   k12 = a(1)*a(2) !+ b(1)*b(2)   !
   k13 = a(1)*a(3) !+ b(1)*b(3)   !
   k22 = a(2)*a(2) + b(2)*b(2)
   k23 = a(2)*a(3) + b(2)*b(3)
   k33 = a(3)*a(3) + b(3)*b(3)
   ! middle surface temperatures
   rl(1,1) = k11*t(1,1) + k12*t(2,1) + k13*t(3,1)
   rl(2,1) = k12*t(1,1) + k22*t(2,1) + k23*t(3,1)
   rl(3,1) = k13*t(1,1) + k23*t(2,1) + k33*t(3,1)
   ta = SUM(t(:,1))     ! sum of temperatures at the middle surface
   IF( ndoft > 1 ) THEN         !first hierarchic DOF
     f = 1d0/(3d0*thnew*thnew)         !1/(3t^2)
     tb = SUM(t(:,2))     ! sum of first H-Dof
     rl(1,2) = ( k11*t(1,2) + k12*t(2,2) + k13*t(3,2) )/3d0 + f*(t(1,2)+tb)
     rl(2,2) = ( k12*t(1,2) + k22*t(2,2) + k23*t(3,2) )/3d0 + f*(t(2,2)+tb)
     rl(3,2) = ( k13*t(1,2) + k23*t(2,2) + k33*t(3,2) )/3d0 + f*(t(3,2)+tb)
   END IF
   IF( ndoft == 3 ) THEN
     rl(1,1) = rl(1,1)+ (k11*t(1,3) + k12*t(2,3) + k13*t(3,3))/3d0
     rl(2,1) = rl(2,1)+ (k12*t(1,3) + k22*t(2,3) + k23*t(3,3))/3d0
     rl(3,1) = rl(3,1)+ (k13*t(1,3) + k23*t(2,3) + k33*t(3,3))/3d0
     tc = SUM(t(:,3))     ! sum of second H-Dof
     f = f*4d0/3d0
     rl(1,3) = ( k11*t(1,1) + k12*t(2,1) + k13*t(3,1) )/3d0 &
             + ( k11*t(1,3) + k12*t(2,3) + k13*t(3,3) )/5d0 + f*(t(1,3)+tc)
     rl(2,3) = ( k12*t(1,1) + k22*t(2,1) + k23*t(3,1) )/3d0 &
             + ( k12*t(1,3) + k22*t(2,3) + k23*t(3,3) )/5d0 + f*(t(2,3)+tc)
     rl(3,3) = ( k13*t(1,1) + k23*t(2,1) + k33*t(3,1) )/3d0 &
             + ( k13*t(1,3) + k23*t(2,3) + k33*t(3,3) )/5d0 + f*(t(3,3)+tc)
   END IF
   k = condu*area2*thnew/2d0          !2 because is AREA2
   rl = k*rl                          !conduction contribution

   !the next two contributions are associated with changes
   !from previous configuration.
   !(beware TDTIME must be greater than 0)

   ! plastic work contribution
   f = e%area1*thick/3d0
   k = f*e%pwork(1)/td
   rl(:,1) = rl(:,1) - k
   IF( ndoft > 1 )THEN
     k = 2d0*f*e%pwork(2)/thnew/td
     rl(:,2) = rl(:,2) - k
   END IF
   IF( ndoft == 3 )THEN
     k = 4d0*f*e%pwork(3)/thnew/thnew/td
     rl(:,3) = rl(:,3) - k
   END IF

   ! volumetric change contribution (use temperature increment at the element center)
   k = bulk*alpha*e%area1*thick/12d0/(1d0+alpha*ta)*alpha*dt /td
   ta = ta + 4d0*tref
   rl(1,1) = rl(1,1) + k*(ta+t(1,1))
   rl(2,1) = rl(2,1) + k*(ta+t(2,1))
   rl(3,1) = rl(3,1) + k*(ta+t(3,1))
   IF( ndoft > 1 )THEN
     k = k/3d0
     rl(1,2) = rl(1,2) + k*(tb+t(1,2))
     rl(2,2) = rl(2,2) + k*(tb+t(2,2))
     rl(3,2) = rl(3,2) + k*(tb+t(3,2))
     IF( ndoft == 3 )THEN
       f = k/5d0*3d0
       rl(1,1) = rl(1,1) + k*(tc+t(1,3))
       rl(2,1) = rl(2,1) + k*(tc+t(2,3))
       rl(3,1) = rl(3,1) + k*(tc+t(3,3))
       rl(1,3) = rl(1,3) + k*(ta+t(1,1)) + f*(tc+t(1,3))
       rl(2,3) = rl(2,3) + k*(ta+t(2,1)) + f*(tc+t(2,3))
       rl(3,3) = rl(3,3) + k*(ta+t(3,1)) + f*(tc+t(3,3))
     END IF
   END IF

   DO i=1,3             !for each node in the element
     n = e%lnods(i)     !node
     SELECT CASE (ndoft)
     CASE (1)
       t1(1) = rl(i,1)
     CASE (2)
       t1(1) = (rl(i,1)-rl(i,2))/2d0
       t1(2) = (rl(i,1)+rl(i,2))/2d0
     CASE (3)
       t1(1) =  rl(i,1)-rl(i,3)
       t1(2) =(-rl(i,2)+rl(i,3))/2d0
       t1(3) =( rl(i,2)+rl(i,3))/2d0
     END SELECT
     DO j=1,ndoft
       l = iftmp(j,n)
       IF( l > 0 ) tresi(l) = tresi(l)+ t1(j) !sum internal force
     END DO
   END DO
   e => e%next          !point to next node in the list
 END DO

 RETURN
 END SUBROUTINE tres24
