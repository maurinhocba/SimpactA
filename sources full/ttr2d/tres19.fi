 SUBROUTINE tres19( heade, coora, tresi)
 !-------------------------------------------------------------------
 !
 !          EVALUATION OF HEAT RESIDUAL for LINEAR 2D TRIANGLE
 !
 !-------------------------------------------------------------------
 IMPLICIT NONE

 ! dummy arguments
 REAL(kind=8),INTENT (IN) :: coora(:,:)       !present coordinates
 REAL(kind=8),INTENT (IN OUT) :: tresi(:)   !internal nodal forces
 TYPE (ele19), POINTER  :: heade              !pointer to first element

 ! local variables
 LOGICAL :: newmt    ! TRUE if material constant computation necessary

 INTEGER (kind=4) ielem,   & !element number
                  isec,    & !associated material
                  i,j,n,l, & !different indexes
                  ic,      & !sector of curve for properties interpolation 
                  mtype,   & !associated material type
                  ncurv(4),& !associated curves for therm dependant properties
                  numct,   & !number of curves in therm depend problem
                  osec       !associated material of previous element

 REAL (kind=8) x(2,3),      & !nodal coordinates of the patch
               young,       & !young elastic modulus
               nu,          & !poisson ratio
               alpha,       & !Thermal dilatation parameter
               bulk,        & !Bulk modulus
               condu,       & !conductivity
               a(3),b(3),   & !side projections
               t1(2),t2(2), & !sides 1 & 2
               area2,       & !present area (twice)
               t(3),        & !nodal temperatures
               dt(3),       & !nodal increment in temperatures
               rl(3),       & !nodal internal equivalent forces
               r1(3),r2(3), & !auxiliar for internal forces
               k11,k12,k13,k22,k23,k33,k,f,ta,td

 TYPE (ele19), POINTER :: e    !pointer to an element data
 TYPE (section), POINTER :: sec    !pointer to a section data
 TYPE (mater), POINTER :: mat    !pointer to a section data
 TYPE (curve), POINTER :: crvy,crvn,crva,crvc

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
     !thick = sec%rprop(1)           !thickness
     young = mat%prope(1)           !young modulus
     nu    = mat%prope(2)           !nu
     bulk  = mat%prope(4)           !Bulk modulus
     alpha = mat%prope(6)           !Thermical Dilatation Coeff
     condu = mat%prope(8)           !conductivity
     !associated curves for thermal dependant properties
     ncurv(1) = INT(sec%mtbas%props(1))     !young modulus
     ncurv(2) = INT(sec%mtbas%props(2))     !poisson ratio
     ncurv(3) = INT(sec%mtbas%props(6))     !thermal dilatation coeff
     ncurv(4) = INT(sec%mtbas%props(8))     !heat conductiv
     CALL cur_point(ncurv(1),mat%chead,crvy) !pointer to young curve
     CALL cur_point(ncurv(2),mat%chead,crvn) !pointer to poisson curve
     CALL cur_point(ncurv(3),mat%chead,crva) !pointer to alpha curve
     CALL cur_point(ncurv(4),mat%chead,crvc) !pointer to heat conduc. curve
   END IF
   numct = SUM(ncurv)

   !*** computes actual configuration
   x(:,1:3)  = coora(:,e%lnods(1:3)) ! get coordinates of the main triangle
   ! get and transform temperatures in to hierachic DOFs
   dt = 0d0
   DO i=1,3
     n     = e%lnods(i)
     t(i)  = tempe(1,n)  !temp in nodes
     dt(i) = dtemp(1,n)  !temp increment in nodes
   END DO

   !  for temperature dependency of thermal parameters
   IF( numct > 0 )THEN !compute properties as a function of temperature
     f = SUM(t(:))/3d0 ! temperature at element center
     ic = 1
     IF(ncurv(1) > 0) young = inte_cr(crvy%val,crvy%np,f,ic)
     IF(ncurv(2) > 0) nu    = inte_cr(crvn%val,crvn%np,f,ic)
     IF(ncurv(3) > 0) alpha = inte_cr(crva%val,crva%np,f,ic)
     IF(ncurv(4) > 0) condu = inte_cr(crvc%val,crvc%np,f,ic)
     IF(ANY(ncurv(1:2) > 0)) bulk = young/(1d0-2d0*nu)/3d0 !bulk modulus
   END IF

   !*** evaluate the side projections
   !derivatives of the in-plane mapping for a triangular element
   t1 = x(:,3) - x(:,2)  !first  side 2-3
   t2 = x(:,1) - x(:,3)  !second side 3-1
   area2 = t1(1)*t2(2)-t2(1)*t1(2) !2A
   DO i=1,3
     j = kk(1,i) ! secnd node for side i
     k = kk(2,i) ! first node for side i
     a(i) = -x(2,j)+x(2,k) ! dN(i)/dx = ( -y(k) + y(j) ) / 2A
     b(i) =  x(1,j)-x(1,k) ! dN(i)/dy = (  x(k) - x(j) ) / 2A
   END DO
   a = a/area2  ! cartesian derivatives dN/dx
   b = b/area2  ! cartesian derivatives dN/dy

   ! residual heat due to conductivity
   ! (a(i)*a(j) + b(i)*b(j) )/(2A)^2
   k11 = a(1)*a(1) + b(1)*b(1)
   k12 = a(1)*a(2) + b(1)*b(2)
   k13 = a(1)*a(3) + b(1)*b(3)
   k22 = a(2)*a(2) + b(2)*b(2)
   k23 = a(2)*a(3) + b(2)*b(3)
   k33 = a(3)*a(3) + b(3)*b(3)
   ! solid temperatures
   rl(1) = k11*t(1) + k12*t(2) + k13*t(3)
   rl(2) = k12*t(1) + k22*t(2) + k23*t(3)
   rl(3) = k13*t(1) + k23*t(2) + k33*t(3)
   k = condu*e%area1          !e%area1 is volume of element 19
   rl = k*rl                  !conduction contribution

   !*** next contributions are associated with changes from previous configuration
   ! plastic work contribution
   k = e%pwork(1)*e%area1/3d0/td
   DO i=1,3
     rl(i) = rl(i) - k
   END DO

   ! elastic-volumetric change contribution (Gough-Joule effect)
   k  = bulk*alpha*e%pwork(2)*e%area1/12d0/td
   r1 = tref + t   !absolute temperature
   ta = SUM(t1)    !average
   r2 = ta + r1    !consistent aprox matrix
   rl(1) = rl(1) + k*r2(1)
   rl(2) = rl(2) + k*r2(2)
   rl(3) = rl(3) + k*r2(3)

   ! change in plastic work due to change in temperature
   k = e%pwork(3)*e%area1/12d0/td
   ta = SUM(dt)   !average
   r1 = ta + dt   !consistent aprox matrix
   rl(1) = rl(1) - k*(r1(1))
   rl(2) = rl(2) - k*(r1(2))
   rl(3) = rl(3) - k*(r1(3))

   DO i=1,3             !for each node in the element
     n = e%lnods(i)     !node
     l = iftmp(1,n)
     IF( l > 0 ) tresi(l) = tresi(l)+ rl(i) !sum internal force
   END DO

   e => e%next          !point to next element in the list
 END DO

 RETURN
 END SUBROUTINE tres19
