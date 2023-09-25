 SUBROUTINE gaus05(heade,coord,istop,angdf,locax,nnode,quad,isg,psg,ngaus)
 !-------------------------------------------------------------------
 !--- evaluates gauss points constants for solid element 05 PRISM
 !--------------------------------------------------------------------
 IMPLICIT NONE

 INTEGER (kind=4), INTENT(OUT) :: istop,isg(:,:)
 INTEGER (kind=4), INTENT(IN) :: locax, &
                                 ngaus,nnode
 LOGICAL, INTENT(IN) :: quad
 REAL (kind=8), INTENT(IN) :: coord(:,:), & !initial coordinates
                              angdf         !default angle local axes - Ort_Axes
 REAL (kind=8), INTENT(OUT) :: psg(:,:)     !factor for output at Gauss points
 TYPE (ele05), POINTER :: heade       !pointer to first element in the set

 INTEGER (kind=4) i,ielem
 REAL (kind=8) :: shape(nnb,ngaud),deriv(nnb,3,ngaud)
 REAL (kind=8) :: ang,factor
 REAL (kind=8) :: x(3,nnode),lb(3,3)
 REAL (kind=8), PARAMETER :: a=0.333333333333333d0            !1/3
 TYPE (ele05), POINTER :: e           !pointer to an element
 LOGICAL :: sides(3)                  !side element existence for quad approach
 REAL (kind=8) :: thf(ngaus),   & !G-P positions
                  wei(ngaus)      !G-P weigths

 !     ---------------------------------------------------------
 IF( ngaus == 2 )THEN
   isg=RESHAPE((/1,2,1,2/),(/2,2/))
   !psg = RESHAPE((/ 1d0,0d0,0d0,1d0/),(/2,2/)) !STP must extrapolate to external surfaces
   psg = RESHAPE((/ 1.366025403784438,-0.366025403784438, &  !extrapolates directly (STP must do nothing)
                   -0.366025403784438, 1.366025403784438/),(/2,2/))
 ELSE
   CALL gaussq(ngaus,thf(1),wei(1)) !integration points through the thickness
   !extrapolates directly (STP must do nothing)
   isg = RESHAPE((/ 1,2, ngaus-1,ngaus /),(/2,2/))
   factor = 1.d0+(1d0-thf(ngaus))/(thf(ngaus)-thf(ngaus-1))
   psg = RESHAPE((/ factor, 1d0-factor, 1d0-factor, factor /),(/2,2/))
   WRITE(58,"(4i12,/,4e12.4)")isg,psg
 END IF
 !     ---------------------------------------------------------

 CALL shape4(nnb,shape(1,1),deriv(1,1,1),a,a,-1d0)   ! bottom face
 CALL shape4(nnb,shape(1,2),deriv(1,1,2),a,a, 0d0)   ! element center
 CALL shape4(nnb,shape(1,3),deriv(1,1,3),a,a,+1d0)   ! top face

 !     ---------------------------------------------------------

 factor = ASIN(1d0)/90d0         !pi/180

 e => heade       !point to first element
 DO               !loop over the elements in the set
   IF (.NOT.ASSOCIATED (e) ) EXIT  !if all element processed, exit

   x(:,1:nnb) = coord(:,e%lnods(1:nnb))     !get original coordinates (element nodes)
   IF( quad ) THEN
     DO i=1,3        !get original coordinates (extra patch nodes)
       sides(i) = ( e%lnods(i+6) /= 0 )         !side element exist
       IF( sides(i) )THEN
         x(:,i+6:i+9:3) = coord(:,e%lnods(i+6:i+9:3))
       ELSE
         x(:,i+6:i+9:3) = 0d0                   !may be unnecessary
       END IF
     END DO
   END IF
   ang = (e%angle+angdf)*factor      ! in radians
   ielem = e%numel
   CALL lcsy05(lb,x,locax,ang)   !compute local shell system
   ! compute element local systems, nodal functions and derivatives
   CALL jacob05(deriv,e%cartd,e%dvol,x,ielem,istop, &
                e%cdq,lb,sides,e%nfdas,e%jacin,quad)

   e%stint = 0d0            !Initializes stresses

   e => e%next          !point to next element
 END DO

 RETURN
 END SUBROUTINE gaus05
