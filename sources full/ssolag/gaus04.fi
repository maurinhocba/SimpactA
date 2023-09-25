 SUBROUTINE gaus04(heade,coord,istop,angdf,locax,isg,psg,ngaus)
 !-------------------------------------------------------------------
 !--- evaluates gauss points constants for solid element 05 PRISM
 !--------------------------------------------------------------------
 IMPLICIT NONE

 INTEGER (kind=4), INTENT(OUT) :: istop,isg(:,:)
 INTEGER (kind=4), INTENT(IN) :: locax,ngaus
 REAL (kind=8), INTENT(IN) :: coord(:,:), & !initial coordinates
                              angdf         !default angle local axes - Ort_Axes
 REAL (kind=8), INTENT(OUT) :: psg(:,:)     !factor for output at Gauss points
 TYPE (ele04), POINTER :: heade       !pointer to first element in the set

 INTEGER (kind=4) ielem,g
 REAL (kind=8) :: ang,factor
 REAL (kind=8) :: x(3,nnode),lb(3,3)
 TYPE (ele04), POINTER :: e           !pointer to an element
 REAL (kind=8) :: thf(ngaus),   & !G-P positions
                  wei(ngaus)      !G-P weigths

 !     ---------------------------------------------------------
 IF( ngaus == 2 )THEN
   isg=RESHAPE((/1,2,1,2/),(/2,2/))
   !psg = RESHAPE((/ 1d0,0d0,0d0,1d0/),(/2,2/)) !SPT must extrapolate to external surfaces
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

 factor = ASIN(1d0)/90d0         !pi/180

 e => heade       !point to first element
 DO               !loop over the elements in the set
   IF (.NOT.ASSOCIATED (e) ) EXIT  !if all element processed, exit

   x(:,1:nnode) = coord(:,e%lnods(1:nnode))     !get original coordinates (element nodes)
   ang = (e%angle+angdf)*factor      ! in radians
   ielem = e%numel
   CALL lcsy04(lb,x,locax,ang)   !compute local shell system
   ! compute element local systems, nodal functions and derivatives
   CALL jacob04(e%cartd,e%dvol,x,ielem,istop, &
                e%cdq,e%nfdas,lb,e%jacin,e%area)

   e%stint = 0d0            !Initializes stresses

   e => e%next          !point to next element
 END DO

 RETURN
 END SUBROUTINE gaus04
