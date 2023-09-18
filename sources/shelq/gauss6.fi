 SUBROUTINE gauss6(ndime,coord,eule0,istop, &
                   heade,gauss,angdf,locax,zigzag)
 !********************************************************************
 !
 !*** evaluates gauss points constants for shell element 6
 !    valid only for 4-noded cuadrilateral element
 !********************************************************************
 IMPLICIT NONE

 INTEGER (kind=4), INTENT(IN) :: ndime,locax
 INTEGER (kind=4), INTENT(OUT) :: istop
 REAL (kind=8), INTENT(IN) :: coord(:,:), & !initial coordinates
&                             eule0(:,:), & !original nodal system
&                             angdf         !default angle X_1 - Ort_1
 TYPE (ele06), POINTER :: heade   !pointer to first element in the set
 LOGICAL, INTENT (IN) :: zigzag       !TRUE if additional in-plane displacement
 LOGICAL, INTENT (IN OUT) :: gauss    !TRUE if initial constant exists

 INTEGER (kind=4) g,i
 REAL    (kind=8) deriv(nnode,2,ngaus),weigp(ngaus),ang,     &
&                 x(ndime,nnode),t(ndime,nnode),rm(9),lb(3,2,nnode)
 TYPE (ele06), POINTER :: e
 REAL (kind=8),PARAMETER :: factor = 1.74532925199433d-02  !pi/180

 INTERFACE
   INCLUDE 'inrotm.h'
 END INTERFACE

 !    gauss points  weigths

 weigp = 1d0

 ! gauss points derivatives of nodal functions

 DO g=1,ngaus
   deriv(1,1,g) = (posgp(2,g)-1d0)/4d0
   deriv(2,1,g) = -deriv(1,1,g)
   deriv(3,1,g) =  (posgp(2,g)+1d0)/4d0
   deriv(4,1,g) = -deriv(3,1,g)
   deriv(1,2,g) =  (posgp(1,g)-1d0)/4d0
   deriv(2,2,g) = -(posgp(1,g)+1d0)/4d0
   deriv(3,2,g) = -deriv(2,2,g)
   deriv(4,2,g) = -deriv(1,2,g)
 END DO

 IF( .NOT.gauss )THEN
 e => heade       !point to first element
 DO               !loop over the elements in the set
   IF (.NOT.ASSOCIATED (e) ) EXIT  !if all element processed, exit
   ! gathers the initial coordinates and directors
   x = coord(:,e%lnods)
   DO i=1,nnode
     CALL inrotm(eule0(:,e%lnods(i)),rm(:))
     t(:,i) = rm(7:9)
     IF( zigzag)THEN
        lb(:,1,i) =  rm(1:3)
        lb(:,2,i) =  rm(4:6)
     END IF
   END DO
   ang = (e%angle+angdf)*factor
   CALL setga6(nnode,ngaus,ngamm,    x,    t,          &
             & shape,deriv,weigp,e%dvolu,              &
             & e%cartd,e%jacin,e%stra0,e%gamm0,e%ambda,ang,locax)
   DO g=1,ngaus
     IF(e%dvolu(g) <= 0) THEN
       WRITE(lures,"(' Negative or null Jacobian at shell elmnt',         &
       & i8,/,' at Gauss point ',i2,/,' Jacobian =',e13.4)",ERR=9999)     &
       & e%numel,g,e%dvolu(g)
       istop = 1
     END IF
   END DO
   IF( zigzag ) CALL comp_ang06(nnode,x,t,lb,e%nangl,e%jas,ang,locax)
   e => e%next
 END DO
 gauss = .TRUE.
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE gauss6
