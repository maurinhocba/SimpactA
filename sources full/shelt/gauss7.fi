 SUBROUTINE gauss7(stype,coord,eule0,istop,ngaus,nnode,ngamm, &
                   heade,gauss,angdf,locax,posgp,shape,weigp,zigzag)
 !********************************************************************
 !
 !*** evaluates gauss points constants for shell element 7
 !    valid only for 6-noded triangular  element
 !********************************************************************
 IMPLICIT NONE

 INTEGER (kind=4), INTENT(IN) :: stype,locax,ngaus,nnode,ngamm
 INTEGER (kind=4), INTENT(OUT) :: istop
 REAL (kind=8), POINTER    :: coord(:,:), & !initial coordinates
                              eule0(:,:), & !original nodal system
                              posgp(:,:), & !gauss point local coord
                              weigp(:)      !gauss point local coord
 REAL (kind=8), INTENT(IN) :: angdf         !default angle X_1 - Ort_1
 REAL (kind=8), POINTER :: shape(:,:)    !gauss point shape functions
 TYPE (ele07), POINTER :: heade   !pointer to first element in the set
 LOGICAL, INTENT (IN) :: zigzag       !TRUE if additional in-plane displacement
 LOGICAL, INTENT (IN OUT) :: gauss    !TRUE if initial constant exists

 INTEGER (kind=4) i,ni,m,g,typep
 REAL    (kind=8) deriv(6,2,ngaus),ang,x(3,6),t(3,6),rm(9),lb(3,2,6)
 LOGICAL :: mems(3),quad
 TYPE (ele07), POINTER :: e
 REAL (kind=8),PARAMETER ::  factor = 1.74532925199433d-02  !pi/180

 INTERFACE
   INCLUDE 'inrotm.h'
 END INTERFACE

  ! gauss points derivatives of nodal functions
 IF( ngaus /= 1)THEN !TQQQ
   DO g=1,ngaus
     CALL shape7(posgp(1,g),posgp(2,g),shape(1,g),deriv(1,1,g))
   END DO
 END IF

 IF( .NOT.gauss )THEN
 typep = 4-stype
 e => heade       !point to first element
 DO               !loop over the elements in the set
   IF (.NOT.ASSOCIATED (e) ) EXIT  !if all element processed, exit
   IF( ngaus == 1 )THEN          !TLLL
     ! gathers the initial coordinates and directors
     x(:,1:3) = coord(:,e%lnods(1:3))
     DO i=1,3
       CALL inrotm(eule0(:,e%lnods(i+3)),rm(:))
       t(:,i) = rm(7:9)
       IF( zigzag)THEN
          CALL inrotm(eule0(:,e%lnods(i)),rm(:))
          lb(:,1,i) =  rm(1:3)
          lb(:,2,i) =  rm(4:6)
       END IF
     END DO
     ang = (e%angle+angdf)*factor
     quad = nnode == 9
     IF(quad)THEN
       DO m =1,3                       !for the side nodes in the patch
         i = m+3                       !associated node
         ni = e%lnods(i+3)               !node number
         IF(ni > 0) THEN
           x(1:3,i) = coord(1:3,ni)    !get original coordinates of adjacent node
           mems(m) = .TRUE.               !side element exists
         ELSE
           mems(m) = .FALSE.              !side element does not exist
           x(1:3,i) = 0d0            !null values (unnecessary perhaps)
         END IF
       END DO
     END IF
     CALL setg07(x,t,e%dvolu(1),e%cartd(:,:,1),e%jacin(:,:,1),e%stra0(:,1),e%gamm0, &
                 e%ambda(:,1),ang,locax,quad,mems,e%cd)
     IF( ASSOCIATED(e%nangl)) CALL comp_ang07(x,t,lb,e%nangl,ang,locax,e%jas)
   ELSE                     !TQQQ
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
     CALL setga7(stype,dn(:,:,typep:4),x,t,deriv,e%dvolu,e%cartd,e%jacin, &
                 e%stra0,e%gamm0,e%ambda,ang,locax,shape,weigp,ngaus,nnode,ngamm)
     IF( ASSOCIATED(e%nangl)) CALL comp_angl7(nnode,x,t,lb,e%nangl,ang,locax,e%jas)
   END IF
   DO g=1,ngaus
     IF(e%dvolu(g) <= 0) THEN
       WRITE(lures,"(' Negative or null Jacobian at shell elmnt',         &
       & i8,/,' at Gauss point ',i2,/,' Jacobian =',e13.4)",ERR=9999)     &
       & e%numel,g,e%dvolu(g)
       istop = 1
     END IF
   END DO
   e => e%next
 END DO
 gauss = .TRUE.
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE gauss7
