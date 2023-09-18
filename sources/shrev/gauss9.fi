 SUBROUTINE gauss9(ndime,ntype,nstre,nnode,ngaus,axesc,     &
                   coord,eule0,istop,                       &
                   heade,gauss,posgp,weigh,shape,deriv)
 !********************************************************************
 !
 !*** evaluates gauss points constants for shell/beam element
 !
 !********************************************************************
 IMPLICIT NONE

 INTEGER (kind=4), INTENT(IN) :: ndime,ntype,nstre,nnode,ngaus,axesc
 INTEGER (kind=4), INTENT(OUT) :: istop
 REAL (kind=8), INTENT(IN) :: coord(:,:), & !initial coordinates
                              eule0(:,:)    !original nodal system
 TYPE (ele09), POINTER :: heade   !pointer to first element in the set
 LOGICAL, INTENT (IN OUT) :: gauss    !TRUE if initial constant exists
 REAL (kind=8), INTENT(OUT) :: posgp(:), & !Gauss points positions
                               weigh(:), & !Gauss points weights
                               shape(:,:), & !shape functions at Gauss points
                               deriv(:,:)    !shape functions derivatives at Gauss points

 INTEGER (kind=4) g,ielem
 REAL    (kind=8) x(ndime,nnode),phin(nnode)
 TYPE (ele09), POINTER :: e
 REAL (kind=8),PARAMETER :: gaus1 = 0.577350269189626,   & !1/SQRT(3)
                            factor = 1.74532925199433d-02  !pi/180

 INTERFACE
   INCLUDE 'shape9.h'
 END INTERFACE

 !       gauss points in local coordinates and weigths
 ! gauss points shape and derivatives of nodal functions
 CALL shape9(weigh,posgp,shape,deriv,nnode,ngaus)

 IF( .NOT.gauss )THEN
 e => heade       !point to first element
 DO               !loop over the elements in the set
   IF (.NOT.ASSOCIATED (e) ) EXIT  !if all element processed, exit
   ! gathers the initial coordinates and directors
   x  = coord(:,e%lnods)
   phin = eule0(1,e%lnods)
   ielem = e%numel
   CALL setga9(ndime,nstre,ntype,nnode,ngaus,axesc,shape,deriv,x,   &
               phin,e%phil,e%jac,e%r0,e%stra0,e%sede0,ielem)
   DO g=1,ngaus
     IF(e%jac(g) <= 0) THEN
       WRITE(lures,"(' Negative or null Jacobian at beam/shell element (2) ',      &
            &    i5, /' at Gauss Point ',i2,/,' Jacobian = ', e12.4 )",ERR=9999)   &
            &    e%numel,g,e%jac(g)
       istop = 1
     END IF
   END DO
   e => e%next
 END DO
 gauss = .TRUE.
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE gauss9
