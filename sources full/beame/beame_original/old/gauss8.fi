 SUBROUTINE gauss8(ndime,nelem,nnode,ngaus,axesc,heade,coord,     &
                   eule0,posgp,shapf,deriv,weigh,flag,istop)
 !********************************************************************
 !
 !*** evaluates gauss points constants for beam element
 !
 !********************************************************************
 IMPLICIT NONE

 INTEGER (kind=4) :: ndime,nelem,nnode,ngaus,axesc,istop
 REAL    (kind=8) :: coord(:,:),eule0(:,:),posgp(:),              &
                     weigh(:),shapf(:,:),deriv(:,:)
 LOGICAL flag
 TYPE (ele08), POINTER :: heade

 !LOCAL variables
 INTEGER, PARAMETER ::  ndofe = 6, nstre = 6
 REAL    (kind=8) gpcod(ndime,nnode),x(ndime,nnode),lbd(9,nnode)

 INTEGER (kind=4) g,l,ielem,iprnt
 TYPE (ele08), POINTER :: e

 INTERFACE
   INCLUDE 'shape9.h'
   INCLUDE 'inrotm.h'
 END INTERFACE

 IF(flag) THEN
   CALL shape9(weigh,posgp,shapf,deriv,nnode,ngaus)
   WRITE(lures,1500,ERR=9999)
 END IF
 e => heade
 DO ielem=1,nelem
   !       gathers the initial coordinates and local coordinate transf.
   x(1:ndime,1:nnode) = coord(1:ndime,e%lnods(1:nnode))
   DO l=1,nnode
     CALL inrotm(eule0(:,e%lnods(l)),lbd(:,l))
   END DO
   IF(flag .AND. e%jac(1) == 0d0) THEN

     CALL setga8(ndime,nstre,nnode,ngaus,axesc,shapf,deriv,x,     &
                 RESHAPE(lbd,(/3,3,nnode/)),gpcod,e%llbd,e%jac,e%stran,e%sedef)
     DO g=1,ngaus
       IF(e%jac(g) <= 0) THEN
         WRITE(lures,1000,ERR=9999) e%numel,g,e%jac(g)
         istop = 1
         RETURN
       END IF
       iprnt = 0
       IF(ABS(e%stran(1,g)-1d0) > 1.e-5) iprnt = 1
       DO l=2,nstre
         IF(ABS(e%stran(l,g)) > 1.e-5) iprnt = 1
       END DO
       IF(iprnt == 1)THEN
         WRITE(lures,2000,ERR=9999) e%numel,g,e%stran(1,g)-1d0, &
                                    (e%stran(l,g),l=2,nstre )
       END IF
     END DO
   ELSE
     CALL proma1(gpcod,x,shapf,ndime,ngaus,nnode)
   END IF
   e => e%next
 END DO
 IF(flag)   axesc = ABS( axesc )
 RETURN
 1000 FORMAT(' Negative or null jacobian at beam element (8) ',i5,/ &
        &    ' at gauss point ',i2,/,' Jacobian = ', e12.4 )
 1500 FORMAT('  Initial strains in beams (if exist) ')
 2000 FORMAT('  Elemt.',i5,' gauss point ',i2,/,6e12.4)
 9999 CALL runen2('')
 END SUBROUTINE gauss8
