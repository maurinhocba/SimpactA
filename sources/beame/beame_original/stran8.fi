 SUBROUTINE stran8(nn,ng,ndime,nstre,x,dx,lbd,stran,shape, &
                   deriv,jac,lbdg)
 !***********************************************************************
 !
 !*****this routine computes strain measures for beam elements
 !
 !***********************************************************************
 IMPLICIT NONE

 !              routine parameters
 INTEGER (kind=4) nn,ng,ndime,nstre
 REAL    (kind=8) deriv(nn,ng),x(ndime,nn),lbd(3,3,nn),         &
                  stran(nstre,ng),dx(ndime,ng),shape(nn,ng),    &
                  jac(ng),lbdg(3,3,ng)

 !              local variables
 INTEGER (kind=4) i,g,n
 REAL    (kind=8) dn,t(3,3),tg(3,3),dt(3),tb(3),dtb(3),mt,tanm2, &
                  aux,auxm(3,3)

 t(1:3,1) = 0.0D0 !initializes first node rotation vector
 DO n=2,nn !for each other node
   auxm = MATMUL(TRANSPOSE(lbd(:,:,1)),lbd(:,:,n)) !relative matrix rotation respecto to 1st node
   CALL cuate8(auxm(1,1),t(1,n)) !compute rotation vector at other nodes relatives to 1st node
 END DO
 dx = 0d0 !initializes tangent vector dx/ds at each Gauss point

 CALL proma1(tg,t,shape,3,ng,nn) !compute rotation vector at Gauss point
 DO g=1,ng
   dt = 0d0  !initializes derivative of rotation vector
   CALL expot8(tg(1,g),auxm(1,1))         !compute relative rotation matrix
   lbdg(:,:,g) = MATMUL(lbd(:,:,1),auxm)  !compute local system at Gauss point
   DO n=1,nn  !for each node
     dn=deriv(n,g)/jac(g)  !nodal function derivative dN/ds
     DO i=1,3     !for each coordinate
       dx(i,g) = dx(i,g) + dn*x(i,n)  !dx/ds
       dt(i)   = dt(i)   + dn*t(i,n)  !dt/ds
     END DO
   END DO
   !       strains
   stran(1:3,g) = MATMUL( TRANSPOSE(lbdg(:,:,g)), dx(1:3,g) )
   !       curvatures
   mt = SQRT(tg(1,g)**2+tg(2,g)**2+tg(3,g)**2) !rotation vector modulus
   IF(mt <= 0.1e-4) THEN !if rotation is "small"
     tb  = tg(1:3,g)/2d0 !
     dtb = dt/2d0
   ELSE
     tanm2 = TAN(mt/2d0)/mt
     tb  =  tg(1:3,g)*tanm2
     aux = -(1d0-mt/SIN(mt))/mt/mt*DOT_PRODUCT(tg(:,g),dt)
     dtb(1) = tanm2*(dt(1)+aux*tg(1,g))
     dtb(2) = tanm2*(dt(2)+aux*tg(2,g))
     dtb(3) = tanm2*(dt(3)+aux*tg(3,g))
   END IF
   aux = 2d0/(1d0+TAN(mt/2d0)**2)
   stran(4,g) = aux*(dtb(1)+dtb(2)*tb(3)-dtb(3)*tb(2)) !torsion
   stran(5,g) = aux*(dtb(2)+dtb(3)*tb(1)-dtb(1)*tb(3)) !1st bending curvature
   stran(6,g) = aux*(dtb(3)+dtb(1)*tb(2)-dtb(2)*tb(1)) !2nd bending curvature
 END DO
 RETURN
 END SUBROUTINE stran8
