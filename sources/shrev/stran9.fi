 SUBROUTINE stran9(nn,ng,nstre,ntype,x,dx,phin,stran,shape,  &
                   deriv,jac,r0,tg)
 !***********************************************************************
 !
 !*****this routine computes strain measures for 2-d beam/shell elements
 !
 !***********************************************************************
 IMPLICIT NONE
 !              routine parameters
 INTEGER (kind=4) nn,ng,nstre,ntype
 REAL (kind=8), INTENT(IN) :: deriv(:,:),x(:,:),phin(:),shape(:,:),jac(:),r0(:)
 REAL (kind=8), INTENT(OUT) :: stran(:,:),dx(:,:),tg(:,:)

 !              local variables
 INTEGER (kind=4) g
 REAL    (kind=8) dn(nn),phi,r,xi,t1,t2

 DO g = 1,ng
   r    = DOT_PRODUCT(x(1,:),shape(:,g))
   phi  = DOT_PRODUCT(phin,shape(:,g))
   dn   = deriv(:,g)/jac(g)
   xi   = DOT_PRODUCT(phin,dn)
   dx(:,g) = MATMUL(x,dn)
   t1 = COS(phi)
   t2 = SIN(phi)
   tg(1,g) = -t2
   tg(2,g) =  t1
   !       strains
   SELECT CASE (ntype)
   CASE (1)
     stran(1,g) = dx(1,g)*t1 + dx(2,g)*t2
     stran(2,g) = xi
   CASE (2)
     stran(1,g) = dx(1,g)*dx(1,g) + dx(2,g)*dx(2,g)
     stran(2,g) = 1d0
     stran(3,g) = -xi*(dx(1,g)*t1+dx(2,g)*t2)
     stran(4,g) = 0d0
   CASE (3)
     stran(1,g) = dx(1,g)*dx(1,g) + dx(2,g)*dx(2,g)
     stran(2,g) = (r/r0(g))**2
     stran(3,g) = -xi*(dx(1,g)*t1+dx(2,g)*t2)
     stran(4,g) = - r*t2/r0(g)**2
   END SELECT
   stran(nstre,g) = -t2*dx(1,g) + t1*dx(2,g)
 END DO
 RETURN
 END SUBROUTINE stran9
