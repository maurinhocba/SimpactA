 SUBROUTINE btstr9(nn,ntype,deriv,shape,phin,x,r0,bts,dx,jac,s)
 !******************************************************************************
 !
 !     computes the contribution of gauss point to internal nodal forces
 !     Tranpose(B) S
 !     for beam/shell elements
 !
 !******************************************************************************
 IMPLICIT NONE
 !*** dummy arguments
 INTEGER (kind=4), INTENT(IN) :: nn,   & !number of nodes
                                 ntype   !problem tiype
 REAL (kind=8), INTENT(IN) :: deriv(:), & !(nn)Nodal function derivatives (natural)
                              shape(:), & !(nn)nodal functions
                              phin(:),  & !(nn)nodal angles
                              x(:,:),   & !(2,nn)nodal coordinates
                              r0,       & !radius at element center
                              dx(:),    & !(2)configuration derivative
                              jac,      & !jacobian
                              s(:)        !(nstre)generalized stress x facto
 REAL (kind=8), INTENT(OUT) :: bts(:,:)   !(2,nn) internal forces

 !*** local variables
 INTEGER (kind=4) n
 REAL    (kind=8) t1,t2,phi,dphi,r,derxi(nn),aux,n1,n2,m1,m2,q

 IF(ntype == 3) r = DOT_PRODUCT(x(1,:),shape) !radius at Gauss point
 phi  = DOT_PRODUCT(phin,shape)               !angle at Gauss Point
 derxi= deriv/jac                             !Nodal function derivatives (arc-length)
 dphi = DOT_PRODUCT(phin,derxi)               !angle derivative (arc-length)
 t1 = COS(phi)                                !tangent vector (r or x1 component)
 t2 = SIN(phi)                                !tangent vector (z or x2 component)
 n1 = s(1)                                    !axial force o membrane 1 force
 IF(ntype == 1) THEN             !plane stress problem
   m1 = s(2)                       !bending moment
   q  = s(3)                       !shear force
 ELSE                            !plane strain or axilsymmetric problem
   m1 = s(3)                       !bending moment - 1
   q  = s(5)                       !shear force
 END IF
 DO n=1,nn                       !for each node
   SELECT CASE (ntype)             !plane stress problem
   CASE (1)
     !       axial force contribution
     bts(1,n)= bts(1,n)+ n1*derxi(n)*t1
     bts(2,n)= bts(2,n)+ n1*derxi(n)*t2
     bts(3,n)= bts(3,n)+ n1*shape(n)*(dx(1)*(-t2)+dx(2)*t1)
     !       bending contribution
     bts(3,n)= bts(3,n)+ m1*derxi(n)
     !       transverse shear contribution
     bts(1,n)= bts(1,n)- q*derxi(n)*t2
     bts(2,n)= bts(2,n)+ q*derxi(n)*t1
     bts(3,n)= bts(3,n)- q*shape(n)*(dx(1)*t1+dx(2)*t2)
   CASE (2:3)                      !plane strain or axilsymmetric problem
     !       axial force contribution
     bts(1,n)= bts(1,n)+ n1*derxi(n)*dx(1)
     bts(2,n)= bts(2,n)+ n1*derxi(n)*dx(2)
     !       bending contribution
     bts(1,n)= bts(1,n)- m1*t1*dphi*derxi(n)
     bts(2,n)= bts(2,n)- m1*t2*dphi*derxi(n)
     bts(3,n)= bts(3,n)- m1*(dphi*shape(n)*(-t2*dx(1)+t1*dx(2))+ &
                                  derxi(n)*( t1*dx(1)+t2*dx(2)))
     IF(ntype == 3) THEN  !axilsymmetric
       n2 = s(2)                !membrane force 2 (hoop)
       m2 = s(4)                !bending moment 2 (hoop)
       aux = shape(n)/r0**2     !N/r^2
       !         hoop stress contribution
       bts(1,n) = bts(1,n) + n2 *aux*r
       !         hoop moment contribution
       bts(1,n) = bts(1,n) - m2 *aux*t2
       bts(3,n) = bts(3,n) - m2 *aux*t1*r
     END IF
     !       transverse shear contribution
     bts(1,n)= bts(1,n)+ q*derxi(n)*(-t2)
     bts(2,n)= bts(2,n)+ q*derxi(n)*t1
     bts(3,n)= bts(3,n)- q*shape(n)*(t1*dx(1)+t2*dx(2))
   END SELECT
 END DO
 RETURN
 END SUBROUTINE btstr9
