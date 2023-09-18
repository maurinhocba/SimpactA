 SUBROUTINE nodnor( nelem,nnode,lnods,coord,eule0,euler)
 !******************************************************************************
 !
 !****this routine compute the position angles of the nodal coordinate systems
 !
 !     input:   nelem:         number of elements
 !              nnode:         number of nodes per element
 !              coord:         nodal coordinates
 !              lnods:         element conectivities
 !              euler:  /= 0   angles of the nodes. not modified
 !                      == 0   nodes where to compute euler angles
 !
 !     output:  euler:  angles of the nodes.
 !
 !******************************************************************************
 IMPLICIT NONE

 INTEGER (kind=4), PARAMETER :: nd = 3
 INTEGER (kind=4), INTENT(IN) :: nelem,nnode, lnods(:,:)
 REAL (kind=8), INTENT(IN) :: coord(:,:)
 REAL (kind=8), INTENT(IN OUT) :: eule0(:,:),euler(:,:)


 INTEGER (kind=4)ipoin,ielem,inode,npoin
 REAL (kind=8)  auxil,lbd(nd,nd),eltan(nd,nnode),                  &
                elcod(nd,nnode),elnor(nd,nnode),s(nnode),t(nnode), &
                deriv(nnode,2,nnode),shape(nnode,nnode)
 REAL (kind=8), ALLOCATABLE :: vecnx(:,:),vecnz(:,:)
 INTERFACE
   INCLUDE 'angeul.h'
 END INTERFACE
 !     initialize normal z' and x' vectors
 npoin = MAXVAL(lnods)
 ALLOCATE ( vecnx(nd,npoin), vecnz(nd,npoin) )
 vecnx = 0d0
 vecnz = 0d0

 !   assign local nodal coordinates
 SELECT CASE (nnode)
 CASE (3)
    s(1) = 0d0      !node 1
    t(1) = 0d0
    s(2) = 1d0      !node 2
    t(2) = 0d0
    s(3) = 0d0      !node 3
    t(3) = 1d0
 CASE (4)
    s(1) = -1d0     !node 1
    t(1) = -1d0
    s(2) =  1d0     !node 2
    t(2) = -1d0
    s(3) =  1d0     !node 3
    t(3) =  1d0
    s(4) = -1d0     !node 4
    t(4) =  1d0
 CASE (6)
    s(1) = 0d0      !node 1
    t(1) = 0d0
    s(2) = 1d0      !node 2
    t(2) = 0d0
    s(3) = 0d0      !node 3
    t(3) = 1d0
    s(4) = 0.5d0    !node 4
    t(4) = 0d0
    s(5) = 0.5d0    !node 5
    t(5) = 0.5d0
    s(6) = 0d0      !node 6
    t(6) = 0.5d0
 END SELECT
 DO inode=1,nnode
   !          compute shape functions and derivatives
   CALL shape2(deriv(1,1,inode),s(inode),t(inode),nnode,shape(1,inode))
 END DO
 !     summ elemental z' and x'

 DO ielem=1,nelem
   elcod(1:nd,1:nnode) = coord(1:nd,lnods(1:nnode,ielem))
   CALL norelm(nnode,elcod,elnor,eltan,deriv)
   CALL scater(elnor,nd,nnode,vecnz,nd,npoin,lnods(1,ielem))
   CALL scater(eltan,nd,nnode,vecnx,nd,npoin,lnods(1,ielem))
 END DO

 !     compute euler's angles asociated to nodal coordinate system

 DO ipoin=1,npoin
   !   check if Euler angles exist
   IF( .NOT.ALL( eule0(1:3,ipoin) == 0d0) ) CYCLE
   !   average normal vector
   lbd(1:nd,3) = vecnz(1:nd,ipoin)
   CALL vecuni(nd,lbd(1,3),auxil)
   !   if node have do NOT have associated elements
   IF(auxil == 0d0) CYCLE
   lbd(1:nd,1) = vecnx(1:nd,ipoin)
   CALL vecuni(nd,lbd(1,1),auxil)
   IF( ABS(auxil) > 1d-8 ) THEN
     ! IF it is a border node
     ! make surface normal perpendicular to direction X
     auxil = DOT_PRODUCT(lbd(1:3,1),lbd(1:3,3))
     lbd(1:nd,3) = lbd(1:nd,3) - auxil*lbd(1:nd,1)
     CALL vecuni(nd,lbd(1,3),auxil)
     CALL vecpro(lbd(1,3),lbd(1,1),lbd(1,2))
   ELSE
     ! IF it is an internal node
     !         SELECT local x=t1 in the global xy plane
     lbd(1:nd,1) = (/ -lbd(2,3), lbd(1,3) , 0d0 /)
     !         of course local y = t2 = t3 x t1
     lbd(1:nd,2) = (/ -lbd(1,3)*lbd(3,3), -lbd(2,3)*lbd(3,3), &
                      (lbd(1,3)*lbd(1,3) + lbd(2,3)*lbd(2,3)) /)
     !         IF t3 is just global Z
     IF(ABS(lbd(3,2)) < 1.0d-5) THEN
       lbd(1:nd,1) = (/  1d0, 0d0, 0d0 /)
       lbd(1:nd,2) = (/  0d0, 1d0, 0d0 /)
     ELSE
       !           normalizes t1 & t2
       CALL vecuni(nd,lbd(1,1),auxil)
       CALL vecuni(nd,lbd(1,2),auxil)
     END IF
   END IF
   CALL angeul(lbd,eule0(1:3,ipoin),.TRUE.)
   euler(1:3,ipoin) = eule0(1:3,ipoin)
 END DO
 DEALLOCATE ( vecnx, vecnz )

 RETURN
 END SUBROUTINE nodnor
