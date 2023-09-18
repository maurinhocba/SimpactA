 SUBROUTINE norelm( nnode,elcod,elnor,eltan,deriv)
 !*************************************************************************
 !
 !****this routine calculate the vz' and vy' of local shell system
 !    in the nodes
 !
 !*************************************************************************
 IMPLICIT NONE
 !     constants
 INTEGER (kind=4), PARAMETER :: nd=3, ndiel=2
 !     routine parameters
 INTEGER (kind=4), INTENT(IN) :: nnode
 REAL (kind=8), INTENT(IN) :: elcod(nd,nnode),deriv(nnode,ndiel,nnode)
 REAL (kind=8), INTENT(OUT) :: eltan(nd,nnode),elnor(nd,nnode)

 !     local variables
 INTEGER inode,iside,nvert,nside,inodl,jnodl
 REAL (kind=8)  v(nd,ndiel)

 SELECT CASE (nnode)
 CASE (3)
    nvert = 3
    nside = 0

 CASE (4)
    nvert = 4
    nside = 0

 CASE (6)
    nvert = 3
    nside = 3

 END SELECT

 DO inode=1,nnode
   !          compute tangent vectors on the nodes
   v = MATMUL(elcod,deriv(1:nnode,1:ndiel,inode))
   !          compute normal at the nodes
   CALL vecpro(v(1,1),v(1,2),elnor(1,inode))
 END DO

 !  compute tangential vectors on the vertex
 DO inode=1,nvert
   ! inodl = previous node         jnodl = next node
   IF(nside == 0)THEN
     inodl = MODULO(inode-2,nvert) + 1      !notice MODULO definition
     jnodl = MOD(inode,nvert) + 1
   ELSE
     inodl = MODULO(inode-2,nvert) + 1 + nside
     jnodl = inode + nside
   END IF
   eltan(:,inode) =  elcod(:,jnodl) - elcod(:,inodl) 
 END DO

 !     compute tangential vectors on the sides
 DO iside=1,nside
   ! inodl = previous node         jnodl = next node
   inode = nvert+iside
   inodl = iside
   jnodl = MOD(iside,nside) + 1
   eltan(:,inode) =  elcod(:,jnodl) - elcod(:,inodl) 
 END DO

 RETURN
 END SUBROUTINE norelm
