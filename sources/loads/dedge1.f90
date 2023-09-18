 SUBROUTINE dedge1(nedge,ntype,ndime,ndofn,npoin,iwrit,            &
                   coord,loadv,heade)

 !     Apply loads distributed over edges

 USE lispa0
 USE loa_db
 IMPLICIT NONE

 INTEGER (kind=4), INTENT (IN) :: nedge,ntype,ndime,ndofn,         &
                                  npoin,iwrit
 REAL (kind=8), INTENT(IN) :: coord(ndime,npoin)
 REAL (kind=8), INTENT(IN OUT) :: loadv(ndofn,npoin)
 TYPE (edg_nod), POINTER :: heade

 INTEGER (kind=4), PARAMETER :: nnmax = 6
 REAL (kind=8), PARAMETER  :: twopi=6.283185307179586
 INTEGER (kind=4) :: iedge,i,n,l,ig,nn,noprs(nnmax)
 REAL (kind=8) :: x(ndime,nnmax),press(ndime,nnmax),&
                  xa(ndime),dx(ndime),dy(ndime),dz(ndime),&
                  rload(ndime,nnmax),px(ndime),pg(ndime),&
                  dvolu,long,auxv
 REAL (kind=8) weigh(3),posgp(3),shape(3,3),deriv(3,3),&
               gpcod(ndime)
!      INTEGER (kind=4) :: ii(4)= (/6,4,5,6/)  , kk(4)= (/4,5,6,5/)

 CHARACTER :: systm
 INTEGER (kind=4) chnode
 TYPE (edg_nod),POINTER :: edg
 INTERFACE
   INCLUDE 'shape9.h'
 END INTERFACE

 !--------------------------------------------------------------

 IF(iwrit == 1) WRITE(lures,"(//' Number of loaded edges ',i8,//&
&   5x,'List of Loaded Edges and Applied Loads'/)",ERR=9999) nedge

 !*** loop over each loaded edge

 edg => heade   ! point to first loaded edge
 DO iedge = 1,nedge  !for each loaded edge

   systm = edg%systm2   !system L: local  G: Global
   nn = edg%nn2         !number of nodes defining the edge
   noprs(1:nn) = edg%lnod2(1:nn)   !nodes defining the edge
   xa = 0d0                        !local load vector
   IF( ndime == 3 .AND. systm == 'L') xa(1:ndime) = edg%xa(1:ndime)
   DO i=1,nn            !internal nodes defining the edge
     noprs(i) = chnode(noprs(i))
   END DO
   press(1:ndime,1:nn) = edg%press2(1:ndime,1:nn)  !pressure in each direction at each node
   IF(iwrit == 1)THEN
     IF(nn == 2)WRITE(lures,"(5x,' nodes ',2i7,6f10.3)",ERR=9999) &
                               noprs(1:2),press(1:ndime,1:2)
     IF(nn == 3)WRITE(lures,"(5x,' nodes ',3i7,9f10.3)",ERR=9999) &
                               noprs(1:3),press(1:ndime,1:3)
   END IF
   x(1:ndime,1:nn) = coord(1:ndime,noprs(1:nn)) ! coordinates of the nodes of the edge
   rload = 0d0                                  !initialize the equivalent load vector
   !*** enter loop for linear numerical integration
   CALL shape9(weigh,posgp,shape,deriv,nn,nn)   !shape functions and derivatives
   DO ig=1,nn     !number of integration points == nn
     dx = MATMUL(x(1:ndime,1:nn),deriv(1:nn,ig)) !local tangent vector
     CALL vecuni(ndime,dx,long)                  !unit direction and jacobian
     pg = MATMUL(press(1:ndime,1:nn),shape(1:nn,ig)) !pressure at gauss point (1:ndime)
     gpcod = MATMUL(x(1:ndime,1:nn),shape(1:nn,ig))  !Gauss point coordinates
     dvolu = weigh(ig)*long                          !Gauss point associated length
     IF( ndime == 2 )THEN                        !For 2D
       IF(systm == 'L')THEN           !for local system
         px(1) = dx(1)*pg(2) - dx(2)*pg(1)  !global x1
         px(2) = dx(1)*pg(1) + dx(2)*pg(2)  !global x2
       ELSE
         px = pg
       END IF
       IF(ntype == 3) dvolu=dvolu*twopi*gpcod(1)
     ELSE !IF(ndime == 3)THEN
       IF(systm == 'L')THEN
         dz = xa - gpcod
         dz = dz - DOT_PRODUCT(dz,dx) * dx
         CALL vecuni(ndime,dz,auxv)
         CALL vecpro(dz,dx,dy)
         px(1) = DOT_PRODUCT(pg,dx)
         px(2) = DOT_PRODUCT(pg,dy)
         px(3) = DOT_PRODUCT(pg,dz)
       ELSE
         px = pg
       END IF
     END IF

     DO n = 1,nn
       rload(1:ndime,n) = rload(1:ndime,n)+shape(n,ig)*px*dvolu
     END DO
   END DO

   DO n=1,nn
     l = noprs(n)
     loadv(1:ndime,l) = loadv(1:ndime,l) + rload(1:ndime,n)
   END DO
   edg => edg%next

 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE dedge1
