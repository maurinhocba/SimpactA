 SUBROUTINE dsurf1(nsurf,ndime,ndofn,npoin,iwrit,coord,loadv,heads)
!     Apply load distributed over surface
 USE lispa0, ONLY : lures
 USE loa_db
 IMPLICIT NONE
 INTEGER (kind=4), INTENT (IN) :: nsurf,ndime,ndofn,npoin,iwrit
 REAL (kind=8), INTENT(IN) :: coord(ndime,npoin)
 REAL (kind=8), INTENT(IN OUT) :: loadv(ndofn,npoin)
 TYPE (srf_nod),POINTER :: heads

 INTEGER (kind=4), PARAMETER :: nnmax = 6
 INTEGER (kind=4) :: i,n,l,nn,isurf,nvert,&
                     noprs(nnmax),g,ll(3)
 REAL (kind=8) :: x(ndime,nnmax),direc(ndime),&
                  dx(ndime),dy(ndime),dz(ndime),&
                  rload(ndime,nnmax),px(ndime),p(nnmax),&
                  auxv,area,fc
 REAL (kind=8) weigh(3,3),posgp(9),shape(21),deriv(21),gldir(ndime)
 INTEGER (kind=4) :: ii(4)= (/6,4,5,6/), kk(4)= (/4,5,6,5/)
 TYPE (srf_nod),POINTER :: srf

 CHARACTER :: systm
 INTEGER (kind=4) chnode

 IF(iwrit == 1) WRITE(lures,"(//6x,'No. of Loaded Surfaces =',i5,&
&                            /)",ERR=9999) nsurf

 !*** distributed surface loads section

 IF(iwrit == 1) WRITE(lures,"(/6x,'list of loaded Faces  and ',&
&  'Applied loads'/)",ERR=9999)

 srf => heads
 DO isurf = 1,nsurf

   systm = srf%systm3
   nn = srf%nn3
   noprs(1:nn) = srf%lnod3(1:nn)
   nvert = nn
   IF( nn > 4) nvert = nn/2
   DO i=1,nn
     noprs(i) = chnode(noprs(i))
   END DO
   direc = 0d0
   direc(1:ndime) = srf%direc3(1:ndime)
   CALL vecuni(ndime,direc,area)
   IF (area == 0d0) direc(3) = 1d+0
   p(1:nn) = srf%press3(1:nn)

   !       prints actual DATA
   IF(iwrit == 1) THEN
     WRITE(lures,"(' FACE',i3,' Sys= ',a1,' NN=',i2,' DIR'&
&         ,3f6.3,' Nodes',6i6)",ERR=9999) isurf,systm,nn,direc,noprs(1:NN)
     WRITE(lures,"(' Nodal p ',6g12.3)",ERR=9999) p(1:nn)
   END IF

   rload = 0d0
   x(1:ndime,1:nn) = coord(1:ndime,noprs(1:nn))

   SELECT CASE (nn)
   CASE (3)
     dx = x(1:ndime,2) - x(1:ndime,1)
     dy = x(1:ndime,3) - x(1:ndime,1)
     CALL vecpro(dx,dy,dz)
     CALL vecuni(ndime,dz,area)
     area = area/24d0    !area/12

     !         direction of the load

     IF ( systm == 'L') THEN
       !         SELECT local x=t1 in the global xy plane
       dx = (/ -dz(2), dz(1) , 0d0 /)
       !         of course local y = t2 = t3 x t1
       dy=(/-dz(1)*dz(3),-dz(2)*dz(3),(dz(1)*dz(1)+dz(2)*dz(2))/)
       IF(ABS(dy(3)) < 1.0d-5) THEN
         dx = (/  1d0, 0d0, 0d0 /)
         dy = (/  0d0, 1d0, 0d0 /)
       ELSE
         !           normalizes t1 & t2
         CALL vecuni(ndime,dx,auxv)
         CALL vecuni(ndime,dy,auxv)
       END IF
       !  transform direc from local to global
       gldir = direc(1)*dx+direc(2)*dy+direc(3)*dz

     ELSE

       gldir = direc

     END IF
     px(1) = 2*p(1)+  p(2)+  p(3)
     px(2) =   p(1)+2*p(2)+  p(3)
     px(3) =   p(1)+  p(2)+2*p(3)

     DO g=1,nvert
       auxv = area*px(g)
       rload(1:ndime,g) = rload(1:ndime,g) + gldir*auxv
     END DO

   CASE (4)
     CALL gaussq(2 ,posgp ,weigh )
     DO g = 1,4
       CALL shape3(deriv,shape,posgp(MOD(g-1,2)+1),posgp((g+1)/2),     4)
       dx(1:ndime) = MATMUL(x(1:ndime,1:4),deriv(1:4))
       dy(1:ndime) = MATMUL(x(1:ndime,1:4),deriv(5:8))
       CALL vecpro(dx,dy,dz)
       CALL vecuni(ndime,dz,area)

       ! direction of the load
       IF (systm == 'L') THEN
         ! transform direc from local to global
         !         SELECT local x=t1 in the global xy plane
         dx = (/ -dz(2), dz(1) , 0d0 /)
         !         of course local y = t2 = t3 x t1
         dy = (/-dz(1)*dz(3), -dz(2)*dz(3),(dz(1)*dz(1)+dz(2)*dz(2))/)
         IF(ABS(dy(3)) < 1.0d-5) THEN
           dx = (/  1d0, 0d0, 0d0 /)
           dy = (/  0d0, 1d0, 0d0 /)
         ELSE
           !           normalizes t1 & t2
           CALL vecuni(ndime,dx,auxv)
           CALL vecuni(ndime,dy,auxv)
         END IF
         gldir = direc(1)*dx+direc(2)*dy+direc(3)*dz

       ELSE

         gldir = direc

       END IF

       auxv = DOT_PRODUCT(shape(1:4),p(1:4))*area
       DO n = 1,nn
         fc = auxv*shape(n)
         rload(1:ndime,n) = rload(1:ndime,n) + fc*gldir
       END DO

     END DO
   CASE (6)
     ! for each sub-triangle
     DO n=1,4
       ll(1) = n
       ll(2) = ii(n)
       ll(3) = kk(n)
       dx = x(1:ndime,kk(n)) - x(1:ndime,n)
       dy = x(1:ndime,ii(n)) - x(1:ndime,n)
       CALL vecpro(dx,dy,dz)
       CALL vecuni(ndime,dz,area)
       area = area/24d0

       !           direction of the load

       IF ( systm == 'L') THEN
         !         SELECT local x=t1 in the global xy plane
         dx = (/ -dz(2), dz(1) , 0d0 /)
         !         of course local y = t2 = t3 x t1
         dy = (/-dz(1)*dz(3), -dz(2)*dz(3),(dz(1)*dz(1)+dz(2)*dz(2))/)
         IF(ABS(dy(3)) < 1.0d-5) THEN
           dx = (/  1d0, 0d0, 0d0 /)
           dy = (/  0d0, 1d0, 0d0 /)
         ELSE
           !           normalizes t1 & t2
           CALL vecuni(ndime,dx,auxv)
           CALL vecuni(ndime,dy,auxv)
         END IF
         !  transform direc from local to global
         gldir = direc(1)*dx+direc(2)*dy+direc(3)*dz

       ELSE

         gldir = direc

       END IF
       px(1) = 2*p(ll(1))+  p(ll(2))+  p(ll(3))
       px(2) =   p(ll(1))+2*p(ll(2))+  p(ll(3))
       px(3) =   p(ll(1))+  p(ll(2))+2*p(ll(3))

       DO g=1,nvert
         auxv = area*px(g)
         rload(1:ndime,ll(g))= rload(1:ndime,ll(g)) + gldir*auxv
       END DO
     END DO
   CASE DEFAULT
     CALL runend('DSURF1:ERROR NN can be only 3 4 6 ')
   END SELECT

   DO n=1,nn
     l=noprs(n)
     loadv(1:ndime,l) = loadv(1:ndime,l) + rload(1:ndime,n)
   END DO
   srf => srf%next
 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE dsurf1
