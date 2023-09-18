 SUBROUTINE chckar(nsegm,nnseg,lnods,x,label)
 !
 !     Checks aspect ratio of triangles defining a master surface
 !     writes warnings when a specified value is exceeded
 !
 IMPLICIT NONE

 INTEGER (KIND=4), INTENT(IN) :: nsegm,nnseg,lnods(:,:),label(:)
 REAL (KIND=8), POINTER :: x(:,:)

 INTEGER i,j,k,jm(1),jj,nar
 REAL (KIND=8) :: y(3,nnseg),s(3,nnseg),l(nnseg),lmax,h(3),ar,proy,maxar

 IF(nnseg /= 3) RETURN    !only for triangles

 maxar = 0d0 !initializes

 nar = 0
 DO i=1,nsegm                        !for each segment
   y=x(1:3,lnods(1:nnseg,i))         !nodal coordinates
   DO j=1,nnseg                      !for each node of the triangle
     k = MOD(j,nnseg) + 1               !next node
     s(1:3,j) = y(1:3,k) - y(1:3,j)     !vector j->k
     l(j) = SQRT(DOT_PRODUCT(s(1:3,j),s(1:3,j)))  !length of side j->k
   END DO
   jm   = MAXLOC(l) !position of longest side of the triangle
   jj = jm(1)       !
   k = MOD(jj,nnseg) + 1  !node k
   lmax = l(jj)           !length of longest side
   proy = DOT_PRODUCT(s(1:3,jj),s(1:3,k))/lmax**2 !Proyection of side k
   h = s(1:3,k) - proy*s(1:3,jj)      !Orthogonal projection
   ar = lmax/SQRT(DOT_PRODUCT(h,h))   !aspect ratio
   IF( ar > 25d0 )THEN
      WRITE(55,"(' WARNING: Aspect ratio excessive in segment',          &
              &  4i7,' Value: ',f7.2)",ERR=9999) i,label(lnods(1:3,i)),ar
      WRITE(55,"(20X,'Internal ',3i7)",ERR=9999) lnods(1:3,i)  !internal numeration
      IF( ar > maxar ) maxar = ar
      nar = nar + 1
   END IF
 END DO
 IF( nar > 0 ) WRITE(*,"(' WARNING: Aspect ratio excessive in',i5,' segments',/,       &
               & ' Maximum Aspect Ratio detected',f6.1,' See file xxx.rep for details')")nar,maxar

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE chckar
