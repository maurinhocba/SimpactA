 SUBROUTINE rare18(sigma,props,epdef,kp,toler,young,k0,yield,flux, &
                   istop,ielem)
 !***********************************************************************
 !
 !     backward euler scheme for plasticity in beams (model 1)
 !
 !***********************************************************************
 IMPLICIT NONE
 !     routine parameters
 INTEGER (kind=4) istop,ielem
 REAL    (kind=8) sigma(6),props(12),epdef,kp,toler,young,k0, &
                  yield(16),flux(6)
 !     local variables
 INTEGER (kind=4) i,j,l,ii,jj,kk,ll,nf,ny,step
 REAL    (kind=8) s(6),dsig(6),sigy,error,de,toly
 !     normalizes stresses for yield FUNCTION
 DO i=1,6
   s(i)=sigma(i)/props(i+6)
 END DO
 !     square of the normalized yield stress
 sigy = (1d0+kp*epdef)**2
 !     verifies yield condition
 ii = 2
 jj = 2
 kk = 2
 ll = 2
 IF(ABS(s(5)*s(6)) < toler) ii = 1
 IF(abs(s(2)*s(4)) < toler) jj = 1
 IF(abs(s(1)*(abs(s(5))+abs(s(6)))) < toler) kk = 1
 IF(abs(s(4)*s(3)) < toler) ll = 1
 nf = ii*jj*kk*ll
 CALL yiel18(ii,jj,kk,ll,s,sigy,yield)
 ny = 0
 de = 0
 DO i=1,nf
   IF(yield(i) > de) THEN
     ny = i
     de = yield(i)
   END IF
 END DO
 IF(ny == 0)RETURN
 !     SAVE trial stresses
 sigma = s
 step = 0
 toly = SQRT(toler)
 !     begins iterative loop,  repeat until (convergence or step > 100)
 100 IF(step > 100)GO TO 200
   step = step+1
   IF(step == 2) THEN
     i = i+1
   END IF
   CALL flux18(ii,jj,kk,ll,ny,sigma,flux,toly)
   dsig = sigma - s
   !       solves the value of lambda
   de = 0.0
   error = 0.0
   DO l=1,6
     de = de + flux(l)*props(l)*flux(l)/props(6+l)**2
   END DO
   error = yield(ny) + DOT_PRODUCT(flux,dsig)
   de    = error/de

   DO j=1,6
     sigma(j)=s(j)-de*props(j)*flux(j)/props(6+j)**2
   END DO

   CALL yiel18(ii,jj,kk,ll,sigma,sigy,yield)
   !       check convergence
 IF(yield(ny) > toler) GO TO 100
 !     elasto-plastic deformation
 de = 0
 DO i=1,6
   dsig(i) = (s(i)-sigma(i))/props(i)*props(i+6) !delta plastic strain
   sigma(i) = sigma(i)*props(i+6)                !real stress resultant
   de = de + dsig(i)*sigma(i)                    !plastic work
 END DO
 sigy  = k0*(1d0+kp*epdef)
 epdef = epdef + de*young/props(1)/sigy
 RETURN
 200 WRITE(lures,1000,ERR=9999) ielem,SQRT(sigy),(s(i),i=1,6)
 istop = 1

 RETURN
 1000 FORMAT(' PROGRAM halted in radial return algorithm for beams ', &
       'model 1',/,'  Element number ',i5,/                           &
       ' input DATA is as follows: sigy =',f10.6,/,                   &
       ' and the normalized stresses are ', 6(t40,e12.4,/))
 9999 CALL runen2('')
 END SUBROUTINE rare18

 !***********************************************************************
 SUBROUTINE yiel18(ii,jj,kk,ll,s,sigy,yield)
 !***********************************************************************
 !
 !     yield functions for beams (model 1)
 !
 !***********************************************************************
 IMPLICIT NONE
 INTEGER (kind=4) ii,jj,kk,ll
 REAL    (kind=8) s(6),sigy,yield(16)

 INTEGER (kind=4) i,j,k,l,surf
 REAL    (kind=8) a,b,c,d,me

 REAL (kind=8), PARAMETER :: signo(2) = (/ 0.288675134,-0.288675134 /)

 a = s(1)**2+s(2)**2+s(3)**2+s(4)**2-sigy
 surf = 0
 DO i=1,ii
   me = s(5)**2+signo(i)*s(5)*s(6)+s(6)**2
   b  = a + me
   me = SQRT(me)
   DO j=1,jj
     c  = b + signo(j)*s(2)*s(4)
     DO k=1,kk
       d = c + signo(k)*s(1)*me
       DO l=1,ll
         surf = surf+1
         yield(surf) = d + signo(l)*s(4)*s(3)
       END DO
     END DO
   END DO
 END DO

 END SUBROUTINE yiel18
 !***********************************************************************
 SUBROUTINE flux18(ii,jj,kk,ll,ny,s,flux,toly)
 !***********************************************************************
 !
 !     flux vectors for assoc. plasticity for beams (model 1)
 !
 !***********************************************************************
 IMPLICIT NONE
 INTEGER (kind=4) ii,jj,kk,ll,ny
 REAL    (kind=8) s(6),flux(6),toly

 INTEGER (kind=4) i,j,k,l,surf
 REAL    (kind=8) me,md,b
 REAL (kind=8), PARAMETER :: signo(2) = (/ 0.288675134,-0.288675134 /)

 surf = 0
 DO i=1,ii
   me = SQRT(s(5)**2+signo(i)*s(5)*s(6)+s(6)**2)
   md = 1d0
   IF(me > 0) md = 0.5d0/me
   DO j=1,jj
     DO k=1,kk
       DO l=1,ll
         surf = surf+1
         IF(ny == surf) THEN
           b = 1d0 + signo(k)*s(1)*md
           flux(1) = 2d0*s(1)
           flux(2) = 2d0*s(2)
           flux(3) = 2d0*s(3)
           flux(4) = 2d0*s(4)
           flux(5) = 2d0*s(5)*b
           flux(6) = 2d0*s(6)*b
           IF(ABS(s(1)) > toly)flux(1) = flux(1) + signo(k)*me
           IF(ABS(s(2)) > toly)flux(2) = flux(2)+signo(j)*s(4)
           IF(ABS(s(3)) > toly)flux(3) = flux(3)+signo(l)*s(4)
           IF(ABS(s(4)) > toly)flux(4) = flux(4)+signo(j)*s(2)+signo(l)*s(3)
           IF(ABS(s(5)) > toly)flux(5) = flux(5)+signo(i)*s(6)
           IF(ABS(s(6)) > toly)flux(6) = flux(6)+signo(i)*s(5)
           RETURN
         END IF
       END DO
     END DO
   END DO
 END DO

 END SUBROUTINE flux18
