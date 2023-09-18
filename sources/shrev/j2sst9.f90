 SUBROUTINE j2sst9(nstre,nn,yistr,ishar,khard,toler,f1,coefn,coefm,&
&                  coenm,ddmat,gkh,pdmat,zeta,strsg,calfa,eqvst,   &
&                  conpa,ambda,stria,pymat,istop)
 !******************************************************************************
 !
 !     this routine compute the stresses for the j2 shell plasticity
 !
 !  input
 !     strsg: trial stress (n+1)
 !  output
 !     strsg: actual stress (n+1)
 !     stria: increment in plastic strains (n+1)
 !     calfa: hardening PARAMETER
 !     eqvst: equivalent stress
 !     conpa(2): consistency parameters
 !******************************************************************************
 USE lispa0
 IMPLICIT NONE
 !***  parameters
 INTEGER (kind=4) nstre,nn,istop
 REAL (kind=8) strsg(nstre),calfa,eqvst,ddmat(nn),pdmat(nn),f1,    &
&              yistr,ishar,khard,toler,coefn,coefm,coenm,gkh,zeta, &
&              ambda,stria(nstre),conpa(2),pymat(nn)
 !***  local variables
 LOGICAL surf(2)
 INTEGER (kind=4) i,j,iyiel,jyiel,kyiel,lyiel,icont,nq
 REAL    (kind=8) yield(2),pcons(2),jmati(2,2),zima1(3),zima2(3),  &
&                 zima3(3),zima4(3),zeta1,actpc(2),h,              &
&                 value,auxi1,auxi2,auxi3,auxi4,actyi,spcon,dpcon, &
&                 grads(nstre,2),flows(nstre,2),stres(nstre)
 REAL (kind=8), SAVE :: signo(2) = (/ 1d0, -1d0 /)

 !     SAVE trial stress, transform to diagonalized base

 nq = 2*nn+1
 stria(1:nstre) = strsg(1:nstre)
 CALL j2sdia(nn,stria)

 !     compute yields FUNCTION

 icont=0
 CALL j2syie(nstre,nn,nq,yield,stria,pymat,coefn,coefm,coenm,yistr)
 eqvst = 3d0*max(yield(1),yield(2))+yistr**2
 conpa = 0d0

 !     check for plastic step

 IF (yield(1) <= 0) THEN
   IF (yield(2) <= 0) THEN
     RETURN
   ELSE
     kyiel=2
   END IF
 ELSE IF (yield(2) <= 0) THEN
   kyiel=1
 ELSE IF (yield(1) > yield(2)) THEN
   kyiel=3
 ELSE
   kyiel=4
 END IF
 lyiel=0

 !     plastic step : look for consistency

 !     SAVE initial hardening parameters ,initialize consistency parameters

 actyi = yistr
 pcons = 0d0
 surf  = (yield > 0d0)
 DO i = 1,nn
   zima1(i) = 1d0
   zima2(i) = 0d0
   zima3(i) = 0d0
   zima4(i) = 1d0
 END DO
 zeta1 = 1d0
 stres(1:nstre) = stria(1:nstre)
 h = SQRT(1d0/coefn)

 !     iterative loop

 icont = 0
 DO WHILE (icont < 100)
   icont = icont + 1

 !       compute actual flow rule

   CALL j2sflw(nstre,nn,nq,flows,stres,pymat,coefn,coefm,coenm,    &
&              signo)

 !       gradients of efective stress respect consistency parameters

   auxi3 = zeta*zeta1
   DO iyiel = 1,2
     IF(surf(iyiel)) THEN
       DO i=1,nn
         j=i+nn
         auxi1 =   coefn*pdmat(i)*stres(i) +                       &
&                  coenm*pdmat(i)*stres(j)*signo(iyiel)
         auxi2 =  (coenm*pdmat(i)*stres(i)*signo(iyiel) +          &
&                  coefm*pdmat(i)*stres(j))*f1
         grads(i,iyiel) = zima1(i)*auxi1 + zima2(i)*auxi2
         grads(j,iyiel) = zima3(i)*auxi1 + zima4(i)*auxi2
       END DO
       DO i = nq,nstre
         grads(i,iyiel) = stres(i)*auxi3
       END DO
     END IF
   END DO

 !       consider isotropic hardening IF exists
   IF(ishar > 0)  value = 4d0/9d0*actyi*2*ishar/h

 !       compute actual jacobian matrix  "jmati"
   DO iyiel = 1,2
     IF(surf(iyiel)) THEN
       DO jyiel = 1,2
         IF(surf(jyiel)) THEN
           jmati(iyiel,jyiel) = - DOT_PRODUCT(flows(1:nstre,iyiel),&
&                                             grads(1:nstre,jyiel))
 !               INCLUDE the effect of isotropic hardening
          IF(ishar > 0) jmati(iyiel,jyiel)=jmati(iyiel,jyiel)-value
         END IF
       END DO
     END IF
   END DO

 !       solve consistency equation and actualize consistency parameters

   IF (ALL(surf)) THEN
 !                                               both surface active
     value = jmati(1,1)*jmati(2,2)-jmati(1,2)*jmati(2,1)
     auxi1 =   jmati(2,2)/value
     auxi2 = - jmati(1,2)/value
     auxi3 = - jmati(2,1)/value
     auxi4 =   jmati(1,1)/value
     actpc(1) = pcons(1) - auxi1*yield(1) - auxi2*yield(2)
     actpc(2) = pcons(2) - auxi3*yield(1) - auxi4*yield(2)
   ELSE IF (surf(1)) THEN
 !                                               only surface 1
     actpc(1) = pcons(1) - yield(1)/jmati(1,1)
     actpc(2) = 0.0d0
   ELSE IF (surf(2)) THEN
 !                                               only surface 2
     actpc(1) = 0.0d0
     actpc(2) = pcons(2) - yield(2)/jmati(2,2)
   END IF
   IF (actpc(1) < 0) THEN
     IF (actpc(2) < 0) THEN
       IF(ALL(surf)) THEN
         IF(kyiel == 3) THEN
 !               set aside surface two
           lyiel = 2
           surf(2) = .FALSE.
         ELSE IF(kyiel == 4) THEN
 !               set aside surface one
           lyiel = 1
           surf(1) = .FALSE.
         END IF
         actpc = 0d0
       ELSE
         write(lures,"('both surfaces desactivated')",ERR=9999)
         istop = 1
         RETURN
       END IF
     ELSE
 !                               deactivate surface 1
       actpc(1) = 0.0d0
       surf(1)  = .FALSE.
 !                               new proy of surf 2
       actpc(2) = pcons(2) -yield(2)/jmati(2,2)
     END IF
   ELSE
      IF (actpc(2) < 0) THEN
 !                               desactive surface 2
         actpc(2) = 0.0d0
         surf(2)  = .FALSE.
 !                               new proy of surf 1
         actpc(1) = pcons(1) -yield(1)/jmati(1,1)
      END IF
   END IF
   pcons = actpc
   spcon = pcons(1) + pcons(2)
   dpcon = pcons(1) - pcons(2)

 !       compute actual efective stress

   DO i=1,nn
     auxi1 = 1d0 + spcon*coefn*pdmat(i)
     auxi2 =       dpcon*coenm*pdmat(i)
     auxi3 =       dpcon*coenm*pdmat(i)*f1
     auxi4 = 1d0 + spcon*coefm*pdmat(i)*f1
     value = 1d0/(auxi1*auxi4-auxi2*auxi3)
     zima1(i) =   auxi4*value
     zima2(i) = - auxi2*value
     zima3(i) = - auxi3*value
     zima4(i) =   auxi1*value
   END DO
   DO i=1,nn
     j=i+nn
     stres(i) = zima1(i)*stria(i) + zima2(i)*stria(j)
     stres(j) = zima3(i)*stria(i) + zima4(i)*stria(j)
   END DO
   zeta1 = 1d0/(1d0+spcon*zeta)
   stres(nq:nstre) = stria(nq:nstre)*zeta1

 !       compute actual hardening PARAMETER

   calfa = 2d0/3d0*actyi/h*(pcons(1)+pcons(2))

 !       compute yields FUNCTION

   actyi = yistr + ishar*calfa
   CALL j2syie(nstre,nn,nq,yield,stres,pymat,coefn,coefm,coenm,    &
&              actyi)

 !       check

   IF (yield(1) > 0 .AND. lyiel /= 1) surf(1)=.TRUE.
   IF (yield(2) > 0 .AND. lyiel /= 2) surf(2)=.TRUE.
 !       EXIT condition error < toler & lyiel=0
   value=0d0
   IF (surf(1)) value = value + abs(yield(1))/actyi/actyi
   IF (surf(2)) value = value + abs(yield(2))/actyi/actyi
   IF (value < toler) THEN
     IF(lyiel /= 0 )THEN
       IF( yield(lyiel) > 0) THEN
        surf(lyiel) = .TRUE.
        lyiel = 0
       ELSE 
         GO TO 200
       END IF
     ELSE
       GO TO 200
     END IF
   END IF
 END DO
 !     STOP due to non convergence
 WRITE(lures,101,ERR=9999) nstre,nn,f1/ambda**2,yistr,ishar,khard,toler,     &
&      coefn,coefm*ambda**2,coenm*ambda,gkh,zeta,(ddmat(i),i=1,nn)
 WRITE(lures,102,ERR=9999) pdmat(1:nn)
 WRITE(lures,1020,ERR=9999) pymat(1:nn)
  101 FORMAT(5x,' nstre =',i3,' nn =',i3,/,5x,  &
     &       'f1     = ',e24.16,/,5x,           &
     &       'yistr  = ',e24.16,/,5x,           &
     &       'ishar  = ',e24.16,/,5x,           &
     &       'khard  = ',e24.16,/,5x,           &
     &       'toler  = ',e24.16,/,5x,           &
     &       'coefn  = ',e24.16,/,5x,           &
     &       'coefm  = ',e24.16,/,5x,           &
     &       'coenm  = ',e24.16,/,5x,           &
     &       'gkh    = ',e24.16,/,5x,           &
     &       'zeta   = ',e24.16,/,5x,           &
     &       'ddmat  = ',3e24.16)
  102 FORMAT(5x,'pdmat  = ',3e24.16)
 1020 FORMAT(5x,'pymat  = ',3e24.16)
 WRITE(lures,103,ERR=9999) ambda,strsg(1:nstre)
  103 FORMAT(' no convergence in radial RETURN algorithm '/ &
     &  ' input values were ',/,' lambda coef. ',e24.16,/,  &
     &  ' stresses      ',3e24.16,/,5e24.16)
 istop = 1
 RETURN
 !     convergence achieved !!
  200 CONTINUE
 !      write(55,"(8e12.4)",ERR=9999) strsg(1:nstre)
 strsg(1:nstre) = stres(1:nstre)

 !     compute increment in plastic strains

 stres(1:nstre) = stria(1:nstre) - stres(1:nstre)
  DO i = 1,nn
    j = i+nn
    stria(i) = stres(i)/ddmat(i)
    stria(j) = stres(j)/ddmat(i)/f1
  END DO
  stria(nq:nstre) = stres(nq:nstre)/gkh

  !     compute actual stress & plastic strains at non-diagonalized base

  CALL j2scar(nn,strsg)
  CALL j2scar(nn,stria)

  eqvst = 3d0*max(yield(1),yield(2))+actyi**2
  conpa = pcons
  !      write(55,"(8e12.4)",ERR=9999) strsg(1:nstre)
  !      write(55,"(8e12.4)",ERR=9999) stria(1:nstre)
  !      write(55,"(4e12.4)",ERR=9999) conpa(1:2),calfa,SQRT(eqvst)

  RETURN
  9999 CALL runen2('')
  END SUBROUTINE j2sst9
