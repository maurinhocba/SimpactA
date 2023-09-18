 SUBROUTINE istgp9(matty,ng,nstre,plast,stres,stra0,sede0,sedef,ehist, &
&                  strap,cmat,props,prop1,stran,energ,jac,weigh,istop)
 !**********************************************************************
 !
 !     stresses at gauss points for beams/shells (2-d)
 !
 !**********************************************************************
 IMPLICIT NONE

 LOGICAL, INTENT(IN) :: plast
 INTEGER (kind=4), INTENT(IN) :: matty,ng,nstre
 INTEGER (kind=4), INTENT(OUT) :: istop
 REAL(kind=8), INTENT(IN) :: stra0(:,:),sede0(:),jac(:),weigh(:), &
&                            cmat(*),props(*),prop1(*)
 REAL(kind=8), INTENT(IN OUT) :: stres(:,:),sedef(:),stran(:,:),energ(:)
 REAL(kind=8), POINTER :: ehist(:,:),strap(:,:)

 ! local variables
 INTEGER (kind=4) g,n
 REAL (kind=8) def,f1,coefm,coenm,epstr,k,kp,strsg(nstre),khard,expo,rfs, &
               dstrp(nstre),dstra(nstre),ws(nstre)
 REAL (kind=8),PARAMETER :: toler = 1d-4, pmat1=0.6666666666666666d0,  &
                            pmat2(2)= (/0.333333333333333d0,1d0/)

 DO g=1,ng
   IF(plast) THEN
     khard = props(5)
     dstrp = strap( :,g)
     epstr = ehist(1,g)
     expo = props(3)
     rfs  = props(4)
     IF(expo == 0d0) THEN
       k  = props(1) + props(2)*epstr
       kp = props(2)
     ELSE IF( rfs == 0d0 )THEN
       IF(epstr > 0.02d0) THEN
         k  = props(1) * (props(2)+epstr)**expo
         kp = props(1)*expo/(props(2)+epstr)**(1d0-expo)
       ELSE
         kp = props(1)*expo/(props(2)+0.02d0)**(1d0-expo)
         k  = props(1)*(props(2)+0.02d0)**expo - kp*(0.02d0-epstr)
       END IF
     ELSE
       k = props(1)+props(2)*epstr+  &              !linear +
&          (rfs-props(1))*(1d0-1d0/EXP(expo*epstr)) !saturation law hardening
       kp= props(2) + (rfs-props(1))*expo/EXP(expo*epstr)    !derivative
     END IF
     epstr = 0d0
   END IF
   IF(nstre == 3) THEN
     sedef(g) = sede0(g)/(stran(1,g)*(1d0-stran(3,g)**2))**0.25
     ! def = sedef(g)
     def = 1d0
     stran(3,g) = stran(3,g)*def
     !         total strains
     dstra(1) = stran(1,g)-stra0(1,g)
     dstra(2) = stran(2,g)-stra0(2,g)
     dstra(3) = stran(3,g)-stra0(3,g)
     !         elastic strains in CASE of plasticity
     IF(plast) dstra = dstra - dstrp
     IF(matty == 2) THEN
       !           elastic almansi strains
       dstra(1) = dstra(1)/stran(1,g)
       dstra(2) = dstra(2)/stran(1,g)
       dstra(3) = dstra(3)/SQRT(stran(1,g))/def
     END IF
     !         elastic Kirchhoff stresses  (consider thickness change)
     strsg(1) = dstra(1)*cmat(1)
     strsg(2) = dstra(2)*cmat(2)*def**2
     strsg(3) = dstra(3)*cmat(3)
     !         closest point proyection for plasticity
     IF(plast) THEN
       f1 = prop1(1)*def**2
       coefm = prop1(3)/def**2
       coenm = prop1(4)/def
       IF(matty == 2) THEN
         !             plastic almansi strain
         dstrp(1) = dstrp(1)/stran(1,g)
         dstrp(2) = dstrp(2)/stran(1,g)
         dstrp(3) = dstrp(3)/SQRT(stran(1,g))/def
       END IF
       !         IF kinematic hardening
       IF(khard > 0) THEN
         strsg(1) = strsg(1) - khard*dstrp(1)
         strsg(2) = strsg(2) - khard*dstrp(2)*f1
         strsg(3) = strsg(3) - khard*dstrp(3)
       END IF
       CALL j2sst9(nstre,1,k,kp,khard,toler,f1,                  &
&                  prop1(2),coefm,coenm,cmat(1),cmat(3),         & 
&                  prop1(5),prop1(6),strsg,epstr,ehist(2,g),     &
&                  ehist(4,g),def,ws,pmat1,istop)
       IF(epstr > 0d0) dstrp = ws + dstrp
       IF(khard > 0) THEN
         strsg(1) = strsg(1) + khard*dstrp(1)
         strsg(2) = strsg(2) + khard*dstrp(2)*f1
         strsg(3) = strsg(3) + khard*dstrp(3)
       END IF
       IF(epstr > 0d0) THEN
         IF(matty == 2) THEN
           !               Green-Lagrange plastic strains
           strap(1,g) = dstrp(1)*stran(1,g)
           strap(2,g) = dstrp(2)*stran(1,g)
           strap(3,g) = dstrp(3)*SQRT(stran(1,g))*def
         ELSE
           strap(:,g) = dstrp
         END IF
         ehist(1,g) = ehist(1,g) + epstr
       END IF
     END IF
     !         2nd Piola-Kirchhoff stress resultants
     IF(matty == 2) THEN
       stres(1,g) = strsg(1)/stran(1,g)
       stres(2,g) = strsg(2)/stran(1,g)
       stres(3,g) = strsg(3)/SQRT(stran(1,g))/def
     ELSE
       stres(:,g) = strsg
     END IF

   ELSE
     sedef(g) =sede0(g)/SQRT((stran(1,g)-stran(5,g)**2)*stran(2,g))
     def = sedef(g)
     ! def = 1d0
     stran(5,g) = stran(5,g)*def
     !         total strains
     dstra(1) = (stran(1,g) - stra0(1,g) )/2d0
     dstra(2) = (stran(2,g) - stra0(2,g) )/2d0
     dstra(3) =  stran(3,g) - stra0(3,g)
     dstra(4) =  stran(4,g) - stra0(4,g)
     dstra(5) =  stran(5,g) - stra0(5,g)
     !         elastic strains in CASE of plasticity
     IF(plast) dstra = dstra - dstrp
     !         elastic Almansi strains
     IF(matty == 2) THEN
       dstra(1) = dstra(1)/stran(1,g)
       dstra(2) = dstra(2)/stran(2,g)
       dstra(3) = dstra(3)/stran(1,g)
       dstra(4) = dstra(4)/stran(2,g)
       dstra(5) = dstra(5)/SQRT(stran(1,g))/def
     END IF
     !         kirchhoff stresses
     strsg(1) =  cmat(1)*dstra(1)+cmat(2)*dstra(2)
     strsg(2) =  cmat(1)*dstra(2)+cmat(2)*dstra(1)
     strsg(3) = (cmat(3)*dstra(3)+cmat(4)*dstra(4))*def**2
     strsg(4) = (cmat(3)*dstra(4)+cmat(4)*dstra(3))*def**2
     strsg(5) =  cmat(5)*dstra(5)
     IF(plast) THEN
       f1 = prop1(1)*def**2
       coefm = prop1(3)/def**2
       coenm = prop1(4)/def
       !           plastic Almansi strain
       IF(matty == 2) THEN
         dstrp(1) = dstrp(1)/stran(1,g)
         dstrp(2) = dstrp(2)/stran(2,g)
         dstrp(3) = dstrp(3)/stran(1,g)
         dstrp(4) = dstrp(4)/stran(2,g)
         dstrp(5) = dstrp(5)/SQRT(stran(1,g))/def
       END IF
       !           IF kinematic hardening
       khard = props(5)
       IF(khard > 0) THEN
         !             substract "back stress"
         strsg(1) = strsg(1) - khard*dstrp(1)
         strsg(2) = strsg(2) - khard*dstrp(2)
         strsg(3) = strsg(3) - khard*dstrp(3)*f1
         strsg(4) = strsg(4) - khard*dstrp(4)*f1
         strsg(5) = strsg(5) - khard*dstrp(5)
       END IF
       CALL j2sst9(nstre,2,k,kp,khard,toler,f1,               &
&                 prop1(2),coefm,coenm,prop1(5),prop1(7),     &
&                 prop1(8),prop1(10),strsg,epstr,ehist(2,g),  &
&                 ehist(4,g),def,ws,pmat2(1),istop)
       IF(epstr > 0d0) dstrp = ws + dstrp
       IF(khard > 0) THEN    !add back stress
         strsg(1) = strsg(1) + khard*dstrp(1)
         strsg(2) = strsg(2) + khard*dstrp(2)
         strsg(3) = strsg(3) + khard*dstrp(3)*f1
         strsg(4) = strsg(4) + khard*dstrp(4)*f1
         strsg(5) = strsg(5) + khard*dstrp(5)
       END IF
       IF(epstr > 0d0) THEN
         IF(matty == 2) THEN
           !               Green-Lagrange plastic strains
           strap(1,g) = dstrp(1)*stran(1,g)
           strap(2,g) = dstrp(2)*stran(2,g)
           strap(3,g) = dstrp(3)*stran(1,g)
           strap(4,g) = dstrp(4)*stran(2,g)
           strap(5,g) = dstrp(5)*SQRT(stran(1,g))*def
         ELSE
           strap(:,g) = dstrp
         END IF
         ehist(1,g) = ehist(1,g) + epstr
       END IF
     END IF
     !         2nd Piola-Kirchhoff stress resultants
     IF(matty == 2) THEN
       stres(1,g) = strsg(1)/stran(1,g)
       stres(2,g) = strsg(2)/stran(2,g)
       stres(3,g) = strsg(3)/stran(1,g)
       stres(4,g) = strsg(4)/stran(2,g)
       stres(5,g) = strsg(5)/SQRT(stran(1,g))/def
     ELSE
       stres(:,g) = strsg
     END IF
   END IF
   DO n=1,nstre
     energ(n) = energ(n) + stres(n,g)*dstra(n)*jac(g)*weigh(g)
   END DO
 END DO
 RETURN
 END SUBROUTINE istgp9
