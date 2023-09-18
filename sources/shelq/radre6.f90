 SUBROUTINE radre6(matty,a1,a2,amb,tgp,stn,st,stp,ehist,cmat,prop,istop )

 !     performs radial return for j2 plasticity in shell element

 IMPLICIT NONE
 INTEGER (kind=4),INTENT(IN) :: matty
 INTEGER (kind=4),INTENT(OUT) :: istop
 REAL (kind=8),INTENT(IN) :: a1(3),a2(3),amb(2),tgp(3),stn(8),cmat(*),prop(*)
 REAL (kind=8),INTENT(IN OUT) :: st(8),stp(8),ehist(5)

 REAL   (kind=8) ap1(3),ap2(3),ab1(3),ab2(3),deter,l11,l12,l21,l22, &
                 f1,coefm,coenm,k,kp,epstr,khard,str(8),expo,       &
                 dsp(8),t11,t12,t21,t22,eas(8),def,ws(8),rfs

 REAL(kind=8),SAVE :: pmat(3) = (/0.3333333333333333d0,1d0,2d0/),   &
                      toler   = 1.0d-04

 def = amb(2)           !present thicknnes ratio
 f1  = prop(6)*def**2   !h^2/12
 IF(matty == 1) THEN
   ! curvilinear system with unit metric (schmidt)
   ab1 = a1
   CALL vecuni(3,ab1,t22)
   t11 = DOT_PRODUCT(ab1,a2)
   ab2 = a2 - t11*ab1
   CALL vecuni(3,ab2,t22)
   ! evaluates contravariant base
   t11  = ab2(2)*tgp(3) - ab2(3)*tgp(2)
   t22  =-ab1(2)*tgp(3) + ab1(3)*tgp(2)
   t12  = ab1(2)*ab2(3) - ab1(3)*ab2(2)
   deter = ab1(1)*t11 + ab2(1)*t22 + tgp(1)*t12
   ap1(1) = t11/deter
   ap1(2) = (-ab2(1)*tgp(3)  + ab2(3)*tgp(1)) /deter
   ap1(3) = ( ab2(1)*tgp(2)  - ab2(2)*tgp(1)) /deter
   ap2(1) = t22/deter
   ap2(2) = ( ab1(1)*tgp(3)  - ab1(3)*tgp(1)) /deter
   ap2(3) = (-ab1(1)*tgp(2)  + ab1(2)*tgp(1)) /deter

   l11 = DOT_PRODUCT(ap1,a1)
   l12 = DOT_PRODUCT(ap1,a2)
   l21 = DOT_PRODUCT(ap2,a1)
   l22 = DOT_PRODUCT(ap2,a2)

   deter = l11*l22-l12*l21

   t11 =  l22/deter
   t12 = -l12/deter
   t21 = -l21/deter
   t22 =  l11/deter

   ! elastic Almansi strains (paso n+1)
   eas(1)=stn(1)*t11**2  +stn(2)*t21**2  +stn(3)*2d0*t11*t21
   eas(2)=stn(1)*t12**2  +stn(2)*t22**2  +stn(3)*2d0*t12*t22
   eas(3)=stn(1)*t11*t12 +stn(2)*t21*t22 +stn(3)*(t11*t22 +t21*t12)
   eas(4)=stn(4)*t11**2  +stn(5)*t21**2  +stn(6)*2d0*t11*t21
   eas(5)=stn(4)*t12**2  +stn(5)*t22**2  +stn(6)*2d0*t12*t22
   eas(6)=stn(4)*t11*t12 +stn(5)*t21*t22 +stn(6)*(t11*t22 +t21*t12)
   eas(7)=stn(7)*t11/def +stn(8)*t21/def
   eas(8)=stn(7)*t12/def +stn(8)*t22/def
   ! plastic Almansi strain (paso n)
   dsp(1)=stp(1)*t11**2  +stp(2)*t21**2  +stp(3)*2d0*t11*t21
   dsp(2)=stp(1)*t12**2  +stp(2)*t22**2  +stp(3)*2d0*t12*t22
   dsp(3)=stp(1)*t11*t12 +stp(2)*t21*t22 +stp(3)*(t11*t22 +t21*t12)
   dsp(4)=stp(4)*t11**2  +stp(5)*t21**2  +stp(6)*2d0*t11*t21
   dsp(5)=stp(4)*t12**2  +stp(5)*t22**2  +stp(6)*2d0*t12*t22
   dsp(6)=stp(4)*t11*t12 +stp(5)*t21*t22 +stp(6)*(t11*t22 +t21*t12)
   dsp(7)=stp(7)*t11/def +stp(8)*t21/def
   dsp(8)=stp(7)*t12/def +stp(8)*t22/def
 ELSE
   eas = stn
   dsp = stp
 END IF
 ! Kirchhoff stresses (linear relation with Almansi strain)
 CALL istgp6(1,8,str(1),eas(1),cmat(1),amb(1),matty,prop(1),cmat(1))
 ! IF kinematic hardening
 khard = prop(5)    !H'h*2/3
 IF(khard > 0) THEN
   str(1) = str(1) - khard*dsp(1)
   str(2) = str(2) - khard*dsp(2)
   str(3) = str(3) - khard*dsp(3)
   str(4) = str(4) - khard*dsp(4)*f1
   str(5) = str(5) - khard*dsp(5)*f1
   str(6) = str(6) - khard*dsp(6)*f1
   str(7) = str(7) - khard*dsp(7)*def
   str(8) = str(8) - khard*dsp(8)*def
 END IF
 ! actual yield stress and hardening parameter
 epstr= ehist(3)
 expo = prop(3)
 rfs  = prop(4)
 IF(expo == 0d0) THEN
   k  = prop(1) + prop(2)*epstr
   kp = prop(2)
 ELSE IF( rfs == 0d0 )THEN
   IF(epstr > 0.02d0) THEN
     k  = prop(1)*(prop(2)+epstr)**expo
     kp = prop(1)*expo/(prop(2)+epstr)**(1d0-expo)
   ELSE
     kp = prop(1)*expo/(prop(2)+0.02d0)**(1d0-expo)
     k  = prop(1)*(prop(2)+0.02d0)**expo-kp*(0.02d0-epstr)
   END IF
 ELSE
   k  = prop(1)+prop(2)*epstr+(rfs-prop(1))*(1d0-1d0/EXP(expo*epstr)) !linear + saturation law hardening
   kp = prop(2) + (rfs-prop(1))*expo/EXP(expo*epstr)            !derivative
 END IF

 epstr = 0d0
 coefm = prop(8)/def**2
 coenm = prop(9)/def

 CALL j2sst9(8,3,k,kp,khard,toler,f1,prop(7),coefm,coenm,    &
             prop(10),prop(13),prop(14),prop(17),str,epstr,  &
             ehist(2),ehist(4),def,ws,pmat,istop)

 IF(epstr > 0d0) dsp = dsp + ws
 IF(khard > 0) THEN    !add back stress
   str(1) = str(1) + khard*dsp(1)
   str(2) = str(2) + khard*dsp(2)
   str(3) = str(3) + khard*dsp(3)
   str(4) = str(4) + khard*dsp(4)*f1
   str(5) = str(5) + khard*dsp(5)*f1
   str(6) = str(6) + khard*dsp(6)*f1
   str(7) = str(7) + khard*dsp(7)*def
   str(8) = str(8) + khard*dsp(8)*def
 END IF
 IF(epstr > 0d0) THEN
   IF(matty == 1) THEN
     ! Green-Lagrange plastic strains
     stp(1)=dsp(1)*l11**2 +dsp(2)*l21**2 +2*dsp(3)*l11*l21
     stp(2)=dsp(1)*l12**2 +dsp(2)*l22**2 +2*dsp(3)*l12*l22
     stp(3)=dsp(1)*l11*l12+dsp(2)*l21*l22+dsp(3)*(l11*l22 +l21*l12)
     stp(4)=dsp(4)*l11**2 +dsp(5)*l21**2 +2*dsp(6)*l11*l21
     stp(5)=dsp(4)*l12**2 +dsp(5)*l22**2 +2*dsp(6)*l12*l22
     stp(6)=dsp(4)*l11*l12+dsp(5)*l21*l22+dsp(6)*(l11*l22 +l21*l12)
     stp(7)=dsp(7)*l11*def+dsp(8)*l21*def
     stp(8)=dsp(7)*l12*def+dsp(8)*l22*def
   ELSE
     stp = dsp
   END IF
   ehist(1) = ehist(3) + epstr
 END IF
 IF(matty == 1) THEN
   ! 2nd Piola-Kirchhoff stresses
   st(1) =str(1)*t11**2  +str(2)*t12**2  +2*str(3)*t11*t12
   st(2) =str(1)*t21**2  +str(2)*t22**2  +2*str(3)*t21*t22
   st(3) =str(1)*t11*t21 +str(2)*t12*t22 +str(3)*(t11*t22 +t12*t21)
   st(4) =str(4)*t11**2  +str(5)*t12**2  +2*str(6)*t11*t12
   st(5) =str(4)*t21**2  +str(5)*t22**2  +2*str(6)*t21*t22
   st(6) =str(4)*t11*t21 +str(5)*t12*t22 +str(6)*(t11*t22 +t12*t21)
   st(7) =str(7)*t11/def +str(8)*t12/def
   st(8) =str(7)*t21/def +str(8)*t22/def
 ELSE
   st = str
 END IF
 RETURN
 END SUBROUTINE radre6
