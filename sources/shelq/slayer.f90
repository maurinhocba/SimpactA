 SUBROUTINE slayer(ngaus,nstre,strsg,elstr,d,prop,strpl,stres, &
                   ambda,nlayr,thico,istop,ehist,ielem,plast)
 !*****************************************************************************
 !
 !*****evaluates total resultant stresses for layered shell element
 !     for linear elastic isotropic material
 !
 !****************************************************************************
 IMPLICIT NONE

 !                        routine parameters

 INTEGER (kind=4), INTENT(IN) :: ngaus, & !number of integration points
                                 nstre, & !number of stress components
                                 nlayr, & !number of layer in thickness
                                 ielem    !element number
 INTEGER (kind=4), INTENT(OUT) :: istop   !error flag

 REAL (kind=8), INTENT(IN ) ::     &
               elstr(nstre,ngaus), & !1st, 2nd FF and shear strains
               prop(17),           & !plastic material properties
               d(5),               & !elastic material properties
               thico                 !original thickness
 REAL (kind=8), INTENT(IN OUT) ::      &
               strpl(6,nlayr,ngaus),   & !layer plastic strains
               stres(5,nlayr,ngaus),   & !layer stresses
               ehist(5,ngaus),         & !global plasticity parameters
               ambda(2,ngaus)            !thickness ratios
 REAL (kind=8), INTENT(OUT) :: strsg(nstre,ngaus)  !Nij, Mij, Qi
 LOGICAL, INTENT(IN) :: plast

 !local variables
 ! Gauss points throught the thickness, 1 to 5 points t t thickness
 REAL (kind=8), PARAMETER ::  wei(5,5) = RESHAPE( (/  &
 1.0000000000000D+00,0.0000000000000D+00,0.0000000000000D+00,0.0000000000000D+00,0.0000000000000D+00, &
 5.0000000000000D-01,5.0000000000000D-01,0.0000000000000D+00,0.0000000000000D+00,0.0000000000000D+00, &
 2.7777777777778D-01,4.4444444444445D-01,2.7777777777778D-01,0.0000000000000D+00,0.0000000000000D+00, &
 1.7392742256873D-01,3.2607257743127D-01,3.2607257743127D-01,1.7392742256873D-01,0.0000000000000D+00, &
 1.1846344252809D-01,2.3931433524968D-01,2.8444444444445D-01,2.3931433524968D-01,1.1846344252809D-01  &
 /),(/5,5/) )
 REAL (kind=8), PARAMETER ::  thf(5,5) = RESHAPE( (/  &
 0.0000000000000D+00,0.0000000000000D+00,0.0000000000000D+00,0.0000000000000D+00,0.0000000000000D+00, &
 -2.886751345948D-01,2.8867513459481D-01,0.0000000000000D+00,0.0000000000000D+00,0.0000000000000D+00, &
 -3.872983346207D-01,0.0000000000000D+00,3.8729833462074D-01,0.0000000000000D+00,0.0000000000000D+00, &
 -4.305681557970D-01,-1.699905217924D-01,1.6999052179243D-01,4.3056815579703D-01,0.0000000000000D+00, &
 -4.530899228193D-01,-2.692346550528D-01,0.0000000000000D+00,2.6923465505284D-01,4.5308992281933D-01  &
 /),(/5,5/) )

 REAL (kind=8) :: strail(5), & !trial strain
                  thick,     & !present thickness
                  zlayr,     & !layer z coordinate
                  strial(5), & !trial stress
                  strsl(5),  & !final stres
                  efstr,     & !effective stress
                  efpst,     & !effective plastic strain
                  r1,r2,lb(2)  !eigen-pair

 INTEGER (kind=4) g,l,ierr

 strsg = 0.d0                      !initializes integrated stresses
 DO g=1,ngaus                      !for each Gauss point
   thick = thico*ambda(2,g)        !present thickness
   DO l=1,nlayr                    !for each layer
     zlayr = thf(l,nlayr)*thick    ! z coord

     !      compute elastic strains at a layer
     strail(1:3) = elstr(1:3,g) + elstr(4:6,g)*zlayr   !layer U^2
     strail(4:5) = elstr(7:8,g)
     ierr = 0                          !initializes
     CALL lgst14(strail,r1,r2,lb,'SLAYER',ierr)  !Hencky (logarithmic) strains
     IF( ierr == 1 )THEN          !too distorted element
       WRITE(55,"(i5,7e12.4)",ERR=9999) ielem,elstr(1:6,g),zlayr
       istop = ierr
       EXIT
       !CALL runend('SLAYER: negative eigenvalues of U^2')
     END IF
     IF( plast ) strail = strail - strpl(1:5,l,g)  !elastic Hencky strains
     ! compute elastic (trial) Kirchhoff stresses
     strial(1) = d(1) * strail(1) + d(2) * strail(2)
     strial(2) = d(2) * strail(1) + d(3) * strail(2)
     strial(3) = d(4) * strail(3)
     strial(4) = d(5) * strail(4)
     strial(5) = d(5) * strail(5)

     !           elastic-plastic calculations
     !           backward euler for plane stress
     !           strsl <- strial
     IF( plast )THEN
       efpst = strpl(6,l,g)
       CALL corr06(strial(1),strpl(1,l,g),efpst, &
                   prop(1),d(1),prop(6),prop(12),ierr,efstr,ielem)
       IF( ierr == 1 )THEN
         istop = ierr
         CYCLE
       END IF
       strpl(6,l,g) = strpl(6,l,g) + efpst                 !layer equivalent plastic strain
       ehist(1,g) = ehist(1,g) + strpl(6,l,g)*wei(l,nlayr)        !average equivalent Mises stress
       ehist(2,g) = ehist(2,g) + efstr*wei(l,nlayr)        !average equivalent Mises stress
     END IF
     !keep layer stresses for tangent moduli
     stres(1:5,l,g) = strial(1:5)

     ! Computes Hencky stress on the natural Frame

     strsl(1) = strial(1)*r1*r1+strial(2)*r2*r2+2d0*strial(3)*r1*r2
     strsl(2) = strial(1)*r2*r2+strial(2)*r1*r1-2d0*strial(3)*r1*r2
     strsl(3) =(strial(2)-strial(1))*r1*r2+strial(3)*(r1*r1-r2*r2)
     ! Computes 2nd P-K stress on the natural Frame
     strial(1) = strsl(1)/lb(1)**2
     strial(2) = strsl(2)/lb(2)**2
     IF( ABS(lb(1)-lb(2)) > 1.d-6)THEN   !lb(1) /= lb(2)
       strial(3) =strsl(3)*2d0*LOG(lb(1)/lb(2))/(lb(1)**2-lb(2)**2)
     ELSE                                !lb(1) = lb(2)
       strial(3) = strsl(3)/lb(1)/lb(2)
     END IF
    ! Computes 2nd P-K on the Lagrangian Frame
     strsl(1) = strial(1)*r1*r1+strial(2)*r2*r2-2d0*strial(3)*r1*r2
     strsl(2) = strial(1)*r2*r2+strial(2)*r1*r1+2d0*strial(3)*r1*r2
     strsl(3) =(strial(1)-strial(2))*r1*r2+strial(3)*(r1*r1-r2*r2)

     !           integrate stresses in the thickness
     strsg(1:3,g) = strsg(1:3,g) + strsl(1:3)*wei(l,nlayr)       !Nij
     strsg(4:6,g) = strsg(4:6,g) + strsl(1:3)*zlayr*wei(l,nlayr) !Mij
     strsg(7:8,g) = strsg(7:8,g) + strial(4:5)*wei(l,nlayr)      !Qi

   END DO
   ehist(4,g) = strpl(6,1,g)        !equivalent plastic strain at first layer
   ehist(5,g) = strpl(6,nlayr,g)        !equivalent plastic strain at first layer
   ! ***********
 END DO
 strsg = strsg * thico !Original thickness (TLF)
 RETURN
 9999 CALL runen2('')
 END SUBROUTINE slayer

 SUBROUTINE corr06(st,dstpl,efpst,props,c,b,d,ierr,fi,ielem)
 !-------------------------------------------------------------------
 !
 !     Planar and transversal Anisotropy
 !
 !-------------------------------------------------------------------
 USE lispa0
 IMPLICIT NONE

 INTEGER (kind=4), INTENT(IN) :: ielem
 REAL (kind=8),INTENT(IN) :: props(4)     !material properties
 REAL (kind=8),INTENT(IN) :: c(5)         !elasticity matrix (orthotropic)
 REAL (kind=8),INTENT(IN) :: b(6)         !flow rule matrix
 REAL (kind=8),INTENT(IN) :: d(6)         !yield function derivative
 REAL (kind=8),INTENT(IN OUT) :: st(5)    !trial and corrected stresses
 REAL (kind=8),INTENT(IN OUT) :: efpst    !effective plastic strain
                                          !(IN) present   (OUT) increment
 REAL (kind=8),INTENT(IN OUT) :: dstpl(5)    !plastic strains
 REAL (kind=8),INTENT(OUT) :: fi          !equivalent von Mises Stress
 INTEGER (kind=4), INTENT(OUT) :: ierr    !error flag (0: O.K., 1:error)

 !local variables

 REAL (kind=8) :: c0,expn,aprim,e0
 REAL (kind=8) :: s11,s22,s12,s13,s23,yield,epbar,f0,ap,ddl, &
                  a1,a2,a3,a4,a5,r1,r2,r3,r4,r5,rr,det,po,rfs

 REAL (kind=8), PARAMETER :: toler=1d-4, toler1=1d-6
 INTEGER (kind=4), PARAMETER :: miter=3

 INTEGER (kind=4) k
 !INTEGER (kind=4), SAVE :: kk(miter)=0 ,j=0

 !------------ begin

 !     setup initial yield FUNCTION radius

 c0    = props(1)     !Initial yield or C0 constant
 expn  = props(3)     !exponent for non-linear hardening
 rfs   = props(4)     !residual flow stress
 IF( expn == 0d0)THEN
   aprim = props(2)   !linear hardening
 ELSE
   e0 = props(2)      !non-linear hardening
 END IF
 epbar = efpst
 IF( expn == 0d0)THEN
   yield = c0 + aprim*epbar     !linear hardening
 ELSE IF( rfs == 0d0 )THEN
   yield = c0*(e0+epbar)**expn            !non-linear (exponential) hardening
   aprim = expn*c0/(e0+epbar)**(1d0-expn) !derivative
 ELSE
   yield = c0+e0*epbar+(rfs-c0)*(1d0-1d0/EXP(expn*epbar)) !linear + saturation law hardening
   aprim = e0 + (rfs-c0)*expn/EXP(expn*epbar)            !derivative
 END IF
 !     calculate initial effective stress

 s11 = st(1)
 s22 = st(2)
 s12 = st(3)
 s13 = st(4)
 s23 = st(5)

 ! derivative of yield function
 a1 = (d(1)*s11 + d(2)*s22)
 a2 = (d(2)*s11 + d(3)*s22)
 a3 =  d(4)*s12
 a4 =  d(5)*s13
 a5 =  d(6)*s23
 fi = SQRT(a1*s11+a2*s22+a3*s12+a4*s13+a5*s23)

 !     check initial yield FUNCTION

 f0 = fi-yield
 IF((f0-toler*fi) <=  0d0) THEN

    !   IF the point is elastic
   efpst = 0d0
   RETURN               !no plastic flow

 ELSE
   !  initial flow rule

   !   start iteration for satisfying consistency condition
   DO k = 1,miter
     ! derivative of yield function
     a1 = a1/fi
     a2 = a2/fi
     a3 = a3/fi
     a4 = a4/fi
     a5 = a5/fi
     ! compute flow rule
     r1 = (b(1)*s11 + b(2)*s22)
     r2 = (b(2)*s11 + b(3)*s22)
     r3 =  b(4)*s12
     r4 =  b(5)*s13
     r5 =  b(6)*s23
     po = SQRT(r1*s11+r2*s22+r3*s12+r4*s13+r5*s23)   !potencial
     r1 = r1/po
     r2 = r2/po
     r3 = r3/po
     r4 = r4/po
     r5 = r5/po
     rr = po/fi
     ap = rr*aprim
     ddl = f0/(c(1)*r1*a1+2d0*c(2)*r1*a2+c(3)*r2*a2+c(4)*r3*a3+ &
               c(5)*r4*a4+c(5)*r5*a5 +ap)
     ! new stres
     s11 = s11 - ddl*(c(1)*r1+c(2)*r2)
     s22 = s22 - ddl*(c(2)*r1+c(3)*r2)
     s12 = s12 - ddl*c(4)*r3
     s13 = s13 - ddl*c(5)*r4
     s23 = s23 - ddl*c(5)*r5
     !   calculate the effective plastic strain
     epbar = epbar + ddl*rr
     !   update the radius of the yield surface
     IF( expn == 0d0)THEN
       yield = c0 + aprim*epbar     !linear hardening
     ELSE IF( rfs == 0d0 )THEN
       yield = c0*(e0+epbar)**expn            !non-linear (exponential) hardening
       aprim = expn*c0/(e0+epbar)**(1d0-expn) !derivative
     ELSE
       yield = c0+e0*epbar+(rfs-c0)*(1d0-1d0/EXP(expn*epbar)) !linear + saturation law hardening
       aprim = c0 + (rfs-c0)*expn/EXP(expn*epbar)             !derivative
     END IF
     !   calculate the effective stress
     ! derivative of yield function
     a1 = (d(1)*s11 + d(2)*s22)
     a2 = (d(2)*s11 + d(3)*s22)
     a3 =  d(4)*s12
     a4 =  d(5)*s13
     a5 =  d(6)*s23
     fi = SQRT(a1*s11+a2*s22+a3*s12+a4*s13+a5*s23)
     f0 = fi-yield     !yield function
     !  IF consistency condition is satisfied exit loop
     IF( ABS(f0/yield) <= toler1 )EXIT
   END DO

   IF( k <= miter )THEN
     !  stress change
     a1 = st(1) - s11
     a2 = st(2) - s22
     a3 = st(3) - s12
     a4 = st(4) - s13
     a5 = st(5) - s23
     !  assign corrected stresses
     st = (/ s11,s22,s12,s13,s23 /)
     ! increments in plastic strains
     det = c(1)*c(3) - c(2)*c(2)
     dstpl(1) = dstpl(1) + (a1*c(3)-c(2)*a2)/det
     dstpl(2) = dstpl(2) + (a2*c(1)-c(2)*a1)/det
     dstpl(3) = dstpl(3) + a3/c(4)
     dstpl(4) = dstpl(4) + a4/c(5)
     dstpl(5) = dstpl(5) + a5/c(5) !Cxz = Cyz
     efpst = epbar - efpst         !increment in effective plastic strain
     !j = j + 1
     !kk(k) = kk(k) + 1
     !IF( j == 10000000 )THEN
     !  WRITE(55,"(5i10)",ERR=9999) kk
     !  j = 0
     !  kk = 0
     !END IF
   ELSE
     !        IF the iterations are greater than maximum STOP
     ierr=1                    !set flag to error in plastic integration
     WRITE(*,*)' NO convergence in constitutive equation * STOP *'
     WRITE(lures,1000,ERR=9999) ielem    !print element number
     !WRITE(lures,*,ERR=9999) props  !material properties
     !WRITE(lures,*,ERR=9999) c      !elasticity matrix (orthotropic)
     !WRITE(lures,*,ERR=9999) b      !flow rule matrix
     !WRITE(lures,*,ERR=9999) d      !yield function derivative
     !WRITE(lures,*,ERR=9999) st     !trial and corrected stresses
     !WRITE(lures,*,ERR=9999) efpst  !effective plastic strain
     !WRITE(lures,*,ERR=9999) dstpl  !plastic strains
     !WRITE(lures,*,ERR=9999) fi          !equivalent von Mises Stress
     !WRITE(lures,*,ERR=9999) ierr    !error flag (0: O.K., 1:error)
   END IF
 END IF

 RETURN
 1000 FORMAT(' Program will be stopped. No convergence in the return ', &
     &       'algorithm, elem. no:',i8)
 9999 RETURN


 END SUBROUTINE corr06
!******************************************************
      SUBROUTINE backstr(strial,c,a,dlam,sigma)
      IMPLICIT NONE
      REAL (kind=8) sigma(3),strial(3),c(3),a(3),dlam
      sigma(1)=strial(1)-dlam*(a(1)*c(1)+a(2)*c(2))
      sigma(2)=strial(2)-dlam*(a(1)*c(2)+a(2)*c(1))
      sigma(3)=strial(3)-dlam*a(3)*c(3)
      RETURN
      END SUBROUTINE backstr
!******************************************************
      SUBROUTINE flow(sigma,a,efstr)
      IMPLICIT NONE
      REAL (kind=8) sigma(3),efstr,a(3)
      IF(efstr == 0.)efstr=SQRT(sigma(1)*sigma(1)+  &
       sigma(2)*sigma(2)-sigma(1)*sigma(2)+3.*sigma(3)*sigma(3))
      a(1)=(2*sigma(1)-sigma(2))/(2.*efstr)
      a(2)=(2*sigma(2)-sigma(1))/(2.*efstr)
      a(3)=6*sigma(3)/(2.*efstr)
      RETURN
      END SUBROUTINE flow
!******************************************************
      SUBROUTINE yield2(c1,c2,dlamprim,xmu,r,syield,f2,conver,df)
      IMPLICIT NONE
      REAL (kind=8) syield,f2,xmu,r,c1,c2,dlamprim,auxil1,auxil2,df, &
            toler
      LOGICAL conver
      DATA toler/1.d-9/
      auxil1=1.+dlamprim*xmu*r
      f2=c1/auxil1/auxil1
      auxil2=1.+3.*dlamprim*xmu
      f2=f2+c2/auxil2/auxil2
      f2=f2/4.-syield*syield
!
!       convergence check
!
      df=f2/syield/syield
      IF(ABS(df) < toler) THEN
        conver=.TRUE.
        RETURN
      END IF
      df=c1*r/auxil1/auxil1/auxil1
      df=df+c2*3./auxil2/auxil2/auxil2
      df=-df*xmu/2.
      RETURN
      END SUBROUTINE yield2
!********************************************************************
      SUBROUTINE yield(sx,sy,tau,f,syield,efstr)
      IMPLICIT NONE
      REAL (kind=8) efstr,syield,f,sx,sy,tau
      efstr=SQRT(sx*sx+sy*sy-sx*sy+3.*tau*tau)
      f=efstr-syield
      RETURN
      END SUBROUTINE yield
