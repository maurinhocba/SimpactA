
   ! check if trought the thickness integration is necessary

   IF( ASSOCIATED(e%gausv) )THEN
     ttti = .TRUE.
     varin = e%gausv
   ELSE IF( plast .OR. logst )THEN
     aux = ABS(stra1(1)-1d0)+ABS(stra1(2)-1d0)+ABS(stra1(3))+ &
          (ABS(stra1(4))+ABS(stra1(5))+ABS(stra1(6)))*thick/2d0
     ttti = aux > minstr
     IF( ttti ) varin = 0d0
   ELSE
     ttti = .FALSE.
   END IF
   IF( ttti ) THEN
     stine(1:6) = 0d0       !Integrated forces and moments
     eql= SQRT(2d0*area1)   ! characteristic length
   END IF
   stine(9:10) = 0d0
   pflag = .FALSE.                    !initializes flag

!**************************SHELL2*************************************************

   IF( secty == 12 )THEN              !for standard solid section
     IF( ttti )THEN
       !Trought the thickness (t.t.t.) integration loop
       large = natst                     !use log strains
       IF( plast ) large = large .OR. ASSOCIATED(e%gausv)  !use log strains
       DO l=1,nlayr                        !for each layer
         zk = thf(l)*thnew                     !Z coordinate
         stran = stra1(1:3)+stra1(4:6)*zk      !layer U^2
         !IF( alpha > 0 )THEN                 !consider thermal strains
         !  j0 = (t0 + zk*t1)**(2d0/3d0)      !modifies jacobian due to temperature
         !  stran = stran/j0                  !modifies C
         !END IF
         IF( large )THEN
           CALL lgst14(stran,r1,r2,lb,'RESVPL',error)    !Hencky (logarithmic) strains
           IF( error == 1 )THEN
             WRITE(55,"(i7,3i6,7e12.4)") e%numel,e%lnods(1:3),stra1(1:6),zk
             WRITE(55,"(3e15.4)") x
             ierr = 1
             error = 0
           END IF
         ELSE
           stran(1:2) = (stran(1:2) - 1d0 )*0.5d0
         END IF

         SELECT CASE (mtype)

         !-----------------------------------------------------------------------------
         CASE ( 1,8 )     !isotropic material
           IF( shear )THEN
             stres(1:2) = stint(7:8,iel)*shf(l)
           ELSE
             stres(1:2) = 0d0
           END IF
           CALL stre14(stran,stres,c,propi,chib,chid,varin(:,l),ierr,3, &
                       plast,elast,val,numpt,aux,mat%matdef(3),pflag(l))
         !-----------------------------------------------------------------------------
         CASE ( 5 )     !orthotropic material
           IF( ASSOCIATED(e%gausv)) THEN
             IF(varin(4,l) > 0d0) stran(1:3) = stran(1:3) - varin(1:3,l)
           END IF
           !elastic (trial) stresses
           stres(1)= c(1)*stran(1)+c(2)*stran(2)
           stres(2)= c(2)*stran(1)+c(3)*stran(2)
           stres(3)= c(4)*stran(3)
           IF( plast) THEN
             CALL corr05(stres(1),stres(2),stres(3),varin(4,l),c(:),propi(:),ierr, &
                         varin(1:3,l),aux)
             IF ( varin(4,l) > 0.d0 )  pflag(l)=.TRUE.
           END IF

         !-----------------------------------------------------------------------------
         CASE ( 6 )  ! Hyperelastic
           lb(3) = 1d0/lb(1)/lb(2)     !thickness ratio
           IF( elast )THEN     !for elastic problems
             CALL rubberps(chi,lb,mat%matdef(8),stres,r1=r1,r2=r2)
             large = .FALSE.   ! 2ndPK stresses already computed
           ELSE                !    elastic-plastic
             efpst = varin(4,l)    !effect plastic strain
             IF( efpst > 0d0 )THEN !if the material has unrecoverable strains
               !elastic (trial) strains
               strpl(1:3) = varin(1:3,l)              !previous (twice) plastic strains
               dstpl(1:3) = stran(1:3) - strpl(1:3)   !trial Elastic log strains (Twice)
               ! compute eigen-decomposition of elastic part
               CALL eige17(dstpl(1),s1,s2,lc(1))      !compute principal log strains
               lc(3) = -lc(1)-lc(2)                   !sum of log strains must be 0
               lc = EXP(lc*0.5)    !lc are twice the principal shear strains
             ELSE                !if the material is elastic
               lc = lb             !eigenvalues and
               s1 = r1             !eigenvectors are the same
               s2 = r2
             END IF

             !   compute deviatoric principal stresses and check plasticity
             CALL core06ps(chi,propi(1:2),mat%matdef(8),pflag(l),lc,stres,efpst,ierr, &
                           s1=s1,s2=s2,ehs=dstpl)
             IF( pflag(l) )THEN  ! if plastic flow in the step
               varin(4,l) = efpst      !effective plastic strain
               varin(1:3,l) = stran(1:3) - dstpl(1:3)   !plastic  log strains
             END IF

           END IF
           aux = SQRT(stres(1)**2+stres(2)**2-stres(1)*stres(2)+3d0*stres(3)**2) !equiv str

         !-----------------------------------------------------------------------------
         CASE ( 30 )      ! User defined material
           CALL ud_shell_2(stran,stres,varin(:,l),ttime,ielem,ierr, &
                           newmt,mat%matdef(8),mat%props,mat%chead)
           aux = SQRT(stres(1)**2+stres(2)**2-stres(1)*stres(2)+3d0*stres(3)**2)  !Von Mises stress

         END SELECT !mtype

         IF( l == 1     ) stine( 9) = aux  !keep equivalent stress at first
         IF( l == nlayr ) stine(10) = aux  !and last layer

         IF( large )THEN  !Compute 2nd Piola-Kirchhoff stresses
          ! Computes Hencky stress on the natural Frame
           sigma(1) = stres(1)*r1*r1+stres(2)*r2*r2+2d0*stres(3)*r1*r2
           sigma(2) = stres(1)*r2*r2+stres(2)*r1*r1-2d0*stres(3)*r1*r2
           sigma(3) =(stres(2)-stres(1))*r1*r2+stres(3)*(r1*r1-r2*r2)
          ! Computes 2nd P-K stress on the natural Frame
           stres(1) = sigma(1)/lb(1)**2
           stres(2) = sigma(2)/lb(2)**2
           IF( ABS(lb(1)-lb(2)) > 1.d-6)THEN   !lb(1) /= lb(2)
             stres(3) = sigma(3)*2d0*LOG(lb(1)/lb(2))/(lb(1)**2-lb(2)**2)
           ELSE                                !lb(1) = lb(2)
             stres(3) = sigma(3)/lb(1)/lb(2)
           END IF
          ! Computes 2nd P-K on the Lagrangian Frame
           sigma(1) = stres(1)*r1*r1+stres(2)*r2*r2-2d0*stres(3)*r1*r2
           sigma(2) = stres(1)*r2*r2+stres(2)*r1*r1+2d0*stres(3)*r1*r2
           sigma(3) =(stres(1)-stres(2))*r1*r2+stres(3)*(r1*r1-r2*r2)
         ELSE
           sigma = stres  !small strain
         END IF

         !***   compute Int(B**t*sigma) on element level

         stine(1:3) = stine(1:3)+sigma*wei(l)     !t.t.t. integrated forces
         IF(shell)stine(4:6) = stine(4:6)+sigma*zk*wei(l)  !integrated moments
       END DO
       stine(1:6) = stine(1:6)*thick                  !Original thick because it's TLF

       IF(ANY(pflag(1:nl)))THEN
         IF( updiv )THEN
           IF(.NOT.ASSOCIATED(e%gausv))  ALLOCATE(e%gausv(nvar,nl))
           e%gausv = varin
         ELSE
           newiv = .TRUE.
         END IF
       END IF

     ELSE          !no t.t.t.i.
       ! membrane part
       stran(1:2) = (stra1(1:2) - 1d0 )*0.5d0            !GL strains
       stine(1) = dm(1)*stran(1) + dm(2)*stran(2)
       stine(2) = dm(2)*stran(1) + dm(3)*stran(2)
       stine(3) = dm(4)*stra1(3)
       ! bending part
       stran(1:2) = stra1(4:5)*0.5d0                     !curvatures
       stine(4) = dm(5)*stran(1) + dm(6)*stran(2)
       stine(5) = dm(6)*stran(1) + dm(7)*stran(2)
       stine(6) = dm(8)*stra1(6)
       stine(9:10) = 0d0
     END IF

!**************************SHELL3*************************************************

   ELSE   !secty = 13
     IF( ttti )THEN
       !Trought the thickness (t.t.t.) integration loop
       DO l=1,nlayr                        !for each layer
         zk = zposl(l)*e%lb                 !Z coordinate
         stran = stra1(1:3)+stra1(4:6)*zk      !layer U^2
         !IF( alpha > 0 )THEN                 !consider thermal strains
         !  j0 = (t0 + zk*t1)**(2d0/3d0)      !modifies jacobian due to temperature
         !  stran = stran/j0                  !modifies C
         !END IF
         IF( natst )THEN
           CALL lgst14(stran,r1,r2,lb,'RESVPL',ierr)    !Hencky (logarithmic) strains
           IF( ierr == 1 )THEN
             !WRITE(55,"(3i5,7e12.4)",ERR=9999) e%lnods(1:3),stra1(1:6),zk
             !WRITE(55,"(3e15.4)",ERR=9999) x
             ierr = 1
             error = 0
           END IF
         ELSE
           stran(1:2) = (stran(1:2) - 1d0 )*0.5d0       !Green strains
         END IF

             !Rotate strains from element basis to layer basis (if necessary)
         IF( rr(1,l) < 0.9999999999 )THEN  ! Angle is not zero
           stral(1) = rr(1,l)*stran(1)+ rr(2,l)*stran(2)+ rr(3,l)*stran(3)
           stral(2) = rr(2,l)*stran(1)+ rr(1,l)*stran(2)- rr(3,l)*stran(3)
           stral(3) =-rr(4,l)*stran(1)+ rr(4,l)*stran(2)+ rr(5,l)*stran(3)
         ELSE
           stral(1:3) = stran(1:3)
         END IF

         IF (l == nlayr) THEN
           ilv = lvari(l)
           jlv = nvar
         ELSE
           ilv = lvari(l)
           jlv = lvari(l+1)-1
         END IF

         mate  = sec%iprop(2+2*l)
         newmt = mate /= oldm
         oldm = mate
         mat => matsl(l)%p        ! pointer to layer material

         SELECT CASE (mat%mtype)
         !-----------------------------------------------------------------------------
         CASE( 1 )    ! isotropic layer material
           IF (newmt) THEN
             propi(1:4) = mat%propp(1:4)       ! isotropic hardening parameters
             propi(5) = REAL( mat%matdef(4),8) ! isotropic hardening model
             chi    = mat%propp(16:27)         ! hill coefficients
             deatht = mat%propp(5)             !end of plasticity
             IF (mat%matdef(4) == 5 )THEN      !iso-hardening defined by points
               val => mat%chead%val
               numpt = mat%chead%np
             ELSE
               NULLIFY (val)
               numpt = 0
             END IF
             plast = deatht > ttime .AND. propi(1) > 0  !consider plasticity ?
             IF(plast)THEN
               SELECT CASE (mat%matdef(3))
               CASE (2,3) !Mises Hill48
                 ! D matrix, derivative of yield function
                 chid(1) = chi(2)+chi(3)  !g+h
                 chid(2) = -chi(3)        !-h
                 chid(3) = chi(1)+chi(3)  !f+h
                 chid(4) = 2d0*chi(6)     !2n
                 chid(5) = 2d0*chi(5)     !2m
                 chid(6) = 2d0*chi(4)     !2l
                 ! B matrix, flow rule matrix
                 chib(1) = chi(8)+chi(9)  !g+h
                 chib(2) = -chi(9)        !-h
                 chib(3) = chi(7)+chi(9)  !f+h
                 chib(4) = 2d0*chi(12)    !2n
               CASE (4)  !Hill-79
                 chid(1:3) = chi(1:3)
               CASE (5)  !Hill-90
                 chid(1:4) = chi(1:4)
                 chib(1:4) = chi(5:8)
               END SELECT
             END IF
             c(1:4) = mat%prope(7:10)          ! plane stress elasticity matrix
           END IF
           IF( shear )THEN
             strel(1:2) = stint(7:8,iel)*shf(l)
           ELSE
             strel(1:2) = 0d0
           END IF
           CALL stre14(stral,strel,c,propi,chib,chid,varin(ilv:jlv,1),ierr,3, &
                       plast,elast,val,numpt,aux,mat%matdef(3),pflag(1))
           !IF ( varin(jlv,1)>0 )  pflag(1)=.TRUE.
         !-----------------------------------------------------------------------------
         CASE ( 5 )   ! orthotropic material
           !elastic (trial) stresses
           IF (newmt) c(1:4) = mat%prope(16:19)     ! plane stress elasticity matrix
           strel(1)= c(1)*stral(1)+c(2)*stral(2)
           strel(2)= c(2)*stral(1)+c(3)*stral(2)
           strel(3)= c(4)*stral(3)
           IF( plast )THEN
             CALL corr05(strel(1),strel(2),strel(3),varin(jlv,l),c(1:4),mat%propp(17:29),ierr, &
                         varin(ilv:jlv-1,1),aux)
             IF ( varin(jlv,1)>0 )  pflag(1)=.TRUE.
           END IF
         !-----------------------------------------------------------------------------
         END SELECT

         !Rotate stresses from layer basis to element basis (if necessary)
         IF( rr(1,l) < 0.9999999999 )THEN  ! Angle is not zero
           stres(1) = strel(1)*rr(1,l)+ strel(2)*rr(2,l)- strel(3)*rr(4,l)
           stres(2) = strel(1)*rr(2,l)+ strel(2)*rr(1,l)+ strel(3)*rr(4,l)
           stres(3) = strel(1)*rr(3,l)- strel(2)*rr(3,l)+ strel(3)*rr(5,l)
         ELSE
           stres(1:3) = strel(1:3)
         END IF


         IF( l == 1     ) stine( 9) = aux
         IF( l == nlayr ) stine(10) = aux

         IF( natst )THEN
          ! Computes Hencky stress on the natural Frame
           sigma(1) = stres(1)*r1*r1+stres(2)*r2*r2+2d0*stres(3)*r1*r2
           sigma(2) = stres(1)*r2*r2+stres(2)*r1*r1-2d0*stres(3)*r1*r2
           sigma(3) =(stres(2)-stres(1))*r1*r2+stres(3)*(r1*r1-r2*r2)
          ! Computes 2nd P-K stress on the natural Frame
           stres(1) = sigma(1)/lb(1)**2
           stres(2) = sigma(2)/lb(2)**2
           IF( ABS(lb(1)-lb(2)) > 1.d-6)THEN   !lb(1) /= lb(2)
             stres(3) = sigma(3)*2d0*LOG(lb(1)/lb(2))/(lb(1)**2-lb(2)**2)
           ELSE                                !lb(1) = lb(2)
             stres(3) = sigma(3)/lb(1)/lb(2)
           END IF
          ! Computes 2nd P-K on the Lagrangian Frame
           sigma(1) = stres(1)*r1*r1+stres(2)*r2*r2-2d0*stres(3)*r1*r2
           sigma(2) = stres(1)*r2*r2+stres(2)*r1*r1+2d0*stres(3)*r1*r2
           sigma(3) =(stres(1)-stres(2))*r1*r2+stres(3)*(r1*r1-r2*r2)
         ELSE
           sigma = stres  !small strain
         END IF

         !***   compute Int(B**t*sigma) on element level

         stine(1:3) = stine(1:3)+sigma*thickl(l)   !t.t.t. integrated forces
         stine(4:6) = stine(4:6)+sigma*zk*thickl(l)  !integrated moments
       END DO
       IF(ANY(pflag(1:nl)))THEN
         IF( updiv )THEN
           IF(.NOT.ASSOCIATED(e%gausv))  ALLOCATE(e%gausv(nvar,nl))
           e%gausv = varin
         ELSE
           newiv = .TRUE.
         END IF
       END IF

     ELSE
       stra1(1:2) = (stra1(1:2) - 1d0 )*0.5d0
       stra1(4:5) = stra1(4:5)*0.5d0
       stine(1) = dm(1)*stra1(1) + dm(2)*stra1(2) + dm(3)*stra1(3)
       stine(2) = dm(2)*stra1(1) + dm(4)*stra1(2) + dm(5)*stra1(3)
       stine(3) = dm(3)*stra1(1) + dm(5)*stra1(2) + dm(6)*stra1(3)
       stine(4) = dm(7)*stra1(4) + dm(8)*stra1(5) + dm(9)*stra1(6)
       stine(5) = dm(8)*stra1(4) + dm(10)*stra1(5) + dm(11)*stra1(6)
       stine(6) = dm(9)*stra1(4) + dm(11)*stra1(5) + dm(12)*stra1(6)
       IF( coupled )THEN
         stine(1) = stine(1) +dm(13)*stra1(4) + dm(14)*stra1(5) + dm(15)*stra1(6)
         stine(2) = stine(2) +dm(16)*stra1(4) + dm(17)*stra1(5) + dm(18)*stra1(6)
         stine(3) = stine(3) +dm(19)*stra1(4) + dm(20)*stra1(5) + dm(21)*stra1(6)
         stine(4) = stine(4) +dm(13)*stra1(1) + dm(16)*stra1(2) + dm(19)*stra1(3)
         stine(5) = stine(5) +dm(14)*stra1(1) + dm(17)*stra1(2) + dm(20)*stra1(3)
         stine(6) = stine(6) +dm(15)*stra1(1) + dm(18)*stra1(2) + dm(21)*stra1(3)
       END IF
       stine(9:10) = 0d0
     END IF
!!!*********************************************************************************
   END IF
