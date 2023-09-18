     sec => psecs(isec)%p             !point to section
     nlayr = sec%iprop(1)             !number of layers
     nvar  = sec%iprop(2)             !number of internal variables per layer
     secty = sec%secty                !section type
     thick = sec%rprop(1)             !original thickness
     minstr= sec%rprop(4)             !strain threshold to use TTTI
     osec = isec                      !keep present section
     shear = sec%iprop(4) == 1        !consider transverse shear stresses in constitutive relation
     visco = .FALSE.
     !***** Allocate vectors thf, wei and shf *****************************************
     IF (ALLOCATED(thf)) DEALLOCATE(thf, wei, shf, pflag)
     ALLOCATE(thf(nlayr), wei(nlayr), shf(nlayr), pflag(nlayr))
     !**************************SHELL2*************************************************
     IF( secty == 12 )THEN            !standard solid section
       min_tr= sec%rprop(2)             !Minimum Thickness ratio
       max_tr= sec%rprop(3)             !Maximum Thickness ratio
       nl=nlayr
       mat => sec%mtbas                 !point to associated material
       shell = nlayr > 1                !membrane or shell
       IF( .NOT.shell )rl(:,4:6) = 0d0
       mtype = mat%mtype                !type of base material
       natst = logst .OR. mtype == 6 .OR. mtype == 10   !use log strains
       elast = mat%matdef(3) == 1       !elastic
       CALL gaussq(nlayr,thf(1),wei(1)) !integration points through the thickness
       thf(1:nlayr) = thf(1:nlayr)/2d0  !positions
       wei(1:nlayr) = wei(1:nlayr)/2d0  !weights
       alpha = 0d0  !alpha = mat%prope(6)             !Thermical Dilatation Coeff
       !IF( .NOT.itemp ) alpha = 0d0
       dm(1:8) = sec%rprop(5:12)             !integrated elasticity matrix
       poiss = mat%prope(2)                  !poisson ratio
       db = dm(5)                            !bending stiffness

       SELECT CASE (mtype)

       CASE ( 1 )              !standard isotropic material
         IF( elast )THEN
           plast = .FALSE.
         ELSE   !for elasto-plastic mats
           ! e1, nu1, uniaxial, efren, consn, r, exponent m, hill 79
           propi(1:7) = mat%propp(1:7)       ! isotropic & kinematic hardening parameters
           deatht = mat%propp(5)             !end of plasticity
           propi(5) = REAL( mat%matdef(4),8) ! isotropic hardening model
           chi    = mat%propp(16:27)         ! hill coefficients
           IF (mat%matdef(4) == 5 )THEN      !if Isotropic Hardening POINts
             val => mat%chead%val            !point to first curve
             numpt = mat%chead%np            !number of points in the curve
           ELSE
             NULLIFY (val)
             numpt = 0
           END IF
           plast = propi(1) > 0  !consider plasticity ?
           IF(plast)THEN                     !plasticity present
             SELECT CASE (mat%matdef(3))     !according to yield surface
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
             plast = .NOT.elastic
           END IF
         END IF
         c(1:4) = mat%prope(7:10)          ! plane stress elasticity matrix
         gh = 5d0/6d0*mat%prope(3)*sec%rprop(1)  !elastic transverse shear factor 5/6 Gt

       CASE ( 5 )            ! orhthotropic material
         c(1:4) = mat%prope(16:19)         ! plane stress elasticity matrix
         IF( .NOT.elast )THEN              ! if elasto-plastic
           propi(1:13) = mat%propp(17:29)    !orthotropic hardening parameters
           deatht = mat%propp(5)             !end of plasticity
           plast = deatht > ttime .AND. propi(1) > 0  !consider plasticity ?
         ELSE
           plast = .FALSE.
         END IF

       CASE ( 6 )            !hyperelastic isotropic
         chi(1:12) = mat%prope(7:18)           ! elastic energy coefficients
         IF( .NOT.elast )THEN                  ! elastic-plastic material
           propi(1:5) = mat%propp(1:5)         ! plastic properties
           plast = propi(1) > 0 ! .AND. propp(5) > ttime
         END IF
       !-----------------------------------------------------------------------------
       CASE ( 8 )              !visco-elastic isotropic material
         plast = .FALSE.                          !no plasticity
         visco = .TRUE.                           !elastic-viscous behavior
         r1 = mat%prope(11)*thick/3d0             ! eta*t/3
         r2 = r1*thick**2/4d0                     ! eta*t^3/12
         gh = 5d0/6d0*mat%prope(3)*sec%rprop(1)   !elastic transverse shear factor 5/6 Gt
         c(1:4) = mat%prope(7:10)          ! plane stress elasticity matrix
         natst = .FALSE.                          !small strain version only

       !-----------------------------------------------------------------------------
       CASE ( 30 )              !user defined material
         plast = .TRUE.   !default
       END SELECT

       !- Obtain shear factors from the rprop array ---------------------------------
       shf(1:nlayr) = sec%rprop(13:12+nlayr)  !shear factors (isotropic material by now)

     !**************************SHELL3*************************************************
     ELSE  !secty == 13  composite laminae
       min_tr= 0d0                      !Minimum Thickness ratio
       max_tr= 1d3                      !Maximum Thickness ratio
       natst = logst                    !use log strains
       nl=1                               !size of second index in varin
       elast = sec%iprop(3) == 0               !elastic problem
       plast = .NOT.elast                      !consider plasticity ?
       shell = .TRUE.                          !bending included
       dm = sec%rprop(6:26)                 !linear elastic integrated matrix
       poiss = 0d0                      !what can we put here?
       db = dm(16)                      !elastic bending stiffness
       s1 = dm(7)/dm(1)                  !ratio
       stine(1:9) = ABS(dm(13:21)*s1)
       coupled = ANY(stine(1:9) > 1d-8 ) !are coupling terms relevant?
       IF( plast )THEN  !plastic lamina
         IF( ALLOCATED ( rr ) )DEALLOCATE( thickl,zposl,rr,matsl,lvari )
         ALLOCATE( thickl(nlayr),zposl(nlayr),rr(5,nlayr),matsl(nlayr),lvari(nlayr) )
         i = 27                                !pointer to layer thickness
         CALL vecasi(nlayr,sec%rprop(i),thickl(1))        !layer thickness
         i = i+2*nlayr                         !pointer to layer z position
         CALL vecasi(nlayr,sec%rprop(i),zposl(1))         !layer z position
         i = i+nlayr                           !pointer to layer rotation matrix
         CALL vecasi(nlayr*5,sec%rprop(i),rr(1,1))        !layer rotation matrix
         lvari(1:nlayr)= sec%iprop(nlayr+5:2*nlayr+4)
         DO l=1,nlayr
           mate= sec%iprop(4+l)
           mat => pmats(mate)%p        !search layer material
           matsl(l)%p => mat
         END DO
         !- Obtain shear factors from the rprop array ---------------------------------
         IF( shear ) shf(1:nlayr) = sec%rprop(27+8*nlayr:26+9*nlayr)  !shear factors
       END IF
       oldm = -1                              !recompute constants

     !*********************************************************************************
     END IF
     IF( plast .OR. logst )THEN
       IF(ALLOCATED( varin ) )THEN
          IF( SIZE(varin,1) /= nvar .OR. SIZE(varin,2) /= nl ) DEALLOCATE(varin)
       END IF
       IF( .NOT.ALLOCATED(varin)) ALLOCATE( varin(nvar,nl) )
       IF(plast) plast = .NOT.elastic       ! for elastic behavior with previous internal variables
     END IF
     newmt = .FALSE.                        !same material than previous ?
