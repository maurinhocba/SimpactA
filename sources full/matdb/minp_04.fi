   SUBROUTINE minp_04 ( mat )
   !            --------------------------------
   !            input material properties for
   !            temperature dependent
   !            elastic - plastic, isotropic
   !            elastic stress-strain curve: linear
   !            plastic stress-strain curve: exponential
   !            --------------------------------
   ! MATDEF: 1: 3      3-D
   !         2: 2      Isotropic
   !         3: 0--6   Rigid, Elastic, von Mises, Hill-48-79-90, BBC'03
   !         4: 1,2,3  No Hard, Linear, Ludwik-Nadai
   !         5: 1      No kinematic Hardening
   !         6: 1,2    Elastic, Elastic-Plastic
   !         7: 0,1    no-thermal, thermal
   !         8: 0      Associative
   ! PROPE: 1:Young 2:Poiss 3:Shear G  4:Bulk K  5:Density  6:Alpha  7:hcapc 8:hcond
   !        10-13: plane stress C (C11,C12,C22,C33)   [reference values]
   ! PROPS: curve numbers for
   !        1:Young 2:Poiss, 3:Cons0 4:Uniax  5:ConsN (or Hards)  6:Alpha  7:HCapc 8:HCond
   ! PROPP:  1-3: isotropic hardening parameters [reference values]
   !              (1)  Uniax  Cons0
   !              (2)  Hards  Uniax => e-bar
   !              (3)         Consn
   !           4: Deactivation time for plasticity
   !         5-7: Nothing yet
   !        8-27: Anisotropic plastic behavior
   !        8-11: s90,biuni,sxy,s45: strength ratios
   !       12-15: rm,r0,r45,r90    : Lankford parameters
   !       16-21: Hill coefficients for yield function
   !       22-27: Hill coefficients for potential function
   IMPLICIT NONE
   ! dummy arguments
   TYPE (mater), POINTER :: mat
   ! local variables
   INTEGER (kind=4) :: nc,k,i,npar
   REAL(kind=8) :: dens,e,nu,cons0,uniax,consn,alpha,hards,hcapc,hcond,nu2,death
   CHARACTER (len=mnam) :: cname(8)
   TYPE (curve), POINTER :: c
   LOGICAL :: pars(4)

   ! check Material Definitions

   IF(mat%matdef(2) == 0 ) mat%matdef(2) = 2  !isotropic
   IF(mat%matdef(2) > 2  ) CALL runend('MINP04: mtype =04 must be ISOTROPIC')
   IF(mat%matdef(2) == 1 ) mat%matdef(3) = 0  !rigid
   IF(mat%matdef(3)  > 6 ) CALL runend('MINP04: mtype =04 Invalid Yield function') !FR
   IF(mat%matdef(4)  > 3 ) CALL runend('MINP04: mtype =04 must be LINEAR or LUDWIK-NADAI')
   IF(mat%matdef(5)  > 1 ) CALL runend('MINP04: mtype =04 cannot use KIN-HARD ')
   mat%matdef(9)  = 13                    !size of array PROPE
   ! get memory for constants & read elastic properties
   ALLOCATE( mat%prope(mat%matdef(9)) )   !get memory for PROPE
   mat%prope = 0d0                        !initializes PROPE
   cname = ''                             !initializes names for curves
   !read elastic constants
   CALL listen('MINP04') !read a line
   dens  = getrea('DENS  ',1.D0,'!DENSITY ..........................')
   IF(mat%matdef(3)  > 0 ) THEN !not rigid
     e     = getrea('YOUNG ',0.D0,'!YOUNG MODULUS ...................')
     nu    = getrea('POISS ',0.D0,'!POISSON RATIO ...................')
     alpha = getrea('ALPHA ',0.D0,' Thermical Dilatation Coeff.......')
     IF( e     == 0d0 ) cname(1) = get_name('YOUNG ') !read curve name
     IF( nu    == 0d0 ) cname(2) = get_name('POISS ') !read curve name
     IF( alpha == 0d0 .AND. exists('ALPHA')) cname(6) = get_name('ALPHA ') !curve name
   END IF

   ! read plastic properties
   IF( mat%matdef(3) > 1 )THEN !plastic
     mat%matdef(10)  = 27               !size of array PROPP
     ALLOCATE( mat%propp(mat%matdef(10)) ) !get memory for PROPP
     mat%propp = 0d0                      !initializes PROPP
     CALL listen('MINP04') !read a line
     SELECT CASE (mat%matdef(4))
       CASE (1)  ! IHNONE: No hardening
         uniax = getrea('UNIAX ',0d0,' Yield Stress .....................')
         IF( uniax == 0d0 ) cname(4) = get_name('UNIAX ') !read curve name
       CASE (2)  ! IHLINE: linear Hardening
         uniax = getrea('UNIAX ',0d0,' Yield Stress .....................')
         hards = getrea('HARDS ',0d0,' Isotropic Hardening Par. .........')
         IF( uniax == 0d0 ) cname(4) = get_name('UNIAX ') !read curve name
         IF( hards == 0d0 ) cname(5) = get_name('HARDS ') !read curve name
       CASE (3)  ! IHLUDW: Ludwik Nadai
         cons0 = getrea('CONS0 ',1d0,'!Reference Yield Value ...........')
         uniax = getrea('UNIAX ',0d0,'!Elastic limit ...................')
         consn = getrea('CONSN ',1d0,'!Isotropic Hardening Exponent.....')
         IF( cons0 == 0d0 ) cname(3) = get_name('CONS0 ') !read curve name
         IF( uniax == 0d0 ) cname(4) = get_name('UNIAX ') !read curve name
         IF( consn == 0d0 ) cname(5) = get_name('CONSN ') !read curve name
     END SELECT
     death = getrea('DEATHT',HUGE(1d0),' End (Deactivation) Time ..........')

     ! read variables associated to plastic anisotropy
     IF( mat%matdef(3) /= 2 )THEN !not von MISES
       ! check existence of strength ratios
       ! strength ratios
       pars(1) = exists('SXY   ')
       pars(2) = exists('S45   ')
       pars(3) = exists('S90   ')
       pars(4) = exists('BIUNI ')
       npar = COUNT(pars)
       ! initializes
       mat%propp(8:11) = 0d0  !s90,biuni,sxy,s45
       IF( npar > 0 )THEN
         IF( pars(3) ) mat%propp(8)= getrea('S90   ',1d0,' Strength Ratio in local X2 direct.')
         IF( pars(4) ) mat%propp(9)= getrea('BIUNI ',1d0,' Biaxial strength ratio............')
         IF( pars(1) .AND. pars(2) )THEN
           CALL runend ('INMA02:Both SXY and S45 must not appear')
         ELSE IF( pars(1) )THEN
           mat%propp(10) = getrea('SXY   ',1d0,' Strength Ratio for local shear    ')
         ELSE IF( pars(2) ) THEN
           mat%propp(11) = getrea('S45   ',1d0,' Strength Ratio in local skew dir. ')
         END IF
         IF( pars(4) .AND. pars(1) ) &   !BIUNI and SXY
           mat%propp(11) =2d0/SQRT(3d0/mat%propp(10)**2+1d0/mat%propp(8)**2-1d0/mat%propp(8))
       END IF

       !plastic anisotropy in terms of Lankford coeffs.
       pars(1) = exists('RM    ')
       pars(2) = exists('R0    ')
       pars(3) = exists('R45   ')
       pars(4) = exists('R90   ')
       ! initializes
       mat%propp(12:15)  = 0d0

       IF( pars(1) .AND. (pars(2).OR.pars(3).OR.pars(4))) &
          CALL runend ('INMA02:Both RM and R?? must not appear')

       IF( pars(1) )THEN
         mat%propp(12) = getrea('RM    ',1d0,'!LANKFORD COEFFICIENT (average)....')
       ELSE IF( ANY (pars) )THEN
         mat%propp(13) = getrea('R0    ',1d0,'!LANKFORD COEFFICIENT IN ROLLING D.')
         mat%propp(14) = getrea('R45   ',1d0,'!LANKFORD COEFFICIENT DIR. 45 DEG. ')
         mat%propp(15) = getrea('R90   ',1d0,'!LANKFORD COEFFICIENT DIR. 90 DEG. ')
       END IF
       IF( exists('MHILL ')) &
         mat%propp(16) = getrea('MHILL ',2d0,'!HILL exponent ................... ')
     END IF
     CALL cohill( mat%propp, mat%matdef(3), mat%matdef(8) )  !comput HILL coefficients

   ELSE
     mat%matdef(10)  = 1                   !size of array PROPP
     ALLOCATE( mat%propp(mat%matdef(10)) ) !get memory for PROPP
     mat%propp = 0d0                      !initializes PROPP

   END IF
   !
   ! read thermal properties
   IF( mat%matdef(7) == 1 )THEN
     CALL listen('MINP04') !read a line (THERMAL_PROPERTIES should exist)
     hcapc = getrea('HCAPC ',0.D0,' Specific Heat Capacity ..........')
     hcond = getrea('HCOND ',0.D0,' Heat Conductivity ...............')
     IF( hcapc == 0d0 ) cname(7) = get_name('HCAPC ') !curve name
     IF( hcond == 0d0 ) cname(8) = get_name('HCOND ') !curve name
   END IF
   !
   ! read curve data for thermal dependant properties
   mat%matdef(11) = 8                    !size of array PROPS
   ALLOCATE( mat%props(mat%matdef(11)) )
   mat%props = 0d0 !initializes
   CALL listen('MINP04') !read a line
   IF( exists('CURVED') )THEN
     CALL rdmpcr(mat%chead,mat%ctail,nc,'TEMP  ','VALUE')
     mat%matdef(12) = nc
   ELSE
     backs = .TRUE. !BACKSPACE (ludat)
   END IF
   !
   ! Initialize elastic, plastic and thermal parameters
   IF(mat%matdef(3)  > 0 ) THEN !not rigid
     ! compute values for elastic materials
     IF( e == 0d0 )THEN
       CALL find_cr(mat%chead,c,cname(1),nc,k)
       mat%props(1) = k
       i = 1
       e     = inte_cr(c%val,c%np,0d0,i)
     END IF
     IF( nu == 0d0 )THEN
       CALL find_cr(mat%chead,c,cname(2),nc,k)
       mat%props(2) = k
       i = 1
       nu    = inte_cr(c%val,c%np,0d0,i)
     END IF
     mat%prope(1) = e
     mat%prope(2) = nu
     mat%prope(3)  = e/(1d0+nu)/2d0      !shear G
     mat%prope(4)  = e/(1d0-2d0*nu)/3d0  !bulk  K
     ! compute plane stress elasticity matrix
     nu2 = 1d0-nu*nu
     mat%prope(10) = e/nu2
     mat%prope(11) = mat%prope(10)*nu
     mat%prope(12) = mat%prope(10)
     mat%prope(13) = e/(2d0*(1d0+nu))
     IF( alpha == 0d0 )THEN
       CALL find_cr(mat%chead,c,cname(6),nc,k)
       mat%props(6) = k
       i = 1
       alpha = inte_cr(c%val,c%np,0d0,i)
     END IF
     mat%prope(6) = alpha
   ELSE
     mat%prope = 0d0
   END IF
   mat%prope(5)  = dens
   ! compute values for plastic materials
   IF( mat%matdef(3) > 1 )THEN
     IF( uniax == 0d0 )THEN
       CALL find_cr(mat%chead,c,cname(4),nc,k)
       mat%props(4) = k
       i = 1
       uniax = inte_cr(c%val,c%np,0d0,i)
     END IF
     SELECT CASE (mat%matdef(4))
       CASE (1)  ! IHNONE: Perfect plasticity
         mat%propp(1) = uniax
       CASE (2)  ! IHLINE: linear Hardening
         IF( hards == 0d0 )THEN
           CALL find_cr(mat%chead,c,cname(5),nc,k)
           mat%props(5) = k
           i = 1
           hards = inte_cr(c%val,c%np,0d0,i)
         END IF
         mat%propp(1) = uniax
         mat%propp(2) = hards
       CASE (3)  ! IHLUDW: Ludwik Nadai
         IF( cons0 == 0d0 )THEN
           CALL find_cr(mat%chead,c,cname(3),nc,k)
           mat%props(3) = k
           i = 1
           cons0 = inte_cr(c%val,c%np,0d0,i)
         END IF
         IF( consn == 0d0 )THEN
           CALL find_cr(mat%chead,c,cname(5),nc,k)
           mat%props(5) = k
           i = 1
           consn = inte_cr(c%val,c%np,0d0,i)
         END IF
         mat%propp(1) = cons0
         mat%propp(3) = consn
     END SELECT
     mat%propp(4) = death
   END IF
  ! compute values for thermo-elasto/plastic materials
   IF( mat%matdef(7) == 1 )THEN
     IF( hcapc == 0d0 )THEN
       CALL find_cr(mat%chead,c,cname(7),nc,k)
       mat%props(7) = k
       i = 1
       hcapc = inte_cr(c%val,c%np,0d0,i)
     END IF
     IF( hcond == 0d0 )THEN
       CALL find_cr(mat%chead,c,cname(8),nc,k)
       mat%props(8) = k
       i = 1
       hcond = inte_cr(c%val,c%np,0d0,i)
     END IF
     mat%prope(7) = hcapc
     mat%prope(8) = hcond
   END IF
   RETURN
   END SUBROUTINE minp_04
