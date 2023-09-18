 SUBROUTINE minp_06 ( mat)
 !            --------------------------------
 !            input material properties for
 !            hyper-elastic rubber type material
 !            --------------------------------
 ! MATDEF: 1: 2,3    2-D, 3-D
 !         2: 2      Isotropic
 !         3: 1      Elastic
 !         8: 1-17   Model
 ! PROPE: 1:Max Young or 0D0       5:Density   6:Alpha
 !        7-27 : different parameters according to the model
 ! PROPP: nothing
 ! PROPS: nothing
 IMPLICIT NONE
 ! dummy arguments
 TYPE (mater), POINTER :: mat
 REAL(Kind=8) :: nu2

 ! check Material Definitions

 IF(mat%matdef(2) == 0 ) mat%matdef(2) = 2
 IF(mat%matdef(2) /= 2 )CALL runend('MINP06: mtype=6  must be ISOTROPIC   ')
 IF(mat%matdef(3)  > 3 )CALL runend('MINP06: mtype=6  cannot use HILL79/90')
 IF(mat%matdef(4)  > 2 )CALL runend('MINP06: mtype=6  can use IHLINE only ')
 IF(mat%matdef(5)  > 1 )CALL runend('MINP06: mtype=6  cannot use KIN-HARD ')

 ! get memory for constants & read elastic properties
 mat%matdef(9)  = 22                        !size of array PROPE
 ALLOCATE( mat%prope(mat%matdef(9)) )       !get memory for PROPE
 mat%prope = 0d0                            !initializes PROPE
 CALL mat_rubber(mat%prope,mat%matdef)      !read elastic constants
 ! compute plane stress elasticity matrix
 nu2 = 1d0-mat%prope(2)*mat%prope(2)                   !1-nu^2
 mat%prope(19) = mat%prope(1)/nu2                      !E/(1-nu^2)
 mat%prope(20) = mat%prope(19)*mat%prope(2)            !*nu
 mat%prope(21) = mat%prope(19)                         !==C11
 mat%prope(22) = mat%prope(1)/(2d0*(1d0+mat%prope(2))) !E/2(1+nu)
 !  reserve space
 IF( mat%matdef(3) /= 1 )THEN
   mat%matdef(10) = 16
   ALLOCATE( mat%propp(mat%matdef(10)) )
   IF(mat%matdef(2) == 0 )mat%matdef(2) = 2  !isotropic
   IF(mat%matdef(3) == 0 )mat%matdef(3) = 2  !mises
   IF(mat%matdef(4) == 0 )mat%matdef(4) = 2  !no isotropic hardening
   mat%propp = 0d0
   CALL mat_pla(mat%propp,mat%matdef)
   IF(mat%matdef(3) > 2 )CALL runend ('MINP06: rubbers are isotropic ')
   IF(mat%matdef(4) > 2 )CALL runend ('MINP06: Linear isot. hard for rubbers')
 ELSE                            !for elastic material
   mat%matdef(10)  = 1                   !set size to 1
   ALLOCATE( mat%propp(mat%matdef(10)) )
   mat%propp = 0d0
 END IF
 mat%matdef(11) = 1                         !size of array PROPS
 ALLOCATE( mat%props(mat%matdef(11)) )      !get memory for PROPS
 !and compute elastic matrices (nothing)

 RETURN
 END SUBROUTINE minp_06

 SUBROUTINE mat_rubber(epr,def)
 !            --------------------------------
 !            input material properties for
 !            hyper elastic rubber type material
 !            -- material can be used for:
 !               a) 3d continuum elements
 !               b) 3d shell elements     (no yet)
 !               b) 2d continuum elements (no yet)
 !
 IMPLICIT NONE
 ! dummy arguments
 INTEGER (kind=4), INTENT(IN OUT) :: def(:) ! (1)Material dimension
                                            ! (2)Elasticity model
                                            ! (6)Submodel

 REAL (kind=8), INTENT(OUT) :: epr(:)   ! properties array
 ! 1:E  2:mu   3:G  4:K  5:dens 6:alpha 7..NN according to material

 ! local variables
 INTEGER (kind=4) :: model,nterms,i,j
 REAL (kind=8) :: e,nu,g,k,lm


 CALL listen('MATELA') !read a line
 IF( exists('HYPERE') ) THEN
   nterms = getint('NTERMS',   1,' NUMBER OF terms in expansion .....')
   IF( exists('ARRUDA'))THEN
     model = 1
   ELSE IF( exists('MOONEY'))THEN
     model = 2
   ELSE IF( exists('NEOHOO'))THEN
     model = 3
   ELSE IF( exists('OGDEN '))THEN
     model = 3 + nterms
   ELSE IF( exists('POLYNO'))THEN
     model = 6 + nterms
   ELSE IF( exists('REDUCE'))THEN
     model = 9 + nterms
   ELSE IF( exists('YEOH  '))THEN
     model = 14
     nterms = 3
   END IF

   IF(model ==0 )CALL runend('MTYPE_6: a model, e.g. OGDEN, must be selected ')

   CALL listen('MATELA') !read a line
   epr(5) = getrea('DENS  ',1.D0,'!DENSITY ..........................')
   epr(6) = getrea('ALPHAT',0.D0,' Thermical Dilatation Coeff........')

   SELECT CASE (model)

   CASE (1) !Arruda-Boyce
     epr(7) = getrea('MU    ',1.D0,'!Shear Coefficient ................')
     epr(8) = getrea('LAMBDA',1.D0,'!LAMBDA_m Coefficient .............')
     epr(9) = getrea('COMPRE',0.D0,' Compresibility ...................')
     lm = 1d0/epr(8)
     epr(10) = 1d0        /     2d0 * epr(7)
     epr(11) = lm**2      /    10d0 * epr(7)
     epr(12) = lm**4* 11d0/   350d0 * epr(7)
     epr(13) = lm**6* 19d0/  1750d0 * epr(7)
     epr(14) = lm**8*519d0/134750d0 * epr(7)
     g = epr(7)

   CASE (2) !Mooney-Rivlin
     epr(7) = getrea('C10   ',1.D0,'!C10 Coefficient ..................')
     epr(8) = getrea('C01   ',1.D0,'!C01 Coefficient ..................')
     epr(9) = getrea('COMPRE',0.D0,' Compresibility ...................')
     g = 2.d0*(epr(7) + epr(8))

   CASE (3) !Neo-Hooke
     epr(7) = getrea('C10   ',1.D0,'!C10 Coefficient ..................')
     epr(8) = 0d0
     epr(9) = getrea('COMPRE',0.D0,' Compresibility ...................')
     g = 2.d0*epr(7)

   CASE (4:6) !Ogden

     epr(7) = getrea('MU1   ',1.D0,'!Ogden First factor................')!/2d0
     epr(8) = getrea('ALPHA1',1.D0,'!Ogden First exponent..............')
     !epr(16) = getrea('COMPR1',0.D0,' Compresibility (first term).......')
     epr(16) = getrea('COMPRE',0.D0,' Compresibility ...................')
     IF( model > 4 )THEN
       epr( 9)= getrea('MU2   ',0.D0,'!Ogden Second factor...............')!/2d0
       epr(10)= getrea('ALPHA2',1.D0,'!Ogden Second exponent.............')
       !epr(17) = getrea('COMPR2',0.D0,' Compresibility (second term)......')
     END IF
     IF( model == 6 )THEN
       epr(11)= getrea('MU3   ',0.D0,'!Ogden Third factor................')!/2d0
       epr(12)= getrea('ALPHA3',1.D0,'!Ogden Third exponent..............')
       !epr(18)= getrea('COMPR3',0.D0,' Compresibility (third term).......')
     END IF
     g = (epr(7)*epr(8) + epr(9)*epr(10)  + epr(11)*epr(12) )!*2d0

   CASE (7:9,10:12) !Polynomial model or reduced polynomial mode

     epr(7) = getrea('C10   ',1.D0,'!Coefficient i=1, j=0 .............')
     epr(8) = getrea('C01   ',1.D0,' Coefficient i=0, j=1 .............')
     epr(16) = getrea('COMPR1',0.D0,' Compresibility (first term).......')
     IF( nterms > 1 )THEN
       epr( 9)= getrea('C20   ',1.D0,'!Coefficient i=2, j=0 .............')
       epr(10)= getrea('C11   ',1.D0,' Coefficient i=1, j=1 .............')
       epr(11)= getrea('C02   ',1.D0,' Coefficient i=0, j=2 .............')
       epr(17)= getrea('COMPR2',0.D0,' Compresibility (second term)......')
     END IF
     IF(  nterms == 3 )THEN
       epr(12)= getrea('C30   ',1.D0,'!Coefficient i=3, j=0 .............')
       epr(13)= getrea('C21   ',1.D0,' Coefficient i=2, j=1 .............')
       epr(14)= getrea('C12   ',1.D0,' Coefficient i=1, j=2 .............')
       epr(15)= getrea('C03   ',1.D0,' Coefficient i=0, j=3 .............')
       epr(18)= getrea('COMPR3',0.D0,' Compresibility (third term).......')
     END IF
     g = 2.d0*(epr(7)+epr(8))

   CASE (14) !Yeoh

     epr( 7) = getrea('C10   ',1.D0,'!Coefficient i=1, j=0 .............')
     epr( 9) = getrea('C20   ',1.D0,'!Coefficient i=2, j=0 .............')
     epr(12) = getrea('C30   ',1.D0,'!Coefficient i=3, j=0 .............')
     epr(16) = getrea('COMPR1',0.D0,' Compresibility (first term).......')
     epr(17) = getrea('COMPR2',0.D0,' Compresibility (second term)......')
     epr(18) = getrea('COMPR3',0.D0,' Compresibility (third term).......')

     g = 2.d0*epr(7)
   END SELECT

   i = 8
   IF( model > 3 )i=15
   DO j=1,nterms
     IF(ABS(epr(i+j)) <  1.d-10) then
       epr(i+j) = (0.1d0/g)**j
       WRITE(lures,*) '*WARNING in hyperelastics: default value was'
       WRITE(lures,*) '         used for compressibility coefficients'
       WRITE(lures,"('   D',i1,' = ',e11.4)") j,epr(i+j)
     END IF
   END DO

   ! consistent Young modulus and Poisson ratio
   k  = 1d0/epr(i+1)
   e  = 9d0*k*g/(3d0*k+g)
   nu = e/2d0/g - 1d0

   epr(1) = e
   epr(2) = nu
   epr(3) = g
   epr(4) = k

 ELSE IF( exists('HYPERF') )THEN
   nterms = getint('NTERMS',   1,' NUMBER OF terms in expansion .....')
   model  = 14 + nterms

   CALL listen('MATELA') !read a line
   epr(5) = getrea('DENS  ',1.D0,'!DENSITY ..........................')
   epr(6) = getrea('ALPHAT',0.D0,' Thermical Dilatation Coeff........')

   epr(7) = getrea('MU1   ',1.D0,'!First factor .....................')
   epr(8) = getrea('ALPHA1',1.D0,'!First exponent....................')
   epr(13) = getrea('POISS1',0.D0,'!Poisson ratio .(first term).......')
   IF( model > 15)THEN
     epr( 9)= getrea('MU2   ',0.D0,'!Second factor ....................')
     epr(10)= getrea('ALPHA2',1.D0,'!Second exponent...................')
     epr(14)=getrea('POISS2',0.D0,'!Poisson ratio .(second term)......')
   END IF
   IF( model == 17)THEN
     epr(11)= getrea('MU3   ',0.D0,'!Third factor        ..............')
     epr(12)= getrea('ALPHA3',1.D0,'!Third exponent      ..............')
     epr(15)= getrea('POISS3',0.D0,'!Poisson ratio  (third term).......')
   END IF
   g = epr(7)*epr(8) + epr(9)*epr(10)  + epr(11)*epr(12)
   ! consistent Young modulus and Poisson ratio
   nu = epr(i)
   e  = 2d0*g*(1+nu)
   k  = e/(1d0-2d0*nu)/3d0

   epr(1) = e
   epr(2) = nu
   epr(3) = g
   epr(4) = k

 ELSE
   CALL runend('MTYPE_6: HyperElastic or HyperFoam keyword is compulsory')
 END IF

 def(8) = model

 RETURN
 END SUBROUTINE mat_rubber
