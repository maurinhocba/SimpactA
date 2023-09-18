 SUBROUTINE minp_30 ( mat )
 !            --------------------------------
 !            input material properties for
 !            USER DEFINED material
 !            --------------------------------
 ! MATDEF: 1: 1,2,3  1D,2D,3D
 !         2: 2,3    Isotropic, Orthotropic
 !         3: 1,3    Elastic,Plastic-Isotropic, Plastic-Orthotropic
 !         4: 1:     No isotropic Hardemomg, User defined isotropic hardening
 !         5: 1:     No kinematic Hardening, User defined kinematic hard
 !         6: 1,2    Elastic, Elastic-Plastic (it can be viscous yet)
 !         7: 0:8    Aditional Internal variables
 !         8: 30-39  Submodel
 !         9: 10:    size of array PROPE (default = 10)
 !        10: 0:     size of array PROPP (default =  0)
 !        11: 0:     size of array PROPS (default =  0)
 !
 ! PROPE: 1:Young 2:Poiss 3:Shear G  4:Bulk K  5:Density   6:Alpha
 !        7-10: plane stress C (C11,C12,C22,C33)
 ! PROPS: 1:    User defined
 ! CURVES:      User defined
 IMPLICIT NONE
 ! dummy arguments
 TYPE (mater), POINTER :: mat
 ! local variables
 REAL(kind=8) :: e,nu,g,k,nu2
 INTEGER (kind=4) :: i,nprop

 ! check Material Definitions

 CALL listen('MATELA')                  !read definitions

 mat%matdef(2) = 2              !isotropic
 IF( exists('MODEL ',i) )THEN
   IF( words(i+1)(1:midn) == 'ORTHOT' )mat%matdef(2) = 3
 END IF
 IF( mat%matdef(2) == 2 )THEN
   WRITE(lures,"(9x,'Elasticity model ',21('.'),' MODEL  =    ISOTROPIC')",ERR=9999)
 ELSE
   WRITE(lures,"(9x,'Elasticity model ',21('.'),' MODEL  =  ORTHOTROPIC')",ERR=9999)
 END IF
 mat%matdef(3) = 2              !plastic as default
 mat%matdef(4) = 1              !
 mat%matdef(5) = 1              !
 mat%matdef(6) = 1              !
 mat%matdef(7) = getint('NVARE ',1,  ' Additional internal variables ....')
 mat%matdef(8) = getint('SUBMO ',30, ' Sub model ........................')
 nprop         = getint('NPROP ',0,  '!Size of array PROPS...............')

 IF( mat%matdef(8) < 30 .OR. mat%matdef(8) > 39 )mat%matdef(8) = 30
 ! get memory for constants & read elastic properties
 mat%matdef(9) = 10                     !basic properties
 ALLOCATE( mat%prope(mat%matdef(9)) )   !get memory for PROPE
 mat%prope = 0d0                        !initializes PROPE
 CALL listen('MATELA')                  !read elastic properties
 IF( exists('YOUNG ') .AND. exists('POISS ')) THEN
   e  = getrea('YOUNG ',0.D0,'!YOUNG MODULUS ....................')
   nu = getrea('POISS ',0.D0,'!POISSON RATIO ....................')
   g = e/(1d0+nu)/2d0
   k = e/(1d0-2d0*nu)/3d0
 ELSE IF ( exists('SHEAR') .AND. exists('BULK')) THEN
   g  = getrea('SHEAR ',0.D0,'!SHEAR MODULUS ....................')
   k  = getrea('BULK  ',0.D0,'!BULK MODULUS .....................')
   nu = (3d0*k-2d0*g)/(6d0*k+2d0*g)
   e  = 9d0*g*k/(g+3d0*k)
  WRITE (lures,"(9X,'COMPUTED YOUNG MODULUS ............... YOUNG  = ',G14.5)",ERR=9999) e
  WRITE (lures,"(9X,'COMPUTED POISSON RATIO ............... POISS  = ',F14.7)",ERR=9999) nu
 ELSE
   CALL runend('MATELA:YOUNG & POISS or SHEAR & BULK compulsory')
 END IF

 mat%prope(1) = e
 mat%prope(2) = nu
 mat%prope(3) = g
 mat%prope(4) = k

 mat%prope(5) = getrea('DENS  ',1.D0,'!DENSITY ..........................')
 mat%prope(6) = getrea('ALPHA ',0.D0,' Thermical Dilatation Coeff........')

 IF( .NOT.exists('C11   ') )THEN
   ! compute plane stress elasticity matrix
   nu2 = 1d0-nu*nu
   mat%prope(7)  = e/nu2                      !E/(1-nu^2)
   mat%prope(8)  = e/nu2*nu                   !*nu
   mat%prope(9)  = mat%prope(7)               !==C11
   mat%prope(10) = e/(2d0*(1d0+nu))           !E/2(1+nu)
 ELSE
   mat%prope(7)  = getrea('C11   ',0.D0,'!C11 coefficient of plane stress...')
   mat%prope(8)  = getrea('C12   ',0.D0,'!C12 coefficient of plane stress...')
   mat%prope(9)  = getrea('C22   ',mat%prope(7),' C22 coefficient of plane stress...')
   mat%prope(10) = getrea('C33   ',0.d0,'!C33 coefficient of plane stress...')
 END IF

 mat%matdef(10) = 1                    !set to 1
 ALLOCATE( mat%propp(mat%matdef(10)) ) !get memory for PROPP
 mat%propp = 0d0                       !initializes PROPP

 ! get memory for constants & read other properties
 CALL listen('MATPRO')                 !read more properties
 IF( nnpar > nprop )THEN
    WRITE(lures,"(' WARNING! Number of parameters read (',i5,')',/, &
          &       '              is greater that NPROP (',i5,')')",ERR=9999) nnpar,nprop
    WRITE(lures,"('          NPROP increased to  (',i5,')')",ERR=9999) nnpar
    nprop = nnpar
 ELSE IF (nnpar < nprop )THEN
    WRITE(lures,"(' WARNING! Number of parameters read (',i5,')',/, &
          &       '              is smaller that NPROP (',i5,')')",ERR=9999) nnpar,nprop
 END IF
 mat%matdef(11)= nprop
 ALLOCATE( mat%props(nprop) ) !get memory for PROPS
 mat%props = 0d0                       !initializes PROPS
 WRITE(lures,"('   Material properties (PROPS)')",ERR=9999)
 DO i=1,nprop
   mat%props(i) = param(i)
   WRITE(lures,"(' PROPS(',i2,') = ',e15.6)",ERR=9999) i, mat%props(i)
 END DO
 ! Read curves
 CALL listen('MATCUR')                 !read curves
 IF( exists('CURVED') )THEN     !IF curve defined material
   CALL rdmpcr(mat%chead,mat%ctail,i,'VALUEX','VALUEY')  !read curves
   mat%matdef(12) = i           !keep number of curves
 ELSE
   !BACKSPACE (ludat)
   mat%matdef(12) = 0
   backs = .TRUE.
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE minp_30
