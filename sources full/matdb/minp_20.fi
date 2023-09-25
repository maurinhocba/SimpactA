 SUBROUTINE minp_20 ( mat )
 !            --------------------------------
 !            input material properties
 !            for a curve dependent material
 !            Used by SPOT
 IMPLICIT NONE
 ! dummy arguments
 TYPE (mater), POINTER :: mat
 ! local variables
 INTEGER (kind=4) :: nc

 mat%matdef(9)  = 6                    !size of array PROPE
 ! get memory for constants & read elastic properties
 ALLOCATE( mat%prope(mat%matdef(9)) )  !get memory for PROPE
 mat%prope = 0d0                       !initializes PROPE
 !read basic elastic constants if present
 CALL listen('MINP20') !read a line
 IF( exists('YOUNG ')) mat%prope(1) = &
   getrea('YOUNG ',0.D0,' YOUNG MODULUS ...................')
 IF( exists('POISS ')) mat%prope(2) = &
   getrea('POISS ',0.D0,' POISSON RATIO ...................')
 mat%prope(3)= mat%prope(1)/(1d0+mat%prope(2))/2d0
 IF( exists('SHEAR') ) mat%prope(3)=  &
   getrea('SHEAR ',0.D0,' SHEAR MODULUS ...................')
 IF( exists('BULK') ) mat%prope(4)=  &
   getrea('BULK  ',0.D0,' BULK MODULUS ....................')
 mat%prope(4)= mat%prope(1)/(1d0-2d0*mat%prope(2))/3d0
 IF( exists('DENS  ')) mat%prope(5) = &
   getrea('DENS  ',1.D0,' DENSITY .........................')

 mat%matdef(10) = 1                    !size of array PROPP
 mat%matdef(11) = 1                    !size of array PROPS
 ALLOCATE( mat%propp(mat%matdef(10)), mat%props(mat%matdef(11)) )

 CALL listen('MINP20') !read a line
 IF( exists('CURVED') )THEN     !IF curve defined material
   CALL rdmpcr(mat%chead,mat%ctail,nc,'VALUEX','VALUEY')  !read curves
 ELSE
   backs = .TRUE.
   nc = 0
 END IF
 mat%matdef(12) = nc  !number of curves

 RETURN
 END SUBROUTINE minp_20
