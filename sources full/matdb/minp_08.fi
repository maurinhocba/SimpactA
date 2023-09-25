   SUBROUTINE minp_08 ( mat )
   !            --------------------------------
   !            input material properties for
   !            visco-elastic material
   !            elastic: linear isotropic
   !            viscous: linear isotropic
   !            --------------------------------
   ! MATDEF: 1: 3      3-D
   !         2: 2      Isotropic
   !         3: 0--6   Elastic
   !         4: 1      No Hard
   !         5: 1      No kinematic Hardening
   !!         6: 1      Elastic
   !         7: 0,1    NO/YES volumetric viscosity
   !!         8: 0      Associative
   ! PROPE: 1:Young 2:Poiss 3:Shear G  4:Bulk K  5:Density  6:Alpha
   !         7-10: plane stress C (C11,C12,C22,C33)   [reference values]
   !        11:xitas 12:xitav
   ! PROPS: nothing
   ! PROPP: nothing
   IMPLICIT NONE
   ! dummy arguments
   TYPE (mater), POINTER :: mat
   ! local variables
   REAL(kind=8) :: nu2

   ! check Material Definitions

   IF(mat%matdef(2) == 0 ) mat%matdef(2) = 2  !isotropic
   IF(mat%matdef(2) > 2  )CALL runend('MINP08: mtype =08 must be ISOTROPIC')
   IF(mat%matdef(2) == 1 )mat%matdef(3) = 0  !rigid
   IF(mat%matdef(3)  > 1 )CALL runend('MINP08: mtype =08 can be Elastic only')
   IF(mat%matdef(4)  > 1  .OR. mat%matdef(5) > 1 ) &
     CALL runend('MINP08: mtype =08 No Hardening function allowed')
   IF(mat%matdef(6) /= 1 )mat%matdef(6) = 1  !viscous-elastic

   mat%matdef(9)  = 12                    !size of array PROPE
   ! get memory for constants & read elastic properties
   ALLOCATE( mat%prope(mat%matdef(9)) )   !get memory for PROPE
   !read elastic constants
   CALL mat_ela(mat%prope,mat%matdef)     !read elastic constants
   IF( mat%matdef(2) /= 1 )THEN ! not rigid material
     ! compute plane stress elasticity matrix
     nu2 = 1d0-mat%prope(2)*mat%prope(2)                   !1-nu^2
     mat%prope(7)  = mat%prope(1)/nu2                      !E/(1-nu^2)
     mat%prope(8)  = mat%prope(7)*mat%prope(2)             !*nu
     mat%prope(9)  = mat%prope(7)                          !==C11
     mat%prope(10) = mat%prope(1)/(2d0*(1d0+mat%prope(2))) !E/2(1+nu)
   END IF
   ! read viscous  properties
   mat%prope(11) = getrea('VISCO ',0d0,'!Viscosity parameter ...............')
   !IF( mat%prope(11) < 1d0)CALL runend('MINP08: Viscosity can not be too low')
   !mat%prope(12) = mat%prope(11)/mat%prope(3)/32d0  !minimum time increment
   mat%prope(12) = 0d0
   mat%matdef(10)  = 1                   !size of array PROPP
   ALLOCATE( mat%propp(mat%matdef(10)) ) !get memory for PROPP
   mat%propp = 0d0                       !initializes PROPP
   mat%matdef(11) = 1                    !size of array PROPS
   ALLOCATE( mat%props(mat%matdef(11)) )
   mat%props = 0d0 !initializes
   RETURN
   END SUBROUTINE minp_08
