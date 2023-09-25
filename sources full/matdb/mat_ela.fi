 SUBROUTINE mat_ela (epr, def)
 !
 !  read standard elasticity constants
 !  this is a general routine and can be adapted to other needs
 !
 IMPLICIT NONE
 ! dummy arguments
 INTEGER (kind=4), INTENT(IN OUT) :: def(:) ! (1)Material dimension
                                            ! (2)Elasticity model

 REAL (kind=8), INTENT(OUT) :: epr(:)   ! properties array
 ! 1:E     (Young modulus)
 ! 2:mu    (Poisson ratio)
 ! 3:G     (Shear modulus)
 ! 4:K     (Bulk modulus)
 ! 5:dens  (density)
 ! 6:alpha (temperature dilatation coefficient)
 ! 7..NN according to material

 ! local variables
 REAL (kind=8) :: e,nu,g,k,lambda


 CALL listen('MATELA') !read a line

 SELECT CASE ( def(2) )  !according to the type of material

 CASE (1)   ! Rigid material
   epr(5) = getrea('DENS  ',1.D0,'!DENSITY ..........................')

 CASE (2)   ! Linear elastic isotropic material
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
   ELSE IF( exists('YOUNG ') .AND. exists('SHEAR ')) THEN
     e  = getrea('YOUNG ',0.D0,'!YOUNG MODULUS ....................')
     g  = getrea('SHEAR ',0.D0,'!SHEAR MODULUS ....................')
     nu= e/2d0/g-1d0
     k = e/(1d0-2d0*nu)/3d0
     WRITE (lures,"(9X,'COMPUTED POISSON RATIO ............... POISS  = ',F14.7)",ERR=9999) nu
   ELSE IF ( exists('MU') .AND. exists('LAMBDA')) THEN
     g = getrea('MU',0.D0,'!Mu Lame Parameter ................')
     lambda = getrea('LAMBDA',0.D0,'!Lambda Lame Parameter.............')
     k = lambda - 2d0/3d0
     nu = (3d0*k-2d0*g)/(6d0*k+2d0*g)
     e  = 9d0*g*k/(g+3d0*k)
    WRITE (lures,"(9X,'COMPUTED YOUNG MODULUS ............... YOUNG  = ',G14.5)",ERR=9999) e
    WRITE (lures,"(9X,'COMPUTED POISSON RATIO ............... POISS  = ',F14.7)",ERR=9999) nu
   ELSE
     CALL runend('MATELA:YOUNG & POISS or SHEAR & BULK compulsory')
   END IF

   epr(1) = e
   epr(2) = nu
   epr(3) = g
   epr(4) = k

   epr(5) = getrea('DENS  ',1.D0,'!DENSITY ..........................')
   epr(6) = getrea('ALPHAT',0.D0,' Thermical Dilatation Coeff........')


 CASE (3)   ! linear elastic Orthotropic materials

   IF (exists('NU12  ')) THEN
     def(8) = 0                  !set a flag
     epr(7)  = getrea('E1    ',0d0,'!YOUNG MODULUS IN DIR. 1 ..........')
     epr(8)  = getrea('E2    ',0d0,' YOUNG MODULUS IN DIR. 2 ..........')
     epr(9)  = getrea('E3    ',0d0,' YOUNG MODULUS IN DIR. 3 ..........')

     epr(10) = getrea('NU12  ',0d0,'!POISSON RATIO IN PLANE 1-2 -e2/e1.')
     epr(11) = getrea('NU13  ',0d0,' POISSON RATIO IN PLANE 1-3 -e3/e1.')
     epr(12) = getrea('NU23  ',0d0,' POISSON RATIO IN PLANE 2-3 -e3/e2.')

     IF( epr(8) /= 0d0 )THEN   !default value
       nu = epr(8)*epr(10)/epr(7)                        !nu21
       g  = 1d0/((1d0+epr(10))/epr(7) + (1d0+nu)/epr(8)) !g12
     ELSE
       g = 0d0
     END IF
     epr(13) = getrea('G12   ',g,' Shear Modulus in Plane 1-2 .......')
     IF( epr(13) == 0d0 ) epr(13) = g

     IF( epr(9) /= 0d0 )THEN   !default value
       nu = epr(9)*epr(11)/epr(7)                        !nu31
       g  = 1d0/((1d0+epr(11))/epr(7) + (1d0+nu)/epr(9)) !g13
     ELSE
       g = 0d0
     END IF
     epr(14) = getrea('G13   ',g,' Shear Modulus in Plane 1-3 .......')
     IF( epr(14) == 0d0 ) epr(14) = g

     IF( epr(8) /= 0d0 .AND. epr(9) /= 0d0)THEN   !default value
       nu = epr(9)*epr(12)/epr(8)                        !nu32
       g  = 1d0/((1d0+epr(12))/epr(8) + (1d0+nu)/epr(9)) !g23
     ELSE
       g = 0d0
     END IF
     epr(15) = getrea('G23   ',g,' Shear Modulus in Plane 2-3 .......')
     IF( epr(15) == 0d0 ) epr(15) = g

     epr(1) = epr(7)                         !Young modulus
     IF( epr(8) > 0 ) THEN                   !E2 present
       epr(2) = epr(1)*epr(10)/epr(8)        !Poisson ratio 21
       epr(3) = epr(1)/(1d0+epr(2))/2d0
       IF( epr(9) > 0 ) &                    !E3 present
         epr(2) = (epr(2) + epr(1)*epr(11)/epr(9))/2d0  !average between 21 y 31
     ELSE
       epr(2) = epr(10)                      !use 12 as isotropic
     END IF

   ELSE  ! elasticity matrix
     def(8) = 1    !set a flag

     epr ( 7)= getrea('D1111 ',1d0,'!D1111 COMPONENT OF D TENSOR ......')
     epr ( 8)= getrea('D1122 ',0d0,'!D1122 COMPONENT OF D TENSOR ......')
     epr ( 9)= getrea('D2222 ',0d0,'!D2222 COMPONENT OF D TENSOR ......')
     epr (10)= getrea('D1212 ',0d0,'!D1212 COMPONENT OF D TENSOR ......')
     epr (11)= getrea('D3333 ',0d0,' D3333 COMPONENT OF D TENSOR ......')
     epr (12)= getrea('D1133 ',0d0,' D1133 COMPONENT OF D TENSOR ......')
     epr (13)= getrea('D2233 ',0d0,' D2233 COMPONENT OF D TENSOR ......')
     epr (14)= getrea('D1313 ',0d0,' D1313 COMPONENT OF D TENSOR ......')
     epr (15)= getrea('D2323 ',0d0,' D2323 COMPONENT OF D TENSOR ......')

     epr(1) = 0d0  !flag

   END IF

   epr(5) = getrea('DENS  ',1.D0,'!DENSITY ..........................')
   epr(6) = getrea('ALPHAT',0.D0,' Thermical Dilatation Coeff........')

 CASE (5)   ! Elastic curve by points
   ! read curve for elastic stress-strain definition

 CASE DEFAULT
 END SELECT

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE mat_ela
