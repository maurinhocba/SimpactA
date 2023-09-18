 SUBROUTINE sc_shel3( sec )
 !            --------------------------------
 !            input material properties for multi-layered lamina
 !            Total Lagrangean BST elements (NBST,LBST,CBST)
 !  iprop: 1 :nlayr
 !         2 :nvars:  number of internal variables for the whole laminate
 !         3 :0  elastic lamina,  1:inelastic
 !         4 :1 if transverse shear is considered in constitutive equation
 !    5:4+nl :associated material label to each layer
 !5+nl:4+2nl :pointer to each component internal variables
 !
 !  rprop  1 :thick
 !         2 :mean young modulus
 !         3 :mean density
 !         4 :elastic generalized strain threshold
 !         5 :surface density
 !      6-26 :integrated elasticity matrix
 !      For each layer (npt values)
 !          1: layer thickness
 !          2: layer angle
 !          3: layer z coordinate
 !        4-8: rotation coefficients
 !          9: shear factor

 IMPLICIT NONE
 TYPE (section), POINTER :: sec

 !local variables
 INTEGER (kind=4), PARAMETER :: npl = 9     !number of real values per layer
 REAL (kind=8), PARAMETER :: factor = 1.7453292519943295769236907684886d-2
 INTEGER (kind=4), PARAMETER :: matsa(2) = (/1,5/) ! material types allowed for this section
 INTEGER (kind=4) matno,nlayr,ith,il,sizer,sizei,ir,ii,j,k,nvars,kvar,imat
 LOGICAL :: found,zauto,readm
 REAL(kind=8) :: thick,youn,dens,thickl,de(21),angle,zl,ce(6),sg,cg,celoc(6),minstr,aux
 REAL(kind=8), ALLOCATABLE :: rprop(:,:)
 TYPE (mater), POINTER :: mat

 sec%secty = 13         !shell type 3
 sec%mabas = -1         ! flag to indicate that no unique associated material exist

 CALL listen('SHELL3')  !read a line

 readm = exists('READMA')  !read elastic matrix ?
 IF( readm )THEN
   sec%secdef = (/ 4, 26, 0, 0, 0 /)
   ALLOCATE (sec%iprop(4), sec%rprop(26) )
   sec%iprop = 0
   sec%rprop = 0d0
   sec%rprop(1) = getrea('THICK ',0.D0,'!Thickness ........................')
   sec%rprop(4) = 1.D3   ! Threshold
   sec%rprop(5) =getrea('DENSI ',0d0,'!Surface density ..................')  !surface density
   sec%rprop( 6)=getrea('C1111 ',0d0,'!Membrane Coeff dir 1 .............')
   sec%rprop( 7)=getrea('C1122 ',0d0,'!Membrane Coeff dir 1-2............')
   sec%rprop( 8)=getrea('C1112 ',0d0,' Membrane Coupling  11-12..........')
   sec%rprop( 9)=getrea('C2222 ',0d0,'!Membrane Coeff dir 2 .............')
   sec%rprop(10)=getrea('C2212 ',0d0,' Membrane Coupling  22-12..........')
   sec%rprop(11)=getrea('C1212 ',0d0,'!Membrane Shear Coeff 1-2..........')
   sec%rprop(12)=getrea('D1111 ',0d0,'!Bending  Coeff dir 1 .............')
   sec%rprop(13)=getrea('D1122 ',0d0,'!Bending  Coeff dir 1-2............')
   sec%rprop(14)=getrea('D1112 ',0d0,' Bending Coupling  11-12 ..........')
   sec%rprop(15)=getrea('D2222 ',0d0,'!Bending  Coeff dir 2 .............')
   sec%rprop(16)=getrea('D2212 ',0d0,' Bending Coupling  22-12 ..........')
   sec%rprop(17)=getrea('D1212 ',0d0,'!Bending  Shear Coeff 1-2..........')

 ELSE
   thick = getrea('THICK ',0.D0,' Thickness ........................')
   IF( thick > 0d0 ) THEN
     ith=1 ! thickl will be input as a percentage of thick
   ELSE
     ith=0 ! thickl will be input as a REAL DIMENSION
   END IF
   !minth = getrea('MINTHI',0.D0,' Mimimum Thickness Ratio...........')
   !maxth = getrea('MAXTHI',1.D3,' Maximum Thickness Ratio...........')
   minstr = getrea('MINSTR',0D0,' Minimum Strain for TTT integration.')
   nlayr = getint('NLAYR ',   0,'!NUMBER OF LAYERS .................')
   IF ( nlayr <= 0) CALL runend('SHELL3: VERIFY NUMBER OF LAYERS ')
   sizei =  4+nlayr*2                !size of integer array
   ALLOCATE (sec%iprop(sizei), rprop(npl,nlayr) )
   sec%iprop(4) = 0
   IF( exists('TSHEAR'))sec%iprop(4) = 1

   youn=0d0              !Initializes mean young modulus
   dens=0d0              !Initializes density
   IF(ith == 0)thick=0d0 !Initializes thickness
   ii =  5               !pointer to integer array
   kvar=1                !pointer to layer internal variables
   DO il=1,nlayr         !for each layer

     CALL listen('SHELL3')  !read a line
     matno    = getint('NUMAT ', 1 ,'!Associated Material ..............')
     ! search for associated material
     CALL mat_search(matno,found,mat,imat)        !search base material
     IF(.NOT.found) CALL runend ('SHELL3: NOT FOUND ASSOCIATED MATERIAL ')
     IF( .NOT. ANY( mat%mtype == matsa) ) &  !check if base material is correct
          CALL runend('SHELL3: Invalid material association')
     sec%iprop(ii)   = imat         !store associated material index (in PMATS)
     sec%iprop(ii+nlayr) = kvar     !pointer to layer varin (internal variables)
     IF(ith == 0) THEN   !thicknes must be computed
       thickl = getrea('THICKL',1.D0,'!LAYER THICKNESS ..................')
       thick = thick + thickl          ! sum layer thickness
     ELSE                !thickness is known
       thickl = getrea('THICKL',100.D0/nlayr,' LAYER THICKNESS ( % ).............')
       thickl = thick * thickl *0.01   ! compute layer thickness (sum not checked)
     END IF
     angle = getrea('ANGLE ',0.0D0, ' ANGLE FROM COMP COORDS TO LA MAT 1')
     youn = youn + thickl*mat%prope(1)  !mean young modulus
     dens = dens + thickl*mat%prope(5)  !surface density
     rprop(1,il) = thickl               !store layer thickness
     rprop(2,il) = angle                !store layer angle
     rprop(3,il) = getrea('POSIT',0.0d0,' LAYER POSITION ..................')

     angle = angle*factor       !angle between orthotropic axes and global axes
     sg = SIN(angle)            !sine
     cg = COS(angle)            !cosine
     ! matrix to transform global strains into orthotropic strains
     rprop(4,il) =  cg**2                ! [ cos^2   sin^2   cos * sin ]
     rprop(5,il) =  sg**2                ! [ sin^2   cos^2  -cos * sin ]
     rprop(6,il) =  cg*sg                ! [-2 c s   2 c s   c^2 - s^2 ]
     rprop(7,il) =  2d0*cg*sg
     rprop(8,il) =  cg**2 - sg**2
     !get required internal variables increment
     kvar = kvar + get_nvar(mat,1)  !1: plane stress

     ii = ii+1                           !updates pointer to integer array

   END DO
   nvars = kvar-1                         !number of internal variables
   sec%iprop(1) = nlayr                   !keep number of layers
   sec%iprop(2) = nvars                   !set total number of internal variables
   sec%iprop(3) = 0                       !initializes to elastic lamina

   sizer = 26             ! minimum size of real array
   IF( nvars>0 )THEN      ! if internal variables exist
     sec%iprop(3) = 1                     !set to plastic lamina
     sizer = sizer+nlayr*npl              !update size of real array
   END IF
   sec%secdef = (/ sizei,sizer, 0,0,0 /)  !size of arrays
   ALLOCATE (sec%rprop(sizer) )           !get memory
   sec%rprop(1) = thick                   !thickness
   sec%rprop(2) = youn/thick              !mean young modulus
   sec%rprop(3) = dens/thick              !mean density
   sec%rprop(4) = minstr                  !strain threshold  (TO BE COMPUTED!!!)
   sec%rprop(5) = dens                    !surface density

   ! compute elastic integrated matrix
   de = 0d0                               !initializes integrated thickness
   zauto = ( rprop(3,nlayr) - rprop(3,1) ) == 0.0d0 !layer positions not given
   IF (zauto)  zl = -thick/2d0      !initializes Z position
   DO il=1,nlayr                           !for each layer
     imat= sec%iprop(4+il)                   !material index
     mat => pmats(imat)%p                    !point to material
     CALL elastiff(celoc,mat,1)              !get elastic stiffnes matrix
     thickl = rprop(1,il)                          !layer thickness
     IF (zauto)  THEN
       zl = zl+thickl/2d0                          !update Z position
       rprop(3,il) = zl                            !keep
     ELSE
       zl = rprop(3,il)                            !update Z position
     END IF
     aux= thickl*thickl/12d0                       !auxiliar h^2/12
     CALL rotortm(celoc,rprop(4:,il),ce)   !rotate elastic in-plane matrix
     ce = ce*thickl
     de( 1: 6) = de( 1: 6) + ce                    !membrane part
     de( 7:12) = de( 7:12) + ce*(zl*zl+aux)        !bending part
     de(13:15) = de(13:15) + ce(1:3)*zl            !coupled part
     de(17:18) = de(17:18) + ce(4:5)*zl
     de(   21) = de(   21) + ce(  6)*zl
     zl = zl+thickl/2d0                          !update Z position
   END DO
   de(16) = de(14)                                 !symmetry values
   de(19) = de(15)
   de(20) = de(18)
   sec%rprop(6:26) = de                            !keep integrated matrix
   WRITE(58,"(' Integrated D matrices')")
   WRITE(58,"(5e15.8)")de
   IF( nvars > 0 ) THEN      !for inelastic behavior
     ir = 27            !initializes pointer to layer thickness
     ii = ir + nlayr    !initializes pointer to angle
     j  = ii + nlayr    !initializes pointer to z
     k  = j  + nlayr    !initializes pointer to rotation matrices
     DO il=1,nlayr    !for each layer
       sec%rprop(ir) = rprop(1,il)               !layer thickness
       sec%rprop(ii) = rprop(2,il)               !layer angle
       sec%rprop( j) = rprop(3,il)               !layer position z
       CALL vecasi( 5,rprop(4,il),sec%rprop(k))  !rotation matrices
       ! update pointers
       ir = ir+ 1   !layer thickness
       ii = ii+ 1   !angle
       j  = j + 1   !z
       k  = k + 5   !rotation matrices
     END DO
   END IF
   IF( sec%iprop(4) > 0 ) THEN  !If shear stress is considered
     !- Compute shear factors for each layer --------------------------------------
     DO il = 1, nlayr                         ! Shear stress through the section - XMG
       sec%rprop(k) = - ( 3.0d0 / (2.0d0 * sec%rprop(1)) )                & ! Computed as: 1/2I * (h^2/4 - z^2) * V
                      + ( 6.0d0 * (rprop(3,il)**2) / (sec%rprop(1)**3) )    ! Where I = b * h^3 / 2 (with b=1)
       k = k+1                                ! h = section thickness
     END DO                                   ! z = layer midpoint
   END IF
   DEALLOCATE (rprop)                     !release local memory

 !------------------------------------------------------------------------------
 END IF
 RETURN
 END SUBROUTINE sc_shel3
!********************************************************************
 SUBROUTINE rotortm (c,rr,ce,flag)
 IMPLICIT NONE
 !
 !     rotate the stress-strain matrix to global/reference coordinates
 !
 !     set the coordinate transformation for rotation of properties
 !        | c(1)  c(2)      |
 !   Cm = | c(2)  c(3)      |    Cs = | c(5)      !
 !        |            c(4) |         |      c(6) !
 !  r = { cos^2  sin^2  cos*sen  2cos*sin  con^2-sin^2 }
 REAL (kind=8) :: c(:),rr(:),ce(:)
 LOGICAL, OPTIONAL :: flag
 !local
 REAL (kind=8) d(3,3)

 !     rotate the stress-strain matrix
 !     t(transpose) * c(material)

 IF( rr(1) < 0.9999999999 )THEN

!   WHICH IS THE BEST VERSION
!   c2 = rr(1)          !cos^2(t)
!   c4 = c2*c2          !cos^4(t)
!   s2 = rr(2)          !sin^2(t)
!   s4 = s2*s2          !sin^4(t)
!   c4s4 = c4+s4        !cos^4(t)+sin^4(t)
!   sc = rr(3)          !cos(t)*sin(t)
!
!   ce(1) = c(1)*c4 + (c(2) +2*c(4))*2*s2*c2 + c(3)*s4                 !c(1,1)
!   ce(2) = (c(1)+c(3)-4*c(4))*s2*c2 + c(2)*(s4+c4)                    !c(1,2)
!   ce(3) = (c(1)-c(2)-2*c(4))*sc*c2 + (c(2)-c(3)+2*c(4))*s2*sc        !c(1,3)
!   ce(4) = c(1)*s4 + (c(2) +2*c(4))*2*s2*c2 + c(3)*c4                 !c(2,2)
!   ce(5) = (c(1)-c(2)-2*c(4))*sc*s2 + (c(2)-c(3)+2*c(4))*c2*sc        !c(2,3)
!   ce(6) = (c(1)+c(3)-2*c(2)-2*c(4))*s2*c2 + c(4)*(s4+c4)             !c(3,3)
!   IF( PRESENT(flag) )THEN
!     ce(7) = c(5)*c2+c(6)*s2
!     ce(8) = (c(6)-c(5))*sc
!     ce(9) = c(5)*s2+c(6)*c2
!   END IF
   !d = MATMUL(TRANSPOSE(t),c)
   d(1,1) = rr(1)*c(1) + rr(2)*c(2)
   d(1,2) = rr(1)*c(2) + rr(2)*c(3)
   d(1,3) = rr(4)*c(4)
   d(2,1) = rr(2)*c(1) + rr(1)*c(2)
   d(2,2) = rr(2)*c(2) + rr(1)*c(3)
   d(2,3) =-rr(4)*c(4)
   d(3,1) =-rr(3)*c(1) + rr(3)*c(2)
   d(3,2) =-rr(3)*c(2) + rr(3)*c(3)
   d(3,3) = rr(5)*c(4)

   !     t(transpose) * c(material) * t
   !ce = MATMUL(d,t)
   ce(1) = d(1,1)*rr(1)+ d(1,2)*rr(2)+ d(1,3)*rr(4)
   ce(2) = d(1,1)*rr(2)+ d(1,2)*rr(1)- d(1,3)*rr(4)
   ce(3) =-d(1,1)*rr(3)+ d(1,2)*rr(3)+ d(1,3)*rr(5)
   ce(4) = d(2,1)*rr(2)+ d(2,2)*rr(1)- d(2,3)*rr(4)
   ce(5) =-d(2,1)*rr(3)+ d(2,2)*rr(3)+ d(2,3)*rr(5)
   ce(6) =-d(3,1)*rr(3)+ d(3,2)*rr(3)+ d(3,3)*rr(5)

   IF( PRESENT(flag) )THEN
     ce(7) = c(5)*rr(1) + c(6)*rr(2)
     ce(8) = (c(6)-c(5))*rr(3)
     ce(9) = c(6)*rr(1) + c(5)*rr(2)
   END IF
 ELSE
   ce(1) = c(1)
   ce(2) = c(2)
   ce(3) = 0d0
   ce(4) = c(3)
   ce(5) = 0d0
   ce(6) = c(4)
   IF( PRESENT(flag) )THEN
     ce(7:9) = (/ c(5), 0d0, c(6) /)
   END IF
 END IF
 RETURN
 END SUBROUTINE rotortm
!********************************************************************
