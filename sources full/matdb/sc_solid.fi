 SUBROUTINE sc_solid (sec)
 ! For non-layered solid
 !    iprop: 1 :nlayr
 !           2 :nvars:  number of internal variables for the whole laminate
 !           3 :number of associated FLC curve
 !    rprop  1 :thick
 !            --------------------------------
 ! For multi-layered lamina
 !    iprop: 1 :nlayr
 !           2 :nvars:  number of internal variables for the whole laminate
 !           3 :0  elastic lamina,  1:inelastic
 !      4:3+nl :associated material label to each layer
 !  4+nl:3+2nl :pointer to each component internal variables
 !    rprop  1 :thick
 !           2 :mean young modulus
 !           3 :mean density
 !           4 :elastic generalized strain threshold
 !           5 :surface density
 !        6-41 :integrated elasticity square matrix
 !        For each layer (npt values)
 !            1: layer thickness
 !            2: layer angle
 !            3: layer z coordinate
 !          4-8: rotation coefficients
 IMPLICIT NONE
 TYPE (section), POINTER :: sec

 ! local variables

 LOGICAL :: found
 INTEGER (kind=4) matno
 INTEGER (kind=4) :: matsa(7) = (/ 1,4,5,6,8,12,30 /)
 INTEGER (kind=4), PARAMETER :: npl = 8     !number of real values per layer
 REAL (kind=8), PARAMETER :: factor = 0.017453292519943
 INTEGER (kind=4) nlayr,ith,il,sizer,sizei,ir,ii,j,k,nvars,kvar,imat,mtype
 REAL(kind=8) :: thick,youn,dens,thickl,de(36),angle,ce(6),sg,cg, &
                 minstr,celoc(4),ce36(36)
 REAL(kind=8), ALLOCATABLE :: rprop(:,:)
 TYPE (mater), POINTER :: mat
!-----------------------------------------------------------------------------------------------

 CALL listen('SC_SOL')

 sec%secty = 0                ! solid
 thick = getrea('THICK ',1.d0 ,' Thickness for Plane Stress Problem')

 IF (exists('NUMAT ')) THEN   ! S I N G L E   M A T E R I A L (NOT LAYERED)
   matno = getint('NUMAT ', 1,'!Associated Material ..............')
   CALL mat_search(matno,found,mat,imat)        !search base material
   IF(.NOT.found) CALL runend ('SC_SOL: NOT FOUND ASSOCIATED MATER ')
   IF( .NOT. ANY( mat%mtype == matsa) ) &
      CALL runend('SC_SOL: Invalid material association')

   sec%mtbas => mat   !pointer to material db
   sec%mabas = matno  !material label

   sec%secdef = (/ 3,1,0,0,0 /)
   ALLOCATE (sec%iprop(sec%secdef(1)), sec%rprop(sec%secdef(2)) )

   nvars = get_nvar(mat,0)  !0: 3D solid
   sec%iprop(1) = 0         !number of layers
   sec%iprop(2) = nvars     !set total number of internal variables
   sec%iprop(3) = getint('LBLFLC',   0,' FLC LABEL ........................') ! read FLC label
   sec%rprop(1) = thick     !save thickness

   IF( sec%mtbas%mtype == 6 )THEN
     IF( sec%mtbas%matdef(3) > 1 )THEN
       WRITE(55,"(' WARNING: Hyperelastic rubber can only be elastic for solid elements')")
       WRITE(lures,"(' WARNING: Hyperelastic rubber can only be elastic for solid elements')")
     END IF
   END IF
   !read FLC curve
   CALL rdflc(sec%iprop(3))

 ELSE    ! L A Y E R E D   S O L I D

   IF( thick > 0 ) THEN
     ith=1 ! thickl will be input as a percentage of thick
   ELSE
     ith=0 ! thickl will be input as a REAL DIMENSION
   END IF
   minstr = getrea('MINSTR',0.D0,' Minimum Strain for TTT integration.')
   nlayr = getint('NLAYR ',   0,'!NUMBER OF LAYERS .................') !compulsory
   IF ( nlayr < 1) CALL runend('SC_SOL: NUMBER OF LAYERS MUST BE POSITIVE')
   sizei =  4+nlayr*2
   ALLOCATE (sec%iprop(sizei), rprop(npl,nlayr) )
   sec%iprop(4) = 0  !Available memory space for an integer value
                     !In shell sections is used as TSHEAR flag
                     ! IF( exists('TSHEAR'))sec%iprop(4) = 1

   youn=0d0              !Initializes mean young modulus
   dens=0d0              !Initializes density
   IF(ith == 0)thick=0.0 !Initializes thickness
   ii =  5               !pointer to integer array
   kvar=1
   DO il=1,nlayr         !for each layer

     CALL listen('SC_SOL')  !read a line
     matno    = getint('NUMAT ', 1 ,'!Associated Material ..............')
     ! search for associated material
     CALL mat_search(matno,found,mat,imat)        !search base material
     IF(.NOT.found) CALL runend ('SC_SOL: NOT FOUND ASSOCIATED MATERIAL ')
     IF( .NOT. ANY( mat%mtype == matsa) ) &  !check if base material is correct
          CALL runend('SC_SOL: Invalid material association')
     sec%iprop(ii)   = imat         !store associated material index (in PMATS)
     sec%iprop(ii+nlayr) = kvar      !pointer to layer varin (internal variables)
     IF(ith == 0) THEN
       thickl = getrea('THICKL',1.D0,'!LAYER THICKNESS ..................')
       thick = thick + thickl          ! sum layer thickness
     ELSE
       thickl = getrea('THICKL',100.D0/nlayr,' LAYER THICKNESS ( % ).............')
       thickl = thick * thickl *0.01   ! compute layer thickness (sum not checked)
     END IF
     angle = getrea('ANGLE ',0.0D0, ' ANGLE FROM COMP COORDS TO LA MAT 1')
     youn = youn + thickl*mat%prope(1)  !mean young modulus
     dens = dens + thickl*mat%prope(5)  !surface density
     rprop(1,il) = thickl               !store layer thickness
     rprop(2,il) = angle                !store layer angle
     rprop(3,il) = 0.0d0
     ! rprop(3,il) = getrea('POSIT',0.0d0,' LAYER POSITION ..................') !irrelevant for solid sections
     ! NOTE: LAYER POSITION is not taken into account in solid elems. All layers
     !       are considered to be coexistent at each integration point.

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
     kvar = kvar + get_nvar(mat,0)  !0: 3D solid

     ii = ii+1

   END DO
   nvars = kvar-1
   sec%iprop(1) = nlayr                   !keep number of layers
   sec%iprop(2) = nvars                   !set total number of internal variables
   sec%iprop(3) = 0                       !elastic lamina

   sec%mabas = -1 !flag to indicate no single material is associated

!   sizer = 26
   sizer = 41
   IF( nvars>0 )THEN
     sec%iprop(3) = 1                     !plastic lamina
     sizer = sizer+nlayr*npl              !
   END IF
   sec%secdef = (/ sizei,sizer, 0,0,0 /)   !size of arrays
   ALLOCATE (sec%rprop(sizer) )           !get memory
   sec%rprop(1) = thick                   !thickness
   sec%rprop(2) = youn/thick              !mean young modulus
   sec%rprop(3) = dens/thick              !mean density
   sec%rprop(4) = minstr                  !strain threshold           (TO BE COMPUTED!!!)
   sec%rprop(5) = dens                    !surface density

   ! compute elastic integrated matrix
   de = 0d0                        !initializes integrated matrix
   DO il=1,nlayr                   !for each layer
     imat= sec%iprop(4+il)
     mat => pmats(imat)%p
     CALL elastiff(ce36,mat,0)     !get elastic stiffnes matrix  ! 0: 3D solid
     celoc(1)= ce36(1) !c11
     celoc(2)= ce36(2) !c12
     celoc(3)= ce36(8) !c22
     celoc(4)= ce36(15)!c44
     CALL rotortm(celoc,rprop(4:,il),ce)  !rotate elastic in-plane matrix
     ce36( 1)= ce(1)   !c11
     ce36( 2)= ce(2)   !c12
     ce36( 3)= ce(3)   !c14
     ce36( 8)= ce(4)   !c22
     ce36( 9)= ce(5)   !c24
     ce36(15)= ce(6)   !c44
     !symmetries
     ce36( 7)= ce(2)   !c21
     ce36(13)= ce(3)   !c41
     ce36(14)= ce(5)   !c42
     thickl = rprop(1,il)/thick                 !layer thickness fraction
     de(1:36) = de(1:36) + ce36(1:36)*thickl    !global stiffness contribution (rule of mixtures)
   END DO
   sec%rprop(6:41) = de(1:36)                   !keep integrated matrix
   IF( nvars>0 ) THEN
     ir = 42            !initializes pointer to layer thickness
     ii = ir + nlayr    !initializes pointer to angle
     j  = ii + nlayr    !initializes pointer to z
     k  = j  + nlayr    !initializes pointer to rotation matrices
     DO il=1,nlayr    !for each layer
       sec%rprop(ir) = rprop(1,il)               !layer thickness
       sec%rprop(ii) = rprop(2,il)               !layer angle
       sec%rprop( j) = rprop(3,il)               !layer position z
       CALL vecasi( 5,rprop(4,il),sec%rprop(k))  !rotation matrices
       ir = ir+ 1   !layer thickness
       ii = ii+ 1   !angle
       j  = j + 1   !z
       k  = k + 5   !rotation matrices
     END DO
   END IF

   DEALLOCATE (rprop)                     !release local memory

 END IF

 IF (sec%mabas==-1) THEN
   mtype = -1
 ELSE
   mtype = sec%mtbas%mtype
 END IF


RETURN
END SUBROUTINE sc_solid
