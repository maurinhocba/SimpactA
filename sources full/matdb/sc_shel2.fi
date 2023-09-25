 SUBROUTINE sc_shel2( sec )
 !            --------------------------------
 !            input material properties for shell constitutive model
 !            Total Lagrangean BST elements (NBST,LBST,BSQ)
 !            --------------------------------
 !  iprop: 1 :nlayr    : number of integration points across thickness
 !         2 :4  number of internal variables per layer
 !         3 :number of associated FLC curve
 !         4 :1 if transverse shear is to be considered in constitutive equation
 !
 !  rprop  1 :thick
 !         2 :minimum thickness ratio
 !         3 :maximum thickness ratio
 !         4 :
 !        5:8: integrated elasticity matrix (membrane)
 !       9:12: integrated elasticity matrix (bending)
 IMPLICIT NONE
 TYPE (section), POINTER :: sec
 TYPE(postv), POINTER :: postp,post1     !pointer to a post-process variable

 !local variables
 INTEGER (kind=4) :: matsa(6) = (/ 1,4,5,6,8,30 /) ! material types allowed for this section
 INTEGER (kind=4) matno,nvar,nvarn,dim,dim2,nnam,mtype,nlayr,l,k
 LOGICAL :: found
 REAL(kind=8) :: d
 REAL(kind=8), ALLOCATABLE :: thf(:)

 sec%secty = 12         !shell type 2

 CALL listen('SHELL2')  !read a line

 matno    = getint('NUMAT ', 1 ,'!Associated Material ..............')
 ! search for associated material
 CALL mat_search(matno,found,sec%mtbas)  !search base material
 sec%mabas = matno                       !material label
 IF(.NOT.found) CALL runend ('SHELL2: NOT FOUND ASSOCIATED MATERIAL ')
 mtype = sec%mtbas%mtype
 IF( .NOT. ANY( mtype == matsa) ) & !check if base material is correct
      CALL runend('SHELL2: Invalid material association')

 nlayr = getint('NLAYR ',   3,' NUMBER OF LAYERS .................')
 sec%secdef = (/ 4,12+nlayr,0,0,0 /)   !size of arrays
 ALLOCATE (sec%iprop(sec%secdef(1)), sec%rprop(sec%secdef(2)) )
 sec%iprop = 0
 sec%rprop = 0d0
 !  read parameters
 sec%iprop(1) = nlayr
 sec%rprop(1) = getrea('THICK ',1.D0,'!Thickness ........................')
 IF (sec%rprop(1) == 0.D0) CALL runend('SHELL2: VERIFY THICKNESS ')
 sec%rprop(2) = getrea('MINTHI',0.D0,' Mimimum Thickness Ratio...........')
 sec%rprop(3) = getrea('MAXTHI',1.D3,' Maximum Thickness Ratio...........')
 sec%rprop(4) = getrea('MINSTR',0d0,' Minimum Strain for TTT integration.')
 IF (sec%iprop(1)  == 0) CALL runend('SHELL2: VERIFY NUMBER OF LAYERS ')

 sec%iprop(3) = getint('LBLFLC',   0,' FLC LABEL ........................') ! read FLC label
 sec%iprop(4) = 0
 IF( exists('TSHEAR'))sec%iprop(4) = 1

 IF( sec%mtbas%matdef(2)==1 )RETURN  ! RIGID MATERIAL

 IF( mtype == 30 )THEN
   CALL listen('SHELL2')  !read a line
   IF( exists ('POSTPR')) THEN
     sec%secdef(5) = getint('NVAR  ',   0,' NUMBER OF VARIABLES (scalars).....')
     nvar = 0
     nvarn = 0
     ALLOCATE(sec%postp)
     postp => sec%postp
     DO
       CALL listen('SHELL2')  !read a line
       IF( exists('ENDPOS') )EXIT
       nvar = nvar+1
       nnam = -getint('NAME  ', 1,'!-')
       postp%name = '              '
       IF( nnam == 0 ) CALL runend('SHELL2: Variable names missing')
       SELECT CASE (words(1)(1:midn))
       CASE ('SCALAR')
         postp%type = 0
         postp%dim  = 0
         postp%name(1) = names(1) !first 15 characters ?
         postp%name(2:7) = '               '
         nvarn = nvarn + 1
         WRITE(lures,"(9X,'Variable Name ',24('.'),' NAME   = ',A)",ERR=9999) TRIM(postp%name(1))
       CASE ('VECTOR')
         postp%type = 1
         dim = getint('DIM   ', 2,'!Vector dimension .................')
         postp%dim = dim
         postp%name(1:nnam) = names(1:nnam) !first 15 characters ?
         WRITE(lures,"(9X,'Variable Name ',24('.'),' NAME   = ',A)",ERR=9999) TRIM(postp%name(1))
         IF( nnam > 1 ) WRITE(lures,"(57X,A)",ERR=9999) postp%name(2:nnam)
         nvarn = nvarn + dim
       CASE ('TENSOR')
         postp%type = 2
         dim = getint('DIM   ', 2,'!Tensor dimension .................')
         dim2 = dim*(dim+1)/2
         postp%dim = dim
         postp%name(1:nnam) = names(1:nnam) !first 15 characters ?
         WRITE(lures,"(9X,'Variable Name ',24('.'),' NAME   = ',A)",ERR=9999) TRIM(postp%name(1))
         IF( nnam > 1 ) WRITE(lures,"(57X,A)",ERR=9999) postp%name(2:nnam)
         nvarn = nvarn + dim2
       CASE DEFAULT
         CALL runend('Error reading post-process variables')
       END SELECT
       ALLOCATE (post1)                  !get memory for next variable
       postp%next => post1               !put in the list
       postp => postp%next               !go to next
     END DO
     sec%secdef(4) = nvar
     IF( (sec%secdef(5) /= nvarn) .AND. (sec%secdef(5) /= 0))THEN
       CALL runend('SHELL2 , No of vars to post-process is incorrect')
     END IF
     sec%secdef(5) = nvarn
   ELSE
     backs = .TRUE.
   END IF
!---------------------------------------------------------------
 ELSE
   nvar = 0
   IF( sec%mtbas%matdef(3) > 1 ) nvar = 4
   ! kinematic hardening not implemented yet
   IF( sec%mtbas%matdef(5) > 1 )THEN
   !  WRITE(lures,"(' **************** WARNING **************************',//, &
   !              & ' kinematic hardening ignored for section type SHELL2')",ERR=9999)
   !  WRITE(55   ,"(' **************** WARNING **************************',//, &
   !              & ' kinematic hardening ignored for section type SHELL2')",ERR=9999)
     nvar = nvar + 3
   END IF
   sec%iprop(2) = nvar  !number of internal variables
   !
   !compute elastic matrix
   d = sec%rprop(1)       !thickness
   ! membrane matrix
   SELECT CASE (mtype )
   CASE ( 1 , 8)        !isotropic
     IF( sec%mtbas%matdef(2) > 1 ) sec%rprop(5:8) = sec%mtbas%prope(7:10)*d
     IF(sec%rprop(4) == 0d0) sec%rprop(4) = 0.5d-3  !minimum strain to consider full integration
   CASE ( 4 )           !thermal isotropic
     IF( sec%mtbas%matdef(2) > 1 ) sec%rprop(5:8) = sec%mtbas%prope(10:13)*d
   CASE ( 5 )           !orthotropic
     sec%rprop(5:8) = sec%mtbas%prope(16:19)*d
     IF(sec%rprop(4) == 0d0) sec%rprop(4) = 0.5d-3  !minimum strain to consider full integration
   CASE ( 6 )           !isotropic rubber
     sec%rprop(5:8) = sec%mtbas%prope(19:22)*d
     IF(sec%rprop(4) == 0d0) sec%rprop(4) = 0.5d-3  !minimum strain to consider full integration
   END SELECT
   !bending matrix
   d = d**2/12d0
   sec%rprop(9:12) = sec%rprop(5:8)*d
 END IF
 ! read FLC curve
 CALL rdflc(sec%iprop(3))
 ! compute shear distribution
 ALLOCATE( thf(nlayr) )
 k = 13
 CALL gaussq(nlayr,thf(1),sec%rprop(k)) !integration points through the thickness
 thf(1:nlayr) = thf(1:nlayr)/2d0  !positions
 DO l=1,nlayr                       ! Shear stress through the section - XMG
   sec%rprop(k) = (1.5d0/sec%rprop(1))*(4d0*thf(l)**2 - 1d0)  ! Max shear stress: 3/2 * 1/(t*b) * V
   k = k+1                          ! Where b = 1.0 (V is computed by unit length)
 END DO                             ! and   t = shell thickness

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE sc_shel2
