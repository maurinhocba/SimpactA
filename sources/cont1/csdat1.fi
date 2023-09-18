 SUBROUTINE csdat1(maxve,iprint,surf)
 !
 !.... read contact surface data and generate database

 IMPLICIT NONE
 !     arguments
 INTEGER (kind=4), INTENT(IN) :: iprint
 INTEGER (kind=4), INTENT(IN) :: maxve
 TYPE (surf1_db), POINTER :: surf  !INTENT(OUT)
 !     local variables
 INTEGER (kind=4), PARAMETER ::  nnseg = 2
 INTEGER (kind=4) ncnod,nsegm,ncnxx,nsexx,maxcc
 INTEGER (kind=4), ALLOCATABLE :: lcnod(:),lcseg(:,:)

 INTERFACE
   INCLUDE 'csdata.h'
 END INTERFACE

 !.... set control surface DATA

 maxcc  = maxve               !maximum work space available

 ncnod = surf%ncnod           !flag to read nodes
 nsegm = surf%nsegm           !flag to read segments

 !.... READ node numbers and segment connectivities

 ! set maximum values for nodes(ncnxx) and segments (nsexx) to be read
 IF( nsegm > 0 )THEN        !if segments are to be read
   ncnxx = MAX( ncnod, maxcc/2, nsegm/2+100 )
   nsexx = MAX( ncnxx,nsegm)
 ELSE                       !only nodes will be read
   ncnxx = MAX( ncnod, maxcc )
   nsexx = 1
 END IF
 ALLOCATE( lcnod(ncnxx), lcseg(nnseg,nsexx) )  !temporary space

 ! read into lcnod <== nodes   lcseg <== segments
 CALL csdata(ncnod,nsegm,ncnxx,nsexx,nnseg,lcnod,lcseg,iprint, &
             coord,surf%sname)

 !.... set memory pointers

! IF( surf%iscod .OR. surf%curved )THEN  !for slave surface keep nodes
 IF( surf%iscod )THEN  !for slave surface keep nodes
   ALLOCATE (surf%lcnod(ncnod))  !reserve space to list of nodes
   surf%lcnod = lcnod(1:ncnod)   !assign list of nodes
 END IF

 IF( nsegm > 0 )THEN          !if segments read keep them
   ALLOCATE ( surf%lcseg(nnseg,nsegm) ) !reserve space for connectivities
   surf%lcseg = lcseg(1:nnseg,1:nsegm)  !assign connectivities
 END IF
 DEALLOCATE( lcnod, lcseg )

 !.... search for neighbours segments if the surface will act as a
 !.... master surface for some pair ==>  NHSEG

! IF( surf%imcod ) THEN
!   !ALLOCATE( surf%xc(2,nsegm) )
!   !surf%cxc = .FALSE.
! END IF

 surf%ncnod = ncnod     !keep number of nodes
 surf%nsegm = nsegm     !keep number of segments

 RETURN
 END SUBROUTINE csdat1
