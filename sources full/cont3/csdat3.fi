 SUBROUTINE csdat3(maxve,iprint,coord,surf)

 !.... read contact surface data and generate database

 IMPLICIT NONE
 !     arguments
 INTEGER (kind=4), INTENT(IN) :: iprint
 INTEGER (kind=4), INTENT(IN) :: maxve
 REAL (kind=8), POINTER :: coord(:,:)
 TYPE (surf3_db), POINTER :: surf  !INTENT(OUT)
 !     local variables
 INTEGER (kind=4), PARAMETER ::  nnseg = 3
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
   ncnxx = MAX( ncnod, maxcc/7, nsegm/2+1000 )
   nsexx = MAX(ncnxx*2,nsegm)
 ELSE                       !only nodes will be read
   ncnxx = MAX( ncnod, maxcc )
   nsexx = 1
 END IF
 ALLOCATE( lcnod(ncnxx), lcseg(nnseg,nsexx) )  !temporary space

 ! read into lcnod <== nodes   lcseg <== segments
 CALL csdata(ncnod,nsegm,ncnxx,nsexx,nnseg,lcnod,lcseg,iprint, &
             coord,surf%sname)

 !.... set memory pointers

 IF(surf%iscod  .OR. surf%curved )THEN  !for slave surface keep nodes
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

 IF( surf%imcod ) THEN
   ALLOCATE( surf%nhseg(nnseg,nsegm) )
   !   form the master surface connection
   WRITE(55,"('     testing connectivities and aspect ratio for surf.',A)",ERR=9999) TRIM(surf%sname)
   CALL mastd3 (nsegm,nnseg,surf%lcseg,surf%nhseg,surf%confor,coord)
   CALL chckar(nsegm,nnseg,surf%lcseg,coord,label)
   ALLOCATE( surf%xc(3,nsegm) )
   surf%cxc = .FALSE.
 END IF

 surf%ncnod = ncnod     !keep number of nodes
 surf%nsegm = nsegm     !keep number of segments

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE csdat3
