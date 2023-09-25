 SUBROUTINE csdat2(maxve,iprint,coord,surf,check,cnode)

 !.... read contact surface data and generate database

 IMPLICIT NONE
 !     arguments
 INTEGER (kind=4), INTENT(IN) :: iprint
 INTEGER (kind=4), INTENT(IN) :: maxve
 INTEGER (kind=4), INTENT(IN) :: cnode
 LOGICAL, INTENT(IN) :: check
 REAL (kind=8), POINTER :: coord(:,:)
 TYPE (surf2_db), POINTER :: surf  !INTENT(OUT)
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

 IF( surf%iscod .OR. surf%curved )THEN  !for slave surface keep nodes
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
!   CALL mastdb (nsegm,nnseg,surf%lcseg,surf%nhseg,.TRUE.,coord)
   CALL mastdb (nsegm,nnseg,surf%lcseg,surf%nhseg,surf%confor,coord)
   ALLOCATE( surf%xc(2,nsegm) )
   surf%cxc = .FALSE.
   IF(check)THEN
     CALL ccheck2(nsegm,surf%lcseg,surf%nhseg,cnode,coord,surf%xc)
     surf%cxc = .TRUE.
   END IF
 END IF

 surf%ncnod = ncnod     !keep number of nodes
 surf%nsegm = nsegm     !keep number of segments

 RETURN
 END SUBROUTINE csdat2
