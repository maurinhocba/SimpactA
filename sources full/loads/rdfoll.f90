 SUBROUTINE rdfoll (numfl,ndime,iwrit,fltype,flparm,headf,tailf,hydfl,curv_name, &
                    fluid)
 !Reads follower loads
 USE c_input,ONLY: ludat, lures, listen, exists, get_name, getint, getrea, backs,  param
 USE loa_db,ONLY: foll_seg, flpar, hydf_load, ini_hydfl, rdhyd, new_foll, add_foll
 USE npo_db,ONLY: coora
 USE ctrl_db,ONLY: ntype
 USE surf_db
 IMPLICIT NONE

   !--- Dummy variables
   CHARACTER(len=*):: fltype      !len=6
   CHARACTER(len=*):: curv_name   !
   INTEGER(kind=4):: numfl, ndime, iwrit
   TYPE(flpar):: flparm
   TYPE(foll_seg),POINTER:: headf, tailf
   TYPE(hydf_load),POINTER:: hydfl
   LOGICAL :: fluid
   !--- Functions
   INTEGER(kind=4):: chnode
   !--- Local variables
   REAL (kind=8), PARAMETER  :: pi=3.141592653589793
   INTEGER(kind=4):: nn, nnode, i, j, nelem
   REAL(kind=8):: p(4),vol,x0,y0,xl,xr,rr
   TYPE(foll_seg),POINTER:: seg
   LOGICAL :: water,full,elset,found,reverse,axisy
   CHARACTER(len=mnam) :: sname
   TYPE(srf_seg),POINTER:: segs

   INTERFACE
     INCLUDE 'elemnt.h'
     INCLUDE 'rdhyd2d.h'
   END INTERFACE

   CALL listen('RDFOLL')           !reads a line

   numfl = 0                       !initializes number of segments

   IF( .NOT.exists('FOLLOW')) THEN !check key-word
     backs = .TRUE.                !one line back
     NULLIFY(hydfl)                !nullifies associated pointer
   ELSE
     IF (iwrit == 1) WRITE(lures,"(/' ====== FOLLOWER LOAD ======'//)",ERR=9999)

     fltype = 'SIMPLE'          !initializes follower type
     flparm%ltcur = curv_name   !initializes name
     fluid = EXISTS('FLUID')    !pressure is provided by an external program

     ! for inflated structs or tube hydroforming
     CALL listen('RDFOLL')

     IF (exists('INFLAT') .OR. exists('TUBHYD')) THEN
       flparm%nodcen =  getint('NODCEN',0,'!Number of cen. point for vol. calc')
       flparm%timini =  getrea('TIMINI',0d0,'!Start time of flux/temp control...')
       flparm%nodcen = chnode(flparm%nodcen)
       flparm%ltcur  = curv_name
       IF (exists('INFLAT')) THEN
         fltype = 'INFLAT'           !inflatable surface
         flparm%pext = getrea('PEXT  ',0d0,'!Exterior pressure                 ')
       ELSE IF ( exists('TUBHYD') ) THEN
         fltype = 'TUBHYD'           !tube hydroforming
         flparm%kvol = getrea('KVOL  ',0d0,'!Bulk modulus of fluid ........... ')
         flparm%qref = getrea('FLUXRE',0d0,'!Reference flux .................. ')
       END IF
     ELSE IF ( exists('HYDROF') ) THEN
       IF (ndime == 2) THEN
         fltype = 'SIMPLE'   !Not enclosed line considered ==> defined as follower load type
         CALL rdhyd2d(numfl,headf,tailf)
       ELSE
         fltype = 'HYDROF'
         CALL ini_hydfl(hydfl,ndime)
         CALL rdhyd(hydfl)
       END IF
       CALL listen('RDFOLL')              !read Final key-word line
       IF (.not.exists('ENDFOL')) CALL runend('RDFOLL:END_FOLLOWER card expected      ') !stop if not found
       RETURN
     ELSE IF ( exists('WATERB') ) THEN
       fltype = 'WATERB'           !water bags problems
       IF( ntype == 3 )THEN
         flparm%nodcen = 3     !flag
         flparm%timini =  0d0
         axisy = .TRUE.
       ELSE
         axisy = .FALSE.
         IF( exists ('SYMMET') )THEN
           flparm%nodcen = 1     !flag
           flparm%timini =  getrea('XSYMME',0d0,' X coordinate for symmetry.........')
         ELSE
           flparm%timini =  getrea('XREFER',0d0,' X coordinate for reference........')
           flparm%nodcen = 2
         END IF
       END IF
       flparm%qref = getrea('SPECIF',0d0,'!Fluid Specific Weight............ ')
       flparm%kvol = getrea('KVOL  ',0d0,'!Bulk modulus of fluid ........... ')
       IF( exists('FULL') )THEN
         flparm%nodcen = flparm%nodcen + 10
         flparm%press= getrea('YINITI',0d0,' Initial Y Coordinate ............ ')
       ELSE IF( exists('FILLIN') ) THEN
         flparm%nodcen = flparm%nodcen + 20
         flparm%press= getrea('YINITI',0d0,' Initial Y Coordinate ............ ')
       ELSE
         flparm%press= getrea('YINITI',0d0,'!Initial Y Coordinate ............ ')
       END IF
     ELSE
       backs = .TRUE.
     END IF


     IF( TRIM(fltype) == 'WATERB' )THEN
       water = .TRUE.                  !flag
       y0    = flparm%press            !water surface coordinate
       x0    = flparm%timini           !x-coordinate for volume computation
       vol = 0d0             !Initializes volume
       rr = 1d0              !Initializes radius
       IF( flparm%nodcen > 10 ) THEN
         xr=-1e9
         full = .TRUE.
       ELSE
         xl = 0d0
         xr = 0d0  !initializes coordinates for free surfaces
         full = .FALSE.
       END IF
     ELSE
       water = .FALSE.
     END IF

     ! see if surface is defined by a set
     !*** loop over each loaded face
     NULLIFY(headf,tailf)     !nullify data base
     numfl = 0                !initializes number of segments in surface
     DO
       CALL listen('RDFOLL')              !read a line
       IF (exists('ENDFOL')) EXIT !exit loop when Final key-word read
       IF( exists('ELSET',i)) THEN        !check key-word
         sname = get_name(posin=i,stype='SURF')     !read surface name
         CALL new_srf(surfa)                        !initializes surface
         CALL elemnt ('SURFAC',name=sname,flag2=found)   !read surface definition
         IF (.NOT.found) CALL runend('RDFOLL:ELEMENT SET NOT FOUND       ') !stop if not found
         WRITE(lures,"('  Follower load applied over Element set ',a)")TRIM(sname)
         nelem = surfa%nelem              !number of segments in surface
         numfl = numfl + nelem            !number of segments in surface
         p(1) = getrea('PRESS',1d0 ,' Pressure value at element ')  !average pressure
         reverse = exists('REVERS')      !if order must be changed
         IF( reverse ) WRITE(lures,"(' set connectivities are REVERSED to compute follower loads')")
         i = 0                            !initializes index
         segs => surfa%head               !point to firt segment in the list
         ! set number of nodes per segment when connectivities are read from a set
         IF (ndime == 2) THEN
           nnode = 2           !line segment
         ELSE
           nnode = 4           !quadrilateral
           IF( surfa%head%nodes(3) == surfa%head%nodes(4) )nnode = 3  !check if first segment is a triangle
         END IF
         elset = .TRUE.
       ELSE
         elset = .FALSE.
         backs = .TRUE.
       END IF
       DO                       !loop over each segment
         CALL new_foll(seg)       !get memory
         IF( elset )THEN             !get connectivities from set
           i = i+1                      !increase counter
           IF( i > nelem )THEN             !all elements read
             CALL dalloc_srf(surfa)          !release memory
             EXIT                            !and exit loop
           ELSE
             IF( reverse )THEN             !if connectivites must b reversed
               DO j=1,nnode
                 seg%lnofl(j) = segs%nodes(nnode+1-j)
               END DO
             ELSE
               seg%lnofl(1:nnode) = segs%nodes(1:nnode)
             END IF
             segs => segs%next  !point to next segment in conns
           END IF
         ELSE                !read connectivities from data input files
           CALL listen('RDFOLL')           !read a line
           IF (exists('ELSET ') .OR. exists('ENDFOL')) THEN !exit loop when Final key-word read or new element set
             backs = .TRUE.                !back one line in data input file
             EXIT                          !exit loop for connectivities
           END IF

           numfl = numfl + 1    ! increase number of segments
           nnode = ndime        ! line segment in 2D and Triangle in 3D (overrided using NOD1..NOD4)
           IF( exists('NOD1',i) )THEN
             nn = INT(param(i)) ! getint('NOD1  ',0,'!Node 1 of element surface ........')
             seg%lnofl(1) = chnode(nn)
             found = exists('NOD2',i)
             nn = INT(param(i)) ! getint('NOD2  ',0,'!Node 2 of element surface ........')
             seg%lnofl(2) = chnode(nn)
             IF (ndime == 3) THEN
               found = exists('NOD3',i)
               nn = INT(param(i)) ! getint('NOD3  ',0,'!Node 3 of element surface ........')
               seg%lnofl(3) = chnode(nn)
               found = exists('NOD4',i)
               IF (found)THEN
                 nn = INT(param(i)) ! getint('NOD4  ',0,'!Node 4 of element surface ........')
                 seg%lnofl(4) = chnode(nn)
                 nnode = 4
               END IF
             END IF
           ELSE
             DO i=1,nnode
               nn = INT(param(i))
               seg%lnofl(i) = chnode(nn)
             END DO
           END IF

         END IF

         IF( .NOT.water )THEN
           IF( elset )THEN  ! for element sets
             seg%fload = p(1) !average nodal pressure
           ELSE
             IF( exists('PRESS1',i))THEN
               p(1) = param(i)  !getrea('PRESS1',0d0 ,'!Pressure value at element node 1 ')
               IF( exists('PRESS2',i ))THEN
                 p(2) = param(i)  !getrea('PRESS2',P(1),' Pressure value at element node 2 ')
                 IF (nnode > 2) THEN
                   IF( exists('PRESS3',i )) &
                   p(3)= param(i) !getrea('PRESS3',p(1),' Pressure value at element node 3 ')
                   IF (nnode > 3) THEN
                     IF( exists('PRESS4',i )) &
                     p(4)= param(i) !getrea('PRESS4',p(1),' Pressure value at element node 4 ')
                   END IF
                 END IF
               ELSE
                 p(2:nnode) = p(1)
               END IF
             ELSE
               p(1) = param(nnode+1)
               IF( p(1) == 0d0 )p(1) = 1d0 !default value
               p(2:nnode) = p(1)
             END IF

             seg%fload = SUM(p(1:nnode))/nnode  !average nodal pressure

           END IF
           IF( nnode == 4 .AND. fluid )THEN    !loads provided by other program
             nnode = 3                        !use 3-node segments only
             seg%nnode = nnode                !set nnode of the segment
             CALL add_foll(seg, headf, tailf) !add to list with 3 nodes
             ! New segment
             numfl = numfl + 1                !create a new segment
             CALL new_foll(seg)
             seg%lnofl(1) = tailf%lnofl(1)    !connectivities
             seg%lnofl(2) = tailf%lnofl(3)
             seg%lnofl(3) = tailf%lnofl(4)
             seg%fload = tailf%fload          !average pressure
           END IF
         END IF

         seg%nnode = nnode
         IF( iwrit == 1 ) WRITE(lures,"('nodes', 3i8,'average pressure: ',e15.4 )")seg%lnofl(1:nnode),seg%fload
         CALL add_foll(seg, headf, tailf)
       END DO   !loop over SEGMENTS
     END DO  !loop to read ELSET or list of SEGMENTS

   END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE rdfoll
