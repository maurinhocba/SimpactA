 SUBROUTINE inpda7 (task,nelem, eule0,euler,coord, iwrit,elsnam,nelms)
 !******************************************************************
 !
 !*** READ control DATA for 6-node triangular shell elements
 !    TQQL and TQLL
 !
 !******************************************************************

 USE ctrl_db, ONLY : ndofn
 USE ele07_db
 USE gvar_db,ONLY : maxsdv,fimpo,overw

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam
 CHARACTER(len=*),INTENT(IN):: task
 INTEGER (kind=4) :: nelem,nelms,iwrit
 REAL    (kind=8) :: coord(:,:),eule0(:,:),euler(:,:)

 ! local variables
 LOGICAL :: cmpea, &   !CoMPute Euler Angles
            oldset     !
 INTEGER (kind=4) :: nreqs, narch,i,stype,g,ngaus,nnode,nstre
 INTEGER (kind=4), ALLOCATABLE  :: lnods(:,:)
 TYPE (ele07_set), POINTER :: elset,anter
 TYPE (ele07), POINTER :: el
 REAL (kind=8) :: a,b,deriv(6,2)
 CHARACTER(len=mnam) :: sname     ! element set name

 INTERFACE
   INCLUDE 'nodnor.h'
   INCLUDE 'nodn07.h'
 END INTERFACE

 sname = elsnam

 !  ******************************
 !  assumed membrane strain points

 a = (1d0-1d0/SQRT(3d0))/2d0
 b = (1d0+1d0/SQRT(3d0))/2d0
 CALL shape7(a  ,0d0,deriv(1,1),nfdas(1,1,1))
 CALL shape7(b  ,0d0,deriv(1,1),nfdas(1,1,2))
 CALL shape7(b  ,a  ,deriv(1,1),nfdas(1,1,3))
 CALL shape7(a  ,b  ,deriv(1,1),nfdas(1,1,4))
 CALL shape7(0d0,b  ,deriv(1,1),nfdas(1,1,5))
 CALL shape7(0d0,a  ,deriv(1,1),nfdas(1,1,6))


 IF (TRIM(task) == 'INPUT') THEN
   ! check if list of sets and set exists and initializes
   CALL srch_ele07 (head, anter, elset, elsnam, oldset)
   IF (oldset) THEN    !if set exists
     CALL commv7 (1, nelem,  nreqs, narch, stype, ngaus, elsnam, nnode, nstre, elset)
     cmpea = .FALSE.
   ELSE                !set ELSNAM does not exist
     CALL new_ele07(elset)       !reserve memory for set
     nelem  =  0           !new set, initializes number of elements
     CALL listen('INPDA7')
     nreqs = getint('NREQS ',0,' GAUSS PT FOR STRESS TIME HISTORY .')
     ngaus = getint('NGAUS ',3,' Number of integration points .....')
     stype = getint('STYPE ',0,' Formulation type for this set ....')
     IF( ngaus == 3 )THEN
       elset%ngamm = 6
       nnode = 6
     ELSE
       IF( stype > 1) stype = 1 !TLLL element
       IF( exists('QUAD  ') )THEN
         WRITE(lures,"(/,5X,' Quadratic approach for membrane ', &
                        &      'behavior Selected'//)",ERR=9999)
         nnode = 9
       ELSE
         nnode = 6
       END IF
       elset%ngamm = 3
       elset%stabq=getrea('STABQ ',1d0,' Stabilization factor for shear    ')
     END IF
     cmpea = exists('EULER ')
     IF(cmpea)WRITE(lures,"(/,5X,' Automatic Evaluation of Local ', &
                          &      'Nodal Systems Selected'//)",ERR=9999)
     elset%angdf=getrea('ANGLE ',0d0,' Default angle between X1 and Ort_1')
     IF( exists('LOCAL  ',i)) THEN   !local axis definition
       elset%locax =  INT(param(i))
       IF( elset%locax < 1 .OR. elset%locax > 3 )THEN
         WRITE(lures,"(/,5X,'Error in the definition of local system', / &
                         5X,'Invalid axis: ',i3,' Default value used (3)')") elset%locax
         elset%locax = 3
       END IF
     END IF
     ALLOCATE(elset%ap1(2,elset%ngamm,ngaus))

     IF( ngaus == 3 )THEN
!       gauss points in local coordinates and weigths
       ALLOCATE( elset%posgp(2,ngaus), elset%weigp(ngaus), elset%shape(nnode,ngaus))

       ! assumed membrane strain points

       a = 1d0/3d0
       CALL shape7(a  ,a  ,deriv(1,1),dn(1,1,1))
       CALL shape7(0d0,0d0,deriv(1,1),dn(1,1,2))
       CALL shape7(1d0,0d0,deriv(1,1),dn(1,1,3))
       CALL shape7(0d0,1d0,deriv(1,1),dn(1,1,4))

       ! ******************************
       ! gauss points in local coordinates and weigths
       a = 1d0/6d0
       elset%weigp(:) = a
       ! version for Integration points at mid-side nodes
       b = 1d0/2d0                   !coordinate at mid point
       elset%posgp(:,1) = (/ b  ,0d0 /)    !
       elset%posgp(:,2) = (/ b  ,b   /)    !
       elset%posgp(:,3) = (/ 0d0,b   /)    !
       !! version for Integration points at interior points
       !b = 2d0/3d0                   !coordinate at interior point
       !elset%posgp(:,1)= (/ a,a /)
       !elset%posgp(:,2)= (/ b,a /)
       !elset%posgp(:,3)= (/ a,b /)

       ! matrix AP^(-T) for assumed shear at Gauss points
       DO g=1,ngaus
         !CALL shape7(elset%posgp(1,g),elset%posgp(2,g),elset%shape(1,g),deriv)
         CALL ap1tm7(elset%ap1(:,:,g),elset%posgp(1,g),elset%posgp(2,g))
       END DO

     ELSE
       !a = 1d0/3d0
       !elset%posgp(:,1) = (/ a,a /)
       !elset%weigp = 0.5d0
       !elset%shape(1:3,1) = a
       !elset%shape(4:6,1) = 0d0
       NULLIFY(elset%posgp, elset%weigp, elset%shape)
       elset%ap1 =  RESHAPE((/ -1d0/3d0,  1d0/3d0,   &
                               -1d0/3d0, -2d0/3d0,   &
                                2d0/3d0,  1d0/3d0/), (/2,elset%ngamm,1/))

     END IF
     nstre = 8
     IF( exists('ZIGZAG') ) THEN !local axis definition
       IF( ndofn == 8 )THEN
         elset%zigzag = .TRUE.
         nstre  = 14
       END IF
     END IF

   END IF

   CALL elmda7(nelem,nnode,elset%head,elset%tail,stype,iwrit,ngaus,elset%ngamm,nstre)
   CALL nodxy7(nelem,elset%head,coord,eule0,euler,ngaus)

   IF(cmpea)THEN
     ALLOCATE (lnods(6,nelem))
     el => elset%head
     DO i=1,nelem
       lnods(:,i) = el%lnods(1:6)
       el => el%next
     END DO
     IF( ngaus /= 1 )THEN
       CALL nodnor(nelem,nnode,lnods,coord,eule0,euler)
     ELSE
       CALL nodn07(nelem,lnods,coord,eule0,euler)  !locax
     END IF
     DEALLOCATE (lnods)
   END IF

   elset%plstr = 0     ! do not compute plastic strains
   IF (.NOT.oldset) CALL rdreqs ( 1 ,nreqs, elset%ngrqs, iwrit )

   CALL commv7(0, nelem,  nreqs, narch, stype, ngaus, elsnam, nnode, nstre, elset)
   ! add to the list of sets
   IF (.NOT.oldset) THEN
     CALL add_ele07 (elset, head, tail) ! add to the list of sets
     nelms = nelms + 1 ! increased set counter for this element type
   END IF

 ELSE IF (TRIM(task) == 'RESTAR') THEN

   ALLOCATE (elset)          !initializes a list
   NULLIFY(elset%head)       !nullify head pointer
   ! read control parameters
   elset%sname = elsnam
   READ (51) elset%nelem, elset%nreqs, elset%narch, elset%stype, elset%ngaus, elset%lside, &
             elset%gauss, elset%plstr, elset%angdf, elset%stabq, elset%nnode, elset%zigzag, elset%nstre
   ngaus = elset%ngaus
   IF( ngaus /= 1 )THEN
     elset%ngamm = 6
   ELSE
     elset%ngamm = 3
   END IF
   ALLOCATE(elset%ap1(2,elset%ngamm,ngaus))
   READ (51) elset%ap1

   IF( ngaus == 3 )THEN
!     gauss points in local coordinates and weigths
     ALLOCATE( elset%posgp(2,ngaus), elset%weigp(ngaus), elset%shape(nnode,ngaus))
     READ (51) elset%posgp, elset%weigp, elset%shape

   ELSE
     NULLIFY(elset%posgp, elset%weigp, elset%shape)

   END IF

   ! restore list of elements
   CALL rest07 (elset%nelem, elset%nreqs, elset%head, elset%tail, &
                elset%ngrqs, elset%stype, elset%ngaus, elset%ngamm, elset%nnode, elset%nstre)
   ! add to list of elements
   CALL add_ele07 (elset, head, tail)

 ELSE IF (TRIM(task) == 'IMPORT') THEN
   ! check if list of sets and set exists and initializes
   IF( overw )THEN
     CALL srch_ele07 (head, anter, elset, elsnam, oldset)
     IF (oldset) THEN    !if set exists
       CALL commv7(1, nelem,  nreqs, narch, stype, ngaus, elsnam, nnode, nstre, elset)
       !elset%lside = .FALSE.   !initializes flag to compute LSIDE
     ELSE
       CALL runen2(' Old set to overwrite does not exist')
     END IF
     READ (fimpo)
   ELSE
     CALL new_ele07(elset)       !reserve memory for set
     elset%sname = sname
     READ (fimpo) nelem,nnode,ngaus,stype,elset%angdf,elset%stabq,nstre
     nreqs = 0
     narch = 0
     IF( ngaus /= 1 )THEN
       elset%ngamm = 6
     ELSE
       elset%ngamm = 3
     END IF
     elset%locax = 3
   END IF
   ! restore list of elements
   CALL impo07 ( nelem,elset%head, elset%tail, ngaus, nnode, stype, elset%ngamm, nstre)
   CALL commv7 (0, nelem,  nreqs, narch, stype, ngaus, elsnam, nnode, nstre, elset)
   ! add to list of elements
   IF( .NOT.overw )THEN
     ALLOCATE(elset%ap1(2,elset%ngamm,ngaus))
     IF( ngaus == 3 )THEN
       !    gauss points in local coordinates and weigths
       ALLOCATE( elset%posgp(2,ngaus), elset%weigp(ngaus), elset%shape(nnode,ngaus))

       ! assumed membrane strain points

       a = 1d0/3d0
       CALL shape7(a  ,a  ,deriv(1,1),dn(1,1,1))
       CALL shape7(0d0,0d0,deriv(1,1),dn(1,1,2))
       CALL shape7(1d0,0d0,deriv(1,1),dn(1,1,3))
       CALL shape7(0d0,1d0,deriv(1,1),dn(1,1,4))

       ! ******************************
       ! gauss points in local coordinates and weigths
       ! version for Integration points at mid-side nodes
       a = 1d0/6d0                   !a = 1d0/6d0
       b = 1d0/2d0                   !b = 2d0/3d0
       elset%posgp(:,1) = (/ b  ,0d0 /)    !posgp(:,1)= (/ a,a /)
       elset%posgp(:,2) = (/ b  ,b   /)    !posgp(:,2)= (/ b,a /)
       elset%posgp(:,3) = (/ 0d0,b   /)    !posgp(:,3)= (/ a,b /)
       elset%weigp = a

       ! gauss points shape of nodal functions

       DO g=1,ngaus
         !CALL shape7(elset%posgp(1,g),elset%posgp(2,g),elset%shape(1,g),deriv)
         CALL ap1tm7(elset%ap1(:,:,g),elset%posgp(1,g),elset%posgp(2,g))
       END DO

     ELSE
       !a = 1d0/3d0
       !elset%posgp(:,1) = (/ a,a /)
       !elset%weigp = 0.5d0
       !elset%shape(1:3,1) = a
       !elset%shape(4:6,1) = 0d0
       NULLIFY(elset%posgp, elset%weigp, elset%shape)
       elset%ap1 =  RESHAPE((/ -1d0/3d0, -2d0/3d0,   &
                                2d0/3d0,  1d0/3d0,   &
                               -1d0/3d0,  1d0/3d0/), (/2,elset%ngamm,1/))

     END IF
     CALL add_ele07 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF

 ELSE
   CALL runend('INPD07: NON-EXISTENT TASK .        ')
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE inpda7
