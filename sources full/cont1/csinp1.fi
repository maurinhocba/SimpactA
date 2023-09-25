 SUBROUTINE csinp1(maxve,iwrit,oldsr,oldpr, &
                   codes,surname,tsurf,surdi)

 !.... input and generate contact surfaces data

 IMPLICIT NONE
 !     arguments
 INTEGER (kind=4), INTENT(IN) :: iwrit,oldpr,maxve,surdi
 INTEGER (kind=4), INTENT(IN OUT) :: oldsr,tsurf,codes(2,surdi)
 CHARACTER(len=*),INTENT(INOUT):: surname(surdi)

 !     local variables
 LOGICAL :: found,initia
 INTEGER (kind=4) :: j,k,isurf,deflt,ipair
 TYPE (pair1_db), POINTER :: pair
 TYPE (surf1_db), POINTER :: surf, anter, posic, mast
 CHARACTER(len=mnam):: sname

 !.... Check if Old Surfaces will be used by New Pairs

 pair => headp                  !point to first pair
 DO ipair=1,oldpr               !skip OLD pair
   pair => pair%next
 END DO
 DO ipair=oldpr+1,npair         !for each new-pair
   DO j=1,2                     !for each surface
     IF( j == 1 )THEN           ! master surface (j=1)
       sname = pair%master                          !master surface label
     ELSE                       ! slave surface  (j=2)
       sname = pair%slave                           !slave surface label
     END IF
     k = 0                      !initializes surface position in list
     surf => shead              !point to first (old) surface
     DO
       k=k+1                    !increase counter
       IF( k > oldsr)EXIT       !all existing surfaces considered EXIT
       IF( TRIM(surf%sname) == TRIM(sname) )THEN  !if surface found
         IF( j == 1)THEN      !master surface
           IF( surf%imcod )THEN     !surface was supposed master
             pair%imast = k            !O.K., assign internal position
             EXIT  !exit surface search
           ELSE
             WRITE(lures,"('Surface',A,' Cannot be MASTER')",ERR=9999) TRIM(sname)
             CALL runend('CSINP1: Surface Cannot be MASTER   ')
           END IF
         ELSE                 !slave surface
           IF( surf%iscod )THEN      !surface was supposed slave
             pair%islav = k             !O.K., assign internal position
             pair%ncnod = surf%ncnod    !O.K.
             ALLOCATE ( pair%issdb(nisdb,surf%ncnod) ) !slave Int. DB
             ALLOCATE ( pair%rssdb(nrsdb,surf%ncnod) ) !slave Real DB
             pair%issdb = 0                            !initializes
             pair%rssdb = 0d0                          !initializes
             IF( pair%press )THEN                      !compute normal force
               ALLOCATE( pair%presn(surf%ncnod) )      !reserve space
               pair%presn = 0d0                        !initializes
             END IF
             EXIT  !exit surface search
           ELSE
             WRITE(lures,"('Surface',A,' Cannot be SLAVE')",ERR=9999) TRIM(sname)
             CALL runend('CSINP1: Surface Cannot be SLAVE    ')
           END IF
         END IF
       END IF
       surf => surf%next    !point to next surface
     END DO  !search in surface Data Base
   END DO ! j=1,2
   pair => pair%next        !consider next pair
 END DO !ipair=oldpr+1,npair

 !.... loop over new contact surfaces (READ)

 DO
   ! loop over surfaces

   CALL listen('CSINP1')
   IF (.NOT.exists('ISURF ')) THEN
     backs = .TRUE.
     EXIT
   END IF

   CALL new_surf1 (surf)              !allocate new surface
   surf%sname = get_name('ISURF ',found, ' CONTACT SURFACE: ', stype='SURF')

   !check if surface label already used
   CALL srch1_srf (shead, anter, posic, surf%sname, found)
   IF (found) CALL runend ('CSINP1: Surface already defined !!!')

   isurf= 0        !search for surface label in CODES array
   DO
     isurf = isurf+1
     IF( isurf > tsurf .AND. isurf < surdi)THEN   !not defined but possible
       surname(isurf) = surf%sname
       WRITE(lures,"(' WARNING: surface ',A,' not defined in any ', &
              & 'pair')",ERR=9999) TRIM(surf%sname)
     ELSE IF(TRIM(surname(isurf)) == TRIM(surf%sname))THEN    !found
       EXIT
     ELSE IF(isurf == surdi)THEN        !space exhausted
       CALL runend('CSINP1: More surfaces than expected')
     END IF
   END DO

   CALL listen('CSINP1')
 !                          read flag WRITE
   deflt = 0
   IF(codes(1,isurf) > 0 .AND. codes(2,isurf) == 0)deflt = 1
   surf%iwrit =getint('WRITE ',deflt, ' Save surface for Post-Processing..') == 1
   IF(iwrit > 0)THEN
     IF( surf%iwrit )THEN
       WRITE(lures,"(12x,'WILL be saved for post-Processing ')",ERR=9999)
     ELSE
       WRITE(lures,"(12x,'WILL NOT be saved for post-Processing')",ERR=9999)
     END IF
   END IF

 !                          read flag NCNOD
   deflt = 0
   IF(codes(1,isurf) == 0 .AND. codes(2,isurf) /= 0)deflt = 1
   surf%ncnod =getint('NCNOD ',deflt, &
                      ' Nodes in the surface will be read ')
   IF(iwrit > 0)THEN
     IF( surf%ncnod > 0)THEN
       WRITE(lures,"(12x,'List of node will be read ')",ERR=9999)
     ELSE
       WRITE(lures,"(12x,'No list of nodes expected')",ERR=9999)
     END IF
   END IF

 !                          read flag NSEGM
   deflt = 0
   IF( surf%iwrit ) deflt = 1
   IF( codes(1,isurf) > 0 .OR. surf%ncnod == 0) deflt = 1
   surf%nsegm =getint('NSEGM ',deflt,' Segments of surface will be read..')
   IF(iwrit > 0)THEN
     IF( surf%nsegm > 0)THEN
       WRITE(lures,"(12x,'Connectivities will be read ')",ERR=9999)
     ELSE
       WRITE(lures,"(12x,'No connectivities expected')",ERR=9999)
     END IF
   END IF

 !                          read flag IMCOD
   deflt = 0
   IF( codes(1,isurf) > 0 ) deflt = 1
   surf%imcod =getint('IMCOD ',deflt,' Code for master surface ..........') == 1
   IF( surf%imcod .AND. surf%nsegm == 0 ) &
        CALL runend ('CSINPT:IMCOD=1 requires nsegm > 0  ')
   IF(iwrit > 0)THEN
     IF( surf%imcod )THEN
       WRITE(lures,"(12x,'Will act as MASTER for some pair')",ERR=9999)
     ELSE
       WRITE(lures,"(12x,'Will not act as MASTER for any pair')",ERR=9999)
     END IF
   END IF
   IF( .NOT.surf%imcod .AND.  deflt == 1 )THEN
     WRITE(lures,"(' Bad Master Code for surface',A)",ERR=9999) TRIM(surf%sname)
     CALL runend('CSINP1: Master Code Must be 1 !!!   ')
   ELSE IF(surf%imcod .AND. deflt == 0 )THEN
     WRITE(lures,"(' Warning Master Code for surface ',A,' set= 1', &
     & /,t15,'but is not as master in any pair')",ERR=9999) TRIM(surf%sname)
   END IF

 !                          read flag ISCOD
   deflt = 0
   IF( codes(2,isurf) > 0 ) deflt = 1
   surf%iscod =getint('ISCOD ',deflt,' Code for slave surface ...........') == 1
   IF(iwrit > 0)THEN
     IF( surf%iscod )THEN
       WRITE(lures,"(12x,'Will act as SLAVE for some pair')",ERR=9999)
     ELSE
       WRITE(lures,"(12x,'Will not act as SLAVE for any pair')",ERR=9999)
     END IF
   END IF

   IF( .NOT.surf%iscod .AND. deflt == 1 )THEN
     WRITE(lures,"(' Bad SLAVE Code for surface',A)",ERR=9999) TRIM(surf%sname)
     CALL runend('CSINPT: SLAVE Code Must be 1 !!!    ')
   ELSE IF( surf%iscod .AND. deflt == 0 )THEN
     WRITE(lures,"(' Warning SLAVE Code for surface',A,' set= 1', &
     & /,t15,'but is not as SLAVE in any pair')",ERR=9999) TRIM(surf%sname)
   END IF

   ! warning if WRITE = 1 and NSEGM = 0
   IF( surf%iwrit .AND. surf%nsegm == 0)THEN
     WRITE(lures,"(' Sorry, NSEGM must be > 0 to save surface for ' &
                 &  'visualization')",ERR=9999)
     surf%iwrit = .FALSE.
   END IF

   surf%density =getrea('DENSIT', 0d0 ,' line density for mass comput...')
   surf%diam    =getrea('DIAMET', 0d0 ,'!surface diameter ..............')

   !surf%curved = exists('CURVED')
   !IF( surf%curved ) &
   !    WRITE(lures,"(12x,'Surface will be treated as CURVED')",ERR=9999)

 !                          read flag INITIA
   deflt = 0
   !IF( codes(2,isurf) > 0 ) deflt = 1
   initia =getint('INITIA',deflt,' Compute intial beam displacements ') == 1
   IF( initia .AND. surf%nsegm == 0 ) &
        CALL runend ('CSINPT:INITIA=1 requires nsegm > 0  ')
   IF( initia .AND. .NOT.surf%iscod ) &
        CALL runend ('CSINPT:INITIA=1 can only be required for slave surfaces')
   IF(iwrit > 0)THEN
     IF( initia )THEN
       WRITE(lures,"(12x,'Inital displacements will be computed')",ERR=9999)
     ELSE
       WRITE(lures,"(12x,'Inital displacements will be computed')",ERR=9999)
     END IF
   END IF

   ! read surface connectivities and generate data base
   CALL csdat1(maxve,iwrit,surf)

   nsurf = 1 + nsurf    !increase number of surfaces

   !Check if the surface will be used as SLAVE for some pairs

   pair => headp               !point to first pair
   DO ipair=1,codes(2,isurf)   !codes(2) = times used as slave
     DO
       IF( TRIM(pair%slave) == TRIM(surname(isurf)) )EXIT    !found
       pair => pair%next                         !point to next pair
     END DO
     pair%islav = nsurf               !assign internal position
     pair%ncnod = surf%ncnod          !assign number of nodes
     ALLOCATE( pair%issdb(nisdb,surf%ncnod),  &   !reserve space for
               pair%rssdb(nrsdb,surf%ncnod) )     !slave data base
     pair%issdb = 0                   !initializes Int. array
     pair%rssdb = 0d0                 !initializes Real array
     IF( pair%press )THEN
       ALLOCATE( pair%presn(pair%ncnod) )   !reserve space
       pair%presn = 0d0                !initialize
     END IF
     pair => pair%next                !next pair
   END DO

   !Check if the surface will be used as MASTER for some pairs

   pair => headp               !point to first pair
   DO ipair=1,codes(1,isurf)   !codes(1) = times used as master
     DO
       IF( TRIM(pair%master) == TRIM(surname(isurf)) )EXIT   !found
       pair => pair%next                         !point to next pair
     END DO
     pair%imast = nsurf               !assign internal position

     pair => pair%next                   !point to next pair
   END DO
   !IF( surf%curved )CALL curve1(coord,surf,npoin)
   IF( initia )surf%ncnod = -surf%ncnod
   CALL add1_srf (surf, shead, stail)    !add surface to data_base
   !Compute initial displacements if requested
 END DO         !surface reading

 surf => shead
 DO j=1,oldsr
   surf => surf%next
 END DO
 DO j=oldsr+1,nsurf
   initia = surf%ncnod < 0
   IF( initia )THEN
     surf%ncnod = -surf%ncnod
     ! find master surface (skin)
     pair => headp               !point to first pair
     DO
       IF( TRIM(pair%slave) == TRIM(surf%sname) ) THEN
         mast => shead
         DO
           IF( TRIM(pair%master) == TRIM(mast%sname) ) EXIT
           mast => mast%next
         END DO
         EXIT   !found
       END IF
       pair => pair%next                         !point to next pair
     END DO
     CALL beam_ini(surf%nsegm,mast%nsegm,surf%lcseg,mast%lcseg,iwrit)
   END IF
   surf => surf%next
 END DO
 ! check if all the pairs have surface correctly defined
 pair => headp               !point to first pair
 found = .FALSE.
 DO ipair=1,npair            !for each pair
   IF( pair%imast < 1 .OR. pair%imast > nsurf )THEN
     WRITE(lures,"(' master surface ',A,' asssociated to pair ',A,  &
&                  ' does NOT exists')",ERR=9999) TRIM(pair%master), TRIM(pair%pname)
     WRITE(55,"(' master surface',A,' associated to pair',A,     &
&                  ' does NOT exists')",ERR=9999) TRIM(pair%master), TRIM(pair%pname)
     found = .TRUE.
   END IF
   IF( pair%islav < 1 .OR. pair%islav > nsurf )THEN
     WRITE(lures,"(' slave surface ',A,' associated to pair ',A,  &
&                  ' does NOT exists')",ERR=9999) TRIM(pair%slave), TRIM(pair%pname)
     WRITE(55,"(' slave surface ',A,' associated to pair ',A,     &
&                  ' does NOT exists')",ERR=9999) TRIM(pair%slave), TRIM(pair%pname)
     found = .TRUE.
   END IF

   pair => pair%next                   !point to next pair
 END DO
 IF( found )CALL runend('CSINP1: contact surface missing   ')

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE csinp1
