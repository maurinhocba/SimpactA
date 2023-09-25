 SUBROUTINE celmn1p(ttime,coord)

 !.... perform contact searching & output report

 IMPLICIT NONE

 !     Dummy arguments
 REAL (kind=8) :: ttime
 REAL (kind=8), POINTER :: coord(:,:) !INTENT(IN)
 !     local variables
 INTEGER (kind=4) ipair,ims,iss,i,np,l,ncnod
 TYPE (pair1_db), POINTER :: pair
 TYPE (surf1_db), POINTER :: master,slave
 REAL (kind=8), POINTER :: nothing(:,:)
 LOGICAL :: flag
 REAL(kind=8), PARAMETER :: tol = -1e-9

 INTERFACE
   INCLUDE 'getnod.h'
 END INTERFACE

 IF( ASSOCIATED (nothing) )DEALLOCATE (nothing)
 ALLOCATE (nothing(1,1))


 pair => headp           !pointer to first pair
 DO ipair = 1, npair     !for each pair
   ! if pair is active for this time
   IF( pair%start <= ttime .AND. pair%end >= ttime)THEN
     IF( .NOT.pair%prev )THEN
       !.... identify master and slave surfaces numbers
       ims = pair%imast    !master surface order in list
       master => shead     !point to first surface
       DO i=1,ims-1        !loop until pointer is correct
         master => master%next
       END DO
       IF( master%nsegm > 0 )THEN
         pair%prev  = .TRUE.    !set Previous Step as TRUE
         iss = pair%islav    !slave surface order in list
         slave => shead      !point to first surface
         DO i=1,iss-1        !loop until pointer is correct
           slave => slave%next
         END DO
         !.... go into next routine
         CALL csrfc1(ipair,coord,pair,master,slave, &
                      1d0,nothing,nothing,nothing,.FALSE.)

         ! print data associated to this pair

         np = 0
         DO i=1,slave%ncnod
           IF( pair%rssdb(2,i) < tol )np = np+1
         END DO
         IF( np > 0 )THEN
           CALL openfi(57,rest=TRIM(rsfprj)//TRIM(rsfstg))   ! List of node penetrations
           WRITE(57,"(' Initial penetrations in pair:', a,/  &
                      '   Master Surface            :', a,/  &
                      '   Slave Surface             :', a,/  &
                      '   Number of penetrated nodes:',i8,/ &
                    )",ERR=9999) pair%pname,master%sname,slave%sname,np
           WRITE(*,"(/,' WARNING Initial Penetrations in pair: ',a)")TRIM(pair%pname)
           WRITE(57,"('   Node   Penetration'   )",ERR=9999)
           DO i=1,slave%ncnod
             IF( pair%rssdb(2,i) < tol )WRITE(57,*,ERR=9999) label(slave%lcnod(i)),pair%rssdb(2,i)
           END DO
           CLOSE(57)
         END IF
         IF( slave%nsegm > 0 .AND. print_data )THEN
           l = LEN_TRIM(pair%pname)
           OPEN(47,file=TRIM(rsfprj)//TRIM(rsfstg)//'-PAIR-'//pair%pname(1:l)//'.msh',form='FORMATTED')
           OPEN(48,file=TRIM(rsfprj)//TRIM(rsfstg)//'-PAIR-'//pair%pname(1:l)//'.res',form='FORMATTED')

           WRITE(47,*,ERR=9999) 'MESH "',TRIM(slave%sname),'"  dimension = 3 ElemType Linear      Nnode =  2'
           WRITE(47,*,ERR=9999) 'Coordinates'
           DO i=1,slave%ncnod
             l = slave%lcnod(i)
             WRITE(47,'(i10,3e15.6)',ERR=9999) label(l),coord(1:3,l)
           END DO
           WRITE(47,*,ERR=9999) 'End Coordinates'
           WRITE(47,*,ERR=9999) 'Elements'
           IF( ASSOCIATED (slave%lcseg) )THEN
             DO i=1,slave%nsegm
               WRITE(47,'(6i8)',ERR=9999) i,label(slave%lcseg(1:2,i))
             END DO
           ELSE
             DO i=1,slave%ncnod-1
               WRITE(47,'(6i8)',ERR=9999) i,label(slave%lcnod(i:i+1))
             END DO
             WRITE(47,'(6i8)',ERR=9999) i,label(slave%lcnod(slave%ncnod)),label(slave%lcnod(1))
           END IF
           WRITE(47,*,ERR=9999) 'end elements'

           IF( ASSOCIATED(master%lcnod) )THEN
             flag = .FALSE.
             ncnod = master%ncnod
           ELSE
             flag = .TRUE.
             l  = 2*master%nsegm+2
             ALLOCATE ( master%lcnod(l) )
             ncnod = 0
             CALL getnod(ncnod,l,master%nsegm,2,master%lcnod,master%lcseg)
           END IF
           WRITE(47,*,ERR=9999) 'MESH "',TRIM(master%sname),                            &
                                '"  dimension = 2 ElemType Linear      Nnode =  2'
           WRITE(47,*,ERR=9999) 'Coordinates'

           DO i=1,ncnod
             l = master%lcnod(i)
             WRITE(47,'(i10,3e15.6)',ERR=9999) label(l),coord(1:3,l)
           END DO
           WRITE(47,*,ERR=9999) 'End Coordinates'
           WRITE(47,*,ERR=9999) 'Elements'
           DO i=1,master%nsegm
             WRITE(47,'(6i8)',ERR=9999) i,label(master%lcseg(1:2,i))
           END DO
           WRITE(47,*,ERR=9999) 'end elements'
           IF(flag)DEALLOCATE (master%lcnod)
           WRITE(48,*,ERR=9999) 'GiD Post Results File 1.0'
           WRITE(48,*,ERR=9999) 'Result "Penetration" "Time Step "',ttime, ' Scalar OnNodes'
           WRITE(48,*,ERR=9999) 'ComponentNames "Gap"'
           WRITE(48,*,ERR=9999) 'Values'
           DO i=1,slave%ncnod
             l = slave%lcnod(i)
             WRITE(48,'(i8,e12.4)',ERR=9999) label(l),MIN(0d0,pair%rssdb(2,i))
           END DO
           WRITE(48,*,ERR=9999) 'End Values'
           WRITE(48,*,ERR=9999) 'Result "Distance" "Time Step "',ttime, ' Scalar OnNodes'
           WRITE(48,*,ERR=9999) 'ComponentNames "Distance"'
           WRITE(48,*,ERR=9999) 'Values'
           DO i=1,slave%ncnod
             l = slave%lcnod(i)
             WRITE(48,'(i8,e12.4)',ERR=9999) label(l),pair%rssdb(2,i)
           END DO
           WRITE(48,*,ERR=9999) 'End Values'
           CLOSE(47)
           CLOSE(48)
         END IF
       END IF
     END IF
   END IF
   pair => pair%next     !point to next pair
 END DO

 CLOSE(57)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE celmn1p
