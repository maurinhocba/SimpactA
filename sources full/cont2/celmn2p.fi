 SUBROUTINE celmn2p(ttime,coord,coorb,coort)

 !.... perform contact searching & output report

 IMPLICIT NONE

 !     Dummy arguments
 REAL (kind=8) :: ttime
 REAL (kind=8), POINTER :: coord(:,:),coorb(:,:),coort(:,:) !INTENT(IN)
 !     local variables
 INTEGER (kind=4) ipair,ims,iss,i,j,k,l,np,ncnod,maxdof
 REAL (kind=8), POINTER :: coors(:,:),coorm(:,:)
 TYPE (pair2_db), POINTER :: pair
 TYPE (surf2_db), POINTER :: master,slave
 REAL (kind=8), POINTER :: nothing(:,:)
 REAL (kind=8) :: disma,x(2,2),t1(2),rad,leng,area
 REAL (kind=8), PARAMETER :: twopi=6.283185307209586
 LOGICAL :: flag
 REAL(kind=8), PARAMETER :: tol = -1e-9


 INTERFACE
   INCLUDE 'getnod.h'
 END INTERFACE


 ALLOCATE (nothing(1,1))

  print *,' entro a celmn2p'
 pair => headp           !pointer to first pair
 DO ipair = 1, npair     !for each pair
   ! if pair is active for this time
   print *,' pair ',ipair
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
         master%cxc = .TRUE.    !set to recompute coordinates

         SELECT CASE (pair%mtsur) !select coordinates to use for master surface
         CASE(-2)
           coorm => coorb               !bottom
         CASE(-1,1)
           coorm => coord               !central
         CASE( 2)
           coorm => coort               !top
         END SELECT
         IF( pair%cutoff > -1d0 )THEN
           disma =  -10*pair%cutoff
         ELSE
           disma =  1d0
         END IF
         IF( pair%auto )THEN
           CALL csrfc2a(ipair,coorm,disma,1d0,pair,master, &
                        1d0,nothing,nothing,.FALSE.,0d0)

         ELSE
           iss = pair%islav    !slave surface order in list
           slave => shead      !point to first surface
           DO i=1,iss-1        !loop until pointer is correct
             slave => slave%next
           END DO
           IF( slave%ncnod > 0 )THEN
             SELECT CASE (pair%slsur) !select coordinates to use for slave surface
             CASE(-2)
               coors => coorb               !bottom
             CASE(-1,1)
               coors => coord               !central
             CASE( 2)
               coors => coort               !top
             END SELECT
             !.... go into next routine
             CALL csrfc2(ipair,coors,coorm,disma,1d0,pair,master,slave, &
                         1d0,nothing,nothing,.FALSE.,0d0)
           END IF
         END IF

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

           WRITE(47,*,ERR=9999) 'MESH "',TRIM(slave%sname),'"  dimension = 2 ElemType Linear      Nnode =  2'
           WRITE(47,*,ERR=9999) 'Coordinates'
           DO i=1,slave%ncnod
             l = slave%lcnod(i)
             WRITE(47,'(i10,3e15.6)',ERR=9999) label(l),coors(1:2,l)
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

           IF( .NOT.pair%auto )THEN
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
               WRITE(47,'(i10,3e15.6)',ERR=9999) label(l),coorm(1:2,l)
             END DO
             WRITE(47,*,ERR=9999) 'End Coordinates'
             WRITE(47,*,ERR=9999) 'Elements'
             DO i=1,master%nsegm
               WRITE(47,'(6i8)',ERR=9999) i,label(master%lcseg(1:2,i))
             END DO
             WRITE(47,*,ERR=9999) 'end elements'
             IF(flag)DEALLOCATE (master%lcnod)
           END IF
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
         ! for thermal analysis, compute associated area
         ! there is problem here with the START-END times
         IF( therm )THEN
           slave%area = 0d0 !initializes
           DO i=1,slave%nsegm
             x(:,1:2) = coord(:,slave%lcseg(:,i)) !coordinates
             t1 = x(:,2) - x(:,1)                 !segment n1 -> n2
             leng = SQRT(t1(1)*t1(1)+t1(2)*t1(2)) !segment longitude
             area = leng/2d0 ! in plain strain A = (L * 1d0)/2             
             IF( ntype == 3 )THEN
               rad  = (x(1,1) + x(1,2))/2d0 !mid segment radius
               area = area*twopi*rad !in axilsymmetric problem A = (L * 2 pi R)/2
             END IF  
             DO j=1,2
               k = slave%lcseg(j,i)
               DO l=1,slave%ncnod
                 IF( k == slave%lcnod(l))THEN
                   slave%area(l) = slave%area(l) + area
                   EXIT
                 END IF
               END DO
             END DO
           END DO
           DO i=1,slave%ncnod
             IF( slave%area(i) > 0 )slave%area(i) = 1d0/slave%area(i)
           END DO
           !  compute Temperature DOFs for master surface
           maxdof = 1
           DO k=2,ndoft
             !DO i=1,master%nsegm
             !  IF( iftmp(k,master%lcseg(1,i)) /= 0 ) maxdof = MAX(k,maxdof)
             !END DO
             IF( iftmp(k,master%lcseg(1,1)) /= 0 ) maxdof = MAX(k,maxdof)
           END DO
           pair%mtdof = maxdof
           IF( maxdof > 1 .AND. pair%mtsur < 0 ) pair%mtdof = maxdof-1
           !  compute Temperature DOFs for slave surface
           maxdof = 1
           DO k=2,ndoft
             !DO i=1,slave%ncnod
             !  IF( iftmp(k,slave%lcnod(i)) /= 0 ) maxdof = MAX(k,maxdof)
             !END DO
             IF( iftmp(k,slave%lcnod(1)) /= 0 ) maxdof = MAX(k,maxdof)
           END DO
           pair%sldof = maxdof
           IF( maxdof > 1 .AND. pair%slsur < 0 ) pair%sldof = maxdof-1
         END IF
       END IF
     END IF
   END IF
   pair => pair%next     !point to next pair
 END DO

 CLOSE(57)
 DEALLOCATE (nothing)
 RETURN
 9999 CALL runen2('')
 END SUBROUTINE celmn2p
