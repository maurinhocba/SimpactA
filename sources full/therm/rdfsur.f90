 SUBROUTINE rdfsur(iwrit,flag)

 !     Reads free element faces with convection or radiation (segments or surfaces)

 USE c_input
 USE heat_db
 USE surf_db
 USE npo_db, ONLY : iftmp
 USE ctrl_db, ONLY : ndoft
 USE curv_db, ONLY: getcun
 IMPLICIT NONE
 ! dummy arguments
 INTEGER (kind=4) :: iwrit  !echo flag

 TYPE (srf_con),POINTER :: s,anter,posic
 LOGICAL :: found
 INTEGER (kind=4) :: flag   !task flag
 CHARACTER (len=mnam) :: cname
 INTEGER (kind=4) :: i,j,nodes(3),chnode,maxdof
 TYPE (srf_seg),POINTER :: elm

 INTERFACE
   INCLUDE 'elemnt.h'
 END INTERFACE

 !----------------------------------------------------------------------------------
 IF( flag == 1 )THEN !read data for convection and radiation surfaces

   ncsur = 0
   !***    convection/radiation surfaces    ***

   DO  !loop to read the convection/radiation surfaces

     CALL listen ('RDFSUR')

     IF( .NOT.exists('SURF  ')) THEN
       IF( .NOT.exists('ENDCON')) CALL runend ('RDFSUR: ENDCON keyword expected')
       EXIT
     ELSE
       ! read surface convection/radiation coefficients and temperatures
       !----------------------------------------------------------------
       ALLOCATE(s)
       s%nn    = getint('NN',2,' Number of Nodes in surface.......')
       s%sname = get_name('ISURF ',found, ' CONVECT SURFACE: ',stype='SURF')
       s%afreq = getint('AFREQ',100000,' Interval to compute areas........')
       s%solid = .FALSE. !initialize

       !check if surface label already used
       CALL srch_csr (headsc, anter, posic, s%sname, found)
       IF (found) CALL runend ('RDFSUR: Surface already defined !!!')

       s%tcur = 0
       s%lbl  = ''
       s%tf   = 0d0
       CALL listen ('RDFSUR')

       IF(exists('HCONV'))THEN  ! read convection data for solid problems
         s%solid = .TRUE.
         s%hf(1)=getrea('HCONV',0d0,' Convection factor for surface......')
         s%conve= s%hf(1) /= 0d0

         IF(exists('TAMBN',i)) THEN
            IF(words(i+1) == 'CURVE' )THEN
              s%tcur(1) = 1
              s%lbl(1) = get_name(posin=i+1,stype='CURV')         !get curve name
              WRITE (lures,"(9X,A34,'.... ',A6,' = ',(A))",ERR=9999) &
             &        'Ambient temperature defined by curve...','TAMBN',s%lbl(1)
            ELSE
              s%tf(1) =getrea('TAMBN',0d0,' Ambient temperature for surface....')
            END IF
         END IF

       ELSE ! read convection data for plate/shell problems

         s%hf(1)=getrea('HCONC',0d0,' Convection factor Central surface..')
         IF( s%hf(1) /= 0d0) THEN  !read ambient temperature
           IF(exists('TAMBC',i)) THEN
              IF(words(i+1) == 'CURVE' )THEN
                s%tcur(1) = 1
                s%lbl(1) = get_name(posin=i+1,stype='CURV')         !get curve name
                WRITE (lures,"(9X,A34,'.... ',A6,' = ',(A))",ERR=9999) &
               &        'Ambient temp. (Central) defined by curve','TAMBC',s%lbl(1)
              ELSE
                s%tf(1) =getrea('TAMBC',0d0,' Ambient temp. for Central surface..')
              END IF
            END IF
         END IF

         s%hf(2)=getrea('HCONT',0d0,' Convection factor Top surface......')
         IF( s%hf(2) /= 0d0) THEN  !read ambient temperature (Top)
           IF(exists('TAMBT',i)) THEN
              IF(words(i+1) == 'CURVE' )THEN
                s%tcur(2) = 1
                s%lbl(2) = get_name(posin=i+1,stype='CURV')         !get curve name
                WRITE (lures,"(9X,A34,'.... ',A6,' = ',(A))",ERR=9999) &
               &        'Ambient temp. (TOP) defined by curve    ','TAMBT',s%lbl(2)
              ELSE
                s%tf(2) =getrea('TAMBT',0d0,' Ambient temp. for Top surface......')
              END IF
            END IF
         END IF

         s%hf(3)=getrea('HCONB',0d0,' Convection factor Bottom surface...')
         IF( s%hf(3) /= 0d0) THEN  !read ambient temperature (Bottom)
           IF(exists('TAMBB',i)) THEN
              IF(words(i+1) == 'CURVE' )THEN
                s%tcur(3) = 1
                s%lbl(3) = get_name(posin=i+1,stype='CURV')         !get curve name
                WRITE (lures,"(9X,A34,'.... ',A6,' = ',(A))",ERR=9999) &
               &        'Ambient temp. (Bottom) defined by curve ','TAMBB',s%lbl(3)
              ELSE
                s%tf(3) =getrea('TAMBB',0d0,' Ambient temp. for Bottom surface...')
              END IF
            END IF
         END IF

         s%conve= ANY(s%hf(:) /= 0d0)

       END IF

       IF(exists('HRADS'))THEN  !read radiation data for solid problems
         s%solid = .TRUE.
         s%hf(4) = getrea('HRADS',0d0,' Radiation factor for surface.......')
         s%radia = s%hf(4) /= 0d0

         IF(exists('TRADS',i)) THEN
            IF(words(i+1) == 'CURVE' )THEN
              s%tcur(4) = 1
              s%lbl(4) = get_name(posin=i+1,stype='CURV')         !get curve name
              WRITE (lures,"(9X,A34,'.... ',A6,' = ',(A))",ERR=9999) &
             &        'Radiation temperature defined by curve.','TRADS',s%lbl(1)
            ELSE
              s%tf(4) =getrea('TRADS',0d0,' Radiation temperature for surface..')
            END IF
         END IF

       ELSE ! read radiation data for plate/shell problems

         s%hf(4) =getrea('HRADC',0d0,' Radiation factor Central surface...')
         IF( s%hf(4) /= 0d0) THEN  !read radiation temperature
           IF(exists('TRADC',i)) THEN
              IF(words(i+1) == 'CURVE' )THEN
                s%tcur(4) = 1
                s%lbl(4) = get_name(posin=i+1,stype='CURV')         !get curve name
                WRITE (lures,"(9X,A34,'.... ',A6,' = ',(A))",ERR=9999) &
               &        'Radiation Temp.(Central) defined by curve','TRADC',s%lbl(4)
              ELSE
                s%tf(4) =getrea('TRADC',0d0,' Radiation temp. for Central surface..')
              END IF
            END IF
         END IF

         s%hf(5) =getrea('HRADT',0d0,' Radiation factor TOP     surface...')
         IF( s%hf(5) /= 0d0) THEN  !read radiation temperature
           IF(exists('TRADT',i)) THEN
              IF(words(i+1) == 'CURVE' )THEN
                s%tcur(5) = 1
                s%lbl(5) = get_name(posin=i+1,stype='CURV')         !get curve name
                WRITE (lures,"(9X,A34,'.... ',A6,' = ',(A))",ERR=9999) &
               &        'Radiation Temp.(TOP) defined by curve    ','TRADT',s%lbl(5)
              ELSE
                s%tf(5) =getrea('TRADT',0d0,' Radiation temp. for Top surface......')
              END IF
            END IF
         END IF

         s%hf(6) =getrea('HRADB',0d0,' Radiation factor Bottom  surface...')
         IF( s%hf(6) /= 0d0) THEN  !read radiation temperature
           IF(exists('TRADB',i)) THEN
              IF(words(i+1) == 'CURVE' )THEN
                s%tcur(6) = 1
                s%lbl(6) = get_name(posin=i+1,stype='CURV')         !get curve name
                WRITE (lures,"(9X,A34,'.... ',A6,' = ',(A))",ERR=9999) &
               &        'Radiation Temp.(Bottom) defined by curve ','TRADB',s%lbl(6)
              ELSE
                s%tf(6) =getrea('TRADB',0d0,' Radiation temp. for Central surface..')
              END IF
            END IF
         END IF

         s%radia = ANY(s%hf(4:6) /= 0d0)

       END IF

       ! read curve data for temperatures
       !----------------------------------------------------------------
       IF(ANY(s%tcur /= 0 )) THEN  !read curve data
         CALL listen ('RDFSUR')
         IF( exists('CURVED',i))THEN
           cname = get_name(posin=i,stype='CURV')     !get curve name
           found = .FALSE.
           DO i=1,6
             IF( s%tcur(i) /= 1 )CYCLE
             IF( TRIM(cname) == TRIM(s%lbl(i)) )THEN
               found = .TRUE.
               s%tcur(i) = -1
               EXIT
             END IF
           END DO
           WRITE(55,"('Curve ',a,' is not associated to any temperature')")cname
           IF( .NOT.found )CALL runend('RDFSUR: Defined Curved not associated')
           backs = .TRUE.
           CALL rdcurv('FORCE',cname)
         ELSE
           backs = .TRUE.
         END IF
       END IF
       IF( ANY(s%tcur > 0 )) THEN
         DO i=1,6
           IF( s%tcur(i) == 1 )WRITE(55,"('Curve undefined',a)")s%lbl(i)
         END DO
         CALL Runend('RDFSUR: There are undefined curves')
       END IF

       !*****    read surface segments definition    *****
       CALL listen('RDFSUR')
       CALL new_srf(surfa)
       IF (.NOT.exists('ELEMEN'))THEN
         backs = .TRUE.
         ! surface definition based on an element set
         found= .FALSE.
         CALL elemnt ('SURFAC',name=s%sname,flag2=found)
         IF (.NOT.found) CALL runend('RDFSUR:ELEMENT SET NOT FOUND       ')

       ELSE  !surface definition based on an node conectivities list
         IF(iwrit == 1 ) WRITE(lures,"(/6x,'List of Faces read in this strategy'/)",ERR=9999)
         s%nelem = 0
         DO
           CALL listen('RDFSUR')
           IF (exists('ENDELE'))EXIT
           IF( nnpar /= s%nn) THEN
             WRITE(lures,"(' nnpar',i5,' nnseg',i5)",ERR=9999) nnpar,s%nn
             CALL runend('RDFSUR:erroneous number of nodes   ')
           END IF
           DO i=1,nnpar
             nodes(i) = INT(param(i))
             nodes(i) = chnode(nodes(i))
           END DO
           IF( iwrit == 1)WRITE(lures,"(3i5)")nodes(1:nnpar)
           s%nelem= s%nelem+1
           CALL new_seg(elm)
           elm%nodes(1:s%nn) = nodes(1:s%nn)
           CALL add_seg(elm,surfa%head,surfa%tail)
         END DO
         !CALL add_name(s%sname,5)
       END IF
       CALL listen ('RDFSUR')
       IF (.NOT.exists('ENDSUR'))CALL runend('RDFSUR: END_SURFACE expected   ')
       s%nelem = surfa%nelem
       ALLOCATE( s%lnods(s%nn,s%nelem), s%area(s%nelem) )
       CALL store_segs(surfa%head,s%lnods,s%nn,s%nelem)
       CALL dalloc_srf(surfa)
       s%carea = .TRUE.
       CALL add_csr (s, headsc, tailsc)
       ncsur = ncsur + 1

     END IF
   END DO

 !----------------------------------------------------------------------------------
 ELSE   !check surface ONLY for shell/plate problems
   s => headsc
   LS: DO
     IF( .NOT.ASSOCIATED(s) )EXIT
     DO i=1,ndoft  !find the number of DOFs per node involved
       DO j=1,s%nelem
         IF( ANY(iftmp(i,s%lnods(:,j)) > 0 ) )maxdof = i
       END DO
     END DO
     s%mdof = maxdof !maximum number of thermal DOFs
     IF( s%solid )THEN ! in solid problems no checks are necesary
       s => s%next     ! point to next surface
       CYCLE LS        ! cycle loop over convection/radiation surface
     END IF
     ! for convection, check central or top surfaces has adequate values
     IF( s%conve )THEN
       IF( maxdof == 1 )THEN  !only mid-surface
         IF( s%hf(1) == 0d0 ) s%hf(1) = s%hf(2)
         IF( s%tf(1) == 0d0 .AND. s%tcur(1) == 0 )THEN
           s%tf(1)  = s%tf(2)
           s%tcur(1) = s%tcur(2)
           IF( s%tcur(1) /= 0 ) s%lbl(1) = s%lbl(2)
         END IF
       ELSE                   !shells
         IF( s%hf(2) == 0d0 ) s%hf(2) = s%hf(1)
         IF( s%tf(2) == 0d0 .AND. s%tcur(2) == 0 )THEN
           s%tf(2)  = s%tf(1)
           s%tcur(2) = s%tcur(1)
           IF( s%tcur(2) /= 0 ) s%lbl(2) = s%lbl(1)
         END IF
       END IF
     END IF
     ! for radiation, check central or top surfaces has adequate values
     IF( s%radia )THEN
       IF( maxdof == 1 )THEN  !only mid-surface
         IF( s%hf(4) == 0d0 ) s%hf(4) = s%hf(5)
         IF( s%tf(4) == 0d0 .AND. s%tcur(4) == 0 )THEN
           s%tf(4)  = s%tf(5)
           s%tcur(4) = s%tcur(5)
           IF( s%tcur(4) /= 0 ) s%lbl(4) = s%lbl(5)
         END IF
       ELSE                   !shells
         IF( s%hf(5) == 0d0 ) s%hf(5) = s%hf(4)
         IF( s%tf(5) == 0d0 .AND. s%tcur(5) == 0 )THEN
           s%tf(5)  = s%tf(4)
           s%tcur(5) = s%tcur(4)
           IF( s%tcur(5) /= 0 ) s%lbl(5) = s%lbl(4)
         END IF
       END IF
     END IF
     DO i=1,6
       IF(s%tcur(i) == -1) s%tcur(i)=getcun (s%lbl(i))    ! get associated curve number for this set
     END DO
     s => s%next
   END DO LS
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE RDFSUR
