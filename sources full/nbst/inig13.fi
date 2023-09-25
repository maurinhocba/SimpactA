 SUBROUTINE inig13(heade,nelem,iwrit,plstr)

 !read initial internal variables
 !THIS ROUTINE MUST BE REVIEWED FOR NEW LOCAL SYSTEMS AND MATERIAL

 IMPLICIT NONE
 INTEGER(kind=4), INTENT(IN) :: nelem,iwrit
 INTEGER(kind=4), INTENT(OUT) :: plstr
 TYPE (ele13), POINTER  :: heade

 INTEGER (kind=4), PARAMETER :: nstre = 3
 INTEGER (kind=4) :: ie,i1,numel,isec,nlayr,i,j,il,osec,nvare,ib
 TYPE (ele13), POINTER :: e
 TYPE (section), POINTER :: sec

 CALL listen('INIG13')          !read line

 IF (exists('ZEROAL')) THEN     !key word ZERO_ALL found
   !Zero stresses and Eff. Pl Strain
   e => heade                   !point to first element
   DO ie=1,nelem                !for each element in the set
     IF( ASSOCIATED( e%gausv) ) THEN
        DEALLOCATE( e%gausv )
        NULLIFY (e%gausv)
     END IF
     e => e%next                !point to next element
   END DO
   plstr = 1                    ! compute plastic strains
 ELSE IF( exists('INTERN') )THEN  !Read initial values for Int. vars.
   IF( iwrit > 0 )WRITE(lures,"(' ELEMENT LAY     EPXX',11x,'EPYY', &
&                               11x,'EPXY',11x,'EFPST')",ERR=9999)
   osec = -1
   e => heade                            !position at top of list
   ie = 1                                !initializes element counter
   DO
     CALL listen('INIG13')          !read a card
     IF( exists('ENDINT'))EXIT      !key-word END_INTERNAL found, exit loop
     IF( exists('NUMEL ',j))THEN    !key-word NUMEL found
       numel = INT(param(j))        !read associated element number
     ELSE
       numel = INT(param(1))        !first parameter is element number
     END IF
     IF(iwrit > 0)WRITE(lures,"(i8)",ERR=9999) numel

     i1 = ie        !remember first element
     DO             !loop search for element with number NUMEL
       IF( e%numel == numel )EXIT          ! element found
       ie = ie+1                           ! next element
       IF( ie == i1 )THEN                  ! if loop completed => ERROR
         WRITE(lures,"(' element label not found ',i6)",ERR=9999) numel
         CALL runend('INIG13: Element label not found !!.')
       END IF
       IF( ie > nelem )THEN                !last element passed
         ie = 1                            !go to top of the list
         e => heade                        !point to first element
       ELSE
         e => e%next                       !go to next element in the list
       END IF
     END DO
     ie = ie + 1                         !correct element order
     isec = e%matno                     !material type
     IF( isec /= osec )THEN
       !IF(matdb(isec)%mt%conso == 0d0) CYCLE   !only for elasto-plastic
       sec => psecs(isec)%p
       nlayr = sec%iprop(1)                     !number of layers
       nvare = sec%iprop(2)                     !number of components
       osec = isec
     END IF
     IF( .NOT.ASSOCIATED(e%gausv)) ALLOCATE ( e%gausv(nvare,nlayr) )
     ib = 1
     DO i=1,nlayr                        !for each layer
       IF( exists('LAYER ',j))THEN       !key-word LAYER found
         il = INT(param(j))           !read associated layer number
       ELSE
         il = INT(param(1+ib))           !layer number is first parameter
       END IF
       IF( i /= il )      &    !layers must be in order
               CALL runend('INIG13: Invalid layer number.      ')
       IF(exists('EPXX  ',j))THEN    !key-word EPXX found
         e%gausv(1,i) = param(j)     !read associated value
       ELSE
         e%gausv(1,i) = param(2+ib)     !EPXX is 2nd parameter
       END IF
       IF(exists('EPYY  ',j)) THEN   !key-word EPYY found
         e%gausv(2,i) = param(j)     !read associated value
       ELSE
         e%gausv(2,i) = param(3+ib)     !EPYY is 3rd parameter
       END IF
       IF(exists('EPXY  ',j)) THEN   !key-word EPXY found
         e%gausv(3,i) = param(j)     !read associated value
       ELSE
         e%gausv(3,i) = param(4+ib)     !EPXY is 4th parameter
       END IF
       IF(exists('EFPST ',j)) THEN   !key-word EFPST found
         e%gausv(4,i) = param(j)     !read associated value
       ELSE
         e%gausv(4,i) = param(5+ib)     !EFPST is 5th parameter
       END IF
       IF( iwrit > 0) &              !echo
           WRITE(lures,"(8x,i2,4e15.5)",ERR=9999) i,e%gausv(1:4,i)

       IF( i < nlayr )CALL listen('INIG13')             !read line
       ib = 0
     END DO
     e => e%next          !point to next element
   END DO
   plstr = 0                    !do not compute plasti strains

 ELSE IF( exists('STRESS') )THEN  !read initial stresses and EFPST
   IF( iwrit > 0 )WRITE(lures,"(' ELEMENT LAY     SXX',12x,'SYY', &
&                               12x,'SXY',12x,'EFPST')",ERR=9999)
   e => heade                            !position at top of list
   ie = 1                                !initializes
   DO
     CALL listen('INIG13')          !read a card
     IF( exists('ENDSTR'))EXIT      !key-word END_STRESSES found, Exit loop
     IF( exists('NUMEL ',j))THEN    !key-word NUMEL found
       numel = INT(param(j))        !read associated element number
     ELSE
       numel = INT(param(1))        !first parameter is element number
     END IF
     IF(iwrit > 0)WRITE(lures,"(i8)",ERR=9999) numel

     i1 = ie        !remember first element
     DO
       IF( e%numel == numel )EXIT          ! element found
       ie = ie+1                           ! next element
       IF( ie == i1 )THEN                  ! if loop completed => ERROR
         WRITE(lures,"(' element label not found ',i6)",ERR=9999) numel
         CALL runend('INIG13: Element label not found !!.')
       END IF
       IF( ie > nelem )THEN                !last element passed
         ie = 1                            !go to top of the list
         e => heade
       ELSE
         e => e%next                       !go to next element in the list
       END IF
     END DO
     ie = ie + 1                         !element order
     isec = e%matno                     !material type
     IF( isec /= osec )THEN
       sec => psecs(isec)%p
       nlayr = sec%iprop(1)                     !number of layers
       nvare = sec%iprop(2)                     !number of components
       osec = isec
     END IF
     IF( .NOT.ASSOCIATED(e%gausv)) ALLOCATE ( e%gausv(nvare,nlayr) )
     DO i=1,nlayr                        !for each layer
       CALL listen('INIG13')             !read line
       IF( exists('LAYER ',j))THEN       !key-word LAYER found
         il = INT(param(j))           !read associated layer number
       ELSE
         il = INT(param(1))           !layer number is first parameter
       END IF
       IF( i /= il )   &    !layers must be in order
               CALL runend('INIG13: Invalid layer number.      ')
       IF(exists('SXX   ',j))THEN      !key-word SXX found
         e%gausv(1,i) = param(j)       !read associated value
       ELSE
         e%gausv(1,i) = param(2)       !SXX is 2th parameter
       END IF
       IF(exists('SYY   ',j)) THEN     !key-word SYY found
         e%gausv(2,i) = param(j)       !read associated value
       ELSE
         e%gausv(2,i) = param(3)       !SXX is 3rd parameter
       END IF
       IF(exists('SXY   ',j)) THEN     !key-word SXY found
         e%gausv(3,i) = param(j)       !read associated value
       ELSE
         e%gausv(3,i) = param(4)       !SXY is 4th parameter
       END IF
       IF(exists('EFPST ',j)) THEN     !key-word EFPST found
         e%gausv(4,i) = param(j)       !read associated value
       ELSE
         e%gausv(4,i) = param(5)       !EFPST is 5th parameter
       END IF
       IF( iwrit > 0)       &          !echo
           WRITE(lures,"(8x,i2,4e15.5)",ERR=9999) i,e%gausv(1:4,i)

     END DO
     e => e%next                 !point to next element
   END DO
   plstr = 1
 ELSE
    backs = .TRUE.               !no internal variables to read
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE inig13
