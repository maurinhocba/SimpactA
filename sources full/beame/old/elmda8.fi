 SUBROUTINE elmda8(ndime,nelem,nnode,axesc,coord,eule0,euler,  &
                   heade,taile,iwrit,ngaus)
 !******************************************************************
 !
 !*** READ element DATA for beam elements
 !
 !******************************************************************
 IMPLICIT NONE
 INTEGER (kind=4) :: ndime,nelem,nnode,axesc,iwrit,ngaus
 REAL    (kind=8) :: coord(:,:),eule0(:,:),euler(:,:)
 TYPE(ele08), POINTER :: heade,taile
 ! local variables
 INTEGER (kind=4) n,isec,i,nr,mtype,osec,ip,na
 LOGICAL found
 INTEGER (kind=4) chnode
 INTEGER (kind=4) :: secs = 42
 TYPE(section), POINTER :: sec
 TYPE(ele08), POINTER :: elem

 ! *** READ the element nodal connections, and the property numbers.

 IF(iwrit == 1) WRITE (lures,900,ERR=9999)
 osec = -1                        !initializes
 na = MAX(ABS(axesc),1)

 elem => heade                    !point to first

 DO
   CALL listen('ELMD08')          !read a card
   IF (exists('ENDELE')) EXIT     !if key word END_ELEMENT read, exit loop

   nelem = nelem + 1         !increase element counter
   ALLOCATE (elem)           !reserve memory for data
   ALLOCATE (elem%lnods(nnode))  !reserve memory for data
   ALLOCATE (elem%llbd(3,3,na),elem%jac(ngaus),elem%stran(6,ngaus),  &  !reserve memory for data
             elem%stres(6,ngaus),elem%epdef(ngaus),elem%sedef(ngaus))
   elem%llbd = 0d0 ; elem%jac  = 0d0 ; elem%stran= 0d0
   elem%stres= 0d0 ; elem%epdef= 0d0 ; elem%sedef= 0d0
   elem%numel = INT (param(1))    !element label
   isec = INT (param(2))         !associated material
   IF( isec /= osec )THEN
     CALL sect_search (isec,found,sec,ip)
     IF(.NOT.found)CALL runend ('BEAME:  Invalid Section Number used')
     mtype = sec%secty       !check section type for element
     IF( mtype /= secs ) CALL runend('BEAME : Invalid sect. type for BEAME')
     !plast = sec%rprop(17) > 0d0  !plasticiy
     osec = isec
   END IF
   elem%matno = ip      !store section position in array
   elem%lnods(1:nnode) = INT (param(3:2+nnode))  !connectivities

   ! read reference nodes if AXES_CODE < 0
   nr = axesc                             !axes code
   SELECT CASE (nr)                       !according to code
   CASE (:-1)                             !read auxiliar nodes if necessary
     DO n=1,-nr
       i = INT(param(n+2+nnode))
       elem%llbd(n,1,1) = chnode(i)
     END DO
     nr = -nr
   CASE (0)
     ! nothing  !elem%llbd(1,1,1)=  0
   CASE (1:)
     elem%llbd(1,1,1)= -2d0               !flag
     nr = 0
   END SELECT
   IF(iwrit == 1) WRITE(lures,901,ERR=9999) elem%numel,isec,elem%lnods(1:nnode),    &
                                            INT(elem%llbd(1:nr,1,1))

   !     original nodes are changed
   DO i=1,nnode
     elem%lnods(i) = chnode(elem%lnods(i))
   END DO
   !     interpolates mid-side nodes
   IF(nnode == 3)CALL nodxy8(ndime,nnode,elem%lnods,coord,eule0,euler)
   CALL add_ele08e (elem, heade, taile) !add element to list
 END DO

 RETURN
 900 FORMAT(//5X,' Element    Mat       Node Numbers ')
 901 FORMAT(6X,I5,I9,5X,6I5)
 9999 CALL runen2('')
 END SUBROUTINE elmda8
