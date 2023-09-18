 SUBROUTINE cupdl1( )
 !
 !     RENUMBER nodes in surface contact Data_Base
 !
 IMPLICIT NONE
 ! routine arguments

 ! local variables
 INTEGER (kind=4), PARAMETER :: nnseg = 2
 LOGICAL :: found,flag
 INTEGER (kind=4) :: isurf
 TYPE (surf1_db), POINTER :: surf
 INTEGER (kind=4), POINTER :: lcseg(:,:)

 INTERFACE
   INCLUDE 'elemnt.h'
 END INTERFACE

 ! Updates surface connectivities

 surf => shead                               ! initializes
 DO isurf=1,nsurf                            ! for each surface
   CALL new_srf(surfa)
   CALL elemnt ('SURFAC', name=surf%sname, flag2=found)  ! get surface from element database
   IF( found )THEN
     ALLOCATE( lcseg(nnseg,2*surfa%nelem) )                 ! temporary auxiliar space
     CALL store_segs (surfa%head,lcseg,nnseg,surfa%nelem) !surface ==> lcseg
     flag = surfa%nelem /= surf%nsegm
     surf%nsegm = surfa%nelem
   ELSE
     CALL upndc1(surf,oldlb,flag)     ! update nodal numeration
   END IF
   CALL chsur1(surf,found,flag,lcseg)
   CALL dalloc_srf(surfa)
   IF(ASSOCIATED(lcseg)) DEALLOCATE ( lcseg )     ! free auxiliar space
   surf => surf%next                         ! next surface
 END DO

 RETURN
 END SUBROUTINE cupdl1
