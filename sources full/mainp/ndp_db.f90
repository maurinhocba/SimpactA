 MODULE ndp_db
   IMPLICIT NONE
   SAVE

   ! Derived type for a list containing dependent nodes

   ! list
   TYPE ndp_nod
     INTEGER (kind=4) :: slave,master
     TYPE (ndp_nod), POINTER :: next
   END TYPE ndp_nod
   INTEGER (kind=4) :: nndp=0  !number of master-slave pairs
   INTEGER (kind=4), POINTER :: ndpdat(:,:)

 CONTAINS
   SUBROUTINE ini_ndp (head, tail)
     !initialize a dependent nodes list

     !Dummy arguments
     TYPE (ndp_nod), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_ndp

   SUBROUTINE add_ndp (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ndp_nod), POINTER :: new, head, tail

     !Check if a list is empty
     IF (.NOT. ASSOCIATED (head)) THEN
       !list is empty, start it
       head => new
       tail => new
       NULLIFY (tail%next)

     ELSE
       !add data to the list
       tail%next => new
       NULLIFY (new%next)
       tail => new

     ENDIF
   END SUBROUTINE add_ndp

   SUBROUTINE store_ndp (head)
     !transfer data from list into array
     !and release memory of list

     !Dummy argument
     TYPE (ndp_nod), POINTER :: head

     !Local variables
     TYPE (ndp_nod), POINTER :: ptr,prev
     INTEGER :: n

     ptr => head

     n = 0                         !initializes array counter
     DO
       IF (.NOT.ASSOCIATED(ptr)) EXIT
       n = n + 1                   !increase counter
       ndpdat(1,n) = ptr%slave     !transfer data
       ndpdat(2,n) = ptr%master
       prev => ptr                 !remember pointer
       ptr => ptr%next             !point to next
       DEALLOCATE (prev)           !release data already transfered
     END DO

     RETURN
   END SUBROUTINE store_ndp

 END MODULE ndp_db
