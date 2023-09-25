 MODULE nsld_db
   IMPLICIT NONE
   SAVE

   ! Derived type for a list containing dependent nodes

   ! list
   TYPE nsld_nod
     INTEGER (kind=4) :: slave,master
     LOGICAL :: autom,nosld,norot
     TYPE (nsld_nod), POINTER :: next
   END TYPE nsld_nod
   INTEGER (kind=4) :: nnsld=0  !number of master-slave pairs
   INTEGER (kind=4), POINTER :: nsldat(:,:)
   LOGICAL, POINTER :: nsldfg(:,:)

 CONTAINS
   SUBROUTINE ini_sld (head, tail)
     !initialize a dependent nodes list

     !Dummy arguments
     TYPE (nsld_nod), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_sld

   SUBROUTINE add_sld (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (nsld_nod), POINTER :: new, head, tail

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
   END SUBROUTINE add_sld

   SUBROUTINE store_sld (head)
     !transfer data from list into array
     !and release memory of list

     !Dummy argument
     TYPE (nsld_nod), POINTER :: head

     !Local variables
     TYPE (nsld_nod), POINTER :: ptr,prev
     INTEGER :: n

     ptr => head

     n = 0                         !initializes array counter
     DO
       IF (.NOT.ASSOCIATED(ptr)) EXIT
       n = n + 1                   !increase counter
       nsldat(1,n) = ptr%slave     !transfer data
       nsldat(2,n) = ptr%master
       nsldfg(1,n) = ptr%autom
       nsldfg(2,n) = ptr%nosld
       nsldfg(3,n) = ptr%norot
       prev => ptr                 !remember pointer
       ptr => ptr%next             !point to next
       DEALLOCATE (prev)           !release data already transfered
     END DO

     RETURN
   END SUBROUTINE store_sld

 END MODULE nsld_db
