 SUBROUTINE rest19 (nelem, heade, taile)
 !******************************************************************
 !
 !***  restore an element set at restart for 2-D triangular element (TLF)
 !
 !******************************************************************
 IMPLICIT NONE

   !--- Dummy arguments
   INTEGER(kind=4),INTENT(IN):: nelem    !number of elements in the set
   TYPE(ele19),POINTER:: heade,taile  !pointers to first and last element
   !--- Variables variable
   INTEGER(kind=4):: i,ielem,nvare
   TYPE(ele19),POINTER:: e

   DO ielem=1,nelem    !for each element

     CALL new_ele19e(e)      !get memory for an element
     READ(51) e%numel,e%matno,e%lnods,nvare
     READ(51) e%area1,e%angle,e%cd,e%stint,e%beta,e%repla
     IF( nvare > 0  )THEN
       ALLOCATE(e%gausv(nvare))          !get space for internal vars
       READ(51) (e%gausv(i),i=1,nvare)   !read internal variables
     END IF

     CALL add_ele19e(e, heade, taile)           !add element to the list

   END DO

 RETURN
 END SUBROUTINE rest19
