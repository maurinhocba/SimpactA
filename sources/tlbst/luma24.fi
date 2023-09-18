 SUBROUTINE luma24 (heade,ymass,sumat)

 ! compute lumped mass vector for element BST-CST (TLF) shell element

 IMPLICIT NONE
 !  dummy arguments

 REAL (kind=8) :: ymass(:,:),      & !nodal mass in each direction
                  sumat              !total mass of elements
 TYPE (ele24), POINTER  :: heade     !pointer to first element

 ! local variables
 INTEGER (kind=4) isec,lnode,idofn,inode,osec
 REAL (kind=8) :: diagm,voli,rhoel,thick
 TYPE (ele24), POINTER :: elem
 TYPE (section), POINTER :: sec

 osec = -1
 elem => heade    !point to first element
 DO
   IF (.NOT.ASSOCIATED (elem) ) EXIT     !if all elements processed, exit
   isec = elem%matno                    !associated material
   IF( isec /= osec )THEN
     sec => psecs(isec)%p
     !IF( sec%secty == 12 )THEN
       thick = sec%rprop(1)                !shell thickness
       rhoel = sec%mtbas%prope(5)*thick    !surface density
     !ELSE
     !  rhoel = sec%rprop(5)                !surface density
     !END IF
     osec = isec
   END IF
   voli=elem%area1*rhoel           !element mass
   diagm = voli/3.d0                     !nodal mass
   !
   ! ---generates lumped mass matrix proportional to diagonal
   !
   DO inode=1,3                          !for each element node
     lnode=elem%lnods(inode)             !node number
     DO idofn=1,3                        !translational DOFs only
       ymass(idofn,lnode)=ymass(idofn,lnode)+diagm
     END DO
   END DO
   sumat=sumat+voli                      !add to total mass
   elem => elem%next                     !point to next element
 END DO
 RETURN
 END SUBROUTINE luma24
