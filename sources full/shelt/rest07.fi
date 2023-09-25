 SUBROUTINE rest07 (nelem, nreqs, heade, taile, ngrqs, stype, ngaus, ngamm, nnode, nstre)
 !******************************************************************
 !
 !***  restore an element set at restart for QUAD shell
 !
 !******************************************************************
 IMPLICIT NONE
 ! dummy arguments
 INTEGER (kind=4), INTENT(IN) :: nreqs, & !Number of elements for requested output
                                 stype, & !formulation type
                                 nstre, & !number of stress values
                                 ngaus, & !number of integration points
                                 ngamm, & !number of shear values
                                 nnode, & !number of nodes in lnods
                                 nelem    !number of elements in the set
 INTEGER (kind=4), POINTER  :: ngrqs(:)   !elements for requested output
 TYPE (ele07), POINTER  :: heade,taile    !pointers to first and last element

 INTEGER (kind=4) :: i,nlayr,ielem,mtype,isec,osec
 TYPE (ele07), POINTER :: e
 TYPE (section), POINTER :: sec      !auxiliar pointer
 LOGICAL :: plast

 osec = -1
 DO ielem=1,nelem    !for each element

   CALL new_ele07e (e)           !reserve memory for data
   ALLOCATE( e%dvolu(ngaus),e%lnods(nnode),e%jacin(2,2,ngaus), &
             e%stra0(6,ngaus),e%gamm0(ngamm),e%strsg(nstre,ngaus),     &
             e%qbar(ngamm),e%ambda(2,ngaus))
   IF(ngaus == 3)THEN
     ALLOCATE(e%cartd(nnode,2,ngaus))
   ELSE
     ALLOCATE(e%cartd(3,2,ngaus))
   END IF
   READ (51) e%numel,e%matno,e%lnods
   READ (51) e%angle,e%dvolu,e%cartd,e%jacin,e%stra0,e%gamm0,e%strsg, &
             e%qbar,e%ambda
   isec = e%matno
   IF( stype > 1 )THEN
     ALLOCATE( e%beta(nbeta) )
     READ (51) e%beta
   END IF
   IF( nnode == 9 )THEN
     ALLOCATE( e%cd(4,2,3) )
     READ (51) e%cd
   END IF
   IF( isec /= osec )THEN
     sec => psecs(isec)%p
     nlayr = sec%iprop(1)     !number of layers
     !nucom = sec%iprop(2)     !number of material components
     mtype = sec%secdef(3)    !section constitutive type
     plast = sec%mtbas%matdef(3) > 1  .AND. mtype <= 2   !plasticiy
     osec = isec
   END IF
   IF( plast )THEN
     SELECT CASE (mtype)              !according to material model
       !  Plasticity is restricted to models 0-2
     CASE (0)     !layered large strain plasticity
       ALLOCATE( e%ehist(5,ngaus), e%strap(6*nlayr,ngaus))
       ALLOCATE( e%stres(5,nlayr,ngaus) )
       READ (51) e%ehist,e%strap
       READ (51) e%stres
     CASE (1,2)  !solid isotropic
       ALLOCATE( e%ehist(5,ngaus), e%strap(nstre,ngaus) )
       READ (51) e%ehist,e%strap
     END SELECT
   ELSE
     NULLIFY( e%ehist,e%strap )
   END IF

   CALL add_ele07e (e, heade, taile)           !add element to the list

 END DO
 IF ( nreqs > 0 ) READ (51) (ngrqs(i), i=1,nreqs)  !read elements for output

 RETURN
 END SUBROUTINE rest07
