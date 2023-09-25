 SUBROUTINE resta8(nelem, nnode, nreqs, axesc, heade, taile, ngrqs, &
                   posgp, shape, deriv, weigh, ngaus)
 !******************************************************************
 !
 !*** restores an element set at restart
 !
 !******************************************************************
 IMPLICIT NONE

 INTEGER (kind=4) :: nelem, nnode, nreqs, axesc, ngaus
 INTEGER (kind=4), POINTER :: ngrqs(:)
 REAL (kind=8) ::  posgp(:), shape(:,:), deriv(:,:),weigh(:)
 TYPE (ele08), POINTER :: heade,taile

 INTEGER (kind=4) :: i,ielem
 TYPE (ele08), POINTER :: e

 READ (51) posgp, shape, deriv, weigh
 DO ielem=1,nelem
   ALLOCATE (e)      !get memory for an element
   ALLOCATE (e%lnods(nnode))  !reserve memory for data
   ALLOCATE (e%llbd(3,3,MAX(1,axesc)),e%jac(ngaus),e%stra0(6,ngaus),  &  !reserve memory for data
             e%stran(6,ngaus),e%stres(6,ngaus),e%epdef(ngaus),e%sedef(ngaus))
   READ (51) e%lnods(1:nnode),e%matno,e%numel
   READ (51) e%llbd,e%jac,e%stra0,e%stran,e%stres,e%epdef,e%sedef
   CALL add_ele08e (e, heade, taile)           !add element to the list
 END DO

 IF ( nreqs > 0 ) READ (51) (ngrqs(i), i=1,nreqs)
 RETURN
 END SUBROUTINE resta8
