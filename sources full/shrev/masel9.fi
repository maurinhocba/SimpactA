 SUBROUTINE masel9(ntype,nnode,nstre,ngaus,nreqs,nelem, &
                   narch,ngrqs,sname,heade)
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: ntype,nnode,nstre,ngaus,nreqs,nelem
 INTEGER (kind=4), INTENT(IN OUT) :: narch
 INTEGER (kind=4), POINTER :: ngrqs(:)
 CHARACTER (len=*) :: sname
 TYPE(ele09), POINTER :: heade

 INTEGER (kind=4) ielem
 TYPE(ele09), POINTER :: e

 WRITE(17,ERR=9) 9,nelem,sname
 WRITE(17,ERR=9) nnode,nstre+3,ngaus,ntype
 e => heade
 DO ielem=1,nelem
   WRITE(17,ERR=9) e%matno,e%lnods(1:nnode)
   e => e%next
 END DO

 IF(nreqs > 0) THEN
    CALL openfi(nfile=21,narch=narch)
    WRITE(narch,ERR=9) nreqs,nstre+3,ngrqs
    SELECT CASE (ntype)
    CASE(1)
      WRITE(narch,ERR=9) 11,sname,'N_1  ','M_1  ','Q_1  ','EP_B ','EP_T ','Lambd'
    CASE(2:3)
      WRITE(narch,ERR=9) 11,sname,'N_1  ','N_2  ','M_1  ','M_2  ','Q_1  ','Q_2  ', &
            'EP_B ', 'EP_T ','Lambd'
    END SELECT
 END IF

 RETURN
 9 CALL runen2('')
 END SUBROUTINE masel9
