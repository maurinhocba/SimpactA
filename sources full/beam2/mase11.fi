 SUBROUTINE mase11(nstre,nreqs,nelem,narch,ntype,ngaus,ngrqs,sname,heade)
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: nstre,nreqs,nelem,ntype,ngaus
 INTEGER (kind=4), INTENT(IN OUT) :: narch
 INTEGER (kind=4), POINTER :: ngrqs(:)
 CHARACTER(len=*):: sname
 TYPE(ele11), POINTER :: heade

 INTEGER (kind=4) ielem
 TYPE(ele11), POINTER :: e

 WRITE(17,ERR=9999) 11,nelem,sname
 WRITE(17,ERR=9999)  2,nstre+4,ngaus,ntype
 e => heade
 DO ielem=1,nelem
   WRITE(17,ERR=9999) e%matno,e%lnods(2:3)
   e => e%next
 END DO

 IF(nreqs > 0) THEN
    CALL openfi(nfile=21,narch=narch)
    WRITE(narch,ERR=9999) nreqs,nstre+3,ngrqs
    SELECT CASE (ntype)
    CASE(1)
      WRITE(narch,ERR=9999) 11,sname,'N_1  ','M_1  ','EP_B ','EP_T ','Lambd'
    CASE(2:3)
      WRITE(narch,ERR=9999) 11,sname,'N_1  ','N_2  ','M_1  ','M_2  ','Q_1  ','EP_B ', &
                                     'EP_T ','Lambd'
    END SELECT
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE mase11
