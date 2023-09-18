 SUBROUTINE comm19(indic,nelem,nreqs,narch,elsnam,elset)

 ! parameters of present set

 IMPLICIT NONE
 CHARACTER(len=*) :: elsnam
 INTEGER (kind=4) :: indic,nelem,nreqs,narch
 TYPE (ele19_set) :: elset

 IF(indic == 0) THEN
    elset%sname = elsnam         !set name
    elset%nelem = nelem          !number of elements
    elset%nreqs = nreqs          !number of required elements for output
    elset%narch = narch          !file for output
  ELSE
    elsnam = elset%sname
    nelem  = elset%nelem
    nreqs  = elset%nreqs
    narch  = elset%narch
  END IF

 END SUBROUTINE comm19
