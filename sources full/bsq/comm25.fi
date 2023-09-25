 SUBROUTINE comm25(indic,nelem,nreqs,narch,elsnam,elset,logst)

 ! parameters of present set

 IMPLICIT NONE
 CHARACTER(len=*) :: elsnam
 LOGICAL  :: logst
 INTEGER (kind=4) :: indic,nelem,nreqs,narch
 TYPE (ele25_set) :: elset

 IF(indic == 0) THEN
    elset%sname = elsnam         !set name
    elset%nelem = nelem          !number of elements
    elset%nreqs = nreqs          !number of required elements for output
    elset%narch = narch          !file for output
    elset%logst = logst          !use logarithmic strain
  ELSE
    elsnam = elset%sname
    nelem  = elset%nelem
    nreqs  = elset%nreqs
    narch  = elset%narch
    logst  = elset%logst
  END IF

 END SUBROUTINE comm25
