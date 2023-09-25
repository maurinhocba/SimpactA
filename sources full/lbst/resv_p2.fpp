   !check maximum and miminum thickness ratio
   ! uses
   !
   !   INTEGER (kind=4), INTENT(OUT) :: ierr    !flag for error detection
   !   REAL (kind=8) min_tr,    & !Minimum Allowable Thickness ratio
   !                 max_tr       !Maximum Allowable Thickness ratio
   !
   !   CHARACTER(len=1 ) letra  !for console input
   !
   !   TYPE (ele14), POINTER :: e    !pointer to an element data
   !   TYPE (section), POINTER :: sec    !pointer to a section data

   IF( e%lb < min_tr .OR. e%lb > max_tr)THEN      !check
     WRITE(*,"(' thickness ratio limit exceeded, TH_RATIO: ',f7.3, &
           &  /,' Do you want to continue [Y/N]')")e%lb
     DO                            !pause until answer
       READ(*,'(a1)')letra         !read the answer
       SELECT CASE (letra)         !process answer
       CASE ( 'y', 'Y')         !continue
         IF( e%lb < min_tr )THEN    !minimum exceeded
           WRITE(*,"(' New value for minimum thickness ratio')")
           READ(*,*) min_tr        !read value
           sec%rprop(2) = min_tr   !Store Minimum Thickness ratio
         ELSE                      !maximum exceeded
           WRITE(*,"(' New value for maximum thickness ratio')")
           READ(*,*) max_tr        !read value
           sec%rprop(3) = max_tr   !Store Maximum Thickness ratio
         END IF
       CASE ( 'n', 'N')         !stop process
         ierr = 1                  !modify flag
       CASE DEFAULT
         CYCLE                     !answer not valid, ask again
       END SELECT
       EXIT
     END DO
   END IF
