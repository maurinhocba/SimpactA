
MODULE GenUseNR
  
  CONTAINS !=================================================================
  ! Routines from:
  !     Numerical Recipes in Fortran 90
  !     B2 Solution of Linear Algebraic Equations   - p. 1014
  !     B3 Interpolation and Extrapolation          - p. 1043
  
  SUBROUTINE tridag_ser(a,b,c,r,u)
    USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c,r
    REAL(SP), DIMENSION(:), INTENT(OUT) :: u
    ! Solves for a vector u of size N the tridiagonal linear set given by equation (2.4.1) using a
    ! serial algorithm. Input vectors b (diagonal elements) and r (right-hand sides) have size N,
    ! while a and c (off-diagonal elements) are size N - 1.
    REAL(SP), DIMENSION(size(b)) :: gam ! One vector of workspace, gam is needed.
    INTEGER(I4B) :: n,j
    REAL(SP) :: bet
    n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_ser')
    bet=b(1)
    if (bet == 0.0) call nrerror('tridag_ser: Error at code stage 1')
    ! If this happens then you should rewrite your equations as a set of order N - 1, with u2
    ! trivially eliminated.
    u(1)=r(1)/bet
    do j=2,n ! Decomposition and forward substitution.
    gam(j)=c(j-1)/bet
    bet=b(j)-a(j-1)*gam(j)
    if (bet == 0.0) & ! Algorithm fails; see below routine in Vol. 1.
    call nrerror('tridag_ser: Error at code stage 2')
    u(j)=(r(j)-a(j-1)*u(j-1))/bet
    end do
    do j=n-1,1,-1 ! Backsubstitution.
    u(j)=u(j)-gam(j+1)*u(j+1)
    end do
  END SUBROUTINE tridag_ser
  !--------------------------------------------------------------------------
  RECURSIVE SUBROUTINE tridag_par(a,b,c,r,u)
    USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
    !USE nr, ONLY : tridag_ser - COMMENTED OUT BY MAURO MAZA <- OTHERWISE IT WOULD NOT LINK
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c,r
    REAL(SP), DIMENSION(:), INTENT(OUT) :: u
    ! Solves for a vector u of size N the tridiagonal linear set given by equation (2.4.1) using a
    ! parallel algorithm. Input vectors b (diagonal elements) and r (right-hand sides) have size
    ! N, while a and c (off-diagonal elements) are size N - 1.
    INTEGER(I4B), PARAMETER :: NPAR_TRIDAG=4 ! Determines when serial algorithm is invoked.
    INTEGER(I4B) :: n,n2,nm,nx
    REAL(SP), DIMENSION(size(b)/2) :: y,q,piva
    REAL(SP), DIMENSION(size(b)/2-1) :: x,z
    REAL(SP), DIMENSION(size(a)/2) :: pivc
    n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_par')
    if (n < NPAR_TRIDAG) then
    call tridag_ser(a,b,c,r,u)
    else
    if (maxval(abs(b(1:n))) == 0.0) & ! Algorithm fails; see below routine in Vol. 1.
    call nrerror('tridag_par: possible singular matrix')
    n2=size(y)
    nm=size(pivc)
    nx=size(x)
    piva = a(1:n-1:2)/b(1:n-1:2) ! Zero the odd a’s and even c’s, giving x, y, z, q.
    pivc = c(2:n-1:2)/b(3:n:2)
    y(1:nm) = b(2:n-1:2)-piva(1:nm)*c(1:n-2:2)-pivc*a(2:n-1:2)
    q(1:nm) = r(2:n-1:2)-piva(1:nm)*r(1:n-2:2)-pivc*r(3:n:2)
    if (nm < n2) then
    y(n2) = b(n)-piva(n2)*c(n-1)
    q(n2) = r(n)-piva(n2)*r(n-1)
    end if
    x = -piva(2:n2)*a(2:n-2:2)
    z = -pivc(1:nx)*c(3:n-1:2)
    call tridag_par(x,y,z,q,u(2:n:2)) ! Recurse and get even u’s.
    u(1) = (r(1)-c(1)*u(2))/b(1) ! Substitute and get odd u’s.
    u(3:n-1:2) = (r(3:n-1:2)-a(2:n-2:2)*u(2:n-2:2) &
    -c(3:n-1:2)*u(4:n:2))/b(3:n-1:2)
    if (nm == n2) u(n)=(r(n)-a(n-1)*u(n-1))/b(n)
    end if
  END SUBROUTINE tridag_par
  !--------------------------------------------------------------------------
  SUBROUTINE spline(x, y, yp1, ypn, y2)
    USE nrtype; USE nrutil, ONLY : assert_eq
    !USE nr, ONLY : tridag - COMMENTED OUT BY MAURO MAZA <- OTHERWISE IT WOULD NOT LINK
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN)  :: x,y
    REAL(SP),               INTENT(IN)  :: yp1,ypn
    REAL(SP), DIMENSION(:), INTENT(OUT) :: y2
    ! Given arrays x and y of length N containing a tabulated function, i.e., yi = f(xi), with x1 <
    ! x2 < ... < xN, and given values yp1 and ypn for the first derivative of the interpolating
    ! function at points 1 and N, respectively, this routine returns an array y2 of length N
    ! that contains the second derivatives of the interpolating function at the tabulated points
    ! xi. If yp1 and/or ypn are equal to 1 × 1030 or larger, the routine is signaled to set the
    ! corresponding boundary condition for a natural spline, with zero second derivative on that
    ! boundary.
    INTEGER(I4B) :: n
    REAL(SP), DIMENSION(size(x)) :: a,b,c,r
    n=assert_eq(size(x),size(y),size(y2),'spline')
    c(1:n-1)=x(2:n)-x(1:n-1) ! Set up the tridiagonal equations.
    r(1:n-1)=6.0_sp*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1)=r(2:n-1)-r(1:n-2)
    a(2:n-1)=c(1:n-2)
    b(2:n-1)=2.0_sp*(c(2:n-1)+a(2:n-1))
    b(1)=1.0
    b(n)=1.0
    if (yp1 > 0.99e30_sp) then ! The lower boundary condition is set either to be “natrural”
    c(1)=0.0
    else ! or else to have a specified first derivative.
    r(1)=(3.0_sp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    end if ! The upper boundary condition is set either to be
    if (ypn > 0.99e30_sp) then ! “natural”
    r(n)=0.0
    a(n)=0.0
    else ! or else to have a specified first derivative.
    r(n)=(-3.0_sp/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
    a(n)=0.5
    end if
        ! ORIGIAL CODE BELOW
        !call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
        ! NEW CODE - Mauro Maza
        call tridag_par(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
        
  ENDSUBROUTINE spline
  !--------------------------------------------------------------------------
  FUNCTION splintMM(xa, ya, y2a, x)
    ! Mauro S. Maza - 17/04/2014
    ! modified to look for klo using "hunt" subroutine
    USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
    USE nr, ONLY: locate
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
    REAL(SP), INTENT(IN) :: x
    REAL(SP) :: splintMM
    ! Given the arrays xa and ya, which tabulate a function (with the xai ’s in increasing or
    ! decreasing order), and given the array y2a, which is the output from spline above, and
    ! given a value of x, this routine returns a cubic-spline interpolated value. The arrays xa, ya
    ! and y2a are all of the same size.
        ! ORIGINAL CODE
        !INTEGER(I4B) :: khi,klo,n
        ! NEW CODE
        INTEGER(I4B)       :: khi,n
        INTEGER(I4B), SAVE :: klo
    REAL(SP) :: a,b,h
    n=assert_eq(size(xa),size(ya),size(y2a),'splint')
        ! ORIGINAL CODE
        !klo=max(min(locate(xa,x),n-1),1)
        !! We will find the right place in the table by means of locate’s bisection algorithm. This is
        !! optimal if sequential calls to this routine are at random values of x. If sequential calls are in
        !! order, and closely spaced, one would do better to store previous values of klo and khi and
        !! test if they remain appropriate on the next call.
        ! NEW CODE
        CALL hunt(xa,x,klo)
    khi=klo+1 ! klo and khi now bracket the input value of x.
    h=xa(khi)-xa(klo)
    if (h == 0.0) call nrerror('bad xa input in splint') ! The xa’s must be distinct.
    a=(xa(khi)-x)/h ! Cubic spline polynomial is now evaluated.
    b=(x-xa(klo))/h
    splintMM=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_sp
  ENDFUNCTION splintMM
  ! NOTE ON FORTRAN IMPLEMENTATION
  ! klo=max(min(locate(xa,x),n-1),1) In the Fortran 77 version of splint,
  ! there is in-line code to find the location in the table by bisection. Here
  ! we prefer an explicit call to locate, which performs the bisection. On
  ! some massively multiprocessor (MMP) machines, one might substitute a different,
  ! more parallel algorithm (see next note).
  !--------------------------------------------------------------------------
  FUNCTION locate(xx,x)
    USE nrtype
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xx
    REAL(SP), INTENT(IN) :: x
    INTEGER(I4B) :: locate
    ! Given an array xx(1:N), and given a value x, returns a value j such that x is between
    ! xx(j) and xx(j + 1). xx must be monotonic, either increasing or decreasing. j = 0 or
    ! j = N is returned to indicate that x is out of range.
    INTEGER(I4B) :: n,jl,jm,ju
    LOGICAL :: ascnd
    n=size(xx)
    ascnd = (xx(n) >= xx(1)) ! True if ascending order of table, false otherwise.
    jl=0 ! Initialize lower
    ju=n+1 ! and upper limits.
    do
    if (ju-jl <= 1) exit ! Repeat until this condition is satisfied.
    jm=(ju+jl)/2 ! Compute a midpoint,
    if (ascnd .eqv. (x >= xx(jm))) then
    jl=jm ! and replace either the lower limit
    else
    ju=jm ! or the upper limit, as appropriate.
    end if
    end do
    if (x == xx(1)) then ! Then set the output, being careful with the endpoints.
    locate=1
    else if (x == xx(n)) then
    locate=n-1
    else
    locate=jl
    end if
  ENDFUNCTION locate
  ! NOTE ON PARALLEL IMPLEMENTATION
  ! The use of bisection is perhaps a sin on a genuinely parallel machine, but
  ! (since the process takes only logarithmicallymany sequential steps) it is at
  ! most a small sin. One can imagine a “fully parallel” implementation like,
  !   k=iminloc(abs(x-xx))
  !   if ((x < xx(k)) .eqv. (xx(1) < xx(n))) then
  !     locate=k-1
  !   else
  !     locate=k
  !   end if
  ! Problem is, unless the number of physical (not logical) processors participating in
  ! the iminloc is larger than N, the length of the array, this “parallel” code turns a
  ! logN algorithm into one scaling as N, quite an unacceptable inefficiency. So we
  ! prefer to be small sinners and bisect.
  !--------------------------------------------------------------------------
  SUBROUTINE hunt(xx,x,jlo)
    USE nrtype
    IMPLICIT NONE
    INTEGER(I4B), INTENT(INOUT) :: jlo
    REAL(SP), INTENT(IN) :: x
    REAL(SP), DIMENSION(:), INTENT(IN) :: xx
    ! Given an array xx(1:N), and given a value x, returns a value jlo such that x is between
    ! xx(jlo) and xx(jlo+1). xx must be monotonic, either increasing or decreasing. jlo = 0
    ! or jlo = N is returned to indicate that x is out of range. jlo on input is taken as the
    ! initial guess for jlo on output.
    INTEGER(I4B) :: n,inc,jhi,jm
    LOGICAL :: ascnd
    n=size(xx)
    ascnd = (xx(n) >= xx(1)) ! True if ascending order of table, false otherwise.
    if (jlo <= 0 .or. jlo > n) then ! Input guess not useful. Go immediately to bisection.
      jlo=0
      jhi=n+1
    else
      inc=1 ! Set the hunting increment.
      if (x >= xx(jlo) .eqv. ascnd) then ! Hunt up:
        do
          jhi=jlo+inc
          if (jhi > n) then ! Done hunting, since off end of table.
            jhi=n+1
            exit
          else
            if (x < xx(jhi) .eqv. ascnd) exit
            jlo=jhi ! Not done hunting,
            inc=inc+inc ! so double the increment
          end if
        end do ! and try again.
      else ! Hunt down:
        jhi=jlo
        do
          jlo=jhi-inc
          if (jlo < 1) then ! Done hunting, since off end of table.
            jlo=0
            exit
          else
            if (x >= xx(jlo) .eqv. ascnd) exit
            jhi=jlo ! Not done hunting,
            inc=inc+inc !so double the increment
          end if
        end do ! and try again.
      end if
    end if ! Done hunting, value bracketed.
    do ! Hunt is done, so begin the final bisection phase:
      if (jhi-jlo <= 1) then
        if (x == xx(n)) jlo=n-1
        if (x == xx(1)) jlo=1
        exit
      else
        jm=(jhi+jlo)/2
        if (x >= xx(jm) .eqv. ascnd) then
          jlo=jm
        else
          jhi=jm
        end if
      end if
    end do
  END SUBROUTINE hunt
  
ENDMODULE GenUseNR

!============================================================================

!============================================================================

!============================================================================
