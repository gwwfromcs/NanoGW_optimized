subroutine dlinear_solver(n, m, A, B, X)

  use myconstants
  implicit none

  integer, intent(in) :: n, m 
  real(dp),intent(in) :: A(n,n)
  real(dp),intent(in) :: B(n,m)
  real(dp),intent(out) :: X(n,m)
  !temporary variables
  real(dp) :: Af(n,n), r(n), c(n), ferr(m), berr(m), work(4*n)
  integer :: ipiv(n), iwork(n)
  character :: equed, fact, trans
  integer :: lda, ldaf, errinfo
  real(dp) :: rcond
  
  fact = 'E'
  trans = 'N'
  ! syntax: call dgesvx( fact, trans, n, nrhs, a, lda, 
  !           af, ldaf, ipiv, equed, r, c,
  !           b, ldb, x, ldx, rcond, ferr,
  !           berr, work, iwork, info )
  write (*,'(" call dgesvx to solve A*X = B")')
  call dgesvx(fact, trans, n, m, A, n, &
    Af, n, ipiv, equed, r, c, &
    B, n, X, n, rcond, ferr, &
    berr, work, iwork, errinfo)
  if(errinfo .eq. 0) then
     write(*,'(" dgesvx is successfully excecuted. ")')
  elseif (errinfo .lt. 0) then
     write(*,'(" The ", i2, "-th parameter of dgesvx has an illegal value.")') &
       -errinfo
  elseif (errinfo .gt. 0 .and. errinfo .le. n) then 
     write(*,'(" The U matrix is singular. LU factorization cannot be", &
       "completed.")')
  else ! errinfo = n+1
     write(*,'(" THe U matrix is nonsingular, but rcond is less than machine", &
       "precision.")')
     write(*,'(" rcond is the reciprocal of condition number: = ", f12.6)') &
       rcond
  endif

  return
end subroutine dlinear_solver
