subroutine dlinear_solver(n, m, A, B, X, inode, verbose)

  use myconstants
  implicit none

  integer, intent(in) :: n, m, inode
  real(dp),intent(in) :: A(n,n)
  real(dp),intent(in) :: B(n,m)
  real(dp),intent(out) :: X(n,m)
  logical   :: verbose
  ! Temporary local variables
  real(dp)  :: Af(n,n), r(n), c(n), ferr(m), berr(m), work(4*n)
  integer   :: ipiv(n), iwork(n)
  character :: equed, fact, trans
  integer   :: lda, ldaf, errinfo
  real(dp)  :: rcond
  
  fact = 'E'
  trans = 'N'
  ! syntax: call dgesvx( fact, trans, n, nrhs, a, lda, 
  !           af, ldaf, ipiv, equed, r, c,
  !           b, ldb, x, ldx, rcond, ferr,
  !           berr, work, iwork, info )
  if(verbose) write (*,'("Node",i6," call dgesvx to solve A*X = B")') inode
  call dgesvx(fact, trans, n, m, A, n, &
    Af, n, ipiv, equed, r, c, &
    B, n, X, n, rcond, ferr, &
    berr, work, iwork, errinfo)
  if(errinfo .eq. 0) then
     if(verbose) write(*,'(" dgesvx is successfully excecuted. ")')
  elseif (errinfo .lt. 0) then
     write(*,'(" Node",i6," The ", i2, "-th parameter of dgesvx has an illegal value.")') &
       inode, -errinfo
  elseif (errinfo .gt. 0 .and. errinfo .le. n) then 
     write(*,'(" Node",i6," The U matrix is singular. LU factorization cannot be", &
       "completed.")') inode
  else ! errinfo = n+1
     write(*,'(" Node",i6," THe U matrix is nonsingular, but rcond is less than machine", &
       "precision.")') inode
     write(*,'(" rcond is the reciprocal of condition number: = ", f12.6)') &
       rcond
  endif

  return
end subroutine dlinear_solver
