subroutine printmatrix ( mat, n, m, iounit)
 
 real, intent(in) :: mat(n,m)
 integer, intent(in) :: n, m, iounit
 integer, parameter :: w = 9, w2 = 4
 integer :: i, j
 
 write ( iounit, '( A, i12, A, i12, A )') "size: (", n, " * ", m, " )"
 if ( n <= w .and. m <= w ) then
   ! print all the matrix 
   do i = 1, n
      do j = 1, m
        write ( iounit, '( 3x, e12.5)', advance='no' ) mat(i,j)
      enddo
      write ( iounit, * ) 
   enddo
 endif

 if ( n <= w .and. m > w ) then
   ! print part of the matrix
   do i = 1, n
      do j = 1, w2
        write ( iounit, '( 3x, e12.5)', advance='no' ) mat(i,j)
      enddo
      write ( iounit, '( "  ...  " )', advance='no' )
      do j = m-w2+1, m
        write ( iounit, '( 3x, e12.5)', advance='no' ) mat(i,j)
      enddo
      write ( iounit, * )
   enddo
 endif

 if ( n > w .and. m <= w ) then
   ! print part of the matrix
   do i = 1, w2
      do j = 1, m
        write ( iounit, '( 3x, e12.5)', advance='no' ) mat(i,j)
      enddo
      write ( iounit, * )
   enddo
   write ( iounit, '( 6x, "  ... ... ...  " )' )
   do i = n-w2+1, n
      do j = 1, m
        write ( iounit, '( 3x, e12.5)', advance='no' ) mat(i,j)
      enddo
      write ( iounit, * )
   enddo
 endif

 if ( n > w .and. m > w ) then
   ! print part of the matrix
   do i = 1, w2
      do j = 1, w2
         write ( iounit, '( 3x, e12.5)', advance='no' ) mat(i,j)
      enddo
      write ( iounit, '( "  ...  " )', advance='no' )
      do j = m-w2+1, m
         write ( iounit, '( 3x, e12.5)', advance='no' ) mat(i,j)
      enddo
      write ( iounit, * )
   enddo
   write ( iounit, '( 6x, "  ... ... ...  " )' )
   do i = n-w2+1, n
      do j = 1, w2
         write ( iounit, '( 3x, e12.5)', advance='no' ) mat(i,j)
      enddo
      write ( iounit, '( "  ...  " )', advance='no' )
      do j = m-w2+1, m
         write ( iounit, '( 3x, e12.5)', advance='no' ) mat(i,j)
      enddo
      write ( iounit, * )
   enddo
 endif

end subroutine printmatrix
