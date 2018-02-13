! Issue need to be resolved: How do you make sure there are no duplicated
! interpolation points??

subroutine cvt(gvec, rho, nspin, n_intp, intp)

 use typedefs
 USE IFPORT   ! This is required for rand() if intel compiler is used
 implicit none

 ! -- PARAMETERS --
 ! number of maximum iterations to find interpolation points
 integer, parameter :: max_iter = 100000
 ! random seed for initializing random interpolation points 
 ! integer, parameter :: rseed = 12348
 integer :: rseed
 ! convergence threshold 
 real(dp), parameter :: tol_conv = 1e-3

 integer, intent(in) :: nspin
 ! gvec stores info about the real-space grids 
 type (gspace), intent(in) :: gvec 
 ! charge density in real-space grids, stored in irreducible wedge
 real(dp), intent(in) :: rho(gvec%nr, nspin)
 ! input parameter, the number of interpolation points needed
 integer, intent(in) :: n_intp
 ! outputs, the index of intp pts in the full grid points
 integer, intent(out) :: intp(n_intp)
 ! the real-space coordinates (in a.u.) of intp pts
 real(dp) :: pts(3,n_intp), newpts(3,n_intp), oldpts(3,n_intp)
 ! counters for different loops
 integer :: iter, ipt, ipt2, itran, igrid, ii, jj, kk, &
   itmp1, itmp2, iduplicate
 integer :: outdbg = 12345 ! unit of file for debugging
 ! full grid points and charge density on full grid points
 real(dp), allocatable :: fullgrid(:,:), fullrho(:)
 integer,allocatable :: icenter(:)
 real(dp) :: bounds(2,3), dist, weightedpos(3), weights, & 
  diff, vtmp(3), dist_tmp
 integer :: pt_tmp(3)

 open(outdbg, file="cvt_debug.dat", form='formatted', status='unknown')
 ! generate all points in the unfolded real space grid
 allocate(fullgrid(3,gvec%nr * gvec%syms%ntrans))
 ! get the charge density on the full grid
 allocate(fullrho(gvec%nr * gvec%syms%ntrans))
 allocate(icenter(gvec%nr * gvec%syms%ntrans))
 igrid = 0 ! counter for full-grid points
 do ipt = 1, gvec%nr 
    do itran = 1, gvec%syms%ntrans
       igrid = igrid + 1
       call unfold(gvec%r(1,ipt), &
         gvec%syms%trans(1,1,itran),gvec%shift(1),pt_tmp(1))
       do ii = 1,3
         fullgrid(ii,igrid) = (pt_tmp(ii) + gvec%shift(ii)) * gvec%step(ii)
       enddo
       fullrho(igrid) = rho(ipt,1)
       if (nspin .eq. 2) fullrho(igrid) = fullrho(igrid) + rho(ipt,2)
    enddo
 enddo
 
 ! find the lower and higher bounds of full-grid points
 do ii = 1,3
    bounds(1,ii) = minval(fullgrid(ii,:))
    bounds(2,ii) = maxval(fullgrid(ii,:))
 enddo
 ! write(outdbg,*) bounds(1,1:3)
 ! write(outdbg,*) bounds(2,1:3)

 ! generate initial guess of interpolation points
 write(outdbg,*) " Initial guess of interpolation points "
 !rseed = time()
 rseed = 1518543090
 write(outdbg,*) " rseed for random number generator is: ", rseed
 call srand(rseed)
 do ipt = 1, n_intp
    do ii = 1,3
        ! generate some random points in the full grids
        ! multiply by 0.9 to make sure these random points are inside the
        ! boundary
        pts(ii,ipt) = bounds(1,ii) + rand(0)*0.9*(bounds(2,ii)-bounds(1,ii)) 
    enddo
    ! print out the random interpolation points
    write(outdbg,'(3f10.4)') pts(1:3,ipt)
 enddo

 ! perform centroidal voronoi tesselation (CVT) algorithm
 ! For more details, see arXiv: 1711.01531v1 by Kun Dong, Wei Hu and Lin Lin
 ! (2017)
 newpts = pts
 oldpts = pts
 do iter = 1, max_iter
    ! for each point in the full grid, find which interpolation points is 
    ! the closest it. Then put the index of intp point to icenter(:)
    do igrid = 1, gvec%nr * gvec%syms%ntrans
       ! set a very large initial value for dist
       dist = 100.0 * ( maxval(bounds(2,:)) - minval(bounds(1,:)) )**2  
       icenter(igrid) = -1
       ! find which intp point is closest to the current grid point
       do ipt = 1, n_intp
          vtmp = newpts(:,ipt) - fullgrid(:,igrid)
          dist_tmp = dot_product(vtmp, vtmp)
          if (dist_tmp < dist) then
             dist = dist_tmp
             icenter(igrid) = ipt
          endif
       enddo 
    enddo

    ! Now update the interpolation pts
    diff = 0.0
    do ipt = 1, n_intp
       weightedpos(:) = 0.0
       weights = 0.0
       do igrid = 1, gvec%nr * gvec%syms%ntrans
         if (icenter(igrid) .ne. ipt) cycle
         weightedpos = weightedpos + fullgrid(:,igrid)*fullrho(igrid)
         weights = weights + fullrho(igrid)
       enddo
       ! update the new intp points with the centroid
       ! This is just a quick simple trick to avoid the case of weights == 0. It
       ! may be changed later with a better method.
       ! Simple minded fix, just generate a new randome points
       if(weights .lt. 1e-7) then
         do ii = 1,3
             newpts(ii,ipt) = bounds(1,ii) + rand(0)*0.9*(bounds(2,ii)-bounds(1,ii)) 
         enddo
       else
         newpts(:,ipt) = weightedpos/weights
       endif
       vtmp = newpts(:,ipt) - oldpts(:,ipt)
       diff = diff + sqrt(dot_product(vtmp,vtmp))
       !write(outdbg,'(8f8.3)') newpts(1:3,ipt), vtmp(1:3), &
       !  sqrt(dot_product(vtmp,vtmp)),weights
    enddo ! loop ipt
    write(*,'(i8,a,f18.12)') iter, " diff (a.u.) ", diff/n_intp
    if (diff < tol_conv) then ! conv threshold is satisfied, break the loop??
       write(*,*) " enter "
       pts = newpts
       iduplicate = -1
       ! For each pts, find the closest full-grid points, and store
       ! the index of the grid points to intp(:)
       intp = 0
       do ipt = 1, n_intp
          ! initialize dist to a very large variabl
          dist = 100.0 * ( maxval(bounds(2,:)) - minval(bounds(1,:)) )**2  
          ! loop over all the points in full-grid to find out which one is
          ! the closest to the current interpolation point
          do igrid = 1,gvec%nr * gvec%syms%ntrans
              vtmp = pts(:,ipt) - fullgrid(:,igrid)
              dist_tmp = dot_product(vtmp, vtmp)
              if (dist_tmp < dist) then
                  dist = dist_tmp
                  intp(ipt) = igrid
              endif
          enddo
          ! Now find out if there are any duplicated points (symmetry related)
          if(ipt .ge. 2) then
            do ipt2 = 1, ipt-1
               itmp1 = (intp(ipt2)-1)/gvec%syms%ntrans
               itmp2 = (intp(ipt)-1)/gvec%syms%ntrans
               ! If there are duplicated points, initialize it with a random point
               ! and continue 
               if(itmp2 .eq. itmp1) then
                  iduplicate = ipt2
                  do ii = 1,3 
                     newpts(ii,ipt) = bounds(1,ii) + &
                        rand(0)*0.9*(bounds(2,ii)-bounds(1,ii))
                  enddo
                  exit ! ipt2 loop
               endif
            enddo ! ipt2 loop
          endif
       enddo ! ipt loop
       ! if there are no duplicated points, then we break the iteration
       if(iduplicate < 0) exit 
    endif
    oldpts = newpts
 enddo ! iter
 
 do ipt = 1, n_intp
    write(outdbg,'(6f8.3)') pts(1:3,ipt), fullgrid(1:3,intp(ipt))
 enddo

 close(outdbg) ! close the dbg file
end subroutine cvt
