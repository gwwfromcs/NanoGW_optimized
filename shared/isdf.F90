! Weiwei Gao, Feb. 2018
! In the density fitting method, the products of wave functions are
! written as linear combinations of some interpolation vectors zeta(r).
! See equation 5 in J. Chem. Theory. Comput. 2017, 13, 5420-5431:
!
! \phi_i(r)\psi_j(r) = \sum^{Nu}_{u=1} \zeta_u(r) \phi_i(r_u) \psi_j(r_u)
! 
! <==> In the matrix form: Z = Zeta * C
!   _______________________________________________________________________________
!   |Matrix  |    Shape      |      Contents                                      |
!   |________|_______________|____________________________________________________|
!   |Z       | Ng * (Nc*Nv)  |  {..., \phi_i(r)\psi_j(r), ...}                    |
!   |        |               |   i=1...Nv, j=1...Nv, r are full-grid points       |
!   |________|_______________|____________________________________________________|
!   |Zeta    | Ng *  Nu      |  {..., \zeta_u(r), ...}                            |
!   |        |               |   u=1...Nu                                         |
!   |________|_______________|____________________________________________________|
!   |C       | Nu * (Nc*Nv)  |  {..., \phi_i(ru)*\psi_j(ru), ...}                 |
!   |________|_______________|____________________________________________________|
!
! This subroutine calculates the interpolation vectors zeta_u(r)
! 
! n_intp : the number of interpolation vectors or points, or Nu
! intp   : the index of interpolation points in the full grid
! zeta(gvec%nr*gvec%sym, n_intp) : the interpolation functions function
! kpt%wfn(isp,ikp)%dwf(:,:) : store the wavefunctions \phi_i, \psi_j

! For now, this is only partially parallelized
! 1. each processor evaluates part of the matrix Z and part of the matrix C
! 2. the matrix Z and C are collected by the master node
! 3. the master node solve the linear equation to get zeta

subroutine isdf(gvec, pol_in, kpt, n_intp, intp, nspin, ncv, maxncv, Cmtrx, Mmtrx, &
     outputdbg)
 
  use typedefs
  use mpi_module
  use myconstants
  use fft_module
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif

  type(gspace), intent(in) :: gvec
  type(polinfo), intent(in), dimension(2) :: pol_in
  type(kptinfo), intent(in) :: kpt
  integer, intent(in) :: n_intp, intp(n_intp), nspin, ncv(2), maxncv
  real(dp), intent(out) :: Cmtrx(n_intp, maxncv, nspin, kpt%nk)  ! Cmtrx on interpolation points
  real(dp), intent(out) :: Mmtrx(n_intp, n_intp, nspin, kpt%nk)
  logical, intent(in) :: outputdbg ! if TRUE, write debuggin information to isdf_dbg.dat


  ! Zmtrx(gvec%ldn, Nc*Nv, nspin, kpt%nk) should be calculated by every processor in wfn_grp 
  ! Cmtrx(n_intp, Nc*Nv, nspin, kpt%nk)
  ! Issue: if there are N wfn_grp, does it mean these matrix need to be
  ! calculated by all of them at the same time, or can we distribute the workload
  ! later ??
  real(dp), allocatable :: &
   Zmtrx(:,:,:,:), &      ! Zmtrx on full grid
   Zmtrx_distr(:), &      ! distributed copy of Zmtrx
   Zmtrx_loc(:),   &      ! local copy of Zmtrx, on irreducible grid
   tmp_Zmtrx(:,:), &      ! temporary copy of Zmtrx
   tmp_Cmtrx(:,:), &      ! temporary copy of Cmtrx
   zeta(:,:,:,:),  &
   ! matrices and vectors used for solving linear equations
   Amtrx(:,:), Bmtrx(:,:), Xmtrx(:,:), tmpmtrx(:,:,:,:), &
   rho_h(:)
  real(dp) :: diff, weight, qkt(3)
  ! counters and temporary integers
  integer :: ipt, ii, jj, iv, ic, icv, & 
             IVV, ICC, JVV, JCC, itrans, isp, ikp, errinfo, igrid, &
             iproc, tag
  integer :: status(MPI_STATUS_SIZE)
  integer :: workgroup, workproc
  integer, dimension(0:w_grp%npes-1) :: offset, ncount
  integer, dimension(0:w_grp%npes-1) :: idum
  ! the number of grid points in irreducible wedge, ngr = gvec%nr
  integer :: ngr 
  ! the number of full grid point, ngf = ngr * (# of sym operations)
  integer :: ngf 

  ! variables for debug and test of accuracy 
  character(50) :: dbg_filename = "isdf_dbg.dat"
  integer :: dbgunit = 20171130
  ! external functions
  real(dp), external :: ddot

  ngr = gvec%nr
  ngf = gvec%nr * gvec%syms%ntrans
  allocate(Zmtrx_distr(w_grp%mydim))
  Cmtrx = 0
  if(peinf%master) then
    allocate(Zmtrx(ngf, maxncv, nspin, kpt%nk))
    allocate(zeta(ngf,n_intp,nspin,kpt%nk))
    allocate(Zmtrx_loc(ngr))
    allocate(Amtrx(n_intp, n_intp))
    allocate(Bmtrx(n_intp, ngf))
    allocate(Xmtrx(n_intp, ngf))
    Zmtrx = 0
    workgroup = w_grp%mygr 
    workproc  = w_grp%inode
    write(*,*) "intp: ", (intp(ii), ii=1,n_intp)
  endif
  call MPI_BCAST(workgroup, 1, MPI_INTEGER, peinf%masterid, peinf%comm, errinfo)
  call MPI_BCAST(workproc, 1, MPI_INTEGER, peinf%masterid, peinf%comm, errinfo)

  idum = 0 
  idum(w_grp%inode) = w_grp%offset + 1
  call MPI_ALLREDUCE(idum, offset, w_grp%npes, MPI_INTEGER, MPI_SUM, &
    w_grp%comm, errinfo)

  idum = 0
  idum(w_grp%inode) = w_grp%mydim
  call MPI_ALLREDUCE(idum, ncount, w_grp%npes, MPI_INTEGER, MPI_SUM, &
    w_grp%comm, errinfo)
  
  !write(*,*) 'w_grp%mygr', w_grp%mygr, ' w_grp%inode=',w_grp%inode, 'workgroup = ', workgroup
  !write(*,*) 'w_grp%mygr', w_grp%mygr, ' w_grp%inode=',w_grp%inode, 'workproc = ', workproc
  !write(*,*) 'w_grp%mygr', w_grp%mygr, ' w_grp%inode=',w_grp%inode, " offset: ", (offset(ii),ii=0,w_grp%npes-1)
  !write(*,*) 'w_grp%mygr', w_grp%mygr, ' w_grp%inode=',w_grp%inode, " ncount: ", (ncount(ii),ii=0,w_grp%npes-1)

  do isp = 1, nspin
     do ikp = 1, kpt%nk
        icv = 0
        if (peinf%master) write(*,*) ' isp = ', isp, ', ikp = ', ikp
        do iv = 1, pol_in(isp)%nval
           do ic = 1, pol_in(isp)%ncond
              !if (peinf%master) write(*,'(a,i5,a,i5)') "iv=", iv," ic=", ic
              Zmtrx_loc = 0.d0
              Zmtrx_distr = 0.d0
              icv = icv+1
              IVV = pol_in(isp)%vmap(iv)
              JVV = kpt%wfn(isp,ikp)%map(IVV)

              ! copy phi_v to Zmtrx_loc
              call dcopy(w_grp%mydim, &
                kpt%wfn(isp,ikp)%dwf(1,JVV),1,Zmtrx_distr(1),1)
              ICC = pol_in(isp)%cmap(ic)
              JCC = kpt%wfn(isp,ikp)%map(ICC)
              ! calculate phi_v(r)*phi_c(r) and store it in Zmtrx_loc
              call dmultiply_vec(w_grp%mydim, &
                kpt%wfn(isp,ikp)%dwf(1,JCC),Zmtrx_distr(1)) 

              if (w_grp%mygr .eq. workgroup) then
                 ! peinf%master collect Zmtrx from other processors in w_grp
                 do iproc = 0, w_grp%npes-1
                    !write(*,*) "w_grp%inode ", w_grp%inode, " loop iproc = ", iproc
                    tag = iproc
                    if (workproc .eq. iproc) then
                      if (w_grp%inode .eq. iproc) &
                        call dcopy(ncount(iproc), Zmtrx_distr(1),1,&
                          Zmtrx_loc(offset(iproc)),1)
                    else
                      if (w_grp%inode .eq. iproc) &
                        call MPI_send( Zmtrx_distr(1), ncount(iproc), &
                          MPI_DOUBLE, workproc, tag, w_grp%comm, errinfo)
                      if (w_grp%inode .eq. workproc) &
                        call MPI_recv( Zmtrx_loc(offset(iproc)), ncount(iproc), &
                          MPI_DOUBLE, iproc, tag, w_grp%comm, status, errinfo)
                    endif
                 enddo
              endif

              !if(peinf%master) write(*,*) ' d'
              if(peinf%master) then
                 igrid = 0
                 do ipt = 1, gvec%nr
                    ! peinf%master generate the Cmtrx
                    do itrans = 1, gvec%syms%ntrans
                      igrid = igrid + 1
                      Zmtrx(igrid, icv, isp, ikp) = Zmtrx_loc(ipt)
                    enddo ! itrans
                 enddo ! ipt
                 
                 do ipt = 1, n_intp
                    Cmtrx(ipt, icv, isp, ikp) = Zmtrx(intp(ipt), icv, isp, ikp)
                 enddo
              endif
           enddo ! ic loop
        enddo ! iv loop
        
        ! For each spin and ikp, calculate zeta
        ! zeta = Z * C^T * ( C * C^T)^-1  dimension: Ng*Nu
        ! Set A = C * C^T  dimension: Nu*(Nc*Nv) x (Nc*Nv)*Nu  = Nu*Nu
        !     B = C * Z^T  dimension: Nu*(Nc*Nv) * (Nc*Nv)*Ng  = Nu*Ng
        ! we solve the linear equation A*(zeta^T) = B to get zeta^T
        ! zeta^T = A^-1 * B = (C * C^T)^-1 * C * Z^T
        
        !  Matrix dimensions:
        !  Cmtrx(n_intp,maxncv,:,:)
        !  Zmtrx(ngf,maxncv,:,:)
        !  zeta(ngf, n_intp,:,:)
        !  Amtrx(n_intp,n_intp) intermediate variable, store C * C^T
        !  Bmtrx(n_intp,ngf)    intermediate variable, store C * Z^T
        !  Xmtrx(n_intp,ngf)    intermediate variable, store zeta^T
        if (peinf%master) then
           ! calculate A = C * C^T
           call dgemm('n','t',n_intp, n_intp, maxncv, one, &
             Cmtrx(1,1,isp,ikp), n_intp, Cmtrx(1,1,isp,ikp), n_intp, zero, Amtrx(1,1), n_intp)
           write(*,'(" A = C*C^T = ")')
           do ii = 1, n_intp
              do jj = 1, n_intp
                 write(*,'(f12.8)',advance="no") Amtrx(ii,jj)
              enddo
              write(*,'()')
           enddo
           ! calculate B = C * Z^T
           call dgemm('n','t',n_intp, ngf, maxncv, one, &
             Cmtrx(1,1,isp,ikp), n_intp, Zmtrx(1,1,isp,ikp), ngf, zero, Bmtrx(1,1), n_intp)
           write(*,'(" B = C*Z^T = ")')
           do ii = 1, n_intp
              do jj = 1, 12
                 write(*,'(f12.8)',advance="no") Bmtrx(ii,jj)
              enddo
              write(*,'(" ... ")',advance="no")
              do jj = 13, 25
                 write(*,'(f12.8)',advance="no") Bmtrx(ii,jj)
              enddo
              write(*,'()')
           enddo
           ! solver the linear equation A * X = B
           call dlinear_solver(n_intp, ngf, Amtrx, Bmtrx, Xmtrx)
           ! Copy Xmtrx to zeta
           write(*,*) 'zeta', isp,ikp
           do ii = 1, n_intp
              do jj = 1, ngf
                 zeta(jj, ii, isp, ikp) = Xmtrx(ii, jj)
                 if(ii.eq.1 .and. jj <= 80 .and. mod(jj,8)==1) &
                     write(*,*) zeta(jj,ii,isp,ikp)
              enddo
           enddo
        endif
     enddo ! ikp loop
  enddo ! isp loop

  write(*,*) 'peinf%inode', peinf%inode
  ! clean up all the allocatable variables
  if (peinf%master) then
     write(*,*) " deallocate arrays"
     deallocate(Zmtrx_loc)
     deallocate(Amtrx)
     deallocate(Bmtrx)
     deallocate(Xmtrx)
  endif

  if (outputdbg .and. peinf%master) then
     write(*,*) "output to zeta_dbg.dat"
     open(dbgunit, file=dbg_filename, form='formatted',status='replace')
     !do igrid = 1, 5
     !   write(dbgunit,'(2x,3f20.8)') zeta(igrid,1,1,1), zeta(igrid,2,1,1), zeta(igrid,3,1,1)
     !enddo
     !write(dbgunit,'(" ... ")')
     !do igrid = 5, 0, -1
     !   write(dbgunit,'(2x,3f20.8)') zeta(ngf-igrid,1,1,1), zeta(ngf-igrid,2,1,1), zeta(ngf-igrid,3,1,1)
     !enddo
     ! calculate zeta * C, which should in principle be a good approximation of Z
     allocate(Xmtrx(ngf,n_intp))
     allocate(tmp_Zmtrx(ngf,maxncv))
     allocate(tmp_Cmtrx(n_intp,maxncv))
     if (.false.) then
     do isp = 1, nspin
        do ikp = 1, kpt%nk
           write(dbgunit, '(a,i2,a,i6)') " isp ", isp," ikp ", ikp
           Xmtrx(:,:) = zeta(:,:,isp,ikp)        ! Xmtrx(ngf,n_intp)
           tmp_Cmtrx(:,:) = Cmtrx(:,:,isp,ikp)   ! tmp_Cmtrx(n_intp, maxncv)
           ! tmp_Zmtrx(ngf,maxncv) = 
           ! Xmtrx(ngf,n_intp) * tmp_Cmtrx(n_intp,maxncv)
           call dgemm('n','n', ngf, maxncv, n_intp, one, &
             Xmtrx, ngf, tmp_Cmtrx, n_intp, zero, tmp_Zmtrx, ngf)  
           icv = 0
           write(dbgunit,'(a)') "   iv    ic      diff "
           do iv = 1, pol_in(isp)%nval
              do ic = 1, pol_in(isp)%ncond
                 icv = icv + 1
                 diff = 0
                 weight = 0
                 do igrid = 1, ngf 
                    diff = diff + abs(Zmtrx(igrid,icv,isp,ikp)-tmp_Zmtrx(igrid,icv))
                    weight = weight + abs(Zmtrx(igrid,icv,isp,ikp))
                    if ( (igrid < 160 .or. (igrid <600 .and. igrid > 500)) .and. &
                        mod(igrid,8) .eq. 1 ) &
                      write(dbgunit,'(i10,2x,3e15.5)') &
                        igrid, &
                        Zmtrx(igrid,icv,isp,ikp), &
                        tmp_Zmtrx(igrid,icv), &
                        Zmtrx(igrid,icv,isp,ikp)-tmp_Zmtrx(igrid,icv)
                 enddo
                 diff = diff/weight
                 write(dbgunit, '(2i6,e15.5)') iv, ic, diff
              enddo
           enddo
        enddo
     enddo
     endif
     deallocate(Xmtrx)
     deallocate(tmp_Zmtrx)
     deallocate(tmp_Cmtrx)
  endif ! if (outputdbg .and. peinf%master) 

  deallocate(Zmtrx_distr)
  if (peinf%master) then
    deallocate(Zmtrx)
  endif
  call MPI_Barrier(peinf%comm, errinfo)

  ! Now calculate <zeta_u(r)|V(r,r')|zeta_w(r')>, where u, w = 1, ..., n_intp
  ! and store it in Mmtrx(n_intp, n_intp)
  
  Mmtrx = 0 ! Mmtrx dimension: Mmtrx(n_intp, n_intp)
  ! only peinf%master calculate Mmtrx, then use MPI_ALLREDUCE to send Mmtrx to
  ! all the other vectors
  if(peinf%master) then
     ! Only for tests of nonperiodic system. This should be updated later.
     qkt = 0
     ! Generate Coulomb potential
     allocate(rho_h(ngf)) ! note: gvec%nr is equal to w_grp%nr
     call dinitialize_FFT(peinf%inode, fft_box)
     write(dbgunit,*) 'fft_box%scale = ', fft_box%scale
     call dcreate_coul_0D(gvec%bdot,qkt,fft_box)
     do isp = 1, nspin
        do ikp = 1, kpt%nk
           write(dbgunit,*) " isp ", isp, " ikp ", ikp
           write(dbgunit,*) " Mmtrx = "
           do ii = 1, n_intp
              ! copy zeta_ii to rho_h
              call dcopy(ngf,zeta(1,ii,isp,ikp),1,rho_h(1),1)
              if(ii == 1 ) then
                 write(dbgunit,'("zeta(:,",3i4,")")') ii, isp, ikp
                 do igrid = 1, ngf
                    write(dbgunit,'(e15.5)',advance='no') zeta(igrid,ii,isp,ikp)
                    if(mod(igrid,8)==0) write(dbgunit,'()')
                 enddo
                 write(dbgunit,'("rho_h(:)")')
                 do igrid = 1, ngf
                    write(dbgunit,'(e15.5)',advance='no') rho_h(igrid)
                    if(mod(igrid,8)==0) write(dbgunit,'()')
                 enddo
              endif
              ! solve poisson equation to get: rho_h(r) = \int V_c(r,r') zeta_ii(r') dr'
              if(ii == 1) then
                  write(dbgunit,'(3i5)') isp,ikp,ii
                  call dfullpoisson(gvec, rho_h, dbgunit,.FALSE.)
              else 
                  call dfullpoisson(gvec, rho_h, dbgunit,.FALSE.)
              endif
              if(ii == 1 .and. .false.) then
                 write(dbgunit,'("After poisson solver")')
                 write(dbgunit,'("rho_h(:)")')
                 do igrid = 1, ngf
                    write(dbgunit,'(e15.5)',advance='no') rho_h(igrid)
                    if(mod(igrid,8)==0) write(dbgunit,'()')
                 enddo
              endif
              do jj = 1, n_intp
                 ! calculate: \int zeta_jj(r) V_c(r,r') zeta_ii(r') drdr' 
                 !          = \int zeta_jj(r) rho_h(r) dr
                 Mmtrx(ii,jj,isp,ikp) = &
                   ddot(ngf,zeta(1,jj,isp,ikp),1,rho_h(1),1) /gvec%hcub
                 if(jj<=5 .or. jj>n_intp-5) write(dbgunit,'(e14.5)',advance="no") Mmtrx(ii, jj, isp, ikp)
                 if(jj == 5 ) write(dbgunit,'(" ... ")',advance='no')
              enddo ! jj
              write(*,'()')
           enddo ! ii
        enddo ! ikp
     enddo ! isp

     close(dbgunit)
  endif ! if(peinf%master) 
  allocate(tmpmtrx(n_intp, n_intp, nspin, kpt%nk))
  tmpmtrx = 0
  call MPI_ALLREDUCE(Mmtrx, tmpmtrx, n_intp*n_intp*nspin*kpt%nk, &
     MPI_DOUBLE, MPI_SUM, peinf%comm, errinfo)
  Mmtrx = tmpmtrx
  deallocate(tmpmtrx)
  allocate(tmpmtrx(n_intp, maxncv, nspin, kpt%nk))
  tmpmtrx = 0
  call MPI_ALLREDUCE(Cmtrx, tmpmtrx, n_intp*maxncv*nspin*kpt%nk, &
     MPI_DOUBLE, MPI_SUM, peinf%comm, errinfo)
  Cmtrx = tmpmtrx
  deallocate(tmpmtrx)

  return
end subroutine isdf

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


