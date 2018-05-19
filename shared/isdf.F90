!
! Weiwei Gao, Feb. 2018
!
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
! 
! This subroutine consists of three main steps:
! 
! 1. Prepare Zmtrx and Cmtrx
! 2. peinf%master use Zmtrx and Cmtrx to solve a linear quation to get zeta
! 3. Calculate Mmtrx = <zeta(i)| f^LDA(r) + V_coul(r,r') | zeta(j)>

subroutine isdf(gvec, pol_in, kpt, n_intp, intp, nspin, ncv, maxncv, invpairmap, &
      Cmtrx, Mmtrx, outputdbg, verbose)
 
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
  ! n_intp: number of interpolation points
  ! intp(n_intp): the index of interpolation points in full grids
  ! ncv(2): the number of wave function pairs for two spins
  integer, intent(in) :: n_intp, intp(n_intp), nspin, ncv(2), maxncv
  real(dp), intent(out) :: Cmtrx(n_intp, maxncv, nspin, kpt%nk)  
  real(dp), intent(out) :: Mmtrx(n_intp, n_intp, nspin, nspin, kpt%nk)
  ! if TRUE, write debuggin information to isdf_dbg.dat
  logical, intent(in) :: outputdbg 
  ! invpairmap maps the index of pair |cv> to c and v, which are the
  ! index of wave functions
  integer, intent(in) :: invpairmap(2, maxncv, nspin, kpt%nk)
  ! If verbose is true, then print out additional debug information
  logical, intent(in) :: verbose   

  ! --- Local variables ---
  ! Zmtrx(gvec%ldn, Nc*Nv, nspin, kpt%nk) should be calculated by every processor in wfn_grp 
  ! Cmtrx(n_intp, Nc*Nv, nspin, kpt%nk)
  ! Issue need to be: if there are N wfn_grp, does it mean these matrix need to be
  ! calculated by all of them at the same time, or can we distribute the workload
  ! later ??
  real(dp), allocatable :: &
   Zmtrx(:,:,:,:), &      ! Zmtrx on full grid
   Zmtrx_distr(:), &      ! distributed copy of Zmtrx
   Zmtrx_loc(:),   &      ! local copy of Zmtrx, on irreducible grid
   tmp_Zmtrx(:,:), &      ! temporary copy of Zmtrx
   zeta(:,:,:,:),  &
   ! matrices and vectors used for solving linear equations
   Amtrx(:,:), Bmtrx(:,:), Xmtrx(:,:), tmpmtrx(:,:,:,:), &
   rho_h(:)
  real(dp) :: diff, weight, qkt(3), vcharac(gvec%syms%ntrans), ccharac(gvec%syms%ntrans)
  ! counters and temporary integers
  integer :: ipt, ii, jj, iv, ic, icv, & 
             IVV, ICC, JVV, JCC, itrans, isp, ikp, errinfo, igrid, &
             iproc, tag, virep, cirep, ipair1, ipair2, isp1, isp2
  integer :: status, mpi_status(MPI_STATUS_SIZE)
  ! The index of w_grp, r_grp, processor that works on the calculation of
  ! zeta(:,:,nspin,nk)
  integer :: workwgrp, workrgrp, workproc
  ! Each processor in a w_grp store part of the wave function
  ! offset: index of grid point from which a processor start to store the wave function
  ! ncount: number of elements of wave functions that a processor stores
  integer, dimension(0:w_grp%npes-1) :: offset, ncount
  ! temporary dummy variable
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
  
  ! the number of real-space grid points in irreducible domain
  ngr = gvec%nr 
  ! the number of grid points in full domain
  ngf = gvec%nr * gvec%syms%ntrans 
  ! each processor store part of the Zmtrx
  ALLOCATE(Zmtrx_distr(w_grp%mydim)) 
  Cmtrx = zero ! initialize Cmtrx
  if(peinf%master) then
     !
     ALLOCATE(Zmtrx(ngf, maxncv, nspin, kpt%nk))
     ALLOCATE(zeta(ngf, n_intp, nspin, kpt%nk))
     ALLOCATE(Zmtrx_loc(ngr))
     ALLOCATE(Amtrx(n_intp, n_intp))
     ALLOCATE(Bmtrx(n_intp, ngf))
     ALLOCATE(Xmtrx(n_intp, ngf))
     Zmtrx = zero 
     !
     ! we want peinf%master to collect Cmtrx and Zmtrix from other processors in
     ! the w_grp
     !
     workrgrp = r_grp%mygr
     workwgrp = w_grp%mygr
     workproc = w_grp%inode ! the master proc should be zero?
     write(*,*) " Index of interpolation points: ", (intp(ii), ii=1,n_intp)
     !
  endif
  !
  call MPI_BCAST(workrgrp, 1, MPI_INTEGER, peinf%masterid, peinf%comm, errinfo)
  call MPI_BCAST(workwgrp, 1, MPI_INTEGER, peinf%masterid, peinf%comm, errinfo)
  call MPI_BCAST(workproc, 1, MPI_INTEGER, peinf%masterid, peinf%comm, errinfo)

  idum = 0 
  idum(w_grp%inode) = w_grp%offset + 1
  call MPI_ALLREDUCE(idum, offset, w_grp%npes, MPI_INTEGER, MPI_SUM, &
     w_grp%comm, errinfo)

  idum = 0
  idum(w_grp%inode) = w_grp%mydim
  call MPI_ALLREDUCE(idum, ncount, w_grp%npes, MPI_INTEGER, MPI_SUM, &
     w_grp%comm, errinfo)
  
  if(verbose) then
     ! 
     if(peinf%master) write(*,*) " in isdf() "
     write(*,*) 'w_grp%mygr', w_grp%mygr, ' w_grp%inode=',w_grp%inode, 'workgroup = ', workwgrp
     write(*,*) 'w_grp%mygr', w_grp%mygr, ' w_grp%inode=',w_grp%inode, 'workproc = ', workproc
     write(*,*) 'w_grp%mygr', w_grp%mygr, ' w_grp%inode=',w_grp%inode, " offset: ", (offset(ii),ii=0,w_grp%npes-1)
     write(*,*) 'w_grp%mygr', w_grp%mygr, ' w_grp%inode=',w_grp%inode, " ncount: ", (ncount(ii),ii=0,w_grp%npes-1)
     !
  endif

  do isp = 1, nspin
     !
     do ikp = 1, kpt%nk
        !
        if (verbose .and. peinf%master) &
           write(*,*) ' isp = ', isp, ', ikp = ', ikp
        !
        do icv = 1, ncv(isp)
           !
           Zmtrx_loc = 0.d0
           Zmtrx_distr = 0.d0
           !
           IVV = invpairmap(1,icv,isp,ikp)
           ICC = invpairmap(2,icv,isp,ikp)
           if (verbose .and. peinf%master) &
              write(*,'(a,i5,a,i5)') "ivv=", ivv," icc=", icc
           JVV = kpt%wfn(isp,ikp)%map(IVV)
           JCC = kpt%wfn(isp,ikp)%map(ICC)
           !
           ! copy phi_v to Zmtrx_loc
           !
           call dcopy(w_grp%mydim, &
             kpt%wfn(isp,ikp)%dwf(1,JVV),1,Zmtrx_distr(1),1)
           !
           ! calculate phi_v(r)*phi_c(r) and store it in Zmtrx_loc
           !
           call dmultiply_vec(w_grp%mydim, &
             kpt%wfn(isp,ikp)%dwf(1,JCC),Zmtrx_distr(1)) 
           !
           if (w_grp%mygr .eq. workwgrp .and. r_grp%mygr .eq. workrgrp) then
              !
              ! peinf%master collect Zmtrx from other processors in w_grp
              !
              do iproc = 0, w_grp%npes-1
                 !write(*,*) "w_grp%inode ", w_grp%inode, " loop iproc = ", iproc
                 !
                 tag = iproc
                 if (workproc .eq. iproc) then
                   ! peinf%master copy Zmtrx_distr to Zmtrx_loc
                   if (w_grp%inode .eq. iproc) &
                     call dcopy(ncount(iproc), Zmtrx_distr(1),1,&
                       Zmtrx_loc(offset(iproc)),1)
                 else
                   ! Each processor in workwgrp send part of Zmtrx to
                   ! peinf%master
                   if (w_grp%inode .eq. iproc) &
                     call MPI_send( Zmtrx_distr(1), ncount(iproc), &
                       MPI_DOUBLE, workproc, tag, w_grp%comm, errinfo)
                   if (w_grp%inode .eq. workproc) &
                     call MPI_recv( Zmtrx_loc(offset(iproc)), ncount(iproc), &
                       MPI_DOUBLE, iproc, tag, w_grp%comm, mpi_status, errinfo)
                 endif
                 !
              enddo ! iproc loop
              !
           endif
           !
           ! Get the symm group representation and characters of wfn |v> and |c>
           !
           virep = kpt%wfn(isp,ikp)%irep(JVV)
           cirep = kpt%wfn(isp,ikp)%irep(JCC)
           vcharac = gvec%syms%chi(virep,:)
           ccharac = gvec%syms%chi(cirep,:)
           !
           ! Zmtrx_loc store the data of grid points in reduced domain, now unfold it
           ! to full domain
           !
           if(peinf%master) then
              igrid = 0
              do ipt = 1, gvec%nr
                 ! peinf%master generate the Cmtrx
                 do itrans = 1, gvec%syms%ntrans
                   igrid = igrid + 1
                   Zmtrx(igrid, icv, isp, ikp) = Zmtrx_loc(ipt) * &
                     vcharac(itrans) * ccharac(itrans)
                 enddo ! itrans loop
              enddo ! ipt loop
              !
              ! Cmtrx store the data at interpolation points
              !
              do ipt = 1, n_intp 
                 Cmtrx(ipt, icv, isp, ikp) = Zmtrx(intp(ipt), icv, isp, ikp)
              enddo
           endif
        enddo ! icv loop
        
        if(peinf%master) then
           write(*, *) (Cmtrx(ii, 1, isp,ikp), ii = 1, n_intp)
        endif 
        
        ! ---------------------
        ! For each spin and ikp, calculate zeta
        ! zeta = Z * C^T * ( C * C^T)^-1  dimension: Ng*Nu
        ! Set A = C * C^T  dimension: Nu*(Nc*Nv) x (Nc*Nv)*Nu  = Nu*Nu
        !     B = C * Z^T  dimension: Nu*(Nc*Nv) * (Nc*Nv)*Ng  = Nu*Ng
        ! we solve the linear equation A*(zeta^T) = B to get zeta^T
        ! zeta^T = A^-1 * B = (C * C^T)^-1 * C * Z^T
        !
        !  Matrix dimensions:
        !  Cmtrx(n_intp,maxncv,:,:)
        !  Zmtrx(ngf,maxncv,:,:)
        !  zeta(ngf, n_intp,:,:)
        !  Amtrx(n_intp,n_intp) intermediate variable, store C * C^T
        !  Bmtrx(n_intp,ngf)    intermediate variable, store C * Z^T
        !  Xmtrx(n_intp,ngf)    intermediate variable, store zeta^T
        ! ---------------------
        if (peinf%master) then
           !
           ! calculate A = C * C^T
           !
           call dgemm('n','t',n_intp, n_intp, maxncv, one, &
             Cmtrx(1,1,isp,ikp), n_intp, Cmtrx(1,1,isp,ikp), n_intp, zero, Amtrx(1,1), n_intp)
           !
           if ( verbose ) then
              write(*,'(" A = C*C^T = ")')
              do ii = 1, n_intp
                 do jj = 1, n_intp
                    write(*,'(f12.8)',advance="no") Amtrx(ii,jj)
                 enddo
                 write(*,'()')
              enddo
           endif
           !
           ! calculate B = C * Z^T
           !
           call dgemm('n','t',n_intp, ngf, maxncv, one, &
             Cmtrx(1,1,isp,ikp), n_intp, Zmtrx(1,1,isp,ikp), ngf, zero, Bmtrx(1,1), n_intp)
           !
           if ( verbose ) then
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
           endif
           !
           ! solver the linear equation A * X = B
           !
           call dlinear_solver(n_intp, ngf, Amtrx, Bmtrx, Xmtrx)
           !
           ! Copy Xmtrx to zeta
           !
           ! if ( verbose ) write(*,*) 'zeta', isp,ikp
           do ii = 1, n_intp
              do jj = 1, ngf
                 zeta(jj, ii, isp, ikp) = Xmtrx(ii, jj)
                 ! if( verbose .and. ii.eq.1 .and. & 
                 !    jj <= 80 .and. mod(jj,8)==1) &
                 !    write(*,*) zeta(jj,ii,isp,ikp)
              enddo ! jj loop
           enddo ! ii loop
           !
        endif ! if peinf%master
     enddo ! ikp loop
  enddo ! isp loop
  !
  ! clean up all the allocatable variables
  !
  deallocate(Zmtrx_distr)
  if ( peinf%master ) then
     if ( verbose ) write(*,*) " deallocate arrays"
     deallocate( Zmtrx_loc )
     deallocate( Amtrx )
     deallocate( Bmtrx )
     deallocate( Xmtrx )
  endif

  if ( outputdbg .and. peinf%master ) then
     !
     write(*,*) "Output to isdf_dbg.dat"
     !
     open(dbgunit, file=dbg_filename, form='formatted',status='replace')
     !
     ! calculate zeta * C, which should in principle be a good approximation of Z
     !
     allocate( tmp_Zmtrx(ngf, maxncv) )
     if ( verbose ) then
        !
        do isp = 1, nspin
           do ikp = 1, kpt%nk
              !
              write(dbgunit, '(a,i2,a,i6)') " isp ", isp," ikp ", ikp
              !
              ! Calculate: tmp_Zmtrx(ngf,maxncv) = 
              !   sum_{n_intp} zeta(ngf,n_intp) * Cmtrx(n_intp,maxncv)
              !
              call dgemm('n','n', ngf, maxncv, n_intp, one, &
                zeta(1,1,isp,ikp), ngf, &
                Cmtrx(1,1,isp,ikp), n_intp, zero, tmp_Zmtrx, ngf)  
              !
              write(dbgunit,'(a)') "   ivv    icc      diff "
              !
              do icv = 1, ncv(isp)
                 ! calculate the relative difference between
                 ! the interpolated |vc> and the real |vc>
                 diff = zero
                 weight = zero
                 ivv = invpairmap(1,icv,isp,ikp)
                 icc = invpairmap(2,icv,isp,ikp)
                 do igrid = 1, ngf 
                    diff = diff + abs(Zmtrx(igrid,icv,isp,ikp)-tmp_Zmtrx(igrid,icv))
                    weight = weight + abs(Zmtrx(igrid,icv,isp,ikp))
                    !if ( (igrid < 160 .or. (igrid <9500 .and. igrid > 9000)) .and. &
                    !    mod(igrid,8) .eq. 1 ) & 
                    !  write(dbgunit,'(i10,2x,3e15.5)') &
                    !    igrid, &
                    !    Zmtrx(igrid,icv,isp,ikp), &
                    !    tmp_Zmtrx(igrid,icv), &
                    !    Zmtrx(igrid,icv,isp,ikp)-tmp_Zmtrx(igrid,icv)  ! For Debug
                 enddo
                 diff = diff/weight
                 write(dbgunit, '(2i6,e15.5)') ivv, icc, diff
              enddo ! icv loop
           enddo ! ikp loop
        enddo ! isp loop
        !
     endif
     !
     deallocate(tmp_Zmtrx)
     !
  endif ! if (outputdbg .and. peinf%master) 

  if ( peinf%master ) then
    deallocate(Zmtrx)
  endif
  call MPI_Barrier(peinf%comm, errinfo)
  !
  ! Now calculate <zeta_u(r,ispin)|V(r,r')|zeta_w(r',jspin)>, where u, w = 1, ..., n_intp
  ! and store it in Mmtrx(n_intp, n_intp)
  !
  Mmtrx = 0 ! Mmtrx dimension: Mmtrx(n_intp, n_intp, nspin, nspin)
  !
  ! only peinf%master calculate Mmtrx, then peinf%master send Mmtrx to
  ! all the other vectors
  !
  if ( peinf%master ) then
     !
     ! Only for tests of nonperiodic system. This should be updated later.
     qkt = 0
     ! Generate Coulomb potential
     allocate(rho_h(ngf)) ! note: gvec%nr is equal to w_grp%nr
     call dinitialize_FFT(peinf%inode, fft_box)
     call dcreate_coul_0D(gvec%bdot,qkt,fft_box)
     !
     do ikp = 1, kpt%nk
        !
        write( dbgunit, * ) " ikp = ", ikp, " Mmtrx = "
        !
        do ipair1 = 1, n_intp * nspin
           !
           ! copy zeta_ii to rho_h
           !
           if ( ipair1 > n_intp ) then
              isp1 = 2
              ii = ipair1 - n_intp
           else
              isp1 = 1
              ii = ipair1
           endif
           call dcopy( ngf, zeta( 1, ii, isp1, ikp ), 1, rho_h(1), 1 )
           !
           ! print out some debug info
           !
           if ( ii == 1 .and. verbose ) then
              write ( dbgunit, '("zeta(:,",3i4,")")' ) ii, isp1, ikp
              do igrid = 1, ngf
                 write(dbgunit,'(e15.5)',advance='no') zeta(igrid,ii,isp1,ikp)
                 if(mod(igrid,8)==0) write(dbgunit,'()')
              enddo
              write(dbgunit,'("rho_h(:)")')
              do igrid = 1, ngf
                 write(dbgunit,'(e15.5)',advance='no') rho_h(igrid)
                 if(mod(igrid,8)==0) write(dbgunit,'()')
              enddo
           endif
           !
           ! solve poisson equation to get: rho_h(r) = \int V_c(r,r') zeta_ii(r') dr'
           !
           call dfullpoisson(gvec, rho_h, dbgunit,.FALSE.)
           !
           ! print out some debug info
           !
           if ( ii == 1 .and. verbose ) then
              write(dbgunit,'("After poisson solver")')
              write(dbgunit,'("rho_h(:)")')
              do igrid = 1, ngf
                 write(dbgunit,'(e15.5)',advance='no') rho_h(igrid)
                 if(mod(igrid,8)==0) write(dbgunit,'()')
              enddo
           endif
           !
           ! Mmtrx is a symmetric matrix 
           ! So we start the loop from ipair1
           !
           do ipair2 = ipair1, n_intp * nspin
              !
              if ( ipair2 > n_intp) then
                 isp2 = 2
                 jj = ipair2 - n_intp
              else
                 isp2 = 1
                 jj = ipair2
              endif
              !
              ! calculate: \int zeta_jj(r) V_c(r,r') zeta_ii(r') drdr' 
              !          = \int zeta_jj(r) rho_h(r) dr
              !
              Mmtrx( ii, jj, isp1, isp2, ikp ) = &
                ddot( ngf, zeta( 1, jj, isp2, ikp ), 1, rho_h(1), 1) / gvec%hcub
              !
              ! Mmtrx is a symmetric matrix
              !
              Mmtrx( jj, ii, isp2, isp1, ikp ) = & 
                Mmtrx( ii, jj, isp1, isp2, ikp )
              if( verbose .and. jj <= 5 .or. jj > n_intp-5) & 
                 write(dbgunit,'(e14.5)',advance="no") Mmtrx(ii, jj, isp1, isp2, ikp)
              if( verbose .and. jj == 5 ) &
                 write(dbgunit,'(" ... ")',advance='no')
              !
           enddo ! ipair2
           !
           if(verbose) write(dbgunit,'()')
           !
        enddo ! ipair1
        !
     enddo ! ikp
     !
     deallocate(zeta) ! no longer needed
     close(dbgunit)
     !
  endif ! if(peinf%master) 
  
  call MPI_BCAST(Mmtrx, n_intp*n_intp*nspin*nspin*kpt%nk, MPI_DOUBLE, &
     peinf%masterid, MPI_COMM_WORLD, errinfo)
  call MPI_BCAST(Cmtrx, n_intp*maxncv*nspin*kpt%nk, MPI_DOUBLE, &
     peinf%masterid, MPI_COMM_WORLD, errinfo)

  return
end subroutine isdf

