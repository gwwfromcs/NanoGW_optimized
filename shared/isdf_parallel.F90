!
! Weiwei Gao, Feb. 2018
!
! In the density fitting method, the products of wave functions are
! written as linear combinations of some interpolation vectors zeta(r).
! See equation 5 in J. Chem. Theory. Comput. 2017, 13, 5420-5431:
!
! \phi_i(r)\psi_j(r) = \sum^{Nu}_{u=1} \zeta_u(r) \phi_i(r_u) \psi_j(r_u)
! 
! <==> written in a matrix form: Z = Zeta * C
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
!   |        |               |                                                    |
!   |________|_______________|____________________________________________________|
!   ______________________________________________________________________
!   |Matrix  |                  Parallelization                          |
!   |________|___________________________________________________________|
!   |Z       |             Each proc store part of Z                     |
!   |        |         Zmtrx(ngfl, Nc*Nv), ngfl = mydim*ntrans           |
!   |________|___________________________________________________________|
!   |Zeta    |             Each proc store part of zeta                  |
!   |        |         zeta(ngfl, Nu), ngfl = mydim*ntrans               |
!   |________|___________________________________________________________|
!   |C       |             Every proc store a full copy of C             |
!   |        |              Cmtrx ( Nu, Nc*Nv )                          |
!   |________|___________________________________________________________|
!
!
! This subroutine calculates the interpolation vectors zeta_u(r)
! 
! n_intp : the number of interpolation vectors or points, or Nu
! intp   : the index of interpolation points in the full grid
! zeta(gvec%nr*gvec%sym, n_intp) : the interpolation vectors 
! kpt%wfn(isp,ikp)%dwf(:,:) : store the wavefunctions \phi_i, \psi_j
!
! For now, this is only partially parallelized
! 1. each processor evaluates part of the matrix Z and part of the matrix C
! 2. collect matrix C from different procs. distribute matrix Z to procs.
! 3. every proc in w_grp solve part of the linear equations to get zeta
! 
! This subroutine consists of three main steps:
! 
! 1. Prepare Zmtrx and Cmtrx
! 2. solve a linear quation to get zeta
! 3. Calculate Mmtrx(1) = <zeta(i)| V_coul(r,r') | zeta(j)> and
!    Mmtrx(2) = <zeta(i)| f^LDA(r) | zeta(j)>

subroutine isdf_parallel(gvec, pol_in, kpt, n_intp, intp, & 
      nspin, ncv, maxncv, invpairmap, kflag, &
      Cmtrx, Mmtrx, verbose)
 
  use typedefs
  use mpi_module
  use myconstants
  use fft_module
  use xc_functionals
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
  integer, intent(in) :: n_intp, intp(n_intp), nspin, ncv(2), maxncv, kflag
  real(dp), intent(out) :: Cmtrx(n_intp, maxncv, nspin, kpt%nk)  
  real(dp), intent(out) :: Mmtrx(n_intp, n_intp, nspin, nspin, kpt%nk, 2)
  ! invpairmap maps the index of pair |cv> to c and v, which are the
  ! index of wave functions
  integer, intent(in) :: invpairmap(2, maxncv, nspin, kpt%nk)
  ! If verbose is true, then print out additional debug information
  logical, intent(in) :: verbose   

  ! --- Local variables ---
  ! Zmtrx(gvec%mydim*ntrans, Nc*Nv, nspin, kpt%nk)
  ! Cmtrx(n_intp, Nc*Nv, nspin, kpt%nk)
  ! Issue need to be addressed: if there are N wfn_grp, does it mean these matrix need to be
  ! calculated by all processors at the same time, or can we distribute the workload
  ! later ??
  real(dp), allocatable :: &
   Zmtrx(:,:,:,:),   &      ! Zmtrx on full grid
   Zmtrx_distr(:),   &      ! distributed copy of Zmtrx
   Zmtrx_loc(:),     &      ! local copy of Zmtrx, on irreducible grid
   tmp_Zmtrx(:),     &      ! temporary copy of Zmtrx
   zeta(:,:,:,:),    &
   fxc(:,:,:), fxc_loc(:,:,:), fzeta(:,:), &
   ! matrices and vectors used for solving linear equations
   Amtrx(:,:), Bmtrx(:,:), Xmtrx(:,:), tmpmtrx(:,:,:,:), &
   rho_h(:), rho_h_local(:)
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
  integer :: ngf, ngfl, istart, iend, iptf, iptr, ioff1, ioff2

  ! variables for debug and test of accuracy 
  character(50) :: dbg_filename = "isdf_dbg.dat"
  integer :: dbgunit = 20171130
  ! external functions
  real(dp), external :: ddot
  type(xc_type) :: xc_lda
  
  !
  workrgrp = 0
  workwgrp = 0
  call MPI_BCAST( workrgrp, 1, MPI_INTEGER, peinf%masterid, peinf%comm, errinfo )
  call MPI_BCAST( workwgrp, 1, MPI_INTEGER, peinf%masterid, peinf%comm, errinfo )
  !
  if ( workrgrp == r_grp%mygr .and. workwgrp == w_grp%mygr ) then
     !
     if( w_grp%master .and. verbose ) then
        !
        write(*,*) "call isdf_parallel(), write debug info to ", dbg_filename
        !
        open(dbgunit, file=dbg_filename, form='formatted', status='replace')
        !
     endif
     !
     ! the number of real-space grid points in irreducible domain
     ngr = gvec%nr 
     ! the number of grid points in full domain
     ngf = gvec%nr * gvec%syms%ntrans 
     ! the number of grid points stored in each processor
     ngfl = w_grp%mydim * gvec%syms%ntrans
     !
     ! each processor store part of the Zmtrx
     !
     ALLOCATE( Zmtrx( ngfl, maxncv, nspin, kpt%nk) )
     ALLOCATE( Zmtrx_loc( w_grp%mydim ) )
     ALLOCATE( zeta ( ngfl, n_intp, nspin, kpt%nk) )
     ALLOCATE( Amtrx( n_intp, n_intp) )
     ALLOCATE( Bmtrx( n_intp, ngfl) )
     ALLOCATE( Xmtrx( n_intp, ngfl) )
     !
     ! initialize matrices to zero
     !
     Cmtrx = zero
     Zmtrx = zero 
     zeta  = zero
     Amtrx = zero
     Bmtrx = zero
     Xmtrx = zero
     !
     if ( w_grp%master .and. verbose ) then
        write(dbgunit, *) " Index of interpolation points: ", ( intp(ii), ii=1, n_intp )
     endif
     !
     idum = 0 
     idum(w_grp%inode) = w_grp%offset + 1
     call MPI_ALLREDUCE(idum, offset, w_grp%npes, MPI_INTEGER, MPI_SUM, &
        w_grp%comm, errinfo)
     !
     idum = 0
     idum(w_grp%inode) = w_grp%mydim
     call MPI_ALLREDUCE(idum, ncount, w_grp%npes, MPI_INTEGER, MPI_SUM, &
        w_grp%comm, errinfo)
     !
     if( verbose .and. w_grp%master ) then
        !
        write(dbgunit, *) " in isdf() "
        write(dbgunit, *) " w_grp%mygr ", w_grp%mygr, " w_grp%inode = ", w_grp%inode, &
        " workwgrp = ", workwgrp, &
        " workrgrp = ", workrgrp, &
        " offset: ", ( offset(ii), ii=0, w_grp%npes-1), &
        " ncount: ", ( ncount(ii), ii=0, w_grp%npes-1)
        !
     endif
     !
     do isp = 1, nspin
        !
        do ikp = 1, kpt%nk
           !
           if (verbose .and. w_grp%master) &
              write(*, *) ' isp = ', isp, ', ikp = ', ikp
           !
           do icv = 1, ncv(isp)
              !
              IVV = invpairmap(1,icv,isp,ikp)
              ICC = invpairmap(2,icv,isp,ikp)
              if (verbose .and. w_grp%master) &
                 write(*, '(a,i5,a,i5)') " ivv = ", ivv, " icc = ", icc
              JVV = kpt%wfn(isp,ikp)%map(IVV)
              JCC = kpt%wfn(isp,ikp)%map(ICC)
              !
              ! copy phi_v to Zmtrx
              !
              call dcopy(w_grp%mydim, &
                kpt%wfn(isp,ikp)%dwf(1,JVV),1,Zmtrx(1,icv,isp,ikp),1)
              !
              ! calculate phi_v(r)*phi_c(r) and store it in Zmtrx
              !
              call dmultiply_vec(w_grp%mydim, &
                kpt%wfn(isp,ikp)%dwf(1,JCC),Zmtrx(1,icv,isp,ikp)) 
              !
              ! Get the symm group representation and characters of wfn |v> and |c>
              !
              virep = kpt%wfn(isp,ikp)%irep(JVV)
              cirep = kpt%wfn(isp,ikp)%irep(JCC)
              vcharac = gvec%syms%chi(virep,:)
              ccharac = gvec%syms%chi(cirep,:)
              !
              ! Zmtrx(1:w_grp%mydim) store the data of grid points in reduced domain, now unfold it
              ! to full domain
              !
              Zmtrx_loc(1:w_grp%mydim) = Zmtrx(1:w_grp%mydim, icv, isp, ikp)
              igrid = 0
              do ipt = 1, w_grp%mydim
                do itrans = 1, gvec%syms%ntrans
                  igrid = igrid + 1
                  Zmtrx(igrid, icv, isp, ikp) = Zmtrx_loc(ipt) * &
                    vcharac(itrans) * ccharac(itrans)
                enddo
              enddo
              !
              ! Cmtrx store the data at interpolation points
              !
              do ipt = 1, n_intp 
                 iptf = intp(ipt)
                 ioff1 = ( offset(w_grp%inode)-1 ) * gvec%syms%ntrans
                 ioff2 = ngf
                 if ( w_grp%inode + 1 < w_grp%npes ) then
                   ioff2 = ( offset(w_grp%inode+1)-1 ) * gvec%syms%ntrans
                 endif
                 if ( iptf .ge. ioff1 + 1 &
                      .and. &
                      iptf .lt. ioff2 + 1 ) then
                   jj = iptf - ioff1 
                   ! if ( icv .eq. 1 .and. verbose ) write(*, '(4(a,i5))') &
                   !    ", wgrp%inode ", w_grp%inode, ", ipt ", ipt, &
                   !    ", iptf ", iptf, ", jj ", jj
                   Cmtrx(ipt, icv, isp, ikp) = Zmtrx(jj, icv, isp, ikp)
                 endif
              enddo
           enddo ! icv loop
           !
           call MPI_ALLREDUCE ( MPI_IN_PLACE, Cmtrx(1,1,isp,ikp), n_intp * ncv(isp), &
              MPI_DOUBLE, MPI_SUM, w_grp%comm, errinfo )
           if(w_grp%master .and. verbose) then
              write(dbgunit, *) " Cmtrx = "
              call printmatrix ( Cmtrx(1,1,isp,ikp), n_intp, ncv(isp), dbgunit )
           endif
           !
           ! ---------------------
           ! For each spin and ikp, calculate zeta
           ! zeta = Z * C^T * ( C * C^T)^-1  dimension: Ng*Nu
           ! Set A = C * C^T  dimension: Nu*(Nc*Nv) * (Nc*Nv)*Nu  = Nu*Nu
           !     B = C * Z^T  dimension: Nu*(Nc*Nv) * (Nc*Nv)*Ng  = Nu*Ng
           ! we solve the linear equation A*(zeta^T) = B to get zeta^T
           ! zeta^T = A^-1 * B = (C * C^T)^-1 * C * Z^T
           !
           !  Matrix dimensions:
           !  Cmtrx(n_intp,maxncv,:,:)
           !  Zmtrx(ngfl,maxncv,:,:) where ngfl = w_grp%mydim * ntrans
           !  zeta(ngfl, n_intp,:,:)
           !  Amtrx(n_intp,n_intp)  intermediate variable, store C * C^T
           !  Bmtrx(n_intp,ngfl)    intermediate variable, store C * Z^T
           !  Xmtrx(n_intp,ngfl)    intermediate variable, store zeta^T
           ! ---------------------
           !
           ! calculate A = C * C^T
           !
           call dgemm('n','t',n_intp, n_intp, maxncv, one, &
             Cmtrx(1,1,isp,ikp), n_intp, Cmtrx(1,1,isp,ikp), n_intp, zero, Amtrx(1,1), n_intp)
           !
           if ( verbose .and. w_grp%master ) then
              write(dbgunit,'(" A = C*C^T = ")')
              call printmatrix ( Amtrx(1,1), n_intp, n_intp, dbgunit )
           endif
           !
           ! calculate B = C * Z^T
           !
           call dgemm( 'n','t',n_intp, ngfl, maxncv, one, &
             Cmtrx(1,1,isp,ikp), n_intp, Zmtrx(1,1,isp,ikp), ngfl, zero, Bmtrx(1,1), n_intp )
           !
           if ( verbose .and. w_grp%master ) then
              write(dbgunit,'(" B = C*Z^T = ")')
              call printmatrix ( Bmtrx(1,1), n_intp, ngfl, dbgunit )
           endif
           !
           ! solver the linear equation A * X = B
           !
           call dlinear_solver( n_intp, ngfl, Amtrx, Bmtrx, Xmtrx )
           !
           ! Copy Xmtrx to zeta
           !
           ! if ( verbose ) write(*,*) 'zeta', isp,ikp
           do ii = 1, n_intp
              do jj = 1, ngfl
                 zeta( jj, ii, isp, ikp ) = Xmtrx( ii, jj )
                 ! if( verbose .and. ii.eq.1 .and. & 
                 !    jj <= 80 .and. mod(jj,8)==1) &
                 !    write(*,*) zeta(jj,ii,isp,ikp)
              enddo ! jj loop
           enddo ! ii loop
           !
        enddo ! ikp loop
     enddo ! isp loop
     if ( verbose ) then

        !
        ! calculate zeta * C, which should in principle be a good approximation of Z
        !
        ALLOCATE( tmp_Zmtrx(ngfl) )
        !
        do isp = 1, nspin
           do ikp = 1, kpt%nk
              !
              if( w_grp%master ) then
                 write( dbgunit, '(a,i2,a,i6)' ) " isp ", isp," ikp ", ikp
                 write( dbgunit, '(a)' ) "   ivv    icc      diff "
              endif
              !
              do icv = 1, ncv(isp)
                 !
                 ! Calculate: tmp_Zmtrx(ngfl) = 
                 !   sum_{n_intp} zeta(ngfl,n_intp) * Cmtrx(n_intp,icv)
                 !
                 call dgemv('n', ngfl, n_intp, one, &
                    zeta(1,1,isp,ikp), ngfl, &
                    Cmtrx(1,icv,isp,ikp), 1, zero, tmp_Zmtrx, 1)
                 !
                 ! calculate the relative difference between
                 ! the interpolated |vc> and the real |vc>
                 !
                 diff = zero
                 weight = zero
                 ivv = invpairmap(1,icv,isp,ikp)
                 icc = invpairmap(2,icv,isp,ikp)
                 do igrid = 1, ngfl
                    !
                    diff = diff + abs(Zmtrx(igrid,icv,isp,ikp)-tmp_Zmtrx(igrid))
                    weight = weight + abs(Zmtrx(igrid,icv,isp,ikp))
                    !if ( (igrid < 160 .or. (igrid <9500 .and. igrid > 9000)) .and. &
                    !    mod(igrid,8) .eq. 1 ) & 
                    !  write(dbgunit,'(i10,2x,3e15.5)') &
                    !    igrid, &
                    !    Zmtrx(igrid,icv,isp,ikp), &
                    !    tmp_Zmtrx(igrid,icv), &
                    !    Zmtrx(igrid,icv,isp,ikp)-tmp_Zmtrx(igrid,icv)  ! For Debug
                    !
                 enddo
                 call MPI_ALLREDUCE(MPI_IN_PLACE, diff, 1, MPI_DOUBLE, MPI_SUM, & 
                   w_grp%comm, errinfo)
                 call MPI_ALLREDUCE(MPI_IN_PLACE, weight, 1, MPI_DOUBLE, MPI_SUM, & 
                   w_grp%comm, errinfo)
                 diff = diff / weight
                 if( w_grp%master ) write(dbgunit, '(2i6,e15.5)') ivv, icc, diff
              enddo ! icv loop
           enddo ! ikp loop
        enddo ! isp loop
        !
        DEALLOCATE(tmp_Zmtrx)
        !
     endif ! if ( outputdbg ) 
     ! 
     !
     ! clean up all the allocatable variables
     !
     if ( w_grp%master .and. verbose ) then
        write(*,*) " DEALLOCATE arrays"
     endif
     DEALLOCATE( Amtrx )
     DEALLOCATE( Bmtrx )
     DEALLOCATE( Xmtrx )
     DEALLOCATE( Zmtrx )
     !
     call MPI_Barrier( w_grp%comm, errinfo )
     !
     ! Now calculate <zeta_u(r,ispin)|V(r,r')|zeta_w(r',jspin)>, where u, w = 1, ..., n_intp
     ! and store it in Mmtrx(n_intp, n_intp)
     !
     Mmtrx = zero ! Mmtrx dimension: Mmtrx(n_intp, n_intp, nspin, nspin)
     !
     ! only workrgrp, workwgrp calculate Mmtrx, then peinf%master send Mmtrx to
     ! all the other procs  
     !
     ! qkt is set to zero. This is only valid for
     !  tests of nonperiodic system, and should be updated later.
     !
     qkt = 0
     ! ngf is the number of grid points in the full wedge
     ALLOCATE(rho_h(ngf)) ! note: gvec%nr is equal to w_grp%nr
     ! ngfl is the number of grid points distributed in each proc, 
     ! ngfl = w_grp%mydim*gvec%syms%ntrans
     ALLOCATE(rho_h_local(ngfl))
     ! if ( kflag > 0 ) ALLOCATE(rho_h_localsp(ngfl, nspin))
     !
     ! kflag = 0 : calculate kernel K^x  ( Coulomb )
     !         1 :                  K^x + K^f ( Coulomb + 1/2 * F_xc )
     !         2 :                  K^f   ( 1/2 * F_xc )
     !         3 :                  K^t   ( 1/2 * F_xc for spin triplet )
     !
     ! Generate Coulomb potential
     !
     if (kflag < 2 ) then
        call dinitialize_FFT(peinf%inode, fft_box)
        call dcreate_coul_0D(gvec%bdot, qkt, fft_box)
     endif
     !
     ! Calculate LDA kernel
     !
     if ( kflag > 0 ) then
        !
        ALLOCATE( fxc_loc( w_grp%mydim, nspin, nspin ), stat = errinfo )
        ALLOCATE( fxc( ngfl, nspin, nspin ), stat = errinfo )
        ALLOCATE( fzeta(ngfl, nspin), stat = errinfo )
        fxc_loc = zero
        do isp = 1, nspin
           call dcopy( w_grp%mydim, kpt%rho(w_grp%offset+1, isp), 1, &
              fxc_loc(1, isp, 1), 1 )
        enddo 
        call xc_init( nspin, XC_LDA_X, XC_LDA_C_PZ, 0, zero, one, .false., xc_lda )
        xc_lda%has_grad = .false.
        call fxc_get( xc_lda, nspin, w_grp%mydim, kflag, fxc_loc )
        call xc_end( xc_lda )
        !
        igrid = 0
        do ipt = 1, w_grp%mydim
          do itrans = 1, gvec%syms%ntrans
            igrid = igrid + 1
            fxc(igrid, :, :) = fxc_loc(ipt, :, :)
          enddo
        enddo
        DEALLOCATE( fxc_loc )
        !
     endif
     !
     do ikp = 1, kpt%nk
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
           !
           ! initialize rho_h and rho_h_local
           rho_h = 0.d0
           rho_h_local = 0.d0
           !
           if ( kflag < 2 ) then
              ioff1 = ( offset(w_grp%inode)-1 ) * gvec%syms%ntrans
              ! ngfl = w_grp%mydim * gvec%syms%ntrans
              call dcopy( ngfl, &
                 zeta( 1, ii, isp1, ikp ), 1, &
                 rho_h( ioff1 + 1 ), 1 )
              call MPI_ALLREDUCE( MPI_IN_PLACE, rho_h(1), ngf, MPI_DOUBLE, &
                 MPI_SUM, w_grp%comm, errinfo )
              !
              ! print out some debug info
              !
              !if ( w_grp%master .and. ii == 1 .and. verbose ) then
              !   !
              !   write(dbgunit,'("rho_h(:)")')
              !   do igrid = 1, ngf
              !      write(dbgunit,'(e15.5)',advance='no') rho_h(igrid)
              !      if(mod(igrid,8)==0) write(dbgunit,'()')
              !   enddo
              !   !
              !endif
              !
              ! solve poisson equation to get: rho_h(r) = \int V_c(r,r') zeta_ii(r') dr'
              !
              call dfullpoisson(gvec, rho_h, dbgunit,.FALSE.)
              !
              ! print out some debug info
              !
              !if ( w_grp%master .and. ii == 1 .and. verbose ) then
              !   write(dbgunit,'("After poisson solver")')
              !   write(dbgunit,'("rho_h(:)")')
              !   do igrid = 1, ngf
              !      write(dbgunit,'(e15.5)',advance='no') rho_h(igrid)
              !      if(mod(igrid,8)==0) write(dbgunit,'()')
              !   enddo
              !endif
              ! Put rho_h to rho_h_local
              ioff1 = ( offset(w_grp%inode)-1 ) * gvec%syms%ntrans
              ! ngfl = w_grp%mydim * gvec%syms%ntrans
              ! rho_h_local = f_x * zeta
              call dcopy( ngfl, rho_h( ioff1 + 1 ), 1, rho_h_local( 1 ), 1)
              !
           endif ! kflag < 2 
           !
           if ( kflag > 0 ) then
              !
              ! rho_h_localsp = 0.d0
              do isp = 1, nspin
                 call dcopy( ngfl, zeta( 1, ii, isp1, ikp ), 1, &
                   fzeta(1, isp), 1 )
                 ! calculate fzeta = f_lda * zeta
                 call dmultiply_vec( ngfl, fxc(1,isp1,isp), fzeta(1, isp) ) 
                 ! calculate rho_h_localsp = rho_h_local + f_lda * zeta
                 ! call dcopy( ngfl, rho_h_local(1), 1, rho_h_localsp(1,isp), 1 )
                 ! call daxpy( ngfl, one, fzeta, 1, rho_h_localsp(1,isp), 1 ) 
              enddo 
              !
           endif ! kflag > 0
           !
           ! Mmtrx is a symmetric matrix
           ! So we start the do loop from "ipair1" to "n_intp * spin"
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
              if ( kflag > 0 ) then
                Mmtrx( ii, jj, isp1, isp2, ikp, 2 ) = &
                ! ddot( ngfl, zeta( 1, jj, isp2, ikp ), 1, rho_h_localsp(1,isp2), 1 ) / gvec%hcub
                  ddot( ngfl, zeta( 1, jj, isp2, ikp ), 1, fzeta(1,isp2), 1 ) / gvec%hcub
                !
                ! Mmtrx is a symmetric matrix
                !
                Mmtrx( jj, ii, isp2, isp1, ikp, 2 ) = &
                  Mmtrx( ii, jj, isp1, isp2, ikp, 2)
              endif
              if ( kflag < 2 ) then
                Mmtrx( ii, jj, isp1, isp2, ikp, 1 ) = &
                  ddot( ngfl, zeta( 1, jj, isp2, ikp ), 1, rho_h_local(1), 1 ) / gvec%hcub
                !
                ! Mmtrx is a symmetric matrix
                !
                Mmtrx( jj, ii, isp2, isp1, ikp, 1 ) = &
                  Mmtrx( ii, jj, isp1, isp2, ikp, 1 )
              endif

              !
           enddo ! ipair2
           !
        enddo ! ipair1
        !
     enddo ! ikp
     !
     if ( kflag > 0 ) then
       DEALLOCATE(fxc)
       DEALLOCATE(fzeta)
       ! DEALLOCATE(rho_h_localsp)
     endif
     DEALLOCATE(zeta) ! no longer needed
     DEALLOCATE(rho_h)
     DEALLOCATE(rho_h_local)
  else
     Mmtrx = zero
  endif ! workrgrp == r_grp%mygr .and. workwgrp == w_grp%mygr
  !
  if ( kflag > 0 ) then
   call MPI_ALLREDUCE( MPI_IN_PLACE, Mmtrx(1,1,1,1,1,2), n_intp * n_intp * nspin * nspin * kpt%nk, &
     MPI_DOUBLE, MPI_SUM, w_grp%comm, errinfo )
  endif
  if ( kflag < 2 ) then
   call MPI_ALLREDUCE( MPI_IN_PLACE, Mmtrx(1,1,1,1,1,1), n_intp * n_intp * nspin * nspin * kpt%nk, &
     MPI_DOUBLE, MPI_SUM, w_grp%comm, errinfo )
  endif
  !
  if ( workrgrp == r_grp%mygr .and. workwgrp == w_grp%mygr ) then
     if ( w_grp%master .and. verbose ) then
        write( *, * ) " ikp = ", ikp, " Mmtrx (:, :, isp=1, isp=1, ikp=1, 1) = "
        call printmatrix ( Mmtrx (1,1,1,1,1,1), n_intp, n_intp, dbgunit ) 
        close ( dbgunit )
     endif
  endif
  return
  !
end subroutine isdf_parallel
