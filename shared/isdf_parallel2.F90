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
! __________________________________________________________________________________
! |  Matrix      |     Shape      |        Contents                                |
! |______________|________________|________________________________________________|
! |   P(r,r_u)   |  Ng * N_u      |    \sum_i \phi_i(r)\psi_i(r_u)                 |
! |              |                |     i=1...Nv, r are full-grid points           | 
! |______________|________________|________________________________________________|
! |   Q(r,r_u)   |  Ng * N_u      |    \sum_i \phi_i(r)\psi_i(r_u)                 |
! |              |                |     i=1...Nc, r are full-grid points           | 
! |______________|________________|________________________________________________|
! |  Zeta_u(r)   |  Ng *  Nu      |    {..., \zeta_u(r), ...}                      |
! |              |                |     u=1...Nu                                   |
! |______________|________________|________________________________________________|
! |  C(r_u,i,j)  |  Nu * (Nc*Nv)  |    {..., \phi_i(ru)*\psi_j(ru), ...}           |
! |              |                |                                                |
! |______________|________________|________________________________________________|
! _________________________________________________________________________
! |  Matrix      |                  Parallelization                       |
! |______________|________________________________________________________|
! |   P(r,r_u)   |             Each proc store part of Z                  |
! |              |             P(ngfl, Nc*Nv), ngfl = mydim               |
! |______________|________________________________________________________|
! |   Q(r,r_u)   |             Each proc store part of Z                  |
! |              |             Q(ngfl, Nc*Nv), ngfl = mydim               |
! |______________|________________________________________________________|
! |  Zeta_u(r)   |             Each proc store part of zeta               |
! |              |         zeta(ngfl, Nu), ngfl = mydim*ntrans            |
! |______________|________________________________________________________|
! |  C(r_u,i,j)  |             Every proc store a full copy of C          |
! |              |              Cmtrx ( Nu, Nc*Nv )                       |
! |______________|________________________________________________________|
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
! This subroutine consists of the following main steps:
! 
! Step 1.0: Prepare P(r,r_u)
! Step 1.5: Prepare A=C.C^T=P(r_u,r_v)Q(r_u,r_v), B=Z.C^T=P(r,r_u)Q(r,r_v)
! Step 1.6: Deallocate P and Q
! Step 2.0: solve a linear quation A.zeta=B to get zeta
! Step 3.0: Calculate Mmtrx(1) = <zeta(i)| V_coul(r,r') | zeta(j)> and
!       Mmtrx(2) = <zeta(i)| f^LDA(r) | zeta(j)>
! Step 3.5: Calculate C(r_u,i,j)

subroutine isdf_parallel2( gvec, pol_in, kpt, n_intp, intp, &
      nspin, ncv, maxncv, invpairmap, nv, maxnv, ivlist, nc, maxnc, iclist, kflag, &
      Cmtrx, Mmtrx, verbose )
 
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
  integer, intent(in) :: n_intp, intp(n_intp), nspin, ncv(2), maxncv, kflag, &
    maxnv, maxnc
  real(dp), intent(out) :: Cmtrx(n_intp, maxncv, nspin, kpt%nk)  
  real(dp), intent(out) :: Mmtrx(n_intp, n_intp, nspin, nspin, kpt%nk, 2)
  ! invpairmap maps the index of pair |cv> to c and v, which are the
  ! index of wave functions
  integer, intent(in) :: invpairmap(2, maxncv, nspin, kpt%nk)
  integer, intent(in) :: ivlist(maxnv, nspin), iclist(maxnc, nspin), nv(nspin), nc(nspin)
  ! If verbose is true, then print out additional debug information
  logical, intent(in) :: verbose   

  ! --- Local variables ---
  ! P(gvec%mydim, Nu, nspin, kpt%nk)
  ! Q(gvec%mydim, Nu, nspin, kpt%nk)
  ! Cmtrx(n_intp, Nc*Nv, nspin, kpt%nk)
  ! Issue need to be addressed: if there are N wfn_grp, does it mean these matrix need to be
  ! calculated by all processors at the same time, or can we distribute the workload
  ! later ??
  real(dp), allocatable ::      &
   PsiV(:,:), PsiV_intp(:,:),   &
   PsiC(:,:), PsiC_intp(:,:),   &
   P(:,:), P_intp(:,:),         &   ! P on irreducible grid
   Q(:,:), Q_intp(:,:),         &   ! Q on irreducible grid
   zeta(:,:,:,:), tmp_Zmtrx(:), Zmtrx(:), Zmtrx_loc(:), &
   fxc(:,:,:), fxc_loc(:,:,:), fzeta(:,:), &
   ! matrices and vectors used for solving linear equations
   Amtrx(:,:), Bmtrx(:,:), Xmtrx(:,:), tmpmtrx(:,:,:,:), &
   rho_h(:), rho_h_local(:)
  real(dp) :: diff, weight, qkt(3), tsec(2), &
   vcharac(gvec%syms%ntrans), ccharac(gvec%syms%ntrans)
  integer, allocatable :: inv_ivlist(:,:), inv_iclist(:,:)
  ! counters and temporary integers
  integer :: ipt, ii, jj, iv, ic, icv, & 
             IVV, ICC, JVV, JCC, itrans, isp, ikp, errinfo, igrid, &
             iproc, tag, virep, cirep, ipair1, ipair2, isp1, isp2, maxivv, maxicc
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
     call timacc(53,1,tsec)
     if( w_grp%master .and. verbose ) then
        !
        write(*,*) "call isdf_parallel(), write debug info to ", dbg_filename
        !
        open(dbgunit, file=dbg_filename, form='formatted', status='replace')
        !
     endif
     !
     ! the number of real-space grid points in reduced real-space domain
     ngr = gvec%nr 
     ! the number of grid points in full real-space domain, gvec%syms%ntrans is the number of sym ops
     ngf = gvec%nr * gvec%syms%ntrans 
     ! the number of grid points of interpolation vectors zeta(:) stored in each processor
     ngfl = w_grp%mydim * gvec%syms%ntrans
     !
     ! each processor stores part of P, Q and zeta
     !
     ! intermediate variables
     ALLOCATE( P( w_grp%mydim, n_intp ) )
     ALLOCATE( P_intp( n_intp, n_intp ) )
     ALLOCATE( Q( w_grp%mydim, n_intp ) )
     ALLOCATE( Q_intp( n_intp, n_intp ) )
     ALLOCATE( Amtrx( n_intp, n_intp ) )
     ALLOCATE( Bmtrx( n_intp, w_grp%mydim ) )
     ALLOCATE( Xmtrx( n_intp, w_grp%mydim ) )
     ! output quantity
     ALLOCATE( zeta ( ngfl, n_intp, nspin, kpt%nk ) ) ! interpolation vectors
     !
     ! initialize matrices with zero
     !
     Cmtrx = zero
     zeta  = zero
     !
     if ( w_grp%master .and. verbose ) then
        write(dbgunit, '(a)') " Index of interpolation points in full domain: "
        write(dbgunit, '(5i15)') ( intp(ii), ii=1, n_intp )
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
     if ( verbose .and. w_grp%master ) then
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
     maxivv = maxval(ivlist)
     maxicc = maxval(iclist)
     ! ivlist(:) maps the index of valence states used in 
     !   calculation (i.e., stored in memory) to the real index of valence states
     allocate(inv_ivlist(maxivv, nspin))
     allocate(inv_iclist(maxicc, nspin))
     ! inv_ivlist(:) maps the real index of valence states 
     !   to the index of valence states used in the calculation
     inv_ivlist = 0
     inv_iclist = 0
     !
     do isp = 1, nspin
        !
        do iv = 1, nv(isp)
           IVV = ivlist(iv,isp)
           inv_ivlist(IVV,isp) = iv
        enddo
        !
        do ic = 1, nc(isp)
           ICC = iclist(ic,isp)
           inv_iclist(ICC,isp) = ic
        enddo
        !
        allocate(PsiV( ngfl, nv(isp) ))
        allocate(PsiV_intp( n_intp, nv(isp) ))
        allocate(PsiC( ngfl, nc(isp) ))
        allocate(PsiC_intp( n_intp, nc(isp) ))
        !
        do ikp = 1, kpt%nk
           !
           ! initialize matrices with zero
           !
           P = zero 
           Q = zero
           Amtrx = zero
           Bmtrx = zero
           Xmtrx = zero
           PsiV  = zero
           PsiC  = zero
           PsiV_intp = zero
           PsiC_intp = zero
           !
           if (verbose .and. w_grp%master) &
              write(dbgunit, *) ' isp = ', isp, ', ikp = ', ikp
           !
           do iv = 1, nv(isp)
              IVV = ivlist(iv,isp)
              JVV = kpt%wfn(isp,ikp)%map(IVV)
              virep = kpt%wfn(isp,ikp)%irep(JVV)
              vcharac = gvec%syms%chi(virep,:)
              ! PsiV_i(r)
              do itrans = 1, gvec%syms%ntrans
                 call daxpy(w_grp%mydim, vcharac(itrans), kpt%wfn(isp,ikp)%dwf(1,JVV), 1, PsiV(itrans,IV), gvec%syms%ntrans)
              enddo
           enddo
           !
           do ic = 1, nc(isp)
              ICC = iclist(ic,isp)
              JCC = kpt%wfn(isp,ikp)%map(ICC)
              cirep = kpt%wfn(isp,ikp)%irep(JCC)
              ccharac = gvec%syms%chi(cirep,:)
              ! PsiC_i(r)
              do itrans = 1, gvec%syms%ntrans
                 call daxpy(w_grp%mydim, ccharac(itrans), kpt%wfn(isp,ikp)%dwf(1,JCC), 1, PsiC(itrans,IC), gvec%syms%ntrans)
              enddo
           enddo
           !
           do ipt = 1, n_intp 
              iptf = intp(ipt)
              ioff1 = ( offset(w_grp%inode)-1 ) * gvec%syms%ntrans
              if ( w_grp%inode + 1 < w_grp%npes ) then
                ioff2 = ( offset(w_grp%inode+1)-1 ) * gvec%syms%ntrans
              else
                ioff2 = ngf
              endif
              if ( iptf .ge. ioff1 + 1 &
                   .and. &
                   iptf .lt. ioff2 + 1 ) then
                jj = iptf - ioff1 
                ! if ( icv .eq. 1 .and. verbose ) write(*, '(4(a,i5))') &
                !    ", wgrp%inode ", w_grp%inode, ", ipt ", ipt, &
                !    ", iptf ", iptf, ", jj ", jj
                PsiV_intp(ipt, :) = PsiV(jj, :)
                PsiC_intp(ipt, :) = PsiC(jj, :)
              endif
           enddo
           !
           call MPI_ALLREDUCE(MPI_IN_PLACE, PsiV_intp, n_intp*nv(isp), MPI_DOUBLE, MPI_SUM, &
             w_grp%comm, errinfo)
           call MPI_ALLREDUCE(MPI_IN_PLACE, PsiC_intp, n_intp*nc(isp), MPI_DOUBLE, MPI_SUM, &
             w_grp%comm, errinfo)
           !
           ! Prepare Cmtrx
           !
           do icv = 1, ncv(isp)
              IVV = invpairmap(1,icv,isp,ikp)
              ICC = invpairmap(2,icv,isp,ikp)
              Cmtrx(1:n_intp, icv, isp, ikp) = PsiV_intp(1:n_intp, inv_ivlist(IVV,isp)) * &
                PsiC_intp(1:n_intp, inv_iclist(ICC,isp)) ! This is element-wise multiplication
           enddo
           !
           ! To reduce the memory cost, here we only store w_grp%mydim points in P and Q, instead of
           !  w_grp%mydim*syms%ntrans points
           !
           do itrans = 0, gvec%syms%ntrans-1 ! A loop over all symm ops to cover all ngfl points, with w_grp%mydim points each times
              !
              ! ( Note: PsiV and PsiC are both Kohn-Sham wavefunctions,
              !         There are no difference between them except for the band index.
              !         In TDLDA (not in SIGMA calculations), PsiV are valence states, and PsiC are empty states.
              !         Here I use different symbols to distinguish them )
              ! Calculate P(r,r_u) = \sum^Nv_i \Psi_i(r) \Psi_i(r_u) 
              !   PsiV(w_grp%mydim, nv) * PsiV_intp^T(nv, n_intp)
              call dgemm('n', 't', w_grp%mydim, n_intp, nv(isp), one, &
                PsiV(itrans*w_grp%mydim+1,1), ngfl, PsiV_intp(1,1), n_intp, zero, P(1,1), w_grp%mydim)
              !
              ! Calculate Q(r,r_u) = \sum^Nc_i \Psi_i(r) \Psi_i(r_u)
              call dgemm('n', 't', w_grp%mydim, n_intp, nc(isp), one, &
                PsiC(itrans*w_grp%mydim+1,1), ngfl, PsiC_intp(1,1), n_intp, zero, Q(1,1), w_grp%mydim)
              !
              ! Calculate P(r_u, r_u) = \sum^Nv_i \PsiV_i(r_u) \PsiV_i(r_u) and
              !      Q(r_u, r_u) = \sum^Nc_i \Psi_i(r_u) \Psi_i(r_u) by selection
              !
              call dgemm('n', 't', n_intp, n_intp, nv(isp), one, &
                PsiV_intp(1,1), n_intp, PsiV_intp(1,1), n_intp, zero, P_intp(1,1), n_intp)
              !
              call dgemm('n', 't', n_intp, n_intp, nc(isp), one, &
                PsiC_intp(1,1), n_intp, PsiC_intp(1,1), n_intp, zero, Q_intp(1,1), n_intp)
              if ( verbose .and. w_grp%master ) then
                 write(dbgunit, '(" PsiV_intp = ")')
                 call printmatrix ( PsiV_intp(1,1), n_intp, nv(isp), dbgunit )
              endif
              if ( verbose .and. w_grp%master ) then
                 write(dbgunit, '(" PsiC_intp = ")')
                 call printmatrix ( PsiC_intp(1,1), n_intp, nc(isp), dbgunit )
              endif
              !
              ! ---------------------
              ! For each spin and ikp, calculate zeta
              ! zeta dimension: Ng*Nu
              ! Set A = (P_intp(r_u,r_u).Q_intp(r_u,r_u))^T  dimension: n_intp * n_intp
              !     B = (P(r,r_u).Q(r,r_u))^T                dimension: n_intp * w_grp%mydim
              ! A is a symmetric matrix!
              ! we solve the linear equation A*(zeta^T) = B to get zeta^T
              ! zeta^T = A^-1 * B = (C * C^T)^-1 * C * Z^T
              !
              !  Matrix dimensions:
              !  zeta(ngfl, n_intp,:,:)
              !  Amtrx(n_intp,n_intp)       intermediate variable, store C * C^T
              !  Bmtrx(n_intp,w_grp%mydim)  intermediate variable, store C * Z^T
              !  Xmtrx(n_intp,w_grp%mydim)         intermediate variable, store zeta^T
              ! ---------------------
              !
              ! calculate A = P_intp.Q_intp (Note: This is an element-wise multipliation)
              !
              Amtrx = P_intp*Q_intp 
              !
              if ( verbose .and. w_grp%master ) then
                 write(dbgunit,'(" A = C*C^T = P_intp*Q_intp = ")')
                 call printmatrix ( Amtrx(1,1), n_intp, n_intp, dbgunit )
              endif
              !
              ! calculate B = (P.Q)^T (Note: This is an element-wise multiplication)
              Bmtrx = transpose(P*Q)
              !
              if ( verbose .and. w_grp%master ) then
                 write(dbgunit,'(" B = C*Z^T = P*Q = ")')
                 call printmatrix ( Bmtrx(1,1), n_intp, w_grp%mydim, dbgunit )
              endif
              !
              ! solver the linear equation A * X = B
              !
              call dlinear_solver( n_intp, w_grp%mydim, Amtrx, Bmtrx, Xmtrx, w_grp%inode, verbose )
              !
              ! Copy Xmtrx to zeta
              !
              ! if ( verbose ) write(*,*) 'zeta', isp,ikp
              do ii = 1, n_intp
                 do jj = 1, w_grp%mydim
                    zeta( itrans*w_grp%mydim+jj, ii, isp, ikp ) = Xmtrx( ii, jj )
                    ! if( verbose .and. ii.eq.1 .and. & 
                    !    jj <= 80 .and. mod(jj,8)==1) &
                    !    write(*,*) zeta(jj,ii,isp,ikp)
                 enddo ! jj loop
              enddo ! ii loop
           enddo ! itrans
           !
        enddo ! ikp loop
        DEALLOCATE( PsiV )
        DEALLOCATE( PsiV_intp )
        DEALLOCATE( PsiC )
        DEALLOCATE( PsiC_intp )
     enddo ! isp loop
     call timacc(53,2,tsec)
     !
     ! clean up all the allocated variables
     !
     if ( w_grp%master .and. verbose ) then
        write(*,*) " DEALLOCATING arrays"
     endif
     DEALLOCATE( Amtrx )
     DEALLOCATE( Bmtrx )
     DEALLOCATE( Xmtrx )
     DEALLOCATE( P )
     DEALLOCATE( P_intp )
     DEALLOCATE( Q )
     DEALLOCATE( Q_intp )
     DEALLOCATE( inv_ivlist )
     DEALLOCATE( inv_iclist )
     if ( verbose ) then
        !
        ! calculate zeta * C, which should in principle be a good approximation of Z
        !
        ALLOCATE( tmp_Zmtrx(ngfl) )
        ALLOCATE( Zmtrx(ngfl) )
        ALLOCATE( Zmtrx_loc(w_grp%mydim) )
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
                 ! Calculate the exact Zmtrx from wfns
                 !
                 IVV = invpairmap(1,icv,isp,ikp)
                 ICC = invpairmap(2,icv,isp,ikp)
                 if (verbose .and. w_grp%master) &
                    write(dbgunit, '(a,i5,a,i5)') " ivv = ", ivv, " icc = ", icc
                 JVV = kpt%wfn(isp,ikp)%map(IVV)
                 JCC = kpt%wfn(isp,ikp)%map(ICC)
                 !
                 ! copy phi_v to Zmtrx
                 !
                 call dcopy(w_grp%mydim, &
                   kpt%wfn(isp,ikp)%dwf(1,JVV),1,Zmtrx_loc(1),1)
                 !
                 ! calculate phi_v(r)*phi_c(r) and store it in Zmtrx
                 !
                 call dmultiply_vec(w_grp%mydim, &
                   kpt%wfn(isp,ikp)%dwf(1,JCC),Zmtrx_loc(1))
                 !
                 ! Get the symm group representation and characters of wfn |v> and |c>
                 !
                 virep = kpt%wfn(isp,ikp)%irep(JVV)
                 cirep = kpt%wfn(isp,ikp)%irep(JCC)
                 vcharac = gvec%syms%chi(virep,:)
                 ccharac = gvec%syms%chi(cirep,:)
                 !
                 ! Zmtrx_loc(1:w_grp%mydim) store the data of grid points in reduced domain, now unfold it
                 ! to full domain
                 !
                 igrid = 0
                 do ipt = 1, w_grp%mydim
                   do itrans = 1, gvec%syms%ntrans
                     igrid = igrid + 1
                     Zmtrx(igrid) = Zmtrx_loc(ipt) * &
                       vcharac(itrans) * ccharac(itrans)
                   enddo
                 enddo
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
                 do igrid = 1, ngfl
                    !
                    diff = diff + abs( Zmtrx(igrid)-tmp_Zmtrx(igrid) )
                    weight = weight + abs( Zmtrx(igrid) )
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
        DEALLOCATE(Zmtrx)
        DEALLOCATE(Zmtrx_loc)
        !
     endif ! if ( outputdbg ) 
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
        ! 
        ! Copy the charge density to fxc_loc
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
        ! Unfold fxc_loc (stored in irreducible domain) with available symm ops.
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
     call timacc(54,1,tsec)
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
     call timacc(54,2,tsec)
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
end subroutine isdf_parallel2
