!===================================================================
!
! Calculates the self-energy for an electronic system using 
! input matrix elements of the Coulomb potential and
! TDLDA exchange-correlation potential in the basis of one-electron
! transitions. Coulomb/exchange-correlation matrix elements are used 
! to set an eigenvalue equation for the full polarizability function, 
! equivalent to the TDLDA method.
!
! After eigenvalues/eigenvectors of the polarizability are calculated
! (subroutine diag_pol), the various terms included in the self-energy
! are calculated: 
!  sig_x = exchange (energy-independent)
!  sig_c = correlation, Sig_c = 2 * V_coul * Pol * V_coul
!  sig_g = TDLDA vertex,
!     Sig_g = V_coul * Pol * f_lda + f_lda * Pol * V_coul
!
! In order to preserve Phi-integrability, the term sig_g should be
! always included if the polarizability Pol is calculated within TDLDA.
!
! This code follows the numerical methodology presented in:
!    M.L. Tiago and J.R. Chelikowsky, Phys. Rev. B 73, 205334 (2006)
! See also:
!    M.E. Casida, in "Recent Advances in Density-Functional Methods", ed. D.P. Chong (1995)
!    I. Vasiliev, S. Ogut and J.R. Chelikowsky, Phys. Rev. B 65, 115416 (2002)
!    W. G. Aulbur, L. Jonsson, and J.W. Wilkins, Solid State Physics series, vol. 54, p. 1 (2000)
!
! Copyright (C) 2009 Murilo Tiago, Univ. of Texas, Austin, TX, USA
! mtiago@ices.utexas.edu
!
! First version written by Murilo Tiago, Univ. of Minnesota, April 2004
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 1, or (at your option)
!    any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston MA  02110-1301 USA
!
!-------------------------------------------------------------------
program sigma

  use typedefs
  use mpi_module
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
  type (gspace) gvec
  type (kptinfo) :: kpt, kpt_sig
  type (qptinfo) :: qpt
  type (kernelinfo), dimension(:,:), allocatable :: k_c, k_p
  type (polinfo), dimension(2) :: pol_in
  type (polinfo), dimension(:,:), allocatable :: pol
  type (siginfo) :: sig_in
  type (siginfo), dimension(:,:), allocatable :: sig
  type (qpinfo), dimension(:,:,:), allocatable :: q_p

  character (len=40), allocatable :: routnam(:)

  logical :: nolda, tamm_d, snorm, writeqp, readvxc, &
       readocc, cohsex, nooffd, init_gr, hqp_sym, lstop
  integer :: ii, irp, iq, ik, nmap, nspin, nbuff, lcache, isp, it_scf, &
       chkpt_in, n_it, nr_buff, static_type, sig_en, nkpt, dft_code
  real(dp) :: tsec(2), mem1, mem2, mem3, xsum, xmax, rtmp, tdldacut, &
        max_sig, max_conv, xbuff, ecuts, qpmix, sig_cut
  logical, allocatable :: wmap(:)

  ! W Gao ISDF method
  integer :: info
  logical :: doisdf
  real(dp), allocatable :: zeta(:,:,:,:)
  real(dp), allocatable :: Cmtrx(:,:,:,:), Mmtrx(:,:,:,:,:,:)
  ! the number of val and cond states pairs. Since there could be two spins,
  ! ncv(1) correspond to spin up, and ncv(2) correspond to spin down.
  ! Assumption: for different kpts, the number of states with the same 
  !  spin are the same 
  integer :: n_intp, maxnj, maxni, i, j, icc, ij, maxnij, &
           ihomo, ikp, intp_type, isdf_type, kflag
  ! cvt.f90
  integer, allocatable :: intp(:), pairmap(:,:,:,:), invpairmap(:,:,:,:), &
           ni(:), nj(:), nij(:), iilist(:,:), ijlist(:,:), timerlist(:)

  ! WG debug
  integer :: outdbg
  character (len=20) :: dbg_filename

  !-------------------------------------------------------------------
  ! Initialization.
  !
  call header('SIGMA')
  ! W Gao open dbg files
  write(dbg_filename,"(i7)") peinf%inode
  outdbg = peinf%inode+198812
  dbg_filename = "kernel_dbg"//adjustl(dbg_filename)
  open(outdbg,file=dbg_filename,status='unknown',iostat=info) 
  !-------------------------------------------------------------------
  ! Read input parameters from rgwbs.in.
  !
  call input_g(pol_in,qpt,tdldacut,nbuff,lcache,w_grp%npes, &
       nolda,tamm_d,r_grp%num,dft_code,doisdf,n_intp,intp_type,isdf_type,.false.)
  write(6, *) "isdf_type = ", isdf_type
  call MPI_BARRIER(peinf%comm,info)
  call input_s(sig_in,kpt_sig,snorm,writeqp,readvxc,readocc,cohsex, &
       nooffd,hqp_sym,n_it,chkpt_in,static_type,sig_en,max_conv,xbuff,ecuts, &
       qpmix,sig_cut,.false.)
  call MPI_BARRIER(peinf%comm,info)

  !-------------------------------------------------------------------
  ! Determine the set of wavefunctions to read: if n-th wavefunction is
  ! needed, then wmap(n) = .true.; otherwise, wmap(n) = .false.
  !
  nmap = 0
  if (sig_in%xc /= XC_GW) then
     if (sig_in%nmap > 0) then
        nmap = maxval(sig_in%map)
        allocate(wmap(nmap))
        wmap = .false.
        do ii = 1, sig_in%nmap
           wmap(sig_in%map(ii)) = .true.
        enddo
     else
        allocate(wmap(1))
     endif
  else
     if (max(pol_in(1)%ncond,pol_in(1)%nval,pol_in(2)%ncond,pol_in(2)%nval, &
          sig_in%nmax_c,sig_in%nmap) > 0) then
        nmap = max(sig_in%nmax_c,sig_in%nmap)
        if (pol_in(1)%ncond > 0) nmap = max(nmap,maxval(pol_in(1)%cmap))
        if (pol_in(1)%nval > 0) nmap = max(nmap,maxval(pol_in(1)%vmap))
        if (pol_in(2)%ncond > 0) nmap = max(nmap,maxval(pol_in(2)%cmap))
        if (pol_in(2)%nval > 0) nmap = max(nmap,maxval(pol_in(2)%vmap))
        if (sig_in%nmap > 0) nmap = max(nmap,maxval(sig_in%map))
        allocate(wmap(nmap))
        wmap = .false.
        do isp = 1, 2
           do ii = 1, pol_in(isp)%ncond
              wmap(pol_in(isp)%cmap(ii)) = .true.
           enddo
           do ii = 1, pol_in(isp)%nval
              wmap(pol_in(isp)%vmap(ii)) = .true.
           enddo
        enddo
        do ii = 1, sig_in%nmap
           wmap(sig_in%map(ii)) = .true.
        enddo
        do ii = 1, sig_in%nmax_c
           wmap(ii) = .true.
        enddo
     else
        allocate(wmap(1))
     endif
     if (min(pol_in(1)%ncond,pol_in(1)%nval,pol_in(2)%ncond,pol_in(2)%nval, &
          sig_in%nmax_c,sig_in%nmap) < 0) then
        deallocate(wmap)
        allocate(wmap(1))
        nmap = 0
     endif
  endif

  !-------------------------------------------------------------------
  ! Read wave-function file.
  !
  init_gr = .true.
  if ( dft_code == PARATEC ) then
     call paratec_wfn(gvec,kpt,nmap,nspin,wmap,init_gr)
  else
     call parsec_wfn(gvec,kpt,nmap,nspin,wmap,init_gr)
  endif
  deallocate(wmap)

  if (doisdf) then
     ! --- prepare some inputs for the ISDF method ---
     ! W Gao find the index of highest occupied orbital
     ihomo = 1
     do isp = 1, nspin
        do ikp = 1, kpt%nk
           do ii = 1, kpt%wfn(isp,ikp)%nstate
              if ( kpt%wfn(isp,ikp)%occ0(ii) > tol_occ .and. &
                   ihomo < ii) then
                 ihomo = ii
              endif
           enddo ! ii loop
        enddo ! ikp loop
     enddo ! isp loop
     ! if n_intp can not be found in rgwbs.in or invalid (i.e., less than the
     ! number of occupied states), then set it to the default value
     if(n_intp .lt. ihomo) then 
        n_intp = int(2.0 * ihomo)
     endif
     allocate(intp(n_intp))
     ! --- find interpolation points for ISDF method ---
     call stopwatch(peinf%master, "before call cvt")
     call timacc(51,1,tsec)
     if(intp_type .eq. 1) then
        if (peinf%master) then
           write(*,*) " intp_type == 1"
           call cvt(gvec, kpt%rho, nspin, n_intp, intp)
        endif
     elseif(intp_type .eq. 2) then
        if (peinf%master) write(*,*) " intp_type == 2"
        call cvt_wfn(gvec, kpt%wfn, nspin, kpt%nk, n_intp, intp)
     else
        write(*,*) 'Type',intp_type,'method for finding interpolation points is',&
           ' not implememted so far. The default method will be used.'
     endif
     call timacc(51,2,tsec)
     call stopwatch(peinf%master, "after call cvt")
     call MPI_BARRIER(peinf%comm, info)
    
     if(peinf%master) write(*,*) " Finding interpolation points successfully. "
     ! broadcast intp to all processors
     call MPI_bcast(intp(1), n_intp, MPI_INTEGER, peinf%masterid, peinf%comm, info)
     
     allocate(ni(nspin))
     allocate(nj(nspin))
     allocate(nij(nspin))
     do isp = 1, nspin 
       ! For now, I assume all the valance states are included, this is not necessarily true in all cases
       ! We should change this later to save some computational costs.
       ni(isp) = maxval(pol_in(isp)%vmap(:)) 
       do ii = 1, sig_in%ndiag_s
          if(sig_in%diag(ii) > ni(isp)) ni(isp) = sig_in%diag(ii)
       enddo
       if(peinf%master) write(outdbg, '(a,i5,a,i5)') "isp", isp, " ni", ni
     enddo
     do isp = 1, nspin
       nj(isp) = max( maxval(pol_in(isp)%cmap(:)), sig_in%nmax_c )
     enddo
     maxni = maxval(ni(:))
     maxnj = maxval(nj(:))
     ! obtain nv(isp), ivlist(nv(isp),isp), nc(isp) and iclist(nc(isp),isp)
     allocate(iilist(maxni,nspin))
     allocate(ijlist(maxnj,nspin))
     do isp = 1, nspin
       do ii = 1, ni(isp)
         iilist(ii,isp) = ii
       enddo
       do ij = 1, nj(isp)
         ijlist(ij,isp) = ij
       enddo
     enddo 
     allocate(pairmap(maxni, maxnj, nspin, kpt%nk))
     pairmap = 0
     ! For now, we assume pairmap is the same for different kpt
     do isp = 1, nspin
        do ikp = 1, kpt%nk
           ij = 0
           do i = 1, ni(isp)
              if (any( sig_in%diag == i)) then
                 do j = 1, nj(isp) 
                    ! we need to make sure we don't deal
                    ! with duplicated pairs
                    ! i.e., phi_i(r) phi_j(r) = phi_j(r) phi_i(r)
                    if( i <= j) then 
                    ij = ij + 1
                    pairmap(i,j,isp,ikp) = ij
                    endif
                 enddo ! j loop
              else 
                 do j = 1, pol_in(isp)%ncond
                    if( i <= j) then
                    ij = ij + 1
                    icc = pol_in(isp)%cmap(j)
                    pairmap(i,icc,isp,ikp) = ij
                    endif
                 enddo ! j loop
              endif
              !if(i<=ihomo) then
              !   !do j = 1, pol_in(isp)%ncond
              !   !   ij = ij + 1
              !   !   icc = pol_in(isp)%cmap(j)
              !   !   pairmap(i, icc, isp, ikp) = ij
              !   !enddo
              !   do j = 1, nj(isp)
              !      if ( i <= j) then
              !      ij = ij + 1
              !      pairmap(i, j, isp, ikp) = ij
              !      endif
              !   enddo
              !else
              !   !do j = 1, nj(isp)
              !   !   if( i < j) then
              !   !   ij = ij + 1
              !   !   pairmap(i, j, isp, ikp) = ij
              !   !   endif
              !   !enddo
              !endif
           enddo ! i loop
           nij(isp) = ij
        enddo ! ikp loop
     enddo ! isp loop
     if(peinf%master) then
        write(outdbg,*) " isp    i    j    pairmap"
        do isp = 1, nspin
           do i = 1,maxni
              do j = 1,maxnj
                 write(outdbg, '(4I7)') isp, i, j, pairmap(i,j,isp,1)
              enddo
           enddo
        enddo
     endif
     maxnij = maxval(nij(:))
     allocate(Cmtrx(n_intp, maxnij, nspin, kpt%nk))
     allocate(Mmtrx(n_intp, n_intp, nspin, nspin, kpt%nk, 2))
     allocate(invpairmap(2,maxnij, nspin, kpt%nk))
     ! set up invpairmap
     do isp = 1, nspin
        do ikp = 1, kpt%nk
           ij = 0
           do i = 1, ni(isp)
              if (any( sig_in%diag == i)) then
                 do j = 1, nj(isp) 
                    if( i <= j) then 
                    ij = ij + 1
                    invpairmap(1,ij,isp,ikp) = i
                    invpairmap(2,ij,isp,ikp) = j
                    endif
                 enddo
              else 
                 do j = 1, pol_in(isp)%ncond
                    if( i <= j) then
                    ij = ij + 1
                    icc = pol_in(isp)%cmap(j)
                    invpairmap(1,ij,isp,ikp) = i
                    invpairmap(2,ij,isp,ikp) = icc
                    endif
                 enddo
              endif
              !if(i<=ihomo) then
              !   !do j = 1, pol_in(isp)%ncond
              !   !   ij = ij + 1
              !   !   icc = pol_in(isp)%cmap(j)
              !   !   invpairmap(1, ij, isp, ikp) = i
              !   !   invpairmap(2, ij, isp, ikp) = icc
              !   !enddo
              !   do j = 1, nj(isp)
              !      if( i <= j) then
              !         ij = ij + 1
              !         invpairmap(1, ij, isp, ikp) = i
              !         invpairmap(2, ij, isp, ikp) = j
              !      endif
              !   enddo
              !else
              !   !do j = 1, nj(isp)
              !   !   if( i < j) then
              !   !   ij = ij + 1
              !   !   invpairmap(1, ij, isp, ikp) = i
              !   !   invpairmap(2, ij, isp, ikp) = j
              !   !   endif
              !   !enddo
              !endif
           enddo ! i loop 
        enddo ! ikp loop
     enddo ! isp loop
     if(peinf%master) then
        write(outdbg,*) " isp    ij   invpairmap"
        do isp = 1, nspin
           do ij = 1,maxnij
              write(outdbg, '(4I7)') isp, ij, &
                invpairmap(1,ij,isp,1), invpairmap(2,ij,isp,1)
           enddo
        enddo
     endif
  endif ! if (doisdf)

  !-------------------------------------------------------------------
  ! Calculate characters of representations.
  !
  if (kpt%lcplx) then
     call zcharac_group(gvec%syms,gvec,kpt,70,nspin,kpt%wfn(1,1)%nstate)
  else
     call dcharac_group(gvec%syms,gvec,kpt,70,nspin,kpt%wfn(1,1)%nstate)
  endif
  !-------------------------------------------------------------------
  ! Initialize pol and k_p structures.
  !
  nkpt = kpt_sig%nk
  allocate(pol(gvec%syms%ntrans,qpt%nk))
  allocate(k_p(gvec%syms%ntrans,qpt%nk))
  allocate(k_c(gvec%syms%ntrans,qpt%nk))
  allocate(sig(nspin,nkpt))
  if (sig_in%xc == XC_GW) then
     pol(:,:)%lv = .false.
     pol(:,:)%ltv = .false.
     pol(:,:)%ntr = 0
     k_p(:,:)%ncol = 0
     k_p(:,:)%nbuff = nbuff
     k_p(:,:)%lcache = lcache
     k_p(:,:)%isdf = doisdf
     k_c(:,:)%ncol = 0
     k_c(:,:)%nbuff = nbuff
     k_c(:,:)%lcache = lcache
     k_c(:,:)%isdf = doisdf
     !
     ! Define number of grid points that fit in buffer size (xbuff is amount
     ! of space available on disk, measured in MB).
     ! If xbuff < 0 at input, then skip the calculation of static correction.
     !
     if (xbuff < 0) then
        nr_buff = 0
     elseif (xbuff == zero) then
        xsum = sum(kpt%rho(1:gvec%nr,1:nspin)) * two/real(nspin,dp)
        do nr_buff = 1, gvec%nr
           rtmp = sum(kpt%rho(1:nr_buff,1:nspin)) * two/real(nspin,dp)
           if (one - rtmp/xsum < 0.0001d0) exit
        enddo
     else
        nr_buff = min(one*gvec%nr,65536*xbuff/sum(pol(:,:)%ntr)/real(nspin,dp))
     endif
  endif

  call timacc(2,1,tsec)
  call stopwatch(peinf%master,'Calling setup_s')
  if (kpt%lcplx) then
     call zsetup_s(gvec,kpt,qpt,kpt_sig,sig_in,sig,pol_in,pol,k_p,k_c, &
          nspin,tdldacut,readocc,peinf%master)
  else
     call dsetup_s(gvec,kpt,qpt,kpt_sig,sig_in,sig,pol_in,pol,k_p,k_c, &
          nspin,tdldacut,readocc,peinf%master)
  endif
  call timacc(2,2,tsec)

  ! --- perform ISDF method to interpolate pair products of wave functions ---
  if(doisdf) then
    if (peinf%master) write(*,*) 'call isdf subroutine'
    call stopwatch(peinf%master, "before call isdf")
    kflag = 1
    if ( nolda ) kflag = 0
    call timacc(52,1,tsec)
    if ( isdf_type == 1 ) then
      write(6,*) " Call isdf_parallel()"
      call isdf_parallel(gvec, pol_in, kpt, n_intp, intp, &
        nspin, nij, maxnij, invpairmap, kflag, &
        Cmtrx, Mmtrx, .TRUE.) 
    else
      write(6,*) " Call isdf_parallel2()"
      call isdf_parallel2(gvec, pol_in, kpt, n_intp, intp, &
        nspin, nij, maxnij, invpairmap, ni, maxni, iilist, nj, maxnj, ijlist, kflag, &
        Cmtrx, Mmtrx, .TRUE.) 
    endif
    if(peinf%master) write(*,*) 'done isdf'
    call MPI_BARRIER(peinf%comm, info)
    call timacc(52,2,tsec)
    call stopwatch(peinf%master, "after call isdf")
  endif
  ! The outputs are Cmtrx and Mmtrx, which are used by k_integrate_isdf() for
  ! calculation of K(i,j,k,l) later !!
  ! Note: For wpol_v we also need zeta(:)!!!
  ! --- finished ISDF ---
  !
  ! --------------------------------------------------------------------------
  ! Calculate Vxc matrix elements and total energy.
  !
  if (kpt%lcplx) then
     do ik = 1, kpt_sig%nk
        do isp = 1, nspin
           call zvxc(nspin,gvec,kpt,sig(isp,ik),readvxc,ik,isp)
        enddo
     enddo
     call zetotal(gvec,kpt,nspin)
  else
     do ik = 1, kpt_sig%nk
        do isp = 1, nspin
           call dvxc(nspin,gvec,kpt,sig(isp,ik),readvxc,ik,isp)
        enddo
     enddo
     call detotal(gvec,kpt,nspin)
  endif

  do ik = 1, kpt%nk
     do isp = 1, nspin
        allocate(kpt%wfn(isp,ik)%cmapi(kpt%wfn(isp,ik)%nstate))
     enddo
  enddo

  !-------------------------------------------------------------------
  ! Print out warnings, information etc.
  !
  if (kpt%lcplx) tamm_d = .true.
  if (peinf%master) then
     write(6,'(/,a,/,/,a,/,2a,/)') repeat('-',65), &
          ' Self-energy input data: ', ' ',repeat('-',23)
     if (sig_in%xc == XC_GW) then
        write(6,'(2a)') ' Number of transitions per representation ', &
             'in TDLDA polarizabillity:'
        write(6,'(8i8)') ((pol(irp,iq)%ntr,irp=1,gvec%syms%ntrans),iq=1,qpt%nk)
        write(6,'(a,i10,/)') ' total = ', sum(pol(:,:)%ntr)
        if (tdldacut > zero) then
           write(6,*) 'Energy cutoff applied in TDLDA polarizability = ', &
                tdldacut*ryd, ' eV'
        else
           write(6,*) 'No energy cutoff in TDLDA polarizability'
        endif
        if (nolda) write(6,'(/,a,/)') &
             ' LDA kernel is not included in polarizability'
        if (tamm_d) then
           write(6,'(2a,/)') ' Calculating TDLDA ', &
                'polarizability within the Tamm-Dancoff approximation.'
        else
           write(6,'(2a,/)') ' Not using the Tamm-Dancoff ', &
                'approximation in TDLDA polarizability.'
        endif
        if (snorm) then
           write(6,'(a,/)') ' Renormalizing Sum rule '
        else
           write(6,'(a,/)') ' Sum rule not renormalized'
        endif
        write(6,'(2a,i5)') ' Order of highest LDA state included ', & 
             'in Green function = ', sig_in%nmax_c
        write(6,'(2a,i5,/)') ' Order of highest LDA state included ', & 
             'in COHSEX approximation = ', sig_in%nmap
        write(6,'(2a,g12.4,/)') ' Energy resolution in energy poles, ', &
             'self-energy (eV) = ', ecuts*ryd
        write(6,'(2a,f10.4,a,/)') ' Energy range used to calculate ', &
             'self-energy  = ', sig_in%deltae*ryd, ' eV '
        write(6,'(2a,i5,/)') ' Number of data points used to calculate ', &
             'self-energy = ', sig_in%nen
        select case(sig_en)
        case (SIG_LEFT)
           write(6,'(a,a,/)') ' Calculating < n_1 | Sigma | n_2 > at ', &
                'the energy of orbital n_1'
        case (SIG_RIGHT)
           write(6,'(a,a,/)') ' Calculating < n_1 | Sigma | n_2 > at ', &
                'the energy of orbital n_2'
        case (SIG_AV)
           write(6,'(a,a,/)') ' Calculating < n_1 | Sigma | n_2 > at ', &
                'the average of energies E_1, E_2.'
        end select
        if (nr_buff > 0) then
           xsum = sum(kpt%rho(1:gvec%nr,1:nspin)) * two/real(nspin,dp)
           rtmp = sum(kpt%rho(1:nr_buff,1:nspin)) * two/real(nspin,dp)
           write(6,*) 'Calculating W_pol in static limit with nr = ', nr_buff
           write(6,*) 'Fraction of electron density within nr = ', &
                rtmp/xsum * 1.d2, ' %'
           xsum = real(nr_buff,dp)*real(sum(pol(:,:)%ntr))/65536.d0
           write(6,'(a,f10.2,a,/)') ' Disk space used = ', xsum, ' MB'
        endif
        if (cohsex) write(6,'(a,a,/)') ' Using the ', &
             'COHSEX (static) approximation in self-energy.'
     else
        select case (sig_in%xc)
        case (XC_HF)
           write(6,'(2a,/)') ' Using the ', &
                'Hartree-Fock approximation in self-energy.'
        case (XC_B3LYP)
           write(6,'(2a,/)') ' Replacing ', &
                'self-energy with the hybrid B3LYP functional.'
        case (XC_LDA_CA)
           write(6,'(2a,/)') ' Replacing ', &
                'self-energy with the LDA CA-PZ functional.'
        case (XC_GGA_PBE)
           write(6,'(2a,/)') ' Replacing ', &
                'self-energy with the GGA PBE functional.'
        case (XC_GGA_BLYP)
           write(6,'(2a,/)') ' Replacing ', &
                'self-energy with the GGA BLYP functional.'
        end select
     endif
     if (qpmix /= one) write(6,'(a,/,a,g12.4,/)') &
          ' Mixing input and output (QP) wavefunctions.', &
          ' Mixing parameter = ', qpmix
     if (sig_cut > zero) write(6,'(a,f20.10,a,/)') &
          ' Applying cut-off ', sig_cut, ' eV to self-energy.'
     if (writeqp) write(6,'(a,/)') 'Printing out new wavefunctions.'
     if (n_it > 0) write(6,'(a,i5,a,/,a,f10.5,a,/)') &
          ' Performing ', n_it, ' SCGW iterations.', &
          ' Maximum converged potential = ', max_conv, ' eV'
     !
     ! Estimate memory usage.
     !
     mem1 = sum(kpt%wfn(:,:)%nmem)*two/real(nspin,dp) * &
          real(nspin*gvec%nr,dp)/two/131072.d0 / real(w_grp%npes,dp)
     if (kpt%lcplx) mem1 = mem1 * two
     write(6,'(a,f10.2,a)') ' Memory needed to store wavefunctions : ', &
           mem1,' MB/proc.'
     mem1 = real(5*gvec%nr,dp)/131072.d0
     if (kpt%lcplx) mem1 = mem1 * two
     write(6,'(a,a,f10.2,a)') ' Memory needed to calculate kernel ', &
          'matrix elements : ',mem1, ' MB/proc.'
     if (sig_in%xc == XC_GW) then
        ! For diagonalization, we store 4 matrices: hamiltonian/eigenvectors, 
        ! temporary array (in eigensolver), and kernel.
        ! Parallelized diagonalization also uses a second temporary matrix.
        xmax = 0
        do iq = 1, qpt%nk
           do irp = 1, gvec%syms%ntrans
              xsum = real(pol(irp,iq)%ntr * pol(irp,iq)%ntr,dp)
              if (xmax < xsum) xmax = xsum
           enddo
        enddo
        mem1 = xmax/1024.d0*3.d0/128.d0/r_grp%npes
        if (r_grp%npes > 1) mem1 = mem1*four / three
        ! For self-energy calculation, we store eigenvectors, 2 potential 
        ! matrices, and 2 kernel matrices.
        xmax = 0
        do iq = 1, qpt%nk
           do irp = 1, gvec%syms%ntrans
              xsum = real(pol(irp,iq)%ntr * pol(irp,iq)%ntr,dp) + &
                   real(k_c(irp,iq)%nrow * k_c(irp,iq)%ncol,dp) * two + &
                   real(pol(irp,iq)%ntr * k_c(irp,iq)%ncol,dp) + &
                   real(pol(irp,iq)%ntr * pol(irp,iq)%ntr,dp) /  real(r_grp%npes,dp)
              if ( xsum > xmax ) xmax = xsum
           enddo
           mem2 = xmax/1024.d0/128.d0/r_grp%npes
           xmax = 0
           do irp = 1, gvec%syms%ntrans
              xsum = real(pol(irp,iq)%ntr * pol(irp,iq)%ntr,dp) + &
                   real(k_c(irp,iq)%nrow * k_c(irp,iq)%ncol,dp) + &
                   real(pol(irp,iq)%ntr * k_c(irp,iq)%ncol,dp) * two + &
                   real(pol(irp,iq)%ntr * pol(irp,iq)%ntr,dp) / real(r_grp%npes,dp)
              if ( xsum > xmax ) xmax = xsum
           enddo
        enddo
        mem3 = real(xmax,dp)/1024.d0/128.d0/r_grp%npes
        write(6,'(a,f10.2,a)') ' Memory needed for diagonalization : ', &
             mem1, ' MB/proc.'
        write(6,'(a,f10.2,a)') ' Memory needed in self-energy calculation :', &
             max(mem2,mem3), ' MB/proc.'
        write(6,'(a,f10.2,a)') ' Minimum memory needed : ', max(mem1, &
             mem2,mem3), ' MB/proc.'
     endif
     write(6,'(/,a,/)') repeat('-',65)
     call flush(6)
  endif

  !-------------------------------------------------------------------
  ! Allocate structures.
  !
  do ik = 1, nkpt
     do isp = 1, nspin
        if (sig(isp,ik)%nmap > 0) then
           allocate(sig(isp,ik)%sigmai(sig(isp,ik)%nmap))
           sig(isp,ik)%sigmai = ecuts
        endif
        if (sig(isp,ik)%ndiag_s > 0) then
           allocate(sig(isp,ik)%xdiag(sig(isp,ik)%ndiag_s))
           allocate(sig(isp,ik)%scsdiag(sig(isp,ik)%ndiag_s))
        endif
        if (sig(isp,ik)%noffd_s > 0) then
           allocate(sig(isp,ik)%xoffd(sig(isp,ik)%noffd_s))
           allocate(sig(isp,ik)%scsoffd(sig(isp,ik)%noffd_s))
        endif
        if (sig_in%xc == XC_GW) then
           if (sig(isp,ik)%ndiag_s > 0) then
              allocate(sig(isp,ik)%sgsdiag(sig(isp,ik)%ndiag_s))
           endif
           if (sig(isp,ik)%ndiag > 0) then
              allocate(sig(isp,ik)%scdiag(sig_in%nen,2*sig(isp,ik)%ndiag))
              allocate(sig(isp,ik)%sexdiag(sig_in%nen,2*sig(isp,ik)%ndiag))
              allocate(sig(isp,ik)%sgdiag(sig_in%nen,2*sig(isp,ik)%ndiag))
           endif
           if (sig(isp,ik)%noffd_s > 0) then
              allocate(sig(isp,ik)%sgsoffd(sig(isp,ik)%noffd_s))
           endif
           if (sig(isp,ik)%noffd > 0) then
              allocate(sig(isp,ik)%scoffd(2,2*sig(isp,ik)%noffd))
              allocate(sig(isp,ik)%sgoffd(2,2*sig(isp,ik)%noffd))
           endif
        endif
     enddo
  enddo
  allocate(q_p(gvec%syms%ntrans,nspin,nkpt))
  !-------------------------------------------------------------------
  ! Start iterations.
  !
  do it_scf = 0, n_it

     if (kpt%lcplx) then
        call zcalculate_sigma(nspin,kpt_sig%nk,n_it,it_scf,nr_buff,static_type, &
             sig_en,dft_code,chkpt_in,gvec,kpt,qpt,k_c,k_p,pol,sig_in, &
             sig,q_p,nolda,tamm_d,writeqp,snorm,cohsex,hqp_sym,lstop, &
             ecuts,qpmix,sig_cut,max_sig,&
             Cmtrx, Mmtrx, n_intp, maxnj, maxni, maxnij, pairmap)
     else
        call dcalculate_sigma(nspin,kpt_sig%nk,n_it,it_scf,nr_buff,static_type, &
             sig_en,dft_code,chkpt_in,gvec,kpt,qpt,k_c,k_p,pol,sig_in, &
             sig,q_p,nolda,tamm_d,writeqp,snorm,cohsex,hqp_sym,lstop, &
             ecuts,qpmix,sig_cut,max_sig,&
             Cmtrx, Mmtrx, n_intp, maxnj, maxni, maxnij, pairmap)
     endif

     if (lstop) then
        call stopwatch(peinf%master,'WARNING: abort SCGW')
        exit
     endif

     if (sig_in%ndiag_s + sig_in%noffd_s == 0) goto 99
     if (max_sig < max_conv) exit

     if (peinf%master) inquire(file='stop_scgw',exist=lstop)
#ifdef MPI
     call MPI_BCAST(lstop,1,MPI_LOGICAL,peinf%masterid,peinf%comm,ii)
#endif
     if (lstop) then
        call stopwatch(peinf%master,'WARNING: abort SCGW')
        exit
     endif
  enddo

  ! Deallocate arrays for ISDF method
  if(doisdf) then
     deallocate(pairmap)
     deallocate(invpairmap)
     deallocate(Cmtrx)
     deallocate(Mmtrx)
     deallocate(ni)
     deallocate(nj)
     deallocate(nij)
  endif
  close(outdbg)

  if (kpt%lcplx) then
     if ( writeqp .and. dft_code == PARSEC ) &
          call zprint_qp(nspin,gvec%syms%ntrans,66,77,gvec,kpt)
  else
     if ( writeqp .and. dft_code == PARSEC ) &
          call dprint_qp(nspin,gvec%syms%ntrans,66,77,gvec,kpt)
  endif

  if (n_it > 0 .and. max_sig > max_conv .and. peinf%master) then
     write(6,'(/,a,/)') repeat('*',65)
     write(6,*) 'WARNING!!!! self-consistent solution did ', &
          'not converge after ', it_scf - 1, ' iterations!'
     write(6,'(a,f10.5)') ' Maximum potential (eV) = ', max_sig
     write(6,'(/,a,/)') repeat('*',65)
  endif

99 continue

  !-------------------------------------------------------------------
  ! Time accounting.
  !
  ii = 9
  ik = 7
  allocate(routnam(ii+ik))
  allocate(timerlist(ii+ik))
  routnam(1)='SETUP_S:'      ; timerlist(1)=2
  routnam(2)='KERNEL:'       ; timerlist(2)=3
  routnam(3)='DIAG_POL:'     ; timerlist(3)=4
  routnam(4)='WPOL_0:'       ; timerlist(4)=5
  routnam(5)='EXCHANGE:'     ; timerlist(5)=6
  routnam(6)='POTENTIAL:'    ; timerlist(6)=7
  routnam(7)='CORRELATION:'  ; timerlist(7)=8
  routnam(8)='VERTEX:'       ; timerlist(8)=9
  routnam(9)='ROTATE:'       ; timerlist(9)=10

  routnam(10)='POISSON_FFT:'       ; timerlist(10)=11
  routnam(11)='EIGENSOLVER:'       ; timerlist(11)=12
  routnam(12)='INTEGRATION:'       ; timerlist(12)=13
  routnam(13)='Find intp pts:'     ; timerlist(13)=51
  routnam(14)='ISDF_PARALLEL:'     ; timerlist(14)=52
  routnam(15)='Calc intp vectors:' ; timerlist(15)=53
  routnam(16)='Calc <zeta|K|zeta>:'; timerlist(16)=54

  call timacc(1,2,tsec)
  call finalize(peinf%master,peinf%comm,ii,ik,routnam, timerlist)

end program sigma
!===================================================================
