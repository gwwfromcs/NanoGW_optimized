!===================================================================
!
! Macros for the complex pre-processing. All source files with extension
! .F90z are pre-processed twice. Once with CPLX and then again without CPLX.
! Most instances of capital "Z" are replaced by either lower case "z" or
! lower case "d" respectively with or without CPLX.
! BE CAREFUL WITH CAPITAL "Z" IN .F90z FILES !!! IT HAS SPECIAL MEANING!!!
!
! Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
! This file is part of RGWBS. It is distributed under the GPL v1.
!
!-------------------------------------------------------------------
#ifdef CPLX
! Constants
#define Zzero zzero
#define Zone zone
#define MPI_DOUBLE_SCALAR MPI_DOUBLE_COMPLEX
! Implicit functions
#define SCALAR complex(dpc)
#define MYCONJG(x) conjg(x)
! External functions and subroutines
#define Zdot_u zdot_u
#define Zdot_c zdot_c
#define Zcopy zcopy
#define Zaxpy zaxpy
#define Zgemm zgemm
#define Zvem zvem
#define Zscal zscal
#define Zgemv zgemv
! Internal functions and subroutines
#define Zgather zgather
#define Zscatter zscatter
#define Zinitialize_FFT zinitialize_FFT
#define Zfinalize_FFT zfinalize_FFT
#define Zcreate_coul_0D zcreate_coul_0D
#define Zcreate_coul_1D zcreate_coul_1D
#define Zcreate_coul_2D zcreate_coul_2D
#define Zpoisson zpoisson
#define Zfullpoisson zfullpoisson
#define Zdo_FFT zdo_FFT
#define Zmultiply_vec zmultiply_vec
#define Zeigensolver zeigensolver
#define Zdiag_pol zdiag_pol
#define Zwrite_pol zwrite_pol
#define Zread_pol zread_pol
#define Zsave_pol zsave_pol
#define Zget_pol zget_pol
#define Zcalculate_tdlda zcalculate_tdlda
#define Zkernel zkernel
#define Zk_integrate zk_integrate
#define Zk_integrate_isdf zk_integrate_isdf
#define Zmultiply_vec zmultiply_vec
#define Zpsum zpsum
#define Zget_dipole zget_dipole
#define Zsetup_g zsetup_g
#define Zsetup_s zsetup_s
#define Zsetup_b zsetup_b
#define Zdefine_pmap zdefine_pmap
#define Zprint_qp zprint_qp
#define Znonloc znonloc
#define Zget_grad_FFT zget_grad_FFT
#define Zget_grad_fd zget_grad_fd
#define Zget_lap_FFT zget_lap_FFT
#define Zget_lap_fd zget_lap_fd
#define Zmatvec3 zmatvec3
#define Zvxcmtxel zvxcmtxel
#define Zvxc zvxc
#define Zpotential zpotential
#define Zcontract_en zcontract_en
#define Zetotal zetotal
#define Zenergy_mtxel zenergy_mtxel
#define Zv_nloc zv_nloc
#define Zfock_exchange zfock_exchange
#define Zmodel_exchange zmodel_exchange
#define Zcorrelation zcorrelation
#define Zmodel_correlation zmodel_correlation
#define Zvertex zvertex
#define Zwpol0 zwpol0
#define Zwpol_v zwpol_v
#define Zstatic_corr zstatic_corr
#define Zstatic zstatic
#define Zcalculate_sigma zcalculate_sigma
#define Zgw_correlation zgw_correlation
#define Zrotate_qp zrotate_qp
#define Zpotential_mtxel zpotential_mtxel
#define Zcharac_group zcharac_group
#define Zproj_calc zproj_calc
#define Zlocalize zlocalize
#define Zblip_interp zblip_interp
#define Zget_kinetic_energy zget_kinetic_energy
#define Zcalculate_bse zcalculate_bse
#define Zread_bse zread_bse
#define Zwrite_bse zwrite_bse
#define Zdirect_b zdirect_b
#define Zdirect_s zdirect_s
#define Zdirect_mix zdirect_mix
#define Zexchange_b zexchange_b
#define Zdiag_bse zdiag_bse
#define Zdiag_bse_mix zdiag_bse_mix
#define Zgroup_reduce_bse zgroup_reduce_bse
#define Zget_phase zget_phase
#define Zk_print zk_print
! Internal arrays
#define Zbox zbox
#define Zcoul zcoul
#define Zdipole zdipole
#define Zwf zwf
#define Zm zm
#define Zm_f zm_f
#define Zv zv
#define Ztv ztv
#else
! Constants
#define Zzero zero
#define Zone one
#define MPI_DOUBLE_SCALAR MPI_DOUBLE_PRECISION
! Implicit functions
#define SCALAR real(dp)
#define MYCONJG(x) (x)
! External functions and subroutines
#define Zdot_u ddot
#define Zdot_c ddot
#define Zcopy dcopy
#define Zaxpy daxpy
#define Zgemm dgemm
#define Zvem dvem
#define Zscal dscal
#define Zgemv dgemv
! Internal functions and subroutines
#define Zgather dgather
#define Zscatter dscatter
#define Zinitialize_FFT dinitialize_FFT
#define Zfinalize_FFT dfinalize_FFT
#define Zcreate_coul_0D dcreate_coul_0D
#define Zcreate_coul_1D dcreate_coul_1D
#define Zcreate_coul_2D dcreate_coul_2D
#define Zpoisson dpoisson
#define Zfullpoisson dfullpoisson
#define Zdo_FFT ddo_FFT
#define Zmultiply_vec dmultiply_vec
#define Zeigensolver deigensolver
#define Zdiag_pol ddiag_pol
#define Zwrite_pol dwrite_pol
#define Zread_pol dread_pol
#define Zsave_pol dsave_pol
#define Zget_pol dget_pol
#define Zcalculate_tdlda dcalculate_tdlda
#define Zkernel dkernel
#define Zk_integrate dk_integrate
#define Zk_integrate_isdf dk_integrate_isdf
#define Zmultiply_vec dmultiply_vec
#define Zpsum dpsum
#define Zget_dipole dget_dipole
#define Zsetup_g dsetup_g
#define Zsetup_s dsetup_s
#define Zsetup_b dsetup_b
#define Zdefine_pmap ddefine_pmap
#define Zprint_qp dprint_qp
#define Znonloc dnonloc
#define Zget_grad_FFT dget_grad_FFT
#define Zget_grad_fd dget_grad_fd
#define Zget_lap_FFT dget_lap_FFT
#define Zget_lap_fd dget_lap_fd
#define Zmatvec3 dmatvec3
#define Zvxcmtxel dvxcmtxel
#define Zvxc dvxc
#define Zpotential dpotential
#define Zcontract_en dcontract_en
#define Zetotal detotal
#define Zenergy_mtxel denergy_mtxel
#define Zv_nloc dv_nloc
#define Zfock_exchange dfock_exchange
#define Zmodel_exchange dmodel_exchange
#define Zcorrelation dcorrelation
#define Zmodel_correlation dmodel_correlation
#define Zvertex dvertex
#define Zwpol0 dwpol0
#define Zwpol_v dwpol_v
#define Zstatic_corr dstatic_corr
#define Zstatic dstatic
#define Zcalculate_sigma dcalculate_sigma
#define Zgw_correlation dgw_correlation
#define Zrotate_qp drotate_qp
#define Zpotential_mtxel dpotential_mtxel
#define Zcharac_group dcharac_group
#define Zproj_calc dproj_calc
#define Zlocalize dlocalize
#define Zblip_interp dblip_interp
#define Zget_kinetic_energy dget_kinetic_energy
#define Zcalculate_bse dcalculate_bse
#define Zread_bse dread_bse
#define Zwrite_bse dwrite_bse
#define Zdirect_b ddirect_b
#define Zdirect_s ddirect_s
#define Zdirect_mix ddirect_mix
#define Zexchange_b dexchange_b
#define Zdiag_bse ddiag_bse
#define Zdiag_bse_mix ddiag_bse_mix
#define Zgroup_reduce_bse dgroup_reduce_bse
#define Zget_phase dget_phase
#define Zk_print dk_print
! Internal arrays
#define Zbox dbox
#define Zcoul dcoul
#define Zdipole ddipole
#define Zwf dwf
#define Zm dm
#define Zm_f dm_f
#define Zv dv
#define Ztv dtv
#endif
!===================================================================
