# NanoGW_optimized
This code is built from the original NanoGW code written by M. Tiago. New experimental features (such as interpolative separable density fitting method) are added to accelerate calculations. 


### Feb 13. 2018
Implemented the isdf method. 
New files include: 
* isdf.f90
* cvt.f90
* k_integrate_isdf.f90.

Modified files are:
* tdlda.F90
* kernel.F90z
* calculate_tdlda.F90z
* esdf_key_module.f90
* fft_module.F90z
* input_t.f90

### May 23. 2018
The ISDF part is parallelized. It requires less memory, but still not very efficient.
Tests show that the CVT part (finding interpolation points) is very slow. 
I also added an option in utils/chkpt_bin_asc.f90 to choose static_type.

Things to be changed/optimized:
1. cvt.f90
2. input_s.f90, some keywords are not very clear
 Refer to RGWBS-new version modified by Linda
3. the printed out Sigma-Vxc (in "sigma_***__****") is calculated at E_dft, this is confusing.
