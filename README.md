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
