gfortran -mcmodel=large -ffpe-summary=none -I/home/hjh/netcdf/include -I/home/hjh/hdf5/include -L/home/hjh/netcdf/lib -lnetcdff -lnetcdf  -L/home/hjh/hdf5/lib -lhdf5_fortran -c fy3g_mwri_fwd_cloudy.f90 -o ./fy3g_mwri_fwd_cloudy.o
gfortran -I/home/hjh/rttov13/mod -I/home/hjh/rttov13/include -fPIC -O3 -fopenmp -ffree-line-length-none  -c rttov_cloudy_fwd.f90 -o ./rttov_cloudy_fwd.o
gfortran -o ./fy3g_mwri_fwd_cloudy_swath.exe \
./fy3g_mwri_fwd_cloudy.o ./rttov_cloudy_fwd.o \
-L/home/hjh/rttov13/lib -lrttov13_brdf_atlas -lrttov13_emis_atlas -lrttov13_mw_scatt -lrttov13_other -lrttov13_coef_io -lrttov13_hdf -lrttov13_parallel -lrttov13_main  \
-L/home/hjh/netcdf/lib -lnetcdff -L/home/hjh/hdf5/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz -fopenmp 