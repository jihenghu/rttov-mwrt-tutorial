# TELSEM Emissivity Atlas的离线使用
&copy;Jiheng Hu 2023-2030, 禁止转载。

## Intro
在Direct Forward的晴空模拟的章节，我们介绍过使用TELSEM2 Atlas进行RTTOV模拟，这里单开一节对这个数据集的离线使用进行介绍。

TELSEM 由CNRS的F.Aires和C. Priegent 在[《A Tool to Estimate Land-Surface Emissivities at Microwave frequencies (TELSEM) for use in numerical weather prediction》](#References) 中提出。TELSEM基于SSM/I地表发射率数据集([C. Prigent, 1997](#References))的19V&H, 22V, 37V&H, 85V&H等通道的Emissivity进行频率插值，得到月平均的，具有一定空间分辨率的（或固定经纬度）的全球陆面微波发射率，其适用频段为10 – 700 GHz.

本节只介绍单独使用TELSEM进行气候态地表发射率估算的方法，不会涉及如何在RTTOV的模拟中使用TEMSEM。

## 下载TELSEM2 Atlas
下载页：https://nwp-saf.eumetsat.int/site/software/rttov/download/#Emissivity_BRDF_atlas_data

下载telsem2_mw_atlas.tar.bz2解压在`/home/hjh/rttov13/emis_data/`目录下：

```bash
$ cd /home/hjh/rttov13/emis_data
$ tar -xvf telsem2_mw_atlas.tar.bz2
├── ssmi_mean_emis_climato_01_cov_interpol_M2
├── ssmi_mean_emis_climato_02_cov_interpol_M2
├── ssmi_mean_emis_climato_03_cov_interpol_M2
├── ssmi_mean_emis_climato_04_cov_interpol_M2
├── ssmi_mean_emis_climato_05_cov_interpol_M2
├── ssmi_mean_emis_climato_06_cov_interpol_M2
├── ssmi_mean_emis_climato_07_cov_interpol_M2
├── ssmi_mean_emis_climato_08_cov_interpol_M2
├── ssmi_mean_emis_climato_09_cov_interpol_M2
├── ssmi_mean_emis_climato_10_cov_interpol_M2
├── ssmi_mean_emis_climato_11_cov_interpol_M2
├── ssmi_mean_emis_climato_12_cov_interpol_M2
├── telsem2
│   ├── Makefile
│   ├── mod_mwatlas_m2.F90
│   ├── README_TELSEM2.pdf
│   ├── readme.txt
│   └── test_telsem2.F90
└── correlations
```
telsem包含了ASCII方式存储的Atlas和独立的Fortran读取代码，可以在别的地方使用。RTTOV内置了针对TELSEM2的发射率地图和相应的插值工具的支持。

## 离线使用Atlas
`mod_mwatlas_m2.F90`提供了使用telsem的接口，`test_telsem2.F90`提供了示例脚本：

较为核心的设定包括，经纬度，频率，月份和分辨率
```fortran
  resol    = 0.25
  theta    = 15.

  lat      = -30.
  lon      = 302.
  freq     = 30.                          ! For single frequency experiment
  freq2(:) = (/30., 25., 38., 60., 90. /) ! For multiple frequency experiment
  i        = 1                            ! Index of freq in freq2 array
```

运行如下：

```shell
[hjh@node05] ~/rttov13/emis_data/telsem2 $ make -f Makefile 
[hjh@node05] ~/rttov13/emis_data/telsem2 $ ./test_telsem2 
```
stdout>>
```
Reading atlas for month   9
Reading number of data in atlas...
Nb data=    233959
reading classes...
 
Inputs:
lat   =   -30.00
lon   =   302.00
theta =    15.00
freq  =    30.00
 
The first four sets of output are identical:
 
Single freq, no spatial averaging
Emis V-pol, H-pol   =   0.954829  0.954446
Stddev V-pol, H-pol =   0.020591  0.021068
Covariance V-/H-pol =   0.000372
 
Multiple freq, no spatial averaging
Emis V-pol, H-pol   =   0.954829  0.954446
Stddev V-pol, H-pol =   0.020591  0.021068
Covariance V-/H-pol =   0.000372
 
Single freq, with spatial averaging, native resol.
Emis V-pol, H-pol   =   0.954829  0.954446
Stddev V-pol, H-pol =   0.020591  0.021068
Covariance V-/H-pol =   0.000372
 
Multiple freq, with spatial averaging, native resol.
Emis V-pol, H-pol   =   0.954829  0.954446
Stddev V-pol, H-pol =   0.020591  0.021068
Covariance V-/H-pol =   0.000372
 
Now the spatial averaging is active and results are different to those above:
 
Single freq, with spatial averaging, non-native resol.
Emis V-pol, H-pol   =   0.952811  0.952665
Stddev V-pol, H-pol =   0.020790  0.023314
Covariance V-/H-pol =   0.000406
 
Multiple freq, with spatial averaging, non-native resol.
Emis V-pol, H-pol   =   0.952811  0.952665
Stddev V-pol, H-pol =   0.020790  0.023314
Covariance V-/H-pol =   0.000406
 
 
Now only return emissivities and print them for all frequencies:
 
Multiple freq, with spatial averaging, non-native resol.
Freq (GHz) :     30.000    25.000    38.000    60.000    90.000
Emis V-pol =   0.952811  0.954156  0.951043  0.953123  0.955559
Emis H-pol =   0.952665  0.954065  0.950784  0.952975  0.955488
```


## Global 0.25° MLSE
按照示例脚本，我们可以对经纬度进行循环，制作全球格点平均的地表发射率产品。

这里分享一下清扬师弟写的脚本：

```fortran
PROGRAM test_telsem
  USE mod_mwatlas_m2
  use HDF5
  IMPLICIT NONE
  INTEGER       :: ilon,ilat
  CHARACTER(2) mm

  INTEGER(jpim) :: error_status

  LOGICAL(jplm) :: verbose   ! For atlas reading subroutine
  INTEGER(jpim) :: verb      !=1 for TRUE and 0 for FALSE - for emissivity routines

  INTEGER(jpim) :: month     ! (1->12)
  CHARACTER(LEN=256) :: dir  ! directory of emis database

  REAL(jprb)    :: resol     ! horizontal resolution for the user
  REAL(jprb)    :: lat       ! (-90->90)
  REAL(jprb)    :: lon       ! (0->360)
  REAL(jprb)    :: theta     ! (0->60?)

  ! For multiple freq interpolations
  INTEGER(jpim), PARAMETER :: nchan = 5,nlon=1440,nlat=720
  REAL(jprb)    ::  ev2(nchan), eh2(nchan), std(2*nchan,2*nchan)
  REAL		    ::  emiss(nlon,nlat,nchan*2), longitude(nlon,nlat),latitude(nlon,nlat)
  REAL		    ::  MLSE(nlon,nlat,nchan*2)
  REAL(jprb)    :: freq(nchan)

  ! Structure containing atlas data
  TYPE(telsem2_atlas_data) :: atlas

	character*255 hdfname
	INTEGER(HID_T) :: file_id       ! File identifier
	INTEGER(HID_T) :: dset_id     
	INTEGER(HID_T) :: dspace_id    
	INTEGER :: rank ,status               ! Dataset rank
	INTEGER(HSIZE_T), DIMENSION(3)  :: data_dims3
	INTEGER(HSIZE_T), DIMENSION(2)  :: data_dims2

	INTEGER(HID_T) :: attr_id     ! Attribute identifier
	INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
	INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
	INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
	INTEGER     ::   arank = 1                      ! Attribure rank
	INTEGER(SIZE_T) :: attrlen=80    ! Length of the attribute string
  !--- End of header ----------------------------------


  verbose = .TRUE.  ! Verbose output for reading subroutine
  verb = 0          ! No verbose output for emissivity subroutines

  !====================================================
  ! Read the atlas
  !====================================================

	dir = '../'
	call getarg(1, mm)! = 12
	read(mm,*) month
	WRITE(0,'(a,i3)') 'Reading atlas for month ',month
	CALL rttov_readmw_atlas(TRIM(dir), month, atlas, verbose, error_status)

	IF (error_status /= 0) THEN
	WRITE(0,'(a)') 'Error reading atlas'
	STOP 1
	ENDIF

	! Multiple frequencies, spatially averaged, return emissivities only
	freq(:) = (/10.65, 18.7, 23.8, 36.64, 89.0/) ! For multiple frequency experiment
	resol    = 0.25
	theta    = 52.8

	do ilat = 1,nlat
		lat=89.875-(ilat-1)*0.25
		do ilon = 1,nlon
		  lon =0.125+(ilon-1)*0.25
		  CALL emis_interp_int_mult(lat, lon, resol, theta, freq, nchan, atlas,&
						ev2(:), eh2(:), verb = verb)	  
			emiss(ilon,ilat,1::2)=ev2(:)
			emiss(ilon,ilat,2::2)=eh2(:)
		end do
	end do
  
	MLSE(721:1440,:,:)=emiss(1:720,:,:)
	MLSE(1:720,:,:)=emiss(721:1440,:,:)
  
	do ilon=1,nlon
		longitude(ilon,:)=-179.875+(ilon-1)*0.25
	end do
	do ilat=1,nlat
		latitude(:,ilat)=89.875-(ilat-1)*0.25
	end do

	! Close the atlas
	CALL rttov_closemw_atlas(atlas)

!! output HDF
	! write(mm,'(0.2i)') month
	hdfname="GMI_TELSEM_"//mm//".HDF5"

	data_dims3 =(/nlon,nlat,10/)
	! Initialize Fortran interface
	CALL h5open_f(status) 
	CALL h5fcreate_f(trim(adjustl(hdfname)), H5F_ACC_TRUNC_F, file_id, status)	    
	
	CALL h5screate_simple_f(3, data_dims3, dspace_id, status)  	  	  
	CALL h5dcreate_f(file_id,"emissivity", H5T_NATIVE_REAL, dspace_id,dset_id, status)
	CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, MLSE, data_dims3, status)
	CALL h5dclose_f(dset_id, status)	
	CALL h5sclose_f(dspace_id, status)
	
	data_dims2=(/nlon,nlat/)
	CALL h5screate_simple_f(2, data_dims2, dspace_id, status)
	
	CALL h5dcreate_f(file_id, "latitude", H5T_NATIVE_REAL, dspace_id,dset_id, status)
	CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, latitude, data_dims2, status)
		  CALL h5screate_simple_f(arank, adims, aspace_id, status)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  attrlen=12
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "units", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "degree_north", adims, status)
		  CALL h5aclose_f(attr_id, status)
		    
		  CALL h5tclose_f(atype_id, status)
		  CALL h5sclose_f(aspace_id, status)
	CALL h5dclose_f(dset_id, status)

	CALL h5dcreate_f(file_id, "longitude", H5T_NATIVE_REAL, dspace_id,dset_id, status)
	CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, longitude, data_dims2, status)
		  CALL h5screate_simple_f(arank, adims, aspace_id, status)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  attrlen=11
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "units", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "degree_east", adims, status)
		  CALL h5aclose_f(attr_id, status)
	CALL h5dclose_f(dset_id, status)	  	
	CALL h5sclose_f(dspace_id, status)
	
	CALL h5fclose_f(file_id, status)
	CALL h5close_f(status)

END PROGRAM test_telsem
```

修改一下Makefile:
```bash
# FC=gfortran
FC=h5fc
FCFLAGS=-g -fPIC -fbounds-check -finit-integer=-1 -finit-real=nan -finit-character=127 #-std=f95
OBJS=mod_mwatlas_m2.o test_telsem2.o
MODS=mod_mwatlas_m2.mod
...

```
`sh Makefile` 即可生成`./test_telsem2`

TELSEM2 6月份地表发射率估算：

![地表发射率](./rttov132-mw-scat/mlse_in_telsem2_mlse_SSMI_06.png)

## References
1. Aires, F., Prigent, C., Bernardo, F., Jiménez, C., Saunders, R. and Brunel, P. (2011), A Tool to Estimate Land-Surface Emissivities at Microwave frequencies (TELSEM) for use in numerical weather prediction. Q.J.R. Meteorol. Soc., 137: 690-699. https://doi.org/10.1002/qj.803.
2. Prigent, C., Rossow, W. B., and Matthews, E. (1997), Microwave land surface emissivities estimated from SSM/I observations, J. Geophys. Res., 102(D18), 21867–21890. https://doi.org/10.1029/97JD01360.