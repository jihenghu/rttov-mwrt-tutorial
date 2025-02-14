
include 'read_FY_mersi.f90'

	CHARACTER*200  mersi,cfrname
	INTEGER nscan_clm, npix_clm
	INTEGER, allocatable,dimension(:,:):: clm
	REAL*4, allocatable,dimension(:,:):: cfr
	REAL*8, allocatable,dimension(:,:):: lon_clm,lat_clm
	INTEGER retcode

	mersi='/data04/0/MWRI/FY3G_IOT/PRDMERSI_L2_CLM/20230801/FY3G_MERSI_GRAN_L2_CLM_MLT_NUL_20230801_0055_0500M_V0.HDF'
	! mersi='/data04/0/MWRI/FY3G_IOT/PRDMERSI_L2_CLM/20230801/FY3G_MERSI_GRAN_L2_CLM_MLT_NUL_20230801_0140_0500M_V0.HDF'

	CALL get_mersi_swath_nscan(mersi,nscan_clm,npix_clm,retcode)	
	ALLOCATE(clm(npix_clm,nscan_clm))
	ALLOCATE(lon_clm(npix_clm/5,nscan_clm/5))
	ALLOCATE(lat_clm(npix_clm/5,nscan_clm/5))
	CALL read_mersi_clm(mersi,nscan_clm,npix_clm,clm,lon_clm,lat_clm,retcode)

	ALLOCATE(cfr(npix_clm/5,nscan_clm/5))
	cfr=-999.9
	CALL calc_cfr_5grids(clm,cfr,npix_clm/5,nscan_clm/5)
	
	WHERE(lon_clm.lt.-200) cfr=-999.9
	WHERE(lat_clm.lt.-200) cfr=-999.9
	cfrname="./cfr.20230801_0055_0500M_V0.HDF"
	CALL write_real_swath_hdf(trim(cfrname),cfr,lon_clm,lat_clm,npix_clm/5,nscan_clm/5)
	
	DEALLOCATE(clm)
	DEALLOCATE(cfr)
	DEALLOCATE(lon_clm)
	DEALLOCATE(lat_clm)
END


