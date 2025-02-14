SUBROUTINE write_hdf5_var(filename,nscan,npix,nch,tbs,emiss,lon,lat,retcode)

      USE HDF5   
      IMPLICIT NONE
	  Integer :: nscan,npix,nch
      CHARACTER(*):: filename
      CHARACTER*20:: dsetname != "TB_89VH_GRID25" ! Dataset name   
      ! Identifiers
      INTEGER(HID_T) :: file_id       ! File identifier
      ! INTEGER(HID_T) :: group_id      ! Group identifier
      INTEGER(HID_T) :: dset_id      ! Dataset 1 identifier
      INTEGER(HID_T) :: dspace_id    ! Dataspace 1 identifier
	  
      ! FP array
      INTEGER :: rank ,retcode                ! Dataset rank
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims3
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims2
      REAL, DIMENSION(nch,npix,nscan) :: tbs,emiss  
      REAL, DIMENSION(npix,nscan) :: lon, lat  
      INTEGER :: error 

  ! =====================================================================
 
      ! Initialize the dset_data array 
      data_dims3 = (/nch,npix,nscan/)
      data_dims2 = (/npix,nscan/)
	  
      ! Initialize Fortran interface
      CALL h5open_f(error) 
      CALL h5fcreate_f(trim(adjustl(filename)), H5F_ACC_TRUNC_F, file_id, error)
	  

	  rank = 3	  
      CALL h5screate_simple_f(rank, data_dims3, dspace_id, error)
      dsetname="TB_FY3G_simu"		  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, tbs, data_dims3, error)
      CALL h5dclose_f(dset_id, error)
	  
      dsetname="Emissivity_simu"		  
	  CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, emiss, data_dims3, error)
      CALL h5dclose_f(dset_id, error)  
      CALL h5sclose_f(dspace_id, error)
 
	  rank = 2	  	  
      CALL h5screate_simple_f(rank, data_dims2, dspace_id, error)
      dsetname="Latitude"	  	  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, lat, data_dims2, error)
      CALL h5dclose_f(dset_id, error)
	  
      dsetname="Longitude"	  	  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, lon, data_dims2, error)
      CALL h5dclose_f(dset_id, error)	  
	   
      CALL h5sclose_f(dspace_id, error)

     ! Close the file
     CALL h5fclose_f(file_id, error)
     ! Close FORTRAN interface
     CALL h5close_f(error)
	 
END SUBROUTINE write_hdf5_var

SUBROUTINE write_hdf5_var2(filename,nscan,npix,nch,tbs,emiss,obs,emissivity,lon,lat,retcode)

      USE HDF5   
      IMPLICIT NONE
	  Integer :: nscan,npix,nch
      CHARACTER(*):: filename
      CHARACTER*20:: dsetname != "TB_89VH_GRID25" ! Dataset name   
      ! Identifiers
      INTEGER(HID_T) :: file_id       ! File identifier
      ! INTEGER(HID_T) :: group_id      ! Group identifier
      INTEGER(HID_T) :: dset_id      ! Dataset 1 identifier
      INTEGER(HID_T) :: dspace_id    ! Dataspace 1 identifier
	  
      ! FP array
      INTEGER :: rank ,retcode                ! Dataset rank
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims3
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims2
      REAL, DIMENSION(nch,npix,nscan) :: tbs,emiss,obs,emissivity  
      REAL, DIMENSION(npix,nscan) :: lon, lat  
      INTEGER :: error 

  ! =====================================================================
 
      ! Initialize the dset_data array 
      data_dims3 = (/nch,npix,nscan/)
      data_dims2 = (/npix,nscan/)
	  
      ! Initialize Fortran interface
      CALL h5open_f(error) 
      CALL h5fcreate_f(trim(adjustl(filename)), H5F_ACC_TRUNC_F, file_id, error)
	  

	  rank = 3	  
      CALL h5screate_simple_f(rank, data_dims3, dspace_id, error)
      dsetname="TB_FY3G_simu"		  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, tbs, data_dims3, error)
      CALL h5dclose_f(dset_id, error)
	  
      dsetname="Emissivity_simu"		  
	  CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, emiss, data_dims3, error)
      CALL h5dclose_f(dset_id, error)  
	  
	  dsetname="TB_FY3G_obs"		  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, obs, data_dims3, error)
      CALL h5dclose_f(dset_id, error)
	  
      dsetname="Emissivity_retrievals"		  
	  CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, emissivity, data_dims3, error)
      CALL h5dclose_f(dset_id, error)  
	  
	  
      CALL h5sclose_f(dspace_id, error)
 
	  rank = 2	  	  
      CALL h5screate_simple_f(rank, data_dims2, dspace_id, error)
      dsetname="Latitude"	  	  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, lat, data_dims2, error)
      CALL h5dclose_f(dset_id, error)
	  
      dsetname="Longitude"	  	  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, lon, data_dims2, error)
      CALL h5dclose_f(dset_id, error)	  
	   
      CALL h5sclose_f(dspace_id, error)

     ! Close the file
     CALL h5fclose_f(file_id, error)
     ! Close FORTRAN interface
     CALL h5close_f(error)
	 
END SUBROUTINE write_hdf5_var2