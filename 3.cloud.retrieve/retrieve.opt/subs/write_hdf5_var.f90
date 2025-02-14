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

SUBROUTINE write_hdf5_var3(filename,nscan,npix,nch,obs,emissivity,telsem2, &
							cfr,lsmask,stime,lon,lat,retcode)

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
      INTEGER :: rank ,retcode                ! Dataset rank
      INTEGER(HSIZE_T), DIMENSION(3)  :: data_dims3
      INTEGER(HSIZE_T), DIMENSION(2)  :: data_dims2
      INTEGER(HSIZE_T), DIMENSION(1)  :: data_dims1
	
	  !! attribute
	  INTEGER(HID_T) :: attr_id, attr_id0       ! Attribute identifier
	  INTEGER(HID_T) :: aspace_id,  aspace_id0     ! Attribute Dataspace identifier
	  INTEGER(HID_T) :: atype_id,   atype_id0      ! Attribute Dataspace identifier
	  INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
	  INTEGER     ::   arank = 1                      ! Attribure rank
	  INTEGER(SIZE_T) :: attrlen=80    ! Length of the attribute string
	
      REAL, DIMENSION(nch,npix,nscan) :: obs,emissivity,telsem2,cfr 
      INTEGER, DIMENSION(npix,nscan)  :: lsmask
      REAL, DIMENSION(nscan) 		  :: stime
      REAL, DIMENSION(npix,nscan) 	  :: lon, lat  
      INTEGER :: error 

  ! =====================================================================
 
      ! Initialize the dset_data array 
      data_dims3 = (/nch,npix,nscan/)
      data_dims2 = (/npix,nscan/)
      data_dims1 = (/nscan/)
	  
      ! Initialize Fortran interface
      CALL h5open_f(error) 
      CALL h5fcreate_f(trim(adjustl(filename)), H5F_ACC_TRUNC_F, file_id, error)
	  
	  rank = 3	  
      CALL h5screate_simple_f(rank, data_dims3, dspace_id, error)  	  
      dsetname="MERSI_Cloud_Fraction"		  
	  CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, cfr, data_dims3, error)
	  	  CALL h5screate_simple_f(arank, adims, aspace_id, error)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)

		  attrlen=1
		  CALL h5tset_size_f(atype_id, attrlen, error)
		  CALL h5acreate_f(dset_id, "units", atype_id, aspace_id, attr_id, error)
		  CALL h5awrite_f(attr_id, atype_id, "%", adims, error)
		  CALL h5aclose_f(attr_id, error)
		  attrlen=7
		  CALL h5tset_size_f(atype_id, attrlen, error)		  
		  CALL h5acreate_f(dset_id, "valid_range", atype_id, aspace_id, attr_id, error)
		  CALL h5awrite_f(attr_id, atype_id, "[0:100]", adims, error)
		  CALL h5aclose_f(attr_id, error)
		  attrlen=14
		  CALL h5tset_size_f(atype_id, attrlen, error)		  
		  CALL h5acreate_f(dset_id, "long_name", atype_id, aspace_id, attr_id, error)
		  CALL h5awrite_f(attr_id, atype_id, "Cloud Fraction", adims, error)
		  CALL h5aclose_f(attr_id, error)
		  
		  attrlen=95
		  CALL h5tset_size_f(atype_id, attrlen, error)
		  CALL h5acreate_f(dset_id, "description", atype_id, aspace_id, attr_id, error)
		  CALL h5awrite_f(attr_id, atype_id, &
 "Footprint-averged Cloud Fraction for each MWRI channel derived from MERSI L1 Cloud Mask product"&
 , adims, error)
		  CALL h5aclose_f(attr_id, error)		  
		  
		  CALL h5tclose_f(atype_id, error)
		  CALL h5sclose_f(aspace_id, error)	
      CALL h5dclose_f(dset_id, error)  
	  

	  
	  dsetname="TBs"		  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, obs, data_dims3, error)
	  	  CALL h5screate_simple_f(arank, adims, aspace_id, error)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)

		  attrlen=1
		  CALL h5tset_size_f(atype_id, attrlen, error)
		  CALL h5acreate_f(dset_id, "units", atype_id, aspace_id, attr_id, error)
		  CALL h5awrite_f(attr_id, atype_id, "K", adims, error)
		  CALL h5aclose_f(attr_id, error)

		  attrlen=35
		  CALL h5tset_size_f(atype_id, attrlen, error)		  
		  CALL h5acreate_f(dset_id, "long_name", atype_id, aspace_id, attr_id, error)
		  CALL h5awrite_f(attr_id, atype_id, "MWRI L1 TOA Brightness Temperatures", adims, error)
		  CALL h5aclose_f(attr_id, error)
		  attrlen=48
		  CALL h5tset_size_f(atype_id, attrlen, error)
		  CALL h5acreate_f(dset_id, "description", atype_id, aspace_id, attr_id, error)
		  CALL h5awrite_f(attr_id, atype_id, "channel: 10V,10H,18V,18H,23V,23H,36V,36H,89V,89H", adims, error)
		  CALL h5aclose_f(attr_id, error)		  
		  
		  CALL h5tclose_f(atype_id, error)
		  CALL h5sclose_f(aspace_id, error)	 
	  
      CALL h5dclose_f(dset_id, error)

      dsetname="Emissivity"		  
	  CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, emissivity, data_dims3, error)
	  	  CALL h5screate_simple_f(arank, adims, aspace_id, error)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
		  attrlen=5
		  CALL h5tset_size_f(atype_id, attrlen, error)
		  CALL h5acreate_f(dset_id, "valid_range", atype_id, aspace_id, attr_id, error)
		  CALL h5awrite_f(attr_id, atype_id, "[0:1]", adims, error)
		  CALL h5aclose_f(attr_id, error)
		  attrlen=28
		  CALL h5tset_size_f(atype_id, attrlen, error)		  
		  CALL h5acreate_f(dset_id, "long_name", atype_id, aspace_id, attr_id, error)
		  CALL h5awrite_f(attr_id, atype_id, "Microwave Surface Emissivity", adims, error)
		  CALL h5aclose_f(attr_id, error)
		  attrlen=39
		  CALL h5tset_size_f(atype_id, attrlen, error)
		  CALL h5acreate_f(dset_id, "description", atype_id, aspace_id, attr_id, error)
		  CALL h5awrite_f(attr_id, atype_id, "Retrieved using RTTOV13,ERA5,MERSI,MWRI", adims, error)
		  CALL h5aclose_f(attr_id, error)		  
		  
		  CALL h5tclose_f(atype_id, error)
		  CALL h5sclose_f(aspace_id, error)	  	  
      CALL h5dclose_f(dset_id, error)  
	  
      dsetname="TELSEM2_Emissivity"		  
	  CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, telsem2, data_dims3, error)
	  	  CALL h5screate_simple_f(arank, adims, aspace_id, error)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
		  attrlen=5
		  CALL h5tset_size_f(atype_id, attrlen, error)
		  CALL h5acreate_f(dset_id, "valid_range", atype_id, aspace_id, attr_id, error)
		  CALL h5awrite_f(attr_id, atype_id, "[0:1]", adims, error)
		  CALL h5aclose_f(attr_id, error)
		  attrlen=39
		  CALL h5tset_size_f(atype_id, attrlen, error)		  
		  CALL h5acreate_f(dset_id, "long_name", atype_id, aspace_id, attr_id, error)
		  CALL h5awrite_f(attr_id, atype_id, "TELSEM2 Monthly Land Surface Emissivity", adims, error)
		  CALL h5aclose_f(attr_id, error)
		  attrlen=47
		  CALL h5tset_size_f(atype_id, attrlen, error)
		  CALL h5acreate_f(dset_id, "description", atype_id, aspace_id, attr_id, error)
		  CALL h5awrite_f(attr_id, atype_id, "TELSEM2 Modeled Monthly Emissivity in RTTOV13.2", adims, error)
		  CALL h5aclose_f(attr_id, error)		  
		  
		  CALL h5tclose_f(atype_id, error)
		  CALL h5sclose_f(aspace_id, error)	  	  
      CALL h5dclose_f(dset_id, error)  
  	  
      CALL h5sclose_f(dspace_id, error)
 
 
	  rank = 2	  	  
      CALL h5screate_simple_f(rank, data_dims2, dspace_id, error)
      dsetname="Latitude"	  	  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, lat, data_dims2, error)
	  	  CALL h5screate_simple_f(arank, adims, aspace_id, error)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
		  CALL h5tset_size_f(atype_id, attrlen, error)
		  attrlen=12
		  CALL h5tset_size_f(atype_id, attrlen, error)
		  CALL h5acreate_f(dset_id, "units", atype_id, aspace_id, attr_id, error)
		  CALL h5awrite_f(attr_id, atype_id, "degree_north", adims, error)
		  CALL h5aclose_f(attr_id, error)
		  
		  attrlen=31
		  CALL h5tset_size_f(atype_id, attrlen, error)		  
		  CALL h5acreate_f(dset_id, "description", atype_id, aspace_id, attr_id, error)
		  CALL h5awrite_f(attr_id, atype_id, "Latitude of each pixel in WGS84", adims, error)
		  CALL h5aclose_f(attr_id, error)		  
		  
		  CALL h5tclose_f(atype_id, error)
		  CALL h5sclose_f(aspace_id, error)
      CALL h5dclose_f(dset_id, error)

      dsetname="Longitude"	  	  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, lon, data_dims2, error)
	  	  CALL h5screate_simple_f(arank, adims, aspace_id, error)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)

		  attrlen=11
		  CALL h5tset_size_f(atype_id, attrlen, error)
		  CALL h5acreate_f(dset_id, "units", atype_id, aspace_id, attr_id, error)
		  CALL h5awrite_f(attr_id, atype_id, "degree_east", adims, error)
		  CALL h5aclose_f(attr_id, error)
		  
		  attrlen=32
		  CALL h5tset_size_f(atype_id, attrlen, error)		  
		  CALL h5acreate_f(dset_id, "description", atype_id, aspace_id, attr_id, error)
		  CALL h5awrite_f(attr_id, atype_id, "Longitude of each pixel in WGS84", adims, error)
		  CALL h5aclose_f(attr_id, error)		  
		  
		  CALL h5tclose_f(atype_id, error)
		  CALL h5sclose_f(aspace_id, error)
      CALL h5dclose_f(dset_id, error)	  
	   	  
      dsetname="Land_Sea_Mask"	  	  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_INTEGER, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, lsmask, data_dims2, error)
	  	  CALL h5screate_simple_f(arank, adims, aspace_id, error)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
		  

		  attrlen=13
		  CALL h5tset_size_f(atype_id, attrlen, error)
		  CALL h5acreate_f(dset_id, "long_name", atype_id, aspace_id, attr_id, error)
		  CALL h5awrite_f(attr_id, atype_id, "Land Sea Mask", adims, error)
		  CALL h5aclose_f(attr_id, error)
		  
		  attrlen=74
		  CALL h5tset_size_f(atype_id, attrlen, error)
		  CALL h5acreate_f(dset_id, "description", atype_id, aspace_id, attr_id, error)
		  CALL h5awrite_f(attr_id, atype_id, &
		  "The type of earth surface, 1 land, 2 continental water, 3 sea, 5  boundary", adims, error)
		  CALL h5aclose_f(attr_id, error)		  
		  
		  CALL h5tclose_f(atype_id, error)
		  CALL h5sclose_f(aspace_id, error)
      CALL h5dclose_f(dset_id, error)	  
	   
      CALL h5sclose_f(dspace_id, error)
	  	  	 

	  rank = 1	  	  
      CALL h5screate_simple_f(rank, data_dims1, dspace_id, error)
      dsetname="Scan_Time_UTC"	  	  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, stime, data_dims1, error)
		
		  CALL h5screate_simple_f(arank, adims, aspace_id, error)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
		  
		  attrlen=18
		  CALL h5tset_size_f(atype_id, attrlen, error)
		  CALL h5acreate_f(dset_id, "long_name", atype_id, aspace_id, attr_id, error)
		  CALL h5awrite_f(attr_id, atype_id, "FY3G Scan Time UTC", adims, error)
		  CALL h5aclose_f(attr_id, error)
		  
		  attrlen=28
		  CALL h5tset_size_f(atype_id, attrlen, error)	  
		  CALL h5acreate_f(dset_id, "description", atype_id, aspace_id, attr_id, error)
		  CALL h5awrite_f(attr_id, atype_id, "Format: HH+MM/60+SS/3600+...", adims, error)
		  CALL h5aclose_f(attr_id, error)
		  
		  CALL h5tclose_f(atype_id, error)
		  CALL h5sclose_f(aspace_id, error)
		  
      CALL h5dclose_f(dset_id, error)  
      CALL h5sclose_f(dspace_id, error)

     ! Close the file
     CALL h5fclose_f(file_id, error)
     ! Close FORTRAN interface
     CALL h5close_f(error)
	 
END SUBROUTINE write_hdf5_var3