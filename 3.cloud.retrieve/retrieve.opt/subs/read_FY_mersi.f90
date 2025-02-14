
subroutine read_mersi_clm(filename,nscan,npix,clm,lon,lat,retcode)

	USE HDF5   
	CHARACTER*(*) :: filename

	INTEGER(HID_T) :: file_id    
	INTEGER(HID_T) :: dset_id     
	INTEGER(HID_T) :: grp_id      
	INTEGER(HID_T) :: dspace_id  
	INTEGER(HSIZE_T), DIMENSION(2) :: dims_2D,maxdims_2D	
	INTEGER        :: error ,retcode     

	INTEGER :: nscan,npix
	INTEGER:: clm(npix,nscan)
	Real*8 :: lon(npix_clm/5,nscan_clm/5),lat(npix_clm/5,nscan_clm/5)
 
	! print*,trim(adjustl(filename))
	CALL h5open_f(error) 
	CALL h5fopen_f (trim(adjustl(filename)), H5F_ACC_RDONLY_F, file_id, error) 
			if (error.ne.0) then
			  retcode = -1
			  return
			end if

		! dims_2D=(/npix,nscan/)
		!!! Cloud_Mask_Index
		CALL h5dopen_f(file_id, "Cloud_Mask_Index", dset_id, error)  
		CALL h5dget_space_f(dset_id,dspace_id,error)  
		CALL h5sget_simple_extent_dims_f(dspace_id, dims_2D, maxdims_2D, error)   
		CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, clm, dims_2D, error)
		CALL h5dclose_f(dset_id, error) 


		!!!!!Scan_Mscnt
		! CALL h5dopen_f(file_id, "Scan_Mscnt", dset_id, error)  
		! CALL h5dget_space_f(dset_id,dspace_id,error)  
		! CALL h5sget_simple_extent_dims_f(dspace_id, dims_1D, maxdims_1D, error)   
		! CALL h5dread_f(dset_id, H5T_IEEE_F64LE, minsec, dims_1D, error)
		! CALL h5dclose_f(dset_id, error) 
		! minsec=minsec*0.1  !! SLOPE: 0.1
		

		CALL h5dopen_f(file_id, "Latitude", dset_id, error)  
		CALL h5dget_space_f(dset_id,dspace_id,error)  
		CALL h5sget_simple_extent_dims_f(dspace_id, dims_2D, maxdims_2D, error)   
		CALL h5dread_f(dset_id, H5T_IEEE_F64LE, lat, dims_2D, error)
		CALL h5dclose_f(dset_id, error) 

		CALL h5dopen_f(file_id, "Longitude", dset_id, error)  
		CALL h5dread_f(dset_id, H5T_IEEE_F64LE, lon, dims_2D, error)
		CALL h5dclose_f(dset_id, error) 		
		! print*,dims_2D
	CALL h5fclose_f(file_id, error)  	
    CALL h5close_f(error)				
end subroutine


subroutine get_mersi_swath_nscan(filein,nscan,npix,retcode)
    USE HDF5 
    implicit none 
    character*(*) :: filein
    character(len=200)::varname
    integer:: error
    integer(HID_T) :: file_id,dset_id,dspace_id
    integer(HSIZE_T) :: data_dims2(2),maxdims2(2)
    
    integer:: retcode 
    integer :: nscan ,npix
    
    retcode = 0
    call h5open_f(error)
    call h5fopen_f(filein, H5F_ACC_RDONLY_F, file_id, error)
    if (error.ne.0) then
      retcode = -1
      return
    end if
    
    varname = "Cloud_Mask_Index"
    CALL h5dopen_f(file_id, trim(varname), dset_id, error)
    if (error.ne.0) then
      retcode = -1
      return
    end if
    CALL h5dget_space_f(dset_id,dspace_id,error)
    CALL h5sget_simple_extent_dims_f(dspace_id, data_dims2, maxdims2, error)
    nscan = data_dims2(2) 
    npix = data_dims2(1) 
    CALL h5sclose_f(dspace_id,error)
    CALL h5dclose_f(dset_id,error)
    
    
    CALL h5fclose_f(file_id,error)
    CALL h5close_f(error)
end subroutine 

subroutine calc_cfr_5grids(clm,cfr,npix,nscan)

	INTEGER:: nscan,npix,iscan,ipix
	INTEGER:: clm(npix*5,nscan*5)
	Real   :: cfr(npix,nscan), ones(5,5)
	
	ones=1.0

	do iscan=1,nscan
		! print*,iscan
		do ipix= 1,npix
			cfr(ipix,iscan)=sum(ones, &
			   (clm((ipix-1)*5+1:ipix*5,(iscan-1)*5+1:iscan*5).lt.1.5))/25.0  
		!! 0:cloud,1:probably cloud,2:probably clear,3:clear,126:space,127:fillvalue	
			! cfr(ipix,iscan)=sum(ones, &
			   ! (clm((ipix-1)*5+1:ipix*5,(iscan-1)*5+1:iscan*5).gt.1.5 .and. &
			    ! clm((ipix-1)*5+1:ipix*5,(iscan-1)*5+1:iscan*5).lt.2.5))/25.0  
		!！ code table 4.217: 	0- clearOcean; 1- clearLand; 2- cloud; 3-nodata ; >4- missing   
		end do
	end do

end subroutine

subroutine write_real_swath_hdf(filename,cfr,lon,lat,npix,nscan)
	USE HDF5   
    IMPLICIT NONE
	
	INTEGER:: nscan,npix
	Real   :: cfr(npix,nscan)
	Real*8 :: lon(npix,nscan),lat(npix,nscan)
	
	CHARACTER(*):: filename
	CHARACTER*20:: dsetname 
	! Identifiers
	INTEGER(HID_T) :: file_id       ! File identifier
	! INTEGER(HID_T) :: group_id      ! Group identifier
	INTEGER(HID_T) :: dset_id      ! Dataset 1 identifier
	INTEGER(HID_T) :: dspace_id    ! Dataspace 1 identifier

	! FP array
	INTEGER :: rank ,retcode                ! Dataset rank
	INTEGER(HSIZE_T), DIMENSION(2) :: data_dims2
	INTEGER :: error 

  ! =====================================================================
 
      ! Initialize the dset_data array 
      data_dims2 = (/npix,nscan/)
	  
      ! Initialize Fortran interface
      CALL h5open_f(error) 
      CALL h5fcreate_f(trim(adjustl(filename)), H5F_ACC_TRUNC_F, file_id, error)
	  
	  rank = 2	  	  
      CALL h5screate_simple_f(rank, data_dims2, dspace_id, error)
      dsetname="Latitude"	  	  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_IEEE_F64LE, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_IEEE_F64LE, lat, data_dims2, error)
      CALL h5dclose_f(dset_id, error)
	  
      dsetname="Longitude"	  	  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_IEEE_F64LE, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_IEEE_F64LE, lon, data_dims2, error)
      CALL h5dclose_f(dset_id, error)	  
	   
	  dsetname="Cloud_Fraction_2500m"		  
	  CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, cfr, data_dims2, error)
      CALL h5dclose_f(dset_id, error)  
      CALL h5sclose_f(dspace_id, error)
	   
     ! Close the file
     CALL h5fclose_f(file_id, error)
     ! Close FORTRAN interface
     CALL h5close_f(error)
	
	
end subroutine

subroutine write_real_grid_hdf(filename,cfr,lon,lat,nlon,nlat)
	USE HDF5   
    IMPLICIT NONE
	
	INTEGER:: nlon,nlat
	Real   :: cfr(nlon,nlat)
	Real :: lon(nlon),lat(nlat)
	
	CHARACTER(*):: filename
	CHARACTER*20:: dsetname 
	! Identifiers
	INTEGER(HID_T) :: file_id       ! File identifier
	! INTEGER(HID_T) :: group_id      ! Group identifier
	INTEGER(HID_T) :: dset_id      ! Dataset 1 identifier
	INTEGER(HID_T) :: dspace_id    ! Dataspace 1 identifier

	! FP array
	INTEGER :: rank ,retcode                ! Dataset rank
	INTEGER(HSIZE_T), DIMENSION(2) :: data_dims2
	INTEGER(HSIZE_T), DIMENSION(1) :: data_dims1
	INTEGER :: error 

  ! =====================================================================
	! print*,nlon,nlat
	! print*,lon
      ! Initialize Fortran interface
      CALL h5open_f(error) 
      CALL h5fcreate_f(trim(adjustl(filename)), H5F_ACC_TRUNC_F, file_id, error)
	  
	  rank = 1		  
      data_dims1 = (/nlat/)	   	  
      CALL h5screate_simple_f(rank, data_dims1, dspace_id, error)
      dsetname="Latitude"	  	  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, lat, data_dims1, error)
      CALL h5dclose_f(dset_id, error)
	  CALL h5sclose_f(dspace_id, error)
	  
	  data_dims1 = (/nlon/)
	  CALL h5screate_simple_f(rank, data_dims1, dspace_id, error)
      dsetname="Longitude"	  	  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, lon, data_dims1, error)
      CALL h5dclose_f(dset_id, error)	
	  CALL h5sclose_f(dspace_id, error)
	  
	  rank = 2		  
	  data_dims2 = (/nlon,nlat/)
	  CALL h5screate_simple_f(rank, data_dims2, dspace_id, error)	  
	  dsetname="Cloud_Fraction_5km"		  
	  CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, cfr, data_dims2, error)
      CALL h5dclose_f(dset_id, error)  
      CALL h5sclose_f(dspace_id, error)
	   
     ! Close the file
     CALL h5fclose_f(file_id, error)
     ! Close FORTRAN interface
     CALL h5close_f(error)
	
	
end subroutine