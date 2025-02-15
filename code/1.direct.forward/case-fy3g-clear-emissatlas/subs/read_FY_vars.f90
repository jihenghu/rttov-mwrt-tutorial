! module fysubs   
    ! contains

subroutine read_FY_Vars(filename,nscan,npix,lon,lat,tbs,minsec,igbp,lsmask,dem,retcode)

	USE HDF5   
	CHARACTER*(*) :: filename
	INTEGER, PARAMETER :: nchn = 10
	INTEGER :: nscan,npix

	INTEGER(HID_T) :: file_id       ! File identifier  文件句柄
	INTEGER(HID_T) :: dset_id       ! Dataset identifier 变量句柄
	INTEGER(HID_T) :: grp_id        ! Dataset identifier  group句柄
	INTEGER(HID_T) :: dspace_id     ! Dataset identifier  只可意会，和维数和大小相关
	INTEGER     	 :: error ,retcode          ! Error flag - success：0

	INTEGER:: tb(nchn, npix, nscan) ,dem(npix, nscan)
	Real*4 :: tbs(nchn, npix, nscan) 
	Real*4 :: lon(npix, nscan),lat(npix, nscan)

	Real*8 :: minsec(nscan)

	INTEGER  :: igbp(npix, nscan)   
	INTEGER  :: lsmask(npix, nscan) 

	INTEGER(HSIZE_T), DIMENSION(3) :: dims_3D,maxdims_3D 
	INTEGER(HSIZE_T), DIMENSION(2) :: dims_2D,maxdims_2D
	INTEGER(HSIZE_T), DIMENSION(1) :: dims_1D,maxdims_1D
  
	! print*,trim(adjustl(filename))
	CALL h5open_f(error) 
	CALL h5fopen_f (trim(adjustl(filename)), H5F_ACC_RDONLY_F, file_id, error) 
			if (error.ne.0) then
			  retcode = -1
			  return
			end if
	!! group S1 
	CALL h5gopen_f (file_id, "S1/Data", grp_id, error)  

		!!!!! EARTH_OBSERVE_BT_10_to_89GHz
		CALL h5dopen_f(grp_id, "EARTH_OBSERVE_BT_10_to_89GHz", dset_id, error)  
		CALL h5dget_space_f(dset_id,dspace_id,error)  
		CALL h5sget_simple_extent_dims_f(dspace_id, dims_3D, maxdims_3D, error)   
		CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, tb, dims_3D, error)
		CALL h5dclose_f(dset_id, error) 
		tbs=tb*0.01+327.68

		!!!!! DEM
		CALL h5dopen_f(grp_id, "DEM", dset_id, error)  
		CALL h5dget_space_f(dset_id,dspace_id,error)  
		CALL h5sget_simple_extent_dims_f(dspace_id, dims_2D, maxdims_2D, error)   
		CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, dem, dims_2D, error)
		CALL h5dclose_f(dset_id, error) 

		!!!!!LandSeaMask
		CALL h5dopen_f(grp_id, "LandSeaMask", dset_id, error)  
		CALL h5dget_space_f(dset_id,dspace_id,error)  
		CALL h5sget_simple_extent_dims_f(dspace_id, dims_2D, maxdims_2D, error)   
		CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, lsmask, dims_2D, error)
		CALL h5dclose_f(dset_id, error) 
		
		!!!!!LandCover
		CALL h5dopen_f(grp_id, "LandCover", dset_id, error)  
		CALL h5dget_space_f(dset_id,dspace_id,error)  
		CALL h5sget_simple_extent_dims_f(dspace_id, dims_2D, maxdims_2D, error)   
		CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, igbp, dims_2D, error)
		CALL h5dclose_f(dset_id, error) 
		
		!!!!!Scan_Mscnt
		CALL h5dopen_f(grp_id, "Scan_Mscnt", dset_id, error)  
		CALL h5dget_space_f(dset_id,dspace_id,error)  
		CALL h5sget_simple_extent_dims_f(dspace_id, dims_1D, maxdims_1D, error)   
		CALL h5dread_f(dset_id, H5T_IEEE_F64LE, minsec, dims_1D, error)
		CALL h5dclose_f(dset_id, error) 
		minsec=minsec*0.1  !! SLOPE: 0.1
		
	CALL h5gclose_f(grp_id, error) 

	!! group Geolocation 
	CALL h5gopen_f (file_id, "S1/Geolocation", grp_id, error)  
		CALL h5dopen_f(grp_id, "Latitude", dset_id, error)  
		CALL h5dget_space_f(dset_id,dspace_id,error)  
		CALL h5sget_simple_extent_dims_f(dspace_id, dims_2D, maxdims_2D, error)   
		CALL h5dread_f(dset_id, H5T_IEEE_F32LE, lat, dims_2D, error)
		CALL h5dclose_f(dset_id, error) 

		CALL h5dopen_f(grp_id, "Longitude", dset_id, error)  
		CALL h5dread_f(dset_id, H5T_IEEE_F32LE, lon, dims_2D, error)
		CALL h5dclose_f(dset_id, error) 		
	CALL h5gclose_f(grp_id, error) 		
	CALL h5fclose_f(file_id, error)  	
    CALL h5close_f(error)				
end subroutine
! end module fysubs


subroutine read_FY_Geometry(filename,nscan,npix,sea,sez,soa,soz,retcode)

	USE HDF5   
	CHARACTER*(*) :: filename
	INTEGER :: nscan,npix

	INTEGER(HID_T) :: file_id     
	INTEGER(HID_T) :: dset_id       
	INTEGER(HID_T) :: grp_id       
	INTEGER(HID_T) :: dspace_id     
	INTEGER        :: error ,retcode          ! Error flag - success：0

	Real*4, DIMENSION(npix, nscan) :: sea,sez,soa,soz
	INTEGER, DIMENSION(npix, nscan) :: temp


	INTEGER(HSIZE_T), DIMENSION(2) :: dims_2D,maxdims_2D

  
	! print*,trim(adjustl(filename))
	CALL h5open_f(error) 
	CALL h5fopen_f (trim(adjustl(filename)), H5F_ACC_RDONLY_F, file_id, error) 
			if (error.ne.0) then
			  retcode = -1
			  return
			end if

	!! group Geolocation 
	CALL h5gopen_f (file_id, "S1/Geolocation", grp_id, error) 
	
		CALL h5dopen_f(grp_id, "Sensor_Azimuth", dset_id, error)  
		CALL h5dget_space_f(dset_id,dspace_id,error)  
		CALL h5sget_simple_extent_dims_f(dspace_id, dims_2D, maxdims_2D, error)   
		CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, temp, dims_2D, error)
		CALL h5dclose_f(dset_id, error) 
		sea=temp*0.01
	
		CALL h5dopen_f(grp_id, "Sensor_Zenith", dset_id, error)  
		CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, temp, dims_2D, error)
		CALL h5dclose_f(dset_id, error) 
		sez=temp*0.01
		
		CALL h5dopen_f(grp_id, "Solar_Azimuth", dset_id, error)  
		CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, temp, dims_2D, error)
		CALL h5dclose_f(dset_id, error) 
		soa=temp*0.01
		
		CALL h5dopen_f(grp_id, "Solar_Zenith", dset_id, error)  
		CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, temp, dims_2D, error)
		CALL h5dclose_f(dset_id, error) 
		soz=temp*0.01
	
	CALL h5gclose_f(grp_id, error) 		
	CALL h5fclose_f(file_id, error)  	
    CALL h5close_f(error)				
end subroutine


subroutine get_FY_swath_nscan(filein,nscan,npix,retcode)
    USE HDF5 
    implicit none 
    character*(*) :: filein
    character(len=200)::varname
    integer:: error
    integer(HID_T) :: file_id,dset_id,dspace_id
    integer(HSIZE_T) :: data_dims3(3),maxdims3(3)
    
    integer:: retcode 
    integer :: nscan ,npix
    
    retcode = 0
    call h5open_f(error)
    call h5fopen_f(filein, H5F_ACC_RDONLY_F, file_id, error)
    if (error.ne.0) then
      retcode = -1
      return
    end if
    
    varname = "S1/Data/EARTH_OBSERVE_BT_10_to_89GHz"
    CALL h5dopen_f(file_id, trim(varname), dset_id, error)
    if (error.ne.0) then
      retcode = -1
      return
    end if
    CALL h5dget_space_f(dset_id,dspace_id,error)
    CALL h5sget_simple_extent_dims_f(dspace_id, data_dims3, maxdims3, error)
    nscan = data_dims3(3) 
    npix = data_dims3(2) 
    CALL h5sclose_f(dspace_id,error)
    CALL h5dclose_f(dset_id,error)
    
    
    CALL h5fclose_f(file_id,error)
    CALL h5close_f(error)
end subroutine 