!!! RTTOV forward MW RT modeling under clear-sky
!!! &copy; Jiheng Hu 2023,August
!!! University of Science and Technology of China


!!! follow insturction of : https://cds.climate.copernicus.eu/api-how-to#install-the-cds-api-key
!!! to facilitate  ERA5 reanalysis download

!!  Ensure internet accessing, by 
!!  $ wlt > wlt.log  

	include './subs/read_FY_vars.f90' 
	include './subs/read_ERA5_vars.f90' 
	include './subs/write_hdf5_var.f90' 
	! include './rttov_clearsky_fwd.f90' 

program fy3g_mwri_fwd_clear

!! Header Start
	CHARACTER*200 :: FY_FILENAME,QuerySet
	CHARACTER*200 :: FY_Folder,FY_Name,pwd,filelist,fileout
	CHARACTER*8   :: yyyymmdd

	!!!  FY vars
	INTEGER :: retcode
	INTEGER :: nscan, npix, iscan, ipix
	INTEGER, Parameter :: nchn=10
    REAL*4 , allocatable, DIMENSION(:,:)  :: FYlon,FYlat
	REAL*4 , allocatable, DIMENSION(:,:,:):: FYtbs
	INTEGER, allocatable, DIMENSION(:,:) :: igbp     
	INTEGER, allocatable, DIMENSION(:,:) :: lsmask     
	INTEGER, allocatable, DIMENSION(:,:):: dem     
	REAL*8 , allocatable, DIMENSION(:) :: minsec
	REAL*4 , allocatable, DIMENSION(:) :: fystime
	REAL*4 , allocatable, DIMENSION(:,:) :: sea,sez,soa,soz !! sensor azimuth, sensor zenith, solar azimu, solar zenith
	integer  :: ls,dem0	
	real*4  sea0,sez0,soa0,soz0,flon,flat
	
	!!! ERA5 Vars
	INTEGER :: UTC_hr,fystart
	REAL*4 ::fystart_hr
	CHARACTER :: UTC_str*2,req_era5*200,req_land*200,req_sigle*200
	LOGICAL :: alive1, alive2, alive3
	!! ERA5
	INTEGER,Parameter :: nlon=1440, nlat=521,nlevel=29, ntime=1
	INTEGER           :: ilon,ilat
	real*4,Parameter  :: grid=0.25
	REAL*4, DIMENSION(nlevel) 				  :: plevel !!hpa 	
	REAL*4 									  :: lon(nlon), lat(nlat) 
	REAL*4, DIMENSION(nlon,nlat,nlevel,ntime) :: qw,rh,ta,zg !!kg/kg, %, K, m2/s2 

	!!ERA5 land
	REAL*4, DIMENSION(nlon,nlat,ntime) 		  :: lst,psrf_land,t2m_land !! K, hPa, K
	! REAL*4, DIMENSION(nlon,nlat,ntime)      :: stl1,swvl1,snowc,sde !! K, m3/m3, %, m
	
	!! ERA5 single level
	REAL*4, DIMENSION(nlon,nlat,ntime) 		  :: u10,v10 !!m/s
	REAL*4, DIMENSION(nlon,nlat,ntime) 		  :: sst,psrf_gbl,t2m_gbl !! K, hPa, K
	
	real*4, dimension(nlevel) :: qw0,ta0
	REAL*4   :: skt,psrf,t2m,u,v
	
	!!! RTTOV OUT
	REAL*4,allocatable,DIMENSION(:,:,:) ::  TBsimu
	!!=================================================================================================

	plevel=(/50.0,70.0,100.0,125.0,150.0,175.0,200.0,225.0,250.0,300.0,350.0,400.0,450.0,500.0	&
        ,550.0,600.0,650.0,700.0,750.0,775.0,800.0,825.0,850.0,875.0,900.0,925.0,950.0,975.0,1000.0/)

!! Directly specify a FY3G fullname to process,
!! e.g., /data04/0/MWRI/FY3G_IOT/ASCEND/20230626/FY3G_MWRI-_ORBA_L1_20230626_0007_7000M_V0.HDF

	! call getarg(1,FY_FILENAME)

!! Or, specify YYYYMMDD to looping process all files in a specific day
	call getarg(1,yyyymmdd)
	
	pwd="/data04/1/hjh/rttov.mlse.algor/1.direct.forward/case-fy3g-clearsky/"
	FY_Folder="/data04/0/MWRI/FY3G_IOT/ASCEND/"//yyyymmdd//"/"
	FY_Name="FY3G_MWRI-_ORBA_L1_"//yyyymmdd//"_????_7000M_V0.HDF"
	QuerySet="cd "//trim(adjustl(FY_Folder))//"; ls "//trim(adjustl(FY_Name))
	filelist=trim(adjustl(pwd))//"/FY3G_filelist_"//yyyymmdd//".txt"

	call system(trim(adjustl(QuerySet))//"> "//trim(adjustl(filelist)))

	open(701,file=TRIM(ADJUSTL(filelist)))
702	read(701,*,end=200) FY_Name  
	FY_FILENAME=trim(adjustl(FY_Folder))//trim(adjustl(FY_Name))
	print*,FY_FILENAME

!! extract FY MWRI IOT Variables
    call get_FY_swath_nscan(FY_FILENAME,nscan,npix,retcode)
	
	allocate(FYtbs(nchn, npix, nscan)) 
	allocate(TBsimu(nchn, npix, nscan)) 
	allocate(FYlon(npix, nscan))
	allocate(FYlat(npix, nscan))
	allocate(dem(npix, nscan))
	allocate(minsec(nscan))
	allocate(igbp(npix, nscan))  
	allocate(lsmask(npix, nscan))
	allocate(FYstime(nscan))

	allocate(sea(npix, nscan))	
	allocate(sez(npix, nscan))	
	allocate(soa(npix, nscan))	
	allocate(soz(npix, nscan))	

	call read_FY_Vars(FY_FILENAME,nscan,npix,FYlon,FYlat,FYtbs,minsec,igbp,lsmask,dem,retcode)
	call read_FY_Geometry(FY_FILENAME,nscan,npix,sea,sez,soa,soz,retcode)
	! print*,MinVal(FYtbs),MaxVal(FYtbs),FYtbs(1, 1, 1:10)
	! print*,MinVal(sea),MaxVal(sea)
	! print*,MinVal(sez),MaxVal(sez)
	! print*,MinVal(soa),MaxVal(soa)
	! print*,MinVal(soz),MaxVal(soz)
	
	! stop
	
	FYstime=real(minsec)/1000./60./60.+12.  
	where(FYstime.gt.24) FYstime=FYstime-24            !! UTC	
	where(FYtbs.gt.600 .or. FYtbs.lt.30) FYtbs=-999.9  !! Kelvin
	DEALLOCATE(minsec) 	
!! extract end

!! ERA5 Variables
	UTC_hr=nint(minVal(FYstime)*0.5+maxVal(FYstime)*0.5)
	read(FY_FILENAME(69:72),'(I4)') fystart
	fystart_hr =fystart/100+mod(fystart,100)/60.
	if (abs(UTC_hr-fystart_hr).gt.3) UTC_hr=nint(fystart_hr)
	write(UTC_str,'(I0.2)') UTC_hr	
	
	!!! download Global ERA5 reanaysis:  [65, -180, -65, 180,]
	req_era5=trim(adjustl(pwd))//"ERA5/"//yyyymmdd//"/ERA5-PL-GBL-"//yyyymmdd//"-"//UTC_str//"00.nc"
	req_land=trim(adjustl(pwd))//"ERA5/"//yyyymmdd//"/ERA5-Land-GBL-"//yyyymmdd//"-"//UTC_str//"00.nc"
	req_sigle=trim(adjustl(pwd))//"ERA5/"//yyyymmdd//"/ERA5-Single-GBL-"//yyyymmdd//"-"//UTC_str//"00.nc"

	inquire(file=trim(adjustl(req_era5)), exist=alive1)
	inquire(file=trim(adjustl(req_land)), exist=alive2)
	inquire(file=trim(adjustl(req_sigle)), exist=alive3)
    if (.not.alive1) then !如果不存在 
		print*, "ERA5 file Not Found, calling Python downloading script"
		call system("python3 "//trim(adjustl(pwd))//"ERA5/era5-pl-download.py "//yyyymmdd//" "//UTC_str )	
	end if
    if (.not.alive2) then !如果不存在
		print*, "ERA5-Land file Not Found, calling Python downloading script"	
		call system("python3 "//trim(adjustl(pwd))//"ERA5/era5-land-download.py "//yyyymmdd//" "//UTC_str )	
	end if
    if (.not.alive3) then !如果不存在 
		print*, "ERA5-Single file Not Found, calling Python downloading script"
		call system("python3 "//trim(adjustl(pwd))//"ERA5/era5-single-download.py "//yyyymmdd//" "//UTC_str )		
	end if

!!! extract ERA5 Vars
	! call read_ERA5_dims(req_era5,nlon,nlat,nlevel,ntime)
	call read_ERA5_profiles(req_era5,lon,lat,qw,rh,ta,zg)
	call read_ERA5_surfs(req_land,lst,psrf_land,t2m_land)        !! K, hPa, K     !!  1 (top) -> 29 (bottom)
	call read_ERA5_single(req_sigle,sst,psrf_gbl,t2m_gbl,u10,v10)


!! Call Rttov 
	do iscan=1,nscan
	print*, iscan, ' / ', nscan
		do ipix=1,npix
		
			ilon=nint((FYlon(ipix,iscan)-minval(lon))/grid)+1
			ilat=nint((MaxVal(lat)-FYlat(ipix,iscan))/grid)+1
			if (ilon>nlon.or. ilat>nlat) then
				TBsimu(:,ipix,iscan)=-999.9
				cycle
			else
				! ls=lsmask(ipix,iscan)  !!! 1 land, 2 continental water, 3 sea, 5  boundary";
				if (lsmask(ipix,iscan).eq.1) then		
					skt = lst(ilon,ilat,1)
					ls= 0
				else if(lsmask(ipix,iscan).eq.3) then
				! cycle
					ls= 1
					skt = sst(ilon,ilat,1)
				else
					TBsimu(:,ipix,iscan)=-999.9
					cycle
				end if
				psrf = psrf_gbl(ilon,ilat,1)
				t2m  = t2m_gbl(ilon,ilat,1)
				u  	 = u10(ilon,ilat,1)
				v    = v10(ilon,ilat,1)
				qw0  = qw(ilon,ilat,:,1)
				ta0  = ta(ilon,ilat,:,1)
				
				flon = FYlon(ipix,iscan)
				flat = FYlat(ipix,iscan)				
				dem0 = dem(ipix, iscan)				
				
			    sea0 = sea(ipix, iscan)	
			    sez0 = sez(ipix, iscan)	
			    soa0 = soa(ipix, iscan)	
			    soz0 = soz(ipix, iscan)
				
				call rttov_clearsky_fwd(flon,flat,plevel,qw0,ta0,&
							skt,real(psrf),t2m,u,v,dem0,ls,sea0,sez0,soa0,soz0,&
							TBsimu(:,ipix,iscan)) 
							
				! print*,TBsimu(:,ipix,iscan)
				! stop
			end if
			
		end do
	end do
			
	fileout="RTTOV_CLEAR_SIMU_Swath_"//trim(adjustl(FY_Name))
	call write_hdf5_var(trim(adjustl(fileout)),nscan,npix,nchn,TBsimu,FYlon,FYlat,retcode)
	
	DEALLOCATE(FYtbs) 
	DEALLOCATE(FYlon) 
	DEALLOCATE(FYlat) 
	DEALLOCATE(dem) 
	DEALLOCATE(igbp) 
	DEALLOCATE(lsmask) 
	DEALLOCATE(FYstime) 
	DEALLOCATE(TBsimu)

	deallocate(sea)	
	deallocate(sez)	
	deallocate(soa)	
	deallocate(soz)
	
	goto 702
200  continue
	print*, yyyymmdd//": process done"
end program fy3g_mwri_fwd_clear
