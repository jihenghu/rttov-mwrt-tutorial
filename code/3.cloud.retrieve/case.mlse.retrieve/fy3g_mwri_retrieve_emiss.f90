!!! prototype: RTTOV forward MW RT modeling under clear-sky
!!! &copy; Jiheng Hu 2023,August
!!! University of Science and Technology of China

!!! on how to install the environment, 
!!     visit: http://home.ustc.edu.cn/~hjh18305/space/research/rttov/rttov132-column/

!!! follow insturction of : https://cds.climate.copernicus.eu/api-how-to#install-the-cds-api-key
!!! to facilitate  ERA5 reanalysis download

!!  optioanl: Ensure internet accessing, by $ wlt > wlt.log  

	include './subs/read_FY_vars.f90' 
	include './subs/read_ERA5_vars.f90' 
	include './subs/write_hdf5_var.f90' 
	! include './rttov_clearsky_fwd.f90' 

program fy3g_mwri_fwd_cloudy

!! Header Start
	CHARACTER*200 :: FY_FILENAME,QuerySet
	CHARACTER*200 :: FY_Folder,FY_Name,pwd,filelist,fileout,era5p
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
	!! optional, may fixed to 0,53.1,0,0
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
	REAL*4, DIMENSION(nlevel) 				  :: plevel,phalf !!hpa 	
	REAL*4 									  :: lon(nlon), lat(nlat) 
	REAL*4, DIMENSION(nlon,nlat,nlevel,ntime) :: qw,rh,ta,zg !!kg/kg, %, K, m2/s2 
	!!! hydros
	REAL*4, DIMENSION(nlon,nlat,nlevel,ntime) :: cc,ciwc,clwc,crwc,cswc
	!!ERA5 land
	REAL*4, DIMENSION(nlon,nlat,ntime) 		  :: lst,psrf_land,t2m_land,snowc,smc !! K, hPa, K, 0-1
	! REAL*4, DIMENSION(nlon,nlat,ntime)      :: stl1,swvl1,sde !! K, m3/m3,  m
	
	!! ERA5 single level
	REAL*4, DIMENSION(nlon,nlat,ntime) 		  :: u10,v10 !!m/s
	REAL*4, DIMENSION(nlon,nlat,ntime) 		  :: sst,psrf_gbl,t2m_gbl !! K, hPa, K
	
	
	!! RTTOV Vars	
	real*4,allocatable,DIMENSION(:,:,:) :: qw_3d,ta_3d	
	real*4,allocatable,DIMENSION(:,:,:) :: swc_3d,rwc_3d,lwc_3d,iwc_3d,cc_3d	
	real*4,allocatable,DIMENSION(:,:) :: skt_2d,psrf_2d,t2m_2d,u10_2d,v10_2d,snowc_2d,smc_2d
	integer :: imonth
	
	! real*4, dimension(nlevel) :: qw0,ta0
	REAL*4   :: skt
	
	!!! RTTOV OUT
	REAL*4,allocatable,DIMENSION(:,:,:) ::  TBsimu, Emiss,emissivity
	!!=================================================================================================

	plevel=(/50.0,70.0,100.0,125.0,150.0,175.0,200.0,225.0,250.0,300.0,350.0,400.0,450.0,500.0	&
        ,550.0,600.0,650.0,700.0,750.0,775.0,800.0,825.0,850.0,875.0,900.0,925.0,950.0,975.0,1000.0/)
	
	phalf=(/25.,60.,85.,112.5,137.5,162.5,187.5,212.5,237.5,275.,325.,375.,425.,475.,525.,575. &
		,625.,675.,725.,762.5,787.5,812.5,837.5,862.5,887.5,912.5,937.5,962.5,987.5/)
!! Directly specify a FY3G fullname to process,
!! e.g., /data04/0/MWRI/FY3G_IOT/ASCEND/20230626/FY3G_MWRI-_ORBA_L1_20230626_0007_7000M_V0.HDF

	! call getarg(1,FY_FILENAME)

!! Or, specify YYYYMMDD to looping process all files in a specific day
	call getarg(1,yyyymmdd)
	
	pwd="/data04/1/hjh/rttov.mlse.algor/3.cloud.retrieve/case.mlse.retrieve/"
	era5p="/data04/1/hjh/rttov.mlse.algor/3.cloud.retrieve/case.mlse.retrieve/ERA5/"
	print*,"working directory :", trim(adjustl(pwd))		
	read(yyyymmdd(5:6),'(I2)') imonth	
	
	FY_Folder="/data04/0/MWRI/FY3G_IOT/DESCEND/"//yyyymmdd//"/"
	FY_Name="FY3G_MWRI-_ORBD_L1_"//yyyymmdd//"_????_7000M_V0.HDF"
	QuerySet="cd "//trim(adjustl(FY_Folder))//"; ls "//trim(adjustl(FY_Name))
	filelist=trim(adjustl(pwd))//"/FY3G_filelist_"//yyyymmdd//".txt"

	call system(trim(adjustl(QuerySet))//"> "//trim(adjustl(filelist)))

	open(701,file=TRIM(ADJUSTL(filelist)))
702	read(701,*,end=200) FY_Name 
	FY_FILENAME=trim(adjustl(FY_Folder))//trim(adjustl(FY_Name))
	print*,"├──Start proc. - ",trim(adjustl(FY_FILENAME))
	
	fileout="descend/RTTOV_TELSEM_CLOUDY_SIMU_Swath_"//trim(adjustl(FY_Name))
	inquire(file=trim(adjustl(fileout)), exist=alive1)
    if (alive1) then 
		print*, "│	└── accomplished output detected, skip the file"	
		goto 702	
	end if

!! extract FY MWRI IOT Variables
    call get_FY_swath_nscan(FY_FILENAME,nscan,npix,retcode)
	
	allocate(FYtbs(nchn, npix, nscan)) 
	allocate(TBsimu( nchn,npix,nscan)) 
	allocate(Emiss( nchn,npix,nscan)) 
	allocate(emissivity( nchn,npix,nscan)) 
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

		
	FYstime=real(minsec)/1000./60./60.+12.  
	where(FYstime.gt.24) FYstime=FYstime-24            !! UTC	
	where(FYtbs.gt.600 .or. FYtbs.lt.30) FYtbs=-999.9  !! Kelvin
	DEALLOCATE(minsec) 	
	!! extract end
	print*,"│	├── FY3G vars read done!  "
	
	!! ERA5 Variables
	UTC_hr=nint(minVal(FYstime)*0.5+maxVal(FYstime)*0.5)
	read(FY_Name(29:32),'(I4)') fystart
	fystart_hr =fystart/100+mod(fystart,100)/60.
	if (abs(UTC_hr-fystart_hr).gt.3) UTC_hr=nint(fystart_hr)
	if (UTC_hr.eq.24) UTC_hr=0
	write(UTC_str,'(I0.2)') UTC_hr	
	
!!! download Global ERA5 reanaysis:  [65, -180, -65, 180,]
	req_era5=trim(adjustl(era5p))//yyyymmdd//"/ERA5-PL-GBL-"//yyyymmdd//"-"//UTC_str//"00.nc"
	req_land=trim(adjustl(era5p))//yyyymmdd//"/ERA5-Land-GBL-"//yyyymmdd//"-"//UTC_str//"00.nc"
	req_sigle=trim(adjustl(era5p))//yyyymmdd//"/ERA5-Single-GBL-"//yyyymmdd//"-"//UTC_str//"00.nc"

	inquire(file=trim(adjustl(req_era5)), exist=alive1)
	inquire(file=trim(adjustl(req_land)), exist=alive2)
	inquire(file=trim(adjustl(req_sigle)), exist=alive3)
    if (.not.alive1) then !如果不存在 
		print*, "│	├── ERA5-Plevel file Not Found, calling Python downloading script"
		call system("python3 "//trim(adjustl(pwd))//"ERA5/era5-pl-download.py "//yyyymmdd//" "//UTC_str )	
	end if
    if (.not.alive2) then !如果不存在
		print*, "│	├── ERA5-Land file Not Found, calling Python downloading script"	
		call system("python3 "//trim(adjustl(pwd))//"ERA5/era5-land-download.py "//yyyymmdd//" "//UTC_str )	
	end if
    if (.not.alive3) then !如果不存在 
		print*, "│	├── ERA5-Single file Not Found, calling Python downloading script"
		call system("python3 "//trim(adjustl(pwd))//"ERA5/era5-single-download.py "//yyyymmdd//" "//UTC_str )		
	end if

!!! extract ERA5 Vars
	! call read_ERA5_dims(req_era5,nlon,nlat,nlevel,ntime)
	call read_ERA5_profiles(req_era5,lon,lat,qw,rh,ta,zg)
	call read_ERA5_surfs(req_land,lst,psrf_land,t2m_land,snowc,smc)        !! K, hPa, K     !!  1 (top) -> 29 (bottom)
	call read_ERA5_single(req_sigle,sst,psrf_gbl,t2m_gbl,u10,v10)
	call read_ERA5_hydros(req_era5,cc,ciwc,clwc,crwc,cswc)
	
	where(qw.eq.0) qw=1E-9
	! print*,minval(qw)
	
	! stop
	
!! collocate RTTOV inputs
	!! profile
	allocate(qw_3d(nlevel,npix,nscan))
	allocate(ta_3d(nlevel,npix,nscan))
	!! hydro	
	allocate(cc_3d (nlevel,npix,nscan))
	allocate(iwc_3d (nlevel,npix,nscan))	
	allocate(lwc_3d (nlevel,npix,nscan))	
	allocate(rwc_3d (nlevel,npix,nscan))	
	allocate(swc_3d (nlevel,npix,nscan))	
	!! surface & 2m
	allocate(skt_2d(npix,nscan))	
	allocate(psrf_2d(npix,nscan))	
	allocate(t2m_2d(npix,nscan))	
	allocate(u10_2d(npix,nscan))	
	allocate(v10_2d(npix,nscan))	
	allocate(snowc_2d(npix,nscan))	
	allocate(smc_2d(npix,nscan))	
	
	print*,"│	├── ERA5 vars read done!  "
	do iscan=1,nscan
		do ipix=1,npix		
			ilon=nint((FYlon(ipix,iscan)-minval(lon))/grid)+1
			ilat=nint((MaxVal(lat)-FYlat(ipix,iscan))/grid)+1
			if (ilon>nlon) ilon=nlon
			if (ilat>nlat) ilat=nlat
			
			! ls=lsmask(ipix,iscan)  !!! 1 land, 2 continental water, 3 sea, 5  boundary";
			if (lsmask(ipix,iscan).eq.1) then		
				skt = lst(ilon,ilat,1)
				lsmask(ipix,iscan)= 0  
			else
				lsmask(ipix,iscan)= 1
				skt = sst(ilon,ilat,1)
			end if

			qw_3d(:,ipix,iscan) = qw(ilon,ilat,:,1)
			ta_3d(:,ipix,iscan)	= ta(ilon,ilat,:,1)
			swc_3d(:,ipix,iscan)= cswc(ilon,ilat,:,1)
			rwc_3d(:,ipix,iscan)= crwc(ilon,ilat,:,1)
			lwc_3d(:,ipix,iscan)= clwc(ilon,ilat,:,1)
			iwc_3d(:,ipix,iscan)= ciwc(ilon,ilat,:,1)
			cc_3d(:,ipix,iscan)	= cc(ilon,ilat,:,1)
			
			skt_2d(ipix,iscan)	= skt
			psrf_2d(ipix,iscan)	= psrf_gbl(ilon,ilat,1)
			t2m_2d(ipix,iscan)	= t2m_gbl(ilon,ilat,1)
			u10_2d(ipix,iscan)	= u10(ilon,ilat,1)
			v10_2d(ipix,iscan)  = v10(ilon,ilat,1)		
			snowc_2d(ipix,iscan)= snowc(ilon,ilat,1)		
			smc_2d(ipix,iscan)  = smc(ilon,ilat,1)		
		end do
	end do
	
	print*, "│	├── Start RTTOV simulations + retrieve"

	! Call Rttov 		
	call rttov_retrieve_emiss(imonth,nscan,npix,FYlon,FYlat,plevel,phalf,&
			qw_3d,ta_3d,cc_3d,iwc_3d,lwc_3d,rwc_3d,swc_3d,&
			skt_2d,psrf_2d,t2m_2d,u10_2d,v10_2d,snowc_2d,smc_2d,&
			dem,lsmask,sea,sez,soa,soz,&
			TBsimu,Emiss,FYtbs,emissivity) 
	
	! output TB simulation
	call write_hdf5_var2(trim(adjustl(fileout)),nscan,npix,nchn,&
						TBsimu,Emiss,FYtbs,emissivity,FYlon,FYlat,retcode)
	
	print*,"│	└── Write down: ", trim(adjustl(fileout))
	deallocate(qw_3d)
	deallocate(ta_3d)		
	deallocate(skt_2d)	
	deallocate(psrf_2d)	
	deallocate(t2m_2d)	
	deallocate(snowc_2d)	
	deallocate(smc_2d)	
	deallocate(u10_2d)	
	deallocate(v10_2d)
	
	deallocate(cc_3d)
	deallocate(iwc_3d)	
	deallocate(lwc_3d)	
	deallocate(rwc_3d)	
	deallocate(swc_3d)	
	
	DEALLOCATE(FYtbs) 
	DEALLOCATE(emissivity) 
	DEALLOCATE(FYlon) 
	DEALLOCATE(FYlat) 
	DEALLOCATE(dem) 
	DEALLOCATE(igbp) 
	DEALLOCATE(lsmask) 
	DEALLOCATE(FYstime) 
	DEALLOCATE(TBsimu)
	DEALLOCATE(Emiss)

	deallocate(sea)	
	deallocate(sez)	
	deallocate(soa)	
	deallocate(soz)
	! print*,"==========================================================================================================="
	goto 702
200  continue
	print*,"└── ", yyyymmdd//": process done"
end program fy3g_mwri_fwd_cloudy
