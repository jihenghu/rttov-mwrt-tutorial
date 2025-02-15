!!! retrieve using RTTOV forward MW RT modeling under clear-sky
!!! &copy; Jiheng Hu 2023,August
!!! University of Science and Technology of China


!!! follow insturction of : https://cds.climate.copernicus.eu/api-how-to#install-the-cds-api-key
!!! to facilitate  ERA5 reanalysis download

!!  Ensure internet accessing, by 
!!  $ wlt > wlt.log  

	include './subs/read_FY_vars.f90' 
	include './subs/read_ERA5_vars.f90' 
	include './subs/write_hdf5_var.f90' 
	include './subs/read_FY_mersi.f90'

program fy3g_mwri_fwd_cloudy

!! Header Start
	CHARACTER*200 :: FY_FILENAME,QuerySet,QueryMERSI
	CHARACTER*200 :: MWRI_Folder,FY_Name,pwd,filelist,fileout,era5p
	CHARACTER*8   :: yyyymmdd

	!!!  FY vars
	INTEGER :: retcode
	INTEGER :: nscan, npix, iscan, ipix
	INTEGER, Parameter :: nchn=10
    REAL*4 , allocatable, DIMENSION(:,:)  :: FYlon,FYlat
	REAL*4 , allocatable, DIMENSION(:,:,:):: FYtbs,cldfrs
	INTEGER, allocatable, DIMENSION(:,:) :: igbp     
	INTEGER, allocatable, DIMENSION(:,:) :: lsmask     
	INTEGER, allocatable, DIMENSION(:,:):: dem     
	REAL*8 , allocatable, DIMENSION(:) :: minsec
	REAL*4 , allocatable, DIMENSION(:) :: fystime
	REAL*4 , allocatable, DIMENSION(:,:) :: sea,sez,soa,soz !! sensor azimuth, sensor zenith, solar azimu, solar zenith
	integer  :: ls,dem0	
	real*4  sea0,sez0,soa0,soz0,flon,flat

	!!!  mersi cloud mask
	CHARACTER*200  mersi,cfrname,mersi_Folder,mersilist,mersiname
	INTEGER nscan_clm, npix_clm
	INTEGER, allocatable,dimension(:,:):: clm
	REAL*4, allocatable,dimension(:,:):: cfr
	REAL*8, allocatable,dimension(:,:):: lon_clm,lat_clm
	
	!!! Global cloud griding
	INTEGER, Parameter :: nlo5=7201,nla5=3601
	real*4, dimension(nlo5,nla5):: cfr5km,cfr5km_org
	real*4, dimension(nlo5):: lon5km
	real*4, dimension(nla5):: lat5km
	integer, dimension(nlo5,nla5):: ncfr5km,valid
	integer ilo,ila,isc,ipi
	
	!!! ERA5 Vars
	INTEGER :: UTC_hr,fystart,MERSIstart
	REAL*4 ::fystart_hr,MERSIstart_hr
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
	REAL*4,allocatable,DIMENSION(:,:,:) ::  TBsimu, TELSEM,emissivity
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
	
	
!!!!===================================================================================================
	!! 目前的工作路径，请根据实际修改
	pwd="/data04/1/hjh/rttov.mlse.algor/3.cloud.retrieve/case.mlse.retrieve.mersi/"
	
	! ERA5的存放路径，包括现在脚本，一般不需要修改
	era5p=trim(adjustl(pwd))//"ERA5/"
	!!! MERSI CLM产品的存在路径，目前使用的是IOT测试数据集，请根据实际修改
	mersi_Folder="/data04/0/MWRI/FY3G_IOT/PRDMERSI_L2_CLM/"//yyyymmdd//"/"
	!!! MWRI数据路径，IOT测试集，请指定
	MWRI_Folder="/data04/0/MWRI/FY3G_IOT/ASCEND/"//yyyymmdd//"/"
!!!!=====================================================================================================

	print*,"working directory :", trim(adjustl(pwd))	

	! call system("mkdir "//pwd//'output')
	call system("mkdir "//pwd//'output/'//yyyymmdd)
	
	read(yyyymmdd(5:6),'(I2)') imonth	
	
	

	FY_Name="FY3G_MWRI-_ORBA_L1_"//yyyymmdd//"_????_7000M_V0.HDF"
	QuerySet="cd "//trim(adjustl(MWRI_Folder))//"; ls "//trim(adjustl(FY_Name))
	filelist=trim(adjustl(pwd))//"/FY3G_filelist_"//yyyymmdd//".txt"
	call system(trim(adjustl(QuerySet))//"> "//trim(adjustl(filelist)))

	mersi="FY3G_MERSI_GRAN_L2_CLM_MLT_NUL_"//yyyymmdd//"_????_0500M_V0.HDF"
	QueryMERSI="cd "//trim(adjustl(mersi_Folder))//"; ls "//trim(adjustl(mersi))
	mersilist=trim(adjustl(pwd))//"/MERSI_filelist_"//yyyymmdd//".txt"
	call system(trim(adjustl(QueryMERSI))//"> "//trim(adjustl(mersilist)))


  ! --------------------------------------------------------------------------
  ! 1. process MWRI obser
  ! --------------------------------------------------------------------------

	open(701,file=TRIM(ADJUSTL(filelist)))
702	read(701,*,end=200) FY_Name 
	FY_FILENAME=trim(adjustl(MWRI_Folder))//trim(adjustl(FY_Name))
	print*,"├── Start proc. - ",trim(adjustl(FY_FILENAME))
	
	fileout="output/"//yyyymmdd//"/Emissivity_Clear_FY3G_MWRIA_"//FY_Name(20:32)//"_V0.HDF"
	inquire(file=trim(adjustl(fileout)), exist=alive1)
    if (alive1) then 
		print*, "└──  Accomplished output detected, skip the file"	
		goto 702	
	end if
	! print*,trim(adjustl(fileout))
!! extract FY MWRI IOT Variables
    call get_FY_swath_nscan(FY_FILENAME,nscan,npix,retcode)
	
	allocate(FYtbs(nchn, npix, nscan),TBsimu(nchn,npix,nscan)) 
	allocate(emissivity(nchn,npix,nscan),TELSEM(nchn,npix,nscan)) 
	allocate(FYlon(npix,nscan),FYlat(npix,nscan),dem(npix,nscan))
	allocate(igbp(npix,nscan),lsmask(npix,nscan),minsec(nscan),FYstime(nscan))
	allocate(sea(npix, nscan),sez(npix, nscan),soa(npix, nscan),soz(npix, nscan))	

	call read_FY_Vars(FY_FILENAME,nscan,npix,FYlon,FYlat,FYtbs,minsec,igbp,lsmask,dem,retcode)
	call read_FY_Geometry(FY_FILENAME,nscan,npix,sea,sez,soa,soz,retcode)

		
	FYstime=real(minsec)/1000./60./60.+12.  
	where(FYstime.gt.24) FYstime=FYstime-24            !! UTC	
	where(FYtbs.gt.600 .or. FYtbs.lt.30) FYtbs=-999.9  !! Kelvin
	DEALLOCATE(minsec) 	
	!! extract end
	print*,"│	├── FY3G vars read done!  "


  ! --------------------------------------------------------------------------
  ! 2. process MERSI Cloud Fraction
  ! --------------------------------------------------------------------------

	read(FY_Name(29:32),'(I4)') fystart  !! mwri start scantime of half-orbit
	fystart_hr =fystart/100+mod(fystart,100)/60.
	
	cfr5km=0.0
	ncfr5km=0
	
	!! %% collect 5-min mersi granules for current mwri 50min-half-orbit
	print*,"│	├── Start proc. - MERSI Granules:"
	open(520,file=TRIM(ADJUSTL(mersilist)))	
516	read(520,*,end=533) mersi       !! E.G., FY3G_MERSI_GRAN_L2_CLM_MLT_NUL_20230801_0055_0500M_V0.HDF
	
	read(mersi(41:44),'(I4)') MERSIstart    !! E.G., 0055
	MERSIstart_hr =MERSIstart/100+mod(MERSIstart,100)/60.
	
	!! %% -6 min ~ 54 min
	IF ((MERSIstart_hr-fystart_hr)<-0.1) goto 516  !! 6 min prior MWRI start
	IF ((MERSIstart_hr-fystart_hr)>0.9)  goto 516  !! 1 hour post MWRI start
	! PRINT*,MERSIstart,fystart	
	
	!! %% load mersi filename
	mersiname=TRIM(ADJUSTL(mersi_Folder))//TRIM(ADJUSTL(mersi))
	print*,"│	│	├── Reading ",trim(adjustl(mersiname))
	
	!! %% open and read mersi HDF Var -- Cloud mask in 500m; Lon, Lat in 2500m
	CALL get_mersi_swath_nscan(mersiname,nscan_clm,npix_clm,retcode)	
	ALLOCATE(clm(npix_clm,nscan_clm))
	ALLOCATE(lon_clm(npix_clm/5,nscan_clm/5))
	ALLOCATE(lat_clm(npix_clm/5,nscan_clm/5))
	CALL read_mersi_clm(mersiname,nscan_clm,npix_clm,clm,lon_clm,lat_clm,retcode)

	!! %% calculate cloud fraction [valid:0-1] for each [5pix,5scan], ie, 2.5km swath grid
	ALLOCATE(cfr(npix_clm/5,nscan_clm/5))
	cfr=-999.9
	CALL calc_cfr_5grids(clm,cfr,npix_clm/5,nscan_clm/5) 
	!! %% trans to fraction in [0-1] of 2.5 km footprint
	cfr=cfr*100	
	WHERE(lon_clm.lt.-200) cfr=-999.9
	WHERE(lat_clm.lt.-200) cfr=-999.9

	!! %% WRITE OUT each 5min 2.5km resolution CFR in HDF format for checking, 
	!! %% shoud better be commented in production environ.
	! cfrname="./Cloud_fraction_"//yyyymmdd//"_"//mersi(41:44)//"_2500M_V0.HDF"
	! CALL write_real_swath_hdf(trim(cfrname),cfr,lon_clm,lat_clm,npix_clm/5,nscan_clm/5)

  ! --------------------------------------------------------------------------
  ! 2.2  Grid CFR to 5km global equal-lon-lat grid
  ! --------------------------------------------------------------------------

	!! %% Drop each sample to 5km grid and do sum
	do isc=1,nscan_clm/5
		do ipi=1,npix_clm/5
			ilo=nint((lon_clm(ipi,isc)+180.)/0.05)+1
			ila=nint((90.-lat_clm(ipi,isc))/0.05)+1	
			! print*,lon_clm(ipi,isc),lat_clm(ipi,isc)
			if(ilo.lt.0 .or. ila.lt.0) cycle
			if(cfr(ipi,isc).lt.0 .or. cfr(ipi,isc).gt.100) cycle
			cfr5km(ilo,ila)=cfr5km(ilo,ila)+cfr(ipi,isc)
			ncfr5km(ilo,ila)=ncfr5km(ilo,ila)+1				
		end do
	end do
	
	DEALLOCATE(clm)
	DEALLOCATE(cfr)
	DEALLOCATE(lon_clm)
	DEALLOCATE(lat_clm)

	goto 516
533 continue
	close(520)
	
	print*,"│	│	├── Concating MERSI granules for MWRI half-orbit"
	!! Normalized each grid with its sample nums.
	do ilo=1,nlo5
		do ila=1,nla5		
			if (ncfr5km(ilo,ila).gt.0) then 
			cfr5km(ilo,ila)=cfr5km(ilo,ila)/ncfr5km(ilo,ila)	!!! CFR cloudy fraction in percentage of 5 km grid
			else
			cfr5km(ilo,ila)=-999.9
			end if
		end do
	end do
	
	!!!%%%  WRITE OUT and Check, better be commented after debuging
	! do ilo=1,nlo5 !7201
		! lon5km(ilo)= -180.0+0.05*(ilo-1)
	! end do
	! do ila=1,nla5 !3601
		! lat5km(ila)= 90.0-0.05*(ila-1)
	! end do		
	! cfrname="./Cloud_Fraction_5km_MWRIORB_"//yyyymmdd//"_"//FY_Name(29:32)//".HDF"
	! CALL write_real_grid_hdf(trim(cfrname),cfr5km,lon5km,lat5km,nlo5,nla5)

	! stop 	
	!!! %%%  WRITE END


  ! --------------------------------------------------------------------------
  ! 2.3  Footprint average
  ! --------------------------------------------------------------------------
	valid=0
	where(cfr5km>=0 .and. cfr5km<=100) valid=1
	cfr5km_org=cfr5km
	where(cfr5km<0 .or. cfr5km>100) cfr5km=0.0
	allocate(cldfrs(nchn,npix,nscan)) 
	cldfrs=-999.9
		
  	do iscan=1,nscan
		do ipix=1,npix		
			ilo5=nint((FYlon(ipix,iscan)+180)/0.05)+1
			ila5=nint((90-FYlat(ipix,iscan))/0.05)+1
			if (ilo5>nlo5.or.ilo5<0.or.ila5>nla5.or.ila5<0) cycle
			
			!! 10.65
			cldfrs(1,ipix,iscan)=  sum(cfr5km(ilo5-2:ilo5+2,ila5-2:ila5+2))/&
									 sum(valid(ilo5-2:ilo5+2,ila5-2:ila5+2))
			cldfrs(2,ipix,iscan)=  sum(cfr5km(ilo5-2:ilo5+2,ila5-2:ila5+2))/&
									 sum(valid(ilo5-2:ilo5+2,ila5-2:ila5+2))										 
			!! 18.7/23.8/36.5
			cldfrs(3,ipix,iscan)=  sum(cfr5km(ilo5-1:ilo5+1,ila5-1:ila5+1))/&
									 sum(valid(ilo5-1:ilo5+1,ila5-1:ila5+1))
			cldfrs(4,ipix,iscan)=  sum(cfr5km(ilo5-1:ilo5+1,ila5-1:ila5+1))/&
									 sum(valid(ilo5-1:ilo5+1,ila5-1:ila5+1))
			cldfrs(5,ipix,iscan)=  sum(cfr5km(ilo5-1:ilo5+1,ila5-1:ila5+1))/&
									 sum(valid(ilo5-1:ilo5+1,ila5-1:ila5+1))
			cldfrs(6,ipix,iscan)=  sum(cfr5km(ilo5-1:ilo5+1,ila5-1:ila5+1))/&
									 sum(valid(ilo5-1:ilo5+1,ila5-1:ila5+1))
			cldfrs(7,ipix,iscan)=  sum(cfr5km(ilo5-1:ilo5+1,ila5-1:ila5+1))/&
									 sum(valid(ilo5-1:ilo5+1,ila5-1:ila5+1))
			cldfrs(8,ipix,iscan)=  sum(cfr5km(ilo5-1:ilo5+1,ila5-1:ila5+1))/&
									 sum(valid(ilo5-1:ilo5+1,ila5-1:ila5+1))						
			!! 89.0						 
			cldfrs(9,ipix,iscan)=  cfr5km_org(ilo5,ila5)											 
			cldfrs(10,ipix,iscan)= cfr5km_org(ilo5,ila5)
		end do
	end do
    	print*, "│	│	└── Calculate Footprint CFR for MWRI channels"	
		
	! call write_hdf5_var3(trim(adjustl(fileout)),nscan,npix,nchn,&
					! FYtbs,emissivity,cldfrs,lsmask,FYstime,FYlon,FYlat,retcode)
	! Stop 
  ! --------------------------------------------------------------------------
  ! 3. process ERA5 Variables
  ! --------------------------------------------------------------------------	
	
	UTC_hr=nint(minVal(FYstime)*0.5+maxVal(FYstime)*0.5)
	! read(FY_Name(29:32),'(I4)') fystart
	! fystart_hr =fystart/100+mod(fystart,100)/60.
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
		
  ! --------------------------------------------------------------------------
  ! 4. collocate RTTOV inputs
  ! --------------------------------------------------------------------------
	!! profile
	ALLOCATE(qw_3d(nlevel,npix,nscan),ta_3d(nlevel,npix,nscan))
	!! hydro	
	ALLOCATE(cc_3d (nlevel,npix,nscan),iwc_3d (nlevel,npix,nscan),lwc_3d (nlevel,npix,nscan))	
	ALLOCATE(rwc_3d (nlevel,npix,nscan),swc_3d (nlevel,npix,nscan))	
	!! surface & 2m
	ALLOCATE(skt_2d(npix,nscan),psrf_2d(npix,nscan),t2m_2d(npix,nscan))	
	ALLOCATE(u10_2d(npix,nscan),v10_2d(npix,nscan),snowc_2d(npix,nscan),smc_2d(npix,nscan))	
	
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
			swc_3d(:,ipix,iscan)= 0!cswc(ilon,ilat,:,1)
			rwc_3d(:,ipix,iscan)= 0!crwc(ilon,ilat,:,1)
			lwc_3d(:,ipix,iscan)= 0!clwc(ilon,ilat,:,1)
			iwc_3d(:,ipix,iscan)= 0!ciwc(ilon,ilat,:,1)
			cc_3d(:,ipix,iscan)	= 0!cc(ilon,ilat,:,1)
			
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


  ! --------------------------------------------------------------------------
  ! 4. ! Call Rttov Retrieve subruntn
  ! --------------------------------------------------------------------------
			
	call rttov_retrieve_emiss(imonth,nscan,npix,FYlon,FYlat,plevel,phalf,&
			qw_3d,ta_3d,cc_3d,iwc_3d,lwc_3d,rwc_3d,swc_3d,&
			skt_2d,psrf_2d,t2m_2d,u10_2d,v10_2d,snowc_2d,smc_2d,&
			dem,lsmask,sea,sez,soa,soz,&
			TBsimu,TELSEM,FYtbs,emissivity) 
	
	where(cldfrs.gt.20) emissivity=-999.9


  ! --------------------------------------------------------------------------
  ! 5. output Retrieve result
  ! --------------------------------------------------------------------------
	
	call write_hdf5_var3(trim(adjustl(fileout)),nscan,npix,nchn,&
			FYtbs,emissivity,TELSEM,cldfrs,lsmask,FYstime,FYlon,FYlat,retcode)
						
	print*,"│	└── Write down: ", trim(adjustl(fileout))

  ! --------------------------------------------------------------------------
  ! 6. free memory
  ! --------------------------------------------------------------------------	
	DEALLOCATE(qw_3d,ta_3d,skt_2d,psrf_2d,t2m_2d,snowc_2d,smc_2d,u10_2d,v10_2d)
	DEALLOCATE(cc_3d,iwc_3d,lwc_3d,rwc_3d,swc_3d)	
	DEALLOCATE(FYtbs,emissivity,FYlon,FYlat,dem,igbp,lsmask,FYstime,TBsimu,TELSEM,cldfrs)
	DEALLOCATE(sea,sez,soa,soz)
		
	goto 702
200  continue
	close(701)
	print*,"└── ", yyyymmdd//": process done"
end program fy3g_mwri_fwd_cloudy
