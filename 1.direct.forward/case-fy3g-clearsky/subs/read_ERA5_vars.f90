subroutine read_ERA5_dims(filename,nlon,nlat,nlevel,ntime)
	include 'netcdf.inc'

	character*(*) ::filename
	integer nlon, nlat, nlevel, ntime
	Integer iopen,ioquire,ncid,ioclose
	Integer lon_dimid,lat_dimid,level_dimid,time_dimid
	
	! print*,filename
	iopen=nf_open(trim(adjustl(filename)),nf_nowrite,ncid)
    if (iopen .ne. 0) then
        print*,"ERA5 open failure!"
        stop
    end if
	
	ioquire=nf_inq_dimid (ncid, 'longitude', lon_dimid) 
	ioquire = nf_inq_dimlen(ncid,lon_dimid,nlon)
	
	ioquire=nf_inq_dimid (ncid, 'latitude', lat_dimid) 
	ioquire = nf_inq_dimlen(ncid,lat_dimid,nlat)	
	
	! nlevel=29
	! ntime=1
	ioquire=nf_inq_dimid (ncid, 'level', level_dimid) 
	ioquire = nf_inq_dimlen(ncid,level_dimid,nlevel)
	
	ioquire=nf_inq_dimid (ncid, 'time', time_dimid) 
	ioquire = nf_inq_dimlen(ncid,time_dimid,ntime)
	
	ioclose=nf_close(ncid)
end subroutine

subroutine read_ERA5_profiles(filename,lon,lat,qw,rh,ta,zg)
	include 'netcdf.inc'	
	
	character*(*) ::  filename
	INTEGER,Parameter :: nlon=1440, nlat=521,nlevel=29, ntime=1
	REAL*4 :: lon(nlon), lat(nlat) 
	REAL*4, DIMENSION(nlon,nlat,nlevel,ntime) :: qw,rh,ta,zg !!kg/kg, %, K, m2/s2 
	Integer*2, DIMENSION(nlon,nlat,nlevel,ntime) :: temp  
	
	Integer iopen,ioquire,ncid,varid,ioclose
 
    iopen=nf_open(trim(adjustl(filename)),nf_nowrite,ncid) !打开netcdf文件，获取文件的ID号(ncid)      
    if (iopen .ne. 0) then
        print*,"Geolocation file, open failure!"
        stop
    end if
    ioquire=nf_inq_varid (ncid, 'longitude', varid) 
    iovar=nf_get_var_real(ncid,varid,lon)
	
    ioquire=nf_inq_varid (ncid, 'latitude', varid) 
    iovar=nf_get_var_real(ncid,varid,lat)	

    ioquire=nf_inq_varid (ncid, 'q', varid) 
    iovar=nf_get_var(ncid,varid,temp)
	qw=temp*3.919505253913927E-7+0.012842628851428717
	where(qw.lt.0) qw=0.0
	! print*, MinVal(qw),MaxVal(qw)	
	
    ioquire=nf_inq_varid (ncid, 'r', varid)
    iovar=nf_get_var(ncid,varid,temp)	
	rh=temp*0.0026103162240518557+73.12762287738563	
	where(rh.lt.0) rh=0.0
	
    ioquire=nf_inq_varid (ncid, 't', varid)
    iovar=nf_get_var(ncid,varid,temp)	
	ta=temp*0.0020127642726428144+254.27396462616446
	where(ta.lt.0) ta=180
	
    ioquire=nf_inq_varid (ncid, 'z', varid)
    iovar=nf_get_var(ncid,varid,temp)
	! print*, minVal(temp), maxVal(temp)
	zg=temp*3.2108395688109903+101238.54484877028
	zg=zg/9.80 
	where(zg.lt.0) zg=0.0
	
	ioclose=nf_close(ncid)
END SUBROUTINE

subroutine read_ERA5_surfs(filename,skt,psrf,t2m)
! subroutine read_ERA5_surfs(filename,skt,psrf,t2m,stl1,swvl1,snowc,sde)
	include 'netcdf.inc'	
	
	character*(*) :: filename
	INTEGER,Parameter :: nlon=1440, nlat=521,ntime=1
	REAL*4, DIMENSION(nlon,nlat,ntime) :: skt,psrf,t2m

	! REAL*4, DIMENSION(nlon,nlat,ntime) :: stl1,swvl1,snowc,sde
	Integer*2, DIMENSION(nlon,nlat,ntime) :: temp  
	
	Integer iopen,ioquire,ncid,varid,ioclose
 
    iopen=nf_open(trim(adjustl(filename)),nf_nowrite,ncid) !打开netcdf文件，获取文件的ID号(ncid)      
    if (iopen .ne. 0) then
        print*,"EAR-Land file, open failure!"
        stop
    end if

    ioquire=nf_inq_varid (ncid, 'skt', varid) 
    iovar=nf_get_var(ncid,varid,temp)
	skt=temp*0.0013377682808661286+282.0693842164455  !! K
	where(skt.lt.0) skt=180
	
    ioquire=nf_inq_varid (ncid, 'sp', varid)
    iovar=nf_get_var(ncid,varid,temp)	
	psrf=temp*0.8168708894755314+76597.46265830526	  !! Pa
	psrf=psrf/100.                                    !! hPa
	where(psrf.lt.0) psrf=0.0
	
    ioquire=nf_inq_varid (ncid, 't2m', varid)
    iovar=nf_get_var(ncid,varid,temp)	
	t2m=temp*0.0011038391440190438+281.6018803314046   !! K
	where(t2m.lt.0) t2m=180

	
    ! ioquire=nf_inq_varid (ncid, 'stl1', varid)
    ! iovar=nf_get_var(ncid,varid,temp)
	! stl1=temp*0.0011601681118673034+285.8065945985613   !! K
	! where(stl1.lt.0) stl1=180
		
    ! ioquire=nf_inq_varid (ncid, 'swvl1', varid)
    ! iovar=nf_get_var(ncid,varid,temp)
	! swvl1=temp* 1.1673498689949148E-5+ 0.3824938580748738  !! m3/m3
	! where(swvl1.lt.0) swvl1=0.0
		
    ! ioquire=nf_inq_varid (ncid, 'snowc', varid)
    ! iovar=nf_get_var(ncid,varid,temp)
	! snowc=temp*0.0015259487586406848+49.99923702562068  !! %
	! where(snowc.lt.0) snowc=0.0
		
    ! ioquire=nf_inq_varid (ncid, 'sde', varid)
    ! iovar=nf_get_var(ncid,varid,temp)
	! sde=temp*5.086446189324462E-4+ 16.666249583940534  !! m
	! where(sde.lt.0) sde=0.0
	
	ioclose=nf_close(ncid)
END SUBROUTINE

subroutine read_ERA5_single(filename,sst,psrf,t2m,u10,v10)
	include 'netcdf.inc'	
	
	character*(*) :: filename
	INTEGER,Parameter :: nlon=1440, nlat=521,ntime=1
	REAL*4, DIMENSION(nlon,nlat,ntime) :: sst,psrf,t2m
	REAL*4, DIMENSION(nlon,nlat,ntime) :: u10,v10
	Integer*2, DIMENSION(nlon,nlat,ntime) :: temp  
	
	Integer iopen,ioquire,ncid,varid,ioclose
 
    iopen=nf_open(trim(adjustl(filename)),nf_nowrite,ncid) !打开netcdf文件，获取文件的ID号(ncid)      
    if (iopen .ne. 0) then
        print*,"EAR-Land file, open failure!"
        stop
    end if

    ioquire=nf_inq_varid (ncid, 'sst', varid) 
    iovar=nf_get_var(ncid,varid,temp)
	sst=temp*5.667916116689301E-4+288.8930271510692  !! K
	where(sst.lt.0) sst=180
	
    ioquire=nf_inq_varid (ncid, 'sp', varid)
    iovar=nf_get_var(ncid,varid,temp)	
	psrf=temp*0.8184426166969313+76686.18452869165	  !! Pa
	psrf=psrf/100.                                    !! hPa
	where(psrf.lt.0) psrf=0.0
	
    ioquire=nf_inq_varid (ncid, 't2m', varid)
    iovar=nf_get_var(ncid,varid,temp)	
	t2m=temp*0.0011056273652205759+282.37044938358304   !! K
	where(t2m.lt.0) t2m=180

    ioquire=nf_inq_varid (ncid, 'u10', varid)
    iovar=nf_get_var(ncid,varid,temp)	
	u10=temp*5.828856024827187E-4+1.1313919068081337   !! m/s
	where(u10.lt.-500) u10=0
	
    ioquire=nf_inq_varid (ncid, 'v10', varid)
    iovar=nf_get_var(ncid,varid,temp)	
	v10=temp*5.233527383150474E-4+0.06637345546677997   !! m/s
	where(v10.lt.-500) v10=0
	
	ioclose=nf_close(ncid)
END SUBROUTINE

! subroutine read_ERA5_hydros(filename,lon,lat,cc,ciwc,clwc,qw,rh,ta,z)
	! include 'netcdf.inc'	
	
	! character  filename*100
	! Real*4 longitude(11104),latitude(4320)
	! Real*4 longitude(11104,4320),latitude(11104,4320)

	! character easeGeoFile*90
	! Integer iopen,ioquire,iovar,latid,lonid,ncid,ioclose
  
	! easeGeoFile="./subs/EASE2_T3.125km.geolocation.single.nc"

    ! iopen=nf_open(trim(adjustl(easeGeoFile)),nf_nowrite,ncid) !打开netcdf文件，获取文件的ID号(ncid)      
    ! if (iopen .ne. 0) then
        ! print*,"Geolocation file, open failure!"
        ! stop
    ! end if

    ! ioquire=nf_inq_varid (ncid, 'latitude', latid) 
    ! ioquire=ioquire+nf_inq_varid (ncid, 'longitude', lonid) 
    ! if (ioquire .ne. 0) then
        ! print*,"Geolocation extraction: variables not found!"
        ! stop
    ! end if

    ! iovar=nf_get_var_real(ncid,latid,latitude)
    ! iovar=iovar+nf_get_var_real(ncid,lonid,longitude)
    ! if (iovar .ne. 0) then
        ! print*,"Geolocation extraction: variables extract failed!"
        ! stop
    ! end if	
		! ioclose=nf_close(ncid)
! END SUBROUTINE
