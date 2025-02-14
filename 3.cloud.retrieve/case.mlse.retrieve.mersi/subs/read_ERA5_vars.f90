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
	real*8 scalef,offsetf
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
	iovar=nf_get_att_double(ncid, varid, 'scale_factor', scalef)
	iovar=nf_get_att_double(ncid, varid, 'add_offset', offsetf)
	qw=temp*scalef+offsetf
! print*,minval(qw),maxval(qw)	
	where(qw.lt.0) qw=0.0

    ioquire=nf_inq_varid (ncid, 'r', varid)
    iovar=nf_get_var(ncid,varid,temp)
	iovar=nf_get_att_double(ncid, varid, 'scale_factor', scalef)
	iovar=nf_get_att_double(ncid, varid, 'add_offset', offsetf)
	rh=temp*scalef+offsetf 	
	where(rh.lt.0) rh=0.0
	
    ioquire=nf_inq_varid (ncid, 't', varid)
    iovar=nf_get_var(ncid,varid,temp)
	iovar=nf_get_att_double(ncid, varid, 'scale_factor', scalef)
	iovar=nf_get_att_double(ncid, varid, 'add_offset', offsetf)
	ta=temp*scalef+offsetf 		
	where(ta.lt.0) ta=180
	
    ioquire=nf_inq_varid (ncid, 'z', varid)
    iovar=nf_get_var(ncid,varid,temp)
	iovar=nf_get_att_double(ncid, varid, 'scale_factor', scalef)
	iovar=nf_get_att_double(ncid, varid, 'add_offset', offsetf)
	zg=temp*scalef+offsetf 	
	zg=zg/9.80 
	where(zg.lt.0) zg=0.0
	
	ioclose=nf_close(ncid)
END SUBROUTINE


subroutine read_ERA5_hydros(filename,cc,ciwc,clwc,crwc,cswc)
	include 'netcdf.inc'	
	
	character*(*) ::  filename
	INTEGER,Parameter :: nlon=1440, nlat=521,nlevel=29, ntime=1
	REAL*4, DIMENSION(nlon,nlat,nlevel,ntime) :: cc,ciwc,clwc,crwc,cswc
	Integer*2, DIMENSION(nlon,nlat,nlevel,ntime) :: temp  
	real*8 scalef,offsetf
	Integer iopen,ioquire,ncid,varid,ioclose
 
    iopen=nf_open(trim(adjustl(filename)),nf_nowrite,ncid) !打开netcdf文件，获取文件的ID号(ncid)      
    if (iopen .ne. 0) then
        print*,"Geolocation file, open failure!"
        stop
    end if
	
    ioquire=nf_inq_varid (ncid, 'cc', varid) 
    iovar=nf_get_var(ncid,varid,temp)
	iovar=nf_get_att_double(ncid, varid, 'scale_factor', scalef)
	iovar=nf_get_att_double(ncid, varid, 'add_offset', offsetf)
	cc=temp*scalef+offsetf 	
	where(cc.lt.0) cc=0.0
	
    ioquire=nf_inq_varid (ncid, 'ciwc', varid)
    iovar=nf_get_var(ncid,varid,temp)	
	iovar=nf_get_att_double(ncid, varid, 'scale_factor', scalef)
	iovar=nf_get_att_double(ncid, varid, 'add_offset', offsetf)
	ciwc=temp*scalef+offsetf 	
	where(ciwc.lt.0) ciwc=0.0
	
    ioquire=nf_inq_varid (ncid, 'clwc', varid)
    iovar=nf_get_var(ncid,varid,temp)	
	iovar=nf_get_att_double(ncid, varid, 'scale_factor', scalef)
	iovar=nf_get_att_double(ncid, varid, 'add_offset', offsetf)
	clwc=temp*scalef+offsetf 
	where(clwc.lt.0) clwc=0.0
	
    ioquire=nf_inq_varid (ncid, 'crwc', varid)
    iovar=nf_get_var(ncid,varid,temp)
	iovar=nf_get_att_double(ncid, varid, 'scale_factor', scalef)
	iovar=nf_get_att_double(ncid, varid, 'add_offset', offsetf)
	crwc=temp*scalef+offsetf 	
	where(crwc.lt.0) crwc=0.0	
	
    ioquire=nf_inq_varid (ncid, 'cswc', varid)
    iovar=nf_get_var(ncid,varid,temp)
	iovar=nf_get_att_double(ncid, varid, 'scale_factor', scalef)
	iovar=nf_get_att_double(ncid, varid, 'add_offset', offsetf)
	cswc=temp*scalef+offsetf 
	where(cswc.lt.0) cswc=0.0
	
	ioclose=nf_close(ncid)
END SUBROUTINE

subroutine read_ERA5_surfs(filename,skt,psrf,t2m,snowc,swvl1)
! subroutine read_ERA5_surfs(filename,skt,psrf,t2m,stl1,swvl1,snowc,sde)
	include 'netcdf.inc'	
	
	character*(*) :: filename
	INTEGER,Parameter :: nlon=1440, nlat=521,ntime=1
	REAL*4, DIMENSION(nlon,nlat,ntime) :: skt,psrf,t2m,snowc,swvl1

	! REAL*4, DIMENSION(nlon,nlat,ntime) :: stl1,swvl1,snowc,sde
	Integer*2, DIMENSION(nlon,nlat,ntime) :: temp  
	
	Integer iopen,ioquire,ncid,varid,ioclose
	real*8 scalef,offsetf
 
    iopen=nf_open(trim(adjustl(filename)),nf_nowrite,ncid) !打开netcdf文件，获取文件的ID号(ncid)      
    if (iopen .ne. 0) then
        print*,"EAR-Land file, open failure!"
        stop
    end if

    ioquire=nf_inq_varid (ncid, 'skt', varid) 
    iovar=nf_get_var(ncid,varid,temp)
	iovar=nf_get_att_double(ncid, varid, 'scale_factor', scalef)
	iovar=nf_get_att_double(ncid, varid, 'add_offset', offsetf)
	skt=temp*scalef+offsetf  !! K
	where(skt.lt.0) skt=180
	
    ioquire=nf_inq_varid (ncid, 'sp', varid)
    iovar=nf_get_var(ncid,varid,temp)	
	iovar=nf_get_att_double(ncid, varid, 'scale_factor', scalef)
	iovar=nf_get_att_double(ncid, varid, 'add_offset', offsetf)
	psrf=temp*scalef+offsetf  !! pa
	psrf=psrf/100.                                    !! hPa
	where(psrf.lt.0) psrf=0.0
	
    ioquire=nf_inq_varid (ncid, 't2m', varid)
    iovar=nf_get_var(ncid,varid,temp)
	iovar=nf_get_att_double(ncid, varid, 'scale_factor', scalef)
	iovar=nf_get_att_double(ncid, varid, 'add_offset', offsetf)
	t2m=temp*scalef+offsetf	!! K 
	where(t2m.lt.0) t2m=180

	
    ! ioquire=nf_inq_varid (ncid, 'stl1', varid)
    ! iovar=nf_get_var(ncid,varid,temp)
	! iovar=nf_get_att_double(ncid, varid, 'scale_factor', scalef)
	! iovar=nf_get_att_double(ncid, varid, 'add_offset', offsetf)
	! stl1=temp*scalef+offsetf	 !! K
	! where(stl1.lt.0) stl1=180
		
    ioquire=nf_inq_varid (ncid, 'swvl1', varid)
    iovar=nf_get_var(ncid,varid,temp)
	iovar=nf_get_att_double(ncid, varid, 'scale_factor', scalef)
	iovar=nf_get_att_double(ncid, varid, 'add_offset', offsetf)
	swvl1=temp*scalef+offsetf	
	 !! m3/m3
	where(swvl1.lt.0) swvl1=0.0
		
    ioquire=nf_inq_varid (ncid, 'snowc', varid)
    iovar=nf_get_var(ncid,varid,temp)
	iovar=nf_get_att_double(ncid, varid, 'scale_factor', scalef)
	iovar=nf_get_att_double(ncid, varid, 'add_offset', offsetf)
	snowc=temp*scalef+offsetf	 !! %
	snowc=snowc/100.0  !! 0-1
	where(snowc.lt.0) snowc=0.0
		
    ! ioquire=nf_inq_varid (ncid, 'sde', varid)
    ! iovar=nf_get_var(ncid,varid,temp)
	! iovar=nf_get_att_double(ncid, varid, 'scale_factor', scalef)
	! iovar=nf_get_att_double(ncid, varid, 'add_offset', offsetf)
	! sde=temp*scalef+offsetf !! m
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
	real*8 scalef,offsetf
	Integer iopen,ioquire,ncid,varid,ioclose
 
    iopen=nf_open(trim(adjustl(filename)),nf_nowrite,ncid) !打开netcdf文件，获取文件的ID号(ncid)      
    if (iopen .ne. 0) then
        print*,"EAR-Land file, open failure!"
        stop
    end if

    ioquire=nf_inq_varid (ncid, 'sst', varid) 
    iovar=nf_get_var(ncid,varid,temp)
		iovar=nf_get_att_double(ncid, varid, 'scale_factor', scalef)
		iovar=nf_get_att_double(ncid, varid, 'add_offset', offsetf)
	sst=temp*scalef+offsetf  !! K
	where(sst.lt.0) sst=180
	
    ioquire=nf_inq_varid (ncid, 'sp', varid)
    iovar=nf_get_var(ncid,varid,temp)
		iovar=nf_get_att_double(ncid, varid, 'scale_factor', scalef)
		iovar=nf_get_att_double(ncid, varid, 'add_offset', offsetf)
	psrf=temp*scalef+offsetf	  !! Pa
	psrf=psrf/100.                                    !! hPa
	where(psrf.lt.0) psrf=0.0
	
    ioquire=nf_inq_varid (ncid, 't2m', varid)
    iovar=nf_get_var(ncid,varid,temp)
		iovar=nf_get_att_double(ncid, varid, 'scale_factor', scalef)
		iovar=nf_get_att_double(ncid, varid, 'add_offset', offsetf)
		t2m=temp*scalef+offsetf	  !! K
	where(t2m.lt.0) t2m=180

    ioquire=nf_inq_varid (ncid, 'u10', varid)
    iovar=nf_get_var(ncid,varid,temp)
		iovar=nf_get_att_double(ncid, varid, 'scale_factor', scalef)
		iovar=nf_get_att_double(ncid, varid, 'add_offset', offsetf)
		u10=temp*scalef+offsetf   !! m/s
	where(u10.lt.-500) u10=0
	
    ioquire=nf_inq_varid (ncid, 'v10', varid)
    iovar=nf_get_var(ncid,varid,temp)
		iovar=nf_get_att_double(ncid, varid, 'scale_factor', scalef)
		iovar=nf_get_att_double(ncid, varid, 'add_offset', offsetf)
		v10=temp*scalef+offsetf   !! m/s 
	where(v10.lt.-500) v10=0
	
	ioclose=nf_close(ncid)
END SUBROUTINE

