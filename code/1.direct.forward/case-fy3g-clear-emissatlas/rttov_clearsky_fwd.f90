
! subroutine clearsky_fwd(lon,lat,qw,ta,plevel,skt,psrf,t2m,u10,v10) 
subroutine rttov_clearsky_fwd_em_atlas(imonth,nscan,npix,lon,lat,plevel,qw,ta,skt,psrf,t2m,&
						u10,v10,snowc,smc,dem,ls,sea,sez,soa,soz,TBsimu,Emiss) 

  ! rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY :     &
         errorstatus_success, &
         errorstatus_fatal,   &
         platform_name,       &
         inst_name,           &
         surftype_sea,        &
         sensor_id_mw,        &
         sensor_id_po

  ! rttov_types contains definitions of all RTTOV data types
  USE rttov_types, ONLY :     &
         rttov_options,       &
         rttov_coefs,         &
         rttov_profile,       &
         rttov_transmission,  &
         rttov_radiance,      &
         rttov_chanprof,      &
         rttov_emissivity,    &
         rttov_reflectance

  ! The rttov_emis_atlas_data type must be imported separately
  USE mod_rttov_emis_atlas, ONLY : &
        rttov_emis_atlas_data, &
        atlas_type_ir, atlas_type_mw

  ! The rttov_brdf_atlas_data type must be imported separately
  USE mod_rttov_brdf_atlas, ONLY : rttov_brdf_atlas_data

  ! jpim, jprb and jplm are the RTTOV integer, real and logical KINDs
  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_unix_env, ONLY : rttov_exit

  IMPLICIT NONE
  
	include "rttov_direct.interface"
	include "rttov_parallel_direct.interface"
	include "rttov_read_coefs.interface"
	include "rttov_dealloc_coefs.interface"
	include "rttov_alloc_direct.interface"
	include "rttov_init_emis_refl.interface"
	include "rttov_user_options_checkinput.interface"
	include "rttov_print_opts.interface"
	include "rttov_print_profile.interface"
	include "rttov_skipcommentline.interface"

! Use emissivity atlas
	include "rttov_setup_emis_atlas.interface"
	include "rttov_get_emis.interface"
	include "rttov_deallocate_emis_atlas.interface"

! Use BRDF atlas
	include "rttov_setup_brdf_atlas.interface"
	include "rttov_get_brdf.interface"
	include "rttov_deallocate_brdf_atlas.interface"

!! ===========================================================================
  INTEGER  :: npix,nscan,iscan,ipix
  INTEGER,dimension(npix,nscan)  :: ls,dem
  real*4,dimension(npix,nscan)   :: lon,lat
  real*4,dimension(npix,nscan)   :: sea,sez,soa,soz
  real*4,dimension(npix,nscan)   :: skt,psrf,t2m,u10,v10,snowc,smc
  real*4, dimension(29) :: plevel
  real*4, dimension(29,npix,nscan) :: qw,ta
  real*4, dimension(10,npix,nscan) :: TBsimu,Emiss
!! ===========================================================================


  !--------------------------
  !
  INTEGER(KIND=jpim), PARAMETER :: iup   = 20   ! unit for input profile file
  INTEGER(KIND=jpim), PARAMETER :: ioout = 21   ! unit for output

  ! RTTOV variables/structures
  !====================
  TYPE(rttov_options)              :: opts                     ! Options structure
  TYPE(rttov_coefs)                :: coefs                    ! Coefficients structure
  TYPE(rttov_chanprof),    POINTER :: chanprof(:)    => NULL() ! Input channel/profile list
  LOGICAL(KIND=jplm),      POINTER :: calcemis(:)    => NULL() ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity),  POINTER :: emissivity(:)  => NULL() ! Input/output surface emissivity
  LOGICAL(KIND=jplm),      POINTER :: calcrefl(:)    => NULL() ! Flag to indicate calculation of BRDF within RTTOV
  TYPE(rttov_reflectance), POINTER :: reflectance(:) => NULL() ! Input/output surface BRDF
  TYPE(rttov_profile),     POINTER :: profiles(:)    => NULL() ! Input profiles
  TYPE(rttov_transmission)         :: transmission             ! Output transmittances
  TYPE(rttov_radiance)             :: radiance                 ! Output radiances

  TYPE(rttov_emis_atlas_data)      :: emis_atlas               ! Data structure for emissivity atlas
  TYPE(rttov_brdf_atlas_data)      :: brdf_atlas    

  INTEGER(KIND=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls

  INTEGER(KIND=jpim) :: atlas_type
  INTEGER(KIND=jpim) :: alloc_status
  CHARACTER(LEN=11)  :: NameOfOutput = 'clearsky_fwd_emiss_atlas'

  ! variables for input
  !====================
  CHARACTER(LEN=256) :: coef_filename
  CHARACTER(LEN=256) :: prof_filename
  INTEGER(KIND=jpim) :: nthreads
  INTEGER(KIND=jpim) :: dosolar
  INTEGER(KIND=jpim) :: imonth
  INTEGER(KIND=jpim) :: nlevels
  INTEGER(KIND=jpim) :: nprof
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim), ALLOCATABLE :: channel_list(:)
  REAL(KIND=jprb)    :: trans_out(10)
  ! loop variables
  INTEGER(KIND=jpim) :: j, jch
  INTEGER(KIND=jpim) :: np, nch
  INTEGER(KIND=jpim) :: ilev, nprint
  INTEGER(KIND=jpim) :: iprof, joff
  INTEGER            :: ios

  !- End of header --------------------------------------------------------

  ! The usual steps to take when running RTTOV are as follows:
  !   1. Specify required RTTOV options
  !   2. Read coefficients
  !   3. Allocate RTTOV input and output structures
  !   4. Set up the chanprof array with the channels/profiles to simulate
  !   5. Read input profile(s)
  !   6. Set up surface emissivity and/or reflectance
  !   7. Call rttov_direct and store results
  !   8. Deallocate all structures and arrays

  ! If nthreads is greater than 1 the parallel RTTOV interface is used.
  ! To take advantage of multi-threaded execution you must have compiled
  ! RTTOV with openmp enabled. See the user guide and the compiler flags.

  errorstatus = 0_jpim

  !=====================================================
  !========== Interactive inputs == start ==============


  coef_filename='/home/hjh/rttov13/rtcoef_rttov13/rttov13pred54L/rtcoef_fy3_4_mwri.dat'
  ! prof_filename='./prof.dat'
  nprof=npix*nscan
  nlevels=29
  dosolar=0
  nchannels=10
 
  ALLOCATE(channel_list(nchannels))
  channel_list=(/1,2,3,4,5,6,7,8,9,10/)
  nthreads=1

  ! --------------------------------------------------------------------------
  ! 1. Initialise RTTOV options structure
  ! --------------------------------------------------------------------------

  IF (dosolar == 1) THEN
    opts % rt_ir % addsolar = .TRUE.           ! Include solar radiation
  ELSE
    opts % rt_ir % addsolar = .FALSE.          ! Do not include solar radiation
  ENDIF
  opts % interpolation % addinterp   = .TRUE.  ! Allow interpolation of input profile
  opts % interpolation % interp_mode = 1       ! Set interpolation method
  opts % rt_all % addrefrac          = .TRUE.  ! Include refraction in path calc
  opts % rt_ir % addclouds           = .FALSE. ! Don't include cloud effects
  opts % rt_ir % addaerosl           = .FALSE. ! Don't include aerosol effects

  opts % rt_all % ozone_data         = .FALSE. ! Set the relevant flag to .TRUE.
  opts % rt_all % co2_data           = .FALSE. !   when supplying a profile of the
  opts % rt_all % n2o_data           = .FALSE. !   given trace gas (ensure the
  opts % rt_all % ch4_data           = .FALSE. !   coef file supports the gas)
  opts % rt_all % co_data            = .FALSE. !
  opts % rt_all % so2_data           = .FALSE. !
  opts % rt_mw % clw_data            = .FALSE. !

  opts % config % verbose            = .TRUE.  ! Enable printing of warnings
  opts % config % do_checkinput            = .FALSE.  ! Enable printing of warnings

  !========== Interactive inputs == end ==============
  !===================================================


  ! --------------------------------------------------------------------------
  ! 2. Read coefficients
  ! --------------------------------------------------------------------------
  CALL rttov_read_coefs(errorstatus, coefs, opts, file_coef=coef_filename)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'fatal error reading coefficients'
    CALL rttov_exit(errorstatus)
  ENDIF
  ! print*,coefs % coef % fmv_chn,nchannels

  ! Ensure input number of channels is not higher than number stored in coefficient file
  IF (nchannels > coefs % coef % fmv_chn) THEN
    nchannels = coefs % coef % fmv_chn
  ENDIF

  ! Ensure the options and coefficients are consistent
  CALL rttov_user_options_checkinput(errorstatus, opts, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'error in rttov options'
    CALL rttov_exit(errorstatus)
  ENDIF


  ! --------------------------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! --------------------------------------------------------------------------

  ! Determine the total number of radiances to simulate (nchanprof).
  ! In this example we simulate all specified channels for each profile, but
  ! in general one can simulate a different number of channels for each profile.

  nchanprof = nchannels * nprof

  ! Allocate structures for rttov_direct
  CALL rttov_alloc_direct( &
        errorstatus,             &
        1_jpim,                  &  ! 1 => allocate
        nprof,                   &
        nchanprof,               &
        nlevels,                 &
        chanprof,                &
        opts,                    &
        profiles,                &
        coefs,                   &
        transmission,            &
        radiance,                &
        calcemis=calcemis,       &
        emissivity=emissivity,   &
        calcrefl=calcrefl,       &
        reflectance=reflectance, &
        init=.TRUE._jplm)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'allocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF


  ! Initialise the RTTOV emissivity atlas
  ! (this loads the default IR/MW atlases: use the atlas_id argument to select alternative atlases)
  IF (coefs%coef%id_sensor == sensor_id_mw .OR. &
      coefs%coef%id_sensor == sensor_id_po) THEN
    atlas_type = atlas_type_mw ! MW atlas
  ELSE
    atlas_type = atlas_type_ir ! IR atlas
  ENDIF
  ! print*,imonth
  CALL rttov_setup_emis_atlas(          &
              errorstatus,              &
              opts,                     &
              imonth,                   &
              atlas_type,               & ! Selects MW (1) or IR (2)
              emis_atlas,               &
              path = '/home/hjh/rttov13/emis_data', & ! The default path to atlas data
              coefs = coefs) ! This is mandatory for the CNRM MW atlas, ignored by TELSEM2;
                             ! if supplied for IR atlases they are initialised for this sensor
                             ! and this makes the atlas much faster to access
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'error initialising emissivity atlas'
    CALL rttov_exit(errorstatus)
  ENDIF

  IF (opts % rt_ir % addsolar) THEN

    ! Initialise the RTTOV BRDF atlas
    CALL rttov_setup_brdf_atlas(        &
                errorstatus,            &
                opts,                   &
                imonth,                 &
                brdf_atlas,             &
                path='/home/hjh/rttov13/brdf_data', &  ! The default path to atlas data
                coefs = coefs) ! If supplied the BRDF atlas is initialised for this sensor and
                               ! this makes the atlas much faster to access
    IF (errorstatus /= errorstatus_success) THEN
      WRITE(*,*) 'error initialising BRDF atlas'
      CALL rttov_exit(errorstatus)
    ENDIF

  ENDIF


  ! --------------------------------------------------------------------------
  ! 4. Build the list of profile/channel indices in chanprof
  ! --------------------------------------------------------------------------

  nch = 0_jpim
  DO j = 1, nprof
    DO jch = 1, nchannels
      nch = nch + 1_jpim
      chanprof(nch)%prof = j
      chanprof(nch)%chan = channel_list(jch)
    ENDDO
  ENDDO


  ! --------------------------------------------------------------------------
  ! 5. Read profile data
  ! --------------------------------------------------------------------------

  !===============================================
  !========== Read profiles == start =============


  ! Read gas units for profiles
	! Gas units (must be same for all profiles)
	! 0 => ppmv over dry air
	! 1 => kg/kg over moist air
	! 2 => ppmv over moist air
  profiles(:) % gas_units = 1
 
  ! Loop over all profiles and read data for each one
  DO iscan = 1, nscan
  DO ipix = 1, npix
	iprof=(iscan-1)*npix+ipix
    ! Read pressure (hPa), temp (K), WV, O3 (gas units ppmv or kg/kg - as read above)
    profiles(iprof) % p(:)= plevel
    profiles(iprof) % t(:)= ta(:,ipix,iscan)
    profiles(iprof) % q(:)= qw(:,ipix,iscan)
	! print*,profiles(iprof) % q(:)

    ! 2 meter air variables
    profiles(iprof) % s2m % t= t2m(ipix,iscan)
    profiles(iprof) % s2m % q= qw(nlevels,ipix,iscan)
    profiles(iprof) % s2m % p= psrf(ipix,iscan)
    profiles(iprof) % s2m % u= u10(ipix,iscan) 
    profiles(iprof) % s2m % v= v10(ipix,iscan)
    profiles(iprof) % s2m % wfetc=0


    ! Skin variables
    profiles(iprof) % skin % t= skt(ipix,iscan)       
    profiles(iprof) % skin % salinity= 35.0
    profiles(iprof) % skin % snow_fraction=snowc(ipix,iscan)
    profiles(iprof) % skin % soil_moisture=smc(ipix,iscan)
    profiles(iprof) % skin % fastem= (/3.0, 5.0, 15.0, 0.1, 0.3/)  
   
    ! Surface type and water type
    profiles(iprof) % skin % surftype= ls(ipix,iscan)
    profiles(iprof) % skin % watertype= 1


    ! Elevation, latitude and longitude
    profiles(iprof) % elevation= dem(ipix,iscan)/1000.
    profiles(iprof) % latitude= lat(ipix,iscan)
    profiles(iprof) % longitude= lon(ipix,iscan)

    ! Satellite and solar angles
    profiles(iprof) % zenangle= sez(ipix,iscan) !53.1
    profiles(iprof) % azangle = sea(ipix,iscan)
    profiles(iprof) % sunzenangle= soz(ipix,iscan)
    profiles(iprof) % sunazangle = soa(ipix,iscan)
  
    ! Cloud variables for simple cloud scheme, set cfraction to 0. to turn this off (VIS/IR only)
    profiles(iprof) % ctp = 500
    profiles(iprof) % cfraction= 0.0
   
  ENDDO
  ENDDO
  ! CLOSE(iup)

  !========== Read profiles == end =============
  !=============================================


  ! --------------------------------------------------------------------------
  ! 6. Specify surface emissivity and reflectance
  ! --------------------------------------------------------------------------
  ! Use emissivity atlas
  CALL rttov_get_emis(             &
            errorstatus,           &
            opts,                  &
            chanprof,              &
            profiles,              &
            coefs,                 &
            emis_atlas,            &
            emissivity(:) % emis_in)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'error reading emissivity atlas'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! Calculate emissivity within RTTOV where the atlas emissivity value is
  ! zero or less
  calcemis(:) = (emissivity(:) % emis_in <= 0._jprb)

  IF (opts % rt_ir % addsolar) THEN

    ! Use BRDF atlas
    CALL rttov_get_brdf(              &
              errorstatus,            &
              opts,                   &
              chanprof,               &
              profiles,               &
              coefs,                  &
              brdf_atlas,             &
              reflectance(:) % refl_in)
    IF (errorstatus /= errorstatus_success) THEN
      WRITE(*,*) 'error reading BRDF atlas'
      CALL rttov_exit(errorstatus)
    ENDIF

    ! Calculate BRDF within RTTOV where the atlas BRDF value is zero or less
    calcrefl(:) = (reflectance(:) % refl_in <= 0._jprb)

  ENDIF

  ! Use the RTTOV emissivity and BRDF calculations over sea surfaces
  DO j = 1, SIZE(chanprof)
    IF (profiles(chanprof(j)%prof) % skin % surftype == surftype_sea) THEN
      calcemis(j) = .TRUE.
      calcrefl(j) = .TRUE.
    ENDIF
  ENDDO

  ! Use default cloud top BRDF for simple cloud in VIS/NIR channels
  reflectance(:) % refl_cloud_top = 0._jprb

  ! Let RTTOV provide diffuse surface reflectances
  reflectance(:) % diffuse_refl_in = 0._jprb


  ! --------------------------------------------------------------------------
  ! 7. Call RTTOV forward model
  ! --------------------------------------------------------------------------
  IF (nthreads <= 1) THEN
    CALL rttov_direct(                &
            errorstatus,              &! out   error flag
            chanprof,                 &! in    channel and profile index structure
            opts,                     &! in    options structure
            profiles,                 &! in    profile array
            coefs,                    &! in    coefficients structure
            transmission,             &! inout computed transmittances
            radiance,                 &! inout computed radiances
            calcemis    = calcemis,   &! in    flag for internal emissivity calcs
            emissivity  = emissivity, &! inout input/output emissivities per channel
            calcrefl    = calcrefl,   &! in    flag for internal BRDF calcs
            reflectance = reflectance) ! inout input/output BRDFs per channel
  ELSE
    CALL rttov_parallel_direct(     &
            errorstatus,              &! out   error flag
            chanprof,                 &! in    channel and profile index structure
            opts,                     &! in    options structure
            profiles,                 &! in    profile array
            coefs,                    &! in    coefficients structure
            transmission,             &! inout computed transmittances
            radiance,                 &! inout computed radiances
            calcemis    = calcemis,   &! in    flag for internal emissivity calcs
            emissivity  = emissivity, &! inout input/output emissivities per channel
            calcrefl    = calcrefl,   &! in    flag for internal BRDF calcs
            reflectance = reflectance,&! inout input/output BRDFs per channel
            nthreads    = nthreads)    ! in    number of threads to use
  ENDIF

  IF (errorstatus /= errorstatus_success) THEN
    WRITE (*,*) 'rttov_direct error'
    CALL rttov_exit(errorstatus)
  ENDIF

  !=====================================================
  !============== Output results == start ==============
  DO iscan = 1, nscan
  DO ipix = 1, npix
	iprof=(iscan-1)*npix+ipix
    joff = (iprof-1_jpim) * nchannels
    TBsimu(:,ipix,iscan)=radiance%bt((1+joff):nchannels+joff)
    Emiss(:,ipix,iscan)=emissivity((1+joff):nchannels+joff) % emis_out
  ENDDO
  ENDDO

  !============== Output results == end ==============
  !=====================================================


  ! --------------------------------------------------------------------------
  ! 8. Deallocate all RTTOV arrays and structures
  ! --------------------------------------------------------------------------
  DEALLOCATE (channel_list, stat=alloc_status)
  IF (alloc_status /= 0) THEN
    WRITE(*,*) 'mem dellocation error'
  ENDIF

  ! Deallocate structures for rttov_direct
  CALL rttov_alloc_direct( &
        errorstatus,             &
        0_jpim,                  &  ! 0 => deallocate
        nprof,                   &
        nchanprof,               &
        nlevels,                 &
        chanprof,                &
        opts,                    &
        profiles,                &
        coefs,                   &
        transmission,            &
        radiance,                &
        calcemis=calcemis,       &
        emissivity=emissivity,   &
        calcrefl=calcrefl,       &
        reflectance=reflectance)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'deallocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

  CALL rttov_dealloc_coefs(errorstatus, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'coefs deallocation error'
  ENDIF


! Format definitions for output
111  FORMAT(1X,10I8)
1115 FORMAT(3X,10I8)
222  FORMAT(1X,10F8.2)
444  FORMAT(1X,10F8.3)
4444 FORMAT(1X,10F8.4)
4445 FORMAT(1X,I2,10F8.4)
777  FORMAT(/,A,A9,I3)

END subroutine rttov_clearsky_fwd_em_atlas
