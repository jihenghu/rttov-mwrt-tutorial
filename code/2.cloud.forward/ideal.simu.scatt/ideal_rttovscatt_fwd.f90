! Description:
!> @file
!!   Example calling RTTOV-SCATT direct model simulation.
!
!> @brief
!!   Example calling RTTOV-SCATT direct model simulation.
!!
!! @details
!!   This example is most easily run using the run_example_rttovscatt_fwd.sh
!!   script (see the user guide).
!!
!!   This program requires the following files:
!!     the file containing input profiles (e.g. prof.dat)
!!     the file containing input cloud data
!!     the RTTOV optical depth coefficient file
!!     the RTTOV hydrotable file
!!
!!   The output is written to a file called example_rttovscatt_fwd_output.dat
!!
!!   You may wish to base your own code on this example in which case you
!!   should edit in particular the sections of code marked as follows:
!!       !================================
!!       !======Read =====start===========
!!            code to be modified
!!       !======Read ===== end ===========
!!       !================================
!!
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 7 December 2016, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, DWD and MeteoFrance.
!
!    Copyright 2017, EUMETSAT, All Rights Reserved.
!
PROGRAM ideal_rttovscatt_fwd

  ! rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY :     &
         errorstatus_success, &
         errorstatus_fatal,   &
         platform_name,       &
         inst_name,           &
         hydro_index_rain,    &
         hydro_index_snow,    &
         hydro_index_graupel, &
         hydro_index_clw,     &
         hydro_index_ciw

  ! rttov_types contains definitions of all RTTOV data types
  USE rttov_types, ONLY :     &
         rttov_options,       &
         rttov_options_scatt, &
         rttov_coefs,         &
         rttov_scatt_coef,    &
         rttov_profile,       &
         rttov_profile_cloud, &
         rttov_radiance,      &
         rttov_chanprof,      &
         rttov_emissivity

  ! jpim, jprb and jplm are the RTTOV integer, real and logical KINDs
  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_unix_env, ONLY : rttov_exit

  IMPLICIT NONE

  include "rttov_scatt.interface"
  include "rttov_parallel_scatt.interface"
  include "rttov_read_scattcoeffs.interface"
  include "rttov_dealloc_scattcoeffs.interface"
  include "rttov_scatt_setupindex.interface"

  include "rttov_read_coefs.interface"
  include "rttov_dealloc_coefs.interface"
  include "rttov_alloc_direct.interface"
  include "rttov_init_emis_refl.interface"
  include "rttov_print_opts_scatt.interface"
  include "rttov_print_profile.interface"
  include "rttov_print_cld_profile.interface"
  include "rttov_skipcommentline.interface"

  !--------------------------
  !
  INTEGER(KIND=jpim), PARAMETER :: iup   = 20   ! unit for input profile file
  INTEGER(KIND=jpim), PARAMETER :: ioout = 21   ! unit for output

  INTEGER(KIND=jpim), PARAMETER :: nhydro_frac = 1 ! number of hydrometeor fractions is one, i.e. a single profile of cloud cover

  ! RTTOV variables/structures
  !====================
  TYPE(rttov_options)                :: opts                     ! Options structure - leave this set to defaults
  TYPE(rttov_options_scatt)          :: opts_scatt               ! RTTOV-SCATT options structure
  TYPE(rttov_coefs)                  :: coefs                    ! Coefficients structure
  TYPE(rttov_scatt_coef)             :: coef_scatt               ! RTTOV-SCATT coefficients structure
  TYPE(rttov_chanprof),      POINTER :: chanprof(:)    => NULL() ! Input channel/profile list
  INTEGER(KIND=jpim),        POINTER :: frequencies(:) => NULL() ! Channel indexes for hydrotable lookup
  LOGICAL(KIND=jplm),        POINTER :: use_chan(:,:)  => NULL() ! Flags to specify channels to simulate
  LOGICAL(KIND=jplm),        POINTER :: calcemis(:)    => NULL() ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity),    POINTER :: emissivity(:)  => NULL() ! Input/output surface emissivity
  TYPE(rttov_profile),       POINTER :: profiles(:)    => NULL() ! Input profiles
  TYPE(rttov_profile_cloud), POINTER :: cld_profiles(:)=> NULL() ! Input RTTOV-SCATT cloud/hydrometeor profiles
  TYPE(rttov_radiance)               :: radiance                 ! Output radiances

  INTEGER(KIND=jpim)                 :: errorstatus              ! Return error status of RTTOV subroutine calls

  INTEGER(KIND=jpim) :: alloc_status
  CHARACTER(LEN=22)  :: NameOfRoutine = 'ideal_rttovscatt_fwd'

  ! variables for input
  !====================
  CHARACTER(LEN=256) :: coef_filename
  CHARACTER(LEN=256) :: hydrotable_filename
  CHARACTER(LEN=256) :: prof_filename
  INTEGER(KIND=jpim) :: nthreads
  INTEGER(KIND=jpim) :: nlevels
  INTEGER(KIND=jpim) :: nprof
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim), ALLOCATABLE :: channel_list(:)
  ! loop variables
  INTEGER(KIND=jpim) :: j
  INTEGER(KIND=jpim) :: ilev
  INTEGER(KIND=jpim) :: iprof, joff
  INTEGER            :: ios

  !- End of header --------------------------------------------------------

  ! The usual steps to take when running RTTOV-SCATT are as follows:
  !   1. Specify required RTTOV-SCATT options
  !   2. Read coefficients and hydrotable file
  !   3. Allocate RTTOV input and output structures
  !   4. Set up the chanprof and frequencies arrays by calling rttov_scatt_setupindex
  !   5. Read input profile(s)
  !   6. Set up surface emissivity
  !   7. Call rttov_scatt and store results
  !   8. Deallocate all structures and arrays

  ! If nthreads is greater than 1 the parallel RTTOV-SCATT interface is used.
  ! To take advantage of multi-threaded execution you must have compiled
  ! RTTOV with openmp enabled. See the user guide and the compiler flags.

  errorstatus = 0_jpim

  !=====================================================
  !========== Interactive inputs == start ==============

  ! WRITE(0,*) 'enter path of coefficient file'
  ! READ(*,*) coef_filename

  ! WRITE(0,*) 'enter path of hydrotable file'
  ! READ(*,*) hydrotable_filename
  
  ! WRITE(0,*) 'enter path of file containing profile data'
  ! READ(*,*) prof_filename

  ! WRITE(0,*) 'enter number of profiles'
  ! READ(*,*) nprof
  ! WRITE(0,*) 'enter number of profile levels'
  ! READ(*,*) nlevels
  ! WRITE(0,*) 'enter number of channels to simulate per profile'
  ! READ(*,*) nchannels
  ! ALLOCATE(channel_list(nchannels))
  ! WRITE(0,*) 'enter space-separated channel list'
  ! READ(*,*,iostat=ios) channel_list(:)
  ! WRITE(0,*) 'enter number of threads to use'
  ! READ(*,*) nthreads

  prof_filename='./prof_rttovscatt.dat'
  hydrotable_filename='/home/hjh/rttov13/rtcoef_rttov13/hydrotable/hydrotable_fy3_mwri.dat'
  coef_filename='/home/hjh/rttov13/rtcoef_rttov13/rttov13pred54L/rtcoef_fy3_4_mwri.dat'
  
  WRITE(0,*) 'coefficient file : ' ,trim(adjustl(coef_filename))
  WRITE(0,*) 'hydrotable file : ' ,trim(adjustl(hydrotable_filename))
  WRITE(0,*) 'profile file : ' ,trim(adjustl(prof_filename))

  nprof=3
  nlevels=61
  nchannels=10
  ALLOCATE(channel_list(nchannels))
  channel_list=(/1,2,3,4,5,6,7,8,9,10/)
  nthreads=1
  ! --------------------------------------------------------------------------
  ! 1. Initialise RTTOV-SCATT options structure
  ! --------------------------------------------------------------------------

  ! The rttov_options structure (opts) should be left with its default values.
  ! RTTOV-SCATT only allows access to a limited number of RTTOV options: these
  ! are set in the rttov_options_scatt structure (opts_scatt).

  ! For example:
  opts_scatt % interp_mode = 1                    ! Set interpolation method
  opts_scatt % config % verbose = .TRUE.          ! Enable printing of warnings
  ! opts_scatt % lusercfrac = .TRUE.          ! Enable user input cc
	
  ! See user guide for full list of RTTOV-SCATT options

  

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

  ! Ensure input number of channels is not higher than number stored in coefficient file
  IF (nchannels > coefs % coef % fmv_chn) THEN
    nchannels = coefs % coef % fmv_chn
  ENDIF

  ! Read the RTTOV-SCATT hydrotable file
  CALL rttov_read_scattcoeffs(errorstatus, opts_scatt, coefs, coef_scatt, file_coef=hydrotable_filename)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'fatal error reading RTTOV-SCATT coefficients'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! --------------------------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! --------------------------------------------------------------------------

  ! Determine the total number of radiances to simulate (nchanprof).
  ! In this example we simulate all specified channels for each profile, but
  ! in general one can simulate a different number of channels for each profile.

  nchanprof = nchannels * nprof
 
  ! opts%rt_mw%clw_data=.true.  !!! jiheng to enable CLW absorption
 ! Allocate structures for RTTOV direct model
  CALL rttov_alloc_direct( &
        errorstatus,                 &
        1_jpim,                      &  ! 1 => allocate
        nprof,                       &
        nchanprof,                   &
        nlevels,                     &
        chanprof,                    &
        opts,                        &
        profiles,                    &
        coefs,                       &
        radiance = radiance,         &
        calcemis = calcemis,         &
        emissivity = emissivity,     &
        frequencies = frequencies,   &
        coef_scatt = coef_scatt,     &
        nhydro_frac = nhydro_frac,   &
        cld_profiles = cld_profiles, &
        init = .TRUE._jplm)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'allocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF


  ! --------------------------------------------------------------------------
  ! 4. Populate chanprof and frequencies arrays
  ! --------------------------------------------------------------------------

  ! RTTOV-SCATT requires the frequencies array to be populated by a call to
  ! rttov_scatt_setupindex. This also populates the chanprof array. To specify
  ! only a subset of channels (i.e. those in channel_list) an array of flags is
  ! passed in (use_chan).

  ! use_chan array is dimensioned by the total number of instrument channels
  ALLOCATE(use_chan(nprof,coefs%coef%fmv_chn))

  ! Set use_chan to .TRUE. only for required channels
  use_chan(:,:) = .FALSE._jplm
  DO j = 1, nprof
    use_chan(j,channel_list(1:nchannels)) = .TRUE._jplm
  ENDDO

  ! Populate chanprof and frequencies arrays
  CALL rttov_scatt_setupindex ( &
        errorstatus,        &
        nprof,              &
        coefs%coef%fmv_chn, &
        coefs,              &
        coef_scatt,         & 
        nchanprof,          &
        chanprof,           &
        frequencies,        &
        use_chan)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'error finding channels, frequencies and polarisations'
    CALL rttov_exit(errorstatus)
  ENDIF


  ! --------------------------------------------------------------------------
  ! 5. Read profile data
  ! --------------------------------------------------------------------------

  !===============================================
  !========== Read profiles == start =============

  OPEN(iup, file=TRIM(prof_filename), status='old', iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error opening profile file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF
  CALL rttov_skipcommentline(iup, errorstatus)

  ! Read gas units for profiles
  READ(iup,*) profiles(1) % gas_units
  profiles(:) % gas_units = profiles(1) % gas_units
  CALL rttov_skipcommentline(iup, errorstatus)

  ! Loop over all profiles and read data for each one
  DO iprof = 1, nprof

    ! Read vertical profile data
    ! NB The bottom-most half pressure level is taken as the 2m pressure (see below)
    !    The definitions of hydrometeors in the hydro(:,:) are governed by the
    !    hydrotable file
    DO ilev = 1, nlevels
      READ(iup,*) &
            profiles    (iprof) % p(ilev),                      &  ! full level pressure (hPa)
            cld_profiles(iprof) % ph(ilev),                     &  ! half level pressure (hPa)
            profiles    (iprof) % t(ilev),                      &  ! temperature (K)
            profiles    (iprof) % q(ilev),                      &  ! specific humidity (ppmv or kg/kg - as read above)
            cld_profiles(iprof) % hydro_frac(ilev,1),           &  ! cloud cover (0-1)
            cld_profiles(iprof) % hydro(ilev,hydro_index_clw),  &  ! cloud liquid water (kg/kg)
            cld_profiles(iprof) % hydro(ilev,hydro_index_ciw),  &  ! cloud ice water (kg/kg)
            cld_profiles(iprof) % hydro(ilev,hydro_index_rain), &  ! rain (kg/kg)
            cld_profiles(iprof) % hydro(ilev,hydro_index_snow)     ! snow (kg/kg)

	  ! profiles(iprof) % clw(ilev)= cld_profiles(iprof) % hydro(ilev,hydro_index_clw)  !!!  jiheng to enable CLW absorption
      ! graupel (kg/kg) (not in data file)
      cld_profiles(iprof) % hydro(ilev,hydro_index_graupel) = 0.0_jprb
    ENDDO
    CALL rttov_skipcommentline(iup, errorstatus)

    ! 2 meter air variables
    READ(iup,*) profiles(iprof) % s2m % t, &
                profiles(iprof) % s2m % q, &
                profiles(iprof) % s2m % p, &
                profiles(iprof) % s2m % u, &
                profiles(iprof) % s2m % v
    CALL rttov_skipcommentline(iup, errorstatus)

    ! The bottom-most half pressure level is taken as the 2m pressure
    cld_profiles(iprof) % ph(nlevels+1) = profiles(iprof) % s2m % p

    ! Skin variables
    READ(iup,*) profiles(iprof) % skin % t,        &
                profiles(iprof) % skin % salinity, &
                profiles(iprof) % skin % fastem
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Surface type and water type
    READ(iup,*) profiles(iprof) % skin % surftype, &
                profiles(iprof) % skin % watertype
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Elevation, latitude and longitude
    READ(iup,*) profiles(iprof) % elevation, &
                profiles(iprof) % latitude,  &
                profiles(iprof) % longitude
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Satellite angles
    READ(iup,*) profiles(iprof) % zenangle, &
                profiles(iprof) % azangle
    CALL rttov_skipcommentline(iup, errorstatus)

  ENDDO
  CLOSE(iup)

  !========== Read profiles == end =============
  !=============================================


  ! --------------------------------------------------------------------------
  ! 6. Specify surface emissivity
  ! --------------------------------------------------------------------------

  ! In this example we have no values for input emissivities
  ! so we initialise all inputs to zero
  CALL rttov_init_emis_refl(emissivity)

  ! Calculate emissivity within RTTOV where the input emissivity value is
  ! zero or less (all channels in this case)
  calcemis(:) = (emissivity(:) % emis_in <= 0._jprb)


  ! --------------------------------------------------------------------------
  ! 7. Call RTTOV-SCATT forward model
  ! --------------------------------------------------------------------------
  IF (nthreads <= 1) THEN
    CALL rttov_scatt ( &
        errorstatus,         &! out   error flag
        opts_scatt,          &! in    RTTOV-SCATT options structure
        nlevels,             &! in    number of profile levels
        chanprof,            &! in    channel and profile index structure
        frequencies,         &! in    channel indexes for hydrotable lookup
        profiles,            &! in    profile array
        cld_profiles,        &! in    cloud/hydrometeor profile array
        coefs,               &! in    coefficients structure
        coef_scatt,          &! in    hydrotable structure
        calcemis,            &! in    flag for internal emissivity calcs
        emissivity,          &! inout input/output emissivities per channel
        radiance)             ! inout computed radiances
  ELSE
    CALL rttov_parallel_scatt ( &
        errorstatus,         &! out   error flag
        opts_scatt,          &! in    RTTOV-SCATT options structure
        nlevels,             &! in    number of profile levels
        chanprof,            &! in    channel and profile index structure
        frequencies,         &! in    channel indexes for hydrotable lookup
        profiles,            &! in    profile array
        cld_profiles,        &! in    cloud/hydrometeor profile array
        coefs,               &! in    coefficients structure
        coef_scatt,          &! in    hydrotable structure
        calcemis,            &! in    flag for internal emissivity calcs
        emissivity,          &! inout input/output emissivities per channel
        radiance,            &! inout computed radiances
        nthreads = nthreads)  ! in    number of threads to use
  ENDIF
  IF (errorstatus /= errorstatus_success) THEN
    WRITE (*,*) 'rttov_scatt error'
    CALL rttov_exit(errorstatus)
  ENDIF


  !=====================================================
  !============== Output results == start ==============

  ! Open output file where results are written
  OPEN(ioout, file='output_'//trim(adjustl(NameOfRoutine))//'.dat', status='unknown', form='formatted', iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error opening the output file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF

  WRITE(ioout,*)' -----------------'
  WRITE(ioout,*)' Instrument ', inst_name(coefs % coef % id_inst)
  WRITE(ioout,*)' -----------------'
  WRITE(ioout,*)' '
  CALL rttov_print_opts_scatt(opts_scatt, lu=ioout)

  DO iprof = 1, nprof

    joff = (iprof-1_jpim) * nchannels

    WRITE(ioout,*)' '
    WRITE(ioout,*)' Profile ', iprof

    CALL rttov_print_profile(profiles(iprof), lu=ioout)
    CALL rttov_print_cld_profile(cld_profiles(iprof), lu=ioout)

    WRITE(ioout,777)'CHANNELS PROCESSED FOR SAT ', platform_name(coefs % coef % id_platform), coefs % coef % id_sat
    WRITE(ioout,111) (chanprof(j) % chan, j = 1+joff, nchannels+joff)
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED BRIGHTNESS TEMPERATURES (K):'
    WRITE(ioout,444) (radiance % bt(j), j = 1+joff, nchannels+joff)
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED SURFACE EMISSIVITIES:'
    WRITE(ioout,444) (emissivity(j) % emis_out, j = 1+joff, nchannels+joff)
  ENDDO

  ! Close output file
  CLOSE(ioout, iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error closing the output file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF

  !============== Output results == end ==============
  !=====================================================


  ! --------------------------------------------------------------------------
  ! 8. Deallocate all RTTOV arrays and structures
  ! --------------------------------------------------------------------------
  DEALLOCATE(channel_list, use_chan, stat=alloc_status)
  IF (alloc_status /= 0) THEN
    WRITE(*,*) 'mem dellocation error'
  ENDIF

  ! Deallocate structures for RTTOV direct model
  CALL rttov_alloc_direct( &
        errorstatus,                 &
        0_jpim,                      &  ! 0 => deallocate
        nprof,                       &
        nchanprof,                   &
        nlevels,                     &
        chanprof,                    &
        opts,                        &
        profiles,                    &
        coefs,                       &
        radiance = radiance,         &
        calcemis = calcemis,         &
        emissivity = emissivity,     &
        frequencies = frequencies,   &
        coef_scatt = coef_scatt,     &
        nhydro_frac = nhydro_frac,   &
        cld_profiles = cld_profiles)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'deallocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

  CALL rttov_dealloc_scattcoeffs(coef_scatt)

  CALL rttov_dealloc_coefs(errorstatus, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'coefs deallocation error'
  ENDIF


! Format definitions for output
111  FORMAT(1X,10I8)
444  FORMAT(1X,10F8.3)
777  FORMAT(/,A,A9,I3)

END PROGRAM ideal_rttovscatt_fwd
