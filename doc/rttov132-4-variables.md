# RTTOV模式学习笔记：(四) RTTOV 变量和结构体

&copy;Jiheng Hu 2023-2030, 禁止转载。

理解和熟悉RTTOV的工程结构和变量含义是我们使用模型的第一步。RTTOV的示例和源码可以帮助我们快速搭建起自己的模拟项目。

## RTTOV工程结构
这里介绍一下最重要的源码和示例文件，这是使用RTTOV的基础。
```bash
$ tree 
|- bin/  ## 存放编译产生的可执行文件，其源码来自于/src/test/下的代码示例
|- rtcoef_rttov13/ ##存放大气吸收系数和散射系数的文件夹
|- rttov_test/     ##测试脚本，调用bin/下的可执行文件和预置的相关大气廓线，水凝物和大气系数以完成模式安装的测试；
|- src/  ##源码路径
  |- main/  RTTOV的F90源码路径包括 正向模拟的接口，
  |     |- rttov_k.F90    ## 雅可比矩阵计算
  |     |- rttov_ad.F90   ## 伴随矩阵计算
  |     |- rttov_tl.F90   ## tangent linear计算
  |     |- rttov_direct.F90  ## 正向模拟-晴空
  |     |- rttov_const.F90 ## 常量定义
  |     |- rttov_types.F90 ## 类型定义
  |- mw_scatt/  ## 微波散射模型-主要用于云、降水情形的模拟 
  |     |- rttov_scatt.F90 
  |     |- rttov_scatt_ad.F90 
  |     |- rttov_scatt_tl.F90 
  |     |- rttov_scatt_emis_retrieval.F90 ## Computes all-sky retrieved emissivities (Baordo and Geer, 2016, DOI:10.1002/qj.2873)
  |- test/   ## 模式调用示例
  |     |- example_fwd.F90   			## 调用模式进行正向传输的示例，晴空
  |     |- example_rttovscatt_fwd.F90   ## 调用模式进行正向传输的示例，水凝物
  |     |- Makefile_examples  ##  Makefile example，生成bin/下的可执行文件，按照rttov_test/ 下的shell脚本示例来运行。 
```

## RTTOV变量及其结构
### 类型定义
模式的运行需要提前声明和创建一系列变量来作为输入量和输出量。这些变量和结构体在使用时需要先声明其类型，形式为：

```fortran  src/test/example_fwd.F90
  ! RTTOV variables/structures
  !====================
  TYPE(rttov_options)              :: opts                     ! Options structure
  TYPE(rttov_coefs)                :: coefs                    ! Coefficients structure
  TYPE(rttov_chanprof),    POINTER :: chanprof(:)    ! Input channel/profile list
  LOGICAL(KIND=jplm),      POINTER :: calcemis(:)    ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity),  POINTER :: emissivity(:)  ! Input/output surface emissivity
  LOGICAL(KIND=jplm),      POINTER :: calcrefl(:)    ! Flag to indicate calculation of BRDF within RTTOV
  TYPE(rttov_reflectance), POINTER :: reflectance(:) ! Input/output surface BRDF
  TYPE(rttov_profile),     POINTER :: profiles(:)    ! Input profiles
  TYPE(rttov_transmission)         :: transmission             ! Output transmittances
  TYPE(rttov_radiance)             :: radiance                 ! Output radiances

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
```

这些TYPE在`rttov_types.F90`中被定义，我们特别看一下这些结构体的内容。
首先，`nchanels`和`nprof`分别为模拟的通道数和廓线条数，`nchanprof = nchannels * nprof` 指定了进行模拟的次数。

#### rttov_options RT方案、插值开关

```fortran 
  !> RTTOV options structure
  TYPE rttov_options
    TYPE(rttov_opts_config)  :: config          !< General configuration options
    TYPE(rttov_opts_rt_all)  :: rt_all          !< General RT options
    TYPE(rttov_opts_rt_ir)   :: rt_ir           !< VIS/IR RT options
    TYPE(rttov_opts_rt_mw)   :: rt_mw           !< MW RT options
    TYPE(rttov_opts_interp)  :: interpolation   !< Interpolation options
    TYPE(rttov_opts_htfrtc)  :: htfrtc_opts     !< HTFRTC options
    TYPE(rttov_opts_dev)     :: dev             !< Developer-only options
  END TYPE
```

其中：

- `rt_all`结构体设定了RT的算法大气和地表设定方案开关；(参见 TYPE rttov_opts_rt_all的定义)；

- `rt_mw`结构体设定了MWRT的算法选项，如算法版本，云参数化方案等(参见 TYPE rttov_opts_rt_mw的定义)；

```fortran
  !> MW-only radiative transfer options
  TYPE rttov_opts_rt_mw
    ......
    LOGICAL(jplm) :: clw_data   = .FALSE.   !< Switch to enable input of cloud liquid water profile
    INTEGER(jpim) :: clw_scheme = 2         !< MW CLW scheme: 1 => Liebe, 2 => Rosenkranz, 3 => TKC
    ......
  END TYPE
```

默认的MW模拟时不考虑cloud liquid water的吸收作用的。

- `interpolation`结构体定义了廓线插值的设定，包括插值开关`opts%interpolation%addinterp = .True.` (参见rttov_opts_interp的定义)；

#### rttov_chanprof 模拟量
```fortran rttov_types.F90
  TYPE rttov_chanprof
    INTEGER(jpim) :: chan              !< Channel index
    INTEGER(jpim) :: prof              !< Profile index
  END TYPE
```
包含两个一维integer数组变量：`chan`,`prof`，其size决定于nchanels*nprof，也就是说，对于每个通道和每条廓线，都要定义其channe index和profile index，注意，这里的变量存储的应该是RTTOV内置的卫星表对应的通道表中的通道序号，而不是实际的频率。
赋值示例：
```fortran  example_fwd.F90
  nch = 0_jpim
  DO j = 1, nprof
    DO jch = 1, nchannels
      nch = nch + 1_jpim
      chanprof(nch)%prof = j
      chanprof(nch)%chan = channel_list(jch)
    ENDDO
  ENDDO
```

#### rttov_coefs 大气参数

coefs是用来读取和缓存RTTOV的核心--internal Optical depth coefs的变量类型。
```fortran
!> @internal Optical depth coefs
TYPE rttov_coef
    ! Structure for the storage of RTTOV coefficients
    ! this may differ from what is stored in the coefficient files especially
    ! for the units (ie kg/kg to ppmv)
    ! Gases are separated in MxG WV O3
    ! Number of levels is the same for all gases (taken from MxG).
    INTEGER(jpim)       :: id_platform   ! platform   (see user guide or rttov_const)
    INTEGER(jpim)       :: id_sat        ! satellite  (.....)
    INTEGER(jpim)       :: id_inst       ! instrument (.....)
    INTEGER(jpim)       :: id_sensor     ! 1 = Infrared, 2 = Microwave, 3 = High-resolution, 4 = Polarimeter

    ......

    !FAST_COEFFICIENTS/SOLAR_FAST_COEFFICIENTS section
    ! For non-PW instruments, "solar" will point to "thermal" coefs
    ! For instruments with solar-affected PW channels, both thermal and solar
    ! structures will be populated from coef file
    ! Fast coefficients now grouped per channel rather than per gas (RTTOV12)
    ! so thermal and solar are now arrays with size nchannels.
    TYPE(rttov_fast_coef), POINTER :: thermal(:)      ! FAST_COEFFICIENTS for gas opdep prediction
    TYPE(rttov_fast_coef), POINTER :: thermal_corr(:) ! Coefficients for correction term
``` 
`id_platform`,`id_sat`,`id_inst`确定了用来读取的快速辐射参数的文件名，这些ID见下一篇[《RTTOV辐射传输模式实践：RTTOV V13.2 约定和方案》](../rttov132-conventions)的表格介绍，当然了你可直接看文档。
`rttov_fast_coef`类型用来读取和缓存多种气体的快速辐射传输FAST RT系数廓线：
```fortran
  !> @internal Actual storage for gas fast coefficients 
  TYPE rttov_fast_coef_gas
    REAL(jprb), POINTER :: coef(:,:) => NULL()
  END TYPE rttov_fast_coef_gas

  !> @internal Fast coefficients
  !  Separate structure to allow for non-PW and PW if necessary
  TYPE rttov_fast_coef
    ! SIZE is fmv_gas, i.e. one element per gas.
    TYPE(rttov_fast_coef_gas), POINTER  :: gasarray(:) => NULL()

    ! SHAPE is (ncoef, levels): these point to coefs within gasarray for ease of use
    REAL(jprb), POINTER :: mixedgas(:,:)    => NULL()
    REAL(jprb), POINTER :: watervapour(:,:) => NULL()
    REAL(jprb), POINTER :: ozone(:,:)       => NULL()
    REAL(jprb), POINTER :: wvcont(:,:)      => NULL()
    REAL(jprb), POINTER :: co2(:,:)         => NULL()
    REAL(jprb), POINTER :: n2o(:,:)         => NULL()
    REAL(jprb), POINTER :: co(:,:)          => NULL()
    REAL(jprb), POINTER :: ch4(:,:)         => NULL()
    REAL(jprb), POINTER :: so2(:,:)         => NULL()
  END TYPE rttov_fast_coef
```
所以，到目前为止，你应该明白了，RTTOV的全部模式中是不存在具体的频率信息作为输入的，对于不同仪器不同通道的模拟体现在上述coefs大气散射和消光参数的设定。这些参数文件可以在官网下载，其使用辐射传输和机器学习方法提前计算的大气光学参数数据集，所有的散射和吸收计算已经在应用先行数据集中进行了电磁学计算（介电常数计算，Mie, Reyleigh, Lieb等方案）。那么到RTTOV中，只需要将Coefs廓线和插值以后的大气水汽和水凝物廓线含量信息进行简单的计算（个人理解），计算出各层的optical depth等宏观辐射参数，再次插值回到`user level`，最后再进行辐射值的逐层计算；个人理解，有很大问题，比如辐射传输的变量其实包含了很多因子，不是就简单的两三个光学参数就概括的。

```fortran
    INTEGER(jpim)          :: nlevels             ! number of levels(pres/absorber) same for all gases
    INTEGER(jpim)          :: nlayers             ! number of layers(pres/absorber) nlevels-1
```
这里的Level指的是大气廓线层，对于绝大部分的大气辐射传输模式来说，所有气体的浓度、气压分层、风场、水汽湿度和对应气体的RT参数都是在Level层上；而Layer是Level的夹层，这样nlayer= nlevels-1, Layer一般是云水、云冰、雨、雪、霰等水凝物的分层。

打开`rtcoef_rttov13/rttov13pred54L/rtcoef_fy3_4_mwri.dat`你会发现，最底层存储了三种气体的FASTEM 透射系数廓线 coef(ncoef, levels)
```fortran
FAST_COEFFICIENTS
 ! 
 ! Transmission coefficients
 ! Order of the gases:
 !     Mixed_gases 
 !     Water_vapour
 !     WV_Continuum
Mixed_gases 
  0.22347463E-12  0.97183104E-12  0.11559487E-10 -0.28401033E-10  0.24825281E-12
 -0.13350612E-11  0.11559461E-10  0.11023843E-10 -0.35628344E-11  0.10916160E-09
 -0.15138966E-11 -0.23916756E-09  0.16813143E-09 -0.15779782E-12  0.16958355E-11
  0.12455733E-10 -0.45041561E-10  0.54250922E-11  0.40792076E-09 -0.31569689E-12
 -0.86874801E-09  0.65760704E-09 -0.48280489E-12  0.40331001E-12  0.28221647E-10
 -0.18337354E-09  0.60724938E-12  0.12792379E-08  0.34500769E-12 -0.26868266E-08
 ......
 END
```

#### rttov_profile 大气廓线
这里的profile指的是大气气体和水汽、温压风的大气廓线，用于晴空大气的微波辐射传输（当然了，本笔记只考虑被动微波模拟，所以理解会狭隘一些）。
结构体定义如下：
```fortran
  !> Atmospheric profiles on model pressure levels
  TYPE rttov_profile
    ......
    REAL(jprb), POINTER :: p(:)          => NULL() !< Pressure (hPa)
    REAL(jprb), POINTER :: t(:)          => NULL() !< Temperature (K)
    REAL(jprb), POINTER :: q(:)          => NULL() !< Water vapour (ppmv or kg/kg)
    REAL(jprb), POINTER :: o3(:)         => NULL() !< O3 (ppmv or kg/kg)
    REAL(jprb), POINTER :: co2(:)        => NULL() !< CO2 (ppmv or kg/kg)
    REAL(jprb), POINTER :: n2o(:)        => NULL() !< N2O (ppmv or kg/kg)
    REAL(jprb), POINTER :: co(:)         => NULL() !< CO (ppmv or kg/kg)
    REAL(jprb), POINTER :: ch4(:)        => NULL() !< CH4 (ppmv or kg/kg)
    REAL(jprb), POINTER :: so2(:)        => NULL() !< SO2 (ppmv or kg/kg)
    REAL(jprb), POINTER :: clw(:)        => NULL() !< Cloud liquid water absorption only (kg/kg)

    LOGICAL(jplm)       :: mmr_cldaer              !< Cloud/aerosol units: False => num density cm-3 (aer), g.m-3 (cld);
                                                   !!                      True => kg/kg (cld, aer)
    REAL(jprb), POINTER :: aerosols(:,:) => NULL() !< Aerosol layer concentrations (units: see mmr_cldaer)
    REAL(jprb), POINTER :: cloud(:,:)    => NULL() !< Cloud liquid/ice water layer concentrations (units: see mmr_cldaer)
    REAL(jprb), POINTER :: cfrac(:)      => NULL() !< Layer cloud fraction (0-1)
    REAL(jprb), POINTER :: clwde(:)      => NULL() !< Cloud liquid water particle effective diameter (microns)
    INTEGER(jpim)       :: clwde_param             !< Cloud liquid water effective diameter parameterisation (1)
    INTEGER(jpim)       :: clw_scheme              !< Select liquid water cloud scheme (1-2)
    REAL(jprb), POINTER :: icede(:)      => NULL() !< Ice particle effective diameter (microns)
    INTEGER(jpim)       :: icede_param             !< Ice particle effective diameter parameterisation (1-4)
    INTEGER(jpim)       :: ice_scheme              !< Select ice cloud scheme (1-2)

    TYPE(rttov_skin) :: skin                       !< Surface skin variables
    TYPE(rttov_s2m)  :: s2m                        !< Surface 2m variables
    .......
  END TYPE rttov_profile
```
可以看出，这里不仅考虑了大气物理参数，还提供过了多种气体分子的浓度廓线。包括O3, CO2, NO2, CO, CH4和 SO2。并且考虑了云水/冰和气溶胶。注意，这里的设置兼顾了IR和VIS等波段的模拟，所以考虑的气体较多，实际上很多非极性分子对微波辐射传输不起作用。
这里的云水设置在MW模拟时用于计算云的的Absorption,而散射则由MW_SCAT modelue 负责（见下文）：
```fortran rttov_direct.F90: line 675
!--------------------------------------------------------------------------
! MW CLW absorption optical depths
!--------------------------------------------------------------------------
  IF (sensor_mw .AND. opts%rt_mw%clw_data) THEN
    CALL rttov_mw_clw_absorption( &0
    ......
  ENDIF
```
profile 结构体还包含了地表参数的输入，分别由`rttov_skin`,`rttov_s2m`子结构体来定义，包含的参数有：地表类型、水体盐度、土壤湿度、积雪分数、skin temperature; 2m 温压湿风和抽样浓度。
此外，profile还包含了必要的时间和经纬度变量。
对于profile的输入，`example.fwd.F90`是通过读取给定文件的方式读取的，后面的章节我们会介绍到。


### 四个辐射参数类型
- rttov_emissivity 
``` fortran 
  !> Input/output surface emissivities, declare as array of size nchanprof,
  !! also used to specify per-channel specularity parameter for use with do_lambertian option
  TYPE rttov_emissivity
    REAL(jprb)    :: emis_in     = 0._jprb !< Input emissivity (0-1)
    REAL(jprb)    :: emis_out    = 0._jprb !< Output emissivity, value used by RTTOV (0-1)
    REAL(jprb)    :: specularity = 0._jprb !< Mix of specular/Lambertian for do_lambertian option (0-1, 0=Lambertian/1=specular)
    REAL(jprb)    :: tskin_eff   = 0._jprb !< Effective skin temperature per channel (if use_tskin_eff true)
  END TYPE
```
emis_out：地表模式或atlas输出的Emissivity数据集，作为RTTOV大气辐射传输的地表背景输入。
- rttov_reflectance微波段不使用，可以理解为1-emissivity；
- rttov_transmission 各层的透过率输出；
- rttov_radiance：天顶出射的模拟辐亮度及对应的亮温和BRF
```fortran
  !> Output radiances and corresponding brightness temperatures and reflectances (BRFs)
  !! Radiance unit: mW/m2/sr/cm-1; BT unit: K; BRFs are unitless.
  !! Array sizes are (nchannels) or (nlayers,nchannels)
  TYPE rttov_radiance
    REAL(jprb),    POINTER :: clear(:)              => NULL() !< Clear sky radiance
    REAL(jprb),    POINTER :: total(:)              => NULL() !< Cloudy radiance for given cloud
    REAL(jprb),    POINTER :: bt_clear(:)           => NULL() !< Brightness temp equivalent to clear radiance
    REAL(jprb),    POINTER :: bt(:)                 => NULL() !< Brightness temp equivalent to total radiance
    REAL(jprb),    POINTER :: refl_clear(:)         => NULL() !< Reflectance calculated from clear radiance
    REAL(jprb),    POINTER :: refl(:)               => NULL() !< Reflectance calculated from total radiance
    REAL(jprb),    POINTER :: overcast(:,:)         => NULL() !< Overcast radiance for opaque cloud at level bounding
                                                              !!   bottom of each layer
    REAL(jprb),    POINTER :: cloudy(:)             => NULL() !< 100% cloudy radiance for given cloud (simple cloud scheme)
                                                              !!   or same as total (addclouds/addaerosl true)
    ......
  END TYPE rttov_radiance
```

### MW-SCAT 变量
除了上述的变量类型，RTTOV还包含用于MW-SCAT模块散射计算特有的变量。
#### rttov_scatt_coef hydrotable系数表
定义了用来缓存水凝物辐射传输系数的变量：
```fortran
  !> @internal RTTOV-SCATT coefs
  TYPE rttov_scatt_coef
    ! Structure for the storage of RTTOV_SCATT coefficients (hydrotable files)
    INTEGER(jpim) :: rttov_version ! RTTOV compatibility version
    INTEGER(jpim) :: nhydro        ! Number of hydrometeors in computation
    INTEGER(jpim) :: mtype         ! Number of hydrometeors
    INTEGER(jpim) :: mfreqm        ! Number of frequencies
    INTEGER(jpim) :: mtemp         ! Number of temperature bins
    INTEGER(jpim) :: mwc           ! Number of water bins
    ......
    REAL(jprb), POINTER :: freq(:)      => NULL()   ! list of frequencies in hydrotable
    INTEGER(jpim), POINTER :: mpol(:)   => NULL()   ! list of polarisations in Hydro table (-1 = polarisation ignored)
    REAL(jprb), POINTER :: ext(:,:,:,:) => NULL()   ! extinction coefficent table
    REAL(jprb), POINTER :: ssa(:,:,:,:) => NULL()   ! single scattering albedo table
    REAL(jprb), POINTER :: asp(:,:,:,:) => NULL()   ! asymmetry parameter table
    ......
  END TYPE rttov_scatt_coef
```
这里水凝物Hydrotable参数表考虑的FAST RT 参数包括`消光系数`，`单散反照率`和`不对称因子`。熟悉辐射传输的同学应该知道，这是在一定的形状（如球形）假设和散射理论（Mie）计算下, 水凝物粒子的散射计算的结果。这些参数辅以一定的粒径谱假设（如Gamma分布）就可以进一步计算出光学厚度，散射系数等介质的单层光学系数。不对称因子则描述了散射项函数的不对称性，在二流近似的情形下决定了前向和后向散射强度的比值。

#### MW-SCATT Cloud/hydrometeor profile
MW-SCAT模块的水凝物廓线变量：
```fortran
  !> Additional atmospheric cloud/hydrometeor profile input for RTTOV-SCATT
  TYPE rttov_profile_cloud
    INTEGER(jpim) :: nlevels                         !< Number of atmospheric levels (same as in rttov_profile)
    INTEGER(jpim) :: nhydro                          !< Number of hydrometeor types
    INTEGER(jpim) :: nhydro_frac                     !< Number of hydrometeor fractions (should be 1 or nhydro)
    INTEGER(jpim), POINTER :: flux_conversion(:) => NULL() !< 0 => no flux conv, input units are kg/kg,
                                                           !< 1,2 => input kg/m2/s rain,snow
    REAL(jprb)    :: cfrac                           !< Average cloud fraction (only used if lusercfrac = TRUE, 0-1)
    REAL(jprb), POINTER :: ph(:)           => NULL() !< nlevels+1 of half-level model pressures (hPa)
    REAL(jprb), POINTER :: hydro(:,:)      => NULL() !< nlevels by ntypes of hydrometeor (kg/kg or (deprecated) kg/m2/s)
    REAL(jprb), POINTER :: hydro_frac(:,:) => NULL() !< nlevels of hydrometeor fraction (cloud cover) (0-1)
  END TYPE rttov_profile_cloud
```

### 反演相关变量类型
(Baordo and Geer, 2016)反演方案实例中使用的中间变量类型，其思想是通过才cfrac将dirct和和SCAT的模拟进行加权，以实现对云天的散射和吸收效应的考量。
```fortran
  !> Output variables from RTTOV-SCATT enabling emissivity retrievals
  !! All arrays are size(chanprof)
  TYPE rttov_scatt_emis_retrieval_type
    REAL(jprb), POINTER  :: cfrac(:)    => NULL() !< RTTOV-SCATT effective cloud fraction
                                                  !!  (Tallsky = cfrac * Tcld + (1-cfrac) * Tclr
    REAL(jprb), POINTER  :: bsfc(:)     => NULL() !< Surface black-body radiance (i.e. Planck(Tsfc))
    REAL(jprb), POINTER  :: tau_cld(:)  => NULL() !< Along-path transmittance, surface to space (cloudy column)
    REAL(jprb), POINTER  :: up_cld(:)   => NULL() !< TOA upwelling radiance from atmosphere (not inc. surface
                                                  !!  emission or reflection)
    REAL(jprb), POINTER  :: down_cld(:) => NULL() !< SFC downwelling radiance (inc. cosmic term)
    REAL(jprb), POINTER  :: tau_clr(:)  => NULL() !< Along-path transmittance, surface to space (clear column)
    REAL(jprb), POINTER  :: up_clr(:)   => NULL() !< TOA upwelling radiance from atmosphere (not inc. surface
                                                  !!  emission or reflection)
    REAL(jprb), POINTER  :: down_clr(:) => NULL() !< SFC downwelling radiance (inc. cosmic term)
  END TYPE rttov_scatt_emis_retrieval_type
```

## 变量分配
从`src/test/exmaple_rttovscatt_fwd.F90`和`src/test/example_fwd.F90`两个正向模拟示例来看，上述变量的分配（allocate）任务交给了`rttov_alloc_direct.F90`的同名函数来完成的: `Allocate/deallocate any/all arrays and structures for RTTOV or RTTOV-SCATT direct models`。并且对于大部分变量结构的输入是可选的（optional）。

## 小结
本节对RTTOV direct和mw-scat两个模块涉及到的主要变量类型和结构体进行了介绍，穿插了自己的理解。
后面的章节会分别对于上述的两种module进行介绍，主要着眼与应用，不会有太多源码的解读，
再次声明，本系列笔记是个人理解，希望在记录学习过程的同时能抛砖引玉，共同进步。

---
参考：1. 源码；2. user guide.