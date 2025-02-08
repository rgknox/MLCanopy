module MLCanopyFluxesMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate  multilayer canopy fluxes
  !
  ! !USES:

  use MLCanopyCouplerMod, only : endrun
  use MLCanopyCouplerMod, only : ispval
  use MLCanopyCouplerMod, only : spval
  use MLCanopyCouplerMod, only : nlevgrnd
  use MLCanopyCouplerMod, only : numrad
  use MLCanopyCouplerMod, only : r8
  use MLCanopyCouplerMod, only : iulog

  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PRIVATE TYPES:
  integer, parameter :: nvar1d = 12     ! Number of single-level fluxes to accumulate over sub-time steps
  integer, parameter :: nvar2d = 4      ! Number of multi-level profile fluxes to accumulate over sub-time steps
  integer, parameter :: nvar3d = 10     ! Number of multi-level leaf fluxes to accumulate over sub-time steps
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: MLCanopyFluxes              ! Compute canopy fluxes
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: SubTimeStepFluxIntegration ! Integrate fluxes over model sub-time steps
  private :: CanopyFluxesDiagnostics    ! Sum leaf and soil fluxes and other canopy diagnostics
  !-----------------------------------------------------------------------

  contains



  !-----------------------------------------------------------------------
  subroutine SubTimeStepFluxIntegration (niter, num_sub_steps, num_filter, filter, &
  flux_accumulator, flux_accumulator_profile, flux_accumulator_leaf, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Integrate fluxes over model sub-time steps
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: niter                                ! Current sub-time step
    integer, intent(in) :: num_sub_steps                        ! Number of sub-time steps
    integer, intent(in) :: num_filter                           ! Number of patches in filter
    integer, intent(in) :: filter(:)                            ! Patch filter
    real(r8), intent(inout) :: flux_accumulator(:,:)            ! Single-level flux accumulator variable
    real(r8), intent(inout) :: flux_accumulator_profile(:,:,:)  ! Multi-level profile flux accumulator variable
    real(r8), intent(inout) :: flux_accumulator_leaf(:,:,:,:)   ! Multi-level leaf flux accumulator variable
    type(mlcanopy_type), intent(in) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                                              ! Filter index
    integer  :: p                                               ! Patch index for CLM g/l/c/p hierarchy
    integer  :: i,j,k                                           ! Variable index
    !---------------------------------------------------------------------

    associate ( &
    ustar       => mlcanopy_inst%ustar_canopy         , &  ! Friction velocity (m/s)
    lwup        => mlcanopy_inst%lwup_canopy          , &  ! Upward longwave radiation above canopy (W/m2)
    qflx_intr   => mlcanopy_inst%qflx_intr_canopy     , &  ! Intercepted precipitation (kg H2O/m2/s)
    qflx_tflrain => mlcanopy_inst%qflx_tflrain_canopy , &  ! Total rain throughfall onto ground (kg H2O/m2/s)
    qflx_tflsnow => mlcanopy_inst%qflx_tflsnow_canopy , &  ! Total snow throughfall onto ground (kg H2O/m2/s)
    lwsoi       => mlcanopy_inst%lwsoi_soil           , &  ! Absorbed longwave radiation: ground (W/m2)
    rnsoi       => mlcanopy_inst%rnsoi_soil           , &  ! Net radiation: ground (W/m2)
    shsoi       => mlcanopy_inst%shsoi_soil           , &  ! Sensible heat flux: ground (W/m2)
    lhsoi       => mlcanopy_inst%lhsoi_soil           , &  ! Latent heat flux: ground (W/m2)
    etsoi       => mlcanopy_inst%etsoi_soil           , &  ! Water vapor flux: ground (mol H2O/m2/s)
    gsoi        => mlcanopy_inst%gsoi_soil            , &  ! Soil heat flux (W/m2)
    gac0        => mlcanopy_inst%gac0_soil            , &  ! Aerodynamic conductance for soil fluxes (mol/m2/s)
    shair       => mlcanopy_inst%shair_profile        , &  ! Canopy layer air sensible heat flux (W/m2)
    etair       => mlcanopy_inst%etair_profile        , &  ! Canopy layer air water vapor flux (mol H2O/m2/s)
    stair       => mlcanopy_inst%stair_profile        , &  ! Canopy layer air storage heat flux (W/m2)
    gac         => mlcanopy_inst%gac_profile          , &  ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
    lwleaf      => mlcanopy_inst%lwleaf_leaf          , &  ! Leaf absorbed longwave radiation (W/m2 leaf)
    rnleaf      => mlcanopy_inst%rnleaf_leaf          , &  ! Leaf net radiation (W/m2 leaf)
    shleaf      => mlcanopy_inst%shleaf_leaf          , &  ! Leaf sensible heat flux (W/m2 leaf)
    lhleaf      => mlcanopy_inst%lhleaf_leaf          , &  ! Leaf latent heat flux (W/m2 leaf)
    trleaf      => mlcanopy_inst%trleaf_leaf          , &  ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    evleaf      => mlcanopy_inst%evleaf_leaf          , &  ! Leaf evaporation flux (mol H2O/m2 leaf/s)
    stleaf      => mlcanopy_inst%stleaf_leaf          , &  ! Leaf storage heat flux (W/m2 leaf)
    anet        => mlcanopy_inst%anet_leaf            , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    agross      => mlcanopy_inst%agross_leaf          , &  ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    gs          => mlcanopy_inst%gs_leaf                &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    )

    do fp = 1, num_filter
       p = filter(fp)

       ! Initialize flux variables that are summed over sub-time steps

       if (niter == 1) then
          flux_accumulator(p,:) = 0._r8
          flux_accumulator_profile(p,:,:) = 0._r8
          flux_accumulator_leaf(p,:,:,:) = 0._r8
       end if

       ! Accumulate fluxes over sub-time steps

       i = 0
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + ustar(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + lwup(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + lwsoi(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + rnsoi(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + shsoi(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + lhsoi(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + etsoi(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + gsoi(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + gac0(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + qflx_intr(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + qflx_tflrain(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + qflx_tflsnow(p)

       j = 0
       j = j + 1; flux_accumulator_profile(p,:,j) = flux_accumulator_profile(p,:,j) + shair(p,:)
       j = j + 1; flux_accumulator_profile(p,:,j) = flux_accumulator_profile(p,:,j) + etair(p,:)
       j = j + 1; flux_accumulator_profile(p,:,j) = flux_accumulator_profile(p,:,j) + stair(p,:)
       j = j + 1; flux_accumulator_profile(p,:,j) = flux_accumulator_profile(p,:,j) + gac(p,:)

       k = 0
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + lwleaf(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + rnleaf(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + shleaf(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + lhleaf(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + trleaf(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + evleaf(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + stleaf(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + anet(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + agross(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + gs(p,:,:)

       if (i > nvar1d .or. j > nvar2d .or. k > nvar3d) then
          call endrun (msg=' ERROR: SubTimeStepFluxIntegration: nvar error')
       end if

       ! Average fluxes over sub-timesteps

       if (niter == num_sub_steps) then

          ! Time averaging

          flux_accumulator(p,:) = flux_accumulator(p,:) / float(num_sub_steps)
          flux_accumulator_profile(p,:,:) = flux_accumulator_profile(p,:,:) / float(num_sub_steps)
          flux_accumulator_leaf(p,:,:,:) = flux_accumulator_leaf(p,:,:,:) / float(num_sub_steps)

          ! Map fluxes to variables: variables must be in the same order as above

          i = 0
          i = i + 1; ustar(p) = flux_accumulator(p,i)
          i = i + 1; lwup(p) = flux_accumulator(p,i)
          i = i + 1; lwsoi(p) = flux_accumulator(p,i)
          i = i + 1; rnsoi(p) = flux_accumulator(p,i)
          i = i + 1; shsoi(p) = flux_accumulator(p,i)
          i = i + 1; lhsoi(p) = flux_accumulator(p,i)
          i = i + 1; etsoi(p) = flux_accumulator(p,i)
          i = i + 1; gsoi(p) = flux_accumulator(p,i)
          i = i + 1; gac0(p) = flux_accumulator(p,i)
          i = i + 1; qflx_intr(p) = flux_accumulator(p,i)
          i = i + 1; qflx_tflrain(p) = flux_accumulator(p,i)
          i = i + 1; qflx_tflsnow(p) = flux_accumulator(p,i)

          j = 0
          j = j + 1; shair(p,:) = flux_accumulator_profile(p,:,j)
          j = j + 1; etair(p,:) = flux_accumulator_profile(p,:,j)
          j = j + 1; stair(p,:) = flux_accumulator_profile(p,:,j)
          j = j + 1; gac(p,:) = flux_accumulator_profile(p,:,j)

          k = 0
          k = k + 1; lwleaf(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; rnleaf(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; shleaf(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; lhleaf(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; trleaf(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; evleaf(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; stleaf(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; anet(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; agross(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; gs(p,:,:) = flux_accumulator_leaf(p,:,:,k)

          if (i > nvar1d .or. j > nvar2d .or. k > nvar3d) then
             call endrun (msg=' ERROR: SubTimeStepFluxIntegration: nvar error')
          end if

       end if

    end do

    end associate
  end subroutine SubTimeStepFluxIntegration

  !-----------------------------------------------------------------------
  subroutine CanopyFluxesDiagnostics (num_filter, filter, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Sum leaf and soil fluxes to get canopy fluxes and calculate
    ! other canopy diagnostics
    !
    ! !USES:
    use clm_varpar, only : ivis, inir
    use MLclm_varctl, only : turb_type
    use MLclm_varpar, only : isun, isha
    use MLWaterVaporMod, only : LatVap
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_filter   ! Number of patches in filter
    integer, intent(in) :: filter(:)    ! Patch filter
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                      ! Filter index
    integer  :: p                       ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic                      ! Aboveground layer index
    real(r8) :: err                     ! Energy imbalance (W/m2)
    real(r8) :: radin                   ! Incoming radiation (W/m2)
    real(r8) :: radout                  ! Outgoing radiation (W/m2)
    real(r8) :: avail                   ! Available energy (W/m2)
    real(r8) :: flux                    ! Turbulent fluxes + storage (W/m2)
    real(r8) :: fracgreen               ! Green fraction of plant area index: lai/(lai+sai)
    real(r8) :: minlwp                  ! Minimum leaf water potential for canopy water stress diagnostic (MPa)
    !---------------------------------------------------------------------

    associate ( &
                                                         ! *** Input ***
    tref        => mlcanopy_inst%tref_forcing       , &  ! Air temperature at reference height (K)
    swskyb      => mlcanopy_inst%swskyb_forcing     , &  ! Atmospheric direct beam solar radiation (W/m2)
    swskyd      => mlcanopy_inst%swskyd_forcing     , &  ! Atmospheric diffuse solar radiation (W/m2)
    lwsky       => mlcanopy_inst%lwsky_forcing      , &  ! Atmospheric longwave radiation (W/m2)
    ncan        => mlcanopy_inst%ncan_canopy        , &  ! Number of aboveground layers
    ntop        => mlcanopy_inst%ntop_canopy        , &  ! Index for top leaf layer
    lai         => mlcanopy_inst%lai_canopy         , &  ! Leaf area index of canopy (m2/m2)
    sai         => mlcanopy_inst%sai_canopy         , &  ! Stem area index of canopy (m2/m2)
    swveg       => mlcanopy_inst%swveg_canopy       , &  ! Absorbed solar radiation: vegetation (W/m2)
    albcan      => mlcanopy_inst%albcan_canopy      , &  ! Albedo above canopy (-)
    lwup        => mlcanopy_inst%lwup_canopy        , &  ! Upward longwave radiation above canopy (W/m2)
    shsoi       => mlcanopy_inst%shsoi_soil         , &  ! Sensible heat flux: ground (W/m2)
    lhsoi       => mlcanopy_inst%lhsoi_soil         , &  ! Latent heat flux: ground (W/m2)
    gsoi        => mlcanopy_inst%gsoi_soil          , &  ! Soil heat flux (W/m2)
    swsoi       => mlcanopy_inst%swsoi_soil         , &  ! Absorbed solar radiation: ground (W/m2)
    lwsoi       => mlcanopy_inst%lwsoi_soil         , &  ! Absorbed longwave radiation: ground (W/m2)
    etsoi       => mlcanopy_inst%etsoi_soil         , &  ! Water vapor flux: ground (mol H2O/m2/s)
    dpai        => mlcanopy_inst%dpai_profile       , &  ! Canopy layer plant area index (m2/m2)
    fwet        => mlcanopy_inst%fwet_profile       , &  ! Canopy layer fraction of plant area index that is wet
    fdry        => mlcanopy_inst%fdry_profile       , &  ! Canopy layer fraction of plant area index that is green and dry
    tair        => mlcanopy_inst%tair_profile       , &  ! Canopy layer air temperature (K)
    wind        => mlcanopy_inst%wind_profile       , &  ! Canopy layer wind speed (m/s)
    shair       => mlcanopy_inst%shair_profile      , &  ! Canopy layer air sensible heat flux (W/m2)
    etair       => mlcanopy_inst%etair_profile      , &  ! Canopy layer air water vapor flux (mol H2O/m2/s)
    stair       => mlcanopy_inst%stair_profile      , &  ! Canopy layer air storage heat flux (W/m2)
    fracsun     => mlcanopy_inst%fracsun_profile    , &  ! Canopy layer sunlit fraction (-)
    vcmax25_profile => mlcanopy_inst%vcmax25_profile, &  ! Canopy layer leaf maximum carboxylation rate at 25C (umol/m2/s)
    lwleaf      => mlcanopy_inst%lwleaf_leaf        , &  ! Leaf absorbed longwave radiation (W/m2 leaf)
    rnleaf      => mlcanopy_inst%rnleaf_leaf        , &  ! Leaf net radiation (W/m2 leaf)
    stleaf      => mlcanopy_inst%stleaf_leaf        , &  ! Leaf storage heat flux (W/m2 leaf)
    shleaf      => mlcanopy_inst%shleaf_leaf        , &  ! Leaf sensible heat flux (W/m2 leaf)
    lhleaf      => mlcanopy_inst%lhleaf_leaf        , &  ! Leaf latent heat flux (W/m2 leaf)
    trleaf      => mlcanopy_inst%trleaf_leaf        , &  ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    evleaf      => mlcanopy_inst%evleaf_leaf        , &  ! Leaf evaporation flux (mol H2O/m2 leaf/s)
    swleaf      => mlcanopy_inst%swleaf_leaf        , &  ! Leaf absorbed solar radiation (W/m2 leaf)
    agross      => mlcanopy_inst%agross_leaf        , &  ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    apar        => mlcanopy_inst%apar_leaf          , &  ! Leaf absorbed PAR (umol photon/m2 leaf/s)
    anet        => mlcanopy_inst%anet_leaf          , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    gs          => mlcanopy_inst%gs_leaf            , &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    tleaf       => mlcanopy_inst%tleaf_leaf         , &  ! Leaf temperature (K)
    lwp         => mlcanopy_inst%lwp_leaf           , &  ! Leaf water potential (MPa)
    vcmax25_leaf=> mlcanopy_inst%vcmax25_leaf       , &  ! Leaf maximum carboxylation rate at 25C (umol/m2/s)
                                                         ! *** Output ***
    rnet        => mlcanopy_inst%rnet_canopy        , &  ! Total net radiation, including soil (W/m2)
    stflx       => mlcanopy_inst%stflx_canopy       , &  ! Canopy storage heat flux (W/m2)
    shflx       => mlcanopy_inst%shflx_canopy       , &  ! Total sensible heat flux, including soil (W/m2)
    lhflx       => mlcanopy_inst%lhflx_canopy       , &  ! Total latent heat flux, including soil (W/m2)
    etflx       => mlcanopy_inst%etflx_canopy       , &  ! Total water vapor flux, including soil (mol H2O/m2/s)
    lwveg       => mlcanopy_inst%lwveg_canopy       , &  ! Absorbed longwave radiation: vegetation (W/m2)
    lwvegsun    => mlcanopy_inst%lwvegsun_canopy    , &  ! Absorbed longwave radiation: sunlit canopy (W/m2)
    lwvegsha    => mlcanopy_inst%lwvegsha_canopy    , &  ! Absorbed longwave radiation: shaded canopy (W/m2)
    shveg       => mlcanopy_inst%shveg_canopy       , &  ! Sensible heat flux: vegetation (W/m2)
    shvegsun    => mlcanopy_inst%shvegsun_canopy    , &  ! Sensible heat flux: sunlit canopy (W/m2)
    shvegsha    => mlcanopy_inst%shvegsha_canopy    , &  ! Sensible heat flux: shaded canopy (W/m2)
    lhveg       => mlcanopy_inst%lhveg_canopy       , &  ! Latent heat flux: vegetation (W/m2)
    lhvegsun    => mlcanopy_inst%lhvegsun_canopy    , &  ! Latent heat flux: sunlit canopy (W/m2)
    lhvegsha    => mlcanopy_inst%lhvegsha_canopy    , &  ! Latent heat flux: shaded canopy (W/m2)
    etveg       => mlcanopy_inst%etveg_canopy       , &  ! Water vapor flux: vegetation (mol H2O/m2/s)
    etvegsun    => mlcanopy_inst%etvegsun_canopy    , &  ! Water vapor flux: sunlit canopy (mol H2O/m2/s)
    etvegsha    => mlcanopy_inst%etvegsha_canopy    , &  ! Water vapor flux: shaded canopy (mol H2O/m2/s)
    gppveg      => mlcanopy_inst%gppveg_canopy      , &  ! Gross primary production: vegetation (umol CO2/m2/s)
    gppvegsun   => mlcanopy_inst%gppvegsun_canopy   , &  ! Gross primary production: sunlit canopy (umol CO2/m2/s)
    gppvegsha   => mlcanopy_inst%gppvegsha_canopy   , &  ! Gross primary production: shaded canopy (umol CO2/m2/s)
    vcmax25veg  => mlcanopy_inst%vcmax25veg_canopy  , &  ! Vcmax at 25C: total canopy (umol/m2/s)
    vcmax25sun  => mlcanopy_inst%vcmax25sun_canopy  , &  ! Vcmax at 25C: sunlit canopy (umol/m2/s)
    vcmax25sha  => mlcanopy_inst%vcmax25sha_canopy  , &  ! Vcmax at 25C: shaded canopy (umol/m2/s)
    gsveg       => mlcanopy_inst%gsveg_canopy       , &  ! Stomatal conductance: canopy (mol H2O/m2/s)
    gsvegsun    => mlcanopy_inst%gsvegsun_canopy    , &  ! Stomatal conductance: sunlit canopy (mol H2O/m2/s)
    gsvegsha    => mlcanopy_inst%gsvegsha_canopy    , &  ! Stomatal conductance: shaded canopy (mol H2O/m2/s)
    windveg     => mlcanopy_inst%windveg_canopy     , &  ! Wind speed: canopy (m/s)
    windvegsun  => mlcanopy_inst%windvegsun_canopy  , &  ! Wind speed: sunlit canopy (m/s)
    windvegsha  => mlcanopy_inst%windvegsha_canopy  , &  ! Wind speed: shaded canopy (m/s)
    tlveg       => mlcanopy_inst%tlveg_canopy       , &  ! Leaf temperature: canopy (K)
    tlvegsun    => mlcanopy_inst%tlvegsun_canopy    , &  ! Leaf temperature: sunlit canopy (K)
    tlvegsha    => mlcanopy_inst%tlvegsha_canopy    , &  ! Leaf temperature: shaded canopy (K)
    taveg       => mlcanopy_inst%taveg_canopy       , &  ! Air temperature: canopy (K)
    tavegsun    => mlcanopy_inst%tavegsun_canopy    , &  ! Air temperature: sunlit canopy (K)
    tavegsha    => mlcanopy_inst%tavegsha_canopy    , &  ! Air temperature: shaded canopy (K)
    laisun      => mlcanopy_inst%laisun_canopy      , &  ! Canopy plant area index (lai+sai): sunlit canopy (m2/m2) 
    laisha      => mlcanopy_inst%laisha_canopy      , &  ! Canopy plant area index (lai+sai): shaded canopy (m2/m2) 
    fracminlwp  => mlcanopy_inst%fracminlwp_canopy  , &  ! Fraction of canopy that is water-stressed
    swsrc       => mlcanopy_inst%swsrc_profile      , &  ! Canopy layer source/sink flux: absorbed solar radiation (W/m2)
    lwsrc       => mlcanopy_inst%lwsrc_profile      , &  ! Canopy layer source/sink flux: absorbed longwave radiation (W/m2)
    rnsrc       => mlcanopy_inst%rnsrc_profile      , &  ! Canopy layer source/sink flux: net radiation (W/m2)
    stsrc       => mlcanopy_inst%stsrc_profile      , &  ! Canopy layer source/sink flux: storage heat flux (W/m2)
    shsrc       => mlcanopy_inst%shsrc_profile      , &  ! Canopy layer source/sink flux: sensible heat (W/m2)
    lhsrc       => mlcanopy_inst%lhsrc_profile      , &  ! Canopy layer source/sink flux: latent heat (W/m2)
    etsrc       => mlcanopy_inst%etsrc_profile      , &  ! Canopy layer source/sink flux: water vapor (mol H2O/m2/s)
    fco2src     => mlcanopy_inst%fco2src_profile    , &  ! Canopy layer source/sink flux: CO2 (umol CO2/m2/s)
    swleaf_mean => mlcanopy_inst%swleaf_mean_profile, &  ! Canopy layer weighted mean: leaf absorbed solar radiation (W/m2 leaf)
    lwleaf_mean => mlcanopy_inst%lwleaf_mean_profile, &  ! Canopy layer weighted mean: leaf absorbed longwave radiation (W/m2 leaf)
    rnleaf_mean => mlcanopy_inst%rnleaf_mean_profile, &  ! Canopy layer weighted mean: leaf net radiation (W/m2 leaf)
    stleaf_mean => mlcanopy_inst%stleaf_mean_profile, &  ! Canopy layer weighted mean: leaf storage heat flux (W/m2 leaf)
    shleaf_mean => mlcanopy_inst%shleaf_mean_profile, &  ! Canopy layer weighted mean: leaf sensible heat flux (W/m2 leaf)
    lhleaf_mean => mlcanopy_inst%lhleaf_mean_profile, &  ! Canopy layer weighted mean: leaf latent heat flux (W/m2 leaf)
    etleaf_mean => mlcanopy_inst%etleaf_mean_profile, &  ! Canopy layer weighted mean: leaf water vapor flux (mol H2O/m2 leaf/s)
    fco2_mean   => mlcanopy_inst%fco2_mean_profile  , &  ! Canopy layer weighted mean: leaf net photosynthesis (umol CO2/m2 leaf/s)
    apar_mean   => mlcanopy_inst%apar_mean_profile  , &  ! Canopy layer weighted mean: leaf absorbed PAR (umol photon/m2 leaf/s)
    gs_mean     => mlcanopy_inst%gs_mean_profile    , &  ! Canopy layer weighted mean: leaf stomatal conductance (mol H2O/m2 leaf/s)
    tleaf_mean  => mlcanopy_inst%tleaf_mean_profile , &  ! Canopy layer weighted mean: leaf temperature (K)
    lwp_mean    => mlcanopy_inst%lwp_mean_profile     &  ! Canopy layer weighted mean: leaf water potential (MPa)
    )

    do fp = 1, num_filter
       p = filter(fp)

       ! Leaf flux profiles

       do ic = 1, ncan(p)
          if (dpai(p,ic) > 0._r8) then

             ! Leaf fluxes/states (per unit leaf area)

             lwleaf_mean(p,ic) = lwleaf(p,ic,isun)*fracsun(p,ic) + lwleaf(p,ic,isha)*(1._r8 - fracsun(p,ic))
             swleaf_mean(p,ic,ivis) = swleaf(p,ic,isun,ivis)*fracsun(p,ic) + swleaf(p,ic,isha,ivis)*(1._r8 - fracsun(p,ic))
             swleaf_mean(p,ic,inir) = swleaf(p,ic,isun,inir)*fracsun(p,ic) + swleaf(p,ic,isha,inir)*(1._r8 - fracsun(p,ic))
             rnleaf_mean(p,ic) = rnleaf(p,ic,isun)*fracsun(p,ic) + rnleaf(p,ic,isha)*(1._r8 - fracsun(p,ic))
             stleaf_mean(p,ic) = stleaf(p,ic,isun)*fracsun(p,ic) + stleaf(p,ic,isha)*(1._r8 - fracsun(p,ic))
             shleaf_mean(p,ic) = shleaf(p,ic,isun)*fracsun(p,ic) + shleaf(p,ic,isha)*(1._r8 - fracsun(p,ic))
             lhleaf_mean(p,ic) = lhleaf(p,ic,isun)*fracsun(p,ic) + lhleaf(p,ic,isha)*(1._r8 - fracsun(p,ic))
             etleaf_mean(p,ic) = (evleaf(p,ic,isun) + trleaf(p,ic,isun)) * fracsun(p,ic) &
                               + (evleaf(p,ic,isha) + trleaf(p,ic,isha)) * (1._r8 - fracsun(p,ic))
             fco2_mean(p,ic) = anet(p,ic,isun)*fracsun(p,ic) + anet(p,ic,isha)*(1._r8 - fracsun(p,ic))

             apar_mean(p,ic) = apar(p,ic,isun)*fracsun(p,ic) + apar(p,ic,isha)*(1._r8 - fracsun(p,ic))
             gs_mean(p,ic) = gs(p,ic,isun)*fracsun(p,ic) + gs(p,ic,isha)*(1._r8 - fracsun(p,ic))
             tleaf_mean(p,ic) = tleaf(p,ic,isun)*fracsun(p,ic) + tleaf(p,ic,isha)*(1._r8 - fracsun(p,ic))
             lwp_mean(p,ic) = lwp(p,ic,isun)*fracsun(p,ic) + lwp(p,ic,isha)*(1._r8 - fracsun(p,ic))

             ! Source fluxes (per unit ground area)

             lwsrc(p,ic) = lwleaf_mean(p,ic) * dpai(p,ic)
             swsrc(p,ic,ivis) = swleaf_mean(p,ic,ivis) * dpai(p,ic)
             swsrc(p,ic,inir) = swleaf_mean(p,ic,inir) * dpai(p,ic)
             rnsrc(p,ic) = rnleaf_mean(p,ic) * dpai(p,ic)
             stsrc(p,ic) = stleaf_mean(p,ic) * dpai(p,ic)
             shsrc(p,ic) = shleaf_mean(p,ic) * dpai(p,ic)
             lhsrc(p,ic) = lhleaf_mean(p,ic) * dpai(p,ic)
             etsrc(p,ic) = etleaf_mean(p,ic) * dpai(p,ic)
             fracgreen = fdry(p,ic) / (1._r8 - fwet(p,ic))
             fco2src(p,ic) = (anet(p,ic,isun)*fracsun(p,ic) + anet(p,ic,isha)*(1._r8 - fracsun(p,ic))) * dpai(p,ic) * fracgreen

          else

             lwleaf_mean(p,ic) = 0._r8
             swleaf_mean(p,ic,ivis) = 0._r8
             swleaf_mean(p,ic,inir) = 0._r8
             rnleaf_mean(p,ic) = 0._r8
             stleaf_mean(p,ic) = 0._r8
             shleaf_mean(p,ic) = 0._r8
             lhleaf_mean(p,ic) = 0._r8
             etleaf_mean(p,ic) = 0._r8
             fco2_mean(p,ic) = 0._r8

             apar_mean(p,ic) = 0._r8
             gs_mean(p,ic) = 0._r8
             tleaf_mean(p,ic) = 0._r8
             lwp_mean(p,ic) = 0._r8

             lwsrc(p,ic) = 0._r8
             swsrc(p,ic,ivis) = 0._r8
             swsrc(p,ic,inir) = 0._r8
             rnsrc(p,ic) = 0._r8
             stsrc(p,ic) = 0._r8
             shsrc(p,ic) = 0._r8
             lhsrc(p,ic) = 0._r8
             etsrc(p,ic) = 0._r8
             fco2src(p,ic) = 0._r8
          end if
       end do

       ! Sum leaf fluxes

       lwveg(p) = 0._r8
       stflx(p) = 0._r8
       shveg(p) = 0._r8
       lhveg(p) = 0._r8
       etveg(p) = 0._r8
       gppveg(p) = 0._r8
       vcmax25veg(p) = 0._r8
       gsveg(p) = 0._r8

       do ic = 1, ncan(p)
          lwveg(p) = lwveg(p) + lwsrc(p,ic)
          stflx(p) = stflx(p) + stsrc(p,ic)
          shveg(p) = shveg(p) + shsrc(p,ic)
          lhveg(p) = lhveg(p) + lhsrc(p,ic)
          etveg(p) = etveg(p) + etsrc(p,ic)
          if (dpai(p,ic) > 0._r8) then
             fracgreen = fdry(p,ic) / (1._r8 - fwet(p,ic))
             gppveg(p) = gppveg(p) &
                       + (agross(p,ic,isun)*fracsun(p,ic) + agross(p,ic,isha)*(1._r8 - fracsun(p,ic))) * dpai(p,ic) * fracgreen
             gsveg(p) = gsveg(p) + gs_mean(p,ic) * dpai(p,ic)
          end if
          vcmax25veg(p) = vcmax25veg(p) + vcmax25_profile(p,ic) * dpai(p,ic)
       end do

       ! Check energy balance for conservation

       err = swveg(p,ivis) + swveg(p,inir) + lwveg(p) - shveg(p) - lhveg(p) - stflx(p)
       if (abs(err) >= 1.e-03_r8) then
          call endrun (msg=' ERROR: CanopyFluxesDiagnostics: energy conservation error (1)')
       end if

       ! Turbulent fluxes

       select case (turb_type)
       case (0, -1)
          ! Sum of vegetation and soil fluxes
          shflx(p) = shveg(p) + shsoi(p)
          etflx(p) = etveg(p) + etsoi(p)
          lhflx(p) = lhveg(p) + lhsoi(p)
       case (1)
          ! Turbulent fluxes are at the top of the canopy
          shflx(p) = shair(p,ntop(p))
          etflx(p) = etair(p,ntop(p))
          lhflx(p) = etair(p,ntop(p)) * LatVap(tref(p))
       case default
          call endrun (msg=' ERROR: CanopyFluxesDiagnostics: turbulence type not valid')
       end select

       ! Add canopy air heat storage to storage flux

       do ic = 1, ntop(p)
          stflx(p) = stflx(p) + stair(p,ic)
       end do

       ! Overall energy balance check:
       ! radiation in - radiation out - soil heat = available energy = turbulent flux + canopy storage flux

       rnet(p) = swveg(p,ivis) + swveg(p,inir) + swsoi(p,ivis) + swsoi(p,inir) + lwveg(p) + lwsoi(p)
       radin = swskyb(p,ivis) + swskyd(p,ivis) + swskyb(p,inir) + swskyd(p,inir) + lwsky(p)
       radout = albcan(p,ivis)*(swskyb(p,ivis)+swskyd(p,ivis)) + albcan(p,inir)*(swskyb(p,inir)+swskyd(p,inir)) + lwup(p)

       err = rnet(p) - (radin - radout)
       if (abs(err) > 0.01_r8) then
          call endrun (msg=' ERROR: CanopyFluxesDiagnostics: energy conservation error (2)')
       end if

       avail = radin - radout - gsoi(p)
       flux = shflx(p) + lhflx(p) + stflx(p)
       err = avail - flux
       if (abs(err) > 0.01_r8) then
          call endrun (msg=' ERROR: CanopyFluxesDiagnostics: energy conservation error (3)')
       end if

       ! Sunlit and shaded canopy fluxes

       laisun(p) = 0._r8 ; laisha(p) = 0._r8
       lwvegsun(p) = 0._r8 ; lwvegsha(p) = 0._r8
       shvegsun(p) = 0._r8 ; shvegsha(p) = 0._r8
       lhvegsun(p) = 0._r8 ; lhvegsha(p) = 0._r8
       etvegsun(p) = 0._r8 ; etvegsha(p) = 0._r8
       gppvegsun(p) = 0._r8 ; gppvegsha(p) = 0._r8
       vcmax25sun(p) = 0._r8 ; vcmax25sha(p) = 0._r8
       gsvegsun(p) = 0._r8 ; gsvegsha(p) = 0._r8

       do ic = 1, ncan(p)
          if (dpai(p,ic) > 0._r8) then
             laisun(p) = laisun(p) + fracsun(p,ic) * dpai(p,ic)
             laisha(p) = laisha(p) + (1._r8 - fracsun(p,ic)) * dpai(p,ic)
             lwvegsun(p) = lwvegsun(p) + lwleaf(p,ic,isun) * fracsun(p,ic) * dpai(p,ic)
             lwvegsha(p) = lwvegsha(p) + lwleaf(p,ic,isha) * (1._r8 - fracsun(p,ic)) * dpai(p,ic)
             shvegsun(p) = shvegsun(p) + shleaf(p,ic,isun) * fracsun(p,ic) * dpai(p,ic)
             shvegsha(p) = shvegsha(p) + shleaf(p,ic,isha) * (1._r8 - fracsun(p,ic)) * dpai(p,ic)
             lhvegsun(p) = lhvegsun(p) + lhleaf(p,ic,isun) * fracsun(p,ic) * dpai(p,ic)
             lhvegsha(p) = lhvegsha(p) + lhleaf(p,ic,isha) * (1._r8 - fracsun(p,ic)) * dpai(p,ic)
             etvegsun(p) = etvegsun(p) + (evleaf(p,ic,isun) + trleaf(p,ic,isun)) * fracsun(p,ic) * dpai(p,ic)
             etvegsha(p) = etvegsha(p) + (evleaf(p,ic,isha) + trleaf(p,ic,isha)) * (1._r8 - fracsun(p,ic)) * dpai(p,ic)
             fracgreen = fdry(p,ic) / (1._r8 - fwet(p,ic))
             gppvegsun(p) = gppvegsun(p) + agross(p,ic,isun) * fracsun(p,ic) * dpai(p,ic) * fracgreen
             gppvegsha(p) = gppvegsha(p) + agross(p,ic,isha) * (1._r8 - fracsun(p,ic)) * dpai(p,ic) * fracgreen
             vcmax25sun(p) = vcmax25sun(p) + vcmax25_leaf(p,ic,isun) * fracsun(p,ic) * dpai(p,ic)
             vcmax25sha(p) = vcmax25sha(p) + vcmax25_leaf(p,ic,isha) * (1._r8 - fracsun(p,ic)) * dpai(p,ic)
             gsvegsun(p) = gsvegsun(p) + gs(p,ic,isun) * fracsun(p,ic) * dpai(p,ic)
             gsvegsha(p) = gsvegsha(p) + gs(p,ic,isha) * (1._r8 - fracsun(p,ic)) * dpai(p,ic)
          end if
       end do

       ! Sunlit and shaded canopy temperatures and wind speed are weighted for sun/shade leaf area

       windveg(p) = 0._r8 ; windvegsun(p) = 0._r8 ; windvegsha(p) = 0._r8
       tlveg(p) = 0._r8 ; tlvegsun(p) = 0._r8 ; tlvegsha(p) = 0._r8
       taveg(p) = 0._r8 ; tavegsun(p) = 0._r8 ; tavegsha(p) = 0._r8
       do ic = 1, ncan(p)
          if (dpai(p,ic) > 0._r8) then
             windveg(p) = windveg(p) + wind(p,ic) * dpai(p,ic) / (laisun(p) + laisha(p))
             windvegsun(p) = windvegsun(p) + wind(p,ic) * fracsun(p,ic) * dpai(p,ic) / laisun(p)
             windvegsha(p) = windvegsha(p) + wind(p,ic) * (1._r8 - fracsun(p,ic)) * dpai(p,ic) / laisha(p)

             tlveg(p) = tlveg(p) + tleaf_mean(p,ic) * dpai(p,ic) / (laisun(p) + laisha(p))
             tlvegsun(p) = tlvegsun(p) + tleaf(p,ic,isun) * fracsun(p,ic) * dpai(p,ic) / laisun(p)
             tlvegsha(p) = tlvegsha(p) + tleaf(p,ic,isha) * (1._r8 - fracsun(p,ic)) * dpai(p,ic) / laisha(p)

             taveg(p) = taveg(p) + tair(p,ic) * dpai(p,ic) / (laisun(p) + laisha(p))
             tavegsun(p) = tavegsun(p) + tair(p,ic) * fracsun(p,ic) * dpai(p,ic) / laisun(p)
             tavegsha(p) = tavegsha(p) + tair(p,ic) * (1._r8 - fracsun(p,ic)) * dpai(p,ic) / laisha(p)
          end if
       end do

       ! Diagnose fraction of the canopy that is water stressed

       minlwp = -2._r8
       fracminlwp(p) = 0._r8

       do ic = 1, ncan(p)
          if (dpai(p,ic) > 0._r8) then
             if (lwp_mean(p,ic) <= minlwp) then
                fracminlwp(p) = fracminlwp(p) + dpai(p,ic)
             end if
          end if
       end do

       if ((lai(p) + sai(p)) > 0._r8) then
          fracminlwp(p) = fracminlwp(p) / (lai(p) + sai(p))
       end if

    end do

    end associate
  end subroutine CanopyFluxesDiagnostics

end module MLCanopyFluxesMod
