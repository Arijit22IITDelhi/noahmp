module SoilWaterTranspirationMod

!!! compute soil water transpiration factor that will be used for 
!!! stomata resistance and evapotranspiration calculations
!!! Includes peatland-specific drought and waterlogging stress conditions

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use PeatlandPhysicsMod,                only : ApplyPeatlandPhysics

  implicit none

contains

  subroutine SoilWaterTranspiration(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: None (embedded in ENERGY subroutine)
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactored code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! Modified for Peatland transpiration and waterlogging stress (Chakraborty et al. 2026)
! -------------------------------------------------------------------------

    implicit none

! in & out variables
    type(noahmp_type), intent(inout) :: noahmp

! local variables
    integer                          :: IndSoil       ! loop index
    real(kind=kind_noahmp)           :: SoilWetFac    ! temporary variable
    real(kind=kind_noahmp)           :: MinThr        ! minimum threshold to prevent division by zero
    real(kind=kind_noahmp)           :: F_wilt        ! PEAT wilting stress based on WTD
    real(kind=kind_noahmp)           :: F_log         ! PEAT waterlogging stress factor (FOXY)

! --------------------------------------------------------------------
    associate(                                                                             &
              SurfaceType               => noahmp%config%domain%SurfaceType               ,&
              ThicknessSnowSoilLayer    => noahmp%config%domain%ThicknessSnowSoilLayer    ,&
              DepthSoilLayer            => noahmp%config%domain%DepthSoilLayer            ,&
              OptSoilWaterTranspiration => noahmp%config%nmlist%OptSoilWaterTranspiration ,&
              OptRunoffSubsurface       => noahmp%config%nmlist%OptRunoffSubsurface       ,&
              OptPeatlandPhysics        => noahmp%config%nmlist%OptPeatlandPhysics        ,& 
              NumSoilLayerRoot          => noahmp%water%param%NumSoilLayerRoot            ,&
              SoilMoistureWilt          => noahmp%water%param%SoilMoistureWilt            ,&
              SoilMoistureFieldCap      => noahmp%water%param%SoilMoistureFieldCap        ,&
              SoilMatPotentialWilt      => noahmp%water%param%SoilMatPotentialWilt        ,&
              SoilMatPotentialSat       => noahmp%water%param%SoilMatPotentialSat         ,&
              SoilMoistureSat           => noahmp%water%param%SoilMoistureSat             ,&
              SoilExpCoeffB             => noahmp%water%param%SoilExpCoeffB               ,&
              SoilLiqWater              => noahmp%water%state%SoilLiqWater                ,&
              SoilTranspFac             => noahmp%water%state%SoilTranspFac               ,&
              SoilTranspFacAcc          => noahmp%water%state%SoilTranspFacAcc            ,&
              SoilMatPotential          => noahmp%water%state%SoilMatPotential            ,&
              WaterTableDepth           => noahmp%water%state%WaterTableDepth              &
             )
! ----------------------------------------------------------------------

    MinThr         = 1.0e-6
    SoilTranspFacAcc = 0.0

    ! Set peatland-specific transpiration option if peatland physics is enabled
    if ( OptPeatlandPhysics == 1 ) then
       OptSoilWaterTranspiration = 4
    endif

    if ( SurfaceType == 1 ) then

       do IndSoil = 1, NumSoilLayerRoot
          if ( OptSoilWaterTranspiration == 1 ) then  ! Noah
             SoilWetFac = (SoilLiqWater(IndSoil) - SoilMoistureWilt(IndSoil)) / &
                           (SoilMoistureFieldCap(IndSoil) - SoilMoistureWilt(IndSoil))

          else if ( OptSoilWaterTranspiration == 2 ) then  ! CLM
             SoilMatPotential(IndSoil) = max(SoilMatPotentialWilt, -SoilMatPotentialSat(IndSoil) * &
                                            (max(0.01,SoilLiqWater(IndSoil))/SoilMoistureSat(IndSoil)) ** &
                                            (-SoilExpCoeffB(IndSoil)))
             SoilWetFac = (1.0 - SoilMatPotential(IndSoil)/SoilMatPotentialWilt) / &
                           (1.0 + SoilMatPotentialSat(IndSoil)/SoilMatPotentialWilt)

          else if ( OptSoilWaterTranspiration == 3 ) then  ! SSiB
             SoilMatPotential(IndSoil) = max(SoilMatPotentialWilt, -SoilMatPotentialSat(IndSoil) * &
                                            (max(0.01,SoilLiqWater(IndSoil))/SoilMoistureSat(IndSoil)) ** &
                                            (-SoilExpCoeffB(IndSoil)))
             SoilWetFac = 1.0 - exp(-5.8*(log(SoilMatPotentialWilt/SoilMatPotential(IndSoil))))

          else if ( OptSoilWaterTranspiration == 4 ) then  ! PEAT (drought + waterlogging stress)
             if (WaterTableDepth < 0.3) then
                F_wilt = 0.0
             else if (WaterTableDepth >= 0.3 .and. WaterTableDepth < 1.15) then
                F_wilt = 1.18 * WaterTableDepth - 0.35
             else
                F_wilt = 1.0
             end if
             F_wilt = max(0.0, min(1.0, F_wilt))

             ! Apply waterlogging stress based only on WaterTableDepth
             if (WaterTableDepth >= 0.29) then
                F_log = 1.0
             else if (WaterTableDepth >= -0.35 .and. WaterTableDepth < 0.29) then
                F_log = 1.0 - max(0.0, min(0.95, (0.29 - WaterTableDepth) / 0.64))
             else
                F_log = 0.0
             end if
             F_log = max(0.0, min(1.0, F_log))

             ! Combine both drought and waterlogging stress factors
             SoilWetFac = (1.0 - F_wilt) * F_log

          end if

          SoilWetFac = min(1.0, max(0.0, SoilWetFac))

          SoilTranspFac(IndSoil) = max(MinThr, ThicknessSnowSoilLayer(IndSoil) / &
                                       (-DepthSoilLayer(NumSoilLayerRoot)) * SoilWetFac)
          SoilTranspFacAcc = SoilTranspFacAcc + SoilTranspFac(IndSoil)
       end do

       SoilTranspFacAcc = max(MinThr, SoilTranspFacAcc)
       SoilTranspFac(1:NumSoilLayerRoot) = SoilTranspFac(1:NumSoilLayerRoot) / SoilTranspFacAcc

    end if

    end associate

  end subroutine SoilWaterTranspiration

end module SoilWaterTranspirationMod
