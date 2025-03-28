module RunoffSubSurfacePeatlandMod

!!! Calculate subsurface runoff using Ivanov-based Peatland Runoff Scheme - Enabled by (Chakraborty & Bechtold, 2025)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use WaterTableEquilibriumMod, only : WaterTableEquilibrium

  implicit none

contains

  subroutine RunoffSubSurfacePeatland(noahmp)

! ------------------------ Code history --------------------------------------------------
! Modified to include Peatland-specific Ivanov-based runoff scheme (Chakraborty & Bechtold, 2025)
! ----------------------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

    ! Define double precision kind parameter
    integer, parameter :: dp = kind(1.0d0)

    ! Declare local peatland-specific parameters
    real(dp) :: Ksz_zero, m_Ivanov, v_slope
    real(dp) :: Ta, BFLOW

! --------------------------------------------------------------------
    associate(                                                           &
              SoilImpervFracMax => noahmp%water%state%SoilImpervFracMax ,& ! in,    maximum soil imperviousness fraction
              WaterTableDepth   => noahmp%water%state%WaterTableDepth   ,& ! out,   water table depth [m]
              FSW_change         => noahmp%water%state%FSW_change        ,& ! inout, 
              RunoffSubsurface  => noahmp%water%flux%RunoffSubsurface    & ! out,   subsurface runoff [mm/s] 
             )
! ----------------------------------------------------------------------

    ! Compute equilibrium water table depth
    ! For very shallow water table detph, use PEATCLSM approximation outside of this routine
    if (WaterTableDepth>0.1) then
       call WaterTableEquilibrium(noahmp)
    endif

    ! ------------------------------------------
    ! Option 9: Ivanov-based Peatland Runoff Scheme (Chakraborty & Bechtold, 2025)
    ! ------------------------------------------

    ! Assign parameter values for peatland runoff scheme
    Ksz_zero = 3165.38_dp      ! Saturated hydraulic conductivity [m^2/s]
    m_Ivanov = 2.06_dp         ! Ivanov exponent
    v_slope = 1.5e-08_dp       ! Slope factor for runoff generation [unitless]

    ! Compute transmissivity function (Ta) [m^2/s]
    Ta = (Ksz_zero * (24.5_dp + 100.0_dp * max(-0.244999_dp, WaterTableDepth))**(1.0_dp - m_Ivanov)) / &
         (100.0_dp * (m_Ivanov - 1.0_dp))

    ! Compute baseflow (BFLOW) in mm/s
    BFLOW = v_slope * Ta * 1000.0_dp  ! Convert from m/s to mm/s

    ! Compute subsurface runoff using Peatland-specific equation
    RunoffSubsurface = (1.0_dp - SoilImpervFracMax) * BFLOW
    
    ! Set FSW_change to zero for following calculations in SoilWaterMain and WaterBalanceError Check
    FSW_change = 0.0

    end associate

  end subroutine RunoffSubSurfacePeatland

end module RunoffSubSurfacePeatlandMod
