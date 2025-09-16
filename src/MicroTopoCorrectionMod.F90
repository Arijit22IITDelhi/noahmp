module MicroTopoCorrectionMod

!!! Enabled by (Chakraborty & Bechtold, 2025)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use WaterTableEquilibriumPeatMod,      only : WaterTableEquilibriumPeat

  implicit none

contains

  subroutine MicroTopoCorrection(noahmp)

! ------------------------ Code history --------------------------------------------------
! SY from WaterTableDepth following PEATCLSM routines (Chakraborty & Bechtold, 2025)
! ----------------------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

    ! Define double precision kind parameter
    integer, parameter :: dp = kind(1.0d0)

    ! Declare local peatland-specific parameters
    real(dp) :: SySoil, bf1, bf2, PEATCLSM_ZBARMAX_4_SYSOIL
    real(dp) :: catdef, ars1, ars2, ars3

! --------------------------------------------------------------------
    associate(                                                           &
              InfilRateSfc           => noahmp%water%flux%InfilRateSfc    ,& ! in,   infiltration rate at surface [m/s]
              SoilTimeStep           => noahmp%config%domain%SoilTimeStep ,& ! in,    noahmp soil time step [s]
              f_soil                 => noahmp%water%state%f_soil         ,& ! inout, fraction of flux in and out of soil [-]
              AR1                    => noahmp%water%state%AR1         ,& ! inout, fraction of flux in and out of soil [-]
              WaterTableDepth   => noahmp%water%state%WaterTableDepth      & ! inout,   water table depth [m]
             )
! ----------------------------------------------------------------------

    ! Assign parameter values for peatland
    ! CO NN
    bf1 = 1.8064707e+02
    bf2 = 2.4298242e-01
    ars1 = -8.9250673e-03
    ars2 = 5.7296452e-02
    ars3 = 2.2656331e-03
    PEATCLSM_ZBARMAX_4_SYSOIL = 0.45    ! [m]

    ! Compute transmissivity function (Ta) [m^2/s]
    SySoil = (2.*bf1*min(max(WaterTableDepth,0.),PEATCLSM_ZBARMAX_4_SYSOIL) + 2.*bf1*bf2)/1000.
    catdef = max(((WaterTableDepth + bf2)**2 - 1.0E-20), 0._dp) * bf1
    AR1 = (1.+ars1*(catdef))/(1.+ars2*(catdef)+ars3*(catdef)**2)
    
    if (WaterTableDepth>0.1) then
        f_soil = MAX(MIN(1.0,(1.-AR1)*SySoil/(AR1+(1.-AR1)*SySoil)),0.0)
    else
        f_soil = 0.0
    endif

    end associate

  end subroutine MicroTopoCorrection

end module MicroTopoCorrectionMod
