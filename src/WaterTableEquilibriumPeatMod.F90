module WaterTableEquilibriumPeatMod

!!! Calculate equilibrium water table depth (Bechtold et al., 2019)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  
  implicit none

contains

  subroutine WaterTableEquilibriumPeat(noahmp)
  
! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: ZWTEQ
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactored:    C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! This version:  fix Campbell/hydrostatic formulation + robust bisection (2025)
! ----------------------------------------------------------------------------------------

    implicit none
    
    type(noahmp_type), intent(inout) :: noahmp

    integer, parameter               :: NumSoilFineLy = 4000
    integer                          :: i, iter
    real(kind=kind_noahmp)           :: dz_fine, zmax
    real(kind=kind_noahmp)           :: deficit_target, deficit_mid
    real(kind=kind_noahmp)           :: zwt_lo, zwt_hi, zwt_mid
    real(kind=kind_noahmp)           :: ae, bb, thetas
    real(kind=kind_noahmp), parameter:: tol_def  = 1.0e-6_kind_noahmp
    real(kind=kind_noahmp), parameter:: tol_zwt  = 1.0e-4_kind_noahmp
    integer,          parameter      :: max_iter = 60
    real(kind=kind_noahmp)           :: z, psi_abs, theta_z, deficit_acc  ! (z, theta_z, deficit_acc not strictly needed here)
! -----------------------------------------------------------------------------------------------------------------------------
    associate(                                                                        &
      NumSoilLayer           => noahmp%config%domain%NumSoilLayer           ,& ! in
      DepthSoilLayer         => noahmp%config%domain%DepthSoilLayer         ,& ! in  layer-bottom depths [m], negative downward
      ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! in
      SoilLiqWater           => noahmp%water%state%SoilLiqWater             ,& ! in  [m3/m3]
      SoilMoistureSat        => noahmp%water%param%SoilMoistureSat          ,& ! in  θ_s [m3/m3]
      SoilMatPotentialSat    => noahmp%water%param%SoilMatPotentialSat      ,& ! in  ψ_e (air-entry), typically negative [m]
      SoilExpCoeffB          => noahmp%water%param%SoilExpCoeffB            ,& ! in  Campbell b
      WaterTableDepth        => noahmp%water%state%WaterTableDepth            & ! out z_wt from surface [m], positive downward
    )
! ------------------------------------------------------------------------------------------------------------------------------

      ! constants for the single soil type used
      thetas = SoilMoistureSat(1)
      ae     = abs(SoilMatPotentialSat(1))   ! air-entry suction head [m], positive
      bb     = SoilExpCoeffB(1)

      ! target deficit from current coarse profile
      deficit_target = 0.0_kind_noahmp
      do i = 1, NumSoilLayer
        deficit_target = deficit_target + ( thetas - SoilLiqWater(i) ) * ThicknessSnowSoilLayer(i)
      end do

      if (deficit_target <= 0.0_kind_noahmp) then
        WaterTableDepth = 0.0_kind_noahmp

      else
        ! integration domain: 3x modeled soil depth
        zmax = 3.0_kind_noahmp * ( -DepthSoilLayer(NumSoilLayer) )
        if (zmax <= 0.0_kind_noahmp) then
          WaterTableDepth = 0.0_kind_noahmp

        else
          dz_fine = zmax / real(NumSoilFineLy, kind_noahmp)

          ! bracket the solution
          zwt_lo = 0.0_kind_noahmp
          zwt_hi = zmax

          if ( deficit_from_zwt(zwt_hi) < deficit_target ) then
            ! target exceeds capacity: deepest WT
            WaterTableDepth = zwt_hi
          else
            ! bisection
            do iter = 1, max_iter
              zwt_mid    = 0.5_kind_noahmp * (zwt_lo + zwt_hi)
              deficit_mid = deficit_from_zwt(zwt_mid)

              if ( abs(deficit_mid - deficit_target) <= tol_def .or. (zwt_hi - zwt_lo) <= tol_zwt ) then
                WaterTableDepth = zwt_mid
                exit
              end if

              if (deficit_mid > deficit_target) then
                zwt_hi = zwt_mid   ! too deep -> move shallower
              else
                zwt_lo = zwt_mid   ! too shallow -> move deeper
              end if

              if (iter == max_iter) WaterTableDepth = zwt_mid
            end do
          end if
        end if
      end if
    end associate

  contains

    function deficit_from_zwt(zwt) result(def)
      implicit none
      real(kind=kind_noahmp), intent(in) :: zwt
      real(kind=kind_noahmp)             :: def
      integer                            :: k
      real(kind=kind_noahmp)             :: zz, th, psi_abs_local

      def = 0.0_kind_noahmp
      do k = 1, NumSoilFineLy
        zz = real(k,kind_noahmp) * dz_fine
        if (zz >= zwt - ae) then
          th = thetas
        else
          psi_abs_local = zwt - zz
          th = thetas * ( (psi_abs_local / ae) ** ( -1.0_kind_noahmp / bb ) )
          if (th > thetas) th = thetas
          if (th < 0.0_kind_noahmp) th = 0.0_kind_noahmp
        end if
        def = def + ( thetas - th ) * dz_fine
      end do
    end function deficit_from_zwt

  end subroutine WaterTableEquilibriumPeat

end module WaterTableEquilibriumPeatMod
