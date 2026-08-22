!----------------------------------------------------------------------
! hsubp
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original hsubp.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Calculates the local pressure scale height at the edge of a
! convective region (PD Dec. 1984; 2/91 MHP revised for YREC/Mark6).
! Calls the equation of state to get the mean molecular weight at the
! edge zone, then returns the density-based pressure scale height
! (pscahe).
subroutine hsubp(composition, density, pressure, radius, mass, &
     temperature, edge_zone, pscahe)

      use const_lib
      use eos_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: composition(15,json), density(json), &
           pressure(json), radius(json), mass(json), temperature(json)
      integer, intent(in) :: edge_zone
      double precision, intent(out) :: pscahe


      double precision :: fxion(3)
      save

      logical :: lderiv, latmo
      integer :: ksaha
      double precision :: hydrogen_fraction, metal_fraction, log10_pressure, &
           log10_temperature, log10_density
      double precision :: pressure_out, temperature_out, density_out
      double precision :: beta, betai, beta14, rmu, amu, emu, eta, qdt, &
           qdp, qcp, dela, qdtt, qdtp, qat, qap, qcpt, qcpp
      double precision :: log10_mass, log10_radius
      double precision :: pscap

! *** CALCULATES THE LOCAL PRESSURE SCALE HEIGHT AT THE EDGE OF A
!     CONVECTIVE REGION (PD DEC. 1984) 2/91 MHP REVISED FOR YREC/MARK6***
!
!  CALL EQUATION OF STATE TO DETERMINE THE MEAN MOLECULAR WEIGHT.
      lderiv = .false.
      latmo = .true.
      ksaha = 0
      hydrogen_fraction = composition(1,edge_zone)
      metal_fraction = composition(3,edge_zone)
      log10_pressure = pressure(edge_zone)
      log10_temperature = temperature(edge_zone)
      log10_density = density(edge_zone)
      call eos_get(log10_temperature,temperature_out,log10_pressure, &
           pressure_out,log10_density,density_out,hydrogen_fraction, &
           metal_fraction,beta,betai,beta14,fxion, &
           rmu,amu,emu,eta,qdt,qdp,qcp,dela,qdtt,qdtp, &
           qat,qap,qcpt,qcpp,lderiv,latmo,ksaha, &
           composition_at_zone=composition(:,edge_zone))
!  COMPUTE PRESSURE SCALE HEIGHT.
      log10_mass = mass(edge_zone)
      log10_radius = radius(edge_zone)
      pscap = rmu*exp(ln10*(log10_temperature-cgl-log10_mass+log10_radius+ &
           log10_radius))
      pscahe = exp(ln10*(log10_pressure + 2.0d0*log10_radius - &
           log10_density - cgl - log10_mass))
      return
end subroutine hsubp
