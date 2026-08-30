!----------------------------------------------------------------------
! getmodel2
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original getmodel2.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! GETMODEL2 - Read a Model2 format stellar model file into memory
!
! llp  4/16/03
!
! Dummy-argument names for the fields shared with the YREC7 format
! (the first three lines of the signature) are chosen to match
! write_yrec7.f90, the writer for the sibling YREC7 format that shares
! the same header/envelope/per-shell record layout.
subroutine read_model2(mixing_length, timestep_yr, trial_sign_flag, &
     iread, core_cz_top_index, envelope_cz_bottom_index, &
     use_extended_composition, rotation_active, atm_code, eos_code, &
     hik_code, use_diffusion_y, use_diffusion_z, &
     disk_locking_active, instability_transport_active, &
     use_wind_torque, alok_code, core_overshoot_active, &
     envelope_overshoot_active, lovstm, use_pure_z_table, &
     use_semiconvection, disk_pressure, disk_temperature, &
     wind_saturation_omega)
      use star_info_lib, only: star, i_lum_3alpha, i_lum_cno, i_lum_grav, i_lum_neu, i_lum_pp1, i_lum_pp2, i_lum_pp3, json
! First three lines above are YREC7 inputs
! Last two lines are MODEL2 add-ons
! MHP 4/25 chanted LOK to ALOK to avoid variable name conflicts

! Several of this subroutine's own dummy arguments below (values read
! from the model file -- rotation_active, envelope_overshoot_active,
! instability_transport_active, core_overshoot_active, lovstm, use_semiconvection, etc.) happen
! to share names with unrelated const_lib runtime-config module
! variables, so `use, only:` the one member actually needed here
! rather than a blanket `use const_lib` (which would conflict with the
! dummy-argument declarations below). Same treatment as
! io/read_yrec7.f90, this file's sibling reader.
! (solar_luminosity_cgs now comes from star% -- 2026 phase-A
! eviction; the former `use const_lib, only:` import is gone.)
      implicit none

      double precision, intent(out) :: mixing_length
      double precision, intent(out) :: timestep_yr
      double precision, intent(out) :: trial_sign_flag
      integer, intent(in) :: iread
      integer, intent(out) :: core_cz_top_index, envelope_cz_bottom_index
      logical, intent(out) :: use_extended_composition, rotation_active
! atm_code/eos_code/hik_code/alok_code/star%initial_composition_code are short physics-
! option codes read verbatim from the model header; their exact
! encoding is not interpreted in this file, so names are conservative.
      character*4, intent(out) :: atm_code
      character*6, intent(out) :: eos_code
      character*4, intent(out) :: hik_code
      logical, intent(out) :: use_diffusion_y, use_diffusion_z, &
           disk_locking_active, instability_transport_active, use_wind_torque
      character*4, intent(out) :: alok_code
      logical, intent(out) :: core_overshoot_active, envelope_overshoot_active, lovstm, &
           use_pure_z_table, use_semiconvection
      double precision, intent(out) :: disk_pressure, disk_temperature, &
           wind_saturation_omega


      character*5 :: header_keyword
      integer :: envelope_record_number
      integer :: i, j, k
      double precision :: max_luminosity_component

      rewind(iread)

! GET THE HEADER RECORD
      read(iread,10) header_keyword,star%model_number,star%nz,star%star_mass, &
           star%log_Teff,star%log_L,star%log_total_mass,star%dage, &
           timestep_yr,star%log_mass(1),star%log_mass(star%nz)
 10   format(A5,2I8,5F16.12,1PE18.10,0P2F16.12)

! Get the flags describing the input physics and the information that
! is stored in the model (rotating or not?, extendec compotition, etc.)
      read(iread,30) core_cz_top_index,envelope_cz_bottom_index, &
           mixing_length,eos_code,atm_code,alok_code,hik_code, &
           use_pure_z_table,star%initial_composition_code,use_extended_composition, &
           use_diffusion_y,use_diffusion_z,use_semiconvection,core_overshoot_active, &
           envelope_overshoot_active,lovstm,rotation_active, &
           instability_transport_active,use_wind_torque,disk_locking_active, &
           disk_temperature,disk_pressure,wind_saturation_omega
   30 format(2I8,F16.10,1X,A6,1X,3(A4,1X),L1,1X,A4,1X,11(L1,1X), &
           3(1PE18.10))

! Get the  LUMINOSITIES---PP(1-2-3)-CNO-HE-NU-GRAV
! CONVERT TO SOLAR LUMINOSITIES
      read(iread,40) (star%luminosity_breakdown(j),j=1,7)
!   40 FORMAT('TLUMX',5X,1P7E17.9)
   40 format(10X,1P7E17.9)

! If TLUMX are in ergs, convert to solar units.  Decide by
! comparing to 10**20.  IF larger, divide by CLSUN.
      max_luminosity_component = dmax1(star%luminosity_breakdown(i_lum_pp1), &
           star%luminosity_breakdown(i_lum_pp2),star%luminosity_breakdown(i_lum_pp3), &
           star%luminosity_breakdown(i_lum_cno),star%luminosity_breakdown(i_lum_3alpha), &
           dabs(star%luminosity_breakdown(i_lum_neu)),star%luminosity_breakdown(i_lum_grav))
      if (max_luminosity_component.gt.1.0D20) then
       do j = 1,7
          star%luminosity_breakdown(j) = star%luminosity_breakdown(j)/star%solar_luminosity_cgs
         enddo
      endif

! Get the ENVELOPE DATA
! FTRI is 1,normally.  It is set to -1  if any of the record numbers
! for the envelope triangle records was set to -1 by WRTLST.
      trial_sign_flag = 1D0
      do i = 1,3
         read(iread,50)envelope_record_number,star%trial_log_temperature(i), &
              star%trial_log_luminosity(i),star%fit_point_pressure(i), &
              star%fit_point_temperature(i),star%fit_point_radius(i), &
              (star%envelope_fit_coeffs(3*i-3+j),j=1,3)
         if (envelope_record_number .lt. 0) trial_sign_flag = -1D0
      enddo
!   50 FORMAT('ENV',I2,5F16.12,1P3E20.12)
   50 format(4X,I2,5F16.12,1P3E20.12)

! Skip over column headings for per-shell data
      read(iread,*)
      read(iread,*)
      read(iread,*)

! Read per-shell data (HENYEY POINGS)
      do i = 1,star%nz
         read(iread,60) k,star%log_mass(i),star%logR(i),star%luminosity_lsun(i), &
              star%logP(i),star%logT(i),star%logRho(i), &
              star%omega(i),star%convective_flag(i),(star%xa(j,i),j=1,15)
       enddo
 60   format(I6,0P2F18.14,1PE24.16,0P3F18.14,1PE24.16,1x,L1, &
           3(0PF12.9),12(0PE16.8))

      return
end subroutine read_model2
