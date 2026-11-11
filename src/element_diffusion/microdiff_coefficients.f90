!----------------------------------------------------------------------
! microdiff_coefficients
!----------------------------------------------------------------------
! microdiff's drift and diffusion coefficients for one species at the
! points of one equal grid. The species' diffusion velocity is (Thoul et
! al. 1994)
!     w = A_P d ln P/dr + A_T d ln T/dr + A_H d ln C_H/dr + A_D d ln C_D/dr,
! the A's from the Burgers equations (thoul_solver) with each pair's
! Coulomb logarithm (coulomb_logarithm), for the four Thoul columns:
! hydrogen, helium, the species (column 3, or iron when hydrogen is
! solved for) and the electrons. In Bahcall & Loeb units the
! coefficients of dX/dt = (1/(rho r^2)) d/dr [D1 + D2 dX/dr] are
!     D1 = f r^2 T^(5/2)/ln Lambda (A_P d ln P/dr + A_T d ln T/dr
!          + A_H d ln C_H/dr) X,
!     D2 = f r^2 T^(5/2)/ln Lambda A_D,
! with A_T d ln T/dr = A_T del d ln P/dr (del = d ln T/d ln P) and the
! signs of the Burgers solution's coefficients reversed. f is
! settling_y_factor for hydrogen (column 1) and settling_z_factor for
! column 3, the metals and the light elements alike (docs/bug_sweep.md
! M250). The hydrogen term drops out of hydrogen's own D1 (it is the
! self term, in D2); its gradient, d ln C_H/dr centred on the grid, is
! computed by the hydrogen call and read by the others.
module microdiff_coefficients
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use composition_layout, only: json
      use element_diffusion_configuration, only: diffusion_settings
      use thoul_solver, only: thoul_diffusion
      use coulomb_logarithm, only: compute_concentrations, compute_coulomb_logarithms
      use math_lib
      implicit none
      private
      public :: compute_microdiff_coefficients

! the Thoul columns: hydrogen, helium, the species, the electrons
      integer, parameter :: num_species = 4
! ln Lambda in the Bahcall & Loeb time unit (bahcall_loeb_units): the
! F77's single-precision 2.2 (2.2000000476837158 as a double), which
! microdiff keeps; the original routines divide by 2.2d0 (docs/bug_sweep.md
! M257: making them one moves results)
      real(dp), parameter :: bl_coulomb_logarithm = real(2.2, dp)
! the least mass fraction the Burgers equations take (none may be zero)
      real(dp), parameter :: mass_fraction_floor = 1.0d-24

contains

! ---------------------------------------------------------------
! The drift (drift_coeff, D1) and diffusion (diffusion_coeff, D2)
! coefficients at the num_points points of grid, for the species in
! Thoul column species_col (atomic weight and charge as given for column
! 3), from the columns' mass fractions species_fraction(column, point).
! hydrogen_dlnc_dr: d ln C_H/dr at the points, written by the hydrogen
! call (species_col 1, zero at the grid's ends) and read by the others.
! Where the species is absent at a point and the next (the depleted
! zone, the first such point excepted) both coefficients are zero.
! ierr: ludcmp's status from a singular Burgers system.
subroutine compute_microdiff_coefficients(num_points, species_fraction, grid, &
     drift_coeff, diffusion_coeff, hydrogen_dlnc_dr, atomic_weight_diffused, &
     atomic_charge_diffused, species_col, ierr)
      use microdiff_mapping, only: microdiff_grid
      use species_table, only: thoul_h1, thoul_he4, thoul_electron, &
           thoul_col_h, thoul_col_he, thoul_col_metal, thoul_col_electron
      use bahcall_loeb_units, only: bl_temp_scale
      integer, intent(in) :: num_points
      real(dp), intent(in) :: species_fraction(:,:)
      type(microdiff_grid), intent(in) :: grid
      real(dp), intent(out) :: drift_coeff(:), diffusion_coeff(:)
      real(dp), intent(inout) :: hydrogen_dlnc_dr(:)
      real(dp), intent(in) :: atomic_weight_diffused, atomic_charge_diffused
      integer, intent(in) :: species_col
      integer, intent(out) :: ierr
      real(dp) :: atomic_weight(num_species), atomic_charge(num_species), &
           mass_frac(num_species), concentration(num_species), &
           coulomb_log(num_species,num_species), pressure_coeff(num_species), &
           temp_coeff(num_species), conc_coeff(num_species,num_species)
! per point: f r^2 T^(5/2)/ln Lambda, the four terms of the velocity
! (signs reversed; the temperature term times del), and C_H
      real(dp) :: coeff_scale(json), pressure_term(json), temp_term(json), &
           hydrogen_term(json), self_term(json), hydrogen_concentration(json)
      real(dp) :: scale_factor, radius_bl, temperature_bl, dlnc_dr
      integer :: i

      ierr = 0
      atomic_weight(thoul_col_h) = thoul_h1%weight
      atomic_weight(thoul_col_he) = thoul_he4%weight
      atomic_weight(thoul_col_metal) = atomic_weight_diffused
      atomic_weight(thoul_col_electron) = thoul_electron%weight
      atomic_charge(thoul_col_h) = thoul_h1%charge
      atomic_charge(thoul_col_he) = thoul_he4%charge
      atomic_charge(thoul_col_metal) = atomic_charge_diffused
      atomic_charge(thoul_col_electron) = thoul_electron%charge
      scale_factor = diffusion_settings%fgrz
      if (species_col == thoul_col_h) scale_factor = diffusion_settings%fgry

! The Burgers solution at each point.
      do i = 1, num_points
         mass_frac(thoul_col_h) = species_fraction(thoul_col_h,i)
         mass_frac(thoul_col_he) = species_fraction(thoul_col_he,i)
         mass_frac(thoul_col_metal) = species_fraction(thoul_col_metal,i)
         if (mass_frac(thoul_col_h) < mass_fraction_floor) mass_frac(thoul_col_h) = mass_fraction_floor
         if (mass_frac(thoul_col_he) < mass_fraction_floor) mass_frac(thoul_col_he) = mass_fraction_floor
         if (mass_frac(thoul_col_metal) < mass_fraction_floor) &
              mass_frac(thoul_col_metal) = mass_fraction_floor
         call compute_concentrations(atomic_weight, atomic_charge, mass_frac, concentration)
         if (species_col == thoul_col_h) hydrogen_concentration(i) = concentration(thoul_col_h)
! In the depleted zone, nothing to solve.
         if (species_fraction(species_col,i) == 0.0d0 .and. i /= num_points) then
            if (species_fraction(species_col,i+1) == 0.0d0) then
               coeff_scale(i) = 0.0d0
               pressure_term(i) = 0.0d0
               temp_term(i) = 0.0d0
               hydrogen_term(i) = 0.0d0
               self_term(i) = 0.0d0
               cycle
            end if
         end if
         call compute_coulomb_logarithms(atomic_weight, atomic_charge, concentration, &
              grid%density(i), grid%temperature(i), coulomb_log)
         call thoul_diffusion(atomic_weight, atomic_charge, mass_frac, coulomb_log, &
              pressure_coeff, temp_coeff, conc_coeff, ierr)
         if (ierr /= 0) return
         radius_bl = grid%radius(i)
         temperature_bl = grid%temperature(i)*bl_temp_scale
         coeff_scale(i) = scale_factor*radius_bl**2*pow(temperature_bl, 2.5d0)/bl_coulomb_logarithm
         pressure_term(i) = -pressure_coeff(species_col)
         temp_term(i) = -temp_coeff(species_col)*grid%del_grad(i)
         hydrogen_term(i) = -conc_coeff(species_col,thoul_col_h)
         self_term(i) = -conc_coeff(species_col,species_col)
      end do

! The coefficients, with the hydrogen gradient.
      do i = 1, num_points
         if (species_col == thoul_col_h) then
            hydrogen_dlnc_dr(i) = centred_log_gradient(hydrogen_concentration, &
                 grid%radius, i, num_points)
            dlnc_dr = 0.0d0
         else
            dlnc_dr = hydrogen_dlnc_dr(i)
         end if
         drift_coeff(i) = coeff_scale(i)*(grid%dlnp_dr(i)*(pressure_term(i)+temp_term(i)) &
              + dlnc_dr*hydrogen_term(i))*species_fraction(species_col,i)
         diffusion_coeff(i) = coeff_scale(i)*self_term(i)
      end do
end subroutine compute_microdiff_coefficients

! ---------------------------------------------------------------
! d ln c/dr at point i of an equally spaced grid of n points (radius),
! the mean of the two one-sided differences, each over the spacing
! r(i+1) - r(i); zero at the grid's ends.
pure function centred_log_gradient(c, radius, i, n) result(gradient)
      real(dp), intent(in) :: c(:), radius(:)
      integer, intent(in) :: i, n
      real(dp) :: gradient
      real(dp) :: log_c, log_c_next, log_c_prev, spacing

      if (i == 1 .or. i == n) then
         gradient = 0.0d0
         return
      end if
      log_c = log(c(i))
      log_c_next = log(c(i+1))
      log_c_prev = log(c(i-1))
      spacing = radius(i+1) - radius(i)
      gradient = 0.5d0*((log_c - log_c_prev)/spacing + (log_c_next - log_c)/spacing)
end function centred_log_gradient

end module microdiff_coefficients

! Sources
! - G. Somers (05/2015), after the original routines' coefficients
!   (MHP 8/94); JvS (01/2026) applied settling_y_factor and
!   settling_z_factor.
! - Thoul, A. A., Bahcall, J. N. & Loeb, A. 1994, ApJ 421, 828.
