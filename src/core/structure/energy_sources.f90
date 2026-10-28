!----------------------------------------------------------------------
! energy_sources
!----------------------------------------------------------------------
! The extra energy sources of the luminosity equation
! (docs/TARGET_ARCHITECTURE.md section 4.5a): heating or cooling that is
! not nuclear, neutrino or gravothermal, such as tidal or magnetic
! heating, an injected energy, or a test source. A source is a
! subroutine with the energy_source_eval interface below, registered by
! name with add_energy_source before the run starts. A program that
! embeds YREC registers its sources and then calls run_yrec;
! examples/other_energy is a complete one.
!
! At every shell of the interior model (inside the fitting point), the
! structure equations add the sum of the sources to the luminosity
! equation,
!   dL/dm = eps_nuc - eps_nu + eps_grav - d(e_rot)/dt + eps_extra,
! and the partials d eps_extra/dln rho (at constant T) and
! d eps_extra/dln T (at constant rho) to its Jacobian, as the nuclear
! term's are. The dependence on r and L is not linearized. eps_extra is
! in erg/g/s; positive heats.
!
! What a source sees (energy_point): the shell's enclosed mass and its
! own mass (g), log10 of T, P, rho and r (cgs), L in Lsun, the mass
! fractions (composition_layout's slots), and the step dt (s; zero when
! the model is not being aged).
!
! Limits:
!  * The envelope above the fitting point is integrated at constant L,
!    so a source there has no effect. Keep sources inside the interior
!    (the fitting mass is printed in the run log).
!  * Within an iteration the sources are evaluated on the current
!    trial model; a source with a steep dependence on the structure
!    should supply honest partials, or the Newton iteration slows.
!
! With no source registered nothing is evaluated, and every result is
! the same, bit for bit, as without this module. The registry is the
! process's configuration, not model state: it outlives run_yrec (a
! second run in the same process uses the same sources) until
! clear_energy_sources is called.
module energy_sources
      use composition_layout, only: n_species_extended
      implicit none
      private
      public :: energy_point, energy_source_eval, max_energy_sources, &
           add_energy_source, clear_energy_sources, num_energy_sources, &
           energy_source_name, evaluate_energy_sources

! What a source sees at one shell (see the header).
      type :: energy_point
         double precision :: mass = 0.0d0, shell_mass = 0.0d0
         double precision :: log_temperature = 0.0d0, log_pressure = 0.0d0, &
              log_density = 0.0d0, log_radius = 0.0d0
         double precision :: luminosity_lsun = 0.0d0
         double precision :: xa(n_species_extended) = 0.0d0
         double precision :: dt = 0.0d0
      end type energy_point

      abstract interface
! eps (erg/g/s) at the point, its derivatives with respect to ln rho at
! constant T and ln T at constant rho, and ierr (nonzero stops the
! iteration with an error).
         subroutine energy_source_eval(pt, eps, deps_dlnrho, deps_dlnt, ierr)
            import :: energy_point
            type(energy_point), intent(in) :: pt
            double precision, intent(out) :: eps, deps_dlnrho, deps_dlnt
            integer, intent(out) :: ierr
         end subroutine energy_source_eval
      end interface

      integer, parameter :: max_energy_sources = 8

      type :: energy_source
         character(len=32) :: name = ''
         procedure(energy_source_eval), pointer, nopass :: eval => null()
      end type energy_source

! The registered sources (the process's configuration; see the header).
      type(energy_source), save :: sources(max_energy_sources)
      integer, save :: num_sources = 0

contains

! ---------------------------------------------------------------
! Register eval under name. ierr = 1 when the list is full or the name
! is already registered.
subroutine add_energy_source(name, eval, ierr)
      character(len=*), intent(in) :: name
      procedure(energy_source_eval) :: eval
      integer, intent(out) :: ierr
      integer :: i

      ierr = 0
      do i = 1, num_sources
         if (sources(i)%name == name) then
            write(*,'(a)') 'add_energy_source: '//trim(name)//' is already registered'
            ierr = 1
            return
         end if
      end do
      if (num_sources >= max_energy_sources) then
         write(*,'(a,i0,a)') 'add_energy_source: no room for ', &
              max_energy_sources + 1, ' sources'
         ierr = 1
         return
      end if
      num_sources = num_sources + 1
      sources(num_sources)%name = name
      sources(num_sources)%eval => eval
end subroutine add_energy_source

! ---------------------------------------------------------------
! Forget every registered source.
subroutine clear_energy_sources()
      integer :: i

      do i = 1, num_sources
         sources(i)%name = ''
         sources(i)%eval => null()
      end do
      num_sources = 0
end subroutine clear_energy_sources

! ---------------------------------------------------------------
integer function num_energy_sources()
      num_energy_sources = num_sources
end function num_energy_sources

! ---------------------------------------------------------------
! The name of source i (blank when there is none).
function energy_source_name(i) result(name)
      integer, intent(in) :: i
      character(len=32) :: name
      name = ''
      if (i >= 1 .and. i <= num_sources) name = sources(i)%name
end function energy_source_name

! ---------------------------------------------------------------
! The sum of the registered sources at pt and of their partials.
subroutine evaluate_energy_sources(pt, eps, deps_dlnrho, deps_dlnt, ierr)
      type(energy_point), intent(in) :: pt
      double precision, intent(out) :: eps, deps_dlnrho, deps_dlnt
      integer, intent(out) :: ierr
      double precision :: e, de_drho, de_dt
      integer :: i

      eps = 0.0d0
      deps_dlnrho = 0.0d0
      deps_dlnt = 0.0d0
      ierr = 0
      do i = 1, num_sources
         call sources(i)%eval(pt, e, de_drho, de_dt, ierr)
         if (ierr /= 0) then
            write(*,'(a)') 'energy source '//trim(sources(i)%name)//' failed'
            return
         end if
         eps = eps + e
         deps_dlnrho = deps_dlnrho + de_drho
         deps_dlnt = deps_dlnt + de_dt
      end do
end subroutine evaluate_energy_sources

end module energy_sources
