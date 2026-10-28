!----------------------------------------------------------------------
! yrec_extra_heat
!----------------------------------------------------------------------
! An example of an extra energy source (src/core/structure/energy_sources.f90):
! YREC with a uniform heating rate in a shell of the star, the kind of
! term a tidal or magnetic heating model adds. The program registers the
! source, then runs YREC exactly as the yrec executable does, on the
! deck named on its command line:
!
!   ./yrec_extra_heat solar_heated.inlist
!
! The source's parameters come from the namelist &extra_heat in
! extra_heat.nml (in the working directory):
!   eps_extra     the heating rate (erg/g/s; negative cools)
!   m_inner_msun  the shell's inner edge (Msun)
!   m_outer_msun  the shell's outer edge (Msun)
! Inside the shell the rate is eps_extra; outside it is zero. It does
! not depend on rho or T, so both partials are zero.
!
! To write a source of your own, copy this file and change uniform_heat,
! keeping its argument list (energy_source_eval in energy_sources.f90).
! A source must be a module procedure, as here, not an internal one.

! The source and its parameters.
module extra_heat_source
      use energy_sources, only: energy_point
      implicit none
      private
      public :: read_extra_heat, uniform_heat

! the solar mass (g), YREC's default solar_mass_cgs
      double precision, parameter :: msun_g = 1.9891d33
      double precision, save :: eps_extra = 0.0d0, m_inner_msun = 0.0d0, &
           m_outer_msun = 0.0d0

contains

! The parameters, from &extra_heat in file.
subroutine read_extra_heat(file, ierr)
      character(len=*), intent(in) :: file
      integer, intent(out) :: ierr
      namelist /extra_heat/ eps_extra, m_inner_msun, m_outer_msun
      integer :: unit

      open(newunit=unit, file=file, status='old', action='read', iostat=ierr)
      if (ierr /= 0) then
         write(*,'(a)') 'yrec_extra_heat: cannot open '//trim(file)
         return
      end if
      read(unit, nml=extra_heat, iostat=ierr)
      close(unit)
      if (ierr /= 0) then
         write(*,'(a)') 'yrec_extra_heat: cannot read &extra_heat from '//trim(file)
         return
      end if
      write(*,'(a,1pe11.3,a,0pf8.4,a,f8.4,a)') ' extra heat: ', eps_extra, &
           ' erg/g/s between ', m_inner_msun, ' and ', m_outer_msun, ' Msun'
end subroutine read_extra_heat

! eps_extra inside the shell, zero outside; no dependence on rho or T.
subroutine uniform_heat(pt, eps, deps_dlnrho, deps_dlnt, ierr)
      type(energy_point), intent(in) :: pt
      double precision, intent(out) :: eps, deps_dlnrho, deps_dlnt
      integer, intent(out) :: ierr

      ierr = 0
      eps = 0.0d0
      if (pt%mass >= m_inner_msun*msun_g .and. pt%mass <= m_outer_msun*msun_g) &
           eps = eps_extra
      deps_dlnrho = 0.0d0
      deps_dlnt = 0.0d0
end subroutine uniform_heat

end module extra_heat_source

! Register the source, then run YREC on the deck named on the command
! line (run_yrec reads it, as the yrec executable does).
program yrec_extra_heat
      use run_driver, only: run_yrec
      use energy_sources, only: add_energy_source
      use extra_heat_source, only: read_extra_heat, uniform_heat
      implicit none
      integer :: ierr

      call read_extra_heat('extra_heat.nml', ierr)
      if (ierr /= 0) stop 1
      call add_energy_source('uniform_heat', uniform_heat, ierr)
      if (ierr /= 0) stop 1
      ierr = 0
      call run_yrec(ierr)
      if (ierr > 0) stop 1
end program yrec_extra_heat
