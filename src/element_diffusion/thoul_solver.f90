!----------------------------------------------------------------------
! thoul_solver
!----------------------------------------------------------------------
! The diffusion coefficients of a fully ionized plasma from Burgers'
! equations, solved as Thoul, Bahcall & Loeb (1994, ApJ 421, 828) do (A.
! A. Thoul's routine, Institute for Advanced Study). For M species (the
! last the electrons), the 2M + 2 unknowns are the M diffusion
! velocities w_s, the M heat fluxes, the electric field and the
! gravitational force; the equations are the M momentum equations, the
! M energy equations, current neutrality (sum_s Z_s c_s w_s = 0) and no
! net mass flow (sum_s A_s c_s w_s = 0). Written as
!     alpha_s d ln P + nu_s d ln T + sum_t gamma_st d ln C_t
!         = sum_t delta_st w_t,
! the system is inverted once (ludcmp) and solved for each right-hand
! side (lubksb): alpha gives the coefficient A_P of each species'
! velocity, nu gives A_T and gamma(:,t) gives A_C for the gradient of
! species t (helium and the electrons excepted: their gradients follow
! from the others'). The velocities are in units of p/K_0, so the
! coefficients scale as 1/ln Lambda when every pair shares one Coulomb
! logarithm.
!
! The notation is the paper's: c the concentrations, cc their sum, ac
! sum_s A_s c_s (proportional to the density), xx, y, yy and k the
! collision terms between pairs.
module thoul_solver
      use, intrinsic :: iso_fortran_env, only: dp => real64
      implicit none
      private
      public :: thoul_diffusion

contains

! ---------------------------------------------------------------
! The coefficients pressure_coeff (A_P), temp_coeff (A_T) and
! conc_coeff (A_C, conc_coeff(s,t) for species t's gradient) of the M =
! size(atomic_weight) species: hydrogen first, helium second, the heavy
! elements, and the electrons last (atomic weight m_e/m_u, charge -1),
! from their atomic weights (amu), charges (e), mass fractions (the
! electrons' is set here) and pairwise Coulomb logarithms. ierr: ludcmp's
! status for a singular system (the caller reports it).
pure subroutine thoul_diffusion(atomic_weight, charge, mass_fraction, coulomb_log, &
     pressure_coeff, temp_coeff, conc_coeff, ierr)
      use numerics_lib, only: lubksb, ludcmp
! the most species, and unknowns, the arrays hold
      integer, parameter :: mmax = 20, nmax = 42
      real(dp), intent(in) :: atomic_weight(:), charge(:)
      real(dp), intent(inout) :: mass_fraction(:)
      real(dp), intent(in) :: coulomb_log(:,:)
      real(dp), intent(out) :: pressure_coeff(:), temp_coeff(:), conc_coeff(:,:)
      integer, intent(out) :: ierr
      integer :: indx(nmax)
      real(dp) :: c(mmax), xx(mmax,mmax), y(mmax,mmax), yy(mmax,mmax), k(mmax,mmax)
      real(dp) :: alpha(nmax), nu(nmax), gamma(nmax,nmax), delta(nmax,nmax), ga(nmax)
      real(dp) :: ko, cc, ac, temp, d
      integer :: num_species, i, j, l, n

      ierr = 0
      num_species = size(atomic_weight)
      ko = 2.0d0
      n = 2*num_species + 2
      c(1:num_species) = 0.0d0
      cc = 0.0d0
      ac = 0.0d0

! The concentrations from the mass fractions, cc, ac, and the electrons'
! mass fraction.
      temp = 0.0d0
      do i = 1, num_species-1
         temp = temp + charge(i)*mass_fraction(i)/atomic_weight(i)
      end do
      do i = 1, num_species-1
         c(i) = mass_fraction(i)/atomic_weight(i)/temp
      end do
      c(num_species) = 1.0d0
      do i = 1, num_species
         cc = cc + c(i)
         ac = ac + atomic_weight(i)*c(i)
      end do
      mass_fraction(num_species) = atomic_weight(num_species)/ac

! The collision terms of Burgers' equations.
      do i = 1, num_species
         do j = 1, num_species
            xx(i,j) = atomic_weight(j)/(atomic_weight(i) + atomic_weight(j))
            y(i,j) = atomic_weight(i)/(atomic_weight(i) + atomic_weight(j))
            yy(i,j) = 3.0d0*y(i,j) + 1.3d0*xx(i,j)*atomic_weight(j)/atomic_weight(i)
            k(i,j) = coulomb_log(i,j)* &
                 sqrt(atomic_weight(i)*atomic_weight(j)/ &
                 (atomic_weight(i) + atomic_weight(j)))*c(i)*c(j)* &
                 charge(i)**2*charge(j)**2
         end do
      end do

! The right-hand sides: alpha_s, nu_s and gamma_st (no column for helium
! or the electrons).
      do i = 1, num_species
         alpha(i) = c(i)/cc
         nu(i) = 0.0d0
         gamma(i,1:num_species) = 0.0d0
         do j = 1, num_species
            if (j /= 2 .and. j /= num_species) then
               gamma(i,j) = -c(j)/cc + c(2)/cc*charge(j)*c(j)/charge(2)/c(2)
               if (j == i) gamma(i,j) = gamma(i,j) + 1.0d0
               if (i == 2) gamma(i,j) = gamma(i,j) - charge(j)*c(j)/charge(2)/c(2)
               gamma(i,j) = gamma(i,j)*c(i)/cc
            end if
         end do
         gamma(i,num_species+1:n) = 0.0d0
      end do
      do i = num_species+1, n-2
         alpha(i) = 0.0d0
         nu(i) = 2.5d0*c(i-num_species)/cc
         gamma(i,1:n) = 0.0d0
      end do
      alpha(n-1) = 0.0d0
      nu(n-1) = 0.0d0
      gamma(n-1,1:n) = 0.0d0
      alpha(n) = 0.0d0
      nu(n) = 0.0d0
      gamma(n,1:n) = 0.0d0

! The left-hand side delta_st: the momentum equations' rows...
      delta(1:n,1:n) = 0.0d0
      do i = 1, num_species
         do j = 1, num_species
            if (j == i) then
               do l = 1, num_species
                  if (l /= i) delta(i,j) = delta(i,j) - k(i,l)
               end do
            else
               delta(i,j) = k(i,j)
            end if
         end do
         do j = num_species+1, n-2
            if (j-num_species == i) then
               do l = 1, num_species
                  if (l /= i) delta(i,j) = delta(i,j) + 0.6d0*xx(i,l)*k(i,l)
               end do
            else
               delta(i,j) = -0.6d0*y(i,j-num_species)*k(i,j-num_species)
            end if
         end do
         delta(i,n-1) = c(i)*charge(i)
         delta(i,n) = -c(i)*atomic_weight(i)
      end do
! ... the energy equations' rows ...
      do i = num_species+1, n-2
         do j = 1, num_species
            if (j == i-num_species) then
               do l = 1, num_species
                  if (l /= i-num_species) delta(i,j) = delta(i,j) + 1.5d0*xx(i-num_species,l)* &
                       k(i-num_species,l)
               end do
            else
               delta(i,j) = -1.5d0*xx(i-num_species,j)*k(i-num_species,j)
            end if
         end do
         do j = num_species+1, n-2
            if (j-num_species == i-num_species) then
               do l = 1, num_species
                  if (l /= i-num_species) delta(i,j) = delta(i,j) - y(i-num_species,l)* &
                       k(i-num_species,l)*(1.6d0*xx(i-num_species,l) + yy(i-num_species,l))
               end do
               delta(i,j) = delta(i,j) - 0.8d0*k(i-num_species,i-num_species)
            else
               delta(i,j) = 2.7d0*k(i-num_species,j-num_species)* &
                    xx(i-num_species,j-num_species)*y(i-num_species,j-num_species)
            end if
         end do
         delta(i,n-1) = 0.0d0
         delta(i,n) = 0.0d0
      end do
! ... current neutrality and no net mass flow.
      do j = 1, num_species
         delta(n-1,j) = c(j)*charge(j)
      end do
      delta(n-1,num_species+1:n) = 0.0d0
      do j = 1, num_species
         delta(n,j) = c(j)*atomic_weight(j)
      end do
      delta(n,num_species+1:n) = 0.0d0

! The solve, for each right-hand side: alpha gives A_P, nu A_T, and
! gamma(:,j) A_C for the gradient of species j. Rows 1..M are the
! velocities, M+1..2M the heat fluxes, 2M+1 the electric field, 2M+2
! the gravitational force.
      call ludcmp(delta, n, nmax, indx, d, ierr)
      if (ierr /= 0) return
      call lubksb(delta, n, nmax, indx, alpha)
      call lubksb(delta, n, nmax, indx, nu)
      do j = 1, n
         ga(1:n) = gamma(1:n,j)
         call lubksb(delta, n, nmax, indx, ga)
         gamma(1:n,j) = ga(1:n)
      end do

! The velocities in units of p/K_0.
      do i = 1, num_species
         alpha(i) = alpha(i)*ko*ac*cc
         nu(i) = nu(i)*ko*ac*cc
         do j = 1, num_species
            gamma(i,j) = gamma(i,j)*ko*ac*cc
         end do
      end do
      do i = 1, num_species
         pressure_coeff(i) = alpha(i)
         temp_coeff(i) = nu(i)
         do j = 1, num_species
            conc_coeff(i,j) = gamma(i,j)
         end do
      end do
end subroutine thoul_diffusion

end module thoul_solver

! Sources
! - Thoul, A. A., Bahcall, J. N. & Loeb, A. 1994, ApJ 421, 828: the
!   routine is A. A. Thoul's; ludcmp and lubksb are Numerical Recipes'
!   (numerics_lib).
