!----------------------------------------------------------------------
! thdiff
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original thdiff.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
!*************************************************************
! This routine was written by Anne A. Thoul, at the Institute
! for Advanced Study, Princeton, NJ 08540.
! See Thoul et al., Ap.J. 421, p. 828 (1994)
! The subroutines LUBKSB and LUDCMP are from Numerical Recipes.
!*************************************************************
! This routine inverses the burgers equations.
!
! The system contains N equations with N unknowns.
! The equations are: the M momentum equations,
!                    the M energy equations,
!                    two constraints: the current neutrality
!                                     the zero fluid velocity.
! The unknowns are: the M diffusion velocities,
!                   the M heat fluxes,
!                   the electric field E
!                   the gravitational force g.
!
!**************************************************
! The parameter num_species is the number of species considered.
!
! Fluid 1 is the hydrogen
! Fluid 2 is the helium
! Fluids 3 to num_species-1 are the heavy elements
! Fluid num_species is the electrons
!
! The vectors atomic_weight,charge and mass_fraction contain the atomic
! mass numbers, the charges (ionization), and the mass fractions, of
! the elements.
! NOTE: Since num_species is the electron fluid, its mass and charge
!      must be
!      atomic_weight(num_species)=m_e/m_u
!      charge(num_species)=-1.
!
! The array coulomb_log contains the values of the Coulomb Logarithms.
! The vector pressure_coeff, temp_coeff, and array conc_coeff contains
! the results for the diffusion coefficients.
!
! Internal notation (c, cc, ac, xx, y, yy, k, alpha, nu, gamma, delta,
! ga) is kept as in the original Thoul et al. (1994) Burgers-equation
! formulation rather than renamed -- these are terse algebraic
! quantities from the published derivation (concentrations, mass-
! weighted sums, and the right/left-hand-side arrays of the linear
! system solved below) and a confident, more descriptive rename is not
! attempted here; see the inline comments below for their roles.
subroutine thoul_diffusion(num_species, atomic_weight, charge, mass_fraction, &
     coulomb_log, pressure_coeff, temp_coeff, conc_coeff, ierr)

      use numerics_lib
      implicit none
      integer, parameter :: mmax=20, nmax=42

      integer, intent(in) :: num_species
      double precision, intent(in) :: atomic_weight(num_species), &
           charge(num_species)
      double precision, intent(inout) :: mass_fraction(num_species)
      double precision, intent(in) :: coulomb_log(num_species,num_species)
      double precision, intent(out) :: pressure_coeff(num_species), &
           temp_coeff(num_species)
      double precision, intent(out) :: conc_coeff(num_species,num_species)
      integer, intent(out) :: ierr

      integer :: indx(nmax)
      double precision :: c(mmax), xx(mmax,mmax), y(mmax,mmax), yy(mmax,mmax)
      double precision :: k(mmax,mmax), nu(nmax)
      double precision :: alpha(nmax), gamma(nmax,nmax), delta(nmax,nmax), &
           ga(nmax)
      double precision :: ko
      integer :: i, j, l, n
      double precision :: cc, ac, temp, d

! The vectors c and the concentration gradients contain the
! concentrations and the concentration gradients;
! cc is the total concentration: cc=sum(c_s)
! ac is proportional to the mass density: ac=sum(a_s c_s)
! The arrays xx,y,yy and k are various parameters which appear in
! Burgers equations.
! The vectors and arrays alpha, nu, gamma, delta, and ga represent
! the "right- and left-hand-sides" of Burgers equations, and later
! the diffusion coefficients.


! Initialize parameters:

      ierr = 0
      ko=2
      n=2*num_species+2
      do i=1,num_species
         c(i)=0.d0
      enddo
      cc=0.d0
      ac=0.d0

! Calculate concentrations from mass fractions:

      temp=0.d0
      do i=1,num_species-1
         temp=temp+charge(i)*mass_fraction(i)/atomic_weight(i)
      enddo
      do i=1,num_species-1
         c(i)=mass_fraction(i)/atomic_weight(i)/temp
      enddo
      c(num_species)=1.d0

! Calculate cc and ac:

      do i=1,num_species
         cc=cc+c(i)
         ac=ac+atomic_weight(i)*c(i)
      enddo

! Calculate the mass fraction of electrons:

      mass_fraction(num_species)=atomic_weight(num_species)/ac

! Calculate the coefficients of the burgers equations
      do i=1,num_species
         do j=1,num_species
            xx(i,j)=atomic_weight(j)/(atomic_weight(i)+atomic_weight(j))
            y(i,j)=atomic_weight(i)/(atomic_weight(i)+atomic_weight(j))
            yy(i,j)=3.0d0*y(i,j)+1.3d0*xx(i,j)*atomic_weight(j)/atomic_weight(i)
            k(i,j)=1.d0*coulomb_log(i,j)* &
                 sqrt(atomic_weight(i)*atomic_weight(j)/ &
                 (atomic_weight(i)+atomic_weight(j)))*c(i)*c(j)* &
                 charge(i)**2*charge(j)**2
         enddo
      enddo

! Write the burgers equations and the two constraints as
! alpha_s dp + nu_s dT + sum_t(not 2 or M) gamma_st dC_t
!                     = sum_t delta_st w_t

      do i=1,num_species
         alpha(i)=c(i)/cc
         nu(i)=0.d0
         do j=1,num_species
            gamma(i,j)=0.d0
         enddo
         do j=1,num_species
            if ((j.ne.2).and.(j.ne.num_species)) then
               gamma(i,j)=-c(j)/cc+c(2)/cc*charge(j)*c(j)/charge(2)/c(2)
               if (j.eq.i) then
                  gamma(i,j)=gamma(i,j)+1.d0
               endif
               if (i.eq.2) then
                  gamma(i,j)=gamma(i,j)-charge(j)*c(j)/charge(2)/c(2)
               endif
               gamma(i,j)=gamma(i,j)*c(i)/cc
            endif
         enddo

         do j=num_species+1,n
            gamma(i,j)=0.d0
         enddo
      enddo

      do i=num_species+1,n-2
         alpha(i)=0.d0
         nu(i)=2.5d0*c(i-num_species)/cc
         do j=1,n
            gamma(i,j)=0.d0
         enddo
      enddo

      alpha(n-1)=0.d0
      nu(n-1)=0
      do j=1,n
         gamma(n-1,j)=0.d0
      enddo

      alpha(n)=0.d0
      nu(n)=0
      do j=1,n
         gamma(n,j)=0.d0
      enddo

      do i=1,n
         do j=1,n
            delta(i,j)=0.d0
         enddo
      enddo

      do i=1,num_species
         do j=1,num_species
            if (j.eq.i) then
               do l=1,num_species
                  if(l.ne.i) then
                     delta(i,j)=delta(i,j)-k(i,l)
                  endif
               enddo
            else
               delta(i,j)=k(i,j)
            endif
         enddo

         do j=num_species+1,n-2
            if(j-num_species.eq.i) then
               do l=1,num_species
                  if (l.ne.i) then
                     delta(i,j)=delta(i,j)+0.6d0*xx(i,l)*k(i,l)
                  endif
               enddo
            else
               delta(i,j)=-0.6d0*y(i,j-num_species)*k(i,j-num_species)
            endif
         enddo

         delta(i,n-1)=c(i)*charge(i)

         delta(i,n)=-c(i)*atomic_weight(i)
      enddo

      do i=num_species+1,n-2
         do j=1,num_species
            if (j.eq.i-num_species) then
               do l=1,num_species
                  if (l.ne.i-num_species) then
                     delta(i,j)=delta(i,j)+1.5d0*xx(i-num_species,l)* &
                          k(i-num_species,l)
                  endif
               enddo
            else
               delta(i,j)=-1.5d0*xx(i-num_species,j)*k(i-num_species,j)
            endif
         enddo

         do j=num_species+1,n-2
            if (j-num_species.eq.i-num_species) then
               do l=1,num_species
                  if (l.ne.i-num_species) then
                     delta(i,j)=delta(i,j)-y(i-num_species,l)*k(i-num_species,l)* &
                          (1.6d0*xx(i-num_species,l)+yy(i-num_species,l))
                  endif
               enddo
               delta(i,j)=delta(i,j)-0.8d0*k(i-num_species,i-num_species)
            else
               delta(i,j)=2.7d0*k(i-num_species,j-num_species)* &
                    xx(i-num_species,j-num_species)*y(i-num_species,j-num_species)
            endif
         enddo

         delta(i,n-1)=0.d0

         delta(i,n)=0.d0
      enddo

      do j=1,num_species
         delta(n-1,j)=c(j)*charge(j)
      enddo
      do j=num_species+1,n
         delta(n-1,j)=0.d0
      enddo

      do j=1,num_species
         delta(n,j)=c(j)*atomic_weight(j)
      enddo
      do j=num_species+1,n
         delta(n,j)=0.d0
      enddo

! Inverse the system for each possible right-hand-side, i.e.,
! if alpha is the r.h.s., we obtain the coefficient A_p
! if nu    ---------------------------------------- A_T
! if gamma(i,j) ----------------------------------- A_Cj
!
! If I=1, we obtain the hydrogen diffusion velocity
! If I=2, ------------- helium   ------------------
! If I=3,M-1, --------- heavy element -------------
! If I=M, ------------- electrons -----------------
! For I=M,2M, we get the heat fluxes
! For I=N-1, we get the electric field
! For I=N, we get the gravitational force g

      call ludcmp(delta,n,nmax,indx,d,ierr)
      if (ierr /= 0) return

      call lubksb(delta,n,nmax,indx,alpha)
      call lubksb(delta,n,nmax,indx,nu)
      do j=1,n
         do i=1,n
            ga(i)=gamma(i,j)
         enddo
         call lubksb(delta,n,nmax,indx,ga)
         do i=1,n
            gamma(i,j)=ga(i)
         enddo
      enddo

! The results for the coefficients must be multiplied by p/K_0:

      do i=1,num_species
         alpha(i)=alpha(i)*ko*ac*cc
         nu(i)=nu(i)*ko*ac*cc
         do j=1,num_species
            gamma(i,j)=gamma(i,j)*ko*ac*cc
         enddo
      enddo

      do i=1,num_species
         pressure_coeff(i)=alpha(i)
         temp_coeff(i)=nu(i)
         do j=1,num_species
            conc_coeff(i,j)=gamma(i,j)
         enddo
      enddo


      return

end subroutine thoul_diffusion
