!----------------------------------------------------------------------
! microdiff_cod
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original microdiff_cod.f; only variable names, source form, and
! comment style were updated. Validated against the Stage 0 regression
! suite (examples/run_standard_solar_model).
!
! CRITICAL: several numeric literals below (0.0, 1.0, 0.5, 2.2) are
! written WITHOUT a D0/D-3 exponent suffix in the original, meaning
! Fortran parses them as single precision before widening to double --
! numerically different from a D0-suffixed literal. These are copied
! character-for-character from microdiff_cod.f; do not "fix" them.
!
! G Somers 5/15: THIS IS A SUBROUTINE WRITEN TO CALCULATE THE DIFFUSION
! COEFFICIENT ON THE EQUALLY-SPACED GRID USED FOR LIGHT ELEMENT DIFFUSION.
! THIS CODE TAKES IN THE NECESSARY PARAMETERS ON THE EQUAL GRID (NAMELY
! HCOMP(1,I), HCOMP(2,I), HCOMP(13-15,I), DENSITY, TEMPERATURE, HQPR,
! AND DELR(2,I)
!
! Part of the microdiff.f90 pipeline (see also microdiff_setup.f90,
! microdiff_mte.f90, microdiff_run.f90, microdiff_etm.f90); calls
! thdiff.f90 (the Thoul et al. 1994 Burgers-equation solver) once per
! equally spaced grid point.
subroutine microdiff_cod(num_eq_points, species_fraction, eq_radius, &
     eq_density, eq_temperature, eq_dlnp_dr, eq_del_grad, diffusion_coeff1, &
     diffusion_coeff2, hydrogen_dlnc_dr, atomic_weight_diffused, &
     atomic_charge_diffused, species_col)

      use scrtch_lib
      use const_lib
      implicit none
      integer, parameter :: json = 5000

      integer, intent(in) :: num_eq_points
      double precision, intent(in) :: species_fraction(3,json), &
           eq_radius(json), eq_density(json), eq_temperature(json), &
           eq_dlnp_dr(json), eq_del_grad(json)
      double precision, intent(out) :: diffusion_coeff1(json), &
           diffusion_coeff2(json)
      double precision, intent(inout) :: hydrogen_dlnc_dr(json)
      double precision, intent(in) :: atomic_weight_diffused, &
           atomic_charge_diffused
      integer, intent(in) :: species_col




! common/theage/: not used in this file. Naming matches mix.f90.
      double precision :: dage
      common/theage/ dage

      double precision :: atomic_weight(4), atomic_charge(4), mass_frac(4), &
           coulomb_log(4,4), pressure_coeff(4), temp_coeff(4), &
           conc_coeff(4,4), &
           coeff_scale(json), pressure_term(json), temp_term(json), &
           hydrogen_term(json), diffusion_term(json), &
           hydrogen_concen(json)
      double precision :: zxa, ac, ni, cz, xij, ne, ao, lambda_debye, &
           lambda, concen(4)
      double precision :: ln_lambda
      integer :: num_species
      data num_species/4/
!       DATA FGRLI/1.0,1.0,1.0,1.0/
      save

      integer :: i, ii, jj
      double precision :: bl_radius_scale_local, bl_temp_scale_local, &
           hru_i, htu_i, fac, ap, at, ah, ad, dlncdr, coni, conip1, conim1, &
           dradi, t1, t2, rho, t

! SET UP THE ATOMIC WEIGHT AND CHARGE MATRICIES
      atomic_weight(1) = 1.008d0
      atomic_weight(2) = 4.004d0
      atomic_weight(3) = atomic_weight_diffused
      atomic_weight(4) = 5.486d-4
      atomic_charge(1) = 1.0d0
      atomic_charge(2) = 2.0d0
      atomic_charge(3) = atomic_charge_diffused
      atomic_charge(4) = -1.0d0
! SET LN_LAMBDA, CON_RAD, AND CON_TEMP.
      ln_lambda = 2.2
      bl_radius_scale_local=1.0d0/6.9598d10
      bl_temp_scale_local=1.0d-7
! CALCULATE DIFFUSION COEFFICIENTS FOR EACH LAYER.
      do 5 i = 1,num_eq_points
         mass_frac(1) = species_fraction(1,i)
         mass_frac(2) = species_fraction(2,i)
         mass_frac(3) = species_fraction(3,i)
! make sure XFRAC = 0.0 isn't used for diff coefficients
         if (mass_frac(1).lt.1.0d-24) mass_frac(1) = 1.0d-24
         if (mass_frac(2).lt.1.0d-24) mass_frac(2) = 1.0d-24
         if (mass_frac(3).lt.1.0d-24) mass_frac(3) = 1.0d-24
!        calculate concentrations from mass fractions:
         zxa=0.d0
         do ii=1,num_species-1
             zxa=zxa+atomic_charge(ii)*mass_frac(ii)/atomic_weight(ii)
         enddo
         do ii=1,num_species-1
            concen(ii)=mass_frac(ii)/(atomic_weight(ii)*zxa)
         enddo
         concen(num_species)=1.d0
!        save the hydrogen concentration when X is diffused.
         if(species_col.eq.1) hydrogen_concen(i) = concen(1)
!        now check whether the Thoul routine must be run. if not,
!        write COD1 = COD2 = 0. If its the first shell in the depleted
!        zone, permit the calculations so that AD is correct.
         if(species_fraction(species_col,i).eq.0.0.and.i.ne.num_eq_points)then
            if(species_fraction(species_col,i+1).eq.0.0)then
               diffusion_term(i) = 0.0
               goto 5
            endif
         endif
!        set relevant physical variables.
         rho = eq_density(i)
         t = eq_temperature(i)
!        calculate density of electrons (NE) from mass density (RHO):
         ac=0.d0
         do ii=1,num_species
          ac=ac+atomic_weight(ii)*concen(ii)
         enddo
         ne=rho/(1.6726d-24*ac)
!        calculate interionic distance (AO):
         ni=0.d0
         do ii=1,num_species-1
            ni=ni+concen(ii)*ne
         enddo
         ao=(0.23873d0/ni)**cc13
!        calculate Debye length (LAMBDAD):
         cz=0.d0
         do ii=1,num_species
          cz=cz+concen(ii)*atomic_charge(ii)**2
         enddo
         lambda_debye=6.9010d0*sqrt(t/(ne*cz))
!        calculate LAMBDA to use in Coulomb logarithm:
         lambda=max(lambda_debye,ao)
!        calculate Coulomb logarithms:
         do ii=1,num_species
            do jj=1,num_species
               xij=2.3939d3*t*lambda/abs(atomic_charge(ii)*atomic_charge(jj))
               coulomb_log(ii,jj)=0.81245d0 &
               *log(1.d0+0.18769d0*xij**1.2d0)
          enddo
         enddo
!
!        calculate the diffusion coefficients
!
         call thdiff(num_species,atomic_weight,atomic_charge,mass_frac, &
              coulomb_log,pressure_coeff,temp_coeff,conc_coeff)
!
         hru_i = eq_radius(i)
         htu_i = t*bl_temp_scale_local
!         FAC=FGRLI(KK)*HRU_I**2*HTU_I**2.5D0/LN_LAMBDA
!        JvS 01/26 Added support for FGRY and FGRZ modifications
!        to diffusion coefficients.
!         FAC=HRU_I**2*HTU_I**2.5D0/LN_LAMBDA
         fac=hru_i**2*htu_i**2.5d0/ln_lambda
         if(species_col.eq.1)then
            fac=fgry*hru_i**2*htu_i**2.5d0/ln_lambda
         endif
         if(species_col.eq.3)then
            fac=fgrz*hru_i**2*htu_i**2.5d0/ln_lambda
         endif
!        collect the first diffusion terms for hydroden.
!        collect the third diffusion terms for everything else.
         ap = -pressure_coeff(species_col)
         at = -temp_coeff(species_col)*eq_del_grad(i)
         ah = -conc_coeff(species_col,1)
         ad = -conc_coeff(species_col,species_col)
!        store the numbers so the hydrogen gradient can finish
!        being calculated; then use them later.
         coeff_scale(i) = fac
         pressure_term(i) = ap
         temp_term(i) = at
         hydrogen_term(i) = ah
         diffusion_term(i) = ad
    5 continue
      do 10 i = 1,num_eq_points
         fac = coeff_scale(i)
         ap = pressure_term(i)
         at = temp_term(i)
         ah = hydrogen_term(i)
         ad = diffusion_term(i)
!
!        The general diffusion velocity is written as:
!        Ap(D)*dlnP/dR + At(D)*dlnT/dR + Ax(D)*dlnX/dR + Ad(D)*dlnD/dR
!        The first two terms are pressure and temperature gradients, and
!        the last two terms are the hydrogen gradient and the diffused element
!        gradient. For hydrogren diffusion, we can neglect the third term, as
!        it is the same as the final term. For metal diffusion, we need all
!        four. (see Thoul et al. 1994 for details)
!
!        calculate the hydrogen gradient if diffusing X, then set
!        the local gradient to 0 to remove 2nd self-diff term.
!        look it up otherwise.
!
!        only calculate the gradient on the hydrogen call.
!        set the gradient at first and last points to 0.
!
         if(species_col.eq.1)then
            dlncdr = 1.0
            if (i.eq.1) then
               dlncdr = 0.0
            elseif (i.eq.num_eq_points) then
               dlncdr = 0.0
            else
               coni = log(hydrogen_concen(i))
               conip1 = log(hydrogen_concen(i+1))
               conim1 = log(hydrogen_concen(i-1))
               dradi = eq_radius(i+1)-eq_radius(i)
               t1 = (coni-conim1)/dradi
               t2 = (conip1-coni)/dradi
               dlncdr = 0.5*(t1+t2)
            endif
            hydrogen_dlnc_dr(i) = dlncdr
!           zero out gradient for hydrogen abundance, as stated above.
            dlncdr = 0.0
         else
!           if not the hydrogen call, retrieve the gradient
            dlncdr = hydrogen_dlnc_dr(i)
         endif
         diffusion_coeff1(i) = fac*(eq_dlnp_dr(i)*(ap+at)+dlncdr*ah)* &
              species_fraction(species_col,i)
         diffusion_coeff2(i) = fac*ad
   10 continue
      return
end subroutine microdiff_cod
