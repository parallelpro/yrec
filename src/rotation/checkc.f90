!----------------------------------------------------------------------
! checkc
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original checkc.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! SR CHECKC PERFORMS SEVERAL FUNCTIONS.
! FIRST, IT CHECKS FOR ANOMALOUS COMPOSITIONS.
! IF THEY ARE ENCOUNTERED, THE TIMESTEP IS CUT.
! IT ALSO COMPUTES THE CHANGE IN MEAN MOLECULAR WEIGHT DUE TO
! COMPOSITION DIFFUSION.
!
! INPUT VARIABLES:
! composition :  MASS FRACTION OF ALL THE SPECIES BEING TRACKED AS A
!    FUNCTION OF MASS.
! iteration_number : THE PROGRAM ITERATES FOR THE DIFFUSION
!    COEFFICIENTS; IT IS THE ITERATION NUMBER.
! max_iterations : USER PARAMETER - MAXIMUM NUMBER OF ITERATIONS.
! cut_count : NUMBER OF TIMES DIFFUSION TIMESTEP HAS BEEN CUT.
! converged_flag : FLAG SET T IF DIFFUSION COEFFICEINTS HAVE CONVERGED.
! print_flag : FLAG SET T IF OUTPUT ABOUT CHECMICAL DIFFUSION DESIRED.
! num_zones : NUMBER OF MODEL POINTS.
!
! OUTPUT VARIABLES:
! dt : DIFFUSION TIMESTEP, WHICH CAN BE CUT IF ERRORS IN THE
!    COMPOSITION DIFFUSION ARE DISCOVERED.
! mix_phys%amum (former common/mdphy/) : NEW RUN OF MEAN MOLECULAR WEIGHT.
! cut_count : NUMBER OF TIMES DIFFUSION TIMESTEP HAS BEEN CUT.
! converged_flag : SET F IF ERRORS IN COMPOSITION DIFFUSION DISCOVERED.
! redo_flag : SET T IF ERRORS IN COMPOSITION DIFFUSION DISCOVERED.
subroutine checkc(composition, iteration_number, print_flag, num_zones, &
     dt, cut_count, converged_flag, redo_flag)

      use const_lib
      use mdphy_lib
      use oldmod_lib
      use luout_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(inout) :: composition(15,json)
      integer, intent(in) :: iteration_number
      logical, intent(in) :: print_flag
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: dt
      integer, intent(inout) :: cut_count
      logical, intent(inout) :: converged_flag
      logical, intent(out) :: redo_flag

! common/comp/: envelope_hydrogen_fraction/amuenv are used here. Naming
! matches getopac.f90/hpoint.f90.
      double precision :: envelope_hydrogen_fraction, &
           envelope_metal_fraction, zenvm, amuenv, fxenv(12), xnew, &
           znew, stotal, senv
      common/comp/ envelope_hydrogen_fraction, envelope_metal_fraction, &
           zenvm, amuenv, fxenv, xnew, znew, stotal, senv

! common/comp2/: envelope_helium_fraction/envelope_he3_fraction
! (originally YENV/Y3ENV), both used here. Not referenced in any
! already-converted file.
      double precision :: envelope_helium_fraction, envelope_he3_fraction
      common/comp2/ envelope_helium_fraction, envelope_he3_fraction

! common/difus/: only max_iterations is used here. Naming matches
! getw.f90/dadcoeft.f90.
      double precision :: dtdif, convergence_tolerance
      integer :: itdif1, max_iterations
      common/difus/ dtdif, convergence_tolerance, itdif1, max_iterations





      double precision :: atomic_weight(4)
      data atomic_weight/1.007825d0,4.002603d0,12.0d0,3.01603d0/
      save

! locals
      integer :: num_diffused_species, species_index, zone_index
      double precision :: max_fractional_comp_change
      integer :: max_change_zone, max_change_species
      double precision :: min_comp_for_check, fractional_comp_change
      double precision :: delta_hydrogen, delta_helium, delta_metal, &
           delta_helium3
! amu_calc_temp is reused for two sequential 1/(mean weight) sums, as
! in the original (first for the ion mean weight, then for the
! electron mean weight).
      double precision :: amu_calc_temp, ion_mean_weight_inverse, &
           electron_mean_weight_inverse

!  CHECK FOR ANOMALOUS COMPOSITIONS.
!  PRIOR TO THE LAST ITERATION, ONLY DIFFUSION OF H,HE,HE3 PERFORMED.
!  FIND NUMBER OF SPECIES BEING DIFFUSED.
      if(iteration_number.eq.max_iterations)then
         num_diffused_species = 11
      else
         num_diffused_species = 4
      endif
      redo_flag = .false.
      do 20 species_index = 1,num_diffused_species
!  composition(3,...) IS Z, WHICH IS NOT DIFFUSED AS A UNIT.
         if(species_index.eq.3)goto 20
         do 10 zone_index = 2,num_zones-1
            if(composition(species_index,zone_index).lt.0.0d0.or. &
                 composition(species_index,zone_index).gt.1.0d0)then
!  SOME SPECIES CAN BE MIXED INTO REGIONS WHERE THEY ARE DESTROYED VERY
!  QUICKLY.  THE INTERPLAY BETWEEN DIFFUSION AND NUCLEAR BURNING CAN LEAD
!  TO SMALL NEGATIVE ABUNDANCES.  THEREFORE, ZERO OUT THESE SMALL NUMBERS
!  RATHER THAN STOPPING THE CODE.
               if(abs(composition(species_index,zone_index)).lt. &
                    1.0d-5*composition(species_index,num_zones))then
                  if(zone_index.eq.1.or.zone_index.eq.num_zones)then
                     composition(species_index,zone_index) = 0.0d0
                  else
                     composition(species_index,zone_index) = &
                          max(0.0d0,0.5d0* &
                          (composition(species_index,zone_index+1)+ &
                          composition(species_index,zone_index-1)))
                  endif
                  goto 10
               endif
               cut_count = cut_count + 1
               if(cut_count.gt.3)then
                  write(6,1010) species_index,zone_index, &
                       composition(species_index,zone_index)
                  write(short_file_unit,1010) species_index,zone_index, &
                       composition(species_index,zone_index)
 1010 format(1x,39('>'),39('>')/' ERROR IN SR CHECKC'/ &
              ' ANOMALOUS COMP NUMBER',i2,' IN ZONE',i5,' ABUNDANCE ', &
              1pe12.3/' 3 ATTEMPTS AT TIMESTEP CUTTING FAILED'/ &
              'RUN STOPPED')
                  stop
               else
                  redo_flag = .true.
                  converged_flag = .false.
                  dt = 0.5d0*dt
                  write(6,1015)cut_count,species_index,zone_index, &
                       composition(species_index,zone_index)
                  write(short_file_unit,1015)cut_count,species_index, &
                       zone_index,composition(species_index,zone_index)
 1015 format(' ERROR IN SR CHECKC'/' TIMESTEP CUT NUMBER ',i2, &
              ' DUE TO ANOMALOUS COMP NUMBER',i2,' IN ZONE',i5, &
              ' ABUNDANCE ',1pe12.3)
                  goto 100
               endif
            endif
   10    continue
   20 continue
      if(iteration_number.eq.max_iterations.and.print_flag) then
!  FIND MAXIMUM FRACTIONAL CHANGE IN COMPOSITION AND PRINT IT OUT.
         max_fractional_comp_change = 0.0d0
         max_change_zone = 0
         max_change_species = 0
         do 40 species_index = 1,num_diffused_species
            if(species_index.eq.3)goto 40
! min_comp_for_check IS USED TO GUARD AGAINST DIVISION BY ZERO.
            min_comp_for_check = max(1.0d-6* &
                 composition(species_index,num_zones),1.0d-20)
            do 30 zone_index = 1,num_zones
               if(prev_model%old_composition(species_index,zone_index).lt. &
                    min_comp_for_check)goto 30
               fractional_comp_change = &
                    (composition(species_index,zone_index)- &
                    prev_model%old_composition(species_index,zone_index))/ &
                    prev_model%old_composition(species_index,zone_index)
               if(abs(fractional_comp_change).gt. &
                    abs(max_fractional_comp_change)) then
                  max_fractional_comp_change = fractional_comp_change
                  max_change_zone = zone_index
                  max_change_species = species_index
               endif
   30       continue
   40    continue
         write(*,50)max_fractional_comp_change,max_change_species, &
              max_change_zone
   50 format(' MAX FRAC.COMP.CHANGE',1pe12.3,' SPECIES',i2, &
              ' AT PT.',i5)
         if(use_extended_composition)write(*,60) &
              composition(14,num_zones),prev_model%old_composition(14,num_zones)
   60 format(5x,'NEW SURFACE LI',1pe14.4,'OLD VALUE',e14.4)
      endif
!  FIND NEW RUN OF MEAN MOLECULAR WEIGHT ASSUMING FULLY IONIZED GAS.
!  AMUENV IS(1/MEAN MOLECULAR WEIGHT PER ION OF THE SURFACE MIXTURE.)
!  CORRECTION FOR PARTIAL IONIZATION NEEDED IN MASSIVE STARS.
      if(iteration_number.gt.1)then
         do 90 zone_index = 1,num_zones
            delta_hydrogen = composition(1,zone_index)- &
                 envelope_hydrogen_fraction
            delta_helium = composition(2,zone_index)- &
                 envelope_helium_fraction
            delta_metal = composition(3,zone_index)- &
                 envelope_metal_fraction
            delta_helium3 = composition(4,zone_index)- &
                 envelope_he3_fraction
            amu_calc_temp = amuenv + delta_hydrogen/atomic_weight(1) + &
                 delta_helium/atomic_weight(2) + &
                 delta_metal/atomic_weight(3) + &
                 delta_helium3/atomic_weight(4)
            ion_mean_weight_inverse = 1.0d0/amu_calc_temp
            amu_calc_temp = composition(1,zone_index)/atomic_weight(1)+ &
                 2.0d0*(composition(4,zone_index)/atomic_weight(4) &
                 + composition(2,zone_index)/atomic_weight(2)) + &
                 0.5d0*composition(3,zone_index)
            electron_mean_weight_inverse = 1.0d0/amu_calc_temp
            mix_phys%amum(zone_index) = ion_mean_weight_inverse* &
                 electron_mean_weight_inverse/ &
                 (ion_mean_weight_inverse+electron_mean_weight_inverse)
   90    continue
      endif
  100 continue

      return
end subroutine checkc
