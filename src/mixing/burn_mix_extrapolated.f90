!----------------------------------------------------------------------
! bsrotmix
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original bsrotmix.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! MHP 6/00 SUBROUTINE.  MIXING AND BURNING ARE PERFORMED SEQUENTIALLY
! FOR AN INCREASING NUMBER OF SUBSTEPS FOR A GIVEN TOTAL STEP.  THIS
! ROUTINE EXTRAPOLATES TO ZERO TIMESTEP (Bulirsch-Stoer / Richardson
! extrapolation in step size, using Neville's algorithm), given the
! composition (composition) computed with extrapolation_order
! (originally IEST) substeps out of the sequence substep_counts.
subroutine burn_mix_extrapolated(timestep, composition, extrapolation_order, num_zones, &
     species_begin, species_end, substep_counts, converged)
      use rotation_scratch_lib
      use star_info_lib, only: star, json
      implicit none

      double precision, intent(in) :: timestep
      double precision, intent(inout) :: composition(15,json)
      integer, intent(in) :: extrapolation_order, num_zones, &
           species_begin, species_end
      integer, intent(in) :: substep_counts(11)
      logical, intent(inout) :: converged


      double precision :: step_size_squared(11), current_value(15,json)
      integer :: active_species_id(15)
      logical :: species_active(15), use_cno_ratio_method(json), &
           he3_extrapolate_log(json), use_cno_ratio_species
      integer :: i, j, j1, species_count, num_active_species, k1, &
           max_error_zone, max_error_species
      double precision :: c12_abundance, c13_abundance, n14_abundance, &
           o16_abundance, cno_sum_check, current_step_size_squared, &
           delta, extrap_weight1, extrap_weight2, prev_estimate, &
           max_relative_error

!      SAVE X,D,JJ,LDO,NMAX,LCNO,LCNCHECK
! DETERMINE WHICH SPECIES REQUIRE CALCULATION
      if (extrapolation_order.eq.1) then
         do j = species_begin,species_end
            species_active(j) = .false.
            do i = num_zones,1,-1
               if (composition(j,i).gt.1.0d-14) then
                  species_active(j) = .true.
                  exit
               end if
            end do
         end do
         if (species_active(4)) then
            do i = 1,num_zones
               if (composition(4,i).gt.1.0d-24) then
                  he3_extrapolate_log(i) = .true.
               else
                  he3_extrapolate_log(i) = .false.
               end if
            end do
         end if
! CHECK FOR ALTERNATE METHOD FOR DOING CNO ABUNDANCES; SOLVE FOR
! C13/C12, N14/C12, O16/N14, SUM OF CNO NUCLEI RATHER THAN FOR
! THE INDIVIDUAL ABUNDANCES OF C12,C13,N14,O16.
         if (species_active(5).and.species_active(6).and. &
              species_active(7).and.species_active(9)) then
            use_cno_ratio_species = .true.
            do i = 1,num_zones
               if (composition(5,i).gt.1.0d-24 .and. &
                    composition(6,i).gt.1.0d-24 &
                    .and. composition(7,i).gt.1.0d-24 .and. &
                    composition(9,i).gt.1.0d-24) then
                  use_cno_ratio_method(i) = .true.
               else
                  use_cno_ratio_method(i) = .false.
               end if
            end do
         else
            use_cno_ratio_species = .false.
            do i = 1,num_zones
               use_cno_ratio_method(i) = .false.
            end do
         end if
         species_count = 1
         do i = species_begin,species_end
            if (species_active(i)) then
               active_species_id(species_count) = i
               species_count = species_count+1
            end if
         end do
         num_active_species = species_count - 1
         if (num_active_species.lt.1) stop 999
      end if
! STORE CURRENT RESULTS IN VECTOR OF COMPOSITION (HCOMPA) AND
! ERROR (DCOMPA)
      do j1 = 1,num_active_species
         j = active_species_id(j1)
         do i = 1,num_zones
            rot_scr%bs_extrapolated_composition(j,i) = composition(j,i)
            rot_scr%bs_extrapolation_increment(j,i) = composition(j,i)
         end do
      end do
! EXTRAPOLATE IN THE LOG OF HE3 RATHER THAN ABSOLUTE ABUNDANCE
      if (species_active(4)) then
         do i = 1,num_zones
            if (he3_extrapolate_log(i)) then
               rot_scr%bs_extrapolated_composition(4,i) = log(composition(4,i))
               rot_scr%bs_extrapolation_increment(4,i) = &
                    rot_scr%bs_extrapolated_composition(4,i)
            end if
         end do
      end if
      if (use_cno_ratio_species) then
         do i = 1,num_zones
            if (use_cno_ratio_method(i)) then
               c12_abundance = composition(5,i)
               c13_abundance = composition(6,i)
               n14_abundance = composition(7,i)
               o16_abundance = composition(9,i)
! VECTOR 5: C13/C12
               rot_scr%bs_extrapolated_composition(5,i) = c13_abundance/c12_abundance
               rot_scr%bs_extrapolation_increment(5,i) = &
                    rot_scr%bs_extrapolated_composition(5,i)
! VECTOR 6: N14/C12
               rot_scr%bs_extrapolated_composition(6,i) = n14_abundance/c12_abundance
               rot_scr%bs_extrapolation_increment(6,i) = &
                    rot_scr%bs_extrapolated_composition(6,i)
! VECTOR 7: O16/N14
               rot_scr%bs_extrapolated_composition(7,i) = o16_abundance/n14_abundance
               rot_scr%bs_extrapolation_increment(7,i) = &
                    rot_scr%bs_extrapolated_composition(7,i)
! VECTOR 9: C12/12 + C13/13 + N14/14 + O16/16
               rot_scr%bs_extrapolated_composition(9,i) = c12_abundance/12.0d0 + &
                    c13_abundance/13.0d0 + n14_abundance/14.0d0 &
                    + o16_abundance/16.0d0
               rot_scr%bs_extrapolation_increment(9,i) = &
                    rot_scr%bs_extrapolated_composition(9,i)
            end if
         end do
      end if
! SET CURRENT RESULTS AS THE INITIAL EXTRAPOLATION FOR THE FIRST SET OF
! TIME STEPS AND EXIT.
      step_size_squared(extrapolation_order) = &
           (timestep/dfloat(substep_counts(extrapolation_order)))**2
      if (extrapolation_order.eq.1) then
         do j1 = 1,num_active_species
            j = active_species_id(j1)
            do i = 1,num_zones
               rot_scr%bs_extrapolation_table(extrapolation_order,j,i) = &
                    composition(j,i)
            end do
         end do
         if (species_active(4)) then
            do i = 1,num_zones
               if (he3_extrapolate_log(i)) then
                  rot_scr%bs_extrapolation_table(extrapolation_order,4,i) = &
                       log(composition(4,i))
               end if
            end do
         end if
         if (use_cno_ratio_species) then
            do i = 1,num_zones
               if (use_cno_ratio_method(i)) then
                  rot_scr%bs_extrapolation_table(extrapolation_order,5,i) = &
                       rot_scr%bs_extrapolated_composition(5,i)
                  rot_scr%bs_extrapolation_table(extrapolation_order,6,i) = &
                       rot_scr%bs_extrapolated_composition(6,i)
                  rot_scr%bs_extrapolation_table(extrapolation_order,7,i) = &
                       rot_scr%bs_extrapolated_composition(7,i)
                  rot_scr%bs_extrapolation_table(extrapolation_order,9,i) = &
                       rot_scr%bs_extrapolated_composition(9,i)
               end if
             end do
          end if
        converged = .false.
      else
! POLYNOMIAL FUNCTION EXTRAPOLATON
! CURRENT STEPSIZE
         current_step_size_squared = step_size_squared(extrapolation_order)
         do j1 = 1,num_active_species
            j = active_species_id(j1)
            do i = 1,num_zones
               current_value(j,i) = composition(j,i)
            end do
         end do
         if (species_active(4)) then
            do i = 1,num_zones
               if (he3_extrapolate_log(i)) then
                  current_value(4,i) = log(composition(4,i))
               end if
            end do
         end if
         if (use_cno_ratio_species) then
            do i = 1,num_zones
               if (use_cno_ratio_method(i)) then
                  c12_abundance = composition(5,i)
                  c13_abundance = composition(6,i)
                  n14_abundance = composition(7,i)
                  o16_abundance = composition(9,i)
! VECTOR 5: C13/C12
                  current_value(5,i) = c13_abundance/c12_abundance
! VECTOR 6: N14/C12
                  current_value(6,i) = n14_abundance/c12_abundance
! VECTOR 7: O16/N14
                  current_value(7,i) = o16_abundance/n14_abundance
! VECTOR 9: C12/12 + C13/13 + N14/14 + O16/16
                  current_value(9,i) = c12_abundance/12.0d0 + &
                       c13_abundance/13.0d0 + n14_abundance/14.0d0 &
                       + o16_abundance/16.0d0
               end if
            end do
         end if
         do k1 = 1,extrapolation_order-1
            delta = 1.0d0/(step_size_squared(extrapolation_order-k1)- &
                 current_step_size_squared)
            extrap_weight1 = current_step_size_squared*delta
            extrap_weight2 = step_size_squared(extrapolation_order-k1)*delta
            do j1 = 1,num_active_species
               j = active_species_id(j1)
               do i = 1,num_zones
                  prev_estimate = rot_scr%bs_extrapolation_table(k1,j,i)
                  rot_scr%bs_extrapolation_table(k1,j,i) = &
                       rot_scr%bs_extrapolation_increment(j,i)
                  delta = current_value(j,i) - prev_estimate
                  rot_scr%bs_extrapolation_increment(j,i) = extrap_weight1*delta
                  current_value(j,i) = extrap_weight2*delta
                  rot_scr%bs_extrapolated_composition(j,i) = &
                       rot_scr%bs_extrapolated_composition(j,i)+ &
                       rot_scr%bs_extrapolation_increment(j,i)
               end do
            end do
         end do
         do j1 = 1,num_active_species
            j = active_species_id(j1)
            do i = 1,num_zones
               rot_scr%bs_extrapolation_table(extrapolation_order,j,i) = &
                    rot_scr%bs_extrapolation_increment(j,i)
            end do
         end do
! NOW CHECK IF ERROR IS WITHIN TOLERANCE
         max_relative_error = 1.0d-30
         max_error_zone = 1
         max_error_species = 1
         do j1 = 1,num_active_species
            j = active_species_id(j1)
            do i = 1,num_zones
               if (rot_scr%bs_extrapolated_composition(j,i).gt.1.0d-12) then
                  prev_estimate = rot_scr%bs_extrapolation_increment(j,i)/ &
                       rot_scr%bs_extrapolated_composition(j,i)
                  if (abs(prev_estimate).gt.max_relative_error) then
                     max_error_zone = i
                     max_error_species = j
                     max_relative_error=abs(prev_estimate)
                  end if
               end if
            end do
         end do
!         WRITE(*,*)HCOMP(JMAX,IMAX),HCOMPA(JMAX,IMAX)
         write(*,10) extrapolation_order,max_relative_error, &
              max_error_species,max_error_zone, &
              composition(max_error_species,max_error_zone), &
              rot_scr%bs_extrapolated_composition(max_error_species,max_error_zone)
         if (max_relative_error.lt.1.0d-3) converged=.true.
         if (converged) then
 10      format(5x,'ITER',i3,' MAX ERR ',1pe10.2,' SPECIES ',i3, &
              ' SHELL ',i5,' ABUND',2e12.4)
! NOW CONVERT BACK TO ABSOLUTE ABUNDANCES.
            if (species_active(4)) then
               do i =1,num_zones
                  if (he3_extrapolate_log(i)) then
                     rot_scr%bs_extrapolated_composition(4,i)= &
                          exp(rot_scr%bs_extrapolated_composition(4,i))
                  end if
               end do
            end if
            if (use_cno_ratio_species) then
               do i = 1,num_zones
                  if (use_cno_ratio_method(i)) then
! VECTOR NINE IS C12/12+ C13/13 + N14/14 + O16/16.
! THEREFORE, C12 = SUMCNO/(1/12 + C13/C12*1/13 + N14/C12*1/14 +
! O16/N14*N14/C12*1/16)
                     cno_sum_check = (1.0d0/12.0d0) + &
                          (rot_scr%bs_extrapolated_composition(5,i)/13.0d0) + &
                          (rot_scr%bs_extrapolated_composition(6,i)/14.0d0) &
                          + (rot_scr%bs_extrapolated_composition(6,i)* &
                          rot_scr%bs_extrapolated_composition(7,i)/16.0d0)
                     c12_abundance = rot_scr%bs_extrapolated_composition(9,i)/ &
                          cno_sum_check
                     c13_abundance = rot_scr%bs_extrapolated_composition(5,i)* &
                          c12_abundance
                     n14_abundance = rot_scr%bs_extrapolated_composition(6,i)* &
                          c12_abundance
                     o16_abundance = rot_scr%bs_extrapolated_composition(7,i)* &
                          n14_abundance
!       WRITE(*,911)X12,X13,X14,X16,HCOMPA(5,I),HCOMPA(6,I),
!     *             HCOMPA(7,I),HCOMPA(9,I),SUM
! 911   FORMAT(1P9E10.2)
                     rot_scr%bs_extrapolated_composition(5,i) = c12_abundance
                     rot_scr%bs_extrapolated_composition(6,i) = c13_abundance
                     rot_scr%bs_extrapolated_composition(7,i) = n14_abundance
                     rot_scr%bs_extrapolated_composition(9,i) = o16_abundance
                  end if
               end do
            end if
            do j1 = 1,num_active_species
               j = active_species_id(j1)
               do i = 1,num_zones
                  composition(j,i) = rot_scr%bs_extrapolated_composition(j,i)
                  if (composition(j,i).lt.1.0d-24) composition(j,i)=0.0d0
               end do
            end do
         end if
      end if
      return
end subroutine burn_mix_extrapolated
