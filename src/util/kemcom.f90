!----------------------------------------------------------------------
! kemcom
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original kemcom.f; only variable names, source form, and comment
! style were updated.
!
! kemcom uses the reaction rates computed in engeb to implicitly solve
! for the abundances of H, He3, He4, C12, C13, N14, O16, and O18
! simultaneously. Don VandenBerg's notes explain what the sr is doing.
!
! Input variables:
!
! timestep_years - model time step in years.
! composition - array containing run of mass fractions of different
!  species. composition(1,..)=H; 2=He4; 3=Z; 4=He3; 5=C12; 6=C13;
!  7=N14; 9=O16; 11=O18. Elements 8(N15) and 10(O17) not currently
!  used; 12-15 are light elements whose burning is treated elsewhere.
! old_composition (common/oldmod/) - the array of abundances at the
!  start of the timestep.
! shell_mass - run of mass contained in each shell.
! log_temperature - run of model temperature.
! zone_begin, zone_end - starting and ending shells; different for a
!  convection zone.
! rate_pp..rate_zero13, frac_c12_alpha - reaction rates (excluding
!  terms that depend only on the composition) per gigayear per amu.
!  needs to be multiplied by the mass fractions of the reactants and
!  the atomic weight (in amu) of the product to get the rate of
!  change of the mass fraction.
!
! Numerical parameters in common block /burtol/:
! Abundances below min_abundance are zeroed out.
! absolute_tolerance and relative_tolerance are the absolute and
! relative tolerances for convergence.
! max_burn_iterations is the maximum number of iterations before the
! program will halt.
!
! Output variables: the new run of abundance composition.
subroutine kemcom(log_temperature, zone_begin, zone_end, rate_pp, &
     rate_he3_he3, rate_he3_he4, rate_c12_p, rate_c13_p, rate_n14_p, &
     rate_o16_p, rate_c13_alpha, rate_c12_alpha, rate_n14_alpha, &
     rate_triple_alpha, frac_c12_alpha, shell_mass, composition, &
     timestep_years)

      use const_lib
      use luout_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: log_temperature(json)
      integer, intent(in) :: zone_begin, zone_end
      double precision, intent(in) :: rate_pp(json), rate_he3_he3(json), &
           rate_he3_he4(json), rate_c12_p(json), rate_c13_p(json), &
           rate_n14_p(json), rate_o16_p(json), rate_c13_alpha(json), &
           rate_c12_alpha(json), rate_n14_alpha(json), &
           rate_triple_alpha(json)
      double precision, intent(in) :: frac_c12_alpha(json)
      double precision, intent(in) :: shell_mass(json)
      double precision, intent(out) :: composition(15,json)
      double precision, intent(in) :: timestep_years

! common/burtol/: min_abundance, absolute_tolerance, relative_tolerance,
! max_burn_iterations (originally CMIN,ABSTOL,RELTOL,KEMMAX). First
! appearance of this common block in the converted sources.
      double precision :: min_abundance, absolute_tolerance, relative_tolerance
      integer :: max_burn_iterations
      common/burtol/ min_abundance, absolute_tolerance, relative_tolerance, &
           max_burn_iterations



! common/oldmod/: only old_composition is used here. Naming matches
! eqburn.f90/dburn.f90.
      double precision :: old_pressure(json), old_temperature(json), &
           old_radius(json), old_luminosity(json), old_density(json), &
           old_composition(15,json), old_shell_mass(json), old_teff
      logical :: old_convective_flag(json), old_cz_flag(json)
      integer :: old_num_zones
      common/oldmod/ old_pressure, old_temperature, old_radius, &
           old_luminosity, old_density, old_composition, old_shell_mass, &
           old_convective_flag, old_cz_flag, old_teff, old_num_zones

! system_matrix(56): flattened 7-row x 8-column augmented matrix
! passed to simeqc -- columns 1-7 are the Jacobian of the 7 implicit
! burning equations w.r.t. the 7 solved species, column 8 (elements
! 50-56) is the residual/RHS vector. correction(7) is equivalenced
! onto system_matrix(50:56) so that simeqc's in-place solve leaves the
! Newton correction directly in correction().
      double precision :: system_matrix(56), correction(7)
      equivalence (system_matrix(50),correction(1))
      double precision :: abundance(7), avg_abundance(11)
      save

      double precision :: total_shell_mass
      integer :: species_idx, zone_idx, write_zone_idx
      double precision :: timestep_gyr, timestep_gyr_3, timestep_gyr_4, &
           timestep_gyr_12, timestep_gyr_13, timestep_gyr_14, timestep_gyr_16
      double precision :: min_abundance_local
      integer :: iteration_count
      double precision :: gr_pp, gr_he3_he3, gr_he3_he4, gr_c12_p, gr_c13_p, &
           gr_n14_p, gr_o16_p, gr_c13_alpha, gr_zero9, gr_c12_alpha, &
           gr_n14_alpha, gr_triple_alpha, gr_zero13
      double precision :: branch_frac_c12, branch_frac_o16
      double precision :: x_start, he3_start, y_start, c12_start, &
           c13_start, n14_start, o16_start
      integer :: any_nonzero_flag, rhs_column_idx, mat_idx
      double precision :: max_abs_change, max_relative_change, &
           relative_change
      double precision :: o18_new, new_metal_fraction
      integer :: solved_species_idx

      if(zone_begin.ne.zone_end) then
!  homogenize convection zones.
!  avg_abundance is the mass-weighted average abundance for the cz.
!  total_shell_mass is the total mass of the cz.
!  initialize sums.
         total_shell_mass = 0.0d0
         do 1 species_idx = 1,11
            avg_abundance(species_idx) = 0.0d0
    1    continue
         do 5 zone_idx = zone_begin,zone_end
            total_shell_mass = total_shell_mass + shell_mass(zone_idx)
            do 3 species_idx = 1,11
               avg_abundance(species_idx) = avg_abundance(species_idx)+ &
                    old_composition(species_idx,zone_idx)*shell_mass(zone_idx)
    3       continue
    5    continue
         do 7 species_idx = 1,11
            avg_abundance(species_idx) = avg_abundance(species_idx)/total_shell_mass
    7    continue
      else
         do 9 species_idx = 1,11
            avg_abundance(species_idx) = &
                 old_composition(species_idx,zone_begin)
    9    continue
      endif
!  skip burning calculations if starting shell below t cutoff for reactions.
      if(log_temperature(zone_begin).lt.tcut(1)) then
         do 13 zone_idx = zone_begin,zone_end
            do 11 species_idx = 1,11
               composition(species_idx,zone_idx) = avg_abundance(species_idx)
   11       continue
   13    continue
         goto 200
      endif
!
!  set up numerical parameters.
!
!  timestep in gigayears. the other delts are this timestep multiplied
!   by the atomic mass of the different species (in amu).
      timestep_gyr=timestep_years*1.0d-9
      timestep_gyr_3 = 3.0d0*timestep_gyr
      timestep_gyr_4 = 4.0d0*timestep_gyr
      timestep_gyr_12 = 12.0d0*timestep_gyr
      timestep_gyr_13 = 13.0d0*timestep_gyr
      timestep_gyr_14 = 14.0d0*timestep_gyr
      timestep_gyr_16 = 16.0d0*timestep_gyr
      min_abundance_local=min_abundance
!  counter for the number of iterations.
      iteration_count=0
!
!  nuclear reaction rates.
!
!  these reactions are no longer included and are zeroed out.
      gr_zero9 = 0.0d0
      gr_zero13 = 0.0d0
      if(zone_begin.eq.zone_end) then
!  pp
         gr_pp = rate_pp(zone_begin)
!  he3,he3
         gr_he3_he3 = rate_he3_he3(zone_begin)
!  he3,he4
         gr_he3_he4 = rate_he3_he4(zone_begin)
!  c12,p
         gr_c12_p = rate_c12_p(zone_begin)
!  c13,p
         gr_c13_p = rate_c13_p(zone_begin)
!  n14,p + n15,p
         gr_n14_p = rate_n14_p(zone_begin)
!  o16,p + o17,p.
         gr_o16_p = rate_o16_p(zone_begin)
!  c13,alpha
         gr_c13_alpha = rate_c13_alpha(zone_begin)
!  o16,alpha (not used)
!        gr_zero9 = rate_zero9(zone_begin)
!  c12,alpha
         gr_c12_alpha = rate_c12_alpha(zone_begin)
!  n14,alpha
         gr_n14_alpha = rate_n14_alpha(zone_begin)
!  triple alpha
         gr_triple_alpha = rate_triple_alpha(zone_begin)
!  c12,c12 (not used)
!        gr_zero13 = rate_zero13(zone_begin)
!  branching ratio for n15,p :
!  branch_frac_c12 = fraction going to c12+alpha, 1-branch_frac_c12 = fraction going to o16
         branch_frac_c12 = frac_c12_alpha(zone_begin)
         branch_frac_o16 = 1.0d0 - branch_frac_c12
      else
!  use the mass-weighted average rates for the cz
         gr_pp = 0.0d0
         gr_he3_he3 = 0.0d0
         gr_he3_he4 = 0.0d0
         gr_c12_p = 0.0d0
         gr_c13_p = 0.0d0
         gr_n14_p = 0.0d0
         gr_o16_p = 0.0d0
         gr_c13_alpha = 0.0d0
         gr_c12_alpha = 0.0d0
         gr_n14_alpha = 0.0d0
         gr_triple_alpha = 0.0d0
         branch_frac_c12 = 0.0d0
         do 15 zone_idx = zone_begin,zone_end
            gr_pp = gr_pp + shell_mass(zone_idx)*rate_pp(zone_idx)
            gr_he3_he3 = gr_he3_he3 + shell_mass(zone_idx)*rate_he3_he3(zone_idx)
            gr_he3_he4 = gr_he3_he4 + shell_mass(zone_idx)*rate_he3_he4(zone_idx)
            gr_c12_p = gr_c12_p + shell_mass(zone_idx)*rate_c12_p(zone_idx)
            gr_c13_p = gr_c13_p + shell_mass(zone_idx)*rate_c13_p(zone_idx)
            gr_n14_p = gr_n14_p + shell_mass(zone_idx)*rate_n14_p(zone_idx)
            gr_o16_p = gr_o16_p + shell_mass(zone_idx)*rate_o16_p(zone_idx)
            gr_c13_alpha = gr_c13_alpha + shell_mass(zone_idx)*rate_c13_alpha(zone_idx)
            gr_c12_alpha = gr_c12_alpha + shell_mass(zone_idx)*rate_c12_alpha(zone_idx)
            gr_n14_alpha = gr_n14_alpha + shell_mass(zone_idx)*rate_n14_alpha(zone_idx)
            gr_triple_alpha = gr_triple_alpha + shell_mass(zone_idx)*rate_triple_alpha(zone_idx)
            branch_frac_c12 = branch_frac_c12 + shell_mass(zone_idx)*frac_c12_alpha(zone_idx)
   15    continue
         gr_pp = gr_pp/total_shell_mass
         gr_he3_he3 = gr_he3_he3/total_shell_mass
         gr_he3_he4 = gr_he3_he4/total_shell_mass
         gr_c12_p = gr_c12_p/total_shell_mass
         gr_c13_p = gr_c13_p/total_shell_mass
         gr_n14_p = gr_n14_p/total_shell_mass
         gr_o16_p = gr_o16_p/total_shell_mass
         gr_c13_alpha = gr_c13_alpha/total_shell_mass
         gr_c12_alpha = gr_c12_alpha/total_shell_mass
         gr_n14_alpha = gr_n14_alpha/total_shell_mass
         gr_triple_alpha = gr_triple_alpha/total_shell_mass
         branch_frac_c12 = branch_frac_c12/total_shell_mass
         branch_frac_o16 = 1.0d0 - branch_frac_c12
      endif
!
!  initial abundances of species.
!
      x_start = avg_abundance(1)
      he3_start = avg_abundance(4)
      y_start = avg_abundance(2)
      c12_start = avg_abundance(5)
      c13_start = avg_abundance(6)
      n14_start = avg_abundance(7)
      o16_start = avg_abundance(9)
!  starting guess for final abundances of species :
!  for calls in crrect (itlvl>1) use the previous result.
!  otherwise use the final abundances found for the previous shell
!  unless i = 1.
      if(zone_begin.ne.1) then
         abundance(1) = composition(1,zone_begin-1)
         abundance(2) = composition(4,zone_begin-1)
         abundance(3) = composition(2,zone_begin-1)
         abundance(4) = composition(5,zone_begin-1)
         abundance(5) = composition(6,zone_begin-1)
         abundance(6) = composition(7,zone_begin-1)
         abundance(7) = composition(9,zone_begin-1)
      else
         abundance(1) = composition(1,zone_begin)
         abundance(2) = composition(4,zone_begin)
         abundance(3) = composition(2,zone_begin)
         abundance(4) = composition(5,zone_begin)
         abundance(5) = composition(6,zone_begin)
         abundance(6) = composition(7,zone_begin)
         abundance(7) = composition(9,zone_begin)
      endif
!
!  iteration for new abundances.
!
   10 continue
!  functions to be minimized. these are of the form
!  system_matrix(#)=abundance(end step)-abundance (start step)-timestep*(d(species)/dt)
!  because the burning rates d(species)/dt use the abundances at the
!  end of the timestep abundance(1)-abundance(7), this is an implicit scheme.
!  equation for h (x)
      system_matrix(50) = abundance(1)-x_start-timestep_gyr* &
           (-3.d0*gr_pp*abundance(1)**2+2.d0*gr_he3_he3*abundance(2)**2 &
           -gr_he3_he4*abundance(2)*abundance(3) &
           -gr_c12_p*abundance(1)*abundance(4)-gr_c13_p*abundance(1)*abundance(5) &
           -2.d0*abundance(1)*(gr_n14_p*abundance(6)+gr_o16_p*abundance(7)))
!  he3
      system_matrix(51)= abundance(2)-he3_start-timestep_gyr_3* &
           (gr_pp*abundance(1)**2-2.d0*gr_he3_he3*abundance(2)**2 &
           -gr_he3_he4*abundance(2)*abundance(3))
!  he4 (y)
      system_matrix(52) = abundance(3)-y_start-timestep_gyr_4* &
           (gr_he3_he3*abundance(2)**2+gr_he3_he4*abundance(2)*abundance(3) &
           +branch_frac_c12*gr_n14_p*abundance(6)*abundance(1) &
           +gr_o16_p*abundance(1)*abundance(7)-gr_c12_alpha*abundance(3)*abundance(4) &
           -gr_c13_alpha*abundance(3)*abundance(5)-gr_n14_alpha*abundance(3)*abundance(6) &
           -3.d0*gr_triple_alpha*abundance(3)**3)
!  c12
      system_matrix(53) = abundance(4)-c12_start-timestep_gyr_12* &
           (gr_triple_alpha*abundance(3)**3-gr_c12_p*abundance(1)*abundance(4) &
           -gr_c12_alpha*abundance(3)*abundance(4) &
           +branch_frac_c12*gr_n14_p*abundance(1)*abundance(6))
!  c13
      system_matrix(54)=abundance(5)-c13_start-timestep_gyr_13* &
           (gr_c12_p*abundance(1)*abundance(4)-gr_c13_p*abundance(1)*abundance(5) &
           -gr_c13_alpha*abundance(3)*abundance(5))
!  n14
      system_matrix(55) = abundance(6)-n14_start-timestep_gyr_14* &
           (gr_c13_p*abundance(1)*abundance(5)+gr_o16_p*abundance(1)*abundance(7) &
           -gr_n14_p*abundance(1)*abundance(6) - gr_n14_alpha*abundance(3)*abundance(6))
!  o16*abundance(4)
      system_matrix(56) = abundance(7)-o16_start-timestep_gyr_16* &
           (branch_frac_o16*gr_n14_p*abundance(1)*abundance(6) &
           +gr_c12_alpha*abundance(3)*abundance(4)-gr_o16_p*abundance(1)*abundance(7) &
           +gr_c13_alpha*abundance(3)*abundance(5))
!  note that the o18 abundance is computed explicitly at the end of the
!  subroutine.
!  the remaining values of system_matrix correspond to the derivatives of the above 7
!  functions with respect to the 7 species.
!  d(rates)/dx
      system_matrix(1) = 1.d0+timestep_gyr* &
           (6.d0*gr_pp*abundance(1)+gr_c12_p*abundance(4)+gr_c13_p*abundance(5) &
           +2.d0*(gr_n14_p*abundance(6)+gr_o16_p*abundance(7)))
      system_matrix(2) = -2.d0*timestep_gyr_3*gr_pp*abundance(1)
      system_matrix(3) = -timestep_gyr_4*(branch_frac_c12*gr_n14_p*abundance(6)+gr_o16_p*abundance(7))
      system_matrix(4) = -timestep_gyr_12*(branch_frac_c12*gr_n14_p*abundance(6)-gr_c12_p*abundance(4))
      system_matrix(5) = -timestep_gyr_13*(gr_c12_p*abundance(4)-gr_c13_p*abundance(5))
      system_matrix(6) = -timestep_gyr_14*(gr_c13_p*abundance(5)+gr_o16_p*abundance(7)-gr_n14_p*abundance(6))
      system_matrix(7) = -timestep_gyr_16*(branch_frac_o16*gr_n14_p*abundance(6)-gr_o16_p*abundance(7))
!  d(rates)/dhe3
      system_matrix(8) = -timestep_gyr*(4.d0*gr_he3_he3*abundance(2)-gr_he3_he4*abundance(3))
      system_matrix(9) = 1.d0+timestep_gyr_3*(4.d0*gr_he3_he3*abundance(2)+gr_he3_he4*abundance(3))
      system_matrix(10) = -timestep_gyr_4*(2.d0*gr_he3_he3*abundance(2)+gr_he3_he4*abundance(3))
      system_matrix(11) = 0.d0
      system_matrix(12) = 0.d0
      system_matrix(13) = 0.d0
      system_matrix(14) = 0.d0
!  d(rates)/dhe4
      system_matrix(15) = timestep_gyr*gr_he3_he4*abundance(2)
      system_matrix(16) = timestep_gyr_3*gr_he3_he4*abundance(2)
      system_matrix(17) = 1.d0-timestep_gyr_4* &
           (gr_he3_he4*abundance(2)-gr_c12_alpha*abundance(4)-gr_c13_alpha*abundance(5) &
           -gr_n14_alpha*abundance(6)-9.d0*gr_triple_alpha*abundance(3)**2)
      system_matrix(18) = -timestep_gyr_12*(3.d0*gr_triple_alpha*abundance(3)**2-gr_c12_alpha*abundance(4))
      system_matrix(19) = timestep_gyr_13*gr_c13_alpha*abundance(5)
      system_matrix(20) = timestep_gyr_14*gr_n14_alpha*abundance(6)
      system_matrix(21) = -timestep_gyr_16*(gr_c12_alpha*abundance(4)+gr_c13_alpha*abundance(5))
!  d(rates)/dc12
      system_matrix(22) = timestep_gyr*gr_c12_p*abundance(1)
      system_matrix(23) = 0.d0
      system_matrix(24) = timestep_gyr_4*gr_c12_alpha*abundance(3)
      system_matrix(25) = 1.d0+timestep_gyr_12*(gr_c12_p*abundance(1)+gr_c12_alpha*abundance(3))
      system_matrix(26) = -timestep_gyr_13*gr_c12_p*abundance(1)
      system_matrix(27) = 0.d0
      system_matrix(28) = -timestep_gyr_16*gr_c12_alpha*abundance(3)
!  d(rates)/dc13
      system_matrix(29) = timestep_gyr*gr_c13_p*abundance(1)
      system_matrix(30) = 0.d0
      system_matrix(31) = timestep_gyr_4*gr_c13_alpha*abundance(3)
      system_matrix(32) = 0.d0
      system_matrix(33) = 1.d0+timestep_gyr_13*(gr_c13_p*abundance(1)+gr_c13_alpha*abundance(3))
      system_matrix(34) = -timestep_gyr_14*gr_c13_p*abundance(1)
      system_matrix(35) = -timestep_gyr_16*gr_c13_alpha*abundance(3)
!  d(rates)/dn14
      system_matrix(36) = 2.d0*timestep_gyr*gr_n14_p*abundance(1)
      system_matrix(37) = 0.d0
      system_matrix(38) = -timestep_gyr_4*(branch_frac_c12*gr_n14_p*abundance(1)-gr_n14_alpha*abundance(3))
      system_matrix(39) = -timestep_gyr_12*branch_frac_c12*gr_n14_p*abundance(1)
      system_matrix(40) = 0.d0
      system_matrix(41) = 1.d0+timestep_gyr_14*(gr_n14_p*abundance(1)+gr_n14_alpha*abundance(3))
      system_matrix(42) = -branch_frac_o16*timestep_gyr_16*gr_n14_p*abundance(1)
!  d(rates)/do16
      system_matrix(43) = 2.d0*timestep_gyr*gr_o16_p*abundance(1)
      system_matrix(44) = 0.d0
      system_matrix(45) = -timestep_gyr_4*gr_o16_p*abundance(1)
      system_matrix(46) = 0.d0
      system_matrix(47) = 0.d0
      system_matrix(48) = -timestep_gyr_14*gr_o16_p*abundance(1)
      system_matrix(49) = 1.d0+timestep_gyr_16*gr_o16_p*abundance(1)
!  reduce the 7x8 system to find the corrections to the relative
!  abundances.
      any_nonzero_flag = 0
      rhs_column_idx = 1
      do 20 mat_idx = 1,49
         if(mat_idx.eq.rhs_column_idx) then
            rhs_column_idx = rhs_column_idx+8
         else if(dabs(system_matrix(mat_idx)).lt.1.d-8) then
            system_matrix(mat_idx) = 0.d0
         else
            any_nonzero_flag = 1
         endif
   20 continue
      if(any_nonzero_flag.ne.0) call simeqc(system_matrix,8,7)
!  check to see if the system has converged within the desired tolerances.
      max_abs_change = 0.d0
      max_relative_change = 0.d0
      do 30 solved_species_idx = 1,7
         abundance(solved_species_idx) = abundance(solved_species_idx)-correction(solved_species_idx)
         if(dabs(correction(solved_species_idx)).ge.max_abs_change) max_abs_change = dabs(correction(solved_species_idx))
         if(abundance(solved_species_idx).lt.min_abundance_local) then
            abundance(solved_species_idx) = 0.d0
         else
            relative_change = dabs(correction(solved_species_idx)/abundance(solved_species_idx))
            if(relative_change.ge.max_relative_change) max_relative_change = relative_change
         endif
   30 continue
      if(max_abs_change.ge.absolute_tolerance.or.max_relative_change.ge.relative_tolerance) then
!  system not converged - see if maximum number of iterations exceeded.
         iteration_count = iteration_count+1
         if(iteration_count.ge.max_burn_iterations) then
!  mhp 10/02 iu not defined
!            WRITE (short_file_unit,1000) iu
            write (short_file_unit,1000) zone_begin
 1000       format(1X,39('>'),40('<')/1X,'ERROR IN SUBROUTINE KEMCOM'/ &
            1X,'UNABLE TO SOLVE FOR NEW ABUNDANCES IN SHELL',I4/1X, &
            'RUN STOPPED AFTER 50 ATTEMPTS')
            stop
         else
            goto 10
         endif
      endif
!  system has converged.
!  update composition matrix.
!  update o18.
      o18_new = avg_abundance(11)+18.0d0*timestep_gyr*gr_n14_alpha*abundance(3)*abundance(6)
!  change metal abundance if x<5.0d-7.
      if(old_composition(1,zone_end).lt.5.0d-7) then
         new_metal_fraction = 1.0d0-abundance(1)-abundance(2)-abundance(3)
      else
         new_metal_fraction = avg_abundance(3)
      endif
      do 40 write_zone_idx = zone_begin,zone_end
         composition(1,write_zone_idx) = abundance(1)
         composition(4,write_zone_idx) = abundance(2)
         composition(2,write_zone_idx) = abundance(3)
         composition(5,write_zone_idx) = abundance(4)
         composition(6,write_zone_idx) = abundance(5)
         composition(7,write_zone_idx) = abundance(6)
         composition(9,write_zone_idx) = abundance(7)
         composition(3,write_zone_idx) = new_metal_fraction
         composition(11,write_zone_idx) = o18_new
   40 continue
  200 return
end subroutine kemcom
