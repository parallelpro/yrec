!----------------------------------------------------------------------
! surfbc
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original surfbc.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Maintains the (Teff, L) triangle used to interpolate the envelope
! solution (P, T, R at the fitting point) around the current model's
! (TEFFL, BL): checks whether the requested point still lies inside
! the existing triangle, re-triangulates and recomputes any envelope
! vertices as needed via atm_lib.f90's atm_get, and rebuilds the bilinear
! interpolation coefficients CFENV. Also handles the automatic
! gray-atmosphere fallback when Teff moves outside the range of a
! tabulated (Kurucz/Allard) T-tau relation.
!
! hydrogen_fraction/metal_fraction/pressure_rotation_factor/
! temperature_rotation_factor are intent(inout), not intent(in): they
! are relayed unchanged into atm_get's own same-named dummies, which
! atm_get's own header documents as intent(inout) because they flow
! through to bsstep/mmid's callback machinery. Declaring them
! intent(in) here was a pre-existing bug of the same shape documented
! in atm_get's header (and eos/eqstat.f90's metal_fraction fix) --
! silently tolerated while atm_get (as envint) was a bare external
! subroutine with no interface to check against, surfaced once
! atm_lib.f90 gave it one. core/crrect.f90 (this routine's only
! caller) already passes real local variables for these positions, so
! widening the intent here changes nothing about how it's called.
subroutine surfbc(tri_teffl, tri_logl, envelope_coeffs, &
     vtx_logp, vtx_logt, &
     vtx_logr, tri_orientation, stored_vertex_index, &
     stored_envelope_state, &
     start_new_triangle, reset_triangle, saha_state, env_call_count, &
     atm_call_count, log10_star_mass, luminosity_linear, &
     log10_teff, hydrogen_fraction, metal_fraction, &
     pressure_rotation_factor, temperature_rotation_factor, &
     envelope_recomputed_flag, log10_pressure_limit, convective_flag, &
     zone_index)

! INPUTS   start_new_triangle = .T.    START UP WITH 3 NEW ENVELOPES ABOUT(TEFFL,BL)
! INPUTS   reset_triangle = .T.  REDO ALL 3 ENVELOPES AND RETRIANGULATE IF NEED
! BOTH start_new_triangle AND reset_triangle ARE RESET TO .FALSE.
      use atm_lib
      use envint_lib, only: atm_get
      use atm_table_lib
      use envelope_comp_lib
      use luout_lib
      use const_lib
      implicit none
      integer, parameter :: json=5000

      double precision, intent(inout) :: tri_teffl(3), tri_logl(3), &
           envelope_coeffs(9)
      double precision, intent(inout) :: vtx_logp(3), &
           vtx_logt(3), vtx_logr(3)
      double precision, intent(inout) :: tri_orientation
      integer, intent(inout) :: stored_vertex_index
      double precision, intent(inout) :: stored_envelope_state(4)
      logical, intent(inout) :: start_new_triangle, reset_triangle
      integer, intent(inout) :: saha_state, env_call_count, atm_call_count
      double precision, intent(in) :: log10_star_mass, luminosity_linear
      double precision, intent(inout) :: log10_teff
      double precision, intent(inout) :: hydrogen_fraction, metal_fraction, &
           pressure_rotation_factor, temperature_rotation_factor
      logical, intent(out) :: envelope_recomputed_flag
      double precision, intent(in) :: log10_pressure_limit
      logical, intent(in) :: convective_flag(json)
      integer, intent(in) :: zone_index

      logical :: tri_vertex_valid(3)

      integer :: numenv
      data numenv/0/

! G Somers END

      save

! --- locals ---
      integer :: i1, i2, i3, i, vertex_being_computed, j
      double precision :: temp, temp1, temp2, tri_err
      logical :: envelope_needs_recompute, print_envelope_flag, &
           save_boundary_flag, pulse_print_flag
      double precision :: b, gl, rl, adjusted_teffl

! MHP 9/01
      if (atm_choice.eq.3) then
         if (log10_teff.ge.3.95d0) then
            write(*,*)
            write(*,5) log10_teff
 5          format('LOG TEFF OF ',F7.3,' ABOVE 3.95 - SWITCH' &
                   ,'TO GRAY ATMOSPHERE BOUNDARY CONDITION')
            write(*,*)
            atm_choice = 0
            start_new_triangle = .true.
! MHP 06/13 Remember that flag is switched
            use_ttau_relation = .true.
         endif
      endif
      if (atm_choice.eq.4) then
         if (log10_teff.ge.atm_table%allard_al_teffl_max) then
            write(*,*)
            write(*,7) log10_teff, atm_table%allard_al_teffl_max
 7          format('LOG TEFF OF ',F7.3,' ABOVE Allard Table max ',F7.3 &
                   ,'  - SWITCH TO GRAY ATMOSPHERE BOUNDARY CONDITION')
            write(*,*)
            atm_choice = 0
            start_new_triangle = .true.
! MHP 06/13 Remember that flag is switched
            use_ttau_relation = .true.
         endif
      endif
      if (use_ttau_relation) then
         if (atm_choice_initial.eq.3.and.log10_teff.lt.3.95d0) then
            atm_choice = atm_choice_initial
            use_ttau_relation = .false.
            write(*,9) log10_teff
 9          format('LOG TEFF OF ',F7.3,' BELOW 3.95 - SWITCH' &
           ,' BACK TO KURUCZ ATMOSPHERE BOUNDARY CONDITION')
         else if (atm_choice_initial.eq.4.and.log10_teff.lt.atm_table%allard_al_teffl_max) then
            atm_choice = atm_choice_initial
            use_ttau_relation = .false.
            write(*,11) log10_teff, atm_table%allard_al_teffl_max
 11         format('LOG TEFF OF ',F7.3,' below Allard Table max ',F7.3 &
           ,'  - SWITCH BACK TO ALLARD ATMOSPHERE BOUNDARY CONDITION')
         endif
      endif
      envelope_recomputed_flag = .false.
! IF LNEW0, REDO ALL 3 ENVELOPES EVERY MODEL
      if (start_new_triangle.or.reset_triangle) then
! REDO ALL THREE ENVELOPES
       tri_vertex_valid(1) = .false.
       tri_vertex_valid(2) = .false.
       tri_vertex_valid(3) = .false.
       stored_vertex_index = 0
       tri_err = 0.0d0
       reset_triangle = .false.
      else
       tri_vertex_valid(1) = .true.
       tri_vertex_valid(2) = .true.
       tri_vertex_valid(3) = .true.
       tri_err = 0.0d0
      endif
      if (start_new_triangle) then
! STARTING PROCEDURE
       tri_orientation = +1.0d0
       tri_teffl(3) = log10_teff
       tri_teffl(1) = tri_teffl(3) - 0.5d0*tri_delta_teffl
       tri_teffl(2) = tri_teffl(1) + tri_delta_teffl
       tri_logl(3) = luminosity_linear + 0.5d0*tri_delta_logl
       tri_logl(1) = tri_logl(3) - tri_delta_logl
       tri_logl(2) = tri_logl(1)
       start_new_triangle = .false.
      else
! CHECK TRIANGULATION OF POINT (TEFFL,BL)
 10      continue
       do i1 = 1,3
          i2 = mod(i1,3) + 1
          i3 = mod(i2,3) + 1
          temp = tri_orientation*((tri_logl(i2)-tri_logl(i3))*(log10_teff-tri_teffl(i2)) + &
                 (tri_teffl(i3)-tri_teffl(i2))*(luminosity_linear-tri_logl(i2)) )
          if (temp.lt.-tri_err) then
             tri_orientation = -tri_orientation
             tri_teffl(i1) = tri_teffl(i2) + tri_teffl(i3) - tri_teffl(i1)
             tri_logl(i1) = tri_logl(i2) + tri_logl(i3) - tri_logl(i1)
             tri_vertex_valid(i1) = .false.
             tri_err = 0.0d0
             goto 10
          endif
 20      continue
       end do
      endif
! COMPUTE NEW ENVELOPES IF NECESSARY
      envelope_needs_recompute = .false.
      do i = 1,3
       if (.not.tri_vertex_valid(i)) then
! NEW ENVELOPE NEEDED
          envelope_recomputed_flag = .true.
          log10_teff = tri_teffl(i)
          b = dexp(ln10*tri_logl(i))
          rl = 0.5d0*(tri_logl(i) + log10_solar_luminosity - 4.0d0*log10_teff - c4pil - csigl)
          gl = cgl + log10_star_mass - rl - rl
          numenv = numenv + 1
! G Somers 11/14, LPENV FUNCTIONALITY NOW INCLUDED IN LSTENV.
            print_envelope_flag = .false.
!          LPENV = MOD(NUMENV,NPENV).EQ.0
!          IF(.NOT.LSTORE) LPENV = .FALSE.
!          IF(LPENV) THEN
!             WRITE(ISHORT,30) I,TRIT(I),TRIL(I)
! 30            FORMAT(/,' NEW ENVELOPE NO.',I2,' LOG(TE) =',F9.5,
!     *              '  LOG(L) =',F8.5)
!          ENDIF
! G Somers END
          vertex_being_computed=i
          save_boundary_flag = .true.
! DBG PULSE: DO NOT DO PULSE OUTPUT
            pulse_print_flag = .false.
! G Somers 10/14, FOR SPOTTED RUNS, FIND THE
! PRESSURE AT THE AMBIENT TEMPERATURE ATEFFL
          if (convective_flag(zone_index).and.spot_filling_factor.ne.0.0.and.spot_temp_contrast.ne.1.0) then
               adjusted_teffl = log10_teff - 0.25*log10(spot_filling_factor * spot_temp_contrast**4.0 + 1.0 - spot_filling_factor)
          else
             adjusted_teffl = log10_teff
          endif
          call atm_get(b,pressure_rotation_factor,temperature_rotation_factor,gl, &
                 log10_star_mass,vertex_being_computed,print_envelope_flag, &
                 save_boundary_flag,log10_pressure_limit,rl, &
                 adjusted_teffl,hydrogen_fraction,metal_fraction, &
                 stored_envelope_state,stored_vertex_index,atm_call_count, &
                 env_call_count,saha_state,vtx_logp, &
                 vtx_logr,vtx_logt,pulse_print_flag)
! G Somers END
          envelope_needs_recompute = .true.
       endif
 40   continue
      end do
      if (envelope_needs_recompute) then
! RECOMPUTE COEFFICIENTS
       temp = 1.0d0/(vtx_logt(2)-vtx_logt(1))
       temp1 = (vtx_logt(3)-vtx_logt(1))*temp
       temp2 = 1.0d0/(vtx_logp(3)-vtx_logp(1)-temp1*(vtx_logp(2)-vtx_logp(1)))
       envelope_coeffs(1) = (  vtx_logr(3)-  vtx_logr(1)-temp1*(  vtx_logr(2)-  vtx_logr(1)))*temp2
       envelope_coeffs(4) = (tri_logl(3)-tri_logl(1)-temp1*(tri_logl(2)-tri_logl(1)))*temp2
       envelope_coeffs(7) = (tri_teffl(3)-tri_teffl(1)-temp1*(tri_teffl(2)-tri_teffl(1)))*temp2
       envelope_coeffs(2) = (  vtx_logr(2)-  vtx_logr(1)-envelope_coeffs(1)*(vtx_logp(2)-vtx_logp(1)))*temp
       envelope_coeffs(5) = (tri_logl(2)-tri_logl(1)-envelope_coeffs(4)*(vtx_logp(2)-vtx_logp(1)))*temp
       envelope_coeffs(8) = (tri_teffl(2)-tri_teffl(1)-envelope_coeffs(7)*(vtx_logp(2)-vtx_logp(1)))*temp
       envelope_coeffs(3) =   vtx_logr(1) - envelope_coeffs(1)*vtx_logp(1) - envelope_coeffs(2)*vtx_logt(1)
       envelope_coeffs(6) = tri_logl(1) - envelope_coeffs(4)*vtx_logp(1) - envelope_coeffs(5)*vtx_logt(1)
       envelope_coeffs(9) = tri_teffl(1) - envelope_coeffs(7)*vtx_logp(1) - envelope_coeffs(8)*vtx_logt(1)
       write(short_file_unit,50)(i,tri_vertex_valid(i),tri_teffl(i),tri_logl(i),vtx_logp(i),vtx_logt(i),vtx_logr(i), &
              (envelope_coeffs(i+i+i-3+j),j=1,3), i=1,3)
 50      format(' ENVELOPE TRIANGLE  LOG(TE)  LOG(L)   LOG(P)   LOG(T)', &
              3X,'LOG(R)   COEFFICIENTS'/I14,L1,2X,5F9.5,8X,'LOG(R)', &
              3F10.5/ I14,L1,2X,5F9.5,5X,'LOG(L/L0)',3F10.5/ &
              I14,L1,2X,5F9.5,7X,'LOG(TE)',3F10.5)
! 60      FORMAT(' COUNTS  KATM',I5,'  KENV',I5,'  KSAHA',I5)
      else
       write(short_file_unit,70)
 70      format(' ENVELOPE TRIANGLE NOT CHANGED')
      endif
      return
end subroutine surfbc
