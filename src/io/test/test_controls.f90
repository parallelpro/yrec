!----------------------------------------------------------------------
! test_controls
!----------------------------------------------------------------------
! Reads one deck through read_controls, exactly as a run does, and
! prints every resulting control: each element of star%ctrl and
! star%job (controls_dump.inc, generated from their declarations), the
! star fields read_controls sets directly, the Bulirsch-Stoer settings
! and the unit numbers. Doubles print at full precision, so two dumps
! agree exactly when the controls are bit-identical.
!   test_controls                  the pristine defaults (no deck read)
!   test_controls deck.inlist      a single-file inlist
!   test_controls deck.nml1 deck.nml2   a legacy CONTROL/PHYSICS pair
! The output files a run opens go under $YREC_OUTPUT (test_controls.py
! points it at a scratch directory). test_controls.py runs this on every
! deck in examples/ and testsuite/.
program test_controls
      use star_info_lib, only: star, control_nml_override, physics_nml_override
      use integration_parameters, only: tolerance_fraction, max_stage_index, extrap_order
      use io_units, only: last_model_unit, run_log_unit, terminal_unit
      use controls_reader, only: read_controls
      implicit none
      character(len=256) :: deck1, deck2
      integer :: ierr, u, i1, i2, nargs

      u = 6
      nargs = command_argument_count()
      if (nargs >= 1) then
         call get_command_argument(1, deck1)
         deck2 = deck1
         if (nargs >= 2) call get_command_argument(2, deck2)
         control_nml_override = deck1
         physics_nml_override = deck2
! read_controls writes its banner to standard output; the dump goes to
! unit 6 after a marker so the two can be told apart
         call read_controls(star, ierr)
         write(u,'(a)') "# controls"
         write(u,'(a,i0)') "ierr = ", ierr
      else
         write(u,'(a)') "# controls (pristine defaults)"
      end if
      include 'controls_dump.inc'
! the table paths and mixture weights read_controls stores in star%job
      call path('alex06_table_path', star%job%alex06_table_path)
      call path('allard_table_path', star%job%allard_table_path)
      call path('atm_table_path', star%job%atm_table_path)
      call path('fermi_table_path', star%job%fermi_table_path)
      call path('kurucz_table_path', star%job%kurucz_table_path)
      call path('kurucz_table2_path', star%job%kurucz_table2_path)
      call path('laol_table_path', star%job%laol_table_path)
      call path('laol_table2_path', star%job%laol_table2_path)
      call path('opal95_table_path', star%job%opal95_table_path)
      call path('opal92_table_path', star%job%opal92_table_path)
      call path('zams_a_table_path', star%job%zams_a_table_path)
      call path('zams_b_table_path', star%job%zams_b_table_path)
      call path('zams_c_table_path', star%job%zams_c_table_path)
      call path('centre1_table_path', star%job%centre1_table_path)
      call path('centre2_table_path', star%job%centre2_table_path)
      call path('centre3_table_path', star%job%centre3_table_path)
      call path('centre4_table_path', star%job%centre4_table_path)
      call path('centre5_table_path', star%job%centre5_table_path)
      call path('opal92_table2_path', star%job%opal92_table2_path)
      call path('pure_z_table_path', star%job%pure_z_table_path)
      call path('scv_h_table_path', star%job%scv_h_table_path)
      call path('scv_he_table_path', star%job%scv_he_table_path)
      call path('scv_z_table_path', star%job%scv_z_table_path)
      call path('conductive_table_path', star%job%conductive_table_path)
      call path('opal_eos_table_path', star%job%opal_eos_table_path)
      do i1 = 1, size(star%job%alex95_table_paths)
         call path('alex95_table_paths('//itoa(i1)//')', star%job%alex95_table_paths(i1))
      end do
      do i1 = 1, size(star%job%mixture_weights)
         write(u,'(a," = ",es25.17e3)') 'job%mixture_weights('//itoa(i1)//')', &
              star%job%mixture_weights(i1)
      end do
      write(u,'(a," = ",i0)') 'star%atm_choice_initial', star%atm_choice_initial
      write(u,'(a," = ",l1)') 'star%use_ttau_relation', star%use_ttau_relation
      write(u,'(a," = ",es25.17e3)') 'star%tenv', star%tenv
      do i1 = 1, size(star%cross_section_scale)
         write(u,'(a," = ",es25.17e3)') 'star%cross_section_scale('//itoa(i1)//')', &
              star%cross_section_scale(i1)
      end do
      write(u,'(a," = ",es25.17e3)') 'star%envelope_hydrogen_fraction', star%envelope_hydrogen_fraction
      write(u,'(a," = ",es25.17e3)') 'star%envelope_metal_fraction', star%envelope_metal_fraction
      write(u,'(a," = ",es25.17e3)') 'star%solar_luminosity_cgs', star%solar_luminosity_cgs
      write(u,'(a," = ",es25.17e3)') 'star%solar_radius_cgs', star%solar_radius_cgs
      write(u,'(a," = ",es25.17e3)') 'tolerance_fraction', tolerance_fraction
      write(u,'(a," = ",i0)') 'max_stage_index', max_stage_index
      write(u,'(a," = ",i0)') 'extrap_order', extrap_order
      write(u,'(a," = ",i0)') 'last_model_unit', last_model_unit
      write(u,'(a," = ",i0)') 'run_log_unit', run_log_unit
      write(u,'(a," = ",i0)') 'terminal_unit', terminal_unit

contains

subroutine path(name, value)
      character(len=*), intent(in) :: name, value
      write(u,'(a," = ",a)') 'job%'//name, '"'//trim(value)//'"'
end subroutine path

function itoa(i) result(s)
      integer, intent(in) :: i
      character(len=:), allocatable :: s
      character(len=12) :: buf
      write(buf,'(i0)') i
      s = trim(buf)
end function itoa

end program test_controls
